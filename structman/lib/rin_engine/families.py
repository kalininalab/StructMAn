"""Interaction detectors, one per family, vectorised over neighbour pairs.

Every detector reports atom pairs between different residues where at least one atom belongs to a
kept chain instance, as (family, atom_a, atom_b, distance, score). Ring interactions report the
first atom of each ring, H-bonds the donor heavy atom. Atom pairs within
`max_excluded_bond_separation` bonds of each other (across a peptide bond, for instance) are never
interactions.

The criteria are in constants.THRESHOLDS. pi-pi stacking uses biotite's implementation; pi-cation is
computed here because it has to use the group-level cation typing of atom_types (biotite's only
sees formal charges, so ligand amines would be missed and metal ions taken for cations).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

import biotite.structure as struc

from structman.lib.rin_engine import constants as C
from structman.lib.rin_engine.atom_types import AtomTypes
from structman.lib.rin_engine.contacts import contact_areas

T = C.THRESHOLDS


@dataclass
class AtomInteractions:
    family: np.ndarray
    atom_a: np.ndarray
    atom_b: np.ndarray
    distance: np.ndarray
    score: np.ndarray
    # Solvent accessible area of every atom (0 for hydrogens), from the tessellation.
    sas_area: np.ndarray = field(default_factory = lambda: np.zeros(0))


class _Collector:
    def __init__(self):
        self.parts = []

    def add(self, family: int, atom_a: np.ndarray, atom_b: np.ndarray, distance: np.ndarray, score: np.ndarray | None = None):
        if len(atom_a) == 0:
            return
        if score is None:
            score = np.ones(len(atom_a), dtype = np.float32)
        self.parts.append((np.full(len(atom_a), family, dtype = np.int8), atom_a.astype(np.int64), atom_b.astype(np.int64),
                           distance.astype(np.float32), score.astype(np.float32)))

    def finish(self, sas_area: np.ndarray) -> AtomInteractions:
        if not self.parts:
            empty_i = np.zeros(0, dtype = np.int64)
            empty_f = np.zeros(0, dtype = np.float32)
            return AtomInteractions(np.zeros(0, dtype = np.int8), empty_i, empty_i, empty_f, empty_f, sas_area)
        columns = list(zip(*self.parts))
        return AtomInteractions(*(np.concatenate(column) for column in columns), sas_area = sas_area)


def _angles(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Angle a-b-c at b in degrees, row-wise."""
    v1 = a - b
    v2 = c - b
    norm = np.linalg.norm(v1, axis = 1) * np.linalg.norm(v2, axis = 1)
    with np.errstate(invalid = 'ignore', divide = 'ignore'):
        cos = np.einsum('ij,ij->i', v1, v2) / norm
    return np.degrees(np.arccos(np.clip(np.nan_to_num(cos), -1.0, 1.0)))


class _Exclusions:
    """Inter-residue heavy atom pairs closer than the maximal excluded bond separation."""

    def __init__(self, atoms, types: AtomTypes):
        n = atoms.array_length()
        self.n = n
        bonds = atoms.bonds.as_array() if atoms.bonds is not None else np.zeros((0, 3), dtype = np.int64)
        bonds = bonds[types.heavy[bonds[:, 0]] & types.heavy[bonds[:, 1]]] if len(bonds) else bonds
        if len(bonds) == 0:
            self.keys = np.zeros(0, dtype = np.int64)
            return
        rows = np.concatenate([bonds[:, 0], bonds[:, 1]]).astype(np.int64)
        cols = np.concatenate([bonds[:, 1], bonds[:, 0]]).astype(np.int64)
        # int32: the products count paths, which a narrower type could wrap to zero.
        adjacency = sparse.csr_matrix((np.ones(len(rows), dtype = np.int32), (rows, cols)), shape = (n, n))
        reach = adjacency.copy()
        power = adjacency
        for _ in range(int(T['max_excluded_bond_separation']) - 1):
            power = power @ adjacency
            reach = reach + power
        reach = reach.tocoo()
        cross = types.node[reach.row] != types.node[reach.col]
        low = np.minimum(reach.row[cross], reach.col[cross]).astype(np.int64)
        high = np.maximum(reach.row[cross], reach.col[cross]).astype(np.int64)
        self.keys = np.unique(low * n + high)

    def contains(self, i: np.ndarray, j: np.ndarray) -> np.ndarray:
        if len(self.keys) == 0 or len(i) == 0:
            return np.zeros(len(i), dtype = bool)
        key = np.minimum(i, j).astype(np.int64) * self.n + np.maximum(i, j)
        position = np.searchsorted(self.keys, key)
        position = np.minimum(position, len(self.keys) - 1)
        return self.keys[position] == key


def _relevant(atoms, types: AtomTypes, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    return (atoms.kept[i] | atoms.kept[j]) & (types.node[i] != types.node[j])


def _unique_pairs(i: np.ndarray, j: np.ndarray, *others: np.ndarray):
    """Keep the first occurrence of every (i, j) pair."""
    if len(i) == 0:
        return (i, j, *others)
    key = i.astype(np.int64) * (int(max(i.max(), j.max())) + 1) + j
    _, first = np.unique(key, return_index = True)
    first.sort()
    return (i[first], j[first], *(other[first] for other in others))


def _cross_pairs(points_a: np.ndarray, points_b: np.ndarray, radius: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(index into a, index into b, distance) of all pairs closer than radius."""
    if len(points_a) == 0 or len(points_b) == 0:
        empty = np.zeros(0, dtype = np.int64)
        return empty, empty, np.zeros(0)
    result = cKDTree(points_a).sparse_distance_matrix(cKDTree(points_b), radius, output_type = 'ndarray')
    return result['i'].astype(np.int64), result['j'].astype(np.int64), result['v']


def detect_contacts(atoms, types: AtomTypes, exclusions: _Exclusions, collector: _Collector) -> np.ndarray:
    heavy = np.nonzero(types.heavy)[0]
    sas_area = np.zeros(atoms.array_length(), dtype = np.float64)
    if len(heavy) < 2:
        return sas_area
    a, b, area, sas = contact_areas(atoms, heavy, T['probe_radius'])
    sas_area[heavy] = sas
    keep = _relevant(atoms, types, a, b) & ~exclusions.contains(a, b) & (area > 0)
    a, b, area = a[keep], b[keep], area[keep]
    distance = np.linalg.norm(atoms.coord[a] - atoms.coord[b], axis = 1)
    collector.add(C.CONTACT, a, b, distance, area)
    return sas_area


def detect_atom_pairs(atoms, types: AtomTypes, exclusions: _Exclusions, collector: _Collector) -> None:
    heavy = np.nonzero(types.heavy)[0]
    if len(heavy) < 2:
        return
    pairs = cKDTree(atoms.coord[heavy]).query_pairs(C.MAX_ATOM_PAIR_DISTANCE, output_type = 'ndarray')
    i = heavy[pairs[:, 0]]
    j = heavy[pairs[:, 1]]
    keep = _relevant(atoms, types, i, j)
    i, j = i[keep], j[keep]
    keep = ~exclusions.contains(i, j)
    i, j = i[keep], j[keep]
    d = np.linalg.norm(atoms.coord[i] - atoms.coord[j], axis = 1)

    element = atoms.element
    is_o = element == 'O'
    is_s = element == 'S'
    coordinating = np.isin(element, ('N', 'O', 'S'))

    polar = ((types.donor[i] & types.acceptor[j]) | (types.acceptor[i] & types.donor[j])) & (d <= T['polar'])
    collector.add(C.POLAR, i[polar], j[polar], d[polar])

    weak_polar = ((types.weak_donor[i] & is_o[j]) | (is_o[i] & types.weak_donor[j])) & (d <= T['weak_polar'])
    collector.add(C.WEAK_POLAR, i[weak_polar], j[weak_polar], d[weak_polar])

    ionic = ((types.cation[i] & types.anion[j]) | (types.anion[i] & types.cation[j])) & (d <= T['ionic'])
    collector.add(C.IONIC, i[ionic], j[ionic], d[ionic])

    hydrophobic = types.hydrophobic[i] & types.hydrophobic[j] & (d <= T['hydrophobic'])
    collector.add(C.HYDROPHOBIC, i[hydrophobic], j[hydrophobic], d[hydrophobic])

    disulfide = is_s[i] & is_s[j] & (d <= T['disulfide'])
    collector.add(C.DISULFIDE, i[disulfide], j[disulfide], d[disulfide])

    metal = ((types.metal[i] & coordinating[j]) | (coordinating[i] & types.metal[j])) & (d <= T['metal'])
    collector.add(C.METAL, i[metal], j[metal], d[metal])

    _halogen_atom_acceptors(atoms, types, i, j, d, collector)


def _halogen_atom_acceptors(atoms, types: AtomTypes, i, j, d, collector: _Collector) -> None:
    forward = types.halogen_donor[i] & types.acceptor[j]
    backward = types.halogen_donor[j] & types.acceptor[i]
    x = np.concatenate([i[forward], j[backward]])
    acceptor = np.concatenate([j[forward], i[backward]])
    distance = np.concatenate([d[forward], d[backward]])
    if len(x) == 0:
        return
    radii = np.array([C.BONDI_RADII.get(e, 1.7) for e in atoms.element[x]]) + \
        np.array([C.BONDI_RADII.get(e, 1.7) for e in atoms.element[acceptor]])
    coord = atoms.coord.astype(np.float64)
    theta1 = _angles(coord[types.halogen_carbon[x]], coord[x], coord[acceptor])
    keep = (distance <= radii) & (theta1 >= T['halogen_theta1'])

    # theta2 = X...A-Y for the heavy neighbours Y of the acceptor; an acceptor without heavy
    # neighbours (an ion) has no such angle to violate.
    neighbors = types.neighbors[acceptor]
    theta2_ok = np.zeros(len(x), dtype = bool)
    has_neighbor = np.zeros(len(x), dtype = bool)
    for column in range(neighbors.shape[1]):
        y = neighbors[:, column]
        valid = (y >= 0) & types.heavy[np.maximum(y, 0)]
        has_neighbor |= valid
        theta2 = _angles(coord[x], coord[acceptor], coord[np.maximum(y, 0)])
        theta2_ok |= valid & (theta2 >= T['halogen_theta2_min']) & (theta2 <= T['halogen_theta2_max'])
    keep &= theta2_ok | ~has_neighbor
    collector.add(C.HALOGEN, x[keep], acceptor[keep], distance[keep])


def detect_hydrogen_bonds(atoms, types: AtomTypes, exclusions: _Exclusions, collector: _Collector) -> None:
    coord = atoms.coord.astype(np.float64)
    hydrogens = np.nonzero(types.parent >= 0)[0]
    parent = types.parent[hydrogens]

    # D-H...A
    polar_h = hydrogens[types.donor[parent]]
    acceptors = np.nonzero(types.acceptor & types.heavy)[0]
    hi, ai, _ = _cross_pairs(coord[polar_h], coord[acceptors], T['hbond_h_acceptor'])
    h, acceptor = polar_h[hi], acceptors[ai]
    donor = types.parent[h]
    keep = (donor != acceptor) & _relevant(atoms, types, donor, acceptor)
    h, acceptor, donor = h[keep], acceptor[keep], donor[keep]
    keep = ~exclusions.contains(donor, acceptor) & (_angles(coord[donor], coord[h], coord[acceptor]) >= T['hbond_angle'])
    donor, acceptor = donor[keep], acceptor[keep]
    donor, acceptor = _unique_pairs(donor, acceptor)
    collector.add(C.HBOND, donor, acceptor, np.linalg.norm(coord[donor] - coord[acceptor], axis = 1))

    # C-H...O
    carbon_h = hydrogens[types.weak_donor[parent]]
    oxygens = np.nonzero(atoms.element == 'O')[0]
    hi, oi, _ = _cross_pairs(coord[carbon_h], coord[oxygens], T['weak_hbond_donor_acceptor'])
    h, oxygen = carbon_h[hi], oxygens[oi]
    carbon = types.parent[h]
    keep = _relevant(atoms, types, carbon, oxygen)
    h, oxygen, carbon = h[keep], oxygen[keep], carbon[keep]
    c_o = np.linalg.norm(coord[carbon] - coord[oxygen], axis = 1)
    keep = (~exclusions.contains(carbon, oxygen) & (c_o <= T['weak_hbond_donor_acceptor'])
            & (_angles(coord[carbon], coord[h], coord[oxygen]) >= T['weak_hbond_angle']))
    carbon, oxygen, c_o = carbon[keep], oxygen[keep], c_o[keep]
    carbon, oxygen, c_o = _unique_pairs(carbon, oxygen, c_o)
    collector.add(C.WEAK_HBOND, carbon, oxygen, c_o)


def detect_ring_interactions(atoms, types: AtomTypes, collector: _Collector) -> None:
    if not types.rings:
        return
    coord = atoms.coord.astype(np.float64)
    ring_atom = np.array([ring[0] for ring in types.rings], dtype = np.int64)
    centroid = types.ring_centroid
    normal = types.ring_normal

    # pi-pi stacking (biotite), parallel or T-shaped
    stacks = struc.find_stacking_interactions(
        atoms, centroid_cutoff = T['stacking_centroid'],
        plane_angle_tol = np.deg2rad(T['stacking_plane_tolerance']),
        shift_angle_tol = np.deg2rad(T['stacking_shift_tolerance']),
    )
    if stacks:
        a = np.array([ring_a[0] for ring_a, _, _ in stacks], dtype = np.int64)
        b = np.array([ring_b[0] for _, ring_b, _ in stacks], dtype = np.int64)
        parallel = np.array([kind == struc.PiStacking.PARALLEL for _, _, kind in stacks], dtype = bool)
        centroid_a = np.array([coord[ring_a].mean(axis = 0) for ring_a, _, _ in stacks])
        centroid_b = np.array([coord[ring_b].mean(axis = 0) for _, ring_b, _ in stacks])
        distance = np.linalg.norm(centroid_a - centroid_b, axis = 1)
        keep = _relevant(atoms, types, a, b)
        collector.add(C.PIPI_PARALLEL, a[keep & parallel], b[keep & parallel], distance[keep & parallel])
        collector.add(C.PIPI_T, a[keep & ~parallel], b[keep & ~parallel], distance[keep & ~parallel])

    # pi-cation: cation within the cutoff of the centroid, close to the ring normal
    cations = np.nonzero(types.cation)[0]
    ri, ci, distance = _cross_pairs(centroid, coord[cations], T['pication'])
    if len(ri):
        cation = cations[ci]
        offset = coord[cation] - centroid[ri]
        cos = np.abs(np.einsum('ij,ij->i', offset, normal[ri])) / np.linalg.norm(offset, axis = 1)
        keep = (np.degrees(np.arccos(np.clip(cos, 0.0, 1.0))) <= T['pication_angle_tolerance'])
        keep &= _relevant(atoms, types, ring_atom[ri], cation)
        collector.add(C.PICATION, ring_atom[ri][keep], cation[keep], distance[keep])

    # pi-hydrogen: X-H pointing at the ring face
    hydrogens = np.nonzero(types.parent >= 0)[0]
    ri, hi, distance = _cross_pairs(centroid, coord[hydrogens], T['pihydrogen_h_centroid'])
    if len(ri):
        h = hydrogens[hi]
        x = types.parent[h]
        offset = coord[h] - centroid[ri]
        cos = np.abs(np.einsum('ij,ij->i', offset, normal[ri])) / np.linalg.norm(offset, axis = 1)
        keep = np.degrees(np.arccos(np.clip(cos, 0.0, 1.0))) <= T['pihydrogen_normal_angle']
        keep &= _angles(coord[x], coord[h], centroid[ri]) >= T['pihydrogen_xh_centroid_angle']
        keep &= _relevant(atoms, types, ring_atom[ri], x)
        ring, x, distance = _unique_pairs(ring_atom[ri][keep], x[keep], distance[keep])
        collector.add(C.PIHYDROGEN, ring, x, distance)

    # halogen bond with a pi system as acceptor
    halogens = np.nonzero(types.halogen_donor)[0]
    ri, xi, distance = _cross_pairs(centroid, coord[halogens], T['halogen_pi_centroid'])
    if len(ri):
        x = halogens[xi]
        theta1 = _angles(coord[types.halogen_carbon[x]], coord[x], centroid[ri])
        keep = (theta1 >= T['halogen_theta1']) & _relevant(atoms, types, ring_atom[ri], x)
        collector.add(C.HALOGEN, x[keep], ring_atom[ri][keep], distance[keep])


def detect_all(atoms, types: AtomTypes) -> AtomInteractions:
    collector = _Collector()
    exclusions = _Exclusions(atoms, types)
    sas_area = detect_contacts(atoms, types, exclusions, collector)
    detect_atom_pairs(atoms, types, exclusions, collector)
    detect_hydrogen_bonds(atoms, types, exclusions, collector)
    detect_ring_interactions(atoms, types, collector)
    return collector.finish(sas_area)
