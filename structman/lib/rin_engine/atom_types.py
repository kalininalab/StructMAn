"""What every atom can take part in, derived from element, bonds and formal charge.

Everything is computed for all atoms at once from the bond table (biotite `get_all_bonds`), so
ligands, nucleotides and modified residues are typed by the same rules as amino acids:

- donor: N, O or S carrying a hydrogen (hydrogens are explicit, placed by hydride);
- acceptor: any O; N without hydrogen, fewer than three heavy neighbours and not cationic;
  thioether or thiolate S;
- weak donor: C carrying a hydrogen;
- cationic / anionic: formal charge, extended to whole charged groups at pH 7 - both carboxylate
  oxygens, all terminal phosphate / sulfate oxygens, all guanidinium / amidinium nitrogens,
  aliphatic amines, quaternary nitrogen. Metals are not counted, they have their own family;
- hydrophobic (Arpeggio): C not bonded to N, O or F; thioether S; Cl, Br, I;
- halogen bond donor: Cl, Br or I bonded to C;
- metal: element in constants.METAL_ELEMENTS, whatever residue it is in (heme iron included).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import biotite.structure as struc

from structman.lib.rin_engine import constants as C

_DOUBLE_BOND_TYPES = (2, 6)                 # DOUBLE, AROMATIC_DOUBLE
_MULTIPLE_BOND_TYPES = (2, 3, 4, 6, 7)      # anything above a single bond
_AROMATIC_BOND_TYPES = (5, 6, 7, 9)


@dataclass
class AtomTypes:
    node: np.ndarray            # residue index of every atom
    heavy: np.ndarray
    parent: np.ndarray          # hydrogens: the bonded heavy atom; -1 otherwise
    neighbors: np.ndarray       # (n, k) bonded atoms, -1 padded
    donor: np.ndarray
    acceptor: np.ndarray
    weak_donor: np.ndarray
    cation: np.ndarray
    anion: np.ndarray
    hydrophobic: np.ndarray
    halogen_donor: np.ndarray
    halogen_carbon: np.ndarray  # the C bonded to a halogen donor; -1 otherwise
    metal: np.ndarray
    atom_class: np.ndarray      # C.MAIN_CHAIN ... C.OTHER
    rings: list                 # aromatic rings, arrays of atom indices
    ring_centroid: np.ndarray   # (n_rings, 3)
    ring_normal: np.ndarray     # (n_rings, 3), unit length


def _residue_index(atoms) -> np.ndarray:
    starts = struc.get_residue_starts(atoms)
    index = np.zeros(atoms.array_length(), dtype = np.int64)
    index[starts[1:]] = 1
    return np.cumsum(index)


def ring_geometry(coord: np.ndarray, rings: list) -> tuple[np.ndarray, np.ndarray]:
    centroids = np.zeros((len(rings), 3), dtype = np.float64)
    normals = np.zeros((len(rings), 3), dtype = np.float64)
    for pos, ring in enumerate(rings):
        points = coord[ring].astype(np.float64)
        centroid = points.mean(axis = 0)
        # Normal of the best plane through the ring atoms.
        _, _, vt = np.linalg.svd(points - centroid)
        centroids[pos] = centroid
        normals[pos] = vt[2] / np.linalg.norm(vt[2])
    return centroids, normals


def type_atoms(atoms) -> AtomTypes:
    n = atoms.array_length()
    element = atoms.element
    is_h = element == 'H'
    heavy = ~is_h
    charge = atoms.charge.astype(np.int64)

    neighbors, bond_types = atoms.bonds.get_all_bonds()
    if neighbors.shape[1] == 0:
        neighbors = np.full((n, 1), -1, dtype = np.int64)
        bond_types = np.zeros((n, 1), dtype = np.int64)
    valid = neighbors >= 0
    safe = np.where(valid, neighbors, 0)
    neighbor_element = np.where(valid, element[safe], '')
    neighbor_is_h = valid & (neighbor_element == 'H')
    neighbor_heavy = valid & ~neighbor_is_h

    n_hydrogens = neighbor_is_h.sum(axis = 1)
    n_heavy_neighbors = neighbor_heavy.sum(axis = 1)

    parent = np.full(n, -1, dtype = np.int64)
    first_heavy = np.argmax(neighbor_heavy, axis = 1)
    has_heavy = neighbor_heavy.any(axis = 1)
    parent[is_h & has_heavy] = safe[np.nonzero(is_h & has_heavy)[0], first_heavy[is_h & has_heavy]]

    is_c = element == 'C'
    is_n = element == 'N'
    is_o = element == 'O'
    is_s = element == 'S'
    metal = np.isin(element, list(C.METAL_ELEMENTS))

    bonded_to_s = (valid & (neighbor_element == 'S')).any(axis = 1)
    bonded_to_nof = (valid & np.isin(neighbor_element, ('N', 'O', 'F'))).any(axis = 1)
    only_c_heavy_neighbors = (neighbor_heavy & (neighbor_element != 'C')).sum(axis = 1) == 0

    donor = heavy & (is_n | is_o | is_s) & (n_hydrogens > 0)
    weak_donor = is_c & (n_hydrogens > 0)
    thioether = is_s & (n_hydrogens == 0) & (n_heavy_neighbors == 2) & ~bonded_to_s
    acceptor = (
        is_o
        | (is_n & (n_hydrogens == 0) & (n_heavy_neighbors < 3) & (charge <= 0))
        | thioether
        | (is_s & (charge < 0))
    )

    cation, anion = _charged_groups(atoms, element, charge, neighbors, bond_types, valid, safe,
                                    neighbor_element, n_hydrogens, n_heavy_neighbors, metal)

    hydrophobic = (
        (is_c & ~bonded_to_nof & (charge == 0))
        | (thioether & only_c_heavy_neighbors)
        | (np.isin(element, ('CL', 'BR', 'I')) & (charge == 0))
    )

    halogen_donor = np.zeros(n, dtype = bool)
    halogen_carbon = np.full(n, -1, dtype = np.int64)
    halogens = np.nonzero(np.isin(element, list(C.HALOGEN_DONOR_ELEMENTS)))[0]
    for atom in halogens:
        carbons = [b for b in neighbors[atom] if b >= 0 and element[b] == 'C']
        if carbons:
            halogen_donor[atom] = True
            halogen_carbon[atom] = carbons[0]

    node = _residue_index(atoms)
    atom_class = _atom_classes(atoms, is_h, parent)

    rings = struc.find_aromatic_rings(atoms) if n else []
    rings = [np.asarray(ring, dtype = np.int64) for ring in rings]
    ring_centroid, ring_normal = ring_geometry(atoms.coord, rings)

    return AtomTypes(
        node = node, heavy = heavy, parent = parent, neighbors = neighbors,
        donor = donor, acceptor = acceptor, weak_donor = weak_donor,
        cation = cation, anion = anion, hydrophobic = hydrophobic,
        halogen_donor = halogen_donor, halogen_carbon = halogen_carbon, metal = metal,
        atom_class = atom_class, rings = rings, ring_centroid = ring_centroid, ring_normal = ring_normal,
    )


def _charged_groups(atoms, element, charge, neighbors, bond_types, valid, safe, neighbor_element,
                    n_hydrogens, n_heavy_neighbors, metal) -> tuple[np.ndarray, np.ndarray]:
    n = atoms.array_length()
    cation = (charge > 0) & ~metal
    anion = (charge < 0) & ~metal

    # Canonical amino acids by atom name: the charged group, not only the atom that carries the
    # formal charge.
    for res_name, names in C.AMINO_ACID_CATIONIC.items():
        cation |= (atoms.res_name == res_name) & np.isin(atoms.atom_name, names)
    for res_name, names in C.AMINO_ACID_ANIONIC.items():
        anion |= (atoms.res_name == res_name) & np.isin(atoms.atom_name, names)

    heavy_neighbor = valid & (neighbor_element != 'H')
    is_o = element == 'O'
    is_n = element == 'N'
    # Terminal oxygen: bonded to exactly one heavy atom.
    terminal_o = is_o & (n_heavy_neighbors == 1)
    first_heavy = safe[np.arange(n), np.argmax(heavy_neighbor, axis = 1)]

    # Count of terminal oxygens on every atom, used to find carboxyl / phosphate / sulfate centres.
    terminal_o_count = np.zeros(n, dtype = np.int64)
    np.add.at(terminal_o_count, first_heavy[terminal_o], 1)
    o_count = (valid & (neighbor_element == 'O')).sum(axis = 1)

    centre = first_heavy
    centre_element = element[centre]
    # Carboxylate (incl. C-terminus): C with two terminal oxygens and no N neighbour (no amide).
    centre_has_n = (valid[centre] & (neighbor_element[centre] == 'N')).any(axis = 1)
    carboxyl = terminal_o & (centre_element == 'C') & (terminal_o_count[centre] == 2) & ~centre_has_n
    # Phosphate / phosphonate, sulfate / sulfonate: terminal oxygens on P or S with >= 3 oxygens.
    oxo_acid = terminal_o & np.isin(centre_element, ('P', 'S')) & (o_count[centre] >= 3)
    anion |= carboxyl | oxo_acid

    # Nitrogen: multiple bonds and aromaticity decide between amine, amide, amidine/guanidine,
    # imine and nitrile.
    double_bonded = (valid & np.isin(bond_types, _MULTIPLE_BOND_TYPES)).any(axis = 1)
    aromatic = (valid & np.isin(bond_types, _AROMATIC_BOND_TYPES)).any(axis = 1)
    # Neighbouring carbon with a double bond to O, S or N (amide, thioamide, amidine) or aromatic.
    conjugated_neighbor = np.zeros(n, dtype = bool)
    carbonyl_like = (element == 'C') & (
        valid & np.isin(bond_types, _DOUBLE_BOND_TYPES) & np.isin(neighbor_element, ('O', 'S', 'N'))
    ).any(axis = 1)
    for column in range(neighbors.shape[1]):
        partner = safe[:, column]
        conjugated_neighbor |= valid[:, column] & (carbonyl_like[partner] | aromatic[partner])

    only_c_or_h = ((valid & ~np.isin(neighbor_element, ('C', 'H'))).sum(axis = 1)) == 0
    amine = is_n & ~double_bonded & ~aromatic & ~conjugated_neighbor & only_c_or_h & (n_heavy_neighbors <= 3)
    quaternary = is_n & (n_heavy_neighbors == 4)
    # Guanidine / amidine carbon: three nitrogen neighbours, or two with a C=N double bond.
    n_neighbors = (valid & (neighbor_element == 'N')).sum(axis = 1)
    c_double_n = (element == 'C') & (valid & np.isin(bond_types, _DOUBLE_BOND_TYPES) & (neighbor_element == 'N')).any(axis = 1)
    amidinium_c = (element == 'C') & ~aromatic & ((n_neighbors == 3) | ((n_neighbors == 2) & c_double_n))
    amidinium_n = np.zeros(n, dtype = bool)
    for column in range(neighbors.shape[1]):
        partner = safe[:, column]
        amidinium_n |= is_n & valid[:, column] & amidinium_c[partner]

    is_amino_acid = np.isin(atoms.res_name, list(C.CANONICAL_AMINO_ACIDS))
    # Amino acids are covered by the tables and the estimated charges (backbone N is an amide,
    # the N-terminus has its charge); the generic rules are for everything else.
    generic = ~is_amino_acid
    cation |= generic & (amine | quaternary | amidinium_n) & ~metal
    return cation, anion & ~metal


def _atom_classes(atoms, is_h: np.ndarray, parent: np.ndarray) -> np.ndarray:
    n = atoms.array_length()
    classes = np.full(n, C.OTHER, dtype = np.int8)
    protein = atoms.kind == C.KIND_PROTEIN
    nucleotide = atoms.kind == C.KIND_NUCLEOTIDE

    names = atoms.atom_name
    backbone = np.isin(names, list(C.AMINO_ACID_BACKBONE))
    nucleic_backbone = np.isin(names, list(C.NUCLEOTIDE_BACKBONE))
    classes[protein] = np.where(backbone[protein], C.MAIN_CHAIN, C.SIDE_CHAIN)
    classes[nucleotide] = np.where(nucleic_backbone[nucleotide], C.NUCLEIC_BACKBONE, C.NUCLEIC_BASE)

    # Hydrogens share the class of the atom they sit on.
    hydrogens = np.nonzero(is_h & (parent >= 0))[0]
    classes[hydrogens] = classes[parent[hydrogens]]
    return classes
