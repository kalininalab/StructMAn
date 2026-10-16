"""Atom-atom contact areas from a radical (Laguerre-Voronoi) tessellation, by Voronota-LT.

Every pair of heavy atoms that shares a face of the tessellation within the solvent accessible
surface is in contact; the face area is the contact score. It replaces the dot counts of
MolProbity Probe: defined for any pair of atoms, hydrogens not required, and summed up per
residue pair it measures how much two residues touch.

Radii are united-atom ProtOr radii (biotite.structure.info.vdw_radius_protor) - they account for
the hydrogens, which are left out of the tessellation - and single-element van der Waals radii for
anything ProtOr does not cover.
"""

from __future__ import annotations

import functools

import numpy as np

import biotite.structure.info as info

_FALLBACK_RADIUS = 1.8


@functools.lru_cache(maxsize = None)
def _radius(res_name: str, atom_name: str, element: str) -> float:
    try:
        radius = info.vdw_radius_protor(res_name, atom_name)
    except Exception:
        radius = None
    if radius is None:
        radius = info.vdw_radius_single(element)
    return _FALLBACK_RADIUS if radius is None else float(radius)


def atom_radii(atoms, indices: np.ndarray) -> np.ndarray:
    return np.array(
        [_radius(r, a, e) for r, a, e in zip(atoms.res_name[indices], atoms.atom_name[indices], atoms.element[indices])],
        dtype = np.float64,
    )


def contact_areas(atoms, indices: np.ndarray, probe: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Contacts between the given atoms: (atom_a, atom_b, area) in atom indices, and the solvent
    accessible area of every given atom."""
    import voronotalt

    coords = atoms.coord[indices].astype(np.float64)
    radii = atom_radii(atoms, indices)
    balls = [voronotalt.Ball(x, y, z, r) for (x, y, z), r in zip(coords.tolist(), radii.tolist())]
    tessellation = voronotalt.RadicalTessellation(balls, probe = probe)

    table = tessellation.contacts.to_pandas()
    local_a = table['index_a'].to_numpy(dtype = np.int64)
    local_b = table['index_b'].to_numpy(dtype = np.int64)
    area = table['area'].to_numpy(dtype = np.float64)
    sas_area = np.fromiter((cell.sas_area for cell in tessellation.cells), dtype = np.float64, count = len(indices))
    return indices[local_a], indices[local_b], area, sas_area
