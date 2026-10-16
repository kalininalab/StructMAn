"""Residue interaction networks (RINs) of PDB mmCIF biological assemblies.

Replaces RINerator (MolProbity Reduce + Probe on PDB format files). Built from open components:

- biotite: mmCIF parsing, CCD bonds and charges, aromatic rings, pi-stacking;
- hydride: hydrogen placement;
- Voronota-LT: atom-atom contact areas (the contact score);
- typed interactions after the published RING 4.0 and Arpeggio criteria (constants.THRESHOLDS):
  H-bonds, weak H-bonds, polar / weak polar contacts, ionic, pi-pi (parallel / T), pi-cation,
  pi-hydrogen, disulfide, metal coordination, halogen bonds, hydrophobic contacts.

Interactions are computed where at least one partner belongs to a kept chain instance (one
instance per chain of the deposited entry, as in the search DB), see structure.py. The result of
compute_assembly_rin is the document storage.write_rin stores.

Worker processes import this package; keep it free of structman.base_utils.base_utils.
"""

from __future__ import annotations

import os
import time

import numpy as np

from structman.lib import mmcifParser
from structman.lib.rin_engine import atom_types, constants, families, metrics, network, storage, structure
from structman.lib.rin_engine.constants import ENGINE_VERSION


def structure_id_of(path: str) -> str:
    parsed = mmcifParser.parse_assembly_filename(os.path.basename(path))
    if parsed is None:
        # Not an assembly file name (a model, a test file): the file name without extensions.
        name = os.path.basename(path)
        for suffix in ('.gz', '.cif', '.bcif', '.pdb'):
            if name.endswith(suffix):
                name = name[:-len(suffix)]
        return name
    return mmcifParser.make_assembly_id(*parsed)


def compute_assembly_rin(path: str, add_hydrogens: bool = True) -> dict:
    """The RIN of one assembly file as a storage document."""
    stat = os.stat(path)
    return compute_rin(path = path, structure_id = structure_id_of(path), add_hydrogens = add_hydrogens,
                       source = {'name': os.path.basename(path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns})


def compute_rin(path: str | None = None, cif_text: str | None = None, structure_id: str | None = None,
                source: dict | None = None, add_hydrogens: bool = True) -> dict:
    """The RIN of a structure given as a file or as mmCIF text, as a storage document.

    `source` describes where it came from (file stamp, or a model's stamp in its SQL table); it is
    stored as given and used by the builders to decide what is up to date.
    """
    timings = {}
    clock = time.perf_counter()

    if cif_text is not None:
        import io
        from biotite.structure.io import pdbx
        block = pdbx.CIFFile.read(io.StringIO(cif_text)).block
        prepared = structure.prepare_assembly(block = block, add_hydrogens = add_hydrogens)
    else:
        prepared = structure.prepare_assembly(path, add_hydrogens = add_hydrogens)
    if structure_id is None:
        structure_id = structure_id_of(path) if path else ''
    timings['prepare'] = time.perf_counter() - clock
    clock = time.perf_counter()

    types = atom_types.type_atoms(prepared.atoms)
    timings['typing'] = time.perf_counter() - clock
    clock = time.perf_counter()

    interactions = families.detect_all(prepared.atoms, types)
    timings['interactions'] = time.perf_counter() - clock
    clock = time.perf_counter()

    net = network.build_network(prepared.atoms, types, interactions)
    timings['network'] = time.perf_counter() - clock
    clock = time.perf_counter()

    centralities, centrality_mode = metrics.betweenness(net)
    flows, flow_mode = metrics.flow(net)
    timings['metrics'] = time.perf_counter() - clock

    family_counts = np.bincount(net.edge_family.astype(np.int64), minlength = len(constants.FAMILY_NAMES))
    stats = dict(prepared.stats)
    stats.update({
        'nodes': net.n_nodes,
        'kept_nodes': int(net.kept.sum()),
        'edges': int(len(net.edge_a)),
        'edges_per_family': {name: int(count) for name, count in zip(constants.FAMILY_NAMES, family_counts)},
        'rings': len(types.rings),
        'seconds': {name: round(value, 3) for name, value in timings.items()},
    })

    return {
        'format': storage.FORMAT_VERSION,
        'engine_version': ENGINE_VERSION,
        'structure_id': structure_id,
        'source': dict(source or {}),
        'params': dict(constants.THRESHOLDS),
        'families': list(constants.FAMILY_NAMES),
        'classes': list(constants.CLASS_NAMES),
        'kinds': list(constants.KIND_NAMES),
        'mode': {
            'region_restricted': bool(prepared.stats.get('region_restricted', False)),
            'centrality': centrality_mode,
            'flow': flow_mode,
            'hydrogens': bool(add_hydrogens),
        },
        'stats': stats,
        'nodes': {
            'chain': net.chain, 'res_id': net.res_id, 'ins_code': net.ins_code, 'res_name': net.res_name,
            'label_asym_id': net.label_asym_id, 'label_seq_id': net.label_seq_id, 'entity_id': net.entity_id,
            'kind': net.kind, 'kept': net.kept, 'sas_area': net.sas_area, 'b_factor': net.b_factor,
        },
        'edges': {
            'a': net.edge_a, 'b': net.edge_b, 'family': net.edge_family,
            'class_a': net.edge_class_a, 'class_b': net.edge_class_b,
            'count': net.edge_count, 'score': net.edge_score, 'min_distance': net.edge_min_distance,
        },
        'metrics': {
            'names': list(metrics.CENTRALITY_NAMES) + list(metrics.FLOW_NAMES),
            'values': np.concatenate([centralities, flows], axis = 1),
        },
    }
