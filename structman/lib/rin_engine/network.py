"""From atom-level interactions to the residue interaction network.

Nodes are residues (polymer residues, ligands, metals, ions, sugars), identified like the runtime
will identify them: author chain id (with the '-N' operator suffix of symmetry copies), author
residue number, insertion code, residue name - plus the label ids, the entity, the node kind and
whether it belongs to a kept chain instance.

Edges are undirected residue pairs (node_a < node_b) per family and per pair of atom classes
(class_a belongs to node_a), with the number of atom (or ring) pairs, the summed score (contact
area for CONTACT, the count for the other families) and the shortest distance.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import biotite.structure as struc

from structman.lib.rin_engine import constants as C
from structman.lib.rin_engine.atom_types import AtomTypes
from structman.lib.rin_engine.families import AtomInteractions

_N_FAMILIES = len(C.FAMILY_NAMES)
_N_CLASSES = len(C.CLASS_NAMES)


@dataclass
class ResidueNetwork:
    # nodes
    chain: np.ndarray
    res_id: np.ndarray
    ins_code: np.ndarray
    res_name: np.ndarray
    label_asym_id: np.ndarray
    label_seq_id: np.ndarray        # -1 for residues without one (non-polymers)
    entity_id: np.ndarray
    kind: np.ndarray
    kept: np.ndarray
    sas_area: np.ndarray
    b_factor: np.ndarray            # mean B of the heavy atoms (pLDDT in AlphaFold models)
    # edges
    edge_a: np.ndarray
    edge_b: np.ndarray
    edge_family: np.ndarray
    edge_class_a: np.ndarray
    edge_class_b: np.ndarray
    edge_count: np.ndarray
    edge_score: np.ndarray
    edge_min_distance: np.ndarray

    @property
    def n_nodes(self) -> int:
        return len(self.res_id)

    def contact_graph(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(node_a, node_b, total contact area) per residue pair in contact."""
        contact = self.edge_family == C.CONTACT
        a = self.edge_a[contact].astype(np.int64)
        b = self.edge_b[contact].astype(np.int64)
        if len(a) == 0:
            empty = np.zeros(0, dtype = np.int64)
            return empty, empty, np.zeros(0)
        key = a * self.n_nodes + b
        unique, inverse = np.unique(key, return_inverse = True)
        area = np.bincount(inverse, weights = self.edge_score[contact].astype(np.float64))
        return unique // self.n_nodes, unique % self.n_nodes, area


def _label_seq_ids(values: np.ndarray) -> np.ndarray:
    out = np.full(len(values), -1, dtype = np.int32)
    for pos, value in enumerate(values):
        if value not in ('', '.', '?'):
            try:
                out[pos] = int(value)
            except ValueError:
                pass
    return out


def build_network(atoms, types: AtomTypes, interactions: AtomInteractions) -> ResidueNetwork:
    starts = struc.get_residue_starts(atoms)
    n_nodes = len(starts)
    sas = np.bincount(types.node, weights = interactions.sas_area, minlength = n_nodes) if n_nodes else np.zeros(0)
    if hasattr(atoms, 'b_factor') and n_nodes:
        heavy = types.heavy
        b_sum = np.bincount(types.node[heavy], weights = atoms.b_factor[heavy].astype(np.float64), minlength = n_nodes)
        b_count = np.bincount(types.node[heavy], minlength = n_nodes)
        b_factor = np.where(b_count > 0, b_sum / np.maximum(b_count, 1), np.nan)
    else:
        b_factor = np.full(n_nodes, np.nan)

    node_a = types.node[interactions.atom_a]
    node_b = types.node[interactions.atom_b]
    class_a = types.atom_class[interactions.atom_a].astype(np.int64)
    class_b = types.atom_class[interactions.atom_b].astype(np.int64)
    swap = node_a > node_b
    na = np.where(swap, node_b, node_a).astype(np.int64)
    nb = np.where(swap, node_a, node_b).astype(np.int64)
    ca = np.where(swap, class_b, class_a)
    cb = np.where(swap, class_a, class_b)
    family = interactions.family.astype(np.int64)

    key = (((na * n_nodes + nb) * _N_FAMILIES + family) * _N_CLASSES + ca) * _N_CLASSES + cb
    if len(key):
        unique, inverse = np.unique(key, return_inverse = True)
        count = np.bincount(inverse).astype(np.int32)
        score = np.bincount(inverse, weights = interactions.score.astype(np.float64)).astype(np.float32)
        min_distance = np.full(len(unique), np.inf, dtype = np.float32)
        np.minimum.at(min_distance, inverse, interactions.distance)
    else:
        unique = np.zeros(0, dtype = np.int64)
        count = np.zeros(0, dtype = np.int32)
        score = np.zeros(0, dtype = np.float32)
        min_distance = np.zeros(0, dtype = np.float32)

    cb_out = unique % _N_CLASSES
    rest = unique // _N_CLASSES
    ca_out = rest % _N_CLASSES
    rest //= _N_CLASSES
    family_out = rest % _N_FAMILIES
    rest //= _N_FAMILIES
    nb_out = rest % n_nodes if n_nodes else rest
    na_out = rest // n_nodes if n_nodes else rest

    return ResidueNetwork(
        chain = atoms.chain_id[starts].astype(str),
        res_id = atoms.res_id[starts].astype(np.int32),
        ins_code = atoms.ins_code[starts].astype(str),
        res_name = atoms.res_name[starts].astype(str),
        label_asym_id = atoms.label_asym_id[starts].astype(str),
        label_seq_id = _label_seq_ids(atoms.label_seq_id[starts]),
        entity_id = atoms.label_entity_id[starts].astype(str),
        kind = atoms.kind[starts].astype(np.int8),
        kept = atoms.kept[starts].astype(bool),
        sas_area = sas.astype(np.float32),
        b_factor = b_factor.astype(np.float32),
        edge_a = na_out.astype(np.int32),
        edge_b = nb_out.astype(np.int32),
        edge_family = family_out.astype(np.int8),
        edge_class_a = ca_out.astype(np.int8),
        edge_class_b = cb_out.astype(np.int8),
        edge_count = count,
        edge_score = score,
        edge_min_distance = min_distance,
    )
