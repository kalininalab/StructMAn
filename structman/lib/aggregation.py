import os
import ray
from ray.util.state import summarize_objects
import time
import psutil
import threading
import gc
import math
import numpy as np
import sys
from collections import OrderedDict, deque
from typing import NamedTuple
from itertools import cycle
from numba import njit
from numba.typed import List as NumbaList

from ray.util.queue import Queue

from structman.scripts.createMMDB import mm_lookup, mmdb_entry_path_readonly
from structman.base_utils.base_utils import (
    pack, pack_frame, unpack, aggregate_times, print_times, 
    add_to_times, reset_logger_for_remotes, is_alphafold_model,
    deep_get_size_of, sizeof_fmt, sequence_from_store
    )
from structman.lib.lib_utils import fetch_msa
from structman.base_utils.config_class import Config
from structman.base_utils.ray_utils import (parallel_db_access, ray_init, memory_load_percent,
                                           process_memory_breakdown, read_object_store_stats)
from structman.lib.sdsc import mappings as mappings_package
from structman.lib.sdsc import interface as interface_package
from structman.lib.sdsc import residue as residue_package
from structman.lib.sdsc import protein as protein_package
from structman.lib.rin import locate_code, simple_class_code_kernel as _simple_class_code, CLASS_SCORE_INDICES
from structman.lib.sdsc.sdsc_utils import get_shortest_distances
from structman.lib.database.insertion_lib import remote_insert_interfaces, remote_insert_classifications
from structman.lib.database.retrieval import getStoredResidues, get_alignment_data
from structman.lib.database.retrieval_2 import getBackmaps

def qualityScore(resolution: float, coverage: float, seq_id: float, resolution_wf=0.25, coverage_wf=0.5, seq_id_wf=1.0, pLDDT:float | None = None) -> float:
    """Structure quality score.

    No longer called on the aggregation path - _collect_mappings_for_position
    reproduces it inline, down to the float widths and CPython's compensated sum -
    but kept as the readable definition it reproduces, and as the reference the
    differential test for that kernel checks against.
    """
    if seq_id > 1.0:
        seq_id = seq_id / 100.0
    if coverage > 1.0:
        coverage = coverage / 100.0

    if pLDDT is None or np.isnan(pLDDT):
        pLDDT = 1.0

    # project the criteria to [0,1] via logistic regression
    resolution_value = (1 + math.exp((1.5 * resolution) - 4))**(-1)

    seq_value = (1 + math.exp(10 * (0.4 - seq_id)))**(-1)

    ws = sum((resolution_wf, coverage_wf, seq_id_wf))
    quality = pLDDT * (sum((resolution_value * resolution_wf, coverage * coverage_wf, seq_value * seq_id_wf)) / ws)
    return quality


def _gather_alignment_columns(alignments_for_protein: dict, annotation_list: list[int]):
    """The alignment map's three fields as parallel lists, indexed by anno_index.

    annotation_list is that dict's own key list, so the lookup cannot miss.
    """
    sub_infos = []
    covs = []
    seq_ids = []
    for struct_db_id in annotation_list:
        alignment, cov, seq_id = alignments_for_protein[struct_db_id]
        sub_infos.append(alignment)
        covs.append(cov)
        seq_ids.append(seq_id)
    return sub_infos, covs, seq_ids


def _flatten_alignments(sub_infos: list) -> tuple[np.ndarray, np.ndarray]:
    """The ragged per-structure alignments as one array plus (start, length) rows.

    Same concat-and-offsets shape the backmap store already uses. The kernel needs
    the length because `sub_infos[pos]` used to raise IndexError for a position past
    the end of an alignment, which is one of the skip branches.
    """
    if len(sub_infos) == 0:
        return np.empty(0, dtype=np.int64), np.zeros((0, 2), dtype=np.int64)

    parts = [np.asarray(alignment) for alignment in sub_infos]
    bounds = np.zeros((len(parts), 2), dtype=np.int64)
    offset = 0
    for index, part in enumerate(parts):
        bounds[index, 0] = offset
        bounds[index, 1] = part.shape[0]
        offset += part.shape[0]
    # No dtype forced: np.concatenate settles on the alignments' common type and the
    # kernel specialises on it, rather than every protein paying an int64 widening.
    return np.concatenate(parts), bounds


def _alignment_value_columns(covs: list, seq_ids: list):
    """cov/seq_id as float64 columns, plus the None masks the kernel branches on.

    None and NaN are not the same thing here: a None coverage is what made
    qualityScore raise TypeError and fall back to (0.01, 3.0), while a NaN one
    would have propagated into the score. Both are possible in principle, so the
    None-ness travels separately from the value.

    A None value itself becomes NaN in the column, which is what weight_all's
    `np.asarray(covs, dtype=np.float64)` already turned it into; the one place that
    now sees NaN instead of None is the cov stored on a mapping, which only
    set_recommended_residues reads back.
    """
    n = len(covs)
    cov_column = np.full(n, np.nan, dtype=np.float64)
    seq_id_column = np.full(n, np.nan, dtype=np.float64)
    cov_is_none = np.zeros(n, dtype=np.uint8)
    seq_id_is_none = np.zeros(n, dtype=np.uint8)
    for index in range(n):
        cov = covs[index]
        if cov is None:
            cov_is_none[index] = 1
        else:
            cov_column[index] = cov
        seq_id = seq_ids[index]
        if seq_id is None:
            seq_id_is_none[index] = 1
        else:
            seq_id_column[index] = seq_id
    return cov_column, seq_id_column, cov_is_none, seq_id_is_none


# Passed into the mapping kernel rather than written as a literal, so the pow stays
# a pow. See the comment at its use.
_INVERSE_POWER: float = -1.0


@njit(cache=True)
def _collect_mappings_for_position(
        pos, aa1_code, has_aa1_code,
        sub_infos_concat, sub_infos_bounds,
        cov_column, seq_id_column, cov_is_none, seq_id_is_none,
        store_index_column, res_left_column, res_right_column, resolution_column,
        res_stores, aa_column, plddt_column, sc_rsa_column,
        profile_left, profile_right, class_score_indices,
        surface_threshold, buried_threshold, inverse_power,
        out_anno, out_res_contig, out_row, out_store,
        out_qual, out_seq_id, out_cov, out_aa_id, out_resolution, out_res_aa,
        out_simple_class, out_quality_measure_seen,
        out_no_seq_id, out_bad_contig, counts):
    """One position's worth of the per-mapping aggregation loop.

    Everything the old python loop did per mapping happens here: the alignment
    lookup, the four skip branches, the residue row address, the simple class, and
    qualityScore. Nothing is allocated per mapping - the caller's output columns are
    filled and the row indices come back for a single gather.

    qualityScore is reproduced rather than called, including its float widths: the
    python version multiplies a float32 pLDDT by a python float, which NEP 50
    evaluates in float32, and falls back to (0.01, 3.0) when the coverage is None
    and the multiplication would have raised TypeError.

    Returns (n_mappings, count_continues, n_no_seq_id, n_bad_contig).
    """
    n_annotations = sub_infos_bounds.shape[0]
    n_mappings = 0
    count_continues = 0
    n_no_seq_id = 0
    n_bad_contig = 0

    for anno_index in range(n_annotations):
        start = sub_infos_bounds[anno_index, 0]
        length = sub_infos_bounds[anno_index, 1]

        # `sub_infos[pos]` raised IndexError for a position past the alignment.
        if pos >= length:
            count_continues += 1
            counts[1] += 1
            continue

        sub_info = sub_infos_concat[start + pos]
        if sub_info == -1:
            count_continues += 1
            counts[2] += 1
            continue

        res_contig = np.int64(sub_info)

        if seq_id_is_none[anno_index] == 1:
            out_no_seq_id[n_no_seq_id] = anno_index
            n_no_seq_id += 1
            count_continues += 1
            counts[3] += 1
            continue

        # Set before the store check, exactly where the python loop set it, so a
        # structure with no residue store still registers its quality measures.
        out_quality_measure_seen[anno_index] = 1

        store_index = store_index_column[anno_index]
        if store_index == -1:
            count_continues += 1
            counts[4] += 1
            continue

        res_left = res_left_column[anno_index]
        row_index = res_left + (res_contig - 1)
        # res_right is one past the structure's last row (build_residue_batch counts it that way),
        # so row_index == res_right is already outside the structure. With '>' that one slipped
        # through: numba does not check bounds, so it read the next structure's first residue -
        # or, for a batch's last structure, past the end of the matrix, which the gather after
        # this kernel then turned into an IndexError that ended a whole-PDB run.
        if row_index >= res_right_column[anno_index]:
            out_bad_contig[n_bad_contig] = anno_index
            n_bad_contig += 1
            continue

        encoded_row = res_stores[store_index][row_index]

        res_aa = encoded_row[aa_column]
        pLDDT = encoded_row[plddt_column]

        location_code = locate_code(encoded_row[sc_rsa_column],
                                    surface_threshold, buried_threshold)
        simple_class_code = _simple_class_code(
            encoded_row[profile_left:profile_right], class_score_indices, location_code)

        seq_id = seq_id_column[anno_index]
        cov = cov_column[anno_index]
        resolution = resolution_column[anno_index]

        if cov_is_none[anno_index] == 1:
            # qualityScore raised TypeError on `coverage > 1.0`.
            qual = 0.01
            resolution = 3.0
        else:
            scaled_seq_id = seq_id
            if scaled_seq_id > 1.0:
                scaled_seq_id = scaled_seq_id / 100.0
            scaled_cov = cov
            if scaled_cov > 1.0:
                scaled_cov = scaled_cov / 100.0

            # qualityScore is written as `(1 + exp(...))**(-1)`, and libm's pow
            # disagrees with a plain division in the last bit for about one input in
            # 1700. The exponent arrives as an argument rather than a literal so LLVM
            # cannot strength-reduce the pow back into that division; np.exp itself
            # matches math.exp bit for bit.
            resolution_value = (1.0 + np.exp((1.5 * resolution) - 4.0)) ** inverse_power
            seq_value = (1.0 + np.exp(10.0 * (0.4 - scaled_seq_id))) ** inverse_power

            # qualityScore adds the three weighted terms with sum(), and since 3.12
            # CPython's sum() runs Neumaier compensated summation over floats rather
            # than a plain left-to-right add - it disagrees with `a + b + c` on about
            # one input in seven. Six extra flops on three terms, and the scores come
            # out bit for bit.
            total = 0.0
            compensation = 0.0
            for term in (resolution_value * 0.25, scaled_cov * 0.5, seq_value * 1.0):
                running = total + term
                if abs(total) >= abs(term):
                    compensation += (total - running) + term
                else:
                    compensation += (term - running) + total
                total = running
            weighted = (total + compensation) / 1.75

            if pLDDT != pLDDT:
                # np.isnan(pLDDT) -> pLDDT = 1.0, a python float, so the product
                # stayed float64.
                qual = weighted
            else:
                qual = np.float64(np.float32(pLDDT) * np.float32(weighted))

        out_anno[n_mappings] = anno_index
        out_res_contig[n_mappings] = res_contig
        out_row[n_mappings] = row_index
        out_store[n_mappings] = store_index
        out_qual[n_mappings] = qual
        out_seq_id[n_mappings] = seq_id
        out_cov[n_mappings] = cov
        out_aa_id[n_mappings] = has_aa1_code and res_aa == aa1_code
        out_resolution[n_mappings] = resolution
        out_res_aa[n_mappings] = res_aa
        out_simple_class[n_mappings] = simple_class_code
        n_mappings += 1

    return n_mappings, count_continues, n_no_seq_id, n_bad_contig


def record_bad_contig(bad_contigs: dict, struct_db_id: int, structure: str, pos: int,
                      res_contig: int, n_rows: int) -> None:
    """Count one mapping whose residue lies past the rows stored for its structure.

    Kept per structure for the whole protein and reported once (summarize_bad_contigs). Written
    as one warning per position and structure, with a stack trace each, it grew with sequence
    length - and the off-by-one fix in _collect_mappings_for_position sends more of them here.
    """
    entry = bad_contigs.get(struct_db_id)
    if entry is None:
        bad_contigs[struct_db_id] = [1, structure, pos, res_contig, n_rows]
    else:
        entry[0] += 1


def summarize_bad_contigs(protein_id: str, bad_contigs: dict, limit: int = 5) -> str:
    """One warning text for all of a protein's bad contigs, the worst structures first."""
    total = sum(entry[0] for entry in bad_contigs.values())
    worst = sorted(bad_contigs.items(), key=lambda item: -item[1][0])[:limit]
    examples = '; '.join(
        f'{structure} (db id {struct_db_id}): {count}x, e.g. position {pos} -> residue {res_contig} '
        f'of {n_rows} stored'
        for struct_db_id, (count, structure, pos, res_contig, n_rows) in worst)
    more = f' and {len(bad_contigs) - limit} more' if len(bad_contigs) > limit else ''
    return (f'res contig value does match left,right: {protein_id} has {total} mapping(s) to a '
            f'residue past those stored for the structure, in {len(bad_contigs)} structure(s); '
            f'skipped. {examples}{more}')


def _gather_structure_columns(
        annotation_list: list[int],
        structs_from_db: np.ndarray,
        structure_db_list: np.ndarray,
        concat_res_data_stores: list[np.ndarray],
        resolutions_stores: list[np.ndarray]):
    """The per-structure columns the mapping loop reads, as parallel lists.

    structs_from_db is indexed once as a block rather than row by row, because
    unpacking a row per mapping boxes seven numpy scalars. The residue data store
    and the resolution are resolved here too, and the fused id is decoded into its
    structure id and chain once per structure rather than once per mapping.
    Structures with no residue store (store_index == -1) are still listed, with
    None in the columns that would need one; the loop skips them on store_index.
    """
    if len(annotation_list) == 0:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty, empty, np.empty(0, dtype=np.float64), [], []

    # structs_from_db is the protein's LocalStructureRows, so rows are looked up per
    # id rather than fancy-indexed. int32, the dtype the kernel was compiled against.
    rows = np.array([structs_from_db[struct_db_id] for struct_db_id in annotation_list],
                    dtype=np.int32).reshape(-1, 7)
    store_indices = np.ascontiguousarray(rows[:, 0])
    res_lefts = np.ascontiguousarray(rows[:, 1])
    res_rights = np.ascontiguousarray(rows[:, 2])

    resolutions = np.zeros(len(annotation_list), dtype=np.float64)
    for index in range(len(annotation_list)):
        store_index = store_indices[index]
        if store_index == -1:
            continue
        resolutions[index] = resolutions_stores[store_index][rows[index, 6]]

    structure_ids = []
    chains = []
    for struct_db_id in annotation_list:
        fused_id = structure_db_list[struct_db_id]
        structure_ids.append(fused_id[:-1].decode())
        chains.append(fused_id[-1:].decode())

    return store_indices, res_lefts, res_rights, resolutions, structure_ids, chains


def _structure_rank_column(structure_ids: list[str]) -> np.ndarray:
    """Lexicographic rank of each structure id among this protein's structures.

    set_recommended_residues breaks quality ties on `mapping_id[0] > ...`, a
    comparison of the structure id strings. Ranking them once per protein turns
    that into an integer comparison the kernel can make. Ranked with python's own
    ordering, and equal ids share a rank, so `>` behaves exactly as it did.
    """
    ranks = {name: rank for rank, name in enumerate(sorted(set(structure_ids)))}
    return np.array([ranks[name] for name in structure_ids], dtype=np.int64)


def _residue_store_list(concat_res_data_stores: list[np.ndarray]) -> NumbaList:
    """The residue data stores as a numba list, so a kernel can index them by store_index.

    A reflected python list of arrays is not addressable from nopython code; a typed
    list is, provided every element shares dtype and layout. An empty input still has
    to produce a typed list, because a protein whose structures all have
    store_index == -1 never indexes it but numba still needs the type.

    Number of dimensions is part of that type. A chain whose analysis list was empty
    is stored 1-D by insertion_lib.remote_insertResidues, and numba types the whole
    list from the first element it is given, so a (0,) array next to real
    (n_residues, n_features) stores fails to compile. build_local_stores already
    hands over a 2-D matrix; the guard stays so any other caller gets the same
    promise. Widening an empty store costs nothing: it is only ever typed, never
    indexed.

    The placeholder is ENCODED_ROW_LEN columns wide rather than one, because the
    mapping loop sizes an empty result matrix from residue_stores[0].shape[1].
    """
    empty_store = np.zeros((0, residue_package.ENCODED_ROW_LEN), dtype=np.float32)

    stores = NumbaList()
    for store in concat_res_data_stores:
        store = np.ascontiguousarray(store, dtype=np.float32)
        if store.ndim != 2 and store.size == 0:
            store = empty_store
        stores.append(store)
    if len(stores) == 0:
        stores.append(empty_store)
    return stores

# The row of the per-structure table for a structure that has no retrieved data. Same
# 7 columns as getStoredResidues' total_structures: store_index, res_left, res_right,
# index_left, index_right, compl_db_id, resolution_index.
_ABSENT_STRUCTURE_ROW = (-1, 0, 0, 0, 0, -1, -1)


class LocalStructureRows(dict):
    """struct_db_id -> row, shaped like getStoredResidues' offset table, for one protein.

    Every structure the protein needs is decoded into local stores, so every row it
    holds has store_index 0 and offsets into those stores. An id that was not decoded -
    because it has no stored residue data - answers with store_index -1, exactly what
    the global table said for it, so the kernel and the interface map skip it without
    a special case. Readers index it by id and unpack the row or take [5]; nothing
    fancy-indexes it.
    """
    __slots__ = ()

    def __missing__(self, key):
        return _ABSENT_STRUCTURE_ROW


class DecodedStructureCache:
    """Bounded LRU of decoded structures, shared by the proteins of one task.

    Proteins of one package often annotate the same templates, and decoding a residue
    blob is the per-protein cost the compressed store introduced. A budget in bytes
    rather than entries, because a structure's matrix ranges from a few rows to
    thousands. Kept small by default: every classification task has one.
    """
    __slots__ = ['budget', 'used', 'entries', 'hits', 'misses']

    def __init__(self, budget_bytes: int):
        self.budget = budget_bytes
        self.used = 0
        self.entries: OrderedDict = OrderedDict()
        self.hits = 0
        self.misses = 0

    def get(self, key):
        entry = self.entries.get(key)
        if entry is None:
            self.misses += 1
            return None
        self.hits += 1
        self.entries.move_to_end(key)
        return entry[0]

    def put(self, key, value, nbytes: int) -> None:
        if nbytes > self.budget:
            return
        if key in self.entries:
            self.used -= self.entries.pop(key)[1]
        self.entries[key] = (value, nbytes)
        self.used += nbytes
        while self.used > self.budget:
            _, (_, evicted_bytes) = self.entries.popitem(last=False)
            self.used -= evicted_bytes


def _blob(buffer: np.ndarray, offsets: np.ndarray, blob_index: int) -> bytes | None:
    left = int(offsets[blob_index])
    right = int(offsets[blob_index + 1])
    if left == right:
        return None
    # A copy of one compressed blob, tens of KB: unpack sniffs the zstd magic with a
    # bytes comparison, which a uint8 view would turn into an elementwise one.
    return buffer[left:right].tobytes()


def decode_residue_blob(blob: bytes | None) -> np.ndarray | None:
    """The encoded residue matrix of one structure, as the mapping kernel needs it.

    Mirrors what process_results_chunk and _residue_store_list did to a stored blob:
    rows written before the float32 switch are narrowed, and a chain with no residues
    - stored 1-D and empty - becomes an empty 2-D matrix. None means no residue data.
    """
    if blob is None:
        return None
    residue_data = unpack(blob)
    if residue_data is None:
        return None
    if residue_data.ndim != 2:
        if residue_data.size == 0:
            return np.zeros((0, residue_package.ENCODED_ROW_LEN), dtype=np.float32)
    return np.ascontiguousarray(residue_data, dtype=np.float32)


def decode_interface_blob(blob: bytes | None) -> list:
    """(chain, combined_keys, scores) triples of one structure; empty when none are stored."""
    if blob is None:
        return []
    interfaces = unpack(blob)
    if interfaces is None:
        return []
    return interfaces


def interface_structure_ids(annotation_list, compl_db_id_map, compl_struct_db_id_matrix) -> list[int]:
    """Every chain of every complex the protein's structures belong to.

    Exactly the ids calculate_aggregated_interface_map resolves: it walks
    compl_struct_db_id_matrix[compl_db_id_map[structure]] for each annotated structure.
    """
    seen_complexes: set[int] = set()
    chain_ids: list[int] = []
    seen_ids: set[int] = set()
    for struct_db_id in annotation_list:
        compl_db_id = int(compl_db_id_map[struct_db_id])
        if compl_db_id == -1 or compl_db_id in seen_complexes:
            continue
        seen_complexes.add(compl_db_id)
        for chain_db_id in compl_struct_db_id_matrix[compl_db_id]:
            chain_db_id = int(chain_db_id)
            if chain_db_id == -1 or chain_db_id in seen_ids:
                continue
            seen_ids.add(chain_db_id)
            chain_ids.append(chain_db_id)
    return chain_ids


def interface_annotation_list(
        annotation_list: list[int],
        seq_id_column: np.ndarray,
        compl_db_id_map: np.ndarray,
        compl_struct_db_id_matrix: np.ndarray,
        config: Config) -> tuple[list[int], int]:
    """The protein's annotations that go into the interface calculation, and how many
    of all its annotations are in a complex with more than one protein chain.

    The selection database.prune_interface_backmaps made when the backmaps were
    written, made again on the same inputs: the stored sequence identities (NaN where
    None, as seq_id_column holds them) and the protein chains of each complex, as
    compl_struct_db_id_matrix holds them - Proteins.protein_chain_count is what filled
    it. A structure without a complex counts as single chain;
    calculate_aggregated_interface_map skips those anyway.
    """
    if len(annotation_list) == 0:
        return [], 0
    struct_db_ids = np.asarray(annotation_list, dtype=np.int64)
    compl_db_ids = np.asarray(compl_db_id_map[struct_db_ids], dtype=np.int64)
    has_complex = compl_db_ids >= 0
    multi_chain = np.zeros(len(annotation_list), dtype=bool)
    if has_complex.any():
        chain_counts = (compl_struct_db_id_matrix[compl_db_ids[has_complex]] != -1).sum(axis=1)
        multi_chain[has_complex] = chain_counts > 1
    selected = interface_package.select_interface_annotations(
        struct_db_ids, seq_id_column, multi_chain,
        *interface_package.interface_annotation_limits(config))
    return ([annotation_list[index] for index in np.flatnonzero(selected).tolist()],
            int(multi_chain.sum()))


# The fill buffer starts here and is then sized from the data itself (see
# build_local_stores). Small enough that a wrong projection can never cost much, big
# enough that a protein with a handful of structures never resizes at all.
_INITIAL_RESIDUE_ROWS: int = 64 * 1024

# Non-empty structures decoded before the total is projected. One structure is a noisy
# sample - an unusually sparse chain would project far too small - and sixteen out of the
# tens of thousands a heavy protein carries cost nothing to buffer.
_CALIBRATION_STRUCTURES: int = 16

# How often the projection fell short and the buffer had to be enlarged, and how many
# rows were allocated against how many were used. Reported on the worker done line: many
# grows or a large slack ratio means the projection is not tracking the data.
_residue_buffer_stats = {'grows': 0, 'rows_allocated': 0, 'rows_used': 0}


def residue_buffer_report() -> str:
    stats = _residue_buffer_stats
    if stats['rows_used'] == 0:
        return 'res_buf=n/a'
    return (f'res_buf grows={stats["grows"]} '
            f'slack={stats["rows_allocated"] / stats["rows_used"]:.2f}x')


def residue_compressed_bytes(residue_ids, structs_from_db: np.ndarray,
                             blob_stores: list[tuple[np.ndarray, ...]]) -> int:
    """Total compressed size of the residue blobs of these structures.

    The blobs stay compressed in the object store, so nothing knows a structure's row
    count without decoding it - but the compressed span of each is already in the store's
    offset array, and rows-per-compressed-byte is near constant within a run. Measuring
    that ratio on the first structures decoded turns this sum into a projection of the
    total, which is what the fill buffer is sized from.
    """
    total = 0
    for struct_db_id in residue_ids:
        global_row = structs_from_db[int(struct_db_id)]
        store_index = int(global_row[0])
        if store_index == -1:
            continue
        blob_index = int(global_row[1])
        res_offsets = blob_stores[store_index][1]
        total += int(res_offsets[blob_index + 1]) - int(res_offsets[blob_index])
    return total


def _grow_residue_matrix(matrix: np.ndarray, needed_rows: int) -> np.ndarray:
    """Enlarge the fill buffer in place, doubling so growth is amortised.

    resize rather than a new array plus a copy: it reallocates where it can and never
    holds the old and the new at python level. C-contiguous, so enlarging the first axis
    appends rows and leaves the existing ones where they are.
    """
    matrix.resize((max(2 * matrix.shape[0], needed_rows), matrix.shape[1]), refcheck=False)
    _residue_buffer_stats['grows'] += 1
    return matrix


def build_local_stores(
        residue_ids,
        interface_ids,
        structs_from_db: np.ndarray,
        blob_stores: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
        cache: DecodedStructureCache):
    """Decode the structures one protein needs into stores shaped like the global ones.

    Returns (rows, residue_stores, resolutions_stores, index_list_stores,
    chain_list_stores, combined_keys_stores, scores_stores): each *_stores a one-element
    list, so every consumer keeps addressing data as store[store_index][offset] with
    store_index 0. The residue matrix is one contiguous float32 array, which also means
    the per-position gather never needs its mixed-store fallback.

    residue_ids get their residue rows decoded, interface_ids their interface entries;
    an id in both gets both, in one row.

    para_classify no longer calls this: it builds the interface half once per protein and
    the residue half per batch, so that a task never holds the whole matrix
    (build_interface_stores / build_residue_batch). Kept as the single-shot form both
    halves are checked against - dev_tests/test_compressed_residue_stores.py and
    test_local_store_buffer.py decode through it - and as the reference for what those two
    together have to produce.
    """
    rows = LocalStructureRows()

    residue_set = set(int(i) for i in residue_ids)
    interface_set = set(int(i) for i in interface_ids)

    # One buffer filled in place, not a list of decoded matrices plus the
    # np.concatenate of it. Both used to be alive at the concatenation, and that was the
    # task's whole peak: a p50 whole-PDB protein's matrix is 5.1 GiB and the measured
    # per-task footprint was 10.6 GiB, which is what held the stage to ~115 of 255 cores.
    #
    # It starts small and is sized once from the data: after the first
    # _CALIBRATION_STRUCTURES non-empty structures, rows-per-compressed-byte is known and
    # the remaining blobs' compressed sizes project the total. Projecting from a constant
    # instead was accurate on the measured run (214 B a row against an assumed 210) but
    # over-allocated 3.7x on data that compresses differently, which is the very thing
    # this is here to avoid.
    residue_matrix = np.empty((_INITIAL_RESIDUE_ROWS, residue_package.ENCODED_ROW_LEN),
                              dtype=np.float32)
    res_right = 0
    remaining_compressed = residue_compressed_bytes(residue_ids, structs_from_db, blob_stores)
    calibration_rows = 0
    calibration_bytes = 0
    calibration_structures = 0
    projected = False
    key_parts: list[np.ndarray] = []
    score_parts: list[np.ndarray] = []
    index_list: list[tuple[int, int]] = []
    chain_list: list[str] = []
    entry_right = 0
    resolutions: list[np.float32] = []

    ordered_ids = list(dict.fromkeys([int(i) for i in residue_ids] + [int(i) for i in interface_ids]))

    for struct_db_id in ordered_ids:
        global_row = structs_from_db[struct_db_id]
        store_index = int(global_row[0])
        if store_index == -1:
            continue
        blob_index = int(global_row[1])
        res_buffer, res_offsets, iface_buffer, iface_offsets, store_resolutions = blob_stores[store_index]

        res_left = res_right
        if struct_db_id in residue_set:
            residue_data = cache.get(('r', struct_db_id))
            if residue_data is None:
                residue_data = decode_residue_blob(_blob(res_buffer, res_offsets, blob_index))
                if residue_data is None:
                    continue
                cache.put(('r', struct_db_id), residue_data, residue_data.nbytes)
            blob_bytes = int(res_offsets[blob_index + 1]) - int(res_offsets[blob_index])
            remaining_compressed -= blob_bytes

            n_rows = len(residue_data)
            if n_rows > 0:
                if res_right + n_rows > residue_matrix.shape[0]:
                    residue_matrix = _grow_residue_matrix(residue_matrix, res_right + n_rows)
                residue_matrix[res_right:res_right + n_rows] = residue_data
                res_right += n_rows

                if not projected:
                    calibration_rows += n_rows
                    calibration_bytes += blob_bytes
                    calibration_structures += 1
                    if calibration_structures >= _CALIBRATION_STRUCTURES and calibration_bytes > 0:
                        # Rows still to come, from the ratio these structures showed. The
                        # margin covers the sample being a little sparser than the rest;
                        # being short only costs a doubling, and the shrink at the end
                        # gives back whatever is unused either way.
                        rows_per_byte = calibration_rows / calibration_bytes
                        target = res_right + int(1.05 * max(0, remaining_compressed) * rows_per_byte) + 1
                        if target > residue_matrix.shape[0]:
                            residue_matrix.resize((target, residue_matrix.shape[1]), refcheck=False)
                        projected = True
            del residue_data

        index_left = len(index_list)
        if struct_db_id in interface_set:
            interfaces = cache.get(('i', struct_db_id))
            if interfaces is None:
                interfaces = decode_interface_blob(_blob(iface_buffer, iface_offsets, blob_index))
                cache.put(('i', struct_db_id), interfaces,
                          sum(keys.nbytes + scores.nbytes for _, keys, scores in interfaces))
            for ichain, combined_keys, scores in interfaces:
                key_parts.append(combined_keys)
                score_parts.append(scores)
                index_list.append((entry_right, entry_right + len(combined_keys)))
                entry_right += len(combined_keys)
                chain_list.append(ichain)
        index_right = len(index_list)

        resolution_index = len(resolutions)
        resolutions.append(np.float32(store_resolutions[int(global_row[6])]))

        rows[struct_db_id] = (0, res_left, res_right, index_left, index_right,
                              int(global_row[5]), resolution_index)

    # Give back whatever the estimate over-allocated, in place, so the slack does not
    # travel on to the kernel with the matrix.
    _residue_buffer_stats['rows_allocated'] += residue_matrix.shape[0]
    _residue_buffer_stats['rows_used'] += res_right
    # Give back whatever the projection over-allocated, in place, so the slack does not
    # travel on to the kernel with the matrix.
    if res_right != residue_matrix.shape[0]:
        residue_matrix.resize((res_right, residue_package.ENCODED_ROW_LEN), refcheck=False)

    if len(key_parts) == 0:
        combined_keys = np.empty(0, dtype=np.int64)
        scores = np.empty(0, dtype=np.float32)
    else:
        combined_keys = np.concatenate(key_parts)
        scores = np.concatenate(score_parts)

    return (
        rows,
        _residue_store_list([residue_matrix]),
        [np.array(resolutions, dtype=np.float32)],
        # Same dtypes as the global stores getStoredResidues used to build.
        [np.array(index_list, dtype=np.int32).reshape(-1, 2)],
        [np.array(chain_list, dtype='S')],
        [combined_keys],
        [scores],
    )


def classification_memory_estimate(n_annotations: int, n_interface_chains: int,
                                   config: Config) -> int:
    """Bytes a para_classify task needs for a package with these counts.

    What a task actually holds is what build_local_stores decodes: one residue matrix per
    annotated structure, and the interface entries of every chain of those structures'
    complexes. Neither follows the package cost. cost is positions x annotations, so a
    long sequence annotated onto small structures scores high and decodes little, while a
    chain of a large assembly scores low and decodes every chain of every complex it
    appears in - 9CI3:A (cost 36.5M) peaked at 6.8 GB where 7JWB:D (cost 12.5M) peaked at
    49 GB. Measured peak per 1M cost spanned 186 MB to 4748 MB across the tasks that
    finished, which is why no single per-cost coefficient can be both safe and useful.

    Sizing on the two counts separates them. The coefficients are first-principles
    starting points, not a fit: no run has yet logged interface chain counts next to a
    memory peak. The done line logs both now, so they can be fitted from the next one.
    """
    return int(classification_fixed_task_bytes(config)
               + classification_protein_bytes(n_annotations, n_interface_chains, config))


def classification_fixed_task_bytes(config: Config) -> int:
    """What a worker costs before it decodes anything: interpreter, imports, caches.

    Per task, not per protein - a package of five proteins starts one worker, not five.
    Charging it per protein made a package of light proteins look five times its size and
    let the budget fill after 36 of 255 slots.
    """
    return int(getattr(config, 'classification_min_task_gib', 1.0) * 1024**3)


def classification_protein_bytes(n_annotations: int, n_interface_chains: int,
                                 config: Config) -> int:
    """The part of a task's footprint that one protein adds, so packages can sum it."""
    per_annotation = getattr(config, 'classification_mb_per_annotation', 0.0) * 1024**2
    per_chain = getattr(config, 'classification_mb_per_interface_chain', 0.0) * 1024**2
    return int(per_annotation * max(0, n_annotations) + per_chain * max(0, n_interface_chains))


def build_interface_stores(
        interface_ids,
        structs_from_db: np.ndarray,
        blob_stores: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
        cache: DecodedStructureCache):
    """The interface half of build_local_stores, for a whole protein at once.

    Batching applies to the residue data and not to this: the interface keys of a heavy
    whole-PDB protein measure ~22 MiB (the [interface][loop] keys= counters), and
    calculate_aggregated_interface_map walks every complex of the protein in one pass, so
    splitting it would cost a restructure and save nothing.

    Returns (rows, index_list_stores, chain_list_stores, combined_keys_stores,
    scores_stores) with the same shapes build_local_stores produced, so the interface map
    reads them unchanged. Only columns 0, 3, 4 and 5 of a row are meaningful here - the
    residue columns belong to the per-batch store.
    """
    rows = LocalStructureRows()
    key_parts: list[np.ndarray] = []
    score_parts: list[np.ndarray] = []
    index_list: list[tuple[int, int]] = []
    chain_list: list[str] = []
    entry_right = 0

    for struct_db_id in dict.fromkeys(int(i) for i in interface_ids):
        global_row = structs_from_db[struct_db_id]
        store_index = int(global_row[0])
        if store_index == -1:
            continue
        blob_index = int(global_row[1])
        _res_buffer, _res_offsets, iface_buffer, iface_offsets, _resolutions = blob_stores[store_index]

        interfaces = cache.get(('i', struct_db_id))
        if interfaces is None:
            interfaces = decode_interface_blob(_blob(iface_buffer, iface_offsets, blob_index))
            cache.put(('i', struct_db_id), interfaces,
                      sum(keys.nbytes + scores.nbytes for _, keys, scores in interfaces))

        index_left = len(index_list)
        for ichain, combined_keys, scores in interfaces:
            key_parts.append(combined_keys)
            score_parts.append(scores)
            index_list.append((entry_right, entry_right + len(combined_keys)))
            entry_right += len(combined_keys)
            chain_list.append(ichain)
        index_right = len(index_list)

        rows[struct_db_id] = (0, 0, 0, index_left, index_right, int(global_row[5]), 0)

    if len(key_parts) == 0:
        combined_keys = np.empty(0, dtype=np.int64)
        scores = np.empty(0, dtype=np.float32)
    else:
        combined_keys = np.concatenate(key_parts)
        scores = np.concatenate(score_parts)

    return (
        rows,
        [np.array(index_list, dtype=np.int32).reshape(-1, 2)],
        [np.array(chain_list, dtype='S')],
        [combined_keys],
        [scores],
    )


def build_residue_batch(
        batch_ids,
        structs_from_db: np.ndarray,
        blob_stores: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
        cache: DecodedStructureCache):
    """The residue half of build_local_stores, for one batch of a protein's annotations.

    The same self-sizing fill buffer build_local_stores uses - start small, project the
    total from the compressed sizes once the first structures have been decoded, grow if
    short, shrink at the end - but over a batch rather than the whole annotation list.
    That is the point of batching: a heavy whole-PDB protein's full matrix reaches 8.6 GiB,
    a 1000-structure batch of it about 150 MB.

    Returns (rows, residue_stores, resolutions_stores). Only columns 0, 1, 2 and 6 of a
    row are meaningful; the interface columns belong to the per-protein store.
    """
    rows = LocalStructureRows()
    residue_matrix = np.empty((_INITIAL_RESIDUE_ROWS, residue_package.ENCODED_ROW_LEN),
                              dtype=np.float32)
    res_right = 0
    resolutions: list[np.float32] = []

    remaining_compressed = residue_compressed_bytes(batch_ids, structs_from_db, blob_stores)
    calibration_rows = 0
    calibration_bytes = 0
    calibration_structures = 0
    projected = False

    for struct_db_id in dict.fromkeys(int(i) for i in batch_ids):
        global_row = structs_from_db[struct_db_id]
        store_index = int(global_row[0])
        if store_index == -1:
            continue
        blob_index = int(global_row[1])
        res_buffer, res_offsets, _iface_buffer, _iface_offsets, store_resolutions = blob_stores[store_index]

        residue_data = cache.get(('r', struct_db_id))
        if residue_data is None:
            residue_data = decode_residue_blob(_blob(res_buffer, res_offsets, blob_index))
            if residue_data is None:
                continue
            cache.put(('r', struct_db_id), residue_data, residue_data.nbytes)

        blob_bytes = int(res_offsets[blob_index + 1]) - int(res_offsets[blob_index])
        remaining_compressed -= blob_bytes

        res_left = res_right
        n_rows = len(residue_data)
        if n_rows > 0:
            if res_right + n_rows > residue_matrix.shape[0]:
                residue_matrix = _grow_residue_matrix(residue_matrix, res_right + n_rows)
            residue_matrix[res_right:res_right + n_rows] = residue_data
            res_right += n_rows

            if not projected:
                calibration_rows += n_rows
                calibration_bytes += blob_bytes
                calibration_structures += 1
                if calibration_structures >= _CALIBRATION_STRUCTURES and calibration_bytes > 0:
                    rows_per_byte = calibration_rows / calibration_bytes
                    target = res_right + int(1.05 * max(0, remaining_compressed) * rows_per_byte) + 1
                    if target > residue_matrix.shape[0]:
                        residue_matrix.resize((target, residue_matrix.shape[1]), refcheck=False)
                    projected = True
        del residue_data

        resolution_index = len(resolutions)
        resolutions.append(np.float32(store_resolutions[int(global_row[6])]))

        rows[struct_db_id] = (0, res_left, res_right, 0, 0, int(global_row[5]), resolution_index)

    _residue_buffer_stats['rows_allocated'] += residue_matrix.shape[0]
    _residue_buffer_stats['rows_used'] += res_right
    if res_right != residue_matrix.shape[0]:
        residue_matrix.resize((res_right, residue_package.ENCODED_ROW_LEN), refcheck=False)

    return (rows,
            _residue_store_list([residue_matrix]),
            [np.array(resolutions, dtype=np.float32)])


def _gather_identity_columns(annotation_list: list[int], structure_db_list: np.ndarray):
    """(structure id, chain) per annotation, for the whole protein.

    Kept whole even when the residue data is batched: MappingIds indexes into these by a
    protein-wide annotation index, and they are strings shared by reference, so a heavy
    protein's pair of lists costs a few MB rather than anything that wants batching.
    """
    structure_ids = []
    chains = []
    for struct_db_id in annotation_list:
        fused_id = structure_db_list[struct_db_id]
        structure_ids.append(fused_id[:-1].decode())
        chains.append(fused_id[-1:].decode())
    return structure_ids, chains


def _gather_residue_columns(batch_ids, local_rows, resolutions_stores):
    """The per-structure columns the mapping kernel reads, for one batch.

    Same columns _gather_structure_columns produces, minus the identity pair that is now
    gathered once per protein.
    """
    if len(batch_ids) == 0:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty, empty, np.empty(0, dtype=np.float64)

    rows = np.array([local_rows[struct_db_id] for struct_db_id in batch_ids],
                    dtype=np.int32).reshape(-1, 7)
    store_indices = np.ascontiguousarray(rows[:, 0])
    res_lefts = np.ascontiguousarray(rows[:, 1])
    res_rights = np.ascontiguousarray(rows[:, 2])

    resolutions = np.zeros(len(batch_ids), dtype=np.float64)
    for index in range(len(batch_ids)):
        if store_indices[index] == -1:
            continue
        resolutions[index] = resolutions_stores[store_indices[index]][rows[index, 6]]

    return store_indices, res_lefts, res_rights, resolutions


def classification_memory_request(package_cost: int, config: Config,
                                  estimate_bytes: int | None = None) -> int | None:
    """Bytes a para_classify task reserves from ray's memory resource, or None for none.

    A task's footprint is set by the proteins it carries: build_local_stores decodes the
    residue matrix of every structure a protein maps to and the interface entries of
    every chain of their complexes, so a ribosomal protein mapped to thousands of
    ribosomes costs gigabytes on its own. Every chunk of a whole-PDB input holds such
    proteins - sequence deduplication only fuses within a chunk - and pack_packages
    hands out the heaviest first, one per core. A verbosity 4 run measured the result:
    available memory fell from 827 GiB to 102 GiB within a minute of 255 tasks
    starting, and the stage then crawled at the cgroup limit for hours. Smaller chunks
    changed nothing, because the per-protein cost does not depend on the chunk.

    Reserving by cost makes ray itself admit only as many tasks as fit. The request is
    capped at the whole pool, since a task asking for more would never be scheduled.
    """
    gib_per_mcost = getattr(config, 'classification_gib_per_mcost', 0.0)
    pool = getattr(config, 'ray_task_memory_bytes', None)
    if not pool:
        return None
    if estimate_bytes is not None:
        # The counts were available, so the estimate is the one built from what the task
        # decodes rather than from its cost.
        request = float(estimate_bytes)
    elif gib_per_mcost > 0:
        request = gib_per_mcost * (package_cost / 1_000_000) * 1024**3
    else:
        return None
    request = max(request, getattr(config, 'classification_min_task_gib', 1.0) * 1024**3)
    return int(min(request, pool))


def package_prot_db_id_map(prot_db_id_map_array: np.ndarray,
                           package: list[tuple[str, int, int]]) -> dict[int, tuple[int, int]]:
    """prot_db_id -> (left, right) into pos_db_id_vec, for this package's proteins only.

    The array itself is shared zero-copy through the object store, but a dict built
    from all of it is a private copy in every task - one entry per protein in the whole
    chunk, ~150k of them with boxed numpy scalars on a whole-PDB input, tens of MB per
    worker - while the only lookup it serves asks for ids that are in this package.
    Plain ints rather than numpy scalars for the same reason the annotation columns are.
    """
    wanted = {entry[1] for entry in package}
    return {
        int(row[0]): (int(row[1]), int(row[2]))
        for row in prot_db_id_map_array if int(row[0]) in wanted
    }


def classification_cost_cap(config: Config) -> int | None:
    """Largest package cost that still lets every core run at once, or None for no cap.

    The mirror image of classification_memory_request: that function turns a package's
    cost into the memory ray admits it against, this one turns the memory a task may
    have into the cost a package may carry. Both read classification_gib_per_mcost, so
    packaging cannot outrun admission.

    pack_packages used to size packages purely by cost - two per core over the whole
    stage, capped at 100M - with no idea what that costs in bytes. On a whole-PDB run
    that produced packages of 44M cost, which reserve 44 GiB each, and only 15 of those
    fit the 646 GiB pool: 240 of 255 cores waited while the packaging happily kept
    handing out more of them.
    """
    gib_per_mcost = getattr(config, 'classification_gib_per_mcost', 0.0)
    pool = getattr(config, 'ray_task_memory_bytes', None)
    # Both of these already mean "no reservation is made", so there is nothing to
    # stay below either.
    if gib_per_mcost <= 0 or not pool:
        return None
    procs = max(1, config.proc_n)
    budget = pool / procs
    return max(1, int(budget / (gib_per_mcost * 1024**3 / 1_000_000)))


# Which stage of a classification task is running, for PrivateMemoryPeak to attribute
# its maximum to. A one-element list rather than a global so the setter needs no `global`
# statement and the sampler thread reads it with a single bytecode - it is touched once
# per protein and once per batch, never per annotation. Process-wide because a worker
# runs one task at a time (max_calls=1).
_TASK_PHASE: list[str] = ['start']


def set_task_phase(name: str) -> None:
    _TASK_PHASE[0] = name


# The decode cache lives inside para_classify, so the wrapper that writes the done line
# cannot reach it. Its final counters are parked here instead, the same way
# mm_lookup_cache_stats reports the complex lookup cache.
_decode_cache_stats: dict[str, int] = {'used': 0, 'budget': 0, 'hits': 0, 'misses': 0}


def record_decode_cache(cache: 'DecodedStructureCache') -> None:
    _decode_cache_stats['used'] = cache.used
    _decode_cache_stats['budget'] = cache.budget
    _decode_cache_stats['hits'] = cache.hits
    _decode_cache_stats['misses'] = cache.misses


def decode_cache_report() -> str:
    """How much the decode cache held and how often it answered.

    A low hit rate means the proteins of a package annotate different templates and the
    cache is only costing memory; a high one means batching would otherwise be re-decoding
    the same blob once per batch.
    """
    stats = _decode_cache_stats
    lookups = stats['hits'] + stats['misses']
    rate = f'{100.0 * stats["hits"] / lookups:.0f}%' if lookups else 'n/a'
    return (f'cache_used={stats["used"] / 1024**2:.0f}/{stats["budget"] / 1024**2:.0f} MB '
            f'cache_hits={stats["hits"]} cache_misses={stats["misses"]} ({rate})')


class PrivateMemoryPeak:
    """Largest private (USS) and proportional (PSS) footprint a task reached.

    smaps_rollup is a point-in-time read and the kernel keeps no peak-USS counter, so
    the only way to a peak is to sample. It matters because the one-shot reading taken
    at the end of a task reported 4.6 GiB where the task had actually peaked at 10.6:
    build_local_stores holds the decoded residue matrix twice while concatenating it,
    and by the time the task finishes one of the copies is gone. Calibrating
    classification_gib_per_mcost from the end reading would set it 2x too low.

    A daemon thread, so it can never keep a worker alive, and every read is guarded -
    monitoring must not be the thing that breaks a task.
    """

    __slots__ = ('interval', 'peak_uss', 'peak_pss', 'samples', 'peak_phase',
                 'phase_peaks', '_stop', '_thread')

    def __init__(self, interval: float = 2.0):
        self.interval = interval
        self.peak_uss = 0
        self.peak_pss = 0
        self.samples = 0
        # Which stage held the maximum, and the maximum each stage reached. Knowing that
        # the peak is 18.7 GiB does not say what to shrink; this does.
        self.peak_phase = 'none'
        self.phase_peaks: dict[str, int] = {}
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def sample(self) -> None:
        breakdown = process_memory_breakdown()
        if breakdown is None:
            return
        self.samples += 1
        phase = _TASK_PHASE[0]
        uss = breakdown['uss']
        if uss > self.phase_peaks.get(phase, 0):
            self.phase_peaks[phase] = uss
        if uss > self.peak_uss:
            self.peak_uss = uss
            self.peak_phase = phase
        if breakdown['pss'] > self.peak_pss:
            self.peak_pss = breakdown['pss']

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            try:
                self.sample()
            except Exception:
                pass

    def start(self) -> 'PrivateMemoryPeak':
        self.sample()  # so a task shorter than one interval still reports something
        self._thread = threading.Thread(target=self._run, name='structman-uss-peak', daemon=True)
        self._thread.start()
        return self

    def stop(self) -> None:
        self._stop.set()
        try:
            self.sample()
        except Exception:
            pass

    def report(self) -> str:
        if self.samples == 0:
            return 'peak_uss=n/a peak_pss=n/a'
        # The three heaviest stages, largest first: enough to say where the footprint
        # sits without a line per stage on every task.
        by_size = sorted(self.phase_peaks.items(), key=lambda entry: -entry[1])[:3]
        phases = ' '.join(f'{name}:{value / 1024**2:.0f}' for name, value in by_size)
        return (f'peak_uss={self.peak_uss / 1024**2:.0f} MB '
                f'peak_pss={self.peak_pss / 1024**2:.0f} MB '
                f'samples={self.samples} peak_phase={self.peak_phase} '
                f'phase_uss_mb=[{phases}]')


class PackageBounds(NamedTuple):
    """What the dispatcher is allowed to do, and which rule said so.

    Computed in one place because pack_packages and the wait loop that re-invokes it have
    to agree: if the loop thinks more tasks may be dispatched than pack_packages will
    accept, it calls it once per iteration for nothing, and if it thinks fewer, the stage
    runs below the window. The candidates are carried so the log line can say which bound
    is binding - that is how the memory bound was found to be the one that mattered.
    """

    structs_thresh: int
    per_task_bytes: int
    window: int
    by_stage: int
    ceiling: int
    by_memory: int | None


def classification_package_bounds(config: Config, total_cost: int) -> PackageBounds:
    by_stage = total_cost // (2 * max(1, config.proc_n))
    ceiling = 100_000_000
    by_memory = classification_cost_cap(config)
    structs_thresh = min([by_stage, ceiling] + ([by_memory] if by_memory is not None else []))
    per_task_bytes = classification_memory_request(structs_thresh, config) or 0
    window = classification_dispatch_window(config, per_task_bytes)
    return PackageBounds(structs_thresh, per_task_bytes, window, by_stage, ceiling, by_memory)


# Measured headroom, refreshed at most this often; see classification_dispatch_window.
_HEADROOM_TTL_SECONDS: float = 10.0
_headroom_cache: dict = {'at': None, 'bytes': None}


def classification_dispatch_window(
        config: Config,
        per_task_bytes: int,
        read_store=None,
        read_self=None,
        ) -> int:
    """How many para_classify tasks to keep in flight, from the memory actually left.

    The window used to be config.proc_n, so 255 tasks were offered whatever they cost.
    The 2026-09-24 run showed what that means: tasks of 10.6 GiB each against 1220 GiB of
    real headroom, so the node could hold ~115 and the OOM killer found that out for us.

    Headroom is the budget ray_init sized everything against, less what the object store
    and this process are actually holding - both measurable, and both far below the
    allowances ray_init subtracts, so measuring recovers the difference. Falls back to
    ray's own task pool, and then to proc_n, when either reading fails; an unknown
    headroom must not stall the stage.
    """
    procs = max(1, config.proc_n)
    if per_task_bytes <= 0:
        return procs
    headroom = classification_headroom_bytes(config, read_store, read_self)
    if not headroom or headroom <= 0:
        return procs
    return max(1, min(procs, int(headroom // per_task_bytes)))


def plan_classification_packages(
        costs: list[int],
        estimates: list[int],
        structs_thresh: int,
        budget: int | None,
        in_flight_bytes: int,
        slots: int,
        fixed_bytes: int = 0,
        ceiling_bytes: int | None = None,
        ) -> tuple[list[list[int]], str]:
    """Which queue entries to send now, as packages of indices, and why it stopped.

    'costs' is the cost-descending queue and 'estimates' the memory each entry is
    expected to need, aligned with it. Returns packages of indices into both.

    Taking the queue strictly from the heavy end - what this replaced - offers the
    proc_n heaviest packages first. On a whole-PDB chunk that asked for 2870 GiB against
    a 972 GiB pool: ray refused most of them, 136 of 255 cores idled, and the node still
    went over because the estimate behind a reservation is not exact. The queue is 31k
    proteins averaging 0.8 GiB estimated, so pairing a few heavy packages with many
    light ones fills every core at a fraction of the memory.

    So: build the next package from the heavy end while the heaviest still fits the
    remaining budget, and from the light end when it does not. Stop when nothing fits,
    since only a returning task can free room - except with nothing running at all,
    where the heaviest has to go out over budget, a protein being indivisible.

    ceiling_bytes, when set, is the most a package is charged, whatever its estimates add
    up to (ClassificationMemoryModel.ceiling): a worker runs its proteins one after another,
    and its peak does not grow with the package beyond that.
    """
    def charged(package_bytes):
        return package_bytes if ceiling_bytes is None else min(package_bytes, ceiling_bytes)

    packages: list[list[int]] = []
    head = 0
    tail = len(costs) - 1
    stop_reason = 'queue drained'

    while head <= tail and len(packages) < slots:
        remaining = None if not budget else max(0, budget - in_flight_bytes)

        from_head = True
        if remaining is not None:
            cheapest = charged(fixed_bytes + estimates[tail])
            heaviest = charged(fixed_bytes + estimates[head])
            # Taking the heaviest must not eat the room the still-empty slots need, or
            # the heavy end fills the budget and the remaining cores get nothing - with
            # a 972 GiB budget and genuinely large packages that left 60 of 255 slots
            # used. The light end is the cheapest package available, so that is what the
            # reserve for those slots is priced at.
            slots_after = max(0, slots - len(packages) - 1)
            if heaviest + slots_after * cheapest > remaining:
                from_head = False
            if not from_head and cheapest > remaining:
                if in_flight_bytes > 0:
                    # Only a returning task can free room now.
                    stop_reason = 'budget exhausted'
                    break
                # Nothing is running, so waiting would deadlock. A protein cannot be
                # split, so the heaviest goes out alone and over budget.
                from_head = True

        package: list[int] = []
        package_cost = 0
        # One worker per package, so its fixed cost is charged once here.
        package_bytes = fixed_bytes

        # Packing stops on the cost threshold, which bounds how long a package runs but
        # not how much memory it holds. At the light end costs are tiny, so reaching the
        # threshold swept in hundreds of proteins and built 48 GiB "light" packages -
        # 20 of them filled the whole budget. A package may take its share of what is
        # left, and one protein always fits so packing cannot stall.
        package_cap = None
        if remaining is not None:
            package_cap = remaining // max(1, slots - len(packages))

        while head <= tail:
            prot_index = head if from_head else tail
            # Flush before overshooting, on either bound: a protein that pushes a
            # non-empty package past the cost threshold makes packages of ~2x the
            # intended cost, and one that pushes it past the remaining memory is what
            # the budget exists to prevent.
            if package and (package_cost + costs[prot_index]) > structs_thresh:
                break
            if package and remaining is not None and charged(package_bytes + estimates[prot_index]) > remaining:
                break
            if package and package_cap is not None and charged(package_bytes + estimates[prot_index]) > package_cap:
                break

            package.append(prot_index)
            package_cost += costs[prot_index]
            package_bytes += estimates[prot_index]
            if from_head:
                head += 1
            else:
                tail -= 1

            if package_cost >= structs_thresh:
                break

        if not package:
            break

        packages.append(package)
        in_flight_bytes += charged(package_bytes)

    if len(packages) >= slots and head <= tail:
        stop_reason = 'all slots filled'
    return packages, stop_reason


def classification_dispatch_stalled(sent_before: int, sent_after: int) -> bool:
    """Whether the wait loop should stop asking pack_packages until a task returns.

    A call that dispatched nothing cannot do better until memory is freed, and only a
    returning task frees it; the wait loop polls every 0.1 s, and asking anyway meant 98,082
    fruitless calls and 4.7 h of driver time on one whole-PDB chunk.

    Never with nothing in flight: then no task could return to clear the stall, and with an
    empty cluster the planner sends the heaviest package over budget anyway.
    """
    return sent_after == sent_before and sent_after > 0


class ClassificationMemoryModel:
    """The memory model of a para_classify task, refitted to the peaks the running job reports.

    The static estimate - fixed + per annotation + per interface chain - is a fit to the code of
    09-24. para_classify has since been rewritten to decode residue matrices in batches while
    keeping per-mapping columns across them, and mappings scale with cost (positions x
    annotations), a term the static estimate does not have. On a whole-PDB chunk the costliest
    protein (cost 51.9M) went out estimated at 7.5 GiB while workers reached 19.7 GB.

    A single correction factor made that worse: light tasks finish first, they came in under the
    estimate, and the factor went below 1 and raised the budget - 255 tasks, Ray's OOM killer.

    So the model is refitted instead, on all four drivers (1, annotations, chains, cost), by
    non-negative least squares over the finished tasks, and scaled by a margin: the
    'quantile' of peak / fitted (at least 1). It is only ever used to raise the prior,
    never to lower it - per protein for the variable part, and for the fixed part - so an early
    fit from light tasks can cost concurrency but not memory. It lives for the run, so chunk 1's
    heavy tail teaches the chunks after it.

    What it cannot do is lower an estimate, and the prior overshoots badly since para_classify
    decodes in batches and the interface calculation was limited: on the first chunk of the
    10-06 whole-PDB run the tasks peaked at 0.25x their estimate (median), at no more than
    3.8 GiB whatever their size, while the heaviest proteins went out at ~15 GiB each. The
    budget then held the last of them back until the phase had run two thirds of its course,
    and they set its end. So the whole estimate of a task is also capped, at ceiling(): a
    margin over the highest peak seen, once enough finished tasks that were estimated above it
    peaked below it. That evidence is what makes the cap safe to apply to heavy tasks; a higher
    peak raises it at once, and without the evidence there is no cap.
    """

    MIN_SAMPLES = 100
    REFIT_EVERY = 50
    MAX_SAMPLES = 20000

    def __init__(self, quantile: float = 0.9):
        self.quantile = quantile
        # (annotations, chains, cost, peak_bytes) per finished task.
        self.samples: deque = deque(maxlen=self.MAX_SAMPLES)
        self.coefficients: np.ndarray | None = None   # bytes per (1, annotation, chain, cost unit)
        self.margin = 1.0
        self.since_fit = 0
        self.refits = 0
        # For ceiling(): the highest peak seen, and what the finished tasks were estimated at
        # without the cap. Set from the config by pack_packages, like quantile.
        self.ceiling_margin = 1.5
        self.ceiling_min_samples = 20
        self.max_peak = 0
        self.uncapped_estimates = np.zeros(1024, dtype=np.int64)
        self.n_uncapped_estimates = 0
        self._ceiling_key: tuple | None = None
        self._ceiling: int | None = None

    def observe(self, annotations: int, chains: int, cost: int, peak_bytes: int,
                uncapped_estimate: int = 0) -> bool:
        """Record a finished task; returns whether this triggered a refit.

        uncapped_estimate is what the task was estimated at without ceiling(), the evidence
        the ceiling is activated on.
        """
        if not peak_bytes or peak_bytes <= 0:
            return False
        self.max_peak = max(self.max_peak, int(peak_bytes))
        if uncapped_estimate and uncapped_estimate > 0:
            if self.n_uncapped_estimates == self.uncapped_estimates.shape[0]:
                self.uncapped_estimates = np.concatenate(
                    (self.uncapped_estimates, np.zeros_like(self.uncapped_estimates)))
            self.uncapped_estimates[self.n_uncapped_estimates] = int(uncapped_estimate)
            self.n_uncapped_estimates += 1
        self.samples.append((annotations or 0, chains or 0, cost or 0, peak_bytes))
        self.since_fit += 1
        if len(self.samples) >= self.MIN_SAMPLES and self.since_fit >= self.REFIT_EVERY:
            self.refit()
            return True
        return False

    def refit(self) -> None:
        from scipy.optimize import nnls

        data = np.array(self.samples, dtype=np.float64)
        features = np.column_stack([np.ones(len(data)), data[:, 0], data[:, 1], data[:, 2]])
        peaks = data[:, 3]
        # Columns span 1 to ~1e8; scaled for the solver, unscaled after.
        scale = features.max(axis=0)
        scale[scale == 0] = 1.0
        scaled_coefficients, _residual = nnls(features / scale, peaks)
        coefficients = scaled_coefficients / scale
        fitted = features @ coefficients
        usable = fitted > 0
        if not usable.any():
            return
        ratios = peaks[usable] / fitted[usable]
        self.coefficients = coefficients
        self.margin = max(1.0, float(np.quantile(ratios, self.quantile)))
        self.since_fit = 0
        self.refits += 1

    def fixed_bytes(self, prior_fixed: int) -> int:
        if self.coefficients is None:
            return int(prior_fixed)
        return int(max(prior_fixed, self.margin * self.coefficients[0]))

    def variable_bytes(self, prior_variable, annotations, chains, costs):
        """Per-protein (or per-package) variable estimate; arrays or scalars alike."""
        if self.coefficients is None:
            return prior_variable
        _fixed, per_annotation, per_chain, per_cost = self.coefficients
        fitted = self.margin * (per_annotation * np.asarray(annotations, dtype=np.float64)
                                + per_chain * np.asarray(chains, dtype=np.float64)
                                + per_cost * np.asarray(costs, dtype=np.float64))
        result = np.maximum(np.asarray(prior_variable, dtype=np.float64), fitted)
        if np.ndim(result) == 0:
            return int(result)
        return result.astype(np.int64)

    def ceiling(self) -> int | None:
        """The most any one task is estimated at, or None while there is no cap.

        ceiling_margin times the highest peak seen, once at least ceiling_min_samples finished
        tasks had been estimated at that or more without the cap: those tasks are the evidence
        that estimates above it overshoot. A margin of 0 switches the cap off.
        """
        if self.ceiling_margin <= 0 or self.max_peak <= 0:
            return None
        key = (self.n_uncapped_estimates, self.max_peak, self.ceiling_margin, self.ceiling_min_samples)
        if key != self._ceiling_key:
            candidate = int(self.ceiling_margin * self.max_peak)
            n_above = int(np.count_nonzero(
                self.uncapped_estimates[:self.n_uncapped_estimates] >= candidate))
            self._ceiling = candidate if n_above >= self.ceiling_min_samples else None
            self._ceiling_key = key
        return self._ceiling

    def cap(self, estimate):
        """estimate (bytes of a whole task; array or scalar) under ceiling()."""
        ceiling = self.ceiling()
        if ceiling is None:
            return estimate
        if np.ndim(estimate) == 0:
            return int(min(int(estimate), ceiling))
        return np.minimum(np.asarray(estimate, dtype=np.int64), ceiling)

    def describe_ceiling(self) -> str:
        ceiling = self.ceiling()
        if ceiling is not None:
            return f'ceiling {ceiling / 1024**3:.1f} GiB'
        if self.ceiling_margin <= 0:
            return 'no ceiling (off)'
        return 'no ceiling yet'

    def describe(self) -> str:
        if self.coefficients is None:
            return (f'static prior ({len(self.samples)} samples, refit from {self.MIN_SAMPLES}), '
                    f'{self.describe_ceiling()}')
        fixed, per_annotation, per_chain, per_cost = self.coefficients
        return (f'refit {self.refits} on {len(self.samples)} tasks: '
                f'{fixed / 1024**3:.2f} GiB + {per_annotation / 1024**2:.3f} MB/annotation + '
                f'{per_chain / 1024**2:.4f} MB/chain + {per_cost:.1f} B/cost, margin {self.margin:.2f}, '
                f'{self.describe_ceiling()}')


_memory_model = ClassificationMemoryModel()

# _v2: the run and chunk columns were added, since runs sharing an output folder write into
# one table at the same time; the first table had to be taken apart by the structure of its
# estimates. A new name rather than new columns in the old file, which would mix two layouts.
MEMORY_SAMPLES_FILE = 'classification_memory_samples_v2.tsv'
_MEMORY_SAMPLES_COLUMNS = ('run', 'chunk', 'proteins', 'cost', 'annotations', 'interface_chains',
                           'estimate_bytes', 'uncapped_estimate_bytes', 'peak_uss_bytes',
                           'elapsed_seconds')


def record_memory_sample(config: Config, record: dict, peak_uss: int) -> None:
    """Append one finished task to the output folder's calibration table, at any verbosity.

    The per-task done line needs verbosity 4, so no run since the batched rewrite of para_classify
    had the records a refit of the static prior needs. A few KB per thousand tasks; failures to
    write are ignored - this must never be the reason a task's results are lost.
    """
    try:
        path = os.path.join(config.outfolder, MEMORY_SAMPLES_FILE)
        new_file = not os.path.exists(path)
        with open(path, 'a') as samples:
            if new_file:
                samples.write('\t'.join(_MEMORY_SAMPLES_COLUMNS) + '\n')
            samples.write('\t'.join(str(value) for value in (
                getattr(config, 'run_id', ''), getattr(config, 'current_chunk', ''),
                record.get('n_proteins'), record.get('cost'), record.get('annotations'),
                record.get('chains'), record.get('estimate'),
                record.get('uncapped_estimate', record.get('estimate')), int(peak_uss),
                round(record.get('elapsed') or 0.0, 1))) + '\n')
    except Exception:
        pass


def classification_memory_budget(config: Config, read_store=None, read_self=None) -> int | None:
    """Estimated bytes the dispatcher may keep in flight, or None for no budget.

    The same headroom the window is derived from. Separate name because the dispatcher
    packs against the total rather than dividing it by a per-task figure: dividing was
    what assumed every package costs the same, and the real ones ranged from 3.8 GiB
    (assumed) to 40 GiB.
    """
    return classification_headroom_bytes(config, read_store, read_self)


def classification_headroom_bytes(config: Config, read_store=None, read_self=None) -> int | None:
    """Bytes the para_classify tasks may hold between them, or None when unknown.

    The budget ray_init sized everything against, less what the object store and this
    process actually hold - both measurable, and both far below the allowances ray_init
    subtracts, so measuring recovers the difference. Falls back to ray's own task pool.
    """
    if read_store is None:
        read_store = read_object_store_stats
    if read_self is None:
        read_self = process_memory_breakdown

    headroom = None
    effective = getattr(config, 'ray_effective_memory_bytes', None)
    if effective:
        # read_object_store_stats is an RPC to the raylet and the wait loop asks about ten
        # times a second, so the measured part is cached briefly. Neither figure moves
        # fast: the store is filled before classification starts and the driver's own
        # footprint is flat through the stage.
        now = time.time()
        if _headroom_cache['at'] is None or (now - _headroom_cache['at']) > _HEADROOM_TTL_SECONDS:
            store_stats = read_store()
            driver = read_self()
            if store_stats is not None and driver is not None:
                _headroom_cache['bytes'] = (
                    effective - int(store_stats.get('used', 0)) - int(driver['rss']))
            else:
                _headroom_cache['bytes'] = None
            _headroom_cache['at'] = now
        headroom = _headroom_cache['bytes']
    if headroom is None:
        # ray_init already subtracted generous allowances for exactly these two, so its
        # pool is the conservative version of the same number.
        headroom = getattr(config, 'ray_task_memory_bytes', None)
    if not headroom or headroom <= 0:
        return None
    return int(headroom)


@ray.remote(max_calls=1)
def para_classify_remote_wrapper(
        package: list[tuple[str, int, int]],
        stores: tuple,
        package_cost: int | None = None,
        memory_request: int | None = None,
        annotations: int | None = None,
        interface_chains: int | None = None,
        ):
    """package_cost and memory_request only feed the done line below: peak RSS against
    cost is what classification_gib_per_mcost is calibrated from, and worker lines
    interleave too much to pair them up any other way."""

    t_start = time.time()

    config: Config = ray.get(stores[0])
    prot_db_id_map: tuple[np.ndarray[int], np.ndarray[int]] = ray.get(stores[1])
    structs_from_db: np.ndarray[int] = ray.get(stores[2])
    matrix_stores: list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]] = stores[3]
    backmap_store: tuple[np.ndarray[int], np.ndarray[int], np.ndarray[int]] = ray.get(stores[4])
    structure_db_list: np.ndarray[bytes] = ray.get(stores[5])
    sequence_storage: tuple[np.ndarray, np.ndarray] = ray.get(stores[6])
    compl_struct_db_id_matrix: np.ndarray[int] = ray.get(stores[7])
    compl_db_id_map: np.ndarray[int] = ray.get(stores[8])
    known_empty_complexes: np.ndarray = ray.get(stores[9])

    reset_logger_for_remotes(config)

    # Always sampled, whatever the log level: the driver calibrates its memory estimate from
    # these peaks (ClassificationMemoryModel). Gated on verbosity 4 before, the driver
    # never learned that the batched rewrite of para_classify had moved the peak, and a whole-PDB
    # chunk ran ~15 GB tasks against ~12.5 GiB estimates into Ray's OOM killer.
    peak_tracker = PrivateMemoryPeak().start()
    if config.verbosity >= 4:
        rss0 = psutil.Process().memory_info().rss / 1024**2
        setup_time = time.time() - t_start
        config.logger.info(f'[para_classify][worker] start: {len(package)} proteins, setup={setup_time:.1f}s, rss={rss0:.0f} MB, avail={psutil.virtual_memory().available / 1024**2:.0f} MB')

    result = (para_classify(config,
                          package,
                          prot_db_id_map,
                          structure_db_list,
                          sequence_storage,
                          structs_from_db,
                          matrix_stores,
                          backmap_store,
                          compl_struct_db_id_matrix,
                          compl_db_id_map,
                          known_empty_complexes,
                          para=True))

    elapsed = time.time() - t_start
    peak_tracker.stop()
    if config.verbosity >= 4:
        rss1 = psutil.Process().memory_info().rss / 1024**2
        # ru_maxrss is KiB on Linux. The peak, not the final RSS, is what a
        # reservation has to cover.
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        reserved = f'{memory_request / 1024**2:.0f} MB' if memory_request else 'none'
        protein_ids = ','.join(str(entry[0]) for entry in package[:5]) + (',...' if len(package) > 5 else '')
        # uss is what one more task costs and therefore what the reservation has to
        # cover; pss summed over the workers is what the node pays. rss is neither: a
        # worker maps the backmap and residue stores and its rss grows as it touches
        # them, which is shared and charged once. A run calibrated 0.46 GiB per Mcost
        # from rss and then OOMed at 255 tasks.
        breakdown = process_memory_breakdown()
        if breakdown is None:
            split = 'uss=n/a pss=n/a'
        else:
            split = (f'uss={breakdown["uss"] / 1024**2:.0f} MB '
                     f'pss={breakdown["pss"] / 1024**2:.0f} MB '
                     f'shared={breakdown["shared"] / 1024**2:.0f} MB')
        peaks = peak_tracker.report()
        peaks = f'{peaks}, {residue_buffer_report()}, {decode_cache_report()}'
        config.logger.info(
            f'[para_classify][worker] done: {len(package)} proteins, elapsed={elapsed:.1f}s, '
            f'rss={rss1:.0f} MB (+{rss1-rss0:.0f} MB), peak_rss={peak:.0f} MB, {split}, {peaks}, '
            f'cost={package_cost}, annotations={annotations}, '
            f'interface_chains={interface_chains}, reserved={reserved}, '
            f'proteins={protein_ids}')
    # The peak rides along as a fourth element; the wait loop also accepts the old three.
    # 0 means the sampler could not read smaps_rollup, which the calibration skips.
    return tuple(result) + (int(peak_tracker.peak_uss),)


def para_classify(
        config: Config,
        package: list[tuple[str, int, list[int], list[tuple[int, str]]]],
        prot_db_id_map: tuple[np.ndarray, np.ndarray],
        structure_db_list: np.ndarray[bytes],
        sequence_storage: tuple[np.ndarray, np.ndarray],
        structs_from_db: np.ndarray[int],
        matrix_stores: list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]],
        backmap_store: tuple[np.ndarray[int], np.ndarray[int], np.ndarray[int]],
        compl_struct_db_id_matrix,
        compl_db_id_map,
        known_empty_complexes: np.ndarray | None = None,
        para: bool= False):
    ta = time.time()
    times = []

    # Stage markers for PrivateMemoryPeak. They cost one list store each and are what
    # names the part of the footprint the residue batching does not explain.
    set_task_phase('setup')

    if config.verbosity >= 4:
        config.logger.info(f'Start of para_classify {len(package)=} {para=} {backmap_store is None=}')

    if config.verbosity >= 5:
        for name, size in sorted(((name, deep_get_size_of(value)) for name, value in locals().items()), key=lambda x: -x[1])[:10]:
            config.logger.info(f"Point 1: {name:>30}: {sizeof_fmt(size):>8}")

    prot_db_id_map_array, pos_db_id_vec = prot_db_id_map
    prot_db_id_map = package_prot_db_id_map(prot_db_id_map_array, package)

    ta = add_to_times(times, ta, 'retrieve res data')

    # Compressed per-package blobs (see retrieval.select_package_from_db_shard). ray.get
    # maps them zero-copy, so holding every package here costs no heap; each protein
    # decodes only its own structures in build_local_stores.
    blob_stores: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    for p_res_buffer, p_res_offsets, p_iface_buffer, p_iface_offsets, p_resolutions in matrix_stores:
        blob_stores.append((ray.get(p_res_buffer), ray.get(p_res_offsets),
                            ray.get(p_iface_buffer), ray.get(p_iface_offsets),
                            ray.get(p_resolutions)))

    decode_cache = DecodedStructureCache(
        int(getattr(config, 'aggregation_decode_cache_mb', 256) * 1024**2))
    need_interfaces = config.compute_ppi and not config.fast_pdb_annotation

    ta = add_to_times(times, ta, f'complete structure and interface maps {len(structs_from_db)=}')

    prot_db_ids = [p[1] for p in package]

    prot_structure_alignment_map = get_alignment_data(config, prot_db_ids)

    ta = add_to_times(times, ta, f'retrieve alignment data {len(prot_db_ids)=} {len(prot_structure_alignment_map)=} {config.provide_msas=}')

    total_sub_times = []
    position_values = [[] for _ in range(20)]

    if config.verbosity >= 5:
        for name, size in sorted(((name, deep_get_size_of(value)) for name, value in locals().items()), key=lambda x: -x[1])[:10]:
            config.logger.info(f"Point 2: {name:>30}: {sizeof_fmt(size):>8}")

    total_nr_of_structs = 0

    stats_list = []

    # Hoisted out of the per-mapping loop, which runs billions of times.
    profile_left, profile_right = residue_package.PROFILE_SPAN
    sc_rsa_column = residue_package.SC_RSA_COLUMN
    aa_column = residue_package.AA_COLUMN
    plddt_column = residue_package.PLDDT_COLUMN
    # The three columns integrated_features.generate counts locations over, and whose
    # side-chain one it takes np.std of. Retained per position by the batched path,
    # because that std is two-pass and so cannot be accumulated.
    rsa_column = residue_package.RSA_COLUMN
    mc_rsa_column = residue_package.MC_RSA_COLUMN
    surface_threshold = config.surface_threshold
    buried_threshold = config.buried_threshold
    # Maps an amino acid letter onto the code the encoded residue stores.
    aa_string_encodes = residue_package.AA_STRING_ENCODES

    for prot_index, (protein_id, prot_db_id, seq_index) in enumerate(package):
        if prot_db_id not in prot_structure_alignment_map:
            continue

        t_start = time.time()
        set_task_phase('alignment columns')
        seq = sequence_from_store(sequence_storage, seq_index)
        
        alignments_for_protein = prot_structure_alignment_map[prot_db_id]
        annotation_list: list[int] = list(alignments_for_protein.keys())

        # Everything the per-mapping loop reads that depends only on the structure,
        # gathered once per protein instead of once per (position, structure) pair -
        # that loop runs len(seq) times over this same list. Held as plain python
        # ints and strs, not numpy scalars, for the same reason annotation_list is:
        # boxed numpy scalars slow down every use downstream.
        anno_sub_infos, anno_cov, anno_seq_id = _gather_alignment_columns(
            alignments_for_protein, annotation_list)
        sub_infos_concat, sub_infos_bounds = _flatten_alignments(anno_sub_infos)
        cov_column, seq_id_column, cov_is_none, seq_id_is_none = _alignment_value_columns(
            anno_cov, anno_seq_id)
        # Only these take part in the interface calculation; the per-position
        # classification below keeps using every annotation.
        if need_interfaces:
            interface_annotations, n_multi_chain = interface_annotation_list(
                annotation_list, seq_id_column, compl_db_id_map, compl_struct_db_id_matrix, config)
            if config.verbosity >= 4:
                config.logger.info(
                    f'In para_classify: {protein_id=} interface annotations: '
                    f'{len(interface_annotations)} selected of {n_multi_chain} multi chain, '
                    f'{len(annotation_list)} in total')
        else:
            interface_annotations = []
        # The interface entries of every chain of every complex the protein's
        # interface annotations belong to, decoded once. Only the residue data is
        # batched below: these keys measure ~22 MiB even for a heavy whole-PDB protein,
        # and calculate_aggregated_interface_map walks all of its complexes in one pass.
        set_task_phase('interface stores')
        (local_rows, index_list_stores, chain_list_stores,
         concat_combined_keys_stores, concat_scores_stores) = build_interface_stores(
            interface_structure_ids(interface_annotations, compl_db_id_map, compl_struct_db_id_matrix),
            structs_from_db, blob_stores, decode_cache)

        # Whole-protein, because MappingIds indexes them by an annotation index that
        # spans every batch, and they are strings held by reference.
        anno_structure_id, anno_chain = _gather_identity_columns(annotation_list, structure_db_list)
        anno_structure_rank = _structure_rank_column(anno_structure_id)

        n_annotations_total = len(annotation_list)
        # The annotations are decoded in batches, so a task holds one batch's residue
        # matrix rather than the protein's whole one - 8.6 GiB at the top end, which is
        # what held classification to ~70 of 255 cores.
        batch_size = max(1, int(getattr(config, 'classification_batch_size', 1000)))

        # The kernel's output columns, sized to a batch and refilled for every position.
        out_width = max(1, min(batch_size, n_annotations_total))
        out_anno = np.empty(out_width, dtype=np.int64)
        out_res_contig = np.empty(out_width, dtype=np.int64)
        out_row = np.empty(out_width, dtype=np.int64)
        out_store = np.empty(out_width, dtype=np.int64)
        out_qual = np.empty(out_width, dtype=np.float64)
        out_seq_id = np.empty(out_width, dtype=np.float64)
        out_cov = np.empty(out_width, dtype=np.float64)
        out_aa_id = np.empty(out_width, dtype=np.bool_)
        out_resolution = np.empty(out_width, dtype=np.float64)
        out_res_aa = np.empty(out_width, dtype=np.float32)
        out_simple_class = np.empty(out_width, dtype=np.int64)
        out_no_seq_id = np.empty(out_width, dtype=np.int64)
        out_bad_contig = np.empty(out_width, dtype=np.int64)
        # Accumulated across the whole protein: the (seq_id, cov) pair a structure
        # registers does not depend on which position, or which batch, first reached it.
        # The kernel is handed the batch's slice of this, which writes through.
        quality_measure_seen = np.zeros(n_annotations_total, dtype=np.uint8)
        counts_array = np.zeros(10, dtype=np.int64)

        cost = len(seq)*len(annotation_list)

        ts = time.time()
        sub_times = []
        #annotation_list: list[tuple[str, str]] = protein_obj.get_annotation_list()

        #config.logger.info(f'{protein_id} {prot_db_id=} {prot_db_id%20=} {annotation_list=}')

        total_nr_of_structs += len(annotation_list)

        position_value_dict = {}
        
        recommended_complexes: list[list[tuple[str, str, int | str]]] = [None]

        #prot_db_id: int = protein_obj.database_id
        
        if config.provide_msas:
            compl_fetched = set()
            for struct_db_id in annotation_list:
                fused_id = structure_db_list[struct_db_id] 
                structure_id = fused_id[:-1].decode()
                if structure_id in compl_fetched:
                    continue
                compl_fetched.add(structure_id)
                if is_alphafold_model(structure_id):
                    try:
                        # The map's entries are (alignment, coverage, seq_id); index 1
                        # is the coverage, and fetch_msa wants the sequence identity.
                        seq_id = prot_structure_alignment_map[prot_db_id][struct_db_id][2]
                    except KeyError as err:
                        continue
                        #raise KeyError(f'{err=} {prot_db_ids=} {protein_id=} {prot_db_id=} {structure_id=} {chain=} {struct_db_id=} {prot_structure_alignment_map=}')
                    fetch_msa(config, structure_id, protein_id, prot_db_id, seq_id, 1)

        ts = add_to_times(sub_times, ts, 'fetch_msa')

        #positions: list[int] = protein_obj.get_position_ids()

        structure_quality_measures: dict[int, tuple[float, float]] = {}

        if config.verbosity >= 4:
            config.logger.info(f"In para_classify: {protein_id=} {prot_db_id=} {prot_db_id%20=} {len(annotation_list)=} {len(seq)=}")

        total_pos_loop_times = []

        # Per-protein counters for the aggregation hot loop. Scalars only, updated once
        # per position (never per annotation), so they cost nothing measurable but tell
        # us the iteration counts and skip reasons that the timing labels alone cannot.
        prot_counts: list[int] = [0]*10
        prot_continues: int = 0
        sum_mappings: int = 0
        min_mappings: int = -1
        max_mappings: int = 0

        #decoded_residue_objs = {}

        # One Mappings per position, filled batch by batch and finalised afterwards. The
        # loops are nested the other way round from before - batches outside, positions
        # inside - so that a batch's residue matrix is dropped before the next is decoded.
        # The fixed cost is 13 KB a position (a bare Mappings plus its accumulator,
        # measured), so 22 MB for a 1800-residue protein. What grows with the mappings is
        # the thin per-mapping columns the recommendation needs, at 73 bytes each against
        # the 780 of an encoded row: a mapping exists per (position, residue) pair and a
        # matrix row per residue, so both scale with the same product and the retained
        # columns are a fixed ~9% of the matrix they replace.
        mappings_per_pos: list[mappings_package.Mappings | None] = [
            mappings_package.Mappings() for _ in range(len(seq))]
        # count_continues is per position and now arrives one batch at a time.
        per_pos_continues = np.zeros(len(seq), dtype=np.int64)
        # Mappings actually collected per position, over all batches: 0 means no structure
        # covers the position, which is normal and not worth a warning (see the finalize loop).
        per_pos_collected = np.zeros(len(seq), dtype=np.int64)
        # Mappings to a residue past those stored for their structure, reported once per
        # protein after the batches (record_bad_contig / summarize_bad_contigs).
        bad_contigs: dict = {}

        for batch_lo in range(0, n_annotations_total, batch_size):
            batch_hi = min(batch_lo + batch_size, n_annotations_total)
            batch_ids = annotation_list[batch_lo:batch_hi]

            set_task_phase('residue batch')
            batch_rows, residue_stores, resolutions_stores = build_residue_batch(
                batch_ids, structs_from_db, blob_stores, decode_cache)
            (anno_store_index, anno_res_left, anno_res_right,
             anno_resolution) = _gather_residue_columns(batch_ids, batch_rows, resolutions_stores)

            # Views, so the kernel's anno indices stay 0-based within the batch while
            # quality_measure_seen is still the protein's own array. sub_infos_bounds
            # holds offsets into the whole concatenation, so that is passed unsliced.
            batch_bounds = sub_infos_bounds[batch_lo:batch_hi]
            batch_cov_column = cov_column[batch_lo:batch_hi]
            batch_seq_id_column = seq_id_column[batch_lo:batch_hi]
            batch_cov_is_none = cov_is_none[batch_lo:batch_hi]
            batch_seq_id_is_none = seq_id_is_none[batch_lo:batch_hi]
            batch_seen = quality_measure_seen[batch_lo:batch_hi]
            batch_ranks = anno_structure_rank[batch_lo:batch_hi]

            set_task_phase('position loop')
            for pos, aa1 in enumerate(seq, start = 1):

                t_pos_loop = time.time()
                pos_loop_times = []

                # The encoded residue carries the numeric aa code, not the letter,
                # so the query residue has to be translated into the same code space
                # once per position. Comparing the code against the sequence
                # character made identical_aa constantly False, which left the
                # identical-aa recommendation buckets permanently empty.
                aa1_code = aa_string_encodes.get(aa1)

                # The per-mapping loop runs len(seq) x len(annotation_list) times -
                # billions over a run - so it lives in a kernel and this is the only
                # python step per position and batch. The columns it fills are reused;
                # n_collected says how much of each is live.
                counts_array[:] = 0
                n_collected, count_continues, n_no_seq_id, n_bad_contig = _collect_mappings_for_position(
                    pos, np.int64(-1 if aa1_code is None else aa1_code), aa1_code is not None,
                    sub_infos_concat, batch_bounds,
                    batch_cov_column, batch_seq_id_column, batch_cov_is_none, batch_seq_id_is_none,
                    anno_store_index, anno_res_left, anno_res_right, anno_resolution,
                    residue_stores, aa_column, plddt_column, sc_rsa_column,
                    profile_left, profile_right, CLASS_SCORE_INDICES,
                    surface_threshold, buried_threshold, _INVERSE_POWER,
                    out_anno, out_res_contig, out_row, out_store,
                    out_qual, out_seq_id, out_cov, out_aa_id, out_resolution, out_res_aa,
                    out_simple_class, batch_seen,
                    out_no_seq_id, out_bad_contig, counts_array)

                counts = counts_array.tolist()

                per_pos_collected[pos - 1] += n_collected

                for warning_index in range(n_bad_contig):
                    anno_index = int(out_bad_contig[warning_index])
                    protein_anno_index = batch_lo + anno_index
                    res_contig = int(sub_infos_concat[sub_infos_bounds[protein_anno_index, 0] + pos])
                    record_bad_contig(
                        bad_contigs, annotation_list[protein_anno_index],
                        f'{anno_structure_id[protein_anno_index]}:{anno_chain[protein_anno_index]}',
                        pos, res_contig, int(anno_res_right[anno_index]) - int(anno_res_left[anno_index]))

                if config.verbosity >= 4:
                    for skip_index in range(n_no_seq_id):
                        struct_db_id = annotation_list[batch_lo + int(out_no_seq_id[skip_index])]
                        config.logger.info(f'Skipped classification of {protein_id} due to {struct_db_id} got no sequence identity')

                live = slice(0, n_collected)

                # One gather, instead of the list of row views that weight_all then had
                # to copy into a matrix on every position. Rows almost always come from a
                # single store, which fancy indexing handles in one go; the mixed case
                # falls back to a row-by-row copy.
                if n_collected == 0:
                    encoded_matrix = np.empty((0, residue_stores[0].shape[1]), dtype=np.float32)
                else:
                    first_store = out_store[0]
                    if bool((out_store[live] == first_store).all()):
                        encoded_matrix = residue_stores[first_store][out_row[live]]
                    else:
                        encoded_matrix = np.empty(
                            (n_collected, residue_stores[first_store].shape[1]), dtype=np.float32)
                        for mapping_nr in range(n_collected):
                            encoded_matrix[mapping_nr] = residue_stores[
                                out_store[mapping_nr]][out_row[mapping_nr]]

                # The batch's rows are folded into this position's accumulator and the
                # matrix is dropped with the batch. out_anno is batch-local, so the batch
                # offset goes on before the index is kept.
                mappings_per_pos[pos - 1].accumulate_encoded_mappings(
                    encoded_matrix, out_qual[live].copy(), out_seq_id[live].copy(),
                    out_cov[live].copy(), out_aa_id[live].copy(), out_resolution[live].copy(),
                    out_res_aa[live].copy(), out_simple_class[live].copy(),
                    batch_ranks[out_anno[live]],
                    (out_anno[live] + batch_lo).astype(np.int32),
                    out_res_contig[live].astype(np.int32),
                    np.ascontiguousarray(encoded_matrix[:, rsa_column]),
                    np.ascontiguousarray(encoded_matrix[:, mc_rsa_column]),
                    np.ascontiguousarray(encoded_matrix[:, sc_rsa_column]),
                    config)

                per_pos_continues[pos - 1] += count_continues
                prot_continues += count_continues
                for _count_index in range(10):
                    prot_counts[_count_index] += counts[_count_index]

                t_pos_loop = add_to_times(pos_loop_times, t_pos_loop, 'pos aggregation loop')
                total_pos_loop_times = aggregate_times(total_pos_loop_times, pos_loop_times)

            # The batch's residue matrix is dead here; the next batch allocates its own.
            del residue_stores, resolutions_stores, batch_rows
            del anno_store_index, anno_res_left, anno_res_right, anno_resolution

        if bad_contigs:
            config.errorlog.add_warning(summarize_bad_contigs(protein_id, bad_contigs))
        del bad_contigs

        # Every batch has contributed, so each position's aggregate can be closed out.
        set_task_phase('finalize positions')
        for pos, aa1 in enumerate(seq, start = 1):
            t_pos_loop = time.time()
            pos_loop_times = []

            mappings_obj = mappings_per_pos[pos - 1]
            mappings_per_pos[pos - 1] = None
            count_continues = int(per_pos_continues[pos - 1])

            # Not n_collected: this statistic has always counted every annotation the
            # loop did not `continue` out of via a counted branch, which includes the
            # res_contig ones that never became a mapping.
            n_mappings = len(annotation_list) - count_continues
            sum_mappings += n_mappings
            if n_mappings > max_mappings:
                max_mappings = n_mappings
            if min_mappings < 0 or n_mappings < min_mappings:
                min_mappings = n_mappings

            if len(annotation_list) > 0:
                try:
                    err, weighting_times = mappings_obj.finalize_encoded(
                        config, anno_structure_id, anno_chain)
                except AttributeError as err:
                    weighting_times = []
                    if per_pos_collected[pos - 1] == 0:
                        # No structure covers this position (alignment gaps, structures
                        # without residue data): nothing to weight, so the rin profile stays
                        # unset and weight_all raises - by design, and the position stays
                        # unclassified as it should. Not a warning: one stack trace per
                        # uncovered position made up 95% of a 1.19 GB warnlog.
                        pass
                    elif len(annotation_list) <= count_continues+1:
                        config.errorlog.add_warning(f'All mappings failed: {err=} {protein_id=} {prot_db_id=} {pos=} {len(annotation_list)=} {count_continues=}')
                    else:
                        raise AttributeError(f'weight all failed{err=} {protein_id=} {pos=} {len(annotation_list)=} {count_continues=}')

                pos_loop_times.append(weighting_times)
            else:
                pos_loop_times.append([])

            t_pos_loop = add_to_times(pos_loop_times, t_pos_loop, 'weight_all')

            if config.verbosity >= 7:
                config.logger.info(f"In para_classify: {protein_id=} after weighting of {pos=} {mappings_obj.rin_based_features.profile=}")
                mappings_obj.structural_features.log_attr_sizes(config.logger, print_values=True, label=f'{pos=} ')
                mappings_obj.integrated_features.log_attr_sizes(config.logger, print_values=True, label=f'{pos=} ')

            recommendation_order = mappings_obj.recommendation_order
            mappings_obj.recommendation_order = None

            try:
                rec_res = mappings_obj.recommended_res
                max_res = mappings_obj.max_seq_res
            except AttributeError:
                rec_res = None
                max_res = None

            if config.verbosity >= 6:
                config.logger.info(f"In para_classify: {protein_id=} before packing recommended struct of {pos=} {mappings_obj.mapping_ids=}")
                config.logger.info(f"In para_classify: {protein_id=} before packing recommended struct of {pos=} {rec_res=} {max_res=}")

            packed_recs = pack((rec_res, max_res))

            pos_db_id_vec_left, _ = prot_db_id_map[prot_db_id]
            pos_db_id = pos_db_id_vec[(pos_db_id_vec_left + pos) -1]
            position_value_dict[pos] = [pos_db_id, packed_recs, pack(mappings_obj) ,None]
            del mappings_obj

            recommended_complexes.append(recommendation_order)

            t_pos_loop = add_to_times(pos_loop_times, t_pos_loop, 'pos loop wrap up')

            total_pos_loop_times = aggregate_times(total_pos_loop_times, pos_loop_times)


        # Built from the kernel's per-structure flag rather than inside the loop: the
        # (seq_id, cov) pair a structure registers is the same whichever position
        # first reached it, so only the fact that some position did matters. The
        # values stay the python floats the alignment map held - calc_interface_
        # annotation_score multiplies them into a float32 nansum, and a numpy float64
        # there would promote the whole product to float64.
        for anno_index in range(n_annotations_total):
            if quality_measure_seen[anno_index]:
                structure_quality_measures[annotation_list[anno_index]] = (
                    anno_seq_id[anno_index], anno_cov[anno_index])

        sub_times.append(total_pos_loop_times)

        ts = add_to_times(sub_times, ts, 'main aggregation loop')

        if config.verbosity >= 6:
            config.logger.info(f"In para_classify: {protein_id=} before calc_mm_features")

        mm_db_dict: dict[str, dict[str, dict[str | int, bytes]]]
        set_task_phase('mm features')
        mm_counters: dict = {}
        mm_db_dict, successful_mm_annotations = calc_mm_features(
            config, recommended_complexes, known_empty_complexes, counters=mm_counters)

        if config.verbosity >= 6:
            config.logger.info(f"In para_classify: {protein_id=} after calc_mm_features {len(mm_db_dict)=} {len(successful_mm_annotations)=}")

        ts = add_to_times(sub_times, ts, 'calc mm feats')

        del recommended_complexes

        for pos in successful_mm_annotations:
            complex_id, chain_id, res = successful_mm_annotations[pos]
            if chain_id in mm_db_dict[complex_id]:
                if res in mm_db_dict[complex_id][chain_id]:
                    # Packed here rather than at lookup time: only these positions need it.
                    packed_mm_feats = pack_mm_features(mm_db_dict[complex_id][chain_id][res])
                    if config.verbosity >= 6:
                        config.logger.info(f'Setting mm_feat_vec {pos} {complex_id} {chain_id} {res}')
                    position_value_dict[pos][3] = packed_mm_feats

        ts = add_to_times(sub_times, ts, 'add mm feats to value dicts')

        del mm_db_dict
        del successful_mm_annotations

        for pos in position_value_dict:
            internal_pos_db_id, packed_recs, packed_mappings, packed_mm_feats = position_value_dict[pos]
            pos_db_id = internal_pos_db_id//100
            # pack_frame, not pack: both halves are already zstd frames, so compressing
            # the pair again only cost time and grew the blob.
            position_values[prot_db_id % 20].append((pos_db_id, packed_recs, pack_frame((packed_mappings, packed_mm_feats))))

            #config.logger.info(f'{protein_id=} {prot_db_id=} {pos=} {internal_pos_db_id=} {pos_db_id=} {prot_db_id%20=}')

        ts = add_to_times(sub_times, ts, 'transform to value lists')

        del position_value_dict

        ts = add_to_times(sub_times, ts, f'delete value dict {config.compute_ppi=}')

        if config.compute_ppi and not config.fast_pdb_annotation:
            annotation_dict = {}
            for struct_db_id in interface_annotations:
                fused_id = structure_db_list[struct_db_id]
                structure_id = fused_id[:-1].decode()
                chain = fused_id[-1:].decode()
                if structure_id not in annotation_dict:
                    annotation_dict[structure_id] = []
                annotation_dict[structure_id].append((chain, struct_db_id))

            aggregated_interfaces: list[interface_package.Aggregated_interface]
            set_task_phase('interface map')
            aggregated_interfaces, interface_aggregation_times = interface_package.calculate_aggregated_interface_map(
                config,
                prot_db_id,
                annotation_dict,
                local_rows,
                index_list_stores,
                chain_list_stores,
                concat_combined_keys_stores,
                concat_scores_stores,
                structure_quality_measures,
                compl_struct_db_id_matrix,
                compl_db_id_map,
                structure_db_list,
                backmap_store)
            
            sub_times.append(interface_aggregation_times)

            ts = add_to_times(sub_times, ts, f'calculate agg interface map {len(aggregated_interfaces)=}')

            if config.verbosity >= 6:
                config.logger.info(f'Len of aggregated_interfaces for {protein_id}: {len(aggregated_interfaces)}')

            if config.verbosity >= 5 and prot_index == 0:
                for name, size in sorted(((name, deep_get_size_of(value)) for name, value in locals().items()), key=lambda x: -x[1])[:10]:
                    config.logger.info(f"Point 2.5: {name:>30}: {sizeof_fmt(size):>8}")


            if len(aggregated_interfaces) > 0:
                ins_times = remote_insert_interfaces(prot_db_id, aggregated_interfaces, config)
            else:
                ins_times = []

            sub_times.append(ins_times)
            ts = add_to_times(sub_times, ts, 'insert interfaces')

            del annotation_dict
            del aggregated_interfaces

            ts = add_to_times(sub_times, ts, 'delete interface data')

        n_annotations = len(annotation_list)

        del annotation_list, interface_annotations
        del structure_quality_measures

        # The decoded residue matrix is released batch by batch inside the loop above,
        # together with its per-annotation columns - naming them here again failed every
        # task with an UnboundLocalError. What remains is the protein's interface stores
        # and the columns built once per protein, released here rather than left bound
        # until the next protein's stores are built, which had both alive at once.
        del local_rows, index_list_stores
        del chain_list_stores, concat_combined_keys_stores, concat_scores_stores
        del anno_structure_id, anno_chain, anno_structure_rank
        del sub_infos_concat, sub_infos_bounds, anno_sub_infos
        # get_alignment_data read the whole package's alignments up front; this
        # protein's are done with, so they go now rather than at the end of the task.
        del alignments_for_protein
        prot_structure_alignment_map.pop(prot_db_id, None)

        t_end = time.time()
        n_positions = len(seq)
        stats_list.append((
            protein_id, cost, t_end-t_start,
            n_positions, n_annotations, n_positions * n_annotations,
            prot_continues, sum_mappings, min_mappings, max_mappings,
            tuple(prot_counts),
            tuple(mm_counters.get(name, 0) for name in MM_COUNTER_NAMES),
            ))

        ts = add_to_times(sub_times, ts, 'delete annotation data')

        total_sub_times = aggregate_times(total_sub_times, sub_times)

    if config.verbosity >= 5:
        for name, size in sorted(((name, deep_get_size_of(value)) for name, value in locals().items()), key=lambda x: -x[1])[:10]:
            config.logger.info(f"Point 3: {name:>30}: {sizeof_fmt(size):>8}")

    times.append(total_sub_times)
    ta = add_to_times(times, ta, 'core aggregation loop')

    set_task_phase('insert classifications')
    remote_insert_classifications(position_values, config)

    ta = add_to_times(times, ta, 'insert classifications')

    record_decode_cache(decode_cache)

    if config.verbosity >= 4:
        config.logger.info(f'End of para_classify {len(package)=} {para=}')

    return total_nr_of_structs, times, stats_list


class Task_tracker:
    """Dispatch/return bookkeeping for the para_classify wait loop.

    The wait loop only ever knew how many tasks were in flight, which says that
    occupancy fell but not why. A ray task holds its cpu from dispatch to return, so
    the question "do a few packages run for the whole phase while the rest churn" is
    answered by the lifetime of each task, not by anything per protein: proteins are
    packed by cost until they reach structs_thresh, and a protein whose own cost is
    over the threshold becomes a package of one that nothing can subdivide.

    Cheap enough to leave on from verbosity 3: one dict entry per in-flight task and
    one tuple per finished one, which is a handful of thousands over a whole run.
    """

    def __init__(self, t0: float):
        # Shared with log_parallelism, so every timestamp in the phase has one origin.
        # It starts before the stores are built, so the first dispatch is not at t=0.
        self.t0 = t0
        self.n_dispatched = 0
        self.t_first_dispatch: float | None = None
        # ObjectRef -> (task_nr, t_dispatch, n_proteins, cost, est_bytes, annotations, chains,
        # prior_variable, uncapped_bytes); see dispatched(). Entries are dropped on
        # return, so this never pins an object the wait loop has already consumed.
        self.in_flight: dict[any, tuple[int, float, int, int]] = {}
        self.finished: list[tuple[int, float, float, int, int]] = []
        self.queue_drained_at: float | None = None
        # Estimated bytes of everything dispatched and not yet returned. This is the
        # quantity the dispatcher packs against, and it has to survive across
        # pack_packages calls - the wait loop re-invokes it once per returned task.
        self.in_flight_bytes = 0
        # Why the previous pack_packages call stopped, so a repeated stall is not re-logged.
        self.last_dispatch_reason: str | None = None
        # The package bounds last logged, so pack_packages logs them only when they change.
        self.last_bounds_logged: tuple | None = None

    def dispatched(self, ref, n_proteins: int, cost: int, est_bytes: int = 0,
                   annotations: int = 0, chains: int = 0, prior_variable: int = 0,
                   uncapped_bytes: int | None = None) -> None:
        now = time.time() - self.t0
        if self.t_first_dispatch is None:
            self.t_first_dispatch = now
        # The package's features are kept so its estimate can be recomputed under a model
        # refitted after it went out (in_flight_estimate), and fed back when it returns.
        # uncapped_bytes is the estimate before the memory model's ceiling; it is what the
        # ceiling is later justified by (ClassificationMemoryModel.ceiling).
        self.in_flight[ref] = (self.n_dispatched, now, n_proteins, cost, est_bytes,
                               annotations, chains, prior_variable,
                               est_bytes if uncapped_bytes is None else uncapped_bytes)
        self.in_flight_bytes += est_bytes
        self.n_dispatched += 1

    def returned(self, ref) -> dict | None:
        """Release the ref; returns what it was dispatched with, or None if unknown."""
        entry = self.in_flight.pop(ref, None)
        if entry is None:
            return None
        task_nr, t_dispatch, n_proteins, cost, est_bytes, annotations, chains, _prior, uncapped = entry
        now = time.time() - self.t0
        self.in_flight_bytes = max(0, self.in_flight_bytes - est_bytes)
        self.finished.append((task_nr, t_dispatch, now, n_proteins, cost))
        return {'estimate': est_bytes, 'uncapped_estimate': uncapped, 'annotations': annotations,
                'chains': chains, 'cost': cost, 'n_proteins': n_proteins,
                'elapsed': now - t_dispatch}

    def in_flight_estimate(self, model: 'ClassificationMemoryModel', prior_fixed: int) -> int:
        """Bytes of everything in flight under the current model, not the one it went out under."""
        fixed = model.fixed_bytes(prior_fixed)
        return int(sum(model.cap(fixed + model.variable_bytes(entry[7], entry[5], entry[6], entry[3]))
                       for entry in self.in_flight.values()))

    def note_queue_drained(self, queue_remaining: int) -> None:
        if self.queue_drained_at is None and queue_remaining == 0:
            self.queue_drained_at = time.time() - self.t0

    def oldest_in_flight(self) -> tuple[int, float, int, int] | None:
        if not self.in_flight:
            return None
        return min(self.in_flight.values(), key=lambda entry: entry[1])

    def log_heartbeat(self, config: Config, queue_remaining: int) -> None:
        now = time.time() - self.t0
        oldest = self.oldest_in_flight()
        if oldest is None:
            oldest_str = 'none'
        else:
            task_nr, t_dispatch, n_proteins, cost, est_bytes = oldest[:5]
            oldest_str = (f'#{task_nr} running={now - t_dispatch:.0f}s '
                          f'(dispatched t={t_dispatch:.0f}s, {n_proteins} prot, cost={cost:_}, '
                          f'est={est_bytes / 1024**3:.1f} GiB)')
        config.logger.info(
            f'[para_classify][occupancy] t={now:.0f}s in_flight={len(self.in_flight)}/{config.proc_n} '
            f'queue_remaining={queue_remaining} finished={len(self.finished)} '
            f'in_flight_est={self.in_flight_estimate(_memory_model, classification_fixed_task_bytes(config)) / 1024**3:.0f} GiB '
            f'(memory model: {_memory_model.describe()}) '
            f'oldest={oldest_str}'
        )

    def log_summary(self, config: Config) -> None:
        if not self.finished:
            return
        # First dispatch to last return: the window in which cpus could be busy at
        # all. Measuring from t0 instead would charge the store building before it to
        # the utilisation figure.
        t_last_return = max(ret for _, _, ret, _, _ in self.finished)
        wall = max(t_last_return - self.t_first_dispatch, 1e-9)
        durations = np.array([ret - disp for _, disp, ret, _, _ in self.finished])
        # What the phase actually got out of the cpus it asked for. A number well
        # under 100% with an empty queue is the tail, not a scheduling problem.
        busy = float(durations.sum())
        capacity = max(wall * config.proc_n, 1e-9)

        config.logger.info(
            f'[para_classify][occupancy] {len(self.finished)} tasks over {wall:.0f}s '
            f'(t={self.t_first_dispatch:.0f}s -> {t_last_return:.0f}s), '
            f'cpu utilisation={100.0 * busy / capacity:.1f}% of {config.proc_n} cores'
        )
        config.logger.info(
            f'[para_classify][occupancy] task seconds: '
            f'min={durations.min():.0f} median={np.median(durations):.0f} '
            f'p90={np.percentile(durations, 90):.0f} max={durations.max():.0f}'
        )

        if self.queue_drained_at is not None:
            # Both are on the t0 clock; wall is a duration, so the tail has to be
            # taken between the two timestamps rather than off wall.
            tail = t_last_return - self.queue_drained_at
            still_running = sum(1 for _, disp, ret, _, _ in self.finished
                                if ret > self.queue_drained_at)
            config.logger.info(
                f'[para_classify][occupancy] queue emptied at t={self.queue_drained_at:.0f}s '
                f'({100.0 * (self.queue_drained_at - self.t_first_dispatch) / wall:.0f}% into the phase); '
                f'tail={tail:.0f}s ({100.0 * tail / wall:.0f}% of the phase) with '
                f'{still_running} task(s) still running and nothing left to backfill with'
            )
            if tail > 0:
                # Mean occupancy over the tail, as core-seconds of overlap divided by
                # its length. This is the number the question comes down to: a value
                # near 1 means a single package is holding the phase open on its own,
                # while a value near proc_n means the tail is simply the last full
                # wave finishing together.
                overlap = sum(min(ret, t_last_return) - max(disp, self.queue_drained_at)
                              for _, disp, ret, _, _ in self.finished
                              if ret > self.queue_drained_at)
                config.logger.info(
                    f'[para_classify][occupancy] mean occupancy during the tail: '
                    f'{overlap / tail:.1f} of {config.proc_n} cores'
                )
        else:
            config.logger.info(
                '[para_classify][occupancy] the queue never emptied - occupancy was '
                'limited by dispatch, not by a straggler tail'
            )

        last_return = max(self.finished, key=lambda entry: entry[2])
        task_nr, t_dispatch, t_return, n_proteins, cost = last_return
        config.logger.info(
            f'[para_classify][occupancy] last task to return: #{task_nr} '
            f'dispatched t={t_dispatch:.0f}s returned t={t_return:.0f}s '
            f'({t_return - t_dispatch:.0f}s, {n_proteins} prot, cost={cost:_}) '
            f'- it held a core for {100.0 * (t_return - t_dispatch) / max(wall, 1e-9):.0f}% of the phase'
        )

        for task_nr, t_dispatch, t_return, n_proteins, cost in sorted(
                self.finished, key=lambda entry: entry[2] - entry[1], reverse=True)[:10]:
            config.logger.info(
                f'[para_classify][occupancy]   slowest #{task_nr}: '
                f'{t_return - t_dispatch:.0f}s (t={t_dispatch:.0f}s -> {t_return:.0f}s) '
                f'proteins={n_proteins} cost={cost:_}'
            )


def pack_packages(
        send_packages: int,
        classification_results: list[ray.ObjectRef],
        proteins: protein_package.Proteins,
        config: Config,
        size_sorted: list[tuple[str, int]],
        para: bool,
        config_store: ray.ObjectRef,
        n_procs: int,
        prot_db_id_map: ray.ObjectRef,
        structs_from_db: ray.ObjectRef,
        matrix_stores: list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]],
        backmap_store: ray.ObjectRef[tuple[np.ndarray[int], np.ndarray[int], np.ndarray[int]]],
        current_nr_of_remote_structs: int,
        total_cost: int,
        structure_db_list: ray.ObjectRef[np.ndarray[bytes]],
        compl_struct_db_id_matrix,
        compl_db_id_map,
        known_empty_complexes_store,
        task_tracker: 'Task_tracker | None' = None,
        memory_estimate_of_protein: dict[str, int] | None = None,
        queue_arrays: dict | None = None,
        ):

    if config.verbosity >= 4:
        config.logger.info(f'Call of pack_packages with: {len(size_sorted)=}, {send_packages=}, {n_procs=} {para=} {total_cost:_=}')
        config.logger.info(f'Current CPU load: {psutil.cpu_percent()}')
    if config.verbosity >= 4 and len(size_sorted) == 0:
        config.logger.info('#############\ncalled pack packages with empty task\n##############')

    sequence_storage: ray.ObjectRef[tuple[np.ndarray, np.ndarray]]
    seq_index_dict: dict[int, tuple[int, int]]
    sequence_storage, seq_index_dict = proteins.seq_store

    # Package cost, per-task memory and how many tasks the node can hold, from one place
    # so the wait loop that re-invokes this agrees with it.
    bounds = classification_package_bounds(config, total_cost)
    structs_thresh, per_task_bytes = bounds.structs_thresh, bounds.per_task_bytes
    # The window deliberately no longer caps n_procs. It divides the headroom by one
    # assumed per-task size, and every real package violated that assumption - which is
    # how 255 slots turned into 119 running workers and 136 idle cores. Memory is bounded
    # by the byte budget in the dispatch loop instead, which lets light packages keep the
    # remaining cores busy rather than leaving them for a heavy one that does not fit.

    # Logged when the bounds change, not per call: they are fixed for a chunk, and at one call
    # per poll this line was written 101,702 times in one whole-PDB chunk.
    bounds_key = (bounds.by_stage, bounds.ceiling, bounds.by_memory, structs_thresh, n_procs)
    bounds_changed = task_tracker is None or task_tracker.last_bounds_logged != bounds_key
    if task_tracker is not None:
        task_tracker.last_bounds_logged = bounds_key
    if config.verbosity >= 2 and bounds_changed:
        config.logger.info(
            f'pack_packages thresholds: by_stage={bounds.by_stage:_} '
            f'ceiling={bounds.ceiling:_} '
            f'by_memory={"none" if bounds.by_memory is None else f"{bounds.by_memory:_}"} '
            f'-> {structs_thresh=:_} ({total_cost=:_}, {config.proc_n=}); '
            f'dispatch window {n_procs} of {config.proc_n} '
            f'at {per_task_bytes / 1024**3:.1f} GiB per task'
        )

    current_cost = 0

    total_times = []
    t_total = time.time()

    # Authoritative, so a caller working from a stale window cannot push one more package
    # in on every wait iteration: the loop below only checks after a dispatch.
    if para and send_packages >= n_procs:
        return send_packages, classification_results, size_sorted, total_times, current_nr_of_remote_structs
    package: list[tuple[str, int, int]] = []

    if isinstance(size_sorted, list):
        snap_size_sorted = np.array(size_sorted)
    else:
        snap_size_sorted = size_sorted
    scheduled_prots = []
    total_loop_times = []

    if config.verbosity >= 4:
        config.logger.info(f'Before loop in pack_packages with: {len(snap_size_sorted)=} {snap_size_sorted=}, {structs_thresh=} {send_packages=}, {n_procs=}')

    def dispatch(package_to_send: list[tuple[str, int, int]], package_cost: int,
                 estimate_bytes: int = 0, annotations: int = 0, interface_chains: int = 0,
                 prior_variable: int = 0, uncapped_bytes: int | None = None):
        if config.verbosity >= 4:
            config.logger.info(
                f'Sending package: {len(package_to_send)=} {package_cost=} '
                f'annotations={annotations} interface_chains={interface_chains} '
                f'estimate={estimate_bytes / 1024**3:.2f} GiB, Pending processes: {send_packages}')
        memory_request = classification_memory_request(
            package_cost, config, estimate_bytes=estimate_bytes if estimate_bytes > 0 else None)
        remote = (para_classify_remote_wrapper if memory_request is None
                  else para_classify_remote_wrapper.options(memory=memory_request))
        result_ref = remote.remote(
            package_to_send,
            [
                config_store, prot_db_id_map, structs_from_db,
                matrix_stores, backmap_store, structure_db_list,
                sequence_storage, compl_struct_db_id_matrix, compl_db_id_map,
                known_empty_complexes_store
            ],
            package_cost = package_cost,
            memory_request = memory_request,
            annotations = annotations,
            interface_chains = interface_chains,
            )
        classification_results.append(result_ref)
        if task_tracker is not None:
            task_tracker.dispatched(result_ref, len(package_to_send), package_cost,
                                    estimate_bytes or memory_request or 0,
                                    annotations, interface_chains, prior_variable,
                                    uncapped_bytes)

    # Two-pointer backfill over the cost-sorted queue; see plan_classification_packages.
    budget = classification_memory_budget(config)
    # Every estimate below - in flight and new - is the static prior raised wherever the model
    # refitted on finished tasks says the prior is too low (ClassificationMemoryModel). Tasks
    # dispatched under an older fit are counted at today's.
    _memory_model.quantile = getattr(config, 'classification_estimate_quantile', 0.9)
    _memory_model.ceiling_margin = getattr(config, 'classification_ceiling_margin', 1.5)
    _memory_model.ceiling_min_samples = getattr(config, 'classification_ceiling_min_samples', 20)
    prior_fixed = classification_fixed_task_bytes(config)
    in_flight_bytes = (task_tracker.in_flight_estimate(_memory_model, prior_fixed)
                       if task_tracker is not None else 0)

    def estimate_of(u_ac) -> int:
        if memory_estimate_of_protein is None:
            return 0
        entry = memory_estimate_of_protein.get(u_ac)
        return entry[0] if entry else 0

    def counts_of(u_ac) -> tuple[int, int]:
        if memory_estimate_of_protein is None:
            return (0, 0)
        entry = memory_estimate_of_protein.get(u_ac)
        return (entry[1], entry[2]) if entry else (0, 0)

    # Integer costs and estimates aligned with the queue. classification() builds them once and
    # this function deletes from them in lockstep with the queue, so a call costs what it
    # dispatches. Rebuilding them here meant re-parsing every cost out of the string-typed queue
    # and a dict lookup per protein, on every call, for all ~27k proteins still queued.
    if queue_arrays is not None and len(queue_arrays['costs']) == len(snap_size_sorted):
        costs = queue_arrays['costs']
        prior_estimates = queue_arrays['estimates']
        annotation_counts = queue_arrays['annotations']
        chain_counts = queue_arrays['chains']
    else:
        costs = np.array([int(entry[1]) for entry in snap_size_sorted], dtype=np.int64)
        prior_estimates = np.array([estimate_of(entry[0]) for entry in snap_size_sorted], dtype=np.int64)
        annotation_counts = np.array([counts_of(entry[0])[0] for entry in snap_size_sorted], dtype=np.int64)
        chain_counts = np.array([counts_of(entry[0])[1] for entry in snap_size_sorted], dtype=np.int64)
        queue_arrays = None
    # Vectorised over the whole queue, so a refit reaches every protein still waiting.
    uncapped_estimates = _memory_model.variable_bytes(prior_estimates, annotation_counts, chain_counts, costs)
    fixed_bytes = _memory_model.fixed_bytes(prior_fixed) if memory_estimate_of_protein else 0
    # A protein alone is a task of fixed + its estimate, so that is what the ceiling bounds per
    # protein; the planner bounds whole packages by it as well. The uncapped estimates are kept
    # for the records the ceiling is activated on.
    ceiling = _memory_model.ceiling() if memory_estimate_of_protein else None
    if ceiling is not None:
        estimates = np.minimum(np.asarray(uncapped_estimates, dtype=np.int64), max(0, ceiling - fixed_bytes))
    else:
        estimates = uncapped_estimates
    slots = (n_procs - send_packages) if para else 1
    planned, stop_reason = plan_classification_packages(
        costs, estimates, structs_thresh, budget, in_flight_bytes, slots, fixed_bytes,
        ceiling_bytes=ceiling)
    dispatched_now = 0

    for indices in planned:
        ti = time.time()
        times = []

        package = []
        current_cost = 0
        package_bytes = fixed_bytes
        package_uncapped_bytes = fixed_bytes
        package_prior_variable = 0
        package_annotations = 0
        package_chains = 0
        for prot_index in indices:
            u_ac = snap_size_sorted[prot_index][0]
            protein_obj: protein_package.Protein = proteins[u_ac]
            package.append((u_ac, protein_obj.database_id, seq_index_dict[protein_obj.database_id][0]))
            current_cost += int(costs[prot_index])
            package_bytes += int(estimates[prot_index])
            package_uncapped_bytes += int(uncapped_estimates[prot_index])
            package_prior_variable += int(prior_estimates[prot_index])
            package_annotations += int(annotation_counts[prot_index])
            package_chains += int(chain_counts[prot_index])
            scheduled_prots.append(prot_index)

        current_nr_of_remote_structs += current_cost

        if not para:
            # The serial path hands the one package it built back to the caller.
            break

        if ceiling is not None:
            package_bytes = min(package_bytes, ceiling)
        dispatch(package, current_cost, package_bytes, package_annotations, package_chains,
                 package_prior_variable, package_uncapped_bytes)
        send_packages += 1
        in_flight_bytes += package_bytes
        dispatched_now += 1

        ti = add_to_times(times, ti, label = 'remote call of para_classify')
        ti = add_to_times(times, ti, label = 'del size_sorted[0]')
        total_loop_times = aggregate_times(total_loop_times, times)

        # Backstop that does not depend on the estimate being right.
        if memory_load_percent() >= getattr(config, 'classification_dispatch_memory_percent', 78.0):
            stop_reason = 'cgroup pressure'
            break

    # Only when something went out or the reason for stopping changed: a stalled budget used to
    # repeat the same line on every poll, ~98k times in one chunk.
    last_reason = task_tracker.last_dispatch_reason if task_tracker is not None else None
    if para and config.verbosity >= 3 and (dispatched_now > 0 or stop_reason != last_reason):
        config.logger.info(
            f'[para_classify][dispatch] {send_packages} in flight, '
            f'{in_flight_bytes / 1024**3:.0f} GiB estimated of '
            f'{(budget / 1024**3) if budget else float("nan"):.0f} GiB budget, '
            f'{_memory_model.describe_ceiling()}, '
            f'{dispatched_now} dispatched now, stopped on: {stop_reason}')
    if task_tracker is not None:
        task_tracker.last_dispatch_reason = stop_reason

    total_times.append(total_loop_times)

    # No trailing flush: the loop above dispatches every package it forms before
    # choosing an end for the next one, so nothing is left half-built in para mode.

    t_total = add_to_times(total_times, t_total, label = 'main loop')

    size_sorted = np.delete(snap_size_sorted, scheduled_prots, axis=0)
    if queue_arrays is not None:
        # Same indices as the queue, so the arrays stay aligned with it for the next call.
        for key in ('costs', 'estimates', 'annotations', 'chains'):
            queue_arrays[key] = np.delete(queue_arrays[key], scheduled_prots)

    if para:
        return send_packages, classification_results, size_sorted, total_times, current_nr_of_remote_structs
    else:
        return package, total_times, current_nr_of_remote_structs


# Complexes already known to have no MicroMiner entry. mm_lookup answers that with an
# os.path.isfile - plus up to three isdir/mkdir calls for AlphaFold ids - and the
# [mm_feats] counters showed ~1.9M lookups a run of which ~98% miss, with the same
# complexes re-checked for every protein a worker handles. Once a complex comes back
# empty it stays empty for the run: update_microminer_active defaults to False so
# nothing creates the entry meanwhile, and when it is on mm_lookup has already tried to
# build it and failed.
def _mmdb_header_length(entry_path: str) -> int | None:
    # Every processed entry starts with the same header line for a given parameter set,
    # so a file whose size equals that header holds no residue rows at all.
    try:
        with open(entry_path, 'rb') as entry_file:
            return len(entry_file.readline())
    except OSError:
        return None


def build_known_empty_complexes(config: Config, complex_ids) -> np.ndarray:
    """Complex ids that yield no MicroMiner features, sorted, for the object store.

    mm_lookup answers "are there features" with a stat plus a parse - and three more
    stats and mkdirs for AlphaFold ids - which the [mm_feats] counters showed happening
    ~1.9M times a run. A run can only query the complexes it loaded, so those are probed
    once here instead of repeatedly in every worker. numpy rather than a set so ray hands
    it to the workers zero-copy.

    An entry counts as empty two ways: no file at all, or a file that is exactly its
    header line and so parses to nothing. The second case is the common one - the first
    version of this only checked for a missing file and left the header-only entries,
    about 40% of those present, being re-parsed by every worker.

    Deliberately the *empty* set and not the populated one: a complex missing from this
    array falls through to the normal lookup, so an incomplete or stale array costs time
    and never correctness. The size test is exact rather than a threshold for the same
    reason - one residue row makes a file strictly larger than its header.
    """
    if config.microminer_exe == '' or config.microminer_db == '':
        return np.empty(0, dtype='S1')

    empty_ids: list[str] = []
    n_absent = 0
    n_header_only = 0
    header_length: int | None = None

    for complex_id in complex_ids:
        entry_path = mmdb_entry_path_readonly(
            complex_id, config.microminer_db,
            identity = config.microminer_identity,
            kmer_matching_rate = config.microminer_kmer_matching_rate,
            fragment_length = config.microminer_fragment_length,
            k = config.microminer_k)

        if entry_path is None:
            empty_ids.append(complex_id)
            n_absent += 1
            continue

        try:
            entry_size = os.path.getsize(entry_path)
        except OSError:
            empty_ids.append(complex_id)
            n_absent += 1
            continue

        if header_length is None:
            header_length = _mmdb_header_length(entry_path)

        if header_length is not None and entry_size == header_length:
            empty_ids.append(complex_id)
            n_header_only += 1

    if config.verbosity >= 4:
        config.logger.info(
            f'[mm_feats] probe: {n_absent:_} complexes without an entry, '
            f'{n_header_only:_} with a header-only entry, header={header_length} bytes')

    if len(empty_ids) == 0:
        return np.empty(0, dtype='S1')

    # dtype='S' sizes itself to the longest id, so nothing is truncated into a match it
    # should not have.
    return np.sort(np.array(empty_ids, dtype='S'))


def is_known_empty_complex(known_empty_complexes: np.ndarray, complex_id: str) -> bool:
    if known_empty_complexes.shape[0] == 0:
        return False
    encoded = complex_id.encode()
    position = np.searchsorted(known_empty_complexes, encoded)
    return bool(position < known_empty_complexes.shape[0] and known_empty_complexes[position] == encoded)


_mm_empty_complexes: set[str] = set()
_mm_lookups_attempted: int = 0
_mm_lookups_served_from_cache: int = 0


def mm_lookup_cache_stats() -> tuple[int, int, int]:
    return _mm_lookups_attempted, _mm_lookups_served_from_cache, len(_mm_empty_complexes)


def single_mm_lookup(package: set[str], config: Config, known_empty_complexes: np.ndarray | None = None)-> dict[str, dict[str, dict[str | int, bytes]]]:
    global _mm_lookups_attempted, _mm_lookups_served_from_cache

    complex_feature_dict_dict: dict[str, dict[str, dict[str | int, bytes]]] = {}
    feature_dict: dict[str, dict[str | int, list[int | float | None]]]
    # The raw feature lists are kept as they come out of the lookup. Only the handful
    # of residues that end up in successful_mm_annotations are ever read back, so
    # building a Microminer_features object and packing it here - for every residue of
    # every chain of every complex fetched - was work thrown away for all but a few.
    # The dict shape and keys are unchanged, the membership tests still work.
    for complex_id in package:
        _mm_lookups_attempted += 1

        if complex_id in _mm_empty_complexes or (
                known_empty_complexes is not None
                and is_known_empty_complex(known_empty_complexes, complex_id)):
            _mm_lookups_served_from_cache += 1
            complex_feature_dict_dict[complex_id] = {}
            continue

        feature_dict, _ = mm_lookup(complex_id, config, n_sub_procs = 1)

        if len(feature_dict) == 0:
            _mm_empty_complexes.add(complex_id)

        complex_feature_dict_dict[complex_id] = {}
        for chain in feature_dict:
            if chain not in complex_feature_dict_dict[complex_id]:
                complex_feature_dict_dict[complex_id][chain] = {}
            for res_id in feature_dict[chain]:
                if res_id not in complex_feature_dict_dict[complex_id][chain]:
                    complex_feature_dict_dict[complex_id][chain][res_id] = feature_dict[chain][res_id]
    return complex_feature_dict_dict


def pack_mm_features(raw_feature_values) -> bytes:
    microminer_features = mappings_package.Microminer_features()
    microminer_features.set_values(raw_feature_values)
    return pack(microminer_features)

class GcTimer:
    """Garbage collections and the time they took while started, from gc.callbacks."""

    def __init__(self):
        self.counts = [0, 0, 0]
        self.seconds = [0.0, 0.0, 0.0]
        self._started_at: float | None = None

    def _callback(self, phase: str, info: dict) -> None:
        if phase == 'start':
            self._started_at = time.perf_counter()
        elif self._started_at is not None:
            generation = min(2, max(0, int(info.get('generation', 2))))
            self.counts[generation] += 1
            self.seconds[generation] += time.perf_counter() - self._started_at
            self._started_at = None

    def start(self) -> 'GcTimer':
        gc.callbacks.append(self._callback)
        return self

    def stop(self) -> None:
        try:
            gc.callbacks.remove(self._callback)
        except ValueError:
            pass

    def report(self) -> str:
        return ' '.join(f'gen{generation}={self.counts[generation]}x/{self.seconds[generation]:.1f}s'
                        for generation in range(3))


def classification(
        proteins: protein_package.Proteins,
        config: Config,
        background_insert_residues_process,
        custom_structure_annotations=None):

    t_class = time.time()
    class_times = []

    classification_results: list[ray.ObjectRef] = []

    size_sorted = []

    if custom_structure_annotations is None:
        protein_ids = proteins.get_protein_ids()
    else:
        protein_ids = custom_structure_annotations.keys()

    if len(protein_ids) == 0:
        config.errorlog.add_warning(f'No proteins to process ad start of classifiction(): {custom_structure_annotations=}')
        return

    total_cost = 0

    prot_db_id_map: list[tuple[int, int, int]] = []
    concat_pos_db_ids = []
    right = 0
    seq_index_dict: dict[int, tuple[int, int]]
    _ , seq_index_dict = proteins.seq_store

    for protein_id in protein_ids:
        position_list: list[int] = proteins[protein_id].get_position_ids()
        prot_db_id: int = proteins[protein_id].database_id
        prot_db_id_map.append((prot_db_id, right, right + len(position_list)))
        right += len(position_list)

        pos_db_id_list = []
        for pos in proteins[protein_id].get_position_ids():
            pos_db_id_list.append(proteins[protein_id].positions[pos].database_id)
        pos_db_id_list = np.array(pos_db_id_list)
        concat_pos_db_ids.append(pos_db_id_list)

        if proteins.protein_map[protein_id].stored:
            continue

        # Counted without building the annotation objects; see Protein.structure_annotations.
        size = proteins[protein_id].annotation_count()
        if size == 0:
            continue

        # Since structure_annotations became flat - one entry per (pdb_id, chain),
        # nowadays a plain list - this counts annotations rather than distinct pdb ids. That is what para_classify
        # actually iterates - its inner loop runs n_positions x n_annotations times,
        # which is exactly what the per-protein stats report - so the cost estimate
        # driving the longest-processing-time packaging is now the right quantity.
        # Packaging therefore differs from runs made before this change.
        seq_len = seq_index_dict[prot_db_id][1]
        cost = size * seq_len
        size_sorted.append((protein_id, cost))

        total_cost += cost

    concat_pos_db_ids = np.concatenate(concat_pos_db_ids)
    prot_db_id_map = np.array(prot_db_id_map)

    t_class = add_to_times(class_times, t_class, label = f'build size_sorted proteins map: {total_cost=:_}')

    if config.verbosity >= 5:
        config.logger.info(f'{prot_db_id_map=}')
        if len(proteins.structures.keys()) <= 1000:
            config.logger.info(f'{proteins.structures.keys()=}')
        else:
            config.logger.info(f'Size of protein.structures: {len(proteins.structures)}')

    size_sorted: list[tuple[str, int]] = sorted(size_sorted, reverse=True, key=lambda x: x[1])

    # Always the ray path. The `else` branch below cannot run: every variable it hands
    # to pack_packages - compl_struct_db_id_matrix, compl_db_id_map,
    # known_empty_complexes_store - is built inside the `if para:` block above it, so it
    # raises UnboundLocalError on the first one. Reviving it means hoisting that whole
    # preparation section out, not flipping this flag.
    #
    # Nothing is lost by that: read_only never reaches classification at all (core skips
    # template selection, alignment and structural analysis, and paraAnnotate is what
    # calls this), so the only caller that would have wanted a serial path is proc_n == 1,
    # which keeps the behaviour it always had.
    para = True

    send_packages = 0

    para_t0 = time.time()
    prev_in_flight = None

    # Per-task lifetimes, so the wait loop can say whether falling occupancy is a
    # straggler tail or a dispatch problem. See Task_tracker.
    # Always present in para mode now: it carries in_flight_bytes, which the dispatcher
    # needs to pack against across pack_packages calls, not just the occupancy logging.
    task_tracker = Task_tracker(para_t0) if para else None
    heartbeat_interval = 60.0 if config.verbosity == 3 else 15.0
    next_heartbeat = heartbeat_interval

    # Backstop for the wait loop below, independent of whether the per-task
    # reservation is calibrated: hold dispatch when the cgroup is this full. Well under
    # ray's own 95% kill threshold, so the tasks already in flight still have room to
    # grow into.
    dispatch_memory_percent = getattr(config, 'classification_dispatch_memory_percent', 78.0)
    memory_hold_logged = False
    # Set when a pack_packages call dispatched nothing. Only a returning task frees memory, so
    # until one returns another call can only reach the same answer. Without this the wait
    # loop, which polls every 0.1 s, asked 98,082 times on a whole-PDB chunk and spent 4.7 h
    # of driver time in pack_packages - time in which returned tasks were not replaced.
    dispatch_stalled = False
    # The window pack_packages derives from the measured headroom; it caps itself, and
    # the wait loop has to agree or it would call it once per iteration for nothing.
    dispatch_window = config.proc_n

    def log_parallelism(label):
        nonlocal prev_in_flight
        if config.verbosity >= 4 and send_packages != prev_in_flight:
            config.logger.info(
                f'[para_classify][parallelism] t={time.time()-para_t0:.1f}s in_flight={send_packages} '
                f'queue_remaining={len(size_sorted)} event={label}'
            )
            prev_in_flight = send_packages

    if config.verbosity >= 3:
        config.logger.info(f'Begin of aggregation process, {para=} {len(size_sorted)=}')

    # Filled in the para branch below, where the annotation counts it is built from are
    # walked. Bound here so the wait loop cannot read it unbound if that ever changes.
    memory_estimate_of_protein: dict[str, tuple[int, int, int]] | None = None
    queue_arrays: dict | None = None

    if para:
        # proc_n == 1 skips the eager ray_init in core, and read_only never gets here at
        # all. Without this, the ray.put below falls into ray's auto_init_hook and starts
        # an unconfigured cluster; ray_init returns immediately if ray is already up.
        ray_init(config, stop_existing = False)
        config_store = ray.put(config)
        # The setup below runs on the driver alone, with every core idle; this says how much of
        # it is garbage collection.
        setup_gc = GcTimer().start()
        prot_db_id_map_store = ray.put((prot_db_id_map, concat_pos_db_ids))

        # Probed once here rather than in every worker; see build_known_empty_complexes.
        known_empty_complexes = build_known_empty_complexes(config, proteins.complexes.keys())
        known_empty_complexes_store = ray.put(known_empty_complexes)

        if config.verbosity >= 4:
            config.logger.info(
                f'[mm_feats] probed {len(proteins.complexes):_} complexes, '
                f'{known_empty_complexes.shape[0]:_} have no MicroMiner entry')

        t_class = add_to_times(class_times, t_class, label = 'prepare db_id store')

        structure_store_db_dict: dict[int, tuple[str, str]] = {}
        structure_id_store_indeces: list[dict[int, int]] = [{} for _ in range(20)]

        # Only the ids actually written are kept, instead of a python list with
        # one pointer per possible db id: that list is ~100x sparser than the
        # number of structures, because an internal id is 100 * a sharded row id.
        # The array below is materialised from this once the widest id is known.
        structure_db_list_entries: dict[int, str] = {}
        # Both of these are addressed by db id rather than packed, so they are
        # mostly -1 padding and are sized by the largest id in the run. They hold
        # complex db ids and internal structure db ids (100 * sharded row id),
        # orders of magnitude below int32's range, and both are ray.put and so
        # replicated per node.
        compl_db_id_map: np.ndarray[np.int32] = np.full(proteins.max_internal_struct_db_id+100, -1, dtype=np.int32)
        
        if len(proteins.complexes) > 0:
            max_number_of_chains = max([len(proteins.complexes[structure_id].chains) for structure_id in proteins.complexes])
            try:
                max_compl_db_id = max([proteins.complexes[structure_id].database_id for structure_id in proteins.complexes if proteins.complexes[structure_id].database_id is not None])
            except ValueError:
                max_compl_db_id = 0
        else:
            max_number_of_chains = 0
            max_compl_db_id = 0

        compl_struct_db_id_matrix = np.full((max_compl_db_id+1, max_number_of_chains), -1, dtype=np.int32)

        t_class = add_to_times(class_times, t_class, label = f'prepare structure stores {max_compl_db_id=} {max_number_of_chains=} {proteins.max_internal_struct_db_id=}')

        annotation_counter = 0
        chain_loop_counter = 0

        stored_structures: set[int] = set()

        # The chain loop below produces two kinds of output: complex-level data that
        # is identical for every protein annotating that complex, and the per-db_nr
        # structure id sets that getBackmaps shards its Alignment_{db_nr} queries by.
        # Deduplicating per (db_nr, complex) instead of per (protein, complex) keeps
        # both correct, and the complex-level writes only have to happen on the very
        # first visit to a complex. compl_seen_per_db_nr and chain_ids_of_complex below.

        # What each protein's task will have to decode, for the memory estimate the
        # dispatcher sizes packages by. This cannot be read off the bookkeeping below:
        # compl_seen is deduplicated per db_nr, so a protein whose complexes another
        # protein of the same db_nr already visited skips the chain loop entirely and
        # would be counted as having no interfaces at all. Counted per protein here, the
        # same way interface_structure_ids counts it in the worker.
        interface_chains_of_protein: dict[str, int] = {}
        annotations_of_protein: dict[str, int] = {}

        # A protein's hits are indices into annotated_structures (Protein.structure_idx),
        # so everything this loop needs per hit - the structure's db id and its complex's
        # db id - is a per-structure column, and the per-hit python pass it used to make
        # (~100M hits on a whole-PDB chunk) is a handful of numpy operations per protein.
        # Only the complexes run python code, and each once per protein at most.
        annotated_structures = proteins.annotated_structures
        annotated_structure_objs = proteins.annotated_structure_objects()
        structure_compl_db_id = np.fromiter(
            (-1 if compl_db_id is None else compl_db_id
             for compl_db_id in (proteins.complexes[pdb_id].database_id for pdb_id, _ in annotated_structures)),
            dtype=np.int64, count=len(annotated_structures))
        used_structures = np.zeros(len(annotated_structures), dtype=bool)

        # Per complex db id: interface chain count (-1 until computed), and whether a
        # protein of each db_nr, or of any, has visited it yet. max_compl_db_id bounds
        # every complex db id, since compl_struct_db_id_matrix above is indexed by them.
        chains_of_complex = np.full(max_compl_db_id + 1, -1, dtype=np.int64)
        compl_seen_per_db_nr = [np.zeros(max_compl_db_id + 1, dtype=bool) for _ in range(20)]
        empty_structure_idx = np.empty(0, dtype=np.int32)
        # Per complex db id, filled on its first visit: the database ids of its stored protein
        # chains, in compl_obj.chains order - what compl_struct_db_id_matrix holds for it. A
        # complex is then walked once, rather than once for each of up to 20 db_nr that reach
        # it: 20.9M chain iterations on a whole-PDB chunk.
        chain_ids_of_complex: dict[int, list[int]] = {}

        # This loop took 943s on a whole-PDB chunk with every core idle, far more than its
        # visible work, and it creates millions of small objects next to the driver's huge
        # Proteins graph - collections over that graph are the suspect. Nothing it builds
        # forms a reference cycle, so refcounting still frees all of it; the same guard as in
        # map_hits_to_structures. The timers say where the time goes if that was not it.
        part4_seconds = {'protein columns': 0.0, 'first visits': 0.0, 'db_nr updates': 0.0}
        gc_was_enabled = gc.isenabled()
        gc.disable()
        try:
            for protein_id, _ in size_sorted:
                t_part = time.perf_counter()
                protein = proteins[protein_id]
                # Index-backed hits only: size_sorted holds unstored proteins, and the only
                # annotation objects added by other means are the stored ones
                # database.getAlignments creates.
                structure_idx = protein.structure_idx if protein.structure_idx is not None else empty_structure_idx
                prot_db_id = protein.database_id
                db_nr = prot_db_id%20

                compl_seen = compl_seen_per_db_nr[db_nr]
                store_indeces = structure_id_store_indeces[db_nr]

                annotations_of_protein[protein_id] = protein.annotation_count()
                annotation_counter += len(structure_idx)
                used_structures[structure_idx] = True

                if config.verbosity >= 6:
                    config.logger.info(
                        f'{protein_id=} ' +
                        str([annotated_structures[index] for index in structure_idx.tolist()]))

                # The protein's complexes in the order its hits first reach them, the order
                # the per-hit loop visited them in, which store_indeces and c_index follow.
                # Hits of structures whose complex has no db id are skipped, as before.
                hit_compl_db_ids = structure_compl_db_id[structure_idx]
                hit_positions = np.flatnonzero(hit_compl_db_ids >= 0)
                hit_compl_db_ids = hit_compl_db_ids[hit_positions]
                _, first_hits = np.unique(hit_compl_db_ids, return_index=True)
                first_hits.sort()
                protein_complexes = hit_compl_db_ids[first_hits]
                # For each complex, the annotated_structures index of the hit that reached it
                # first; its pdb id names the complex.
                first_structures = structure_idx[hit_positions[first_hits]]

                # Counted before the per-db_nr dedup below, which would skip complexes
                # another protein already registered and leave this protein looking free of
                # interfaces.
                uncounted = np.flatnonzero(chains_of_complex[protein_complexes] < 0)
                for position in uncounted.tolist():
                    compl_db_id = int(protein_complexes[position])
                    structure_id = annotated_structures[int(first_structures[position])][0]
                    # The chains compl_struct_db_id_matrix ends up holding for this complex,
                    # which is what interface_structure_ids walks.
                    chains_of_complex[compl_db_id] = proteins.protein_chain_count(structure_id)
                interface_chains = int(chains_of_complex[protein_complexes].sum())

                new_for_db_nr = np.flatnonzero(~compl_seen[protein_complexes])
                compl_seen[protein_complexes] = True
                t_now = time.perf_counter()
                part4_seconds['protein columns'] += t_now - t_part
                t_part = t_now

                for position in new_for_db_nr.tolist():
                    compl_db_id = int(protein_complexes[position])
                    chain_ids = chain_ids_of_complex.get(compl_db_id)
                    if chain_ids is None:
                        t_first = time.perf_counter()
                        structure_id = annotated_structures[int(first_structures[position])][0]
                        compl_obj = proteins.complexes[structure_id]
                        struct_chains = proteins.structures[structure_id]

                        if config.verbosity >= 6:
                            config.logger.info(f'{protein_id=} {structure_id=} {compl_db_id=} {compl_obj.chains=}')

                        chain_ids = []
                        for c_chain, c_chain_type in compl_obj.chains.items():
                            chain_loop_counter += 1
                            if c_chain_type != 'Protein':
                                continue
                            c_struct = struct_chains.get(c_chain)
                            if c_struct is None:
                                continue
                            c_structure_database_id = c_struct.database_id
                            if c_structure_database_id is None:
                                continue

                            structure_store_db_dict[c_structure_database_id] = (structure_id, c_chain)
                            compl_struct_db_id_matrix[compl_db_id][len(chain_ids)] = c_structure_database_id
                            structure_db_list_entries[c_structure_database_id] = f'{structure_id}{c_chain}'
                            compl_db_id_map[c_structure_database_id] = compl_db_id
                            chain_ids.append(c_structure_database_id)

                        chain_ids_of_complex[compl_db_id] = chain_ids
                        part4_seconds['first visits'] += time.perf_counter() - t_first

                    # Only the keys are read (getBackmaps), in insertion order.
                    store_indeces.update(dict.fromkeys(chain_ids))

                interface_chains_of_protein[protein_id] = interface_chains
                part4_seconds['db_nr updates'] += time.perf_counter() - t_part
        finally:
            if gc_was_enabled:
                gc.enable()
        # The first-visit work is inside the per-db_nr loop, so it is taken out of that here.
        part4_seconds['db_nr updates'] -= part4_seconds['first visits']
        del chain_ids_of_complex

        # The hits' own structures, once each instead of once per hit. The order the
        # keys arrive in does not matter: getStoredResidues bins and sorts them.
        for index in np.flatnonzero(used_structures).tolist():
            struct_db_id = annotated_structure_objs[index].database_id
            structure_store_db_dict[struct_db_id] = annotated_structures[index]
            stored_structures.add(struct_db_id)
        del annotated_structure_objs, structure_compl_db_id, used_structures
        del chains_of_complex, compl_seen_per_db_nr

        t_class = add_to_times(class_times, t_class, label = (
            f'build structure_store_db_dict {len(structure_store_db_dict)=} {annotation_counter=} '
            f'{chain_loop_counter=} ' + ' '.join(f'{name}={seconds:.0f}s' for name, seconds in part4_seconds.items())))

        if config.verbosity >= 3:
            total_chains = sum(interface_chains_of_protein.values())
            worst = max(interface_chains_of_protein.items(), key=lambda item: item[1],
                        default=(None, 0))
            config.logger.info(
                f'[para_classify][estimate] interface chains: total={total_chains:_} '
                f'over {len(interface_chains_of_protein):_} proteins, '
                f'mean={total_chains / max(1, len(interface_chains_of_protein)):.0f}, '
                f'worst={worst[0]} with {worst[1]:_} chains '
                f'({classification_memory_estimate(annotations_of_protein.get(worst[0], 0), worst[1], config) / 1024**3:.1f} GiB estimated)')

        # Per protein: (estimated bytes, annotations, interface chains). Carried to
        # pack_packages, which sizes both the reservation and how much may be in flight
        # from the bytes, and logs the two counts so the coefficients behind them can be
        # fitted against the peaks the workers report.
        # The variable part only: the fixed per-worker cost is charged once per package
        # by the dispatcher, not once per protein in it.
        memory_estimate_of_protein = {
            protein_id: (
                classification_protein_bytes(
                    annotations_of_protein.get(protein_id, 0), chains, config),
                annotations_of_protein.get(protein_id, 0),
                chains,
            )
            for protein_id, chains in interface_chains_of_protein.items()
        }

        # Same array np.array(list_of_str, dtype='S') produced: unwritten slots
        # are b'' and the width is the longest id written.
        fused_id_width = max((len(fused_id) for fused_id in structure_db_list_entries.values()), default=1)
        structure_db_list = np.zeros(proteins.max_internal_struct_db_id + 100, dtype=f'S{fused_id_width}')
        for c_structure_database_id, fused_id in structure_db_list_entries.items():
            structure_db_list[c_structure_database_id] = fused_id
        del structure_db_list_entries

        t_class = add_to_times(class_times, t_class, label = f'numpyfy structure_db_list {len(structure_db_list)=}')

        if config.verbosity >= 6:
            config.logger.info(f'{structure_db_list=}')
            config.logger.info(f'{structure_id_store_indeces=}')

        compl_struct_db_id_matrix = ray.put(compl_struct_db_id_matrix)
        compl_db_id_map = ray.put(compl_db_id_map)
        structure_db_list: ray.ObjectRef[np.ndarray[bytes]] = ray.put(structure_db_list)

        t_class = add_to_times(class_times, t_class, label = 'ray.put structure stores')

        if config.verbosity >= 5:
            config.logger.info(f'{stored_structures=} {structure_store_db_dict=}')

        structs_from_db: ray.ObjectRef[np.ndarray[int]]
        matrix_stores: list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]]
        _, structs_from_db, matrix_stores, res_retrieval_times = getStoredResidues(None, config, custom_ids = structure_store_db_dict)

        class_times.append(res_retrieval_times)

        t_class = add_to_times(class_times, t_class, label = 'getStoredResidues')

        if config.fast_pdb_annotation:
            structure_backmaps = ([], [], [])
        else:
            structure_backmaps, getBackmap_times = getBackmaps(
                    structure_id_store_indeces,
                    config,
                    proteins.max_internal_struct_db_id,
                    proteins.database_ids)

            if config.verbosity >= 7:
                config.logger.info(f'{structure_backmaps[2]=}')

            class_times.append(getBackmap_times)

        t_class = add_to_times(class_times, t_class, label = f' getBackmaps {len(structure_id_store_indeces)=}')
        if config.verbosity >= 4:
            config.logger.info(f'getBackmaps {len(structure_id_store_indeces)=}')

        backmap_store: ray.ObjectRef[tuple[np.ndarray[int], np.ndarray[int], np.ndarray[int]]] = ray.put(structure_backmaps)
        # getBackmaps concatenates the 20 shards into one array in this process, and
        # ray.put has now copied it into the object store, where the workers read it
        # from. The driver's copy is dead from here - but it used to stay bound for the
        # whole of classification, which on a whole-PDB chunk is 60 GiB (70.4M
        # alignments) of driver memory held alongside the store's copy of the same
        # thing, on top of the ~93 GiB its Proteins graph already holds.
        del structure_backmaps

        t_class = add_to_times(class_times, t_class, label = 'put structure backmaps')
        setup_gc.stop()
        if config.verbosity >= 3:
            config.logger.info(f'[para_classify][setup] garbage collection before dispatch: {setup_gc.report()}')
        if config.verbosity >= 4:
            config.logger.info('put structure backmaps')

        if config.verbosity >= 8:
            config.logger.info(f'{summarize_objects(raise_on_missing_output=False)=}')
        if config.verbosity >= 5:
            for name, size in sorted(((name, deep_get_size_of(value)) for name, value in locals().items()), key=lambda x: -x[1])[:20]:
                config.logger.info(f"Before pack_packages: {name:>30}: {sizeof_fmt(size):>8}")


        #start_with = 0
        total_packaging_times = []
        #for group_number, (size_sorted, _) in enumerate(groups[start_with:]):
        current_nr_of_remote_structs: int = 0

        # Built once, in queue order; pack_packages deletes from them in lockstep with
        # size_sorted so every later call reads them instead of rebuilding them.
        queue_arrays = {
            'costs': np.array([int(cost) for _, cost in size_sorted], dtype=np.int64),
            'estimates': np.array(
                [memory_estimate_of_protein[protein_id][0]
                 if memory_estimate_of_protein and protein_id in memory_estimate_of_protein else 0
                 for protein_id, _ in size_sorted],
                dtype=np.int64),
            # The model refitted during the stage needs these to re-estimate what is still queued.
            'annotations': np.array(
                [memory_estimate_of_protein[protein_id][1]
                 if memory_estimate_of_protein and protein_id in memory_estimate_of_protein else 0
                 for protein_id, _ in size_sorted],
                dtype=np.int64),
            'chains': np.array(
                [memory_estimate_of_protein[protein_id][2]
                 if memory_estimate_of_protein and protein_id in memory_estimate_of_protein else 0
                 for protein_id, _ in size_sorted],
                dtype=np.int64),
        }

        (send_packages, classification_results, size_sorted, packaging_times, current_nr_of_remote_structs) = pack_packages(send_packages, classification_results,
                                                                                proteins, config, size_sorted,
                                                                                para, config_store, config.proc_n,
                                                                                prot_db_id_map_store, structs_from_db, matrix_stores,
                                                                                backmap_store, current_nr_of_remote_structs, total_cost,
                                                                                structure_db_list,
                                                                                compl_struct_db_id_matrix, compl_db_id_map,
                                                                                known_empty_complexes_store,
                                                                                task_tracker = task_tracker,
                                                                                memory_estimate_of_protein = memory_estimate_of_protein,
                                                                                queue_arrays = queue_arrays)
        total_packaging_times = aggregate_times(total_packaging_times, packaging_times)
        #if send_packages >= config.proc_n:
        #    break
        #start_with += 1

        log_parallelism('initial dispatch')

        class_times.append(total_packaging_times)

        t_class = add_to_times(class_times, t_class, label = f'para pack_packages {send_packages=} {len(size_sorted)=}')
        
    else:

        (package, packaging_times, current_nr_of_remote_structs) = pack_packages(send_packages, classification_results,
                                                                                    proteins, config, size_sorted,
                                                                                    False, None, config.proc_n,
                                                                                    (prot_db_id_map, concat_pos_db_ids), None,
                                                                                    None, 0, total_cost,
                                                                                    compl_struct_db_id_matrix, compl_db_id_map,
                                                                                    known_empty_complexes_store)
        
        class_times.append(packaging_times)
        t_class = add_to_times(class_times, t_class, label = 'serialized pack_packages')
                               
        # Unreachable - see the comment on `para` above. Left in the shape it would
        # need (the getBackmaps call here was still on the old two argument signature)
        # rather than deleted, but it is not exercised and the variables it reads are
        # not in scope on this path.
        if config.fast_pdb_annotation:
            structure_backmaps = ([], [], [])
        else:
            structure_backmaps, getBackmap_times = getBackmaps(
                    structure_id_store_indeces,
                    config,
                    proteins.max_internal_struct_db_id,
                    proteins.database_ids)

            class_times.append(getBackmap_times)

        total_nr_of_structs, agg_times = para_classify(config, package, (prot_db_id_map, concat_pos_db_ids), backmap_store=structure_backmaps)

        class_times.append(agg_times)
        
        t_class = add_to_times(class_times, t_class, label = f'Serialized aggregation {total_nr_of_structs=}')

    if para:
        total_loop_times = []

        loop_count = 0

        stats_list = []

        while True:
            ti = time.time()
            loop_times = []

            loop_count += 1
            if loop_count >= 5000 and config.verbosity >= 3:
                config.logger.info(f'5k waits {len(classification_results)=}')
                loop_count = 0
            ready, not_ready = ray.wait(classification_results, timeout = 0.1)

            ret_times = []
            if len(ready) > 0:
                # A return frees memory and a slot, so dispatching may succeed again.
                dispatch_stalled = False

                remote_returns = ray.get(ready)

                # ray.get keeps the order of 'ready', so record i belongs to return i.
                returned_records = []
                if task_tracker is not None:
                    for ready_ref in ready:
                        returned_records.append(task_tracker.returned(ready_ref))

                for return_index, returned_package in enumerate(remote_returns):
                    send_packages -= 1
                    returned_nr_of_structs, agg_times, proc_stats_list = returned_package[:3]
                    record = returned_records[return_index] if return_index < len(returned_records) else None
                    if len(returned_package) > 3 and record is not None:
                        peak_uss = returned_package[3]
                        record_memory_sample(config, record, peak_uss)
                        if _memory_model.observe(record['annotations'], record['chains'],
                                                 record['cost'], peak_uss,
                                                 uncapped_estimate=record.get('uncapped_estimate', 0)) and config.verbosity >= 2:
                            config.logger.info(f'[para_classify][memory model] {_memory_model.describe()}')
                    stats_list += proc_stats_list
                    current_nr_of_remote_structs -= returned_nr_of_structs
                    ret_times = aggregate_times(ret_times, agg_times)
                
                classification_results = not_ready
                loop_times.append(ret_times)

            else:
                classification_results = not_ready
                loop_times.append([])

            ti = add_to_times(loop_times, ti, label = 'wait')

            log_parallelism('task completed')

            if task_tracker is not None:
                # The queue draining is the moment after which occupancy can only
                # fall: nothing is left to backfill a returning task with.
                task_tracker.note_queue_drained(len(size_sorted))
                if (time.time() - para_t0) >= next_heartbeat:
                    task_tracker.log_heartbeat(config, len(size_sorted))
                    next_heartbeat += heartbeat_interval

            # Ray admits tasks against its memory resource, and a reservation that is
            # too small oversubscribes the node: the OOM killer then takes workers and
            # eventually the raylet, losing the whole stage. This is the backstop that
            # does not depend on the reservation being right - the cgroup's actual load,
            # which is what the killer watches. Dispatching is held off above the
            # threshold and resumes as tasks return and memory frees.
            #
            # Never held off with nothing in flight: that would be a deadlock, since
            # only a running task can free anything.
            # Every core may hold a task; how much memory those tasks may hold between
            # them is pack_packages' byte budget, not a count. Capping the count here as
            # well would idle cores that light packages could have used.
            dispatch_window = config.proc_n

            dispatch_blocked = False
            if send_packages > 0:
                current_memory_load = memory_load_percent()
                if current_memory_load >= dispatch_memory_percent:
                    dispatch_blocked = True
                    if not memory_hold_logged and config.verbosity >= 1:
                        memory_hold_logged = True
                        config.logger.info(
                            f'[para_classify] holding dispatch at {current_memory_load:.1f}% '
                            f'of the memory limit (>= {dispatch_memory_percent:.1f}%), '
                            f'{send_packages} tasks in flight'
                        )
                elif memory_hold_logged and current_memory_load < dispatch_memory_percent - 5.0:
                    memory_hold_logged = False
                    if config.verbosity >= 1:
                        config.logger.info(
                            f'[para_classify] resuming dispatch at {current_memory_load:.1f}% '
                            f'of the memory limit'
                        )

            if (len(size_sorted) > 0 and send_packages < dispatch_window
                    and not dispatch_blocked and not dispatch_stalled):
                sent_before = send_packages
                ret_times = []
                #for group_number, (size_sorted, _) in enumerate(groups[start_with:]):
                #if len(size_sorted) > 0:
                    #group_number += sw
                (send_packages, classification_results, size_sorted, packaging_times, current_nr_of_remote_structs) = pack_packages(send_packages, classification_results, proteins,
                                                                                    config, size_sorted, para, config_store,
                                                                                    config.proc_n, prot_db_id_map_store, structs_from_db, matrix_stores,
                                                                                    backmap_store, current_nr_of_remote_structs, total_cost,
                                                                                    structure_db_list,
                                                                                    compl_struct_db_id_matrix, compl_db_id_map,
                                                                                    known_empty_complexes_store,
                                                                                    task_tracker = task_tracker,
                                                                                    memory_estimate_of_protein = memory_estimate_of_protein,
                                                                                queue_arrays = queue_arrays)
                ret_times = aggregate_times(ret_times, packaging_times)
                dispatch_stalled = classification_dispatch_stalled(sent_before, send_packages)
                log_parallelism('re-dispatch')
                #if config.verbosity >= 6:
                #    config.logger.info(f'Packaged new aggregation process {group_number=} {send_packages=} {len(size_sorted)=}')
                    #if send_packages >= config.proc_n:
                    #    break
                #start_with += 1
                loop_times.append(ret_times)
            else:
                loop_times.append([])
                

            ti = add_to_times(loop_times, ti, label = 'in loop pack_packages')
            total_loop_times = aggregate_times(total_loop_times, loop_times)

            if len(classification_results) == 0:
                break
        if config.verbosity >= 4:
            config.logger.info(f'[para_classify][parallelism] t={time.time()-para_t0:.1f}s all tasks finished, event=drained')
        class_times.append(total_loop_times)

        if task_tracker is not None:
            task_tracker.log_summary(config)

        if config.verbosity >= 3 and len(stats_list) > 0:
            # Aggregate shape of the hot loop. The timing labels tell us where seconds
            # go; these tell us how many iterations produced them, which is what decides
            # whether a further optimisation is worth anything.
            skip_labels = (
                'alignment KeyError', 'pos out of sub_infos', 'alignment gap',
                'no seq_id', 'no store_index', '(unused 5)',
                '(unused 6)', '(unused 7)', 'res_info None', 'no res_aa',
                )
            n_prots = len(stats_list)
            tot_pos = sum(s[3] for s in stats_list)
            tot_annos = sum(s[4] for s in stats_list)
            tot_iters = sum(s[5] for s in stats_list)
            tot_continues = sum(s[6] for s in stats_list)
            tot_mappings = sum(s[7] for s in stats_list)
            tot_counts = [sum(s[10][i] for s in stats_list) for i in range(10)]

            config.logger.info(
                f'[para_classify][loop] proteins={n_prots} positions={tot_pos:_} '
                f'annotations={tot_annos:_} inner_iterations={tot_iters:_} '
                f'skipped={tot_continues:_} ({100.0*tot_continues/max(tot_iters,1):.1f}%) '
                f'mappings_per_position={tot_mappings/max(tot_pos,1):.1f} '
                f'(min={min(s[8] for s in stats_list)} max={max(s[9] for s in stats_list)})'
            )
            for label, count in zip(skip_labels, tot_counts):
                if count > 0:
                    config.logger.info(f'[para_classify][loop]   skip {label:>22}: {count:_} ({100.0*count/max(tot_iters,1):.1f}%)')

            mm_line = mm_feats_summary(stats_list)
            if mm_line is not None:
                config.logger.info(mm_line)

        if config.verbosity >= 4:
            stats_list = sorted(stats_list, key=lambda x:x[2], reverse=True)
            for entry in stats_list[:25]:
                (protein_id, cost, elapsed, n_positions, n_annotations, n_iters,
                 continues, sum_mappings, min_mappings, max_mappings, counts) = entry[:11]
                config.logger.info(
                    f'[para_classify][protein] {protein_id} elapsed={elapsed:.1f}s cost={cost:_} '
                    f'positions={n_positions} annotations={n_annotations} iters={n_iters:_} '
                    f'mappings/pos={sum_mappings/max(n_positions,1):.1f} '
                    f'(min={min_mappings} max={max_mappings}) skipped={continues:_} {counts=}'
                )

    else:
        gc.collect()
        class_times.append([])

    t_class = add_to_times(class_times, t_class, label = 'wait loop')

    del config_store
    del prot_db_id_map_store
    del structs_from_db
    del matrix_stores
    del backmap_store
    del structure_db_list
    del compl_struct_db_id_matrix
    del compl_db_id_map
    del proteins.seq_store

    if background_insert_residues_process is not None:
        background_insert_residues_process.join()
        background_insert_residues_process.close()
        background_insert_residues_process = None

    t_class = add_to_times(class_times, t_class, label = 'background joins')

    if config.verbosity >= 2:
        print_times(class_times, label = 'classification', logger =config.logger)

def _depth_summary(success_depths: list[int]) -> str:
    if len(success_depths) == 0:
        return 'none'
    ordered = sorted(success_depths)
    def at(fraction: float) -> int:
        return ordered[min(len(ordered) - 1, int(fraction * len(ordered)))]
    return (f'n={len(ordered)},p50={at(0.5)},p90={at(0.9)},p99={at(0.99)},'
            f'max={ordered[-1]},within10={sum(1 for d in ordered if d < 10)}')


# The counters calc_mm_features fills in, in the order para_classify's per-protein stats carry
# them and the driver sums them into the [para_classify][mm_feats] line.
MM_COUNTER_NAMES = ('positions', 'rounds', 'complexes_looked_up', 'served_from_cache',
                    'files_parsed', 'lookup_seconds', 'total_seconds', 'successful',
                    'depth_0', 'depth_1_2', 'depth_3_9', 'depth_10_plus')


def mm_feats_summary(stats_list: list) -> str | None:
    """The [para_classify][mm_feats] line: calc_mm_features' counters summed over a chunk.

    MicroMiner features were 69% of the classification task time on a whole-PDB chunk, and
    their counters only existed per protein at verbosity 4. None when no protein carries them.
    """
    mm_rows = [entry[11] for entry in stats_list if len(entry) > 11]
    if not mm_rows:
        return None
    mm = dict(zip(MM_COUNTER_NAMES, (sum(row[i] for row in mm_rows)
                                     for i in range(len(MM_COUNTER_NAMES)))))
    protein_seconds = sum(entry[2] for entry in stats_list)
    return (
        f'[para_classify][mm_feats] proteins={len(mm_rows)} '
        f'seconds={mm["total_seconds"]:_.0f} ({100.0*mm["total_seconds"]/max(protein_seconds, 1e-9):.1f}% of protein time), '
        f'lookup={mm["lookup_seconds"]:_.0f}s ({100.0*mm["lookup_seconds"]/max(mm["total_seconds"], 1e-9):.1f}%) '
        f'positions={mm["positions"]:_} complexes_looked_up={mm["complexes_looked_up"]:_} '
        f'files_parsed={mm["files_parsed"]:_} served_from_cache={mm["served_from_cache"]:_} '
        f'rounds={mm["rounds"]:_} successful={mm["successful"]:_} '
        f'({100.0*mm["successful"]/max(mm["positions"], 1):.1f}% of positions) '
        f'depths 0={mm["depth_0"]:_} 1-2={mm["depth_1_2"]:_} 3-9={mm["depth_3_9"]:_} '
        f'10+={mm["depth_10_plus"]:_}'
    )


def calc_mm_features(
        config,
        recommended_complexes: list[list[tuple[str, str, int | str]]],
        known_empty_complexes: np.ndarray | None = None,
        counters: dict | None = None,
        )-> tuple[dict[str, dict[str, dict[str | int, bytes]]], dict[int, tuple[str, str, str | int]]]:
    """counters, when given, is filled with MM_COUNTER_NAMES: what this protein's lookups cost
    and found, at any verbosity, for the per-chunk summary."""
    mm_db_dict: dict[str, dict[str, dict[str | int, bytes]]] = {}
    successful_mm_annotations: dict[int, tuple[str, str, str | int]] = {}
    t_total = time.perf_counter()
    if counters is not None:
        for name in MM_COUNTER_NAMES:
            counters[name] = 0
    if config.microminer_exe != '':
        attempted_before = _mm_lookups_attempted
        served_before = _mm_lookups_served_from_cache
        lookup_seconds = 0.0
        mm_inputs: set[str] = set()
        processed_complexes = set()
        rec_complex_pos_dict = {}
        
        for pos, recommendation_order in enumerate(recommended_complexes):
            if recommendation_order is None:
                continue
            if len(recommendation_order) > 0:
                recommended_structure_triple = recommendation_order[0]
                complex_id = recommended_structure_triple[0]
                mm_inputs.add(complex_id)
                processed_complexes.add(complex_id)
                if complex_id not in rec_complex_pos_dict:
                    rec_complex_pos_dict[complex_id] = []
                rec_complex_pos_dict[complex_id].append((pos, 0))

        # Work counters: how many complexes get looked up, over how many expansion
        # rounds, and how many of those lookups actually yield a used annotation.
        n_lookup_rounds = 0
        n_complexes_looked_up = len(mm_inputs)
        # How far down the quality-ordered recommendation list each position had to go
        # before it resolved. This is what decides whether the walk can be capped, and
        # at what depth.
        success_depths: list[int] = []

        t_lookup = time.perf_counter()
        complex_feature_dict_dict: dict[str, dict[str, dict[str | int, bytes]]] = single_mm_lookup(mm_inputs, config, known_empty_complexes)
        lookup_seconds += time.perf_counter() - t_lookup

        while len(mm_inputs) > 0:
            n_lookup_rounds += 1
            new_mm_inputs = set()
            n_empty = 0
            for complex_id in complex_feature_dict_dict:
                feature_dict = complex_feature_dict_dict[complex_id]

                if len(feature_dict) == 0:
                    n_empty += 1
                    for (pos, current_try) in rec_complex_pos_dict[complex_id]:
                        stay_in_loop = True
                        next_try = current_try + 1
                        while len(recommended_complexes[pos]) > next_try and stay_in_loop:
                            recommended_structure_triple = recommended_complexes[pos][next_try]
                            complex_id, chain_id, res = recommended_structure_triple

                            if complex_id in mm_db_dict:
                                if (chain_id, res) in mm_db_dict[complex_id]:
                                    successful_mm_annotations[pos] = complex_id, chain_id, res
                                    success_depths.append(next_try)
                                    stay_in_loop = False
                                else:
                                    next_try += 1
                            elif complex_id not in processed_complexes:
                                new_mm_inputs.add(complex_id)
                                processed_complexes.add(complex_id)
                                if complex_id not in rec_complex_pos_dict:
                                    rec_complex_pos_dict[complex_id] = []
                                rec_complex_pos_dict[complex_id].append((pos, next_try))
                                stay_in_loop = False
                            else:
                                next_try += 1

                else:
                    mm_db_dict[complex_id] = feature_dict
                    for (pos, current_try) in rec_complex_pos_dict[complex_id]:
                        recommended_structure_triple = recommended_complexes[pos][current_try]
                        complex_id, chain_id, res = recommended_structure_triple
                        if chain_id in feature_dict:
                            if res in feature_dict[chain_id]:
                                successful_mm_annotations[pos] = complex_id, chain_id, res
                                success_depths.append(current_try)
                                success = True
                            else:
                                success = False
                        else:
                            success = False
                        if not success:
                            next_try = current_try + 1
                            stay_in_loop = True
                            while len(recommended_complexes[pos]) > next_try and stay_in_loop:
                                recommended_structure_triple = recommended_complexes[pos][next_try]
                                complex_id, chain_id, res = recommended_structure_triple
                                if complex_id in mm_db_dict:
                                    successful_mm_annotations[pos] = complex_id, chain_id, res
                                    success_depths.append(next_try)
                                    stay_in_loop = False
                                elif complex_id not in processed_complexes:
                                    new_mm_inputs.add(complex_id)
                                    processed_complexes.add(complex_id)
                                    if complex_id not in rec_complex_pos_dict:
                                        rec_complex_pos_dict[complex_id] = []
                                    rec_complex_pos_dict[complex_id].append((pos, next_try))
                                    stay_in_loop = False
                                else:
                                    next_try += 1

            if config.verbosity >= 6:
                config.logger.info(f'MicroMiner-lookup succesful for {len(mm_inputs) - n_empty} Complexes.')
            
            mm_inputs = new_mm_inputs
            n_complexes_looked_up += len(mm_inputs)
            t_lookup = time.perf_counter()
            complex_feature_dict_dict = single_mm_lookup(mm_inputs, config, known_empty_complexes)
            lookup_seconds += time.perf_counter() - t_lookup

        if counters is not None:
            served = _mm_lookups_served_from_cache - served_before
            counters['positions'] = sum(1 for order in recommended_complexes if order)
            counters['rounds'] = n_lookup_rounds
            counters['complexes_looked_up'] = n_complexes_looked_up
            counters['served_from_cache'] = served
            counters['files_parsed'] = (_mm_lookups_attempted - attempted_before) - served
            counters['lookup_seconds'] = lookup_seconds
            counters['successful'] = len(successful_mm_annotations)
            for depth in success_depths:
                if depth == 0:
                    counters['depth_0'] += 1
                elif depth <= 2:
                    counters['depth_1_2'] += 1
                elif depth <= 9:
                    counters['depth_3_9'] += 1
                else:
                    counters['depth_10_plus'] += 1

        if config.verbosity >= 4:
            n_residues_held = sum(len(chains) for complex_chains in mm_db_dict.values() for chains in complex_chains.values())
            config.logger.info(
                f'[mm_feats] positions={len(recommended_complexes):_} rounds={n_lookup_rounds} '
                f'complexes_looked_up={n_complexes_looked_up:_} complexes_with_data={len(mm_db_dict):_} '
                f'residues_held={n_residues_held:_} successful_annotations={len(successful_mm_annotations):_} '
                f'lookup_cache={_mm_lookups_served_from_cache:_}/{_mm_lookups_attempted:_} '
                f'known_empty={len(_mm_empty_complexes):_} '
                f'depths={_depth_summary(success_depths)}'
            )

    if counters is not None:
        counters['total_seconds'] = time.perf_counter() - t_total
    return mm_db_dict, successful_mm_annotations