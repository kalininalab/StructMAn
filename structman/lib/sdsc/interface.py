import os
import pickle
import ray
import time
import numpy as np
from numba import njit
from numba.core import types as numba_types
from numba.typed import Dict as NumbaDict
from structman.lib.sdsc.sdsc_utils import doomsday_protocol
from structman.base_utils.base_utils import Slotted_obj, add_to_times, aggregate_times
from structman.base_utils.config_class import Config

_MASK_64 = (1 << 64) - 1
_MASK_63 = (1 << 63) - 1
_HASH_SEED = 0x9E3779B97F4A7C15

# The splitmix64 finalizer constants, as uint64 so that the vectorised form
# below stays in uint64: with a python int on the other side, numpy 2 weak
# promotion would keep the array's (possibly narrower) dtype and numpy 1 value
# based casting would turn uint64 into float64. Both corrupt the digest.
_MIX_MULTIPLIER_1 = np.uint64(0xBF58476D1CE4E5B9)
_MIX_MULTIPLIER_2 = np.uint64(0x94D049BB133111EB)
_MIX_SHIFT_1 = np.uint64(30)
_MIX_SHIFT_2 = np.uint64(27)
_MIX_SHIFT_3 = np.uint64(31)


def _splitmix64(value: int) -> int:
    """Avalanche a 64 bit value, so that neighbouring inputs land far apart.

    Python ints with explicit masking, for the same reason the constants above
    are typed: the scalar path must not depend on numpy promotion rules.
    """
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9 & _MASK_64
    value = (value ^ (value >> 27)) * 0x94D049BB133111EB & _MASK_64
    return value ^ (value >> 31)


def _splitmix64_sum(values, count: int) -> int:
    """Avalanche every value on its own, sum the results mod 2**64."""
    mixed = np.fromiter(values, dtype=np.uint64, count=count)
    mixed ^= mixed >> _MIX_SHIFT_1
    mixed *= _MIX_MULTIPLIER_1
    mixed ^= mixed >> _MIX_SHIFT_2
    mixed *= _MIX_MULTIPLIER_2
    mixed ^= mixed >> _MIX_SHIFT_3
    return int(np.add.reduce(mixed, dtype=np.uint64))


def calc_interface_annotation_score(seq_id, coverage, interactions):
    if isinstance(interactions, np.ndarray):
        total_score = np.nansum(interactions)
    else:
        total_score = 0.
        for res_a in interactions:
            for res_b in interactions[res_a]:
                total_score += interactions[res_a][res_b]
    return total_score * seq_id * coverage

class Position_Position_Interaction(Slotted_obj):
    __slots__ = ['protein_a', 'position_a', 'protein_b',
                 'position_b', 'recommended_complex', 'recommended_interaction',
                 'interaction_score']

    slot_mask = [False, False, False,
                 False, True, True,
                 True]

    def __init__(self, protein_a = None, position_a = None, protein_b = None, position_b = None):
        self.protein_a = protein_a
        self.position_a = position_a
        self.protein_b = protein_b
        self.position_b = position_b
        self.interaction_score = 0.

    def deconstruct(self):
        del self.protein_a
        del self.position_a
        del self.protein_b
        del self.position_b
        doomsday_protocol(self)

    def set_recommended_complex(self, recommended_complex, recommended_interaction, residue_interaction_score):
        self.recommended_complex = recommended_complex
        self.recommended_interaction = recommended_interaction
        self.interaction_score = residue_interaction_score


# --- the aggregation core ----------------------------------------------------
#
# The residue loops of calculate_aggregated_interface_map used to run one numba
# call per interface key, which put a jit dispatch (and a numpy scalar read, and
# a couple of python dict operations) in front of every one of the hundreds of
# millions of key iterations a run performs. Interface_Store holds the whole
# working state of one protein's aggregation in flat arrays and numba typed
# dicts instead, so a full pass over an interface entry's keys - overlap
# counting, position bookkeeping, and pos-pos interaction insertion - is a
# single kernel call.
#
# Everything the loops touch is therefore integer-keyed:
#   * structure ids and chain identifiers are interned into `strings`, so the
#     kernels only ever see indices. Chain identifiers are not guaranteed to be a
#     single character, which is why they are interned rather than ord()'d.
#   * position_interface_map (position -> set of interface numbers) becomes a
#     chained adjacency list: pim_head maps a position onto the head of a linked
#     list in pim_next/pim_iface, with -1 for an empty set. The lists hold a
#     handful of entries each, so membership is a short walk.
#   * an interface's positions and its pos-pos interactions are rows in shared
#     arrays, addressed by (interface, position) and
#     (interface, protein_b, position_a, position_b).
#
# The Aggregated_interface objects the rest of the pipeline consumes are built
# from the store once, at the end of the call.

# UniTuple keys rather than packed integers on purpose: whole_pdb reaches ~13.7M
# positions, which needs 24 bits, so packing several fields into one int64 would
# overflow silently on exactly the dataset that matters.
_POSITION_KEY = numba_types.UniTuple(numba_types.int64, 2)
_INTERACTION_KEY = numba_types.UniTuple(numba_types.int64, 4)

_INITIAL_CAPACITY = 64
# Upper bound on the interaction rows reserved ahead of a single pass. The exact
# bound (keys x interacting proteins) is used when it is smaller; beyond this the
# kernel reports back and the caller grows and resumes, so a pathologically large
# interface entry cannot force a huge allocation up front.
_MAX_INTERACTION_RESERVE = 1 << 19

# Statuses and state slots of _process_entries, the kernel that owns the whole
# interface-entry loop. It hands control back at three points it cannot handle
# itself: when the interaction arrays are full, when an entry fused with more than
# one known interface (superfusion allocates fresh index dicts), and when the
# quality measures the annotation score needs are missing.
_ENTRIES_DONE = 0
_ENTRIES_NEED_CAPACITY = 1
_ENTRIES_NEED_SUPERFUSION = 2
_ENTRIES_MISSING_QUALITY = 3
# Position and pim rows are bounded by keys x passes, not by keys: an entry that
# fuses with several known interfaces walks the same keys once per interface and
# adds a row for each. Rather than reserve that product, the kernel hands back the
# same way it does for interactions.
_ENTRIES_NEED_ROW_CAPACITY = 4

_S_POSITIONS = 0
_S_INTERACTIONS = 1
_S_PIM_NODES = 2
_S_CREATED = 3
_S_SKIPPED = 4
_S_ENTRY = 5
_S_PHASE = 6            # 0 begin entry, 1 fusion candidates, 2 new interface
_S_CANDIDATE = 7
_S_KEY = 8
_S_PROTEIN = 9          # -1: the key's position bookkeeping has not run yet
_S_INTERFACES = 10
_S_BACK_MAPPABLE = 11
_S_N_CANDIDATES = 12
_S_N_FUSED = 13
_S_TARGET = 14          # interface the in-flight pass writes to, -1 if none
_S_ALWAYS_SET = 15
_S_STAMP = 16
_S_FUSION_PASSES = 17
_S_NEW_PASSES = 18
_S_KEY_ITERATIONS = 19
_S_SKIP_3 = 20
_S_KEYS_TOTAL = 21
_S_LEN = 22

_PHASE_BEGIN = 0
_PHASE_FUSION = 1
_PHASE_NEW = 2


@njit(cache=True)
def _superfuse(fused_with, n_interfaces,
               iface_score, iface_complex, iface_chain, iface_ichain, iface_pos_count,
               pos_index, pos_arrays, n_positions,
               int_index, int_arrays, int_score, int_has_rec, n_interactions,
               pim_head, pim_arrays, n_pim_nodes):
    """Merge every interface in `fused_with` into the first one, then renumber.

    Reproduces Aggregated_interface.fuse: of two merged interfaces the
    higher-scoring one keeps the recommendation and the score, and the
    lower-scoring one's rows win wherever both hold the same key. Merges are
    applied in the given (sorted) order, because that chaining is what the old
    code did.
    """
    pos_position, pos_struct, pos_chain, pos_res, pos_iface = pos_arrays
    (int_iface, int_protein_b, int_pos_a, int_pos_b, int_rec_complex,
     int_rec_chain_a, int_rec_res_a, int_rec_chain_b, int_rec_res_b) = int_arrays
    pim_next, pim_iface = pim_arrays

    removed = np.zeros(n_interfaces, np.uint8)
    target = fused_with[0]

    for fusion_index in range(1, fused_with.shape[0]):
        source = fused_with[fusion_index]
        # `fuse` hands the merge to the higher-scoring interface, so its identity
        # survives and the other one's values overwrite on a key collision.
        source_keeps_identity = iface_score[source] > iface_score[target]

        for row in range(n_positions):
            if pos_iface[row] != source:
                continue
            key = (target, pos_position[row])
            if key not in pos_index:
                pos_index[key] = row
                pos_iface[row] = target
                continue
            target_row = pos_index[key]
            if not source_keeps_identity:
                pos_struct[target_row] = pos_struct[row]
                pos_chain[target_row] = pos_chain[row]
                pos_res[target_row] = pos_res[row]
            pos_iface[row] = -1

        for row in range(n_interactions):
            if int_iface[row] != source:
                continue
            key = (target, int_protein_b[row], int_pos_a[row], int_pos_b[row])
            if key not in int_index:
                int_index[key] = row
                int_iface[row] = target
                continue
            target_row = int_index[key]
            if not source_keeps_identity:
                int_score[target_row] = int_score[row]
                int_has_rec[target_row] = int_has_rec[row]
                if int_has_rec[row] != 0:
                    int_rec_complex[target_row] = int_rec_complex[row]
                    int_rec_chain_a[target_row] = int_rec_chain_a[row]
                    int_rec_res_a[target_row] = int_rec_res_a[row]
                    int_rec_chain_b[target_row] = int_rec_chain_b[row]
                    int_rec_res_b[target_row] = int_rec_res_b[row]
            int_iface[row] = -1

        if source_keeps_identity:
            iface_score[target] = iface_score[source]
            iface_complex[target] = iface_complex[source]
            iface_chain[target] = iface_chain[source]
            iface_ichain[target] = iface_ichain[source]

        removed[source] = 1

    # Renumbering. A kept interface shifts down by the number of interfaces
    # removed before it; a removed one is now the merge target, so everything that
    # referred to it has to be redirected there. The old code assigned it
    # `number - diff_count` with diff_count already counting itself, which is the
    # new number of whichever interface happens to sit below it - so the
    # position_interface_map entries of a superfused interface used to end up
    # pointing at an unrelated one, and its positions then counted as overlap for
    # that one in every later fusion decision.
    number_map = np.empty(n_interfaces, np.int64)
    diff_count = 0
    new_n_interfaces = 0
    for interface_number in range(n_interfaces):
        if removed[interface_number] == 0:
            number_map[interface_number] = interface_number - diff_count
            iface_score[new_n_interfaces] = iface_score[interface_number]
            iface_complex[new_n_interfaces] = iface_complex[interface_number]
            iface_chain[new_n_interfaces] = iface_chain[interface_number]
            iface_ichain[new_n_interfaces] = iface_ichain[interface_number]
            new_n_interfaces += 1
        else:
            diff_count += 1
            number_map[interface_number] = -1

    target_number = number_map[target]
    for interface_number in range(n_interfaces):
        if removed[interface_number] != 0:
            number_map[interface_number] = target_number

    new_pos_index = NumbaDict.empty(_POSITION_KEY, numba_types.int64)
    write = 0
    for row in range(n_positions):
        interface_number = pos_iface[row]
        if interface_number == -1:
            continue
        new_number = number_map[interface_number]
        pos_iface[write] = new_number
        pos_position[write] = pos_position[row]
        pos_struct[write] = pos_struct[row]
        pos_chain[write] = pos_chain[row]
        pos_res[write] = pos_res[row]
        new_pos_index[(new_number, pos_position[write])] = write
        write += 1
    new_n_positions = write

    new_int_index = NumbaDict.empty(_INTERACTION_KEY, numba_types.int64)
    write = 0
    for row in range(n_interactions):
        interface_number = int_iface[row]
        if interface_number == -1:
            continue
        new_number = number_map[interface_number]
        int_iface[write] = new_number
        int_protein_b[write] = int_protein_b[row]
        int_pos_a[write] = int_pos_a[row]
        int_pos_b[write] = int_pos_b[row]
        int_score[write] = int_score[row]
        int_rec_complex[write] = int_rec_complex[row]
        int_rec_chain_a[write] = int_rec_chain_a[row]
        int_rec_res_a[write] = int_rec_res_a[row]
        int_rec_chain_b[write] = int_rec_chain_b[row]
        int_rec_res_b[write] = int_rec_res_b[row]
        int_has_rec[write] = int_has_rec[row]
        new_int_index[(new_number, int_protein_b[write],
                       int_pos_a[write], int_pos_b[write])] = write
        write += 1
    new_n_interactions = write

    # The chains are rebuilt rather than patched in place, because renumbering can
    # map two members of one position's set onto the same interface and the old
    # code deduplicated by rebuilding the set.
    positions = np.empty(len(pim_head), np.int64)
    write = 0
    for position in pim_head:
        positions[write] = position
        write += 1

    rebuilt_next = np.empty(n_pim_nodes, np.int64)
    rebuilt_iface = np.empty(n_pim_nodes, np.int64)
    write = 0
    for position_index in range(positions.shape[0]):
        position = positions[position_index]
        node = pim_head[position]
        head = -1
        while node != -1:
            new_number = number_map[pim_iface[node]]
            duplicate = False
            scan = head
            while scan != -1:
                if rebuilt_iface[scan] == new_number:
                    duplicate = True
                    break
                scan = rebuilt_next[scan]
            if not duplicate:
                rebuilt_next[write] = head
                rebuilt_iface[write] = new_number
                head = write
                write += 1
            node = pim_next[node]
        pim_head[position] = head
    new_n_pim_nodes = write
    for node in range(new_n_pim_nodes):
        pim_next[node] = rebuilt_next[node]
        pim_iface[node] = rebuilt_iface[node]

    for interface_number in range(new_n_interfaces):
        iface_pos_count[interface_number] = 0
    for row in range(new_n_positions):
        iface_pos_count[pos_iface[row]] += 1

    return (new_pos_index, new_int_index, new_n_interfaces,
            new_n_positions, new_n_interactions, new_n_pim_nodes)


@njit(cache=True)
def _run_entry_pass(
        concat_keys, concat_scores, key_left, n_keys, backmap,
        bm_protein, bm_left, bm_length, bm_offset, n_proteins, concat_backmaps,
        interface_number, always_set,
        struct_idx, chain_idx, i_chain_idx, compl_db_id, iface_pos_count,
        pos_index, pos_position, pos_struct, pos_chain, pos_res, pos_iface,
        int_index, int_iface, int_protein_b, int_pos_a, int_pos_b,
        int_rec_complex, int_rec_chain_a, int_rec_res_a, int_rec_chain_b, int_rec_res_b,
        int_score, int_has_rec, int_capacity, pos_capacity, pim_capacity,
        pim_head, pim_next, pim_iface, state):
    """One pass over an entry's keys, for one target interface.

    Covers both residue loops of the original: `always_set` is
    `interface_annotation_score > interface.interface_annotation_score` for a
    fusion pass and True for a freshly created interface, which is what decided
    the conditional add_position and the recommended-complex overwrite there.

    Both halves of each key are unpacked here, and the interacting protein's
    position is read straight out of concat_backmaps through that protein's
    (left, length) slice; an i_res past the end gives -1, which is the branch the
    old IndexError took. state[_S_KEY] and state[_S_PROTEIN] are the resume point,
    with a protein index of -1 meaning the key's position bookkeeping has not run.
    """
    n_positions = state[_S_POSITIONS]
    n_interactions = state[_S_INTERACTIONS]
    n_pim_nodes = state[_S_PIM_NODES]
    created = state[_S_CREATED]
    skipped = state[_S_SKIPPED]
    index = state[_S_KEY]
    protein_start = state[_S_PROTEIN]

    backmap_length = backmap.shape[0]

    while index < n_keys:
        key = concat_keys[key_left + index]
        res = key >> 32
        if res >= backmap_length:
            index += 1
            protein_start = -1
            continue
        position = np.int64(backmap[res])
        if position == -1:
            index += 1
            protein_start = -1
            continue

        if protein_start < 0:
            # At most one position row and one pim node come out of this block.
            if n_positions >= pos_capacity or n_pim_nodes >= pim_capacity:
                state[_S_POSITIONS] = n_positions
                state[_S_INTERACTIONS] = n_interactions
                state[_S_PIM_NODES] = n_pim_nodes
                state[_S_CREATED] = created
                state[_S_SKIPPED] = skipped
                state[_S_KEY] = index
                state[_S_PROTEIN] = -1
                return _ENTRIES_NEED_ROW_CAPACITY

            # position_interface_map[position].add(interface_number)
            head = -1
            if position in pim_head:
                head = pim_head[position]
            node = head
            already_known = False
            while node != -1:
                if pim_iface[node] == interface_number:
                    already_known = True
                    break
                node = pim_next[node]
            if not already_known:
                pim_next[n_pim_nodes] = head
                pim_iface[n_pim_nodes] = interface_number
                pim_head[position] = n_pim_nodes
                n_pim_nodes += 1

            pos_key = (interface_number, position)
            if pos_key in pos_index:
                pos_row = pos_index[pos_key]
                write_position = always_set
            else:
                pos_row = n_positions
                pos_index[pos_key] = pos_row
                pos_position[pos_row] = position
                pos_iface[pos_row] = interface_number
                n_positions += 1
                iface_pos_count[interface_number] += 1
                write_position = True
            if write_position:
                pos_struct[pos_row] = struct_idx
                pos_chain[pos_row] = chain_idx
                pos_res[pos_row] = res

            protein_start = 0

        score = np.float64(concat_scores[key_left + index])
        if score == score:
            i_res = key & 0xFFFFFFFF
            for protein_position in range(protein_start, n_proteins):
                slot = bm_offset + protein_position
                if i_res >= bm_length[slot]:
                    continue
                position_b = np.int64(concat_backmaps[bm_left[slot] + i_res])
                if position_b == -1:
                    continue
                protein_b = bm_protein[slot]
                int_key = (interface_number, protein_b, position, position_b)

                known_row = int_key in int_index
                is_new = False
                if not known_row:
                    if n_interactions >= int_capacity:
                        state[_S_POSITIONS] = n_positions
                        state[_S_INTERACTIONS] = n_interactions
                        state[_S_PIM_NODES] = n_pim_nodes
                        state[_S_CREATED] = created
                        state[_S_SKIPPED] = skipped
                        state[_S_KEY] = index
                        state[_S_PROTEIN] = protein_position
                        return _ENTRIES_NEED_CAPACITY
                    row = n_interactions
                    int_index[int_key] = row
                    int_iface[row] = interface_number
                    int_protein_b[row] = protein_b
                    int_pos_a[row] = position
                    int_pos_b[row] = position_b
                    int_score[row] = 0.0
                    int_has_rec[row] = 0
                    n_interactions += 1
                    is_new = True
                else:
                    row = int_index[int_key]
                    if score > int_score[row]:
                        # A beaten row is reset first, exactly as the old code
                        # swapped in a freshly built interaction object.
                        int_score[row] = 0.0
                        int_has_rec[row] = 0
                        is_new = True

                if is_new:
                    created += 1
                else:
                    skipped += 1

                if is_new or always_set:
                    int_score[row] = score
                    int_rec_complex[row] = compl_db_id
                    int_rec_chain_a[row] = chain_idx
                    int_rec_res_a[row] = res
                    int_rec_chain_b[row] = i_chain_idx
                    int_rec_res_b[row] = i_res
                    int_has_rec[row] = 1

        index += 1
        protein_start = -1

    state[_S_POSITIONS] = n_positions
    state[_S_INTERACTIONS] = n_interactions
    state[_S_PIM_NODES] = n_pim_nodes
    state[_S_CREATED] = created
    state[_S_SKIPPED] = skipped
    state[_S_KEY] = index
    state[_S_PROTEIN] = -1
    return _ENTRIES_DONE


@njit(cache=True)
def _process_entries(
        # the entries of one (structure, chain), as columns
        entry_key_left, entry_key_right, entry_i_chain_idx, entry_score,
        entry_bm_offset, entry_bm_count, entry_has_quality,
        bm_protein, bm_left, bm_length,
        concat_keys, concat_scores, concat_backmaps, backmap,
        struct_idx, chain_idx, compl_db_id,
        # scratch that has to survive a hand-back
        scan_overlap, scan_candidates, fused_with,
        seen_stamp,
        # the store
        iface_score, iface_complex, iface_chain, iface_ichain, iface_pos_count,
        pos_index, pos_arrays,
        int_index, int_arrays, int_score, int_has_rec, int_capacity,
        pos_capacity, pim_capacity,
        pim_head, pim_arrays,
        state):
    """The whole per-interface-entry loop of calculate_aggregated_interface_map.

    One call walks every entry of one (structure, chain): the residue scan, the
    fusion decision, every fusion pass, and the new-interface pass. Nothing here
    is per-entry python any more - the entries arrive as columns addressing the
    key and score stores directly.

    `state` is the resume cursor as well as the counters, so the three hand-backs
    (see _ENTRIES_*) return to exactly the key and protein they stopped on. A
    hand-back for superfusion happens at an entry boundary, because the merge
    replaces the index dicts.
    """
    pos_position, pos_struct, pos_chain, pos_res, pos_iface = pos_arrays
    (int_iface, int_protein_b, int_pos_a, int_pos_b, int_rec_complex,
     int_rec_chain_a, int_rec_res_a, int_rec_chain_b, int_rec_res_b) = int_arrays
    pim_next, pim_iface = pim_arrays

    n_entries = entry_key_left.shape[0]
    backmap_length = backmap.shape[0]

    while state[_S_ENTRY] < n_entries:
        entry = state[_S_ENTRY]
        key_left = entry_key_left[entry]
        n_keys = entry_key_right[entry] - key_left
        bm_offset = entry_bm_offset[entry]
        n_proteins = entry_bm_count[entry]

        if state[_S_PHASE] == _PHASE_BEGIN:
            if entry_i_chain_idx[entry] == -1:
                state[_S_SKIP_3] += 1
                state[_S_ENTRY] += 1
                continue

            state[_S_KEYS_TOTAL] += n_keys

            # --- residue scan: back-mappable positions, and the interfaces they
            # already belong to. A residue counts once however many keys carry it,
            # which is what walking set(res_list) did.
            state[_S_STAMP] += 1
            stamp = state[_S_STAMP]
            back_mappable = 0
            n_candidates = 0
            for key_index in range(n_keys):
                res = concat_keys[key_left + key_index] >> 32
                if res >= backmap_length:
                    continue
                if seen_stamp[res] == stamp:
                    continue
                seen_stamp[res] = stamp
                position = np.int64(backmap[res])
                if position == -1:
                    continue
                back_mappable += 1

                if position not in pim_head:
                    pim_head[position] = -1
                    continue
                node = pim_head[position]
                while node != -1:
                    known_interface = pim_iface[node]
                    if scan_overlap[known_interface] == 0:
                        scan_candidates[n_candidates] = known_interface
                        n_candidates += 1
                    scan_overlap[known_interface] += 1
                    node = pim_next[node]

            if back_mappable == 0:
                for candidate in range(n_candidates):
                    scan_overlap[scan_candidates[candidate]] = 0
                state[_S_ENTRY] += 1
                continue

            if entry_has_quality == 0:
                for candidate in range(n_candidates):
                    scan_overlap[scan_candidates[candidate]] = 0
                return _ENTRIES_MISSING_QUALITY

            # Sorted so the fusion passes run in a defined order; they are
            # independent of one another, and fused_with was sorted anyway.
            scan_candidates[:n_candidates].sort()

            state[_S_BACK_MAPPABLE] = back_mappable
            state[_S_N_CANDIDATES] = n_candidates
            state[_S_N_FUSED] = 0
            state[_S_CANDIDATE] = 0
            state[_S_TARGET] = -1
            state[_S_PHASE] = _PHASE_FUSION

        score_of_entry = entry_score[entry]

        if state[_S_PHASE] == _PHASE_FUSION:
            while True:
                if state[_S_TARGET] == -1:
                    if state[_S_CANDIDATE] >= state[_S_N_CANDIDATES]:
                        break
                    known_interface = scan_candidates[state[_S_CANDIDATE]]
                    overlap = scan_overlap[known_interface]
                    fuses = (overlap / state[_S_BACK_MAPPABLE] > 0.5
                             or overlap / iface_pos_count[known_interface] > 0.5)
                    if not fuses:
                        state[_S_CANDIDATE] += 1
                        continue
                    state[_S_TARGET] = known_interface
                    state[_S_ALWAYS_SET] = 1 if score_of_entry > iface_score[known_interface] else 0
                    state[_S_KEY] = 0
                    state[_S_PROTEIN] = -1
                    state[_S_FUSION_PASSES] += 1
                    state[_S_KEY_ITERATIONS] += n_keys

                status = _run_entry_pass(
                    concat_keys, concat_scores, key_left, n_keys, backmap,
                    bm_protein, bm_left, bm_length, bm_offset, n_proteins,
                    concat_backmaps,
                    state[_S_TARGET], state[_S_ALWAYS_SET] == 1,
                    struct_idx, chain_idx, entry_i_chain_idx[entry], compl_db_id,
                    iface_pos_count,
                    pos_index, pos_position, pos_struct, pos_chain, pos_res, pos_iface,
                    int_index, int_iface, int_protein_b, int_pos_a, int_pos_b,
                    int_rec_complex, int_rec_chain_a, int_rec_res_a,
                    int_rec_chain_b, int_rec_res_b, int_score, int_has_rec,
                    int_capacity, pos_capacity, pim_capacity,
                    pim_head, pim_next, pim_iface, state)
                if status != _ENTRIES_DONE:
                    return status

                known_interface = state[_S_TARGET]
                fused_with[state[_S_N_FUSED]] = known_interface
                state[_S_N_FUSED] += 1
                if state[_S_ALWAYS_SET] == 1:
                    iface_complex[known_interface] = struct_idx
                    iface_chain[known_interface] = chain_idx
                    iface_ichain[known_interface] = entry_i_chain_idx[entry]
                    iface_score[known_interface] = score_of_entry
                state[_S_TARGET] = -1
                state[_S_CANDIDATE] += 1

            for candidate in range(state[_S_N_CANDIDATES]):
                scan_overlap[scan_candidates[candidate]] = 0

            if state[_S_N_FUSED] > 1:
                # Superfusion renumbers the interfaces and rebuilds the index
                # dicts, so it is done outside; the entry is finished either way.
                state[_S_PHASE] = _PHASE_BEGIN
                state[_S_ENTRY] += 1
                return _ENTRIES_NEED_SUPERFUSION
            if state[_S_N_FUSED] > 0:
                state[_S_PHASE] = _PHASE_BEGIN
                state[_S_ENTRY] += 1
                continue

            state[_S_TARGET] = -1
            state[_S_PHASE] = _PHASE_NEW

        # --- no fusion: a new interface, over the same keys.
        if state[_S_TARGET] == -1:
            new_interface = state[_S_INTERFACES]
            iface_score[new_interface] = score_of_entry
            iface_complex[new_interface] = struct_idx
            iface_chain[new_interface] = chain_idx
            iface_ichain[new_interface] = entry_i_chain_idx[entry]
            iface_pos_count[new_interface] = 0
            state[_S_INTERFACES] = new_interface + 1
            state[_S_TARGET] = new_interface
            # back_mappable > 0, so the old code's `trigger` always fired and the
            # interface was always created; its arguments did not depend on which
            # key fired it.
            state[_S_ALWAYS_SET] = 1
            state[_S_KEY] = 0
            state[_S_PROTEIN] = -1
            state[_S_NEW_PASSES] += 1
            state[_S_KEY_ITERATIONS] += n_keys

        status = _run_entry_pass(
            concat_keys, concat_scores, key_left, n_keys, backmap,
            bm_protein, bm_left, bm_length, bm_offset, n_proteins, concat_backmaps,
            state[_S_TARGET], True,
            struct_idx, chain_idx, entry_i_chain_idx[entry], compl_db_id,
            iface_pos_count,
            pos_index, pos_position, pos_struct, pos_chain, pos_res, pos_iface,
            int_index, int_iface, int_protein_b, int_pos_a, int_pos_b,
            int_rec_complex, int_rec_chain_a, int_rec_res_a,
            int_rec_chain_b, int_rec_res_b, int_score, int_has_rec,
            int_capacity, pos_capacity, pim_capacity,
                    pim_head, pim_next, pim_iface, state)
        if status != _ENTRIES_DONE:
            return status

        state[_S_TARGET] = -1
        state[_S_PHASE] = _PHASE_BEGIN
        state[_S_ENTRY] += 1

    return _ENTRIES_DONE


def _complex_backmap_columns(structure_backmaps):
    """Every structure's interacting proteins as three flat columns, plus a span map.

    (protein id, backmap start, backmap length) per protein, all structures of the
    complex concatenated once, with span[structure] giving the (offset, count) into
    them. Built per complex rather than per chain, because the same interacting
    structure is reached from every chain of the complex - and rather than the
    (keys x proteins) matrix of positions the caller used to gather, since the
    interface pass reads the positions itself out of concat_backmaps.
    """
    total = sum(len(proteins) for proteins in structure_backmaps.values())
    protein_ids = np.empty(total, dtype=np.int64)
    lefts = np.empty(total, dtype=np.int64)
    lengths = np.empty(total, dtype=np.int64)

    span: dict[int, tuple[int, int]] = {}
    offset = 0
    for structure_db_id, proteins in structure_backmaps.items():
        span[structure_db_id] = (offset, len(proteins))
        for protein_id, (left, right) in proteins.items():
            protein_ids[offset] = protein_id
            lefts[offset] = left
            lengths[offset] = right - left
            offset += 1
    return (protein_ids, lefts, lengths), span


def _grow(array, capacity):
    grown = np.zeros((capacity,) + array.shape[1:], dtype=array.dtype)
    grown[:array.shape[0]] = array
    return grown


class Interface_Store:
    """The whole aggregation state of one protein, in numba-addressable form."""

    def __init__(self):
        self.strings: list[str] = []
        self.string_index: dict[str, int] = {}

        self.n_interfaces = 0
        self.interface_capacity = _INITIAL_CAPACITY
        self.iface_score = np.zeros(_INITIAL_CAPACITY, dtype=np.float64)
        self.iface_complex = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.iface_chain = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.iface_ichain = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.iface_pos_count = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)

        self.n_positions = 0
        self.position_capacity = _INITIAL_CAPACITY
        self.pos_index = NumbaDict.empty(_POSITION_KEY, numba_types.int64)
        self.pos_position = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.pos_struct = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.pos_chain = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.pos_res = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.pos_iface = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)

        self.n_interactions = 0
        self.interaction_capacity = _INITIAL_CAPACITY
        self.int_index = NumbaDict.empty(_INTERACTION_KEY, numba_types.int64)
        self.int_iface = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.int_protein_b = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.int_pos_a = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.int_pos_b = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.int_score = np.zeros(_INITIAL_CAPACITY, dtype=np.float64)
        self.int_rec_complex = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.int_rec_chain_a = np.full(_INITIAL_CAPACITY, -1, dtype=np.int64)
        self.int_rec_res_a = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.int_rec_chain_b = np.full(_INITIAL_CAPACITY, -1, dtype=np.int64)
        self.int_rec_res_b = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.int_has_rec = np.zeros(_INITIAL_CAPACITY, dtype=np.uint8)

        self.n_pim_nodes = 0
        self.pim_capacity = _INITIAL_CAPACITY
        self.pim_head = NumbaDict.empty(numba_types.int64, numba_types.int64)
        self.pim_next = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)
        self.pim_iface = np.zeros(_INITIAL_CAPACITY, dtype=np.int64)

        self.state = np.zeros(7, dtype=np.int64)

        # Residue-deduplication scratch for scan_residues; grown to the longest
        # backmap this protein uses, and stamped rather than cleared.
        self.seen_stamp = np.zeros(0, dtype=np.int64)
        self.scan_stamp = 0

    # --- interning --------------------------------------------------------

    def intern(self, text) -> int:
        number = self.string_index.get(text)
        if number is None:
            number = len(self.strings)
            self.strings.append(text)
            self.string_index[text] = number
        return number

    # --- capacity ---------------------------------------------------------

    @staticmethod
    def _next_capacity(current: int, needed: int) -> int:
        while current < needed:
            current *= 2
        return current

    def _ensure_interfaces(self, needed: int) -> None:
        if needed <= self.interface_capacity:
            return
        capacity = self._next_capacity(self.interface_capacity, needed)
        self.iface_score = _grow(self.iface_score, capacity)
        self.iface_complex = _grow(self.iface_complex, capacity)
        self.iface_chain = _grow(self.iface_chain, capacity)
        self.iface_ichain = _grow(self.iface_ichain, capacity)
        self.iface_pos_count = _grow(self.iface_pos_count, capacity)
        self.interface_capacity = capacity

    def _ensure_positions(self, needed: int) -> None:
        if needed <= self.position_capacity:
            return
        capacity = self._next_capacity(self.position_capacity, needed)
        self.pos_position = _grow(self.pos_position, capacity)
        self.pos_struct = _grow(self.pos_struct, capacity)
        self.pos_chain = _grow(self.pos_chain, capacity)
        self.pos_res = _grow(self.pos_res, capacity)
        self.pos_iface = _grow(self.pos_iface, capacity)
        self.position_capacity = capacity

    def _ensure_interactions(self, needed: int) -> None:
        if needed <= self.interaction_capacity:
            return
        capacity = self._next_capacity(self.interaction_capacity, needed)
        self.int_iface = _grow(self.int_iface, capacity)
        self.int_protein_b = _grow(self.int_protein_b, capacity)
        self.int_pos_a = _grow(self.int_pos_a, capacity)
        self.int_pos_b = _grow(self.int_pos_b, capacity)
        self.int_score = _grow(self.int_score, capacity)
        self.int_rec_complex = _grow(self.int_rec_complex, capacity)
        self.int_rec_chain_a = _grow(self.int_rec_chain_a, capacity)
        self.int_rec_res_a = _grow(self.int_rec_res_a, capacity)
        self.int_rec_chain_b = _grow(self.int_rec_chain_b, capacity)
        self.int_rec_res_b = _grow(self.int_rec_res_b, capacity)
        self.int_has_rec = _grow(self.int_has_rec, capacity)
        self.interaction_capacity = capacity

    def _ensure_pim_nodes(self, needed: int) -> None:
        if needed <= self.pim_capacity:
            return
        capacity = self._next_capacity(self.pim_capacity, needed)
        self.pim_next = _grow(self.pim_next, capacity)
        self.pim_iface = _grow(self.pim_iface, capacity)
        self.pim_capacity = capacity

    # --- array bundles handed to the kernels ------------------------------

    def _pos_arrays(self):
        return (self.pos_position, self.pos_struct, self.pos_chain,
                self.pos_res, self.pos_iface)

    def _int_arrays(self):
        return (self.int_iface, self.int_protein_b, self.int_pos_a, self.int_pos_b,
                self.int_rec_complex, self.int_rec_chain_a, self.int_rec_res_a,
                self.int_rec_chain_b, self.int_rec_res_b)

    def _pim_arrays(self):
        return (self.pim_next, self.pim_iface)

    # --- the aggregation steps -------------------------------------------

    def new_interface(self, structure_id, chain, i_chain, interface_annotation_score) -> int:
        interface_number = self.n_interfaces
        self._ensure_interfaces(interface_number + 1)
        self.iface_score[interface_number] = interface_annotation_score
        self.iface_complex[interface_number] = self.intern(structure_id)
        self.iface_chain[interface_number] = self.intern(chain)
        self.iface_ichain[interface_number] = self.intern(i_chain)
        self.iface_pos_count[interface_number] = 0
        self.n_interfaces = interface_number + 1
        return interface_number

    def set_recommended_complex(self, interface_number, structure_id, chain, i_chain,
                                interface_annotation_score) -> None:
        self.iface_complex[interface_number] = self.intern(structure_id)
        self.iface_chain[interface_number] = self.intern(chain)
        self.iface_ichain[interface_number] = self.intern(i_chain)
        self.iface_score[interface_number] = interface_annotation_score

    def process_entries(self, entry_columns, backmap_columns, stores,
                        struct_idx, chain_idx, compl_db_id, counters) -> None:
        """Runs every interface entry of one (structure, chain) in the kernel.

        Owns the three hand-backs _process_entries cannot do itself: growing the
        interaction arrays, superfusing, and reporting missing quality measures.
        `counters` collects the work counts the caller logs.
        """
        (entry_key_left, entry_key_right, entry_i_chain_idx, entry_score,
         entry_bm_offset, entry_bm_count, has_quality) = entry_columns
        bm_protein, bm_left, bm_length = backmap_columns
        concat_keys, concat_scores, concat_backmaps, backmap = stores

        n_entries = entry_key_left.shape[0]
        if n_entries == 0:
            return

        # At most one new interface per entry, so the interface arrays and the
        # scan scratch can be sized once instead of grown inside the kernel.
        self._ensure_interfaces(self.n_interfaces + n_entries)
        if self.seen_stamp.shape[0] < backmap.shape[0]:
            self.seen_stamp = np.zeros(backmap.shape[0], dtype=np.int64)
        scratch = self.interface_capacity
        scan_overlap = np.zeros(scratch, dtype=np.int64)
        scan_candidates = np.empty(scratch, dtype=np.int64)
        fused_with = np.empty(scratch, dtype=np.int64)

        state = np.zeros(_S_LEN, dtype=np.int64)
        state[_S_POSITIONS] = self.n_positions
        state[_S_INTERACTIONS] = self.n_interactions
        state[_S_PIM_NODES] = self.n_pim_nodes
        state[_S_INTERFACES] = self.n_interfaces
        state[_S_STAMP] = self.scan_stamp
        state[_S_TARGET] = -1
        state[_S_PROTEIN] = -1

        while True:
            # One pass over the remaining keys is the common case, so it is the
            # opening reservation; an entry that fuses with several interfaces
            # walks them again per interface and hands back for more.
            remaining_keys = int(
                (entry_key_right[state[_S_ENTRY]:] - entry_key_left[state[_S_ENTRY]:]).sum())
            self._ensure_positions(int(state[_S_POSITIONS]) + remaining_keys + 1)
            self._ensure_pim_nodes(int(state[_S_PIM_NODES]) + remaining_keys + 1)
            self._ensure_interactions(
                int(state[_S_INTERACTIONS])
                + min(remaining_keys * max(bm_protein.shape[0], 1), _MAX_INTERACTION_RESERVE) + 1)

            status = _process_entries(
                entry_key_left, entry_key_right, entry_i_chain_idx, entry_score,
                entry_bm_offset, entry_bm_count, np.int64(1 if has_quality else 0),
                bm_protein, bm_left, bm_length,
                concat_keys, concat_scores, concat_backmaps, backmap,
                np.int64(struct_idx), np.int64(chain_idx), np.int64(compl_db_id),
                scan_overlap, scan_candidates, fused_with, self.seen_stamp,
                self.iface_score, self.iface_complex, self.iface_chain,
                self.iface_ichain, self.iface_pos_count,
                self.pos_index, self._pos_arrays(),
                self.int_index, self._int_arrays(), self.int_score, self.int_has_rec,
                np.int64(self.interaction_capacity),
                np.int64(self.position_capacity), np.int64(self.pim_capacity),
                self.pim_head, self._pim_arrays(),
                state)

            self.n_positions = int(state[_S_POSITIONS])
            self.n_interactions = int(state[_S_INTERACTIONS])
            self.n_pim_nodes = int(state[_S_PIM_NODES])
            self.n_interfaces = int(state[_S_INTERFACES])
            self.scan_stamp = int(state[_S_STAMP])

            if status == _ENTRIES_DONE:
                break
            if status == _ENTRIES_MISSING_QUALITY:
                raise KeyError('quality measures missing for the annotated structure')
            if status == _ENTRIES_NEED_SUPERFUSION:
                self.superfuse(fused_with[:int(state[_S_N_FUSED])].tolist())
                # The merge renumbers the interfaces, replaces the index dicts and
                # compacts the position, interaction and pim rows, so every count
                # in the cursor has to come back from the store. The cursor already
                # points at the next entry.
                state[_S_POSITIONS] = self.n_positions
                state[_S_INTERACTIONS] = self.n_interactions
                state[_S_PIM_NODES] = self.n_pim_nodes
                state[_S_INTERFACES] = self.n_interfaces
                self._ensure_interfaces(self.n_interfaces + n_entries)
                scratch = self.interface_capacity
                if scan_overlap.shape[0] < scratch:
                    scan_overlap = np.zeros(scratch, dtype=np.int64)
                    scan_candidates = np.empty(scratch, dtype=np.int64)
                    fused_with = np.empty(scratch, dtype=np.int64)
                continue
            if status == _ENTRIES_NEED_ROW_CAPACITY:
                self._ensure_positions(self.n_positions + 1)
                self._ensure_pim_nodes(self.n_pim_nodes + 1)
                continue
            # _ENTRIES_NEED_CAPACITY: grown at the top of the next round.
            self._ensure_interactions(self.n_interactions + 1)

        counters[0] += int(state[_S_KEYS_TOTAL])
        counters[1] += int(state[_S_FUSION_PASSES])
        counters[2] += int(state[_S_NEW_PASSES])
        counters[3] += int(state[_S_KEY_ITERATIONS])
        counters[4] += int(state[_S_CREATED])
        counters[5] += int(state[_S_SKIPPED])
        counters[6] += int(state[_S_SKIP_3])

    def superfuse(self, fused_with: list[int]) -> None:
        (self.pos_index, self.int_index, self.n_interfaces,
         self.n_positions, self.n_interactions, self.n_pim_nodes) = _superfuse(
            np.asarray(sorted(fused_with), dtype=np.int64), np.int64(self.n_interfaces),
            self.iface_score, self.iface_complex, self.iface_chain, self.iface_ichain,
            self.iface_pos_count,
            self.pos_index, self._pos_arrays(), np.int64(self.n_positions),
            self.int_index, self._int_arrays(), self.int_score, self.int_has_rec,
            np.int64(self.n_interactions),
            self.pim_head, self._pim_arrays(), np.int64(self.n_pim_nodes))
        self.n_interfaces = int(self.n_interfaces)
        self.n_positions = int(self.n_positions)
        self.n_interactions = int(self.n_interactions)
        self.n_pim_nodes = int(self.n_pim_nodes)

    # --- reporting --------------------------------------------------------

    def describe(self, interface_number: int) -> str:
        return (f'interface {interface_number}: '
                f'recommended_complex={self.strings[self.iface_complex[interface_number]]} '
                f'{self.strings[self.iface_chain[interface_number]]} '
                f'{self.strings[self.iface_ichain[interface_number]]}, '
                f'positions={self.iface_pos_count[interface_number]}, '
                f'score={self.iface_score[interface_number]}')

    def materialize(self, protein_db_id) -> list['Aggregated_interface']:
        """The Aggregated_interface list the rest of the pipeline consumes."""
        interfaces: list[Aggregated_interface] = []
        if self.n_interfaces == 0:
            return interfaces

        numbers = np.arange(self.n_interfaces)

        pos_iface = self.pos_iface[:self.n_positions]
        pos_order = np.argsort(pos_iface, kind='stable')
        pos_sorted = pos_iface[pos_order]
        pos_starts = np.searchsorted(pos_sorted, numbers, side='left')
        pos_ends = np.searchsorted(pos_sorted, numbers, side='right')

        int_iface = self.int_iface[:self.n_interactions]
        int_order = np.argsort(int_iface, kind='stable')
        int_sorted = int_iface[int_order]
        int_starts = np.searchsorted(int_sorted, numbers, side='left')
        int_ends = np.searchsorted(int_sorted, numbers, side='right')

        strings = self.strings

        for interface_number in range(self.n_interfaces):
            rows = pos_order[pos_starts[interface_number]:pos_ends[interface_number]]
            # Python ints, not numpy scalars: these keys end up in the
            # Position column of RS_Position_Interface and feed
            # get_interface_hash, both of which want unbounded ints.
            position_ids = self.pos_position[rows].tolist()
            struct_ids = self.pos_struct[rows].tolist()
            chain_ids = self.pos_chain[rows].tolist()
            residues = self.pos_res[rows].tolist()
            positions = {}
            for row_index in range(rows.shape[0]):
                positions[position_ids[row_index]] = (
                    strings[struct_ids[row_index]], strings[chain_ids[row_index]],
                    residues[row_index])

            interface = Aggregated_interface(
                protein_db_id,
                recommended_complex=strings[self.iface_complex[interface_number]],
                chain=strings[self.iface_chain[interface_number]],
                interacting_chain=strings[self.iface_ichain[interface_number]],
                positions=positions,
                interface_annotation_score=self.iface_score[interface_number])

            rows = int_order[int_starts[interface_number]:int_ends[interface_number]]
            interface.interaction_count = int(rows.shape[0])
            interface.interaction_key_array = np.stack(
                (self.int_protein_b[rows], self.int_pos_a[rows], self.int_pos_b[rows]),
                axis=1)
            interface.interaction_scores = self.int_score[rows]
            interface.interaction_rec_complex = self.int_rec_complex[rows]
            interface.interaction_rec_chain_a = self.int_rec_chain_a[rows]
            interface.interaction_rec_res_a = self.int_rec_res_a[rows]
            interface.interaction_rec_chain_b = self.int_rec_chain_b[rows]
            interface.interaction_rec_res_b = self.int_rec_res_b[rows]
            interface.interaction_has_rec = self.int_has_rec[rows]
            interface.chain_ids = strings

            interfaces.append(interface)

        return interfaces


class Aggregated_interface(Slotted_obj):
    """One aggregated interface, as read by the insert and output code.

    Built either by Interface_Store.materialize or, for interfaces coming back
    from the database, directly with a positions dict and no interactions. The
    interaction rows are flat arrays; the recommended interaction is decomposed
    into (chain_a, res_a, chain_b, res_b) with the chains interned in chain_ids,
    and protein_a is not stored because it is always this interface's protein.
    """

    __slots__ = [
        'protein',  'recommended_complex', 'chain',
        'interacting_chain', 'positions', 'interface_annotation_score',
        'interaction_key_array', 'interaction_scores',
        'interaction_rec_complex', 'interaction_rec_chain_a', 'interaction_rec_res_a',
        'interaction_rec_chain_b', 'interaction_rec_res_b', 'interaction_has_rec',
        'interaction_count', 'interaction_protein_groups', 'chain_ids',
        'database_id', 'interface_hash']

    def __init__(self, protein_id = None, recommended_complex = None, chain = None, interacting_chain = None, positions = None, interface_annotation_score = None, database_id = None):
        self.protein = protein_id
        self.recommended_complex = recommended_complex
        self.chain = chain
        self.interacting_chain = interacting_chain
        self.interface_annotation_score = interface_annotation_score
        if positions is None:
            self.positions = {}
        else:
            self.positions = positions

        self.interaction_count = 0
        self.interaction_key_array = np.empty((0, 3), dtype=np.int64)
        self.interaction_scores = np.empty(0, dtype=np.float64)
        self.interaction_rec_complex = np.empty(0, dtype=np.int64)
        self.interaction_rec_chain_a = np.empty(0, dtype=np.int64)
        self.interaction_rec_res_a = np.empty(0, dtype=np.int64)
        self.interaction_rec_chain_b = np.empty(0, dtype=np.int64)
        self.interaction_rec_res_b = np.empty(0, dtype=np.int64)
        self.interaction_has_rec = np.empty(0, dtype=np.uint8)
        self.interaction_protein_groups: dict[int, list[int]] | None = None
        self.chain_ids: list[str] = []

        self.database_id = database_id
        self.interface_hash = None

    def deconstruct(self):
        del self.positions

        doomsday_protocol(self)

    def add_position(self, position, recommended_residue):
        self.positions[position] = recommended_residue

    # --- interaction access -------------------------------------------------

    def has_interactions(self) -> bool:
        return self.interaction_count > 0

    def _protein_groups(self) -> dict:
        # Derived on demand instead of maintained per insert: the aggregation
        # kernels append rows and cannot touch a python dict. Runs once per
        # interface at insert time.
        if self.interaction_protein_groups is None:
            groups: dict[int, list[int]] = {}
            proteins = self.interaction_key_array[:self.interaction_count, 0]
            for row, protein_b in enumerate(proteins.tolist()):
                rows = groups.get(protein_b)
                if rows is None:
                    groups[protein_b] = [row]
                else:
                    rows.append(row)
            self.interaction_protein_groups = groups
        return self.interaction_protein_groups

    def number_of_interacting_proteins(self) -> int:
        return len(self._protein_groups())

    def interacting_protein_ids(self):
        return self._protein_groups().keys()

    def _build_interaction(self, row: int) -> Position_Position_Interaction:
        protein_b, position_a, position_b = self.interaction_key_array[row].tolist()
        interaction = Position_Position_Interaction(self.protein, position_a, protein_b, position_b)
        interaction.interaction_score = float(self.interaction_scores[row])
        # Left unset when nothing was recorded, exactly as an untouched
        # Position_Position_Interaction slot behaved.
        if self.interaction_has_rec[row]:
            interaction.recommended_complex = int(self.interaction_rec_complex[row])
            interaction.recommended_interaction = (
                self.chain_ids[self.interaction_rec_chain_a[row]],
                int(self.interaction_rec_res_a[row]),
                self.chain_ids[self.interaction_rec_chain_b[row]],
                int(self.interaction_rec_res_b[row]),
            )
        return interaction

    def nested_interactions_for(self, protein_b) -> dict:
        # {position_a: {position_b: Position_Position_Interaction}} is the on-disk shape
        # of Pos_Pos_Data, so it is rebuilt here on demand. Once per interacting protein
        # at insert time, not on the hot path.
        nested: dict[int, dict[int, Position_Position_Interaction]] = {}
        for row in self._protein_groups()[protein_b]:
            position_a = int(self.interaction_key_array[row, 1])
            position_b = int(self.interaction_key_array[row, 2])
            positions_for_a = nested.get(position_a)
            if positions_for_a is None:
                positions_for_a = {}
                nested[position_a] = positions_for_a
            positions_for_a[position_b] = self._build_interaction(row)
        return nested

    def get_interaction(self, protein_b, position_a, position_b) -> Position_Position_Interaction:
        for row in self._protein_groups()[protein_b]:
            if (self.interaction_key_array[row, 1] == position_a
                    and self.interaction_key_array[row, 2] == position_b):
                return self._build_interaction(row)
        raise KeyError((protein_b, position_a, position_b))

    def set_recommended_complex(self, recommended_complex, chain, interacting_chain, interface_annotation_score):
        self.recommended_complex = recommended_complex
        self.chain = chain
        self.interacting_chain = interacting_chain
        self.interface_annotation_score = interface_annotation_score

    def get_interface_hash(self) -> int:
        """Re-find an Interface row right after it was inserted blank.

        The keys of self.positions are position *database ids*, and the whole
        set of them identifies the interface: the hash only has to separate the
        interfaces of one protein from each other (see blank_interfaces in
        insertion_lib.remote_insert_interfaces), and aggregation fuses any two
        interfaces that overlap by more than half their positions, so two
        interfaces that reach insertion differ in at least half of theirs.

        Every id is avalanched on its own and the results are added up mod
        2**64, which makes the digest a function of the key *set* and not of
        the order dict iteration happens to hand them over in. The count goes
        in as well, then a final avalanche spreads it over the whole width. The
        result is folded to 63 bits: the column is bigint unsigned in MySQL but
        a signed 64 bit integer in SQLite, which refuses anything above
        2**63 - 1.

        Interfaces of one protein number in the tens, so with 63 bits of
        digest, a collision inside one protein is ~1e-17 likely; insertion
        checks for one anyway, since it cannot be recovered from.
        """
        if self.interface_hash is not None:
            return self.interface_hash
        if len(self.positions) == 0:
            return 0
        digest = _splitmix64_sum(self.positions, len(self.positions))
        digest = _splitmix64(digest ^ _splitmix64(_HASH_SEED ^ len(self.positions)))
        self.interface_hash = digest & _MASK_63
        return self.interface_hash

    def print_interface(self):
        max_prints = 100
        print('==============================================')
        print(f'Aggregated interface object: {self.protein}')
        print(f'Recommended complex: {self.recommended_complex} {self.chain} {self.interacting_chain}')
        if len(self.positions) < max_prints:
            print(f'Positions:\n{self.positions}')
        else:
            print(f'Positions, in total {len(self.positions)}')
        print(f'Interface score: {self.interface_annotation_score}')


        print('==============================================')

class Interface(Slotted_obj):
    __slots__ = ['chain', 'interacting_chain', 'residues', 'support_residues', 'interactions', 'stored']

    def __init__(self, chain = None, interacting_chain = None, stored = False):
        self.chain: str | None = chain
        self.interacting_chain: str | None = interacting_chain
        self.interactions: dict[int, dict[int, float]] = {}
        self.residues: set[int] = set()
        self.support_residues: set[int] = set()
        self.stored: bool = stored

    def deconstruct(self):
        del self.residues
        del self.interactions
        del self.support_residues
        doomsday_protocol(self)

    def add_interaction(self, res_a: int, res_b: int, score: float):
        if res_a not in self.interactions:
            self.interactions[res_a] = {}
        self.interactions[res_a][res_b] = score
        self.residues.add(res_a)

    def add_support(self, res):
        self.support_residues.add(res)
        self.residues.add(res)


def interface_mat_lookup(res_a: int, res_b: int, keys: np.ndarray[np.int64], scores: np.ndarray[np.float32]):
    target = (np.int64(res_a) << 32) | np.int64(res_b)
    idx = np.searchsorted(keys, target)
    if idx < len(keys) and keys[idx] == target:
        return scores[idx]
    return np.nan

def decode_key(combined_key):
    # Shift right to get the high 32 bits
    row = combined_key >> 32

    # Use a 32-bit mask (0xFFFFFFFF) to get the low 32 bits
    col = combined_key & 0xFFFFFFFF

    return int(row), int(col)


# --- replay harness ----------------------------------------------------------
#
# Opt-in capture of calculate_aggregated_interface_map's inputs and outputs, so a
# refactor of the interface data model can be checked by replaying recorded calls
# instead of by a full pipeline run. Off unless STRUCTMAN_INTERFACE_DUMP names a
# directory; STRUCTMAN_INTERFACE_DUMP_N caps how many calls each worker records
# (default 2). The large shared stores are identical for every call in a process,
# so they are written once per pid and the per-call payload stays small.
#
# structman.scripts.synthetic_interface_corpus writes the same format from
# fabricated inputs, for when no pipeline run is at hand.

_INTERFACE_DUMP_DIR: str | None = os.environ.get('STRUCTMAN_INTERFACE_DUMP') or None
_INTERFACE_DUMP_LIMIT: int = int(os.environ.get('STRUCTMAN_INTERFACE_DUMP_N') or 2)
_interface_dump_calls: int = 0
_interface_dump_shared_written: bool = False


def _interface_dump_enabled() -> bool:
    return _INTERFACE_DUMP_DIR is not None and _interface_dump_calls < _INTERFACE_DUMP_LIMIT


def _interface_dump_write(shared: dict, call: dict) -> None:
    global _interface_dump_calls, _interface_dump_shared_written

    os.makedirs(_INTERFACE_DUMP_DIR, exist_ok=True)
    pid = os.getpid()

    if not _interface_dump_shared_written:
        shared_path = os.path.join(_INTERFACE_DUMP_DIR, f'shared_{pid}.pkl')
        tmp = f'{shared_path}.{pid}.tmp'
        with open(tmp, 'wb') as f:
            pickle.dump(shared, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, shared_path)
        _interface_dump_shared_written = True

    call_path = os.path.join(_INTERFACE_DUMP_DIR, f'call_{pid}_{_interface_dump_calls}.pkl')
    tmp = f'{call_path}.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(call, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, call_path)

    _interface_dump_calls += 1


def summarize_interfaces(aggregated_interfaces: list) -> list:
    """Canonical, comparable form of the function's result, for diffing replays."""
    summary = []
    for interface in aggregated_interfaces:
        interactions = []
        for prot_b in sorted(interface.interacting_protein_ids()):
            interactions_for_b = interface.nested_interactions_for(prot_b)
            for pos_a in sorted(interactions_for_b):
                for pos_b in sorted(interactions_for_b[pos_a]):
                    obj = interactions_for_b[pos_a][pos_b]
                    interactions.append((
                        prot_b, pos_a, pos_b,
                        getattr(obj, 'interaction_score', None),
                        getattr(obj, 'recommended_complex', None),
                        getattr(obj, 'recommended_interaction', None),
                    ))
        summary.append({
            'protein': interface.protein,
            'recommended_complex': interface.recommended_complex,
            'chain': interface.chain,
            'interacting_chain': interface.interacting_chain,
            'interface_annotation_score': interface.interface_annotation_score,
            'positions': sorted(interface.positions.items()),
            'interactions': interactions,
        })
    return summary


def interface_annotation_limits(config) -> tuple[int, int, float]:
    """(keep_all, max_annotations, min_seq_id) for select_interface_annotations.

    getattr with the Config defaults, like the other aggregation knobs, so that a
    config not built by Config.__init__ still works.
    """
    return (getattr(config, 'interface_annotations_keep_all', 100),
            getattr(config, 'interface_annotations_max', 1000),
            getattr(config, 'interface_annotations_min_seq_id', 0.5))


def select_interface_annotations(
        struct_db_ids: np.ndarray,
        seq_ids: np.ndarray,
        multi_chain: np.ndarray,
        keep_all: int,
        max_annotations: int,
        min_seq_id: float) -> np.ndarray:
    """Which of a protein's structure annotations take part in the interface calculation.

    This decides both which alignments keep their backmap in the database
    (database.prune_interface_backmaps) and which structures para_classify hands to
    calculate_aggregated_interface_map, so that the two agree. The three columns hold
    one entry per annotation of the protein; the result is a mask over them.

    1. An annotation whose complex has a single protein chain (multi_chain False)
       never takes part, since nothing can interact with that chain, and is not
       counted below either.
    2. Up to keep_all remaining annotations are all taken.
    3. Up to max_annotations, the keep_all best by sequence identity are taken, plus
       every one above min_seq_id.
    4. Beyond that, the max_annotations best.

    Best means the highest sequence identity, NaN (none stored) last, and among equal
    identities the lower structure db id, so the result does not depend on the order
    the rows arrive in.
    """
    selected = np.array(multi_chain, dtype=bool)
    candidates = np.flatnonzero(selected)
    n_candidates = candidates.shape[0]
    if n_candidates <= keep_all:
        return selected

    candidate_seq_ids = np.asarray(seq_ids, dtype=np.float64)[candidates]
    candidate_struct_db_ids = np.asarray(struct_db_ids, dtype=np.int64)[candidates]
    rank_key = np.where(np.isnan(candidate_seq_ids), -np.inf, candidate_seq_ids)
    ranked = candidates[np.lexsort((candidate_struct_db_ids, -rank_key))]

    selected[:] = False
    if n_candidates > max_annotations:
        selected[ranked[:max_annotations]] = True
    else:
        selected[ranked[:keep_all]] = True
        # NaN compares False, so an annotation without identity only gets in by rank.
        selected[candidates[candidate_seq_ids > min_seq_id]] = True
    return selected


def calculate_aggregated_interface_map(
        config: Config,
        protein_db_id: int,
        annotations: dict[str, list[tuple[str, int]]],
        structs_from_db: np.ndarray[int],
        index_list_stores: list[np.ndarray[int]],
        chain_list_stores: list[np.ndarray[bytes]],
        concat_combined_keys_stores: list[np.ndarray[int]],
        concat_scores_stores: list[np.ndarray[np.float32]],
        #interface_map: dict[int, tuple[list[tuple[str, np.ndarray[float]]], int, float]],
        #stored_interfaces: tuple[dict[int, np.ndarray[np.int64]], dict[int, np.ndarray[float]]],
        quality_measures: dict[int, tuple[float, float]],
        compl_struct_db_id_matrix: np.ndarray[int],
        compl_db_id_map: np.ndarray[int],
        structure_db_list: np.ndarray[bytes],
        backmap_store: tuple[np.ndarray, np.ndarray, np.ndarray[int]],
        ) -> list[Aggregated_interface]:
    dump_this_call = _interface_dump_enabled()
    if dump_this_call:
        # Captured before anything runs, so a replay starts from the true pre-state.
        dump_shared = {
            'structs_from_db': structs_from_db,
            'index_list_stores': index_list_stores,
            'chain_list_stores': chain_list_stores,
            'concat_combined_keys_stores': concat_combined_keys_stores,
            'concat_scores_stores': concat_scores_stores,
            'compl_struct_db_id_matrix': compl_struct_db_id_matrix,
            'compl_db_id_map': compl_db_id_map,
            'structure_db_list': structure_db_list,
            'backmap_store': backmap_store,
        }
        dump_call = {
            'protein_db_id': protein_db_id,
            'annotations': annotations,
            'quality_measures': quality_measures,
            'verbosity': config.verbosity,
        }

    store = Interface_Store()

    if config.verbosity >= 4:
        config.logger.info(f'Call of calculate_aggregated_interface_map: {protein_db_id} {len(annotations)=} {len(backmap_store)=}')

    n_interfaces_0 = 0
    n_interfaces_1 = 0
    n_interfaces_2 = 0

    # Work counters for the interface core loop: how many keys it walks, how often the
    # fusion branch re-walks the same key array, and how many pos-pos interactions come
    # out. The entry kernel keeps them in this array and adds them up per chain, so
    # they cost nothing next to the loop itself.
    # [keys, fusion passes, new-interface passes, key iterations, created, skipped,
    #  entries whose interacting chain is not in the complex]
    entry_counters = [0] * 7
    count_skip_0 = 0
    count_skip_1 = 0
    count_skip_2 = 0
    total_loop_times = []
    retrieved_interface_data: dict[int, tuple[int, np.ndarray, list[str]]] = {}
    # int64 view of each key store, taken once rather than per chain. It is a no-op
    # when the store is already int64, which it is; the cast only pins the width the
    # kernel's shift and mask are written for.
    keys_of_store: dict[int, np.ndarray] = {}

    (concat_backmaps, prot_db_id_vec, structure_db_id_vec) = backmap_store

    if config.verbosity >= 7:
        config.logger.info(f'In calculate_aggregated_interface_map: {protein_db_id=} {concat_backmaps=} {prot_db_id_vec=} {structure_db_id_vec=}')

    for structure_id in annotations:
        example_structure_db_id = annotations[structure_id][0][1]

        compl_db_id = compl_db_id_map[example_structure_db_id]
        if compl_db_id == -1:
            continue
        structure_db_ids: np.ndarray[int] = compl_struct_db_id_matrix[compl_db_id]

        loop_times = []
        t_loop = time.time()

        interfaces: dict[int, tuple[int, np.ndarray, list[str]]] = {}
        structure_backmaps: dict[int, dict[int, tuple[int, int]]] = {}

        chain_to_s_db_ids = {}

        for structure_db_id in structure_db_ids:
            if structure_db_id == -1:
                continue
            fused_id = structure_db_list[structure_db_id]
            chain = fused_id[-1:].decode()

            chain_to_s_db_ids[chain] = structure_db_id
            store_index, _, _, index_left, index_right, _, _ = structs_from_db[structure_db_id]
            if store_index == -1:
                continue

            if structure_db_id not in retrieved_interface_data:
                # Bounds, not slices: the entry kernel indexes the store arrays
                # directly, so an entry costs nothing to describe.
                bounds = np.ascontiguousarray(
                    index_list_stores[store_index][index_left:index_right], dtype=np.int64)
                i_chains = [chain_list_stores[store_index][index_left + index].decode()
                            for index in range(bounds.shape[0])]
                retrieved_interface_data[structure_db_id] = (store_index, bounds, i_chains)
            interfaces[structure_db_id] = retrieved_interface_data[structure_db_id]

        t_loop = add_to_times(loop_times, t_loop, 'ray.get interfacemap')

        for structure_db_id in structure_db_ids:
            if structure_db_id == -1:
                continue

            prot_db_id_vec_left, prot_db_id_vec_right = structure_db_id_vec[structure_db_id]
            if prot_db_id_vec_left != -1:
                # The (left, right) bounds rather than the slice: the interface pass
                # reads the interacting proteins' positions straight out of
                # concat_backmaps, and the target protein's slice is taken once per
                # chain below.
                bm_dict: dict[int, tuple[int, int]] = {}
                for prot_triple in prot_db_id_vec[prot_db_id_vec_left:prot_db_id_vec_right]:
                    prot_db_id, backmap_left, backmap_right = prot_triple
                    bm_dict[prot_db_id] = (backmap_left, backmap_right)
                structure_backmaps[structure_db_id] = bm_dict

        if config.verbosity >= 7:
            config.logger.info(f'In calculate_aggregated_interface_map: {protein_db_id=} {structure_id=} {structure_db_ids=} {[(s_db_id, structure_backmaps[s_db_id]) for s_db_id in structure_db_ids if s_db_id in structure_backmaps]=}\n{interfaces=}')

        backmap_columns, backmap_span = _complex_backmap_columns(structure_backmaps)

        t_loop = add_to_times(loop_times, t_loop, 'check backmap_store')

        n_interfaces_0 += len(interfaces)

        total_inner_loop_times = []
        for chain, t_structure_db_id in annotations[structure_id]:

            if t_structure_db_id not in structure_backmaps:
                count_skip_1 += 1
                if config.verbosity >= 6:
                    config.logger.info(f'skip in calculate_aggregated_interface_map: {protein_db_id}, {structure_id}, {chain}, {list(structure_backmaps.keys())=}')
                continue
            if protein_db_id not in structure_backmaps[t_structure_db_id]:
                count_skip_2 += 1
                if config.verbosity >= 6:
                    #This is quite normal for chains without interfaces
                    config.logger.info(f'skip 2 in calculate_aggregated_interface_map: {protein_db_id}, {structure_id}, {chain}, {list(structure_backmaps[t_structure_db_id].keys())=}')
                continue

            backmap_left, backmap_right = structure_backmaps[t_structure_db_id][protein_db_id]
            backmap: np.ndarray = concat_backmaps[backmap_left:backmap_right]

            if config.verbosity >= 6:
                config.logger.info(f'\ncalculate_aggregated_interface_map: {protein_db_id}, {structure_id}, {chain}, {interfaces.keys()=}')

            if interfaces is None:
                config.errorlog.add_error(f'Interfaces is None: {protein_db_id} {structure_id} {chain}')

            n_interfaces_1 += 1
            if t_structure_db_id not in interfaces:
                continue
            if interfaces[t_structure_db_id] is None:
                continue

            compl_db_id = structs_from_db[t_structure_db_id][5]
            store_index, entry_bounds, entry_i_chains = interfaces[t_structure_db_id]
            if store_index not in keys_of_store:
                keys_of_store[store_index] = np.asarray(
                    concat_combined_keys_stores[store_index], dtype=np.int64)
            n_interfaces_2 += len(entry_i_chains)

            inner_loop_times = []
            t_inner = time.time()

            # The whole entry loop runs in one kernel now. What is left here is
            # O(number of entries): resolving each entry's interacting chain, its
            # annotation score, and where its interacting proteins' backmaps sit.
            # The per-entry logging that used to be inside the loop went with it;
            # the counters below and the summary at the end still describe the work.
            n_entries = len(entry_i_chains)
            entry_i_chain_idx = np.empty(n_entries, dtype=np.int64)
            entry_score = np.zeros(n_entries, dtype=np.float64)
            entry_bm_offset = np.zeros(n_entries, dtype=np.int64)
            entry_bm_count = np.zeros(n_entries, dtype=np.int64)

            # calc_interface_annotation_score stays in python: it is one nansum per
            # entry, and numpy's pairwise summation over the float32 scores is not
            # what a plain loop in a kernel would produce.
            quality = quality_measures.get(t_structure_db_id)
            has_quality = quality is not None
            if has_quality:
                seq_id, cov = quality

            entry_scores_store = concat_scores_stores[store_index]
            for entry in range(n_entries):
                i_structure_db_id = chain_to_s_db_ids.get(entry_i_chains[entry], -1)
                if i_structure_db_id == -1:
                    entry_i_chain_idx[entry] = -1
                    continue
                entry_i_chain_idx[entry] = store.intern(entry_i_chains[entry])
                entry_bm_offset[entry], entry_bm_count[entry] = backmap_span.get(
                    i_structure_db_id, (0, 0))

                if has_quality:
                    left, right = entry_bounds[entry]
                    entry_score[entry] = calc_interface_annotation_score(
                        seq_id, cov, entry_scores_store[left:right])

            try:
                store.process_entries(
                    (np.ascontiguousarray(entry_bounds[:, 0]),
                     np.ascontiguousarray(entry_bounds[:, 1]),
                     entry_i_chain_idx, entry_score, entry_bm_offset, entry_bm_count,
                     has_quality),
                    backmap_columns,
                    (keys_of_store[store_index], entry_scores_store,
                     concat_backmaps, backmap),
                    store.intern(structure_id), store.intern(chain), compl_db_id,
                    entry_counters)
            except KeyError as err:
                raise KeyError(f'{err=} {structure_id=} {chain=} {t_structure_db_id=} {protein_db_id=} {len(quality_measures)=}')

            t_inner = add_to_times(inner_loop_times, t_inner, 'inner loop')

            total_inner_loop_times = aggregate_times(total_inner_loop_times, inner_loop_times)

        loop_times.append(total_inner_loop_times)

        t_loop = add_to_times(loop_times, t_loop, 'outer loop')

        total_loop_times = aggregate_times(total_loop_times, loop_times)

    aggregated_interfaces: list[Aggregated_interface] = store.materialize(protein_db_id)

    if config.verbosity >= 4:
        config.logger.info(
            f'[interface][loop] {protein_db_id=} keys={entry_counters[0]:_} '
            f'fusion_passes={entry_counters[1]:_} new_interface_passes={entry_counters[2]:_} '
            f'key_iterations={entry_counters[3]:_} '
            f'pos_pos_created={entry_counters[4]:_} pos_pos_skipped={entry_counters[5]:_} '
            f'interfaces={len(aggregated_interfaces)}'
        )

    if config.verbosity >= 4:
        config.logger.info(f'End of calculate_aggregated_interface_map: {protein_db_id=} {n_interfaces_0=} {n_interfaces_1=} {n_interfaces_2=} {count_skip_0=} {count_skip_1=} {count_skip_2=} count_skip_3={entry_counters[6]} {len(aggregated_interfaces)=}')

    if dump_this_call:
        dump_call['result'] = summarize_interfaces(aggregated_interfaces)
        _interface_dump_write(dump_shared, dump_call)

    return aggregated_interfaces, total_loop_times
