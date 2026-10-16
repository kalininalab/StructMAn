import gc
import os
import random
import ray
import shutil
import string
import subprocess
import time
from array import array
from typing import NamedTuple

import numpy as np

from structman.scripts import createCustomDb
from structman.base_utils.base_utils import build_sequence_store
from structman.base_utils.config_class import Config
from structman.lib.sdsc.protein import Proteins

# The hit columns are grown as array.array and then reinterpreted as numpy without a
# copy, so the typecodes and the dtypes below have to agree on width. 'i' is 4 bytes
# and 'q' is 8 on every platform CPython builds for; dev_tests/test_hit_store.py
# asserts it rather than leaving it to a silent reinterpretation.
_CODE_TYPECODE = 'i'      # -> np.int32
_OFFSET_TYPECODE = 'q'    # -> np.int64
_CHAIN_TYPECODE = 'B'     # -> np.uint8

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


def geneSeqMapToFasta(proteins, outfile, config: Config, filtering_db=None, custom_prot_ids=None):
    lines = []
    n = 0
    m = 0

    if custom_prot_ids is None:
        prot_ids = proteins.get_protein_ids()
    else:
        prot_ids = custom_prot_ids

    n_not_stored = 0

    for prot_id in prot_ids:

        # Previous if-statement:
        # if proteins[prot_id].stored and prot_id not in proteins.indels:
        # Didnt understood why proteins in the indels dict should be forced into the mmseqs search
        # Remove for now, retrieve it, if something comes up, then better comment it
        if proteins[prot_id].stored:
            continue
        n_not_stored += 1

        seq = proteins.get_sequence(prot_id)

        if seq == 0 or seq == 1 or seq == '' or seq is None:
            continue
        if filtering_db is not None:
            folder_key = prot_id.split('-')[0][-2:]
            filename = '%s/%s/%s_ref50_gpw.fasta.gz' % (filtering_db, folder_key, prot_id)
            if os.path.isfile(filename):
                n += 1
                continue

        lines.append('>%s' % prot_id)
        lines.append(seq)
        m += 1

    if n_not_stored == 0:
        return False, m

    if len(lines) > 0:
        if config.verbosity >= 2:
            config.logger.info(f'Filtered {n} Proteins before mmseqs2')
            config.logger.info(f'{m} sequences going into mmseqs2')
        f = open(outfile, 'w')
        f.write('\n'.join(lines))
        f.close()
        return True, m
    else:
        return 'Empty fasta file', m


class HitStore(NamedTuple):
    """The mmseqs hit map as numpy columns instead of nested dicts of tuples of sets.

    On a whole-PDB chunk this map holds hundreds of millions of hits, and it is built
    in ray workers and handed to the driver through ray.get. As python containers that
    cost ~400 bytes a hit and had to be pickled out of the worker and rebuilt object by
    object in the driver, because ray's zero-copy path only applies to numpy buffers.
    As columns it costs ~13 bytes a hit and crosses the object store as raw buffers.

    Rows are grouped by query protein, CSR style: the hits of protein `p` are the rows
    `row_offsets[p]:row_offsets[p + 1]`, and the protein ids themselves are a packed
    utf-8 blob with its own offsets (the same shape base_utils.build_sequence_store
    already uses for sequences). PDB ids repeat enormously while the distinct set is
    tiny, so they live in a string table and every row stores only its int32 code.
    Chain names are a single ASCII byte, since parseHits drops every hit whose chain
    name is longer than one character.

    `aln_length` only exists so merge_parsed_chunks can keep its max-alignment-length
    tie-break; nothing downstream reads it, and search drops it before handing the
    store on (see drop_aln_length).
    """

    # one entry per query protein, in first-appearance order
    prot_blob: np.ndarray       # uint8, packed protein ids
    prot_offsets: np.ndarray    # int64, n_prots + 1
    row_offsets: np.ndarray     # int64, n_prots + 1, CSR into the hit columns

    # the pdb id string table
    pdb_blob: np.ndarray        # uint8, packed pdb ids
    pdb_offsets: np.ndarray     # int64, n_pdbs + 1

    # one entry per hit
    pdb_codes: np.ndarray       # int32, index into the pdb table
    chains: np.ndarray          # uint8, the chain name's single ASCII byte
    target_len: np.ndarray      # int32, the hit structure's sequence length
    aln_length: np.ndarray      # int32, merge-only, see above

    @property
    def n_prots(self) -> int:
        return len(self.row_offsets) - 1

    @property
    def n_hits(self) -> int:
        return len(self.pdb_codes)


def empty_hit_store() -> HitStore:
    return HitStore(
        prot_blob = np.empty(0, dtype = np.uint8),
        prot_offsets = np.zeros(1, dtype = np.int64),
        row_offsets = np.zeros(1, dtype = np.int64),
        pdb_blob = np.empty(0, dtype = np.uint8),
        pdb_offsets = np.zeros(1, dtype = np.int64),
        pdb_codes = np.empty(0, dtype = np.int32),
        chains = np.empty(0, dtype = np.uint8),
        target_len = np.empty(0, dtype = np.int32),
        aln_length = np.empty(0, dtype = np.int32),
    )


def _decode_blob(blob: np.ndarray, offsets: np.ndarray) -> list[str]:
    """The packed ids of one blob as a python list, in code order.

    One tobytes() for the whole blob and then plain bytes slicing: taking a numpy
    slice per id would allocate an array object per id, and these tables are read
    once per store by callers that then index them millions of times.
    """
    raw: bytes = blob.tobytes()
    bounds: list[int] = offsets.tolist()
    return [raw[bounds[i]:bounds[i + 1]].decode('utf-8') for i in range(len(bounds) - 1)]


def decode_prot_ids(store: HitStore) -> list[str]:
    return _decode_blob(store.prot_blob, store.prot_offsets)


def decode_pdb_table(store: HitStore) -> list[str]:
    return _decode_blob(store.pdb_blob, store.pdb_offsets)


def drop_aln_length(store: HitStore) -> HitStore:
    """aln_length is decided by the time the store leaves search, and no consumer reads
    it, so it does not need to stay resident for the rest of the run."""
    return store._replace(aln_length = np.empty(0, dtype = np.int32))


def build_self_hit_store(prot_ids: list[str], seq_lens: list[int]) -> HitStore:
    """PDB chains hitting their own structures, in the shape parseHits produces.

    Both callers used to spell this out themselves and both were left behind when the
    hit map changed shape last time, so there is one definition. A chain name that is
    not a single character contributes no hit - the same rule parseHits applies - but
    the protein still gets its (empty) block, so the membership tests in search see it.
    """
    pdb_table: list[str] = []
    pdb_code_of: dict[str, int] = {}

    pdb_codes = array(_CODE_TYPECODE)
    chains = array(_CHAIN_TYPECODE)
    target_len = array(_CODE_TYPECODE)
    row_offsets = array(_OFFSET_TYPECODE, [0])

    for prot_id, seq_len in zip(prot_ids, seq_lens):
        pdb_id, _, chain = prot_id.partition(':')
        if len(chain) == 1:
            code = pdb_code_of.get(pdb_id)
            if code is None:
                code = len(pdb_table)
                pdb_code_of[pdb_id] = code
                pdb_table.append(pdb_id)
            pdb_codes.append(code)
            chains.append(ord(chain))
            target_len.append(seq_len)
        row_offsets.append(len(pdb_codes))

    pdb_blob, pdb_offsets = build_sequence_store(pdb_table)
    prot_blob, prot_offsets = build_sequence_store(list(prot_ids))

    return HitStore(
        prot_blob = prot_blob,
        prot_offsets = prot_offsets,
        row_offsets = np.frombuffer(row_offsets, dtype = np.int64),
        pdb_blob = pdb_blob,
        pdb_offsets = pdb_offsets,
        pdb_codes = np.frombuffer(pdb_codes, dtype = np.int32),
        chains = np.frombuffer(chains, dtype = np.uint8),
        target_len = np.frombuffer(target_len, dtype = np.int32),
        # A self hit's alignment is the whole chain, which is what the tuple version
        # of this put in both fields.
        aln_length = np.frombuffer(target_len, dtype = np.int32),
    )


def _concat_offsets(offset_arrays: list[np.ndarray]) -> np.ndarray:
    """Several stores' CSR offsets laid end to end.

    Every array but the first drops its leading zero and is shifted by the total of
    everything before it.
    """
    parts: list[np.ndarray] = [offset_arrays[0]]
    shift: int = int(offset_arrays[0][-1])
    for offsets in offset_arrays[1:]:
        parts.append(offsets[1:] + shift)
        shift += int(offsets[-1])
    return np.concatenate(parts)


def concat_stores(stores: list[HitStore]) -> HitStore:
    """The stores' protein blocks one after another, sharing one pdb string table.

    Protein ids are not deduplicated here - a protein appearing in two stores keeps two
    blocks, which is what _resolve_prot_collisions exists to fix. Every input must
    still carry aln_length.
    """
    stores = [store for store in stores if store.n_prots > 0]

    if len(stores) == 0:
        return empty_hit_store()
    if len(stores) == 1:
        return stores[0]

    # One global pdb table. Each store's codes are then translated with a single
    # vectorised gather rather than a python lookup per row.
    pdb_table: list[str] = []
    pdb_code_of: dict[str, int] = {}
    remaps: list[np.ndarray] = []

    for store in stores:
        local_table = decode_pdb_table(store)
        remap = np.empty(len(local_table), dtype = np.int32)
        for local_code, pdb_id in enumerate(local_table):
            code = pdb_code_of.get(pdb_id)
            if code is None:
                code = len(pdb_table)
                pdb_code_of[pdb_id] = code
                pdb_table.append(pdb_id)
            remap[local_code] = code
        remaps.append(remap)

    pdb_blob, pdb_offsets = build_sequence_store(pdb_table)

    return HitStore(
        prot_blob = np.concatenate([store.prot_blob for store in stores]),
        prot_offsets = _concat_offsets([store.prot_offsets for store in stores]),
        row_offsets = _concat_offsets([store.row_offsets for store in stores]),
        pdb_blob = pdb_blob,
        pdb_offsets = pdb_offsets,
        pdb_codes = np.concatenate([remap[store.pdb_codes] for remap, store in zip(remaps, stores)]),
        chains = np.concatenate([store.chains for store in stores]),
        target_len = np.concatenate([store.target_len for store in stores]),
        aln_length = np.concatenate([store.aln_length for store in stores]),
    )


def _resolve_prot_collisions(store: HitStore, prot_ids: list[str]) -> tuple[HitStore, int]:
    """Collapse a protein that got more than one block into a single one.

    Chunks are query-aligned (see _query_aligned_chunk_ranges) and parseHits only ever
    opens a new block when the query id on the line changes, so in the ordinary case
    every protein has exactly one block and this returns the store untouched - no per
    hit work at all. The careful merge below is for when that does not hold, e.g.
    mmseqs splitting its own target database and re-emitting a query's hits once per
    split, which neither the chunking nor the parse loop can see.

    Merging is the same rule parseHits applies within a block: the longest alignment
    for a (pdb, chain) wins, and first-appearance order decides the rest.
    """
    first_block_of: dict[str, int] = {}
    blocks_of: dict[str, list[int]] = {}
    order: list[str] = []

    for block, prot_id in enumerate(prot_ids):
        first_block = first_block_of.get(prot_id)
        if first_block is None:
            first_block_of[prot_id] = block
            order.append(prot_id)
        else:
            collided = blocks_of.get(prot_id)
            if collided is None:
                blocks_of[prot_id] = [first_block, block]
            else:
                collided.append(block)

    if len(blocks_of) == 0:
        return store, 0

    row_offsets = store.row_offsets
    pdb_codes = store.pdb_codes
    chains = store.chains
    aln_length = store.aln_length

    index_parts: list[np.ndarray] = []
    new_row_offsets = np.zeros(len(order) + 1, dtype = np.int64)
    n_rows: int = 0

    for prot_index, prot_id in enumerate(order):
        collided = blocks_of.get(prot_id)

        if collided is None:
            block = first_block_of[prot_id]
            low, high = int(row_offsets[block]), int(row_offsets[block + 1])
            index_parts.append(np.arange(low, high, dtype = np.int64))
            n_rows += high - low
        else:
            best: dict[tuple[int, int], tuple[int, int]] = {}
            for block in collided:
                low, high = int(row_offsets[block]), int(row_offsets[block + 1])
                block_pdb = pdb_codes[low:high].tolist()
                block_chain = chains[low:high].tolist()
                block_aln = aln_length[low:high].tolist()
                for offset in range(high - low):
                    entry_key = (block_pdb[offset], block_chain[offset])
                    kept = best.get(entry_key)
                    if kept is None or block_aln[offset] > kept[1]:
                        best[entry_key] = (low + offset, block_aln[offset])
            index_parts.append(np.fromiter(
                (row for row, _ in best.values()), dtype = np.int64, count = len(best)))
            n_rows += len(best)

        new_row_offsets[prot_index + 1] = n_rows

    index = np.concatenate(index_parts) if len(index_parts) > 0 else np.empty(0, dtype = np.int64)
    prot_blob, prot_offsets = build_sequence_store(order)

    return store._replace(
        prot_blob = prot_blob,
        prot_offsets = prot_offsets,
        row_offsets = new_row_offsets,
        pdb_codes = pdb_codes[index],
        chains = chains[index],
        target_len = store.target_len[index],
        aln_length = aln_length[index],
    ), len(blocks_of)


def parseHits(
        temp_outfile: str,
        option_seq_thresh: float,
        small_genes: set[str],
        id_add: None | str,
        min_seq_id: float = 0.0,
        start: int = 0,
        end: int | None = None,
        ) -> tuple[HitStore, list[float]]:

    t0 = time.time()
    t1 = t0

    # One row per stored hit, grown as typed C arrays and reinterpreted as numpy at the
    # end without a copy. The dict of 5-tuples this used to build cost ~400 bytes a hit
    # and had to be pickled object by object on the way back through ray's object store.
    pdb_code_column = array(_CODE_TYPECODE)
    chain_column = array(_CHAIN_TYPECODE)
    target_len_column = array(_CODE_TYPECODE)
    aln_length_column = array(_CODE_TYPECODE)

    # One entry per query protein, in first-appearance order. row_offsets is CSR: it
    # gains the end of a block whenever a block is closed, so it always leads with 0.
    prot_ids: list[str] = []
    row_offsets = array(_OFFSET_TYPECODE, [0])

    # The hit lines stay as bytes; encoding the small-gene set once lets the
    # length filter below run without decoding every line first.
    small_genes_bytes: set[bytes] = {gene_id.encode('ascii') for gene_id in small_genes}
    small_gene_seq_thresh: float = option_seq_thresh * 2.0

    n_lines: int = 0

    # The same PDB codes recur enormously across the file (within one line's
    # comma-separated target list and across many lines), while the universe of
    # distinct values is tiny. They go into a string table, so the bytes->str decode
    # and the id_add formatting happen once per distinct id and every row only carries
    # its int32 code.
    pdb_code_cache: dict[bytes, int] = {}
    pdb_table: list[str] = []

    # The hits of the query currently being read: (pdb_code, chain byte) -> row index,
    # so a second, better alignment for the same chain overwrites its row in place.
    #
    # Query blocks are grouped consecutively in mmseqs output (all hits for a query
    # appear on adjacent lines), so clearing this whenever the query changes means only
    # one query's hits are ever held as python objects, instead of the whole file's - on
    # a whole-PDB chunk that is the difference between a few thousand entries and a few
    # hundred million. A query that does reappear later (mmseqs re-emitting it once per
    # internal target-database split) simply opens a second block with the same id,
    # which is exactly the collision _resolve_prot_collisions exists to resolve.
    cur_rows: dict[tuple[int, int], int] = {}
    last_query_bytes: bytes | None = None

    # This loop already allocated tens of millions of small containers while the process
    # holds a very large live heap (>1TB for whole-PDB inputs). Each generation-2
    # collection triggered along the way has to walk that whole heap, so the collections
    # cost far more here than they can reclaim: nothing built in this loop takes part in
    # reference cycles. The columns themselves are not even gc tracked any more, but the
    # per-query map, the pdb cache and the line-local maps still are.
    gc_was_enabled: bool = gc.isenabled()
    gc.disable()

    try:
        # Streaming beats read().split(b'\n'): the outfile reaches tens of GB for
        # whole-PDB inputs, and materialising every line up front doubles peak RSS
        # on top of an already huge process.
        with open(temp_outfile, 'rb', buffering = 4 * 1024 * 1024) as f:

            if start > 0:
                f.seek(start)

            pos: int = start

            t1 = time.time()

            for line in f:
                pos += len(line)

                words: list[bytes] = line.split()

                if len(words) < 6:
                    if end is not None and pos >= end:
                        break
                    continue

                n_lines += 1

                seq_id: float = 100.0 * float(words[2])

                if seq_id < min_seq_id:
                    if end is not None and pos >= end:
                        break
                    continue

                aln_length: int = int(words[3])

                if aln_length < 50 and (words[0] not in small_genes_bytes or seq_id < small_gene_seq_thresh):
                    if end is not None and pos >= end:
                        break
                    continue

                target_len: int = int(words[4])

                if words[0] != last_query_bytes:
                    if last_query_bytes is not None:
                        prot_ids.append(last_query_bytes.decode('ascii'))
                        row_offsets.append(len(pdb_code_column))
                        cur_rows.clear()
                    last_query_bytes = words[0]

                # A target in the clustered search DBs names several chains, comma
                # separated. Chains belonging to the same PDB entry on one line collapse
                # into a single entry, keyed by the first chain seen -- so only one chain
                # per complex goes on to be aligned.
                line_hits: dict[int, bytes] = {}

                for hit in words[1].split(b','):
                    pdb_bytes, _, chain_bytes = hit.rpartition(b'-')

                    pdb_code = pdb_code_cache.get(pdb_bytes)
                    if pdb_code is None:
                        pdb = pdb_bytes.decode('ascii')
                        if id_add is not None:
                            pdb = f'{id_add}{pdb}'
                        pdb_code = len(pdb_table)
                        pdb_table.append(pdb)
                        pdb_code_cache[pdb_bytes] = pdb_code

                    if pdb_code not in line_hits:
                        line_hits[pdb_code] = chain_bytes

                for pdb_code, chain_bytes in line_hits.items():

                    # Chain names longer than a single character are not handled
                    # downstream, and the column holds the chain as its one ASCII byte.
                    # This used to be tested against the loop variable left over from the
                    # hit loop above, so a single long chain name discarded every hit on
                    # the line instead of just its own.
                    if len(chain_bytes) != 1:
                        continue

                    chain_byte: int = chain_bytes[0]
                    entry_key = (pdb_code, chain_byte)

                    row = cur_rows.get(entry_key)
                    if row is None:
                        cur_rows[entry_key] = len(pdb_code_column)
                        pdb_code_column.append(pdb_code)
                        chain_column.append(chain_byte)
                        target_len_column.append(target_len)
                        aln_length_column.append(aln_length)
                    elif aln_length > aln_length_column[row]:
                        aln_length_column[row] = aln_length
                        target_len_column[row] = target_len

                if end is not None and pos >= end:
                    break

        if last_query_bytes is not None:
            prot_ids.append(last_query_bytes.decode('ascii'))
            row_offsets.append(len(pdb_code_column))

    finally:
        if gc_was_enabled:
            gc.enable()

    n_entries: int = len(pdb_code_column)

    t2 = time.time()

    pdb_blob, pdb_offsets = build_sequence_store(pdb_table)
    prot_blob, prot_offsets = build_sequence_store(prot_ids)

    # np.frombuffer wraps each array's own buffer - no copy, and the array stays alive
    # as the view's exporter. Nothing may append to them from here on: an exported
    # buffer makes array.array unresizable.
    store = HitStore(
        prot_blob = prot_blob,
        prot_offsets = prot_offsets,
        row_offsets = np.frombuffer(row_offsets, dtype = np.int64),
        pdb_blob = pdb_blob,
        pdb_offsets = pdb_offsets,
        pdb_codes = np.frombuffer(pdb_code_column, dtype = np.int32),
        chains = np.frombuffer(chain_column, dtype = np.uint8),
        target_len = np.frombuffer(target_len_column, dtype = np.int32),
        aln_length = np.frombuffer(aln_length_column, dtype = np.int32),
    )

    # cur_rows only ever deduplicates within the run of adjacent lines sharing a query
    # id, so a query mmseqs emitted twice (once per internal target-database split)
    # leaves two blocks with the same id here. Resolving them before returning keeps the
    # "one block per protein" invariant true of every store, including the one the
    # sequential fallback in parseHits_parallel returns without going through
    # merge_parsed_chunks. Costs a dict of protein ids and nothing per hit when, as
    # usual, there are none.
    store, _ = _resolve_prot_collisions(store, prot_ids)

    return store, [t1 - t0, t2 - t1, n_lines, n_entries]


# Below this size, splitting into chunks just adds ray scheduling/serialisation
# overhead without giving the workers enough to do.
MIN_BYTES_PER_PARSE_CHUNK: int = 512 * 1024 * 1024


def _query_aligned_chunk_ranges(temp_outfile: str, n_chunks: int) -> list[tuple[int, int]]:
    """Split temp_outfile into up to n_chunks non-overlapping byte ranges, each
    starting and ending on a query-block boundary (never mid-line, and never in
    the middle of the run of consecutive lines sharing one query id). mmseqs
    output groups every hit for a query onto adjacent lines, so this keeps a
    single query's hits -- and the cross-line max-alignment-length merge state
    that parseHits builds for it -- entirely inside one chunk. That in turn
    means the per-chunk `entries` dicts have disjoint keys and can be combined
    with a plain dict update instead of replaying the (order-dependent) merge.
    """
    file_size = os.path.getsize(temp_outfile)

    if n_chunks <= 1 or file_size == 0:
        return [(0, file_size)]

    boundaries: list[int] = [0]

    with open(temp_outfile, 'rb', buffering = 4 * 1024 * 1024) as f:
        for i in range(1, n_chunks):
            target = file_size * i // n_chunks
            if target <= boundaries[-1]:
                continue

            f.seek(target)
            f.readline()  # discard the partial line straddling target

            if f.tell() >= file_size:
                break

            line = f.readline()
            if not line:
                break

            tokens = line.split(maxsplit=1)
            prev_query = tokens[0] if tokens else None

            boundary_pos = f.tell()
            while True:
                save = f.tell()
                nxt = f.readline()
                if not nxt:
                    boundary_pos = save
                    break
                tokens = nxt.split(maxsplit=1)
                nxt_query = tokens[0] if tokens else None
                if nxt_query != prev_query:
                    boundary_pos = save
                    break
                boundary_pos = f.tell()

            if boundaries[-1] < boundary_pos < file_size:
                boundaries.append(boundary_pos)

    boundaries.append(file_size)

    return [(s, e) for s, e in zip(boundaries, boundaries[1:]) if e > s]


@ray.remote
def _parseHits_remote(
        temp_outfile: str,
        option_seq_thresh: float,
        small_genes: set[str],
        id_add: None | str,
        min_seq_id: float,
        start: int,
        end: int,
        ):
    return parseHits(temp_outfile, option_seq_thresh, small_genes, id_add, min_seq_id = min_seq_id, start = start, end = end)


def merge_parsed_chunks(
        chunk_results: list[tuple[HitStore, list[float]]],
        config: Config,
        ) -> tuple[HitStore, list[float]]:
    """Combine the per-chunk hit stores parseHits_parallel's workers returned.

    Split out from parseHits_parallel so it can be exercised without ray: chunk_results
    is already a plain list here, exactly what ray.get(result_refs) hands back.

    Chunks are query-aligned (see _query_aligned_chunk_ranges), so their protein blocks
    are disjoint in the ordinary case - mmseqs groups every hit for a query onto
    adjacent lines - and concatenating the columns reproduces what a single sequential
    pass over the whole file would give, at memcpy speed instead of a python loop over
    every one of hundreds of millions of hits. That loop, entry by entry, was measured
    taking ~1470-1960s of a ~1520-2005s parse step on whole-PDB inputs, against a
    slowest worker of ~45s - almost all of it spent re-confirming what the
    query-alignment already guaranteed.

    _resolve_prot_collisions is what handles the case where the guarantee does not hold
    - e.g. mmseqs splitting its own target database and re-emitting a query's hits once
    per split. It finds the colliding protein ids from the (far smaller) id list and
    only those get the careful per-hit merge.
    """
    total_n_lines = 0
    total_n_entries = 0
    max_parse_t = 0.0
    total_open_t = 0.0

    # This merge is serial in the main process. Unlike the version that walked every
    # hit it no longer builds python containers per hit at all, but the pdb table
    # remapping and the collision bookkeeping still allocate while the driver holds a
    # very large live heap, and every gen-2 collection walks that heap and reclaims
    # nothing. Same reasoning as the guard in parseHits itself.
    gc_was_enabled: bool = gc.isenabled()
    gc.disable()

    try:
        t_merge0 = time.time()

        stores: list[HitStore] = []
        for chunk_store, chunk_times in chunk_results:
            stores.append(chunk_store)

            open_t, parse_t, n_lines, n_entries = chunk_times
            total_open_t += open_t
            max_parse_t = max(max_parse_t, parse_t)
            total_n_lines += n_lines
            total_n_entries += n_entries

        combined = concat_stores(stores)
        del stores

        prot_ids = decode_prot_ids(combined)
        total_chunk_prot_ids = len(prot_ids)

        combined, total_collided_prot_ids = _resolve_prot_collisions(combined, prot_ids)
        del prot_ids

        t_merge = time.time() - t_merge0
    finally:
        if gc_was_enabled:
            gc.enable()

    if config.verbosity >= 2:
        # merge_t and the collision count are what tells this step apart from ray.get()
        # itself (timed by the caller, since only it holds the wait) and from
        # open_t/parse_t, which only ever described the workers' own parsing.
        config.logger.info(
            f'parseHits_parallel merge: merge_t={t_merge:.1f}s '
            f'{total_collided_prot_ids=:_} of {total_chunk_prot_ids=:_} chunk-level prot_ids collided, '
            f'{combined.n_prots=:_} distinct, {combined.n_hits=:_} hits'
        )

    return combined, [total_open_t, max_parse_t, total_n_lines, total_n_entries]


def parseHits_parallel(
        temp_outfile: str,
        option_seq_thresh: float,
        small_genes: set[str],
        id_add: None | str,
        config: Config,
        min_seq_id: float = 0.0,
        ) -> tuple[HitStore, list[float]]:

    file_size = os.path.getsize(temp_outfile)
    n_chunks = min(config.proc_n, max(1, file_size // MIN_BYTES_PER_PARSE_CHUNK))

    if n_chunks <= 1 or not ray.is_initialized():
        return parseHits(temp_outfile, option_seq_thresh, small_genes, id_add, min_seq_id = min_seq_id)

    ranges = _query_aligned_chunk_ranges(temp_outfile, n_chunks)

    if len(ranges) <= 1:
        return parseHits(temp_outfile, option_seq_thresh, small_genes, id_add, min_seq_id = min_seq_id)

    small_genes_ref = ray.put(small_genes)

    result_refs = [
        _parseHits_remote.remote(temp_outfile, option_seq_thresh, small_genes_ref, id_add, min_seq_id, start, end)
        for start, end in ranges
    ]

    # ray.get() is the slowest worker's task time plus ray's own transfer of the results
    # out of the object store. The nested python dicts/tuples/sets this used to return
    # had to be pickled and rebuilt object by object, which could dominate on its own;
    # HitStore's columns take the same zero-copy path as the numpy-backed transfers
    # elsewhere in this pipeline. Still timed separately from the merge below so a slow
    # part 1.2 can be told apart: worker-side (open_t/parse_t, logged by the caller),
    # ray transfer (here), or the merge (merge_parsed_chunks).
    t_wait0 = time.time()
    chunk_results = ray.get(result_refs)
    t_get = time.time() - t_wait0

    combined, times = merge_parsed_chunks(chunk_results, config)

    if config.verbosity >= 2:
        config.logger.info(f'parseHits_parallel ray.get: get_t={t_get:.1f}s')

    return combined, times


def apply_mmseqs(
        mmseqs_tmp_folder,
        mmseqs2_path,
        temp_fasta,
        search_db,
        config: Config,
        small_proteins: set[str] | None = None,
        number_of_returned_hits = '999999',
        sensitivity = 7.5,
        id_add = None,
        ):

    if small_proteins is None:
        small_proteins: set[str] = set()

    temp_outfile = f'{mmseqs_tmp_folder}/tmp_outfile_{randomString()}.fasta'

    if config.verbosity >= 2:
        t0 = time.time()
        config.logger.info(f'\nApply MMseqs2: {mmseqs2_path}, easy-search, {temp_fasta=}, {search_db=}, {temp_outfile=}, {mmseqs_tmp_folder=} {number_of_returned_hits=}\n')

    out_format_str = 'query,target,fident,alnlen,tlen,qcov'

    debug_store = False

    # Hits below this identity are dropped again after the global alignment (see
    # the option_seq_thresh check in serializedPipeline.paraAlignment), so letting
    # mmseqs discard them keeps them out of the outfile, out of the annotation
    # objects and out of the alignment stage entirely.
    min_seq_id: float = getattr(config, 'mmseqs_min_seq_id', 0.0)

    cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder,
            '--format-output', out_format_str,
            '--max-seqs', number_of_returned_hits,
            '-s', f'{sensitivity}',
            '--split-memory-limit', '%1.0fG' % (0.5 * config.gigs_of_ram)]

    # Alignments shorter than 50 residues are filtered in parseHits anyway, but
    # small proteins are exempted there (at a doubled identity threshold), so the
    # filter can only be pushed down to mmseqs when there are none.
    if len(small_proteins) == 0:
        cmds += ['--min-aln-len', '50']

    if min_seq_id > 0.0:
        cmds += ['--min-seq-id', f'{min_seq_id}']

    if config.verbosity >= 5:
        p = subprocess.Popen(cmds)
        p.wait()
    else:
        FNULL = open(os.devnull, 'w')
        p = subprocess.Popen(cmds, stdout=FNULL)
        p.wait()

        if not os.path.exists(temp_outfile):
            #Something went wrong, run again in not silenced mode
            config.errorlog.add_error(f'Error in MMSEQS2 call, see command line prints for possible bug log.\nCommand that called mmseqs:\n{cmds}\nDeletion of temp file: {temp_fasta} was suppressed')
            debug_store = True
            p = subprocess.Popen(cmds)


    if config.verbosity >= 2:
        t1 = time.time()
        config.logger.info(f'apply_mmseqs part 1: {t1-t0}')

    if config.verbosity >= 3:
        config.logger.info(f'MMseqs2 returned: {temp_outfile}')

    if config.verbosity >= 8:
        f = open(temp_outfile, 'r')
        config.logger.info(f.read())
        f.close()

    store, parse_times = parseHits_parallel(temp_outfile, config.option_seq_thresh, small_proteins, id_add, config, min_seq_id = 100.0 * min_seq_id)

    if config.verbosity >= 2:
        t2 = time.time()
        open_t, parse_t, n_lines, n_entries = parse_times
        config.logger.info(f'apply_mmseqs part 1.2: {t2-t1}, {open_t=} {parse_t=} {n_lines=:_} {n_entries=:_} {min_seq_id=}')

    if config.verbosity >= 3:
        config.logger.info(f'MMseqs2 results parsed: size of hit map: {store.n_prots} proteins, {store.n_hits} hits')

    os.remove(temp_outfile)

    return store, debug_store


def wipe_folder(config, folder_path):
    for fn in os.listdir(folder_path):
        subfolder_path = '%s/%s' % (folder_path, fn)
        if os.path.exists(subfolder_path):
            if os.path.getmtime(subfolder_path) > config.prog_start_time:
                try:
                    shutil.rmtree(subfolder_path)
                except:
                    if config.verbosity >= 6:
                        config.errorlog.add_warning('Tmp folder wipe failed for: %s' % subfolder_path)


# called by serializePipeline
def search(proteins: Proteins, config: Config) -> list[tuple[HitStore, bool]] | None:

    # accepts one more argument, custom db

    mmseqs_tmp_folder = f'{config.mmseqs_tmp_folder}/mmseqs_tmp'
    if not os.path.exists(mmseqs_tmp_folder):
        os.mkdir(mmseqs_tmp_folder)

    mmseqs2_path = config.mmseqs2_path

    # custom db mmseqs search processed here
    if config.custom_db_path:
        createCustomDb.main(config)

        mmseqs_custom_db_path = config.custom_db_path + "/custom_db"

        if not os.path.exists(mmseqs_custom_db_path):
            cmds = ["mkdir", "custom_db"]
            p = subprocess.Popen(cmds, cwd=config.custom_db_path)

        cmds = [mmseqs2_path, 'createdb', config.custom_db_fasta_path, "./structman_custom_db"]

        p = subprocess.Popen(cmds, cwd=mmseqs_custom_db_path)
        p.wait()

    max_sensi = 7.5

    search_dbs = [config.mmseqs2_db_path]
    numbers_of_returned_hits = [None]
    is_model_dbs = [False]
    sensitivities = [(max_sensi,)]
    id_adds = [None]

    # checks if custom db path is provided in config file
    if (config.custom_db_path is not None) and (config.custom_db_path != ''):
        search_dbs.append(f"{mmseqs_custom_db_path}/structman_custom_db")
        numbers_of_returned_hits.append(None)
        is_model_dbs.append(True)
        sensitivities.append((max_sensi,))
        id_adds.append(None)

    if config.model_db_active:
        search_dbs.append(config.mmseqs2_model_db_path)
        numbers_of_returned_hits.append('1')
        is_model_dbs.append(True)
        sensitivities.append((1.0, max_sensi))
        id_adds.append(None)

    if (config.mmseqs2_complex_model_db_path is not None) and (config.mmseqs2_complex_model_db_path != ''):
        search_dbs.append(config.mmseqs2_complex_model_db_path)
        numbers_of_returned_hits.append(None)
        is_model_dbs.append(True)
        sensitivities.append((max_sensi,))
        id_adds.append(None)

    if (config.afdb_hetero_mmseqs_path is not None) and (config.afdb_hetero_mmseqs_path != ''):
        search_dbs.append(config.afdb_hetero_mmseqs_path)
        numbers_of_returned_hits.append(None)
        is_model_dbs.append(True)
        sensitivities.append((max_sensi,))
        id_adds.append('afdb_hetero_')

    if (config.afdb_homo_mmseqs_path is not None) and (config.afdb_homo_mmseqs_path != ''):
        search_dbs.append(config.afdb_homo_mmseqs_path)
        numbers_of_returned_hits.append('1')
        is_model_dbs.append(True)
        sensitivities.append((1.0, max_sensi))
        id_adds.append('afdb_homo_')


    small_proteins: set[str] = set()

    u_acs = proteins.get_protein_ids()
    for u_ac in u_acs:
        try:
            if len(proteins.get_sequence(u_ac)) < 100:
                small_proteins.add(u_ac)
        except:
            continue

    t0 = time.time()

    wipe_folder(config, mmseqs_tmp_folder)

    temp_fasta = f'{mmseqs_tmp_folder}/tmp_{randomString()}.fasta'
    to_fasta_out, num_of_prots = geneSeqMapToFasta(proteins, temp_fasta, config)

    if isinstance(to_fasta_out, str):
        config.errorlog.add_warning(f'{to_fasta_out} , mmseqs2 skipped, {list(u_acs)[:10]!s}')
        return None
    if not to_fasta_out:
        # All proteins are stored, no need for mmseqs
        if config.verbosity >= 3:
            config.logger.info('No need for sequence similarity search, all proteins are stored already.')
        return None

    default_number_of_returned_hits = '999999'

    search_results: list[tuple[HitStore, bool]] = []

    t1 = time.time()

    for pos, search_db in enumerate(search_dbs):

        if numbers_of_returned_hits[pos] is None:
            number_of_returned_hits = default_number_of_returned_hits
        else:
            number_of_returned_hits = numbers_of_returned_hits[pos]

        is_model_db = is_model_dbs[pos]

        internal_temp_fasta = temp_fasta

        id_add = id_adds[pos]

        for sensi in sensitivities[pos]:
            t10 = time.time()

            store: HitStore

            store, debug_store = apply_mmseqs(
                mmseqs_tmp_folder, mmseqs2_path, internal_temp_fasta, search_db, config, 
                small_proteins = small_proteins, number_of_returned_hits = number_of_returned_hits, 
                sensitivity = sensi, id_add = id_add)

            t11 = time.time()
            if config.verbosity >= 2:
                config.logger.info(f'MMseqs2 {search_db} Part 2.1: {t11 - t10}')

            # The store answers "which proteins were hit" by decoding its id blob, which
            # is one pass over a few thousand ids; the three membership tests below then
            # read this set rather than the store.
            hit_prot_ids: set[str] = set(decode_prot_ids(store))

            self_hit_ids: list[str] = []
            self_hit_lens: list[int] = []
            for u_ac in small_proteins:
                if u_ac.count(':') == 1 and u_ac not in hit_prot_ids:
                    # A chain mmseqs returned nothing for still gets its own structure.
                    self_hit_ids.append(u_ac)
                    self_hit_lens.append(len(proteins.get_sequence(u_ac)))

            if len(self_hit_ids) > 0:
                store = concat_stores([store, build_self_hit_store(self_hit_ids, self_hit_lens)])
                hit_prot_ids.update(self_hit_ids)

            t12 = time.time()
            if config.verbosity >= 2:
                config.logger.info(f'MMseqs2 {search_db} Part 2.2: {t12 - t11}')

            wipe_folder(config, mmseqs_tmp_folder)

            t13 = time.time()
            if config.verbosity >= 2:
                config.logger.info(f'MMseqs2 {search_db} Part 2.3: {t13 - t12}')

            # aln_length only ever existed for the chunk merge's tie-break, and this
            # store is final now, so it does not need to stay resident for the run.
            search_results.append((drop_aln_length(store), is_model_db))

            # The undropped store is the only thing still holding aln_length - hundreds
            # of millions of int32s on a whole-PDB search - and nothing reads it again.
            del store
    
            t14 = time.time()
            if config.verbosity >= 2:
                config.logger.info(f'MMseqs2 {search_db} Part 2.4: {t14 - t13}')

            if len(hit_prot_ids) >= num_of_prots:
                break

            not_hit_yet = []
            for prot_id in proteins.get_protein_ids():
                if prot_id not in hit_prot_ids:
                    not_hit_yet.append(prot_id)

            internal_temp_fasta = f'{mmseqs_tmp_folder}/tmp_{randomString()}.fasta'
            to_fasta_out, num_of_prots = geneSeqMapToFasta(proteins, internal_temp_fasta, config, custom_prot_ids=not_hit_yet)



    if (not debug_store) and not (config.verbosity >= 6):
        os.remove(temp_fasta)

    t2 = time.time()

    if config.verbosity >= 2:
        config.logger.info(f"MMseqs2 Part 1: {t1 - t0}")
        config.logger.info(f"MMseqs2 Part 2: {t2 - t1}")

    return search_results


def fake_search(proteins: Proteins, config: Config) -> list[tuple[HitStore, bool]]:
    #This gets only called in fast_pdb_annotation, all protein IDs are PDB Ids
    structure_ids = proteins.get_protein_ids()

    # dict.fromkeys rather than a set: the map this used to build kept the input's
    # order, and that order decides the order Structures are created in downstream.
    unique_ids: list[str] = list(dict.fromkeys(structure_ids))

    store = build_self_hit_store(unique_ids, [400] * len(unique_ids))

    return [(drop_aln_length(store), False)]