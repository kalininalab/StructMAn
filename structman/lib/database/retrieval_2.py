import numpy as np
import ray
import time
from structman.lib.database.database_core_functions import binningSelect, select
from structman.base_utils.base_utils import unpack, add_to_times, aggregate_times, reset_logger_for_remotes
from structman.base_utils.config_class import Config
from structman.base_utils.ray_utils import ray_init, parallel_db_access

# Starting capacity of a shard's backmap buffer, in int32 entries (4 MiB). Big enough
# that a normal shard doubles a handful of times, small enough not to matter when a
# shard turns out to be empty.
_INITIAL_BACKMAP_CAPACITY = 1024 * 1024


@ray.remote(max_calls=1)
def remote_getBackmaps(*args, **kwargs):
    # Thin wrapper so the body can also be called directly, without ray, by the
    # serial branch of its dispatcher (see ray_utils.parallel_db_access).
    return get_backmaps_shard(*args, **kwargs)


def get_backmaps_shard(
        db_nr: int,
        config: Config,
        structure_db_ids: list[int],
        database_ids: set[int] | None,
        ) -> tuple[np.ndarray[int], np.ndarray[int], list]:
    """
    Selects the alignments of one shard and returns them as flat numpy arrays.

    The backmaps are unpacked here rather than in the parent: it spreads the
    unpacking over the workers, and it keeps the packed blobs out of the result
    that has to travel back through the object store.

    Returns the shard's concatenated backmaps together with one row per
    (structure, protein) pair holding
    (structure_db_id, prot_db_id, backmap_left, backmap_right), where the backmap
    bounds index into that concatenated array. Those bounds are local to the
    shard; getBackmaps shifts them when it joins the shards together.
    """

    times = []
    t0 = time.time()

    table = f'Alignment_{db_nr}'
    columns = ['Structure', 'Protein', 'Backmap']
    results, sel_times = binningSelect(structure_db_ids, columns, table, config, locked=True, return_run_times=True, full_select=True)

    times.append(sel_times)
    t0 = add_to_times(times, t0, 'binningSelect Alignments')

    isolate_ppi: bool = config.isolate_ppi

    entry_list: list[tuple[int, int, int, int]] = []
    bm_right: int = 0

    # One buffer that grows, rather than a list of unpacked backmaps plus the
    # np.concatenate of it. Holding all three - the packed rows, the unpacked list and
    # the concatenation - put three copies of a shard's backmaps in a worker at once,
    # and on a whole-PDB input the shards are asked about every structure. The rows are
    # dropped as they are consumed, so the packed blobs do not survive the loop either.
    # Backmaps hold position db ids, which stay well below 2**31; the dtype is forced
    # rather than inherited so that rows written before db_backmap was narrowed are
    # brought down on read too.
    concat_backmaps: np.ndarray[np.int32] = np.empty(_INITIAL_BACKMAP_CAPACITY, dtype=np.int32)

    for index in range(len(results)):
        row = results[index]
        results[index] = None

        prot_db_id: int = row[1]

        if isolate_ppi and prot_db_id not in database_ids:
            continue

        backmap: np.ndarray[np.int32] = unpack(row[2])
        if backmap is None:
            # NULL: the annotation takes no part in the interface calculation - its
            # complex has a single protein chain, or select_interface_annotations
            # left it out (database.prune_interface_backmaps) - or the database was
            # written with --fast_pdb_annotation. Not listed at all, which the
            # interface calculation treats exactly as it did the empty backmap this
            # used to become: a target protein without an entry is skipped before its
            # keys are walked rather than failing the bounds check on every one of
            # them, and a missing partner adds no slot instead of one nothing reads.
            continue

        length = len(backmap)
        if bm_right + length > concat_backmaps.shape[0]:
            # resize rather than a new array plus a copy: it reallocates in place where
            # it can, and never holds old and new at Python level.
            concat_backmaps.resize(max(2 * concat_backmaps.shape[0], bm_right + length),
                                   refcheck=False)
        if length > 0:
            concat_backmaps[bm_right:bm_right + length] = backmap

        entry_list.append((row[0], prot_db_id, bm_right, bm_right + length))
        bm_right += length

    del results

    # Shrinks in place, so the slack from the last doubling is given back here rather
    # than travelling to the parent through the object store.
    concat_backmaps.resize(bm_right, refcheck=False)

    t0 = add_to_times(times, t0, 'unpack backmaps loop')

    # Deliberately int64: columns 2 and 3 are offsets into the concatenation
    # above, whose total length is not bounded the way a single db id is. This
    # array has one row per alignment, so it is small next to the backmaps.
    entries: np.ndarray[np.int64] = np.array(entry_list, dtype=np.int64).reshape(-1, 4)

    t0 = add_to_times(times, t0, 'concat shard arrays')

    return concat_backmaps, entries, times

def getBackmaps(
        structure_db_ids: list[dict[int, tuple[str, str]]],
        config: Config,
        max_struct_db_id: int,
        database_ids: set[int]
        ) -> (
            tuple[
                tuple[
                    np.ndarray[int],
                    np.ndarray[int],
                    np.ndarray[int]
                    ],
                list
                ]
            ):
    
    times = []
    ta = time.time()
    ta_0 = ta

    parallel = parallel_db_access(config)

    if parallel:
        # See prepare_struct_base_info_store: ray can be down when the output stage
        # calls in, and ray.put would otherwise auto-init an unconfigured cluster.
        ray_init(config, stop_existing = False)
        # One entry per protein, so worth putting into the object store once rather
        # than serialising it into every one of the 20 remote calls.
        database_ids_store = ray.put(database_ids) if config.isolate_ppi else None
    else:
        database_ids_store = database_ids if config.isolate_ppi else None

    # The Alignment tables are sharded by protein db id, but the caller bins these
    # structure ids by the db id of the *session* protein that annotated them. The
    # other alignments of a structure - an interaction partner that a previous
    # session stored - live in the shard of that partner's own protein id, which is
    # only visited when it happens to agree modulo 20 with a session protein. That is
    # why an interaction with an already stored protein was invisible: its backmap
    # was in an Alignment_N table nobody asked. When interfaces are being computed
    # every shard is asked about every structure instead. Structure is indexed on all
    # 20 tables, so that is an index lookup per shard rather than a scan, and a
    # protein lives in exactly one shard so no row can come back twice.
    if config.compute_ppi:
        all_structure_db_ids: set[int] = set()
        for id_bin in structure_db_ids:
            all_structure_db_ids.update(id_bin.keys())
        # Put once rather than serialising the whole id list into all 20 task specs.
        if len(all_structure_db_ids) == 0:
            shared_ids = None
        elif parallel:
            shared_ids = ray.put(list(all_structure_db_ids))
        else:
            shared_ids = list(all_structure_db_ids)
        shard_slices = [shared_ids] * 20
    else:
        shard_slices = [(list(id_bin.keys()) or None) for id_bin in structure_db_ids]

    ta = add_to_times(times, ta, 'put structure_db_ids')

    shards = [db_nr for db_nr in range(20) if shard_slices[db_nr] is not None]

    if parallel:
        process_ids = [remote_getBackmaps.remote(db_nr, config, shard_slices[db_nr], database_ids_store) for db_nr in shards]
        remote_results = ray.get(process_ids)
    else:
        remote_results = [get_backmaps_shard(db_nr, config, shard_slices[db_nr], database_ids_store) for db_nr in shards]

    ta = add_to_times(times, ta, 'remote getBackmaps')

    concat_parts: list[np.ndarray[int]] = []
    entry_parts: list[np.ndarray[int]] = []
    bm_offset: int = 0

    total_loop_times = []
    total_alignments = 0
    total_backmap_bytes = 0
    for shard_index, (concat_backmaps, entries, remote_times) in enumerate(remote_results):
        total_loop_times = aggregate_times(total_loop_times, remote_times)

        total_alignments += len(entries)
        total_backmap_bytes += concat_backmaps.nbytes

        # No run has ever reached the end of getBackmaps on a whole-PDB input, so the
        # volume it moves has only ever been an estimate. One line per shard, in the
        # format the getStoredResidues packages use.
        if config.verbosity >= 2:
            config.logger.info(
                f'getBackmaps shard {shard_index + 1}/{len(remote_results)}: '
                f'{len(entries)} alignments, backmaps {concat_backmaps.nbytes / 1024**2:.0f} MiB '
                f'| total {total_alignments} alignments, {total_backmap_bytes / 1024**3:.2f} GiB, '
                f'elapsed {time.time() - ta_0:.0f}s'
            )

        if len(entries) == 0:
            continue

        # Each shard numbers its backmap bounds from zero; shift them onto the
        # concatenation of all shards.
        if bm_offset > 0:
            entries = entries.copy()
            entries[:, 2:4] += bm_offset

        entry_parts.append(entries)
        concat_parts.append(concat_backmaps)
        bm_offset += len(concat_backmaps)

    times.append(total_loop_times)

    ta = add_to_times(times, ta, 'proccess remote results')

    # Addressed by structure db id, so mostly -1 padding. Its two columns are row
    # bounds into the alignment entries, one row per (structure, protein) pair,
    # which is nowhere near int32's range.
    structure_db_id_vec: np.ndarray[np.int32] = np.full((max_struct_db_id+100, 2), -1, dtype=np.int32)

    if len(entry_parts) > 0:
        entries = np.concatenate(entry_parts) if len(entry_parts) > 1 else entry_parts[0]
        concat_backmaps_list = np.concatenate(concat_parts) if len(concat_parts) > 1 else concat_parts[0]

        # structure_db_id_vec addresses each structure's proteins as a single
        # slice of prot_db_id_vec, so those rows have to sit together. The
        # Alignment tables are sharded by protein id, so one structure can turn
        # up in several shards and the grouping has to span all of them. The sort
        # is stable, which keeps the ordering the serial version produced.
        entries = entries[np.argsort(entries[:, 0], kind='stable')]

        struct_db_ids, starts, counts = np.unique(entries[:, 0], return_index=True, return_counts=True)
        structure_db_id_vec[struct_db_ids, 0] = starts
        structure_db_id_vec[struct_db_ids, 1] = starts + counts

        # Carries the backmap offsets, so int64 like 'entries' itself.
        prot_db_id_vec: np.ndarray[np.int64] = entries[:, 1:4]
    else:
        concat_backmaps_list = np.empty(0, dtype=np.int32)
        prot_db_id_vec = np.empty((0, 3), dtype=np.int64)

    ta = add_to_times(times, ta, 'proccess backmaps')

    return (concat_backmaps_list, prot_db_id_vec, structure_db_id_vec), times

def retrieve_stored_ligands(config):
    stored_ligands = {}
    results = select(config, ['Ligand_Id', 'Name'], 'Ligand')
    for row in results:
        ligand_name = row[1]
        stored_ligands[ligand_name] = row[0]
    return stored_ligands

@ray.remote(max_calls=1)
def prepare_struct_base_info_store_subprocess(*args, **kwargs):
    # Thin wrapper so the body can also be called directly, without ray, by the
    # serial branch of its dispatcher (see ray_utils.parallel_db_access).
    return prepare_struct_base_info_shard(*args, **kwargs)


def prepare_struct_base_info_shard(
    config, structure_db_ids, db_nr
    ) -> (
        tuple[
            np.ndarray[int],
            ray.ObjectRef,
            ray.ObjectRef,
            ray.ObjectRef,
            int,
            list
        ]
    ):
    #reset_logger_for_remotes(config)
    outer_times = []
    to = time.time()

    columns = ['Structure_Id', 'PDB', 'Chain', 'Parse_Data', 'Resolution']
    table = f'Structure_{db_nr}'

    results, sel_times = binningSelect(structure_db_ids, columns, table, config, return_run_times=True)#, full_select=True)

    outer_times.append(sel_times)

    to = add_to_times(outer_times, to, 'binningSelect Structures')

    db_id_vector: np.ndarray[int] = np.empty((len(results), 4), dtype=int)
    id_vector: list[str] = []
    bytes_matrix: list[bytes] = []
    resolution_vec: np.ndarray[float] = np.empty(len(results), dtype=float)

    for row_nr, row in enumerate(results):

        parse_data: bytes = row[3]
        struct_db_id: int = row[0]
        struct_id: str = row[1]
        chain: str = row[2]
        resolution: float = row[4]

        internal_struct_db_id: int = 100 * struct_db_id + db_nr
        fused_id: str = f'{struct_id}{chain}'

        trailing_count = len(parse_data) - len(parse_data.rstrip(b'\x00'))
        db_id_vector[row_nr] = (internal_struct_db_id, db_nr, row_nr, trailing_count)
        id_vector.append(fused_id)
        
        bytes_matrix.append(parse_data)
        resolution_vec[row_nr] = resolution

        #config.logger.info(f'{internal_struct_db_id=} {fused_id=} {parse_data=}')

    to = add_to_times(outer_times, to, 'inner loop')

    # Handed back as ObjectRefs when running under ray, so these arrays stay in the
    # object store instead of travelling back through the driver. Running serially
    # there is nowhere to put them, so the arrays themselves are returned and the
    # consumer unwraps either form with ray_utils.deref.
    if parallel_db_access(config):
        id_vector_store = ray.put(np.array(id_vector, dtype='S'))
        bytes_matrix_store = ray.put(np.array(bytes_matrix, dtype='S'))
        resolution_store = ray.put(resolution_vec)
    else:
        id_vector_store = np.array(id_vector, dtype='S')
        bytes_matrix_store = np.array(bytes_matrix, dtype='S')
        resolution_store = resolution_vec

    to = add_to_times(outer_times, to, 'ray.puts')

    return db_id_vector, id_vector_store, bytes_matrix_store, resolution_store, db_nr, outer_times


def prepare_struct_base_info_store(
        structure_db_id_mask: np.ndarray,
        config
    ):
    """Base infos of the structures the alignment needs.

    structure_db_id_mask is a bool array indexed by internal structure db id
    (100 * Structure_Id + db_nr): a set of ints the size of the structure count used
    to be built one hit at a time. Entries below 100 stand for structures without a
    database id and are skipped, as the 0 of Structure.get_database_id always was.
    """
    #columns = ['Structure_Id', 'PDB', 'Chain', 'Parse_Data', 'Resolution']

    times = []
    t0 = time.time()

    internal_struct_db_ids = np.flatnonzero(structure_db_id_mask)
    db_nrs = internal_struct_db_ids % 100
    struct_db_ids = internal_struct_db_ids // 100
    has_db_id = struct_db_ids != 0

    # Sorted per shard, which is what binningSelect's binning wants anyway. Plain ints
    # rather than numpy scalars, since they end up as SQL parameters.
    shards = []
    for db_nr in range(20):
        ids = struct_db_ids[has_db_id & (db_nrs == db_nr)]
        if len(ids) > 0:
            shards.append((db_nr, ids.tolist()))

    if parallel_db_access(config):
        # ray may well be down here: the pipeline shuts it down when it finishes and is
        # skipped entirely when the session already exists, so the output stage arrives
        # with no instance. Without this, ray.put falls into ray's auto_init_hook and a
        # bare ray.init() starts an unconfigured cluster that times out under load.
        ray_init(config, stop_existing = False)
        conf_store = ray.put(config)
        process_ids = [prepare_struct_base_info_store_subprocess.remote(conf_store, ids, db_nr) for db_nr, ids in shards]
    else:
        process_ids = None
        subproc_results = [prepare_struct_base_info_shard(config, ids, db_nr) for db_nr, ids in shards]

    t0 = add_to_times(times, t0, 'prepare subprocesses')

    if process_ids is not None:
        subproc_results = ray.get(process_ids)

    base_info_stores: list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]] = [None for _ in range(20)]
    concat_db_id_vectors = []

    total_sub_proc_times = []
    for db_id_vector, id_vector_store, bytes_matrix_store, resolution_store, db_nr, sub_proc_times in subproc_results:
        concat_db_id_vectors.append(db_id_vector)
        base_info_stores[db_nr] = (id_vector_store, bytes_matrix_store, resolution_store)

        total_sub_proc_times = aggregate_times(total_sub_proc_times, sub_proc_times)

    times.append(total_sub_proc_times)

    if len(concat_db_id_vectors) > 0:
        concat_db_id_vectors: np.ndarray[int] = np.concatenate(concat_db_id_vectors, dtype=int)
    else:
        concat_db_id_vectors: np.ndarray[int] = np.empty((0, 4), dtype=int)

    t0 = add_to_times(times, t0, 'subprocesses')

    return concat_db_id_vectors, base_info_stores, times


def get_struct_base_infos(
        structure_db_ids: list[int],
        config,
        locked = True
    ):
    columns = ['Structure_Id', 'PDB', 'Chain', 'Parse_Data', 'Resolution']

    struct_db_bins = [[] for _ in range(20)]

    for internal_struct_db_id in structure_db_ids:
        db_nr = internal_struct_db_id % 100
        struct_db_id = internal_struct_db_id//100

        struct_db_bins[db_nr].append(struct_db_id)

    total_outer_times = []

    structure_info_dict = {}
    for db_nr, structure_db_ids in enumerate(struct_db_bins):
        if len(structure_db_ids) == 0:
            continue
        outer_times = []
        to = time.time()
        table = f'Structure_{db_nr}'

        results, sel_times = binningSelect(structure_db_ids, columns, table, config, locked=locked, return_run_times=True)#, full_select=True)

        outer_times.append(sel_times)

        to = add_to_times(outer_times, to, 'binningSelect Structures')

        total_inner_times = []
        for row in results:
            inner_times = []
            ti = time.time()

            parse_data = unpack(row[3])

            ti = add_to_times(inner_times, ti, 'unpack parse_data')

            struct_db_id = row[0]
            struct_id = row[1]
            chain = row[2]
            
            resolution = row[4]
            internal_struct_db_id = 100 * struct_db_id + db_nr
            structure_info_dict[internal_struct_db_id] = (struct_id, chain, parse_data, resolution)

            ti = add_to_times(inner_times, ti, 'build dict')

            total_inner_times = aggregate_times(total_inner_times, inner_times)
        outer_times.append(total_inner_times)

        to = add_to_times(outer_times, to, 'inner loop')

        total_outer_times = aggregate_times(total_outer_times, outer_times)

    return structure_info_dict, total_outer_times
