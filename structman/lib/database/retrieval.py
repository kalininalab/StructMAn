import shutil
import time
import ray
import numpy as np
from structman.base_utils.base_utils import unpack, aggregate_times, reset_logger_for_remotes, add_to_times
from structman.base_utils.config_class import Config
from structman.base_utils.ray_utils import ray_init, check_store_budget, object_store_capacity_bytes
from structman.lib.sdsc.protein import Proteins
from structman.lib.database.database_core_functions import binningSelect, select

@ray.remote(max_calls=1)
def remote_select_package_from_db(*args, **kwargs):
    # Thin wrapper so the body can also be called directly, without ray, by the
    # serial branch of its dispatcher (see ray_utils.parallel_db_access).
    return select_package_from_db_shard(*args, **kwargs)


def select_package_from_db_shard(
    id_package: list[int],
    config: Config,
    columns: list[str],
    db_nr: int,
    store_index: int,
    skip_interfaces: bool,
    pr_is_none: bool = False,
    ) -> (
            tuple[list, int, int] |
            tuple[
                dict,
                np.ndarray[int],
                tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef],
                int
                ]
    ):
    # One task handles exactly one bounded package of structure ids. It used to
    # loop over sub-packages of the whole Structure_{db_nr} bin while
    # accumulating every result, so a task's heap grew with the bin instead of
    # with the package size, and the 20 bin tasks ran unconditionally in
    # parallel. On whole-PDB sized inputs that put the sum of all retrieved
    # residue data into worker heaps at once. Packaging and the parallelism
    # limit now live in getStoredResidues, where they can be capped.
    reset_logger_for_remotes(config)

    table = f'Structure_{db_nr}'
    t0 = time.time()
    results = binningSelect(id_package, columns, table, config, picky_select=True)
    t1 = time.time()
    if config.verbosity >= 4:
        config.logger.info(f'Time for structure select: {t1-t0} {db_nr=} {store_index=} {len(results)=}')

    if not pr_is_none:
        return results, db_nr, store_index

    # The blobs go into the object store still compressed, and aggregation decodes
    # only the structures of the protein it is working on (see
    # aggregation.build_local_stores). Unpacking here meant 780 bytes per residue for
    # every structure of the run held in the store at once - hundreds of GiB for a
    # whole-PDB input, which is what made it spill - against roughly a sixth of that
    # compressed. It also takes the decoding and concatenation out of this task.
    res_blobs: list[bytes | None] = []
    iface_blobs: list[bytes | None] = []
    resolutions: list[np.float32] = []
    packed_structures: list[tuple[int, int, int, int, int, int, int, int]] = []

    for row_index in range(len(results)):
        row = results[row_index]
        results[row_index] = None

        residue_blob = row[1]
        if residue_blob is None:
            # No residue data: the structure keeps store_index -1, as it did when a
            # blob unpacked to nothing.
            continue

        internal_structure_db_id = 100 * row[0] + db_nr
        blob_index = len(res_blobs)
        res_blobs.append(residue_blob)

        if skip_interfaces:
            compl_db_id = -1
            resolution = 2.0
            iface_blobs.append(None)
        else:
            compl_db_id = row[2]
            resolution = row[3]
            if resolution is None or np.isnan(resolution):
                resolution = 2.0
            iface_blobs.append(row[4])

        resolutions.append(np.float32(resolution))

        # Columns of total_structures: store_index, blob_index, three unused offset
        # slots, compl_db_id, resolution_index. The residue and interface blobs of a
        # structure share its blob_index, and its resolution sits at the same index.
        packed_structures.append((internal_structure_db_id, store_index, blob_index, -1, -1, -1, compl_db_id, blob_index))

    del results

    residue_bytes = sum(len(blob) for blob in res_blobs)
    interface_bytes = sum(len(blob) for blob in iface_blobs if blob is not None)

    res_buffer, res_offsets = concatenate_blobs(res_blobs)
    iface_buffer, iface_offsets = concatenate_blobs(iface_blobs)

    t2 = time.time()

    res_buffer = ray.put(res_buffer)
    res_offsets = ray.put(res_offsets)
    iface_buffer = ray.put(iface_buffer)
    iface_offsets = ray.put(iface_offsets)
    resolutions = ray.put(np.array(resolutions, dtype=np.float32))

    t3 = time.time()

    packed_structures = np.array(packed_structures, dtype=np.int32).reshape(-1, 8)

    package_stats = {
        'n_structures': int(packed_structures.shape[0]),
        'select_seconds': t1 - t0,
        'pack_seconds': t2 - t1,
        'put_seconds': t3 - t2,
        'residue_bytes': residue_bytes,
        'interface_bytes': interface_bytes,
    }

    return (package_stats, packed_structures,
            (res_buffer, res_offsets, iface_buffer, iface_offsets, resolutions),
            store_index)


def concatenate_blobs(blobs: list[bytes | None]) -> tuple[np.ndarray, np.ndarray]:
    """One uint8 buffer holding every blob, and int64 offsets with len(blobs) + 1 entries.

    Blob i is buffer[offsets[i]:offsets[i + 1]]; a missing blob is an empty span. The
    list is emptied while it is copied, so the blobs and the buffer are not both fully
    alive at the peak.
    """
    lengths = np.fromiter((0 if blob is None else len(blob) for blob in blobs),
                          dtype=np.int64, count=len(blobs))
    offsets = np.zeros(len(blobs) + 1, dtype=np.int64)
    np.cumsum(lengths, out=offsets[1:])

    buffer = np.empty(int(offsets[-1]), dtype=np.uint8)
    for index in range(len(blobs)):
        blob = blobs[index]
        blobs[index] = None
        if blob:
            buffer[offsets[index]:offsets[index + 1]] = np.frombuffer(blob, dtype=np.uint8)

    return buffer, offsets


def process_results_chunk(
        config: Config,
        results: list[list[int, bytes] | list[int, bytes, int, float, bytes]],
        stored_ids,
        db_nr: int,
        proteins: None | Proteins,
        skip_interfaces: bool
    ) -> (
        list | 
        tuple[
            list,
            list[tuple[int, np.ndarray]],
            list[tuple[int, list[tuple[str, np.ndarray]], int, float]]
            ]
        ):
    times = []
    t_init = time.time()

    if config.verbosity >= 4:
        config.logger.info(f'process resutls chunk {len(results)=} {db_nr=} {(proteins is None)=} {skip_interfaces=}')

    if proteins is None:
        structures: list[tuple[int, np.ndarray]] = []
        complexes: list[tuple[int, list[tuple[str, np.ndarray]], int, float]] = []

    # These matrices dominate the ray object store during aggregation and hold
    # nothing that needs 64 bit floats. They are written as float32 now, so this
    # only fires on rows stored before that change; the check below makes it a
    # no-op for everything else.
    narrow_residue_data = getattr(config, 'residue_data_float32', True)

    total_subtimes = []

    # The rows are dropped from the list while they are consumed. Each one holds
    # a compressed residue blob and a compressed interface blob, so keeping the
    # full result set alive next to the unpacked arrays roughly doubles the peak
    # of this step. The caller must treat 'results' as consumed afterwards.
    for row_index in range(len(results)):
        row = results[row_index]
        results[row_index] = None

        subtimes = []
        ti = time.time()

        s_id = row[0]
        internal_structure_db_id = 100*s_id + db_nr
        
        residue_data: np.ndarray[float] = unpack(row[1])

        ti = add_to_times(subtimes, ti, label = 'unpack res data')

        if residue_data is None:
            total_subtimes = aggregate_times(total_subtimes, subtimes)
            continue

        if narrow_residue_data and residue_data.dtype == np.float64:
            residue_data = residue_data.astype(np.float32)

        if proteins is None:
            """
            if pdb_id not in structures:
                structures[pdb_id] = {}
            if chain not in structures[pdb_id]:
                structures[pdb_id][chain] = residue_data
            """
            structures.append((internal_structure_db_id, residue_data))
        else:
            pdb_id, chain = stored_ids[internal_structure_db_id]
            proteins.structures[pdb_id][chain] = residue_data
        
        ti = add_to_times(subtimes, ti, label = 'build dict')

        if not skip_interfaces:
            compl_db_id = row[2]
            resolution = row[3]
            if resolution is None or np.isnan(resolution):
                resolution = 2.0
            interfaces: list[tuple[str, np.ndarray, np.ndarray]] = unpack(row[4])
            ti = add_to_times(subtimes, ti, label = 'unpack interface data')

            if proteins is None:
                """
                if pdb_id not in complexes:
                    complexes[pdb_id] = {}
                if chain not in complexes[pdb_id]:
                    complexes[pdb_id][chain] = interfaces
                """
                #if not unpacked:
                #    interfaces = ray.put(interfaces)
                complexes.append((internal_structure_db_id, interfaces, compl_db_id, resolution))
            else:
                if chain not in proteins.complexes[pdb_id].interfaces:
                    proteins.complexes[pdb_id].interfaces[chain] = interfaces
                

            ti = add_to_times(subtimes, ti, label = 'build interface dicts')
        else:
            ti = add_to_times(subtimes, ti, label = 'unpack interface data')
            ti = add_to_times(subtimes, ti, label = 'build interface dicts')

        total_subtimes = aggregate_times(total_subtimes, subtimes)

    times.append(total_subtimes)

    t_init = add_to_times(times, t_init, 'process residue table data core loop')
        
    if proteins is None:
        return times, structures, complexes
    return times

def split_bins_into_jobs(struct_db_bins: list[list[int]], max_structures_per_package: int) -> list[tuple[int, list[int]]]:
    """Turn the 20 Structure_{db_nr} bins into select packages of bounded size.

    Empty bins produce no job at all: nothing references a store that holds no
    structure, so spending a task and a set of empty object store entries on
    them is pointless.
    """
    jobs: list[tuple[int, list[int]]] = []

    for db_nr, bin_ids in enumerate(struct_db_bins):
        if len(bin_ids) == 0:
            continue

        if len(bin_ids) <= max_structures_per_package:
            jobs.append((db_nr, bin_ids))
            continue

        number_of_packages = -(-len(bin_ids) // max_structures_per_package)
        package_size = -(-len(bin_ids) // number_of_packages)

        for left in range(0, len(bin_ids), package_size):
            jobs.append((db_nr, bin_ids[left:left + package_size]))

    return jobs


def shared_memory_used_gib() -> float | None:
    """How much of /dev/shm is in use, which on Linux is where the object store lives.

    A plain statfs, so it is cheap enough to call per package. None where there is no
    /dev/shm.
    """
    try:
        return shutil.disk_usage('/dev/shm').used / 1024**3
    except OSError:
        return None


def log_package_progress(
        config: Config,
        done: int,
        total: int,
        package_stats: dict,
        package_totals: dict,
        t_dispatch: float,
        store_budget_state: dict | None = None) -> None:
    """One line per finished select package.

    A whole-PDB run used to go silent here for its entire duration: the per-task
    timings only existed at verbosity 4, so a run that took 19 hours left nothing to
    tell the database, the packing and a filling object store apart.

    The store capacity is part of the line because /dev/shm usage on its own never said
    how close the run was to the ceiling - the number that decides whether the stages
    after this one still have room.
    """
    shm = shared_memory_used_gib()
    shm_str = 'n/a' if shm is None else f'{shm:.1f} GiB'

    if store_budget_state is not None:
        if 'capacity' not in store_budget_state:
            store_budget_state['capacity'] = object_store_capacity_bytes()
        capacity = store_budget_state['capacity']
        if capacity:
            blob_bytes = package_totals['residue_bytes'] + package_totals['interface_bytes']
            shm_str = (f'{shm_str}, store {blob_bytes / 1024**3:.1f}/'
                       f'{capacity / 1024**3:.1f} GiB ({blob_bytes / capacity:.0%})')
    config.logger.info(
        f'getStoredResidues package {done}/{total}: '
        f'{package_stats["n_structures"]} structures, '
        f'select {package_stats["select_seconds"]:.1f}s, '
        f'pack {package_stats["pack_seconds"]:.1f}s, '
        f'put {package_stats["put_seconds"]:.1f}s, '
        f'blobs {(package_stats["residue_bytes"] + package_stats["interface_bytes"]) / 1024**2:.0f} MiB | '
        f'total {package_totals["n_structures"]} structures, '
        f'{(package_totals["residue_bytes"] + package_totals["interface_bytes"]) / 1024**3:.2f} GiB, '
        f'elapsed {time.time() - t_dispatch:.0f}s, /dev/shm used {shm_str}'
    )


def dispatch_selects(
        pending_jobs: list[tuple[int, int, list[int]]],
        in_flight: list[ray.ObjectRef],
        max_parallel: int,
        conf_dump: ray.ObjectRef,
        columns: list[str],
        skip_interfaces: bool,
        pr_is_none: bool) -> None:
    """Start select tasks until max_parallel of them are running. Mutates both lists."""
    while len(pending_jobs) > 0 and len(in_flight) < max_parallel:
        store_index, db_nr, id_package = pending_jobs.pop()
        in_flight.append(
            remote_select_package_from_db.remote(
                id_package, conf_dump, columns, db_nr, store_index, skip_interfaces, pr_is_none = pr_is_none
            )
        )


def getStoredResidues(
        proteins: Proteins | None,
        config: Config,
        custom_ids: None | dict[int, tuple[str, str]] = None,
        exclude_interacting_chains: bool = False,
        skip_interfaces: bool = False
        ) -> (
            None | 
            tuple[
                int,
                ray.ObjectRef[np.ndarray[int]],
                list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]],
                list
                ]
        ):
    t_total = time.time()
    total_times = []

    if custom_ids is None:
        stored_ids: dict[int, tuple[str, str]] = proteins.getStoredStructureIds(exclude_interacting_chains = exclude_interacting_chains)
        if config.verbosity >= 3:
            config.logger.info(f'Call of getStoredResidues, exclude_interacting_chains {exclude_interacting_chains}: {len(stored_ids)}')
    else:
        stored_ids = custom_ids
        if config.verbosity >= 3:
            config.logger.info(f'Call of getStoredResidues, exclude_interacting_chains {exclude_interacting_chains} with custom_ids: {len(custom_ids)}')
    
    filtered_stored_keys = [x for x in stored_ids.keys() if x is not None]
    if len(filtered_stored_keys) > 0:
        max_db_id = max(filtered_stored_keys)
    else:
        max_db_id = 0

    t_total = add_to_times(total_times, t_total, label = f'getStoredStructureIds {max_db_id}')
    
    # Sized by the largest structure db id rather than by the number of stored
    # structures, so it is mostly -1 padding. Columns: store_index, blob_index, three
    # unused slots, compl_db_id, resolution_index - all far inside int32. The unused
    # slots once held decoded row offsets; aggregation.build_local_stores now fills
    # them per protein, in a table of the same shape.
    total_structures: np.ndarray[np.int32] = np.full((max_db_id+1, 7), -1, dtype=np.int32)
    
    matrix_stores: list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]] = []

    if len(stored_ids) > 0:

        # ray may well be down here: the pipeline shuts it down when it finishes and is
        # skipped entirely when the session already exists, so the output stage arrives
        # with no instance. Without this, ray.put falls into ray's auto_init_hook and a
        # bare ray.init() starts an unconfigured cluster that times out under load.
        ray_init(config, stop_existing = False)

        conf_dump = ray.put(config)

        if skip_interfaces:
            columns = ['Structure_Id', 'Residue_Data']
        else:
            columns = ['Structure_Id', 'Residue_Data', 'Complex', 'Resolution', 'Interacting_Residue_Data']

        total_loop_times = []

        struct_db_bins = [[] for _ in range(20)]

        for internal_struct_db_id in stored_ids:
            if internal_struct_db_id is None:
                continue
            db_nr = internal_struct_db_id % 100
            struct_db_id = internal_struct_db_id//100

            struct_db_bins[db_nr].append(struct_db_id)

        # Sorted before packaging. The ids arrive in protein order, which is random
        # with respect to Structure_Id, so every package used to span its whole shard:
        # picky_select's density test failed and the package became one IN() over
        # 20000 scattered primary keys, each row's longblob stored off-page. Sorted,
        # a package is a contiguous id range and goes out as a single BETWEEN - on a
        # database filled by the same run, almost exactly the rows it asks for.
        # Nothing depends on the order: total_structures maps every id to its store.
        for bin_ids in struct_db_bins:
            bin_ids.sort()

        jobs: list[tuple[int, list[int]]] = split_bins_into_jobs(struct_db_bins, config.max_structures_per_select_package)

        # One matrix store per job. The store index is what gets written into
        # total_structures, so a bin split over several jobs simply produces
        # several stores; every offset in a packed structure is relative to the
        # store of its own job.
        matrix_stores = [None for _ in range(len(jobs))]

        max_parallel_selects = min(config.max_parallel_structure_selects, config.proc_n, len(jobs))
        if max_parallel_selects < 1:
            max_parallel_selects = 1

        if config.verbosity >= 3:
            config.logger.info(
                f'getStoredResidues dispatching {len(jobs)} select packages '
                f'(max {config.max_structures_per_select_package} structures each), '
                f'at most {max_parallel_selects} in parallel'
            )

        # Running totals for the per-package progress lines.
        package_totals = {'residue_bytes': 0, 'interface_bytes': 0, 'n_structures': 0}
        # Carries the object store capacity between the per-package budget checks.
        store_budget_state: dict = {}
        t_dispatch = time.time()

        # Dispatched back to front so that pop() takes them in order.
        pending_jobs: list[tuple[int, int, list[int]]] = [
            (store_index, db_nr, id_package) for store_index, (db_nr, id_package) in reversed(list(enumerate(jobs)))
        ]
        remote_selects_process_ids = []

        tloop = time.time()
        times = []
        dispatch_selects(pending_jobs, remote_selects_process_ids, max_parallel_selects, conf_dump, columns, skip_interfaces, (proteins is None))
        tloop = add_to_times(times, tloop, label = "remote call")
        total_loop_times = aggregate_times(total_loop_times, times)

        total_times.append(total_loop_times)
        t_total = add_to_times(total_times, t_total, label='call loop')

        if len(remote_selects_process_ids) > 0:

            returned_pkgs = set()
            total_loop_times = []
            while True:
                tloop = time.time()
                loop_times = []
                ready, not_ready = ray.wait(remote_selects_process_ids, timeout = 0.1)

                tloop = add_to_times(loop_times, tloop, label = "wait")

                if len(ready) > 0:
                    result_chunks = ray.get(ready)
                    total_inner_loop_times = []
                    for results_chunk in result_chunks:
                        inner_loop_times = []

                        if proteins is not None:

                            results, db_nr, pkg_number = results_chunk
                            proc_times = process_results_chunk(config, results, stored_ids, db_nr, proteins, skip_interfaces)
                            del results
                        else:

                            package_stats, packed_structures, blob_store, pkg_number = results_chunk
                            matrix_stores[pkg_number] = blob_store

                            if len(packed_structures) > 0:
                                total_structures[packed_structures[:, 0]] = packed_structures[:, 1:]

                            proc_times = []
                            package_totals['residue_bytes'] += package_stats['residue_bytes']
                            package_totals['interface_bytes'] += package_stats['interface_bytes']
                            package_totals['n_structures'] += package_stats['n_structures']

                            if config.verbosity >= 2:
                                log_package_progress(config, len(returned_pkgs) + 1, len(jobs),
                                                     package_stats, package_totals, t_dispatch,
                                                     store_budget_state)

                            # The blobs stay in the store for the whole of aggregation,
                            # so this is the running claim of this stage, not a peak.
                            check_store_budget(
                                config,
                                package_totals['residue_bytes'] + package_totals['interface_bytes'],
                                store_budget_state)

                        returned_pkgs.add(pkg_number)
                        inner_loop_times.append(proc_times)
                        total_inner_loop_times = aggregate_times(total_inner_loop_times, inner_loop_times)
                    del result_chunks
                    loop_times.append(inner_loop_times)
                else:
                    loop_times.append([])

                tloop = add_to_times(loop_times, tloop, label = "results procession")
                total_loop_times = aggregate_times(total_loop_times, loop_times)

                remote_selects_process_ids = list(not_ready)
                # Only top up once results have been taken off the previous batch,
                # so the number of tasks holding retrieved residue data at the same
                # time never exceeds max_parallel_selects.
                dispatch_selects(pending_jobs, remote_selects_process_ids, max_parallel_selects, conf_dump, columns, skip_interfaces, (proteins is None))

                if len(remote_selects_process_ids) == 0 and len(returned_pkgs) == len(jobs):
                    break

            total_times.append(total_loop_times)
            t_total = add_to_times(total_times, t_total, label='wait loop')

        else:
            # Nothing to select (every stored id was None). Keep the shape of
            # total_times identical to the dispatching case.
            total_times.append([])
            t_total = add_to_times(total_times, t_total, label='wait loop')

    else:
        t_total = add_to_times(total_times, t_total, label=f'loop prep {len(stored_ids)=}')
        total_times.append([])
        t_total = add_to_times(total_times, t_total, label='call loop')
        total_times.append([])
        t_total = add_to_times(total_times, t_total, label='wait loop')

                    
    t_total = add_to_times(total_times, t_total, label = 'main proc')

    if proteins is None:
        total_structures = ray.put(total_structures)

        return max_db_id, total_structures, matrix_stores, total_times
       


def retrieve_structure_info_from_struct_ids(config: Config, struct_ids, db_nr: int, cols = None) -> list[tuple[str, int, str]]:
    if cols is None:
        columns = ['PDB', 'Structure_Id', 'Chain']
    else:
        columns = cols
    table = f'Structure_{db_nr}'
    if len(struct_ids) < 5000:
        results = select(config, columns, table, in_rows={'PDB':struct_ids}, locked=True)
    else:
        results = select(config, columns, table, locked=True)
    return results

def get_alignment_data(config: Config, prot_db_ids: list[int]):
    columns = ['Protein', 'Structure', 'Sequence_Identity', 'Coverage', 'Alignment']

    binned_prot_db_ids = [[] for _ in range(20)]

    for prot_db_id in prot_db_ids:
        binned_prot_db_ids[prot_db_id%20].append(prot_db_id)

    prot_structure_alignment_map = {}

    for db_nr, prot_db_ids in enumerate(binned_prot_db_ids):
        table = f'Alignment_{db_nr}'
        results = binningSelect(prot_db_ids, columns, table, config)

        for row in results:

            prot_db_id = row[0]

            structure_id = row[1]
            seq_id = row[2]
            coverage = row[3]
            packed_alignment = row[4]

            if prot_db_id not in prot_structure_alignment_map:
                prot_structure_alignment_map[prot_db_id] = {}

            prot_structure_alignment_map[prot_db_id][structure_id] = (unpack(packed_alignment), coverage, seq_id)
    return prot_structure_alignment_map

def get_structure_alignment_context(
        config: Config,
        internal_structure_db_ids
        ) -> dict[int, tuple[str, str, list, str, np.ndarray | None]]:
    """Everything an alignment matrix cell needs, without touching Residue_Data.

    Returns {fused structure db id: (pdb_id, chain, seq_res_map, template_seq,
    simple_class_codes)}. seq_res_map and template_seq come out of the single
    Parse_Data blob and are both indexed by residue contig index - 1, the same
    space Alignment_N.Alignment points into. simple_class_codes is indexed the
    same way, and is None for a structure annotated before Simple_Class_Data
    existed.

    Residue_Data is deliberately not selected: it is ~780 bytes per residue and
    the classification is the only thing in it this output needs.
    """
    # Structure ids arrive fused as Structure_Id * 100 + db_nr, so the shard is
    # already known and no hashing of the pdb id is needed.
    binned_struct_db_ids: list[list[int]] = [[] for _ in range(20)]
    for internal_structure_db_id in internal_structure_db_ids:
        if internal_structure_db_id is None:
            continue
        binned_struct_db_ids[internal_structure_db_id % 100].append(internal_structure_db_id // 100)

    columns = ['Structure_Id', 'PDB', 'Chain', 'Parse_Data', 'Simple_Class_Data']

    context: dict[int, tuple[str, str, list, str, np.ndarray | None]] = {}

    for db_nr, struct_db_ids in enumerate(binned_struct_db_ids):
        if len(struct_db_ids) == 0:
            continue

        results = binningSelect(struct_db_ids, columns, f'Structure_{db_nr}', config)

        for row in results:
            internal_structure_db_id = 100 * row[0] + db_nr

            packed_parse_data = row[3]
            if packed_parse_data is None:
                continue
            (seq_res_map, template_seq) = unpack(packed_parse_data)

            packed_simple_class_data = row[4]
            if packed_simple_class_data is None:
                simple_class_codes = None
            else:
                simple_class_codes = unpack(packed_simple_class_data)

            context[internal_structure_db_id] = (
                row[1], row[2], seq_res_map, template_seq, simple_class_codes
            )

    if config.verbosity >= 3:
        config.logger.info(
            f'get_structure_alignment_context: {len(context)} of '
            f'{sum(len(b) for b in binned_struct_db_ids)} structures'
        )

    return context
