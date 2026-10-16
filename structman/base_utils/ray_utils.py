import ray
from ray._private import ray_constants
import json
import os
import shutil
import threading
from structman import settings
from structman.base_utils.config_class import Config
import subprocess
import time
import sys
import traceback
import psutil


def _read_cgroup_limit_file(path):
    """Read one cgroup limit file, returning bytes or None when it is unlimited."""
    try:
        with open(path) as f:
            val = f.read().strip()
    except (FileNotFoundError, PermissionError, NotADirectoryError, OSError):
        return None

    if val == "max":
        return None
    try:
        limit = int(val)
    except ValueError:
        return None

    # Values near INT64_MAX mean unlimited in cgroup v1
    if limit >= (2**63 - 4096):
        return None
    return limit


def _get_cgroup_memory_limit_gigs():
    """Return the effective cgroup memory limit in GiB for this process, or None if unlimited.

    The limit that actually kills us is rarely on the leaf cgroup. Under systemd
    a login session sits in .../user.slice/user-<uid>.slice/session-<n>.scope and
    the memory limit is set on the user slice, while the scope itself reads
    'max'. Only inspecting the leaf therefore reported "no limit" and let Ray
    size its object store against total system memory, well above what the
    cgroup would actually allow. Every ancestor constrains us, so take the
    smallest limit found along the whole path.

    memory.high is included: it is a soft limit, but exceeding it puts the
    cgroup into aggressive reclaim, which stalls processes badly enough for the
    raylet to miss its heartbeats.
    """
    try:
        cgroup_path = None
        v1_memory_path = None
        with open("/proc/self/cgroup") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(":", 2)
                if len(parts) < 3:
                    continue
                # cgroup v2 unified hierarchy: id=0, subsystems=""
                if parts[0] == "0" and parts[1] == "":
                    cgroup_path = parts[2]
                # cgroup v1: look for the memory subsystem
                elif "memory" in parts[1].split(",") and v1_memory_path is None:
                    v1_memory_path = parts[2]

        limits = []

        for base, rel_path, file_names in (
            ("/sys/fs/cgroup", cgroup_path, ("memory.max", "memory.high")),
            ("/sys/fs/cgroup/memory", v1_memory_path, ("memory.limit_in_bytes",)),
        ):
            if rel_path is None:
                continue

            # Walk from the leaf up to the cgroup root.
            parts = [p for p in rel_path.split("/") if p]
            while True:
                directory = "/".join([base] + parts)
                for file_name in file_names:
                    limit = _read_cgroup_limit_file(f"{directory}/{file_name}")
                    if limit is not None:
                        limits.append(limit)
                if not parts:
                    break
                parts.pop()

        if limits:
            return min(limits) / 1024**3

    except Exception:
        pass

    return None

def _cgroup_v2_path() -> str | None:
    """This process's cgroup v2 path, or None on cgroup v1 or without /proc."""
    try:
        with open("/proc/self/cgroup") as f:
            for line in f:
                parts = line.strip().split(":", 2)
                if len(parts) == 3 and parts[0] == "0" and parts[1] == "":
                    return parts[2]
    except OSError:
        pass
    return None


def memory_load_percent() -> float:
    """Memory in use as a percentage of the limit that actually applies to this run.

    psutil.virtual_memory() measures the whole machine. On the compute server the
    binding limit is the user slice's cgroup, 1435.4 of 1510.9 GiB - 95.0% of the
    machine - so a "stop dispatching at 95%" check on psutil could only fire once the
    cgroup was already full. On a verbosity 4 whole-PDB run, para_classify drew
    available memory from 827 GiB to 102 GiB in one minute without that check ever
    triggering.

    Walks from the leaf cgroup to the root and reports the fullest level, since any
    level with a limit can be the one that throttles. Falls back to psutil where there
    is no cgroup v2 hierarchy to read.
    """
    rel_path = _cgroup_v2_path()
    worst = None
    if rel_path is not None:
        parts = [p for p in rel_path.split("/") if p]
        while True:
            directory = "/".join(["/sys/fs/cgroup"] + parts)
            limits = [limit for limit in (
                _read_cgroup_limit_file(f"{directory}/memory.max"),
                _read_cgroup_limit_file(f"{directory}/memory.high"),
            ) if limit is not None]
            if limits:
                try:
                    with open(f"{directory}/memory.current") as f:
                        current = int(f.read().strip())
                    share = 100.0 * current / min(limits)
                    worst = share if worst is None else max(worst, share)
                except (OSError, ValueError):
                    pass
            if not parts:
                break
            parts.pop()
    if worst is not None:
        return worst
    return psutil.virtual_memory()[2]


def process_memory_breakdown(pid: str | int = 'self') -> dict[str, int] | None:
    """This process's RSS split into private and shared, from /proc/<pid>/smaps_rollup.

    RSS on its own was what misled the classification memory reservation. A
    para_classify worker maps the whole backmap and residue stores out of the object
    store, and its RSS grows as it touches them - 50 GB of growth that correlated with
    elapsed time (+0.78) and not with the package's cost (-0.16). Those pages are
    shared: 255 workers touching them cost the node one copy, not 255.

    Three numbers, because they answer three different questions:

    * uss (Private_Clean + Private_Dirty) - what one more worker costs, so what a
      reservation has to cover;
    * pss - the proportional share, so summing it over the workers gives what the node
      actually pays;
    * rss - what it used to report, kept so old and new lines can be compared.

    One read of one file, so it is cheap enough to call per task. None where the kernel
    does not provide smaps_rollup (before 4.14) or the read fails, and callers then just
    do not report the split.
    """
    fields = {
        'Rss': 'rss', 'Pss': 'pss',
        'Private_Clean': 'private_clean', 'Private_Dirty': 'private_dirty',
        'Shared_Clean': 'shared_clean', 'Shared_Dirty': 'shared_dirty',
    }
    out: dict[str, int] = {}
    try:
        with open(f'/proc/{pid}/smaps_rollup') as rollup:
            for line in rollup:
                name, _, rest = line.partition(':')
                key = fields.get(name)
                if key is None:
                    continue
                # 'Rss:  2268 kB'
                out[key] = int(rest.split()[0]) * 1024
    except (OSError, ValueError, IndexError):
        return None

    if 'private_clean' not in out or 'private_dirty' not in out:
        return None

    out['uss'] = out['private_clean'] + out['private_dirty']
    out['shared'] = out.get('shared_clean', 0) + out.get('shared_dirty', 0)
    return out


def parallel_db_access(config: Config) -> bool:
    """Whether the 20-way sharded database helpers should fan out over ray.

    False means the caller runs them in a plain loop and ray is never started.

    read_only runs can never reach a parallel stage: protCheck deletes every protein
    that is not already stored, and size_sorted skips stored proteins, so
    para_classify is never handed a package; paraAlignment skips completely stored
    proteins and paraAnnotate finds no unstored complex. All that is left is the
    sharded selects and the output stage, and ray.init cost ~33s of a ~60s read_only
    run - it forks one python worker per cpu (services.py passes num_cpus straight to
    --num_prestart_python_workers), so a big machine pays the most.

    proc_n == 1 takes the serial path too. Those helpers used to require ray even
    then, which is why a single-core run still had to start a cluster.

    defer_ray_init is the small-input case (see serializedPipeline.main): the whole
    run up to and including protCheck stays serial, and only if protCheck finds a
    protein that is not in the database yet does core() clear the flag and start ray.
    An input of a handful of already stored proteins therefore never starts a cluster.
    """
    if getattr(config, 'read_only_mode', False):
        return False
    if getattr(config, 'defer_ray_init', False):
        return False
    return config.proc_n > 1


def deref(value):
    """ray.get an ObjectRef, pass anything else straight through.

    The sharded database helpers hand back ObjectRefs when they fan out over ray and
    the plain arrays when they run serially (see parallel_db_access), so everything
    downstream has to accept either.
    """
    if isinstance(value, ray.ObjectRef):
        return ray.get(value)
    return value


# Share of the free /dev/shm the object store may claim. Ray refuses to start when the
# requested store exceeds it, and the rest has to stay free for other shm users.
SHM_STORE_SHARE = 0.95


def shared_memory_available_bytes() -> int | None:
    """Free bytes in /dev/shm, computed the way ray itself does before it starts the store.

    None where there is no /dev/shm, in which case ray picks another location itself.
    """
    try:
        stats = os.statvfs('/dev/shm')
    except OSError:
        return None
    return stats.f_bsize * stats.f_bavail


def clamp_object_store_to_shm(requested_bytes: int, shm_available_bytes: int | None) -> int:
    """The object store size ray will actually accept.

    Ray places the store in /dev/shm and raises a ValueError at startup when the request
    is larger than what /dev/shm has free (unless RAY_OBJECT_STORE_ALLOW_SLOW_STORAGE
    puts it on disk). /dev/shm defaults to half of RAM, so raising
    ray_object_store_fraction past that used to turn into a failed ray.init instead of
    a bigger store. Clamped here, the fraction can be raised for a big input and gets
    as much as the machine can give.
    """
    if shm_available_bytes is None:
        return requested_bytes
    return min(requested_bytes, int(SHM_STORE_SHARE * shm_available_bytes))


# How often the overflow watchdog asks the raylet for its store statistics. A cheap
# local gRPC call; at whole-PDB fill rates a few GiB can reach disk between two polls.
SPILL_CHECK_SECONDS = 10.0

# Exit status of a run stopped by the watchdog, distinct from the pipeline's own codes.
OBJECT_STORE_OVERFLOW_EXIT_CODE = 3

# Share of the store taken by primary copies - objects nothing can evict - at which the
# log gets a one-time warning that an overflow is coming.
STORE_PRESSURE_WARNING_SHARE = 0.9

_watchdog_lock = threading.Lock()
_watchdog_thread: threading.Thread | None = None

# Whether the cluster this process talks to was started by ray_init, i.e. carries the
# settings below. False while a cluster is up that ray started on its own through
# auto_init_hook - then none of _system_config applies and ray_init can only say so.
_cluster_configured = False
# Waited on instead of time.sleep: setting it stops the thread, and it is independent of
# the module's time binding, which tests replace.
_watchdog_stop = threading.Event()


def read_object_store_stats() -> dict | None:
    """Store statistics of this node's raylet, or None when they cannot be read.

    Uses ray's internal node_stats call, the one `ray memory --stats-only` is built on.
    It is not a stable API, so every failure - ray down, attribute renamed, RPC error -
    comes back as None and the watchdog just skips the poll rather than stopping a
    healthy run over a monitoring problem.
    """
    try:
        if not ray.is_initialized():
            return None
        from ray._private import internal_api
        node = ray._private.worker.global_worker.node
        reply = internal_api.node_stats(
            node.node_ip_address, node.node_manager_port, include_memory_info=False)
        stats = reply.store_stats
        return {
            'used': int(stats.object_store_bytes_used),
            'capacity': int(stats.object_store_bytes_avail),
            'primary': int(stats.object_store_bytes_primary_copy),
            'fallback': int(stats.object_store_bytes_fallback),
            'spilled': int(stats.spilled_bytes_total),
        }
    except Exception:
        return None


def abort_on_object_store_overflow(config: Config, stats: dict) -> None:
    """Stop the whole process: the object store overflowed and spilling is not allowed.

    os._exit rather than an exception. core() ends in a bare except that logs and moves
    on to the next chunk, so any exception - ray's own ObjectStoreFullError included -
    would be swallowed, and this runs on a watchdog thread while the main thread is
    usually blocked inside ray.get, where no exception could reach it anyway.
    """
    GiB = 1024**3
    message = (
        f'Object store overflow, stopping the run (ray_allow_spilling is False): '
        f'{stats["spilled"] / GiB:.2f} GiB spilled, {stats["fallback"] / GiB:.2f} GiB fallback-allocated '
        f'to disk, store {stats["used"] / GiB:.1f} of {stats["capacity"] / GiB:.1f} GiB in use, '
        f'{stats["primary"] / GiB:.1f} GiB of it primary copies that cannot be evicted. '
        f'Give the store more room (--ray_object_store_fraction, capped by /dev/shm) or '
        f'split the input; set ray_allow_spilling=True to let ray spill instead.'
    )
    try:
        config.errorlog.add_error(message)
    except Exception:
        pass
    try:
        config.logger.error(message)
        for handler in config.logger.handlers:
            handler.flush()
    except Exception:
        pass
    print(message, file=sys.stderr, flush=True)
    os._exit(OBJECT_STORE_OVERFLOW_EXIT_CODE)


def object_store_capacity_bytes(read_stats=read_object_store_stats) -> int | None:
    """Capacity of the running cluster's object store, or None when it cannot be read."""
    stats = read_stats()
    if stats is None or stats.get('capacity', 0) <= 0:
        return None
    return int(stats['capacity'])


def check_store_budget(config: Config, used_bytes: int, state: dict,
                       capacity_reader=object_store_capacity_bytes) -> None:
    """Stop the run when one stage has claimed more of the object store than its share.

    getStoredResidues filled the store to 99% on a whole-PDB input and left nothing for
    getBackmaps and the classification packages that come after it, so the run died
    later, in a place that said nothing about the cause. This fails in the stage that
    actually consumed the store, and says how much of what.

    'state' carries the capacity between calls, so the stats are read once per stage.
    """
    share = getattr(config, 'ray_residue_store_share', 0.7)
    if share <= 0:
        return

    if 'capacity' not in state:
        state['capacity'] = capacity_reader()
    capacity = state['capacity']
    if not capacity:
        return

    budget = share * capacity
    if used_bytes <= budget:
        return

    GiB = 1024**3
    message = (
        f'Object store budget exceeded: this stage holds {used_bytes / GiB:.1f} GiB of the '
        f'{capacity / GiB:.1f} GiB object store, more than the {share:.0%} it may take '
        f'(ray_residue_store_share). The stages after it - getBackmaps and the '
        f'classification packages - need the rest, so the run would overflow the store '
        f'there instead. Give the store more room (ray_object_store_fraction, capped by '
        f'the free /dev/shm), split the input into smaller chunks, or raise '
        f'ray_residue_store_share if you know the rest fits.'
    )
    try:
        config.errorlog.add_error(message)
    except Exception:
        pass
    try:
        config.logger.error(message)
        for handler in config.logger.handlers:
            handler.flush()
    except Exception:
        pass
    print(message, file=sys.stderr, flush=True)
    # Same reasoning as abort_on_object_store_overflow: core() ends in a bare except
    # that would swallow an exception raised here and carry on to the next chunk.
    os._exit(OBJECT_STORE_OVERFLOW_EXIT_CODE)


def spill_watchdog_tick(config: Config, state: dict, read_stats=read_object_store_stats,
                        abort=abort_on_object_store_overflow) -> None:
    """One poll: abort on any spilled or fallback-allocated byte, warn once under pressure.

    Split out from the thread so it can be exercised without a running cluster.
    """
    stats = read_stats()
    if stats is None:
        return

    if stats['spilled'] > 0 or stats['fallback'] > 0:
        abort(config, stats)
        return

    if (not state.get('warned')
            and stats['capacity'] > 0
            and stats['primary'] >= STORE_PRESSURE_WARNING_SHARE * stats['capacity']):
        state['warned'] = True
        try:
            config.logger.info(
                f'Object store under pressure: {stats["primary"] / 1024**3:.1f} of '
                f'{stats["capacity"] / 1024**3:.1f} GiB held by primary copies. The run stops '
                f'if anything has to spill.')
        except Exception:
            pass


def stop_spill_watchdog() -> None:
    """Stop the watchdog thread; mainly for tests and embedding callers."""
    _watchdog_stop.set()


def warn_about_foreign_cluster(config: Config) -> None:
    """Say that the cluster ray_init found is not one of ours, and guard it anyway.

    A cluster ray started through its auto_init_hook has none of ray_init's settings:
    spilling is on and aimed at /tmp, the object store is sized by ray's defaults
    rather than by ray_object_store_fraction and the free /dev/shm, and num_cpus is
    every core on the node. The pipeline cannot fix that without throwing away whatever
    is already in the store, so it reports it and starts the watchdog, which turns the
    eventual overflow into a clean stop instead of an ENOSPC in the raylet.
    """
    message = (
        'ray was already initialized when ray_init was called, by something other than '
        'ray_init - most likely a ray.put/ray.get/remote call that ran before it and '
        'fell into ray auto_init_hook. This run therefore uses ray default '
        'configuration: object spilling enabled into /tmp, default object store size, '
        'all cores. Find that call and give it a ray_init(config) first.'
    )
    try:
        config.errorlog.add_warning(message)
    except Exception:
        pass
    try:
        config.logger.warning(message)
    except Exception:
        pass

    if not getattr(config, 'ray_allow_spilling', False):
        start_spill_watchdog(config)


def start_spill_watchdog(config: Config) -> None:
    """Start the overflow watchdog once per process; later calls are no-ops.

    It outlives ray.shutdown and re-inits - the pipeline starts the cluster lazily more
    than once - and simply skips polls while ray is down. Daemon, so it never keeps the
    process alive.
    """
    global _watchdog_thread
    with _watchdog_lock:
        if _watchdog_thread is not None and _watchdog_thread.is_alive():
            return

        state: dict = {}
        _watchdog_stop.clear()

        def _run():
            while not _watchdog_stop.wait(SPILL_CHECK_SECONDS):
                try:
                    spill_watchdog_tick(config, state)
                except Exception:
                    # Monitoring must never be the thing that breaks a run.
                    pass

        _watchdog_thread = threading.Thread(target=_run, name='structman-spill-watchdog', daemon=True)
        _watchdog_thread.start()

    if config.verbosity >= 2:
        config.logger.info(
            f'Object store spilling disabled; the run stops if the store overflows '
            f'(checked every {SPILL_CHECK_SECONDS:.0f}s)')


def worker_startup_system_config(config: Config) -> dict:
    """Ray settings for how its worker processes start, or {} to leave Ray's defaults.

    Most stages run their tasks with max_calls=1, so every task starts a fresh Python process,
    and Ray lets as many start at once as there are CPUs. On a whole-PDB chunk that meant
    ~15,000 starts in structural analysis alone, 255 at a time, each importing the environment
    from NFS - and a start that has not registered with the raylet within
    worker_register_timeout_seconds (Ray's default: 60) counts as failed. 512 of them failed
    on one chunk. The initial burst of para_parse_pdb was already throttled to proc_n // 4
    against the same thundering herd; this applies the throttle to every worker start, and
    gives a slow but healthy start under load time to register.

    Both keys exist in Ray's core config. If a Ray version rejects them, ray_init drops them
    for the rest of the run (ray_worker_startup_config_rejected) instead of failing.
    """
    if getattr(config, 'ray_worker_startup_config_rejected', False):
        return {}
    settings = {}
    # -1 (the default) means proc_n // 4, the throttle para_parse_pdb already used; 0 leaves Ray's.
    concurrency = getattr(config, 'ray_worker_startup_concurrency', -1)
    if concurrency is not None and concurrency < 0:
        concurrency = max(8, config.proc_n // 4)
    if concurrency:
        settings['worker_maximum_startup_concurrency'] = int(concurrency)
    timeout = getattr(config, 'ray_worker_register_timeout_seconds', 300)
    if timeout:
        settings['worker_register_timeout_seconds'] = int(timeout)
    return settings


RAY_MEMORY_THRESHOLD_ENV = 'RAY_memory_usage_threshold'
RAY_DEFAULT_MEMORY_THRESHOLD = 0.95


def ray_memory_usage_threshold(cgroup_limit_gigs: float | None, system_gigs: float | None) -> float | None:
    """Ray's OOM-kill threshold that makes it fire at 95% of the cgroup limit, or None to leave it.

    Ray's memory monitor kills a task - which Ray then retries - once memory use passes
    memory_usage_threshold (0.95). Measured against the whole machine, 0.95 of 1510.9 GiB is
    1435 GiB: exactly the user slice's cgroup limit (1435.4 GiB), where the kernel's OOM killer
    acts instead - and on 2026-09-24 it took the raylet and with it the whole run. Scaling the
    threshold by limit / machine puts Ray's recoverable kill in front of the kernel's. Should Ray
    measure the cgroup itself, it fires somewhat early, never late.
    """
    if not cgroup_limit_gigs or not system_gigs or cgroup_limit_gigs >= system_gigs:
        return None
    return round(RAY_DEFAULT_MEMORY_THRESHOLD * cgroup_limit_gigs / system_gigs, 3)


def _drop_rejected_startup_config(config: Config, startup_config: dict, err) -> bool:
    """After a failed ray.init that carried the worker startup settings: leave them out from now on.

    Returns whether they were dropped. A version of Ray that does not know one of the keys fails
    ray.init over it, and one unknown tuning key must not cost the run its cluster.
    """
    if not startup_config or getattr(config, 'ray_worker_startup_config_rejected', False):
        return False
    config.ray_worker_startup_config_rejected = True
    message = (f'ray.init failed with the worker startup settings {startup_config}; '
               f'retrying without them. {err}')
    try:
        config.errorlog.add_warning(message)
    except Exception:
        pass
    try:
        config.logger.warning(message)
    except Exception:
        pass
    return True


def ray_init(config: Config, n_try = 0, redis_mem_default = False, overwrite_logging_level = None, total_memory_quantile = 0.9, debug = False, num_gpus=0, stop_existing = True):
    global _cluster_configured

    if ray.is_initialized():
        # Normal case: an earlier ray_init in this process brought the cluster up and
        # this call has nothing to do. The other case is a cluster ray started itself,
        # which used to pass unnoticed and left the run with spilling enabled and no
        # watchdog - say so, and at least put the watchdog on it.
        if not _cluster_configured:
            warn_about_foreign_cluster(config)
        return

    if not config.python_env_expanded_for_ray:
        os.environ["PYTHONPATH"] = f'{settings.ROOT_DIR}:{os.environ.get("PYTHONPATH", "")}'
        os.environ["PYTHONPATH"] = f'{settings.LIB_DIR}:{os.environ.get("PYTHONPATH", "")}'
        os.environ["PYTHONPATH"] = f'{settings.RINERATOR_DIR}:{os.environ.get("PYTHONPATH", "")}'
        os.environ["PYTHONPATH"] = f'{settings.OUTPUT_DIR}:{os.environ.get("PYTHONPATH", "")}'
        
        config.python_env_expanded_for_ray = True

    logging_level = 20
    if config.verbosity <= 1:
        logging_level = 0

    if overwrite_logging_level is not None:
        logging_level = overwrite_logging_level

    if config.verbosity <= 8 and not debug:
        log_to_driver = False
    else:
        log_to_driver = True
        config.logger.info(f'Setting log_to_driver in ray.init to True, local mode: {config.ray_local_mode}')

    configured_limit_gigs = getattr(config, 'ray_memory_limit_gigs', 0.0)
    cgroup_limit_gigs = _get_cgroup_memory_limit_gigs()

    if configured_limit_gigs and configured_limit_gigs > 0:
        effective_gigs = configured_limit_gigs
        if config.verbosity >= 2:
            config.logger.info(
                f'Using configured ray_memory_limit_gigs: {effective_gigs:.1f} GiB '
                f'(system total: {config.gigs_of_ram:.1f} GiB, '
                f'detected cgroup limit: {cgroup_limit_gigs if cgroup_limit_gigs is None else round(cgroup_limit_gigs, 1)}).'
            )
    elif cgroup_limit_gigs is not None and cgroup_limit_gigs < config.gigs_of_ram:
        effective_gigs = cgroup_limit_gigs
        if config.verbosity >= 2:
            config.logger.info(
                f'Cgroup memory limit detected: {cgroup_limit_gigs:.1f} GiB (system total: {config.gigs_of_ram:.1f} GiB). '
                f'Using cgroup limit for Ray memory sizing.'
            )
    else:
        effective_gigs = config.gigs_of_ram
        if config.verbosity >= 2:
            config.logger.info(
                f'No cgroup memory limit below the system total found, sizing Ray against '
                f'{config.gigs_of_ram:.1f} GiB. If the process is nevertheless killed or stalled '
                f'by a cgroup limit, check memory.max/memory.high (cgroup v2) or '
                f'memory.limit_in_bytes (cgroup v1) on the parent slices, or set '
                f'ray_memory_limit_gigs explicitly.'
            )

    store_fraction = getattr(config, 'ray_object_store_fraction', 0.5)
    redis_mem = int(total_memory_quantile * store_fraction * effective_gigs * 1024 * 1024 * 1024)

    # read_only mode cannot reach anything that fills the object store, so sizing it
    # from total RAM reserves a large block of /dev/shm that the run will never touch.
    #   - protCheck deletes every protein that is not already stored, so every protein
    #     left is stored; size_sorted then skips stored proteins, so para_classify is
    #     never handed a package (a real read_only run logs total_cost=0).
    #   - paraAlignment skips completely stored proteins, which read_only marks them
    #     all as, and paraAnnotate finds no unstored complex - 0 structures analysed.
    #   - the output stage, which is where a read_only run does its work, calls
    #     getPositions and plain selects; it never calls getStoredResidues or
    #     getBackmaps, the two that build the residue matrices and the backmap
    #     concatenation.
    # What remains moves megabytes, so the store is capped rather than proportional.
    if getattr(config, 'read_only_mode', False):
        read_only_cap = int(getattr(config, 'ray_read_only_store_gigs', 32.0) * 1024**3)
        if redis_mem > read_only_cap:
            if config.verbosity >= 2:
                config.logger.info(
                    f'read_only mode: capping object store at '
                    f'{read_only_cap / 1024**3:.1f} GiB instead of {redis_mem / 1024**3:.1f} GiB'
                )
            redis_mem = read_only_cap

    requested_store = redis_mem
    redis_mem = clamp_object_store_to_shm(redis_mem, shared_memory_available_bytes())
    if redis_mem < requested_store and config.verbosity >= 1:
        config.logger.info(
            f'Object store clamped to {redis_mem / 1024**3:.1f} GiB, 95% of the free /dev/shm, '
            f'instead of the requested {requested_store / 1024**3:.1f} GiB '
            f'(ray_object_store_fraction={store_fraction}). A bigger store needs a bigger /dev/shm.'
        )

    if redis_mem > 200 * 1024**3:
        ray_constants.DEFAULT_OBJECT_STORE_MAX_MEMORY_BYTES = redis_mem

    # The pool ray admits memory= requests against, i.e. how much worker heap the run
    # believes it has. This used to be a hardcoded half of the node -
    # (total_memory_quantile*5/10) * effective_gigs - which was not derived from
    # store_fraction at all; the two agreeing at ~646 GiB on a 1435 GiB machine was a
    # coincidence of the constants, and together they claimed 1292 of 1435 GiB.
    #
    # The store's *reservation* costs nothing (the /dev/shm file is sparse, and a
    # whole-PDB chunk was measured using 12.2 of 645.9 GiB of it), but the pool being
    # half the machine caps concurrency for no reason. So the pool is what is actually
    # left: the effective limit, less the store's real appetite, less what the driver
    # holds - on a whole-PDB run its Proteins and annotation graph reach ~93 GiB, and
    # that is not available to workers whatever ray thinks.
    effective_bytes = int(total_memory_quantile * effective_gigs * 1024**3)
    # Both allowances are ceilings that also scale down with the machine. The absolute
    # figures come from a 1435 GiB node - the store holding 72.6 GiB of residue blobs
    # plus 60 GiB of backmaps, the driver reaching ~93 GiB - and taking them literally
    # on a 100 GiB machine subtracted more than the machine has, leaving only the floor.
    store_allowance = min(redis_mem, int(getattr(config, 'ray_store_allowance_gigs', 200.0) * 1024**3))
    driver_allowance = min(
        int(getattr(config, 'ray_driver_allowance_gigs', 120.0) * 1024**3),
        int(getattr(config, 'ray_driver_allowance_share', 0.15) * effective_bytes),
    )
    # Recorded so the classification dispatcher can size its window against the same
    # budget ray_init used, rather than re-deriving the cgroup limit itself.
    config.ray_effective_memory_bytes = effective_bytes
    mem = effective_bytes - store_allowance - driver_allowance
    # Never hand ray a non-positive pool, or nothing that reserves memory could ever be
    # scheduled.
    floor = int(getattr(config, 'ray_task_memory_floor_gigs', 8.0) * 1024**3)
    mem = max(mem, min(floor, effective_bytes))

    if config.verbosity >= 2:
        config.logger.info(
            f'Call of ray.init: num_cpus = {config.proc_n}, object_store_memory = {redis_mem:_} '
            f'({redis_mem / 1024**3:.1f} GiB, {store_fraction:.2f} of effective), '
            f'effective_gigs = {effective_gigs:.1f}'
        )
        config.logger.info(
            f'Ray memory split: effective {effective_bytes / 1024**3:.1f} GiB '
            f'- store allowance {store_allowance / 1024**3:.1f} GiB '
            f'- driver allowance {driver_allowance / 1024**3:.1f} GiB '
            f'= task pool {mem / 1024**3:.1f} GiB '
            f'({mem / 1024**3 / max(1, config.proc_n):.2f} GiB per core at proc_n={config.proc_n})'
        )

    # 'ray stop --force' walks every process it can see and SIGKILLs anything matching
    # ray's process names (ray/scripts/scripts.py). It is not scoped to this run, so on
    # a shared machine it takes down a concurrently starting run's GCS and raylet too,
    # leaving that run's ray.init waiting for a node that no longer exists. Only the
    # pipeline's own cold start has any business doing that; the lazy ray_init calls
    # that bring ray back up for the output stage pass stop_existing=False.
    if stop_existing:
        # The backoff used to be 10**loops seconds - 100, then 1000, then 10000 - so a
        # stderr message that never cleared parked the pipeline for over three hours,
        # indistinguishable from a freeze. communicate() had no timeout either, so a
        # wedged 'ray stop' blocked here forever. Both are bounded now: at most
        # STOP_ATTEMPTS tries, each capped at STOP_TIMEOUT, with a linear backoff.
        STOP_ATTEMPTS: int = 4
        STOP_TIMEOUT: int = 120

        for attempt in range(1, STOP_ATTEMPTS + 1):
            outs = b''
            errs = b''
            try:
                p = subprocess.Popen(['ray', 'stop', '--force'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                outs, errs = p.communicate(timeout=STOP_TIMEOUT)
            except subprocess.TimeoutExpired:
                p.kill()
                outs, errs = p.communicate()
                config.errorlog.add_warning(
                    f'ray stop did not finish within {STOP_TIMEOUT}s and was killed')
            except Exception:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                config.errorlog.add_warning(f'Ray stop failed:\n{e}\n{f}\n{g}')

            if config.verbosity >= 3:
                config.logger.info(f'Returns of ray stop:\n{outs}\n{errs}')

            # Compared as bytes: the old decode('ascii') would itself raise on any
            # non-ascii byte ray happened to print.
            if outs.count(b'Stopped') > 0:
                # Something really was killed; let its ports come free.
                time.sleep(5)

            if len(errs) == 0:
                break

            if attempt < STOP_ATTEMPTS:
                time.sleep(5 * attempt)

    if config.verbosity >= 8:
        include_dashboard = True
    else:
        include_dashboard = False

    _system_config={
        # Reduce the time the driver thread blocks when the store approaches capacity
        "object_store_full_delay_ms": 100,
        
        # Parallelize Ray's internal filesystem/spilling metadata operations
        "max_io_workers": 4,
        
        # Ensure Ray group evicts smaller or temporary objects first
        "min_spilling_size": 100 * 1024 * 1024
    }

    startup_config = worker_startup_system_config(config)
    _system_config.update(startup_config)
    if startup_config and config.verbosity >= 2:
        config.logger.info(f'Ray worker startup settings: {startup_config}')

    # Read by the raylet when it starts, so it has to be in the environment before ray.init.
    # A value the user set is left alone.
    if RAY_MEMORY_THRESHOLD_ENV not in os.environ:
        threshold = ray_memory_usage_threshold(cgroup_limit_gigs, getattr(config, 'gigs_of_ram', None))
        if threshold is not None:
            os.environ[RAY_MEMORY_THRESHOLD_ENV] = str(threshold)
            if config.verbosity >= 2:
                config.logger.info(
                    f'Ray OOM-kill threshold {threshold} ({RAY_MEMORY_THRESHOLD_ENV}): 95% of the '
                    f'{cgroup_limit_gigs:.1f} GiB cgroup limit on a {config.gigs_of_ram:.1f} GiB machine')

    allow_spilling = getattr(config, 'ray_allow_spilling', False)
    if not allow_spilling:
        # No spilling. A full store now blocks puts for ray's grace period and then
        # falls back to disk-backed allocations - which cannot be switched off - so the
        # watchdog started below is what actually stops the run.
        _system_config["automatic_object_spilling_enabled"] = False


    # Without this Ray spills objects and writes its logs into /tmp, which on most
    # compute nodes is far too small for the spill volume of a large run. Running
    # out of space there kills the IO workers with an ENOSPC that looks unrelated
    # to the actual problem.
    extra_init_kwargs = {}
    ray_temp_dir = getattr(config, 'ray_temp_dir', None)
    if ray_temp_dir:
        extra_init_kwargs['_temp_dir'] = ray_temp_dir
        if config.verbosity >= 2:
            config.logger.info(f'Setting ray temp/spill directory to {ray_temp_dir}')

    try:
        if not redis_mem_default:
            ray.init(
                num_cpus=config.proc_n, include_dashboard=include_dashboard,
                ignore_reinit_error=True, logging_level = logging_level,
                log_to_driver = log_to_driver,
                local_mode = config.ray_local_mode,
                object_store_memory = redis_mem, _memory = mem, num_gpus=num_gpus,
                _system_config=_system_config,
                **extra_init_kwargs
                )
            # The pool that tasks declaring memory= are admitted against. Recorded so
            # a request can be capped at it: a task asking for more than the whole
            # pool would never be scheduled and the stage would wait forever.
            config.ray_task_memory_bytes = mem
        else:
            # Ray sizes its memory resource itself here, and the size is not known
            # without a GCS round trip; callers then skip memory requests.
            config.ray_task_memory_bytes = None
            ray.init(
                num_cpus=config.proc_n, include_dashboard=include_dashboard,
                ignore_reinit_error=True, logging_level = logging_level,
                log_to_driver = log_to_driver,
                local_mode = config.ray_local_mode, num_gpus=num_gpus,
                _system_config=_system_config,
                **extra_init_kwargs
                )
    except ray.exceptions.RaySystemError as err:
        # The retry below then starts without the startup settings, in case they were the cause.
        _drop_rejected_startup_config(config, startup_config, err)
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        if n_try == 3:
            config.errorlog.add_error(f'Ray init failed: {err=}\n{e}\n{f}\n{g}')
            return

        config.errorlog.add_warning(f'Ray init failed, retry {n_try}... \n{e}\n{f}\n{g}')
        
        if stop_existing:
            outs = b''
            errs = b''
            try:
                p = subprocess.Popen(['ray', 'stop', '--force'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                outs, errs = p.communicate(timeout=120)
            except subprocess.TimeoutExpired:
                p.kill()
                outs, errs = p.communicate()
                config.errorlog.add_warning('ray stop did not finish within 120s and was killed')
            except Exception:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                config.errorlog.add_warning(f'Ray stop failed:\n{e}\n{f}\n{g}')

            if config.verbosity >= 3:
                config.logger.info(f'Returns of ray stop:\n{outs}\n{errs}')

        time.sleep(10*n_try)
        # stop_existing is carried through: a caller that asked not to kill other
        # runs' ray must not start doing so just because the first attempt failed.
        if n_try <= 4:
           ray_init(config, n_try = (n_try + 1), total_memory_quantile = total_memory_quantile - 0.05, num_gpus=num_gpus, stop_existing = stop_existing)
        else:
           ray_init(config, n_try = (n_try + 1), redis_mem_default = True, num_gpus=num_gpus, stop_existing = stop_existing)
    except Exception as err:
        # Anything else is retried only when the worker startup settings may be the reason - a
        # Ray version that does not know a _system_config key refuses to start over it.
        if not _drop_rejected_startup_config(config, startup_config, err):
            raise
        try:
            ray.shutdown()
        except Exception:
            pass
        return ray_init(config, n_try = n_try + 1, redis_mem_default = redis_mem_default,
                        overwrite_logging_level = overwrite_logging_level,
                        total_memory_quantile = total_memory_quantile, debug = debug,
                        num_gpus = num_gpus, stop_existing = stop_existing)

    # The cluster now up is one ray_init configured - either here or in the recursive
    # retry above. Read from ray rather than assumed, so a run of retries that ended
    # without a cluster does not leave the flag claiming one. It goes stale only if
    # something shuts ray down and lets ray auto-init it again, which cannot happen
    # while RAY_ENABLE_AUTO_CONNECT=0 (see structman/__init__.py).
    _cluster_configured = ray.is_initialized()

    if not allow_spilling:
        start_spill_watchdog(config)

    # Deliberately no ray.cluster_resources() here. It only reported the object store
    # size, which the 'Call of ray.init' line above already logs from redis_mem, and it
    # reaches GCS through GlobalStateAccessor.connect() - a second GCS connection,
    # independent of the core worker's, and a blocking call with no deadline. A run was
    # seen sitting in it for three hours with a healthy idle cluster and not a single
    # task submitted. Nothing here is worth an unbounded blocking call.

def ray_hack():
    # hack proposed by the devs of ray to prevent too many processes being spawned
    return
    #resources = ray.ray.get_resource_ids()
    #cpus = [v[0] for v in resources['CPU']]
    #psutil.Process().cpu_affinity(cpus)

def monitor_ray_store():
    p = subprocess.Popen(['ray','memory'])

    page, err = p.communicate()
    print(page)

    return
