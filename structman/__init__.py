import os

# Ray starts a cluster of its own the moment anything calls ray.put, ray.get or a
# remote function while no cluster is up: ray/_private/auto_init_hook.py wraps those
# APIs and calls a bare ray.init(). That cluster has none of the settings ray_utils.
# ray_init applies - object spilling stays on and writes into /tmp, the object store is
# sized by ray's defaults rather than by ray_object_store_fraction and the free
# /dev/shm, num_cpus is every core on the node, and the spill watchdog never starts.
# Worse, it is silent: every later ray_init returns at its 'if ray.is_initialized()'
# guard, so a single unguarded ray.put anywhere before the pipeline's own init decides
# the configuration of the whole run. A whole-PDB run was lost that way, filling /tmp
# with spilled objects until the raylet died of ENOSPC.
#
# With auto-connect off, such a call raises 'Ray has not been started yet' instead,
# which names the file and line that has to call ray_init first. Set here, at package
# import, because auto_init_hook reads the flag once when ray is imported. Workers
# inherit it and do not care: ray.is_initialized() is true inside a task.
os.environ.setdefault('RAY_ENABLE_AUTO_CONNECT', '0')

# numba's kernel cache is keyed on the numba version but not on the llvmlite/LLVM that compiled a
# kernel, so a reinstall that swaps llvmlite - every switch between the alex and fetch_alignments
# branches does - left kernels an older LLVM could not read, and the annotate workers died on
# 'LLVM bitcode parsing error'. Each numba/llvmlite pair gets a cache directory of its own; see
# _numba_cache.py. Set here for the same reason as above: numba reads NUMBA_CACHE_DIR once, when it
# is imported, and ray workers inherit the driver's environment.
from structman._numba_cache import apply as _apply_numba_cache
_apply_numba_cache()
