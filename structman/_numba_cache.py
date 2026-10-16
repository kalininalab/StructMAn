"""Give every numba/llvmlite pair its own on-disk kernel cache.

numba keys a cached kernel on the numba version, the source file's mtime and size, the bytecode
and the CPU - but not on the llvmlite / LLVM version that compiled it. The `alex` branch requires
numba>=0.67 (llvmlite 0.49, LLVM 22) and `fetch_alignments` numba<0.66 (LLVM 20), so a checkout
plus reinstall swaps the whole stack while the kernels stay in __pycache__ next to the installed
modules - pip only removes files it installed itself. The first run after switching back to
`fetch_alignments` died in every annotate worker with

    RuntimeError: LLVM bitcode parsing error
    Unknown attribute kind (105) (Producer: 'LLVM22.1.0' Reader: 'LLVM 20.1.8')

Pointing NUMBA_CACHE_DIR at a directory named after the pair makes such a load impossible: a
different stack finds an empty directory and compiles, once. numba reads NUMBA_CACHE_DIR when
numba.config is imported, and ray workers inherit the driver's environment, so structman sets it
at package import, before anything imports numba.

Pure standard library, and it never raises: a package import must not fail over a cache location.
"""
import os
import sys
import tempfile

ENV = 'NUMBA_CACHE_DIR'
DEFAULT_BASE = os.path.join('~', '.cache', 'structman', 'numba')


def stack_tag(version_of=None) -> str | None:
    """'numba-<version>_llvmlite-<version>', or None when either version cannot be read.

    Read from the installed distributions' metadata, so numba itself is not imported - that would
    read NUMBA_CACHE_DIR before it is set.
    """
    if version_of is None:
        from importlib.metadata import version as version_of
    try:
        return f'numba-{version_of("numba")}_llvmlite-{version_of("llvmlite")}'
    except Exception:
        return None


def cache_dir_for(tag: str, current: str | None) -> str:
    """The tagged cache directory under 'current' (a user's NUMBA_CACHE_DIR) or the default.

    Idempotent: a ray worker inherits the driver's already-tagged value, and appending the tag
    again would nest it and give every worker generation a cache of its own.
    """
    base = current or os.path.expanduser(DEFAULT_BASE)
    if os.path.basename(os.path.normpath(base)) == tag:
        return base
    return os.path.join(base, tag)


def writable_dir(path: str) -> bool:
    try:
        os.makedirs(path, exist_ok=True)
        with tempfile.TemporaryFile(dir=path):
            pass
        return True
    except OSError:
        return False


def apply(environ=None, modules=None, version_of=None, writable=writable_dir) -> str | None:
    """Point numba's cache at the directory of the installed numba/llvmlite pair.

    Returns that directory, or None when it was left alone: versions unreadable, or the directory
    cannot be written (a container may have no writable home). numba then keeps its default,
    next to the sources. If numba was already imported, its config is updated as well, since the
    environment variable alone would come too late for it.
    """
    try:
        environ = os.environ if environ is None else environ
        modules = sys.modules if modules is None else modules
        tag = stack_tag(version_of)
        if tag is None:
            return None
        path = cache_dir_for(tag, environ.get(ENV))
        if not writable(path):
            return None
        environ[ENV] = path
        numba = modules.get('numba')
        if numba is not None:
            try:
                numba.config.CACHE_DIR = path
            except Exception:
                pass
        return path
    except Exception:
        return None
