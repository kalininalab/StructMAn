"""Record which commit an install was built from, and say so when a run starts.

install.sh does a regular `pip install`, so a `git checkout` changes nothing on the server until
the next install - and _version.py is a static string that no log printed. On 2026-09-25 a
checkout of `alex` was installed two minutes after committing on `fetch_alignments`, and the next
three whole-PDB runs silently ran the older branch; their crashes were then read as failures of
code they never contained. Every install now carries its branch and commit, and every run and
errorlog names them.

Pure standard library: setup.py runs this before the package's dependencies are installed, the
same way it reads _version.py.
"""
import datetime
import os
import subprocess

UNKNOWN = 'unknown'
FIELDS = ('BRANCH', 'COMMIT', 'DIRTY', 'BUILT_AT')


def collect(repo_dir: str, run=subprocess.run) -> dict:
    """Branch, commit and dirty state of the checkout being installed.

    Every git failure degrades to UNKNOWN instead of raising: a public install from a tarball has
    no git, and stamping must never be the reason an install fails. 'run' is injectable so this
    can be tested without a repository.
    """
    def git(*args, timeout=30):
        try:
            completed = run(['git', *args], cwd=repo_dir, capture_output=True, text=True, timeout=timeout)
        except (OSError, subprocess.SubprocessError):
            return None
        if completed.returncode != 0:
            return None
        return completed.stdout.strip()

    branch = git('rev-parse', '--abbrev-ref', 'HEAD')
    commit = git('rev-parse', '--short', 'HEAD')
    # core.fileMode=false: the checkout lives on a share where the executable bit differs between
    # the Windows and the Linux view, and without this every file would count as modified.
    # A generous timeout: status walks the whole tree, which took 21 s over the share.
    status = git('-c', 'core.fileMode=false', 'status', '--porcelain', '--untracked-files=no',
                 timeout=120)

    return {
        'BRANCH': branch or UNKNOWN,
        'COMMIT': commit or UNKNOWN,
        'DIRTY': UNKNOWN if status is None else bool(status),
        'BUILT_AT': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
    }


def write(path: str, info: dict) -> None:
    with open(path, 'w') as stamp_file:
        stamp_file.write('# Written by setup.py at install time (see _build_stamp.py). Not tracked by git.\n')
        for field in FIELDS:
            stamp_file.write(f'{field} = {info.get(field, UNKNOWN)!r}\n')


def read() -> dict:
    """The stamp of the running install, all UNKNOWN when there is none.

    No stamp means an install from before stamping existed, or one not made through setup.py.
    """
    try:
        from structman import _build_info
    except ImportError:
        return {field: UNKNOWN for field in FIELDS}
    return {field: getattr(_build_info, field, UNKNOWN) for field in FIELDS}


def describe(version: str, info: dict, location: str) -> str:
    """One line naming the code a run is using, for the main log and the errorlog header."""
    dirty = info.get('DIRTY', UNKNOWN)
    if dirty is True:
        state = 'with uncommitted changes'
    elif dirty is False:
        state = 'clean'
    else:
        state = 'dirty state unknown'
    return (f'structman {version} built from {info.get("BRANCH", UNKNOWN)}@{info.get("COMMIT", UNKNOWN)} '
            f'({state}) at {info.get("BUILT_AT", UNKNOWN)}, running from {location}')


def stack_description(version_of=None, environ=None) -> str:
    """The numba stack and where its kernel cache lives.

    Checking out the other branch swaps numba and llvmlite on reinstall, and a kernel cache left
    by the previous stack is what killed the annotate workers with 'LLVM bitcode parsing error'.
    Named in the first lines of every run, such a swap is visible before anything fails.
    """
    if version_of is None:
        from importlib.metadata import version as version_of
    environ = os.environ if environ is None else environ

    def installed(name):
        try:
            return version_of(name)
        except Exception:
            return UNKNOWN

    cache = environ.get('NUMBA_CACHE_DIR') or "numba's default, next to the sources"
    return f'numba {installed("numba")}, llvmlite {installed("llvmlite")}, kernel cache {cache}'


def describe_running() -> str:
    """describe() for the install this process imported, with its numba stack."""
    try:
        from structman._version import __version__ as version
    except ImportError:
        version = UNKNOWN
    return (describe(version, read(), os.path.dirname(os.path.abspath(__file__)))
            + '; ' + stack_description())
