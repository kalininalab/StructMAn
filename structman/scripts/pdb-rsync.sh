#!/bin/bash
# Deliberately no `set -e`: one unreachable directory should not stop the rest of
# the mirror from updating. Every rsync's exit status is collected instead and
# reported together at the end, and the script exits non-zero if any of them
# failed - without that, a transfer dying half way left a partially populated tree
# behind, the script walked on to the next directory and still exited 0, so
# update.py printed 'Update PDB done' over an incomplete mirror.
set -o pipefail

trap 'echo >&2; echo "INTERRUPTED - the mirror is incomplete, re-run to resume." >&2; exit 130' INT

FAILED_DIRS=()
SYNCED_DIRS=0

# Attempts per directory, and the base for the linear backoff between them
# (60s, then 120s). A dropped connection part way through a multi-hour transfer
# used to cost the whole directory; a retry picks up where it stopped.
MAX_ATTEMPTS=3
RETRY_BACKOFF=60

BASE_DIR="/scratch/SCRATCH_NVME/pdb"
PORT=33444

# Candidates tried in order, first one that answers wins. The hostname is
# listed first since it always tracks wherever wwPDB currently points;
# the IPs are fallbacks in case DNS resolution isn't available.
CANDIDATE_SERVERS=(
    "rsync.wwpdb.org::ftp"
    "128.6.158.50::ftp"
    "132.249.213.147::ftp"
    "128.6.159.100::ftp"
)

SERVER=""
for candidate in "${CANDIDATE_SERVERS[@]}"; do
    echo "Testing connection to ${candidate} (port ${PORT})..."
    if rsync --contimeout=10 --port=${PORT} --list-only "${candidate}/data/status/" > /dev/null 2>&1; then
        SERVER="${candidate}"
        echo "Connected to ${SERVER}"
        break
    else
        echo "Could not reach ${candidate}"
    fi
done

if [[ -z "${SERVER}" ]]; then
    echo "ERROR: could not connect to any known PDB rsync server, aborting." >&2
    exit 1
fi

if [[ ! -d "${BASE_DIR}" ]]
then
    mkdir -p "${BASE_DIR}"
    chmod 2775 "${BASE_DIR}"
fi

# rsync exit codes that do not mean the transfer failed:
#   0   everything transferred
#   24  some source files vanished mid-transfer. wwPDB rewrites entries
#       continuously, so on a mirror this size that is routine, and the next run
#       picks up whatever moved.
rsync_succeeded () {
    [[ "$1" -eq 0 || "$1" -eq 24 ]]
}

# Transport problems rather than anything wrong with the request, so the same
# command can simply be run again. Retrying resumes rather than restarts: whatever
# already landed matches on size and mtime and is skipped.
#   10  socket I/O error - the connection dropped. This is what cut the two mmCIF
#       trees short; they are the big ones, so they are the ones that stay
#       connected long enough to be hit by it.
#   11  file I/O error
#   12  error in the rsync protocol data stream
#   23  partial transfer due to error
#   30  timeout in data send/receive
#   35  timeout waiting for daemon connection
#   255 unexplained transport failure
# Deliberately NOT retried: 20 (killed by signal - the user meant it), and the
# 1-5 range, which are usage and protocol-negotiation errors that a retry repeats.
rsync_retryable () {
    case "$1" in
        10|11|12|23|30|35|255) return 0 ;;
        *) return 1 ;;
    esac
}

# Mirrors one directory of the wwPDB tree into the same relative location
# below BASE_DIR, creating (and chmod'ing) all missing parent directories.
# Returns non-zero and records the directory if its rsync did not complete; note
# that the target directory is created before the transfer, so a failure here
# leaves an empty or half-filled tree that looks exactly like a finished one.
sync_dir () {
    local rel_path="$1"
    local target="${BASE_DIR}/${rel_path}"
    local status

    if [[ ! -d "${target}" ]]
    then
        mkdir -p "${target}"
        local -a components
        local component
        local partial="${BASE_DIR}"
        IFS='/' read -r -a components <<< "${rel_path}"
        for component in "${components[@]}"
        do
            partial="${partial}/${component}"
            chmod 2775 "${partial}"
        done
    fi

    local attempt=1
    local backoff
    while true
    do
        if [[ "${attempt}" -eq 1 ]]
        then
            echo "Syncing ${rel_path} ..."
        else
            echo "Syncing ${rel_path} (attempt ${attempt}/${MAX_ATTEMPTS}) ..."
        fi

        rsync --info=progress2 --no-inc-recursive -rt -v -z --perms --no-g --chmod=D2775,F664 --delete --port=${PORT} ${SERVER}/${rel_path}/ "${target}/"
        status=$?

        if rsync_succeeded "${status}"
        then
            SYNCED_DIRS=$((SYNCED_DIRS + 1))
            echo "OK: ${rel_path} (rsync exit ${status}, attempt ${attempt})"
            return 0
        fi

        if ! rsync_retryable "${status}" || [[ "${attempt}" -ge "${MAX_ATTEMPTS}" ]]
        then
            break
        fi

        backoff=$((attempt * RETRY_BACKOFF))
        echo "WARNING: ${rel_path} failed (rsync exit ${status}), retrying in ${backoff}s ..." >&2
        sleep "${backoff}"
        attempt=$((attempt + 1))
    done

    echo "ERROR: ${rel_path} did NOT complete after ${attempt} attempt(s) (rsync exit ${status})" >&2
    FAILED_DIRS+=("${rel_path} (rsync exit ${status} after ${attempt} attempt(s))")
    return 1
}

# Asymmetric units, both in the legacy PDB format and as mmCIF. The mmCIF
# copies also cover the entries that are too large for the PDB format and
# therefore have no .ent.gz counterpart.
sync_dir "data/structures/divided/pdb"
sync_dir "data/structures/divided/mmCIF"

# Biological assemblies, again in both formats. The mmCIF assemblies are not
# below data/biounit/, wwPDB keeps them in their own tree and names them
# <pdb_id>-assembly<n>.cif.gz instead of <pdb_id>.pdb<n>.gz.
sync_dir "data/biounit/PDB/divided"
sync_dir "data/assemblies/mmCIF/divided"

# Obsolete entries; this directory holds both the pdb/ and the mmCIF/ subtree.
sync_dir "data/structures/obsolete"

sync_dir "data/status"

echo
if [[ ${#FAILED_DIRS[@]} -gt 0 ]]
then
    echo "===============================================================" >&2
    echo "ERROR: ${#FAILED_DIRS[@]} directory/directories did not finish:" >&2
    for failed in "${FAILED_DIRS[@]}"
    do
        echo "  - ${failed}" >&2
    done
    echo >&2
    echo "The local PDB mirror is INCOMPLETE. Re-run this script to resume:" >&2
    echo "rsync skips what is already there, so it continues where it stopped" >&2
    echo "rather than starting over." >&2
    echo "===============================================================" >&2
    exit 1
fi

echo "All ${SYNCED_DIRS} directories synced."
