#!/bin/bash
trap "exit" INT

BASE_DIR="/scratch/SCRATCH_NVME/pdb"
SERVER=rsync.wwpdb.org::ftp
PORT=33444

if [[ ! -d "${BASE_DIR}" ]]
then
    mkdir -p "${BASE_DIR}"
    chmod 2775 "${BASE_DIR}"
fi

if [[ ! -d "${BASE_DIR}/data/structures/divided/pdb" ]]
then
    mkdir -p "${BASE_DIR}/data/structures/divided/pdb"
    chmod 2775 "${BASE_DIR}/data"
    chmod 2775 "${BASE_DIR}/data/structures"
    chmod 2775 "${BASE_DIR}/data/structures/divided"
    chmod 2775 "${BASE_DIR}/data/structures/divided/pdb"
    
fi
rsync --info=progress2 --no-inc-recursive -rt -v -z --perms --no-g --chmod=D2775,F664 --delete --port=${PORT} ${SERVER}/data/structures/divided/pdb/ "${BASE_DIR}/data/structures/divided/pdb/"

if [[ ! -d "${BASE_DIR}/data/biounit/PDB/divided" ]]
then
    mkdir -p "${BASE_DIR}/data/biounit/PDB/divided"
    chmod 2775 "${BASE_DIR}/data"
    chmod 2775 "${BASE_DIR}/data/biounit"
    chmod 2775 "${BASE_DIR}/data/biounit/PDB"
    chmod 2775 "${BASE_DIR}/data/biounit/PDB/divided"
fi
rsync --info=progress2 --no-inc-recursive -rt -v -z --perms --no-g --chmod=D2775,F664 --delete --port=${PORT} ${SERVER}/data/biounit/PDB/divided/ "${BASE_DIR}/data/biounit/PDB/divided/"

if [[ ! -d "${BASE_DIR}/data/structures/obsolete" ]]
then
    mkdir -p "${BASE_DIR}/data/structures/obsolete"
    chmod 2775 "${BASE_DIR}/data"
    chmod 2775 "${BASE_DIR}/data/structures"
    chmod 2775 "${BASE_DIR}/data/structures/obsolete"
fi
rsync --info=progress2 --no-inc-recursive -rt -v -z --perms --no-g --chmod=D2775,F664 --delete --port=${PORT} ${SERVER}/data/structures/obsolete/ "${BASE_DIR}/data/structures/obsolete/"

if [[ ! -d "${BASE_DIR}/data/status" ]]
then
    mkdir -p "${BASE_DIR}/data/status"
    chmod 2775 "${BASE_DIR}/data"
    chmod 2775 "${BASE_DIR}/data/status"
fi
rsync --info=progress2 --no-inc-recursive -rt -v -z --perms --no-g --chmod=D2775,F664 --delete --port=${PORT} ${SERVER}/data/status/ "${BASE_DIR}/data/status/"
