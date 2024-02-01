#!/usr/bin/python3
# curl --metalink ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz > uniprot_trembl.fasta.gz
# curl --metalink ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz > uniprot_sprot_varsplic.fasta.gz
# curl --metalink ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz > uniprot_sprot.fasta.gz

import getopt
import gzip
import os
import sys

import pymysql as MySQLdb

import structman
from structman.lib.database import database
from structman.base_utils import resolve_path


def parseSeqFile(seq_file_path, config):
    f = gzip.open(seq_file_path, 'r')
    lines = f.read().decode('ascii').split('\n')
    f.close()

    seq_map = {}
    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            u_ac = line.split('|')[1]
            seq_map[u_ac] = ''
        else:
            seq_map[u_ac] += (line)

    values = []
    for u_ac in seq_map:
        values.append((u_ac, seq_map[u_ac]))

    db = MySQLdb.connect(host=config.db_address, user=config.db_user_name, password=config.db_password, database=config.mapping_db)
    cursor = db.cursor()

    database.insert(db, cursor, 'Sequences', ['Uniprot_Ac', 'Sequence'], values)

    db.close()


if __name__ == "__main__":

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "c:i:f:", [])
    except getopt.GetoptError:
        print("Illegal Input")
        sys.exit(2)

    infolder = None
    infile = None
    for opt, arg in opts:
        if opt == '-c':
            config = arg
        if opt == '-f':
            infolder = arg
        if opt == '-i':
            infile = arg

    if infolder is None:
        #infolder = resolve_path(sys.argv[0]).rsplit('/',1)[0]
        pass  # for now
    else:
        infolder = resolve_path(infolder)

    config = resolve_path(config)

    config = structman.Config(config)

    if infile is not None:
        parseSeqFile(infile, config)

    if infolder is not None:
        seq_files = ['uniprot_sprot_varsplic.fasta.gz', 'uniprot_sprot.fasta.gz', 'uniprot_trembl.fasta.gz', ]

        for seq_file in seq_files:
            seq_file_path = f'{infolder}/{seq_file}'
            parseSeqFile(seq_file_path, config)
