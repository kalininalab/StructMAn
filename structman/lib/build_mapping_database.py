#!/usr/bin/python3
import getopt
import gzip
import os
import sys

import pymysql as MySQLdb
from structman.base_utils import resolve_path


def parseConfig(config):

    global db_address
    global db_user_name
    global db_password

    f = open(config, "r")
    lines = f.read().split('\n')
    f.close()
    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '#':
            continue
        words = line.split("=")
        if len(words) > 1:
            opt = words[0]
            arg = words[1].replace("\n", "")

            if opt == 'db_address':
                db_address = arg
            elif opt == 'db_user_name':
                db_user_name = arg
            elif opt == 'db_password':
                db_password = arg


if __name__ == "__main__":

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "c:i:d:", [])
    except getopt.GetoptError:
        print("Illegal Input")
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-c':
            config = arg
        if opt == '-i':
            infile = arg
        if opt == '-d':
            db_name = arg

    infile = resolve_path(infile)
    config = resolve_path(config)

    parseConfig(config)

    if infile[-3:] == '.gz':
        f = gzip.open(infile, 'r')
        lines = f.readlines()
        f.close()
    else:
        f = open(infile, 'r')
        lines = f.readlines()
        f.close()

    ac_refs = []
    ac_ids = []
    for line in lines:
        words = line[:-1].split()
        u_ac = words[0]
        id_name = words[1]
        id_value = words[2]

        if id_name == 'UniProtKB-ID':
            ac_ids.append("('%s','%s','%s','%s')" % (u_ac, u_ac[-2:], id_value, id_value[:2]))
        if id_name == 'RefSeq':
            ac_refs.append("('%s','%s','%s','%s')" % (u_ac, u_ac.split('-')[0][-2:], id_value, id_value[:2]))

    current_ids = {'UniProtKB-ID': 'Uniprot_Id', 'RefSeq': 'Refseq_Prot', 'RefSeq_NT': 'Refseq_Nucl'}
    current_id_list = ['UniProtKB-ID', 'RefSeq', 'RefSeq_NT']

    db = MySQLdb.connect(db_address, db_user_name, db_password, db_name)
    cursor = db.cursor()

    size = len(ac_refs)
    print(size)
    if size > 0:
        bound = 100000
        if size > bound:
            m = size / bound
            rest = size % bound
            if rest == 0:
                part_size = size / m
            else:
                m += 1
                part_size = size / m
            parts = []
            for i in range(0, m - 1):
                parts.append(ac_refs[(i * part_size):((i + 1) * part_size)])
            parts.append(ac_refs[(m - 1) * part_size:])
        else:
            parts = [ac_refs]
        for part in parts:

            sql = ("""INSERT IGNORE INTO AC_Refseq(Uniprot_Ac,Ac_Hash,Refseq,Refseq_Hash)
                    VALUES %s""") % (','.join(part))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Couldn't insert Mapping: %s,%s" % (e, f))

    size = len(ac_ids)
    print(size)
    if size > 0:
        bound = 300000
        if size > bound:
            m = size / bound
            rest = size % bound
            if rest == 0:
                part_size = size / m
            else:
                m += 1
                part_size = size / m
            parts = []
            for i in range(0, m - 1):
                parts.append(ac_ids[(i * part_size):((i + 1) * part_size)])
            parts.append(ac_ids[(m - 1) * part_size:])
        else:
            parts = [ac_ids]
        for part in parts:

            sql = ("""INSERT IGNORE INTO AC_ID(Uniprot_Ac,Ac_Hash,Uniprot_Id,Id_Hash)
                    VALUES %s""") % (','.join(part))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Couldn't insert Mapping: %s,%s" % (e, f))

    db.close()
