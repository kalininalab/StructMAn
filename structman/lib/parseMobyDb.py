#!/usr/bin/python3
import getopt
import gzip
import json
import os
import sys

import pymysql as MySQLdb

import structman


def parseDB(db_path):

    f = gzip.open(db_path, 'r')
    lines = f.read().split(b'\n')
    f.close()

    entries = []
    for line in lines:
        str_line = line.decode('ascii')
        # print(str_line)
        if str_line == '':
            continue
        entries.append(json.loads(str_line))

    return entries


if __name__ == "__main__":

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "c:i:", [])
    except getopt.GetoptError:
        print("Illegal Input")
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-c':
            config = arg
        if opt == '-i':
            infile = arg

    infile = resolve_path(infile)
    config = resolve_path(config)

    mobiDB = parseDB(infile)

    config = resolve_path(config)
    config = structman.Config(config)

    score_values = []
    region_values = []

    db = MySQLdb.connect(host=config.db_address, user=config.db_user_name, password=config.db_password, database=config.mapping_db)
    cursor = db.cursor()

    u_acs = []
    n = 0
    m = 0

    # could be replaced with the update() function in database

    for entry in mobiDB:
        u_ac = entry['acc']
        u_acs.append(u_ac)
        predictors = entry['mobidb_consensus']['disorder']['predictors']
        for pred in predictors:
            method = pred['method']
            if method == 'mobidb-lite':
                scores = pred['scores']
                regions = pred['regions']

                score_values.append("WHEN '%s' THEN '%s'" % (u_ac, str(scores)))
                region_values.append("WHEN '%s' THEN '%s'" % (u_ac, str(regions).replace("'", "")))
        n += 1
        m += 1
        if n == 10000:
            print(m)
            sql = ("""UPDATE Sequences SET Disorder_Scores = CASE Uniprot_Ac %s ELSE Disorder_Scores END WHERE Uniprot_Ac IN ('%s')""") % (' '.join(score_values), "','".join(u_acs))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Couldn't insert Mapping: %s,%s" % (e, f))

            sql = ("""UPDATE Sequences SET Disorder_Regions = CASE Uniprot_Ac %s ELSE Disorder_Regions END WHERE Uniprot_Ac IN ('%s')""") % (' '.join(region_values), "','".join(u_acs))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Couldn't insert Mapping: %s,%s" % (e, f))
            u_acs = []
            n = 0
            score_values = []
            region_values = []

    sql = ("""UPDATE Sequences SET Disorder_Scores = CASE Uniprot_Ac %s ELSE Disorder_Scores END WHERE Uniprot_Ac IN ('%s')""") % (' '.join(score_values), "','".join(u_acs))
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Couldn't insert Mapping: %s,%s" % (e, f))

    sql = ("""UPDATE Sequences SET Disorder_Regions = CASE Uniprot_Ac %s ELSE Disorder_Regions END WHERE Uniprot_Ac IN ('%s')""") % (' '.join(region_values), "','".join(u_acs))
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Couldn't insert Mapping: %s,%s" % (e, f))

    db.close()
