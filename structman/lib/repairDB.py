#!/usr/bin/python3
import getopt
import gzip
import os
import subprocess
import sys
import traceback
import sqlite3

import structman
from structman.lib.database.database_core_functions import binningSelect, select
from structman import _version

# Tries to reclassify all positions with no classification in the database
# Usefull when,
#   - There was an error in insertClassifications
def reclass_null(config):
    db = MySQLdb.connect(config.db_address, config.db_user_name, config.db_password, config.db_name)
    cursor = db.cursor()
    # get all positions without classifications (together with their protein id)
    columns = ['Position_Id', 'Gene']
    table = 'Position'
    null_columns = set(['Class'])
    results = select(db, cursor, columns, table, null_columns=null_columns)

    prot_ids = set()
    for row in results:
        prot_ids.add(row[1])

    # for all proteins search for alignments in the database

    rows = ['Gene', 'Structure', 'Alignment']
    table = 'Alignment'
    results = binningSelect(prot_ids, rows, table, db, cursor)

    # try to map the positions into the alignments
    # for all positions, which could be mapped, reclassify


def remove_sessions(config):
    db, cursor = config.getDB()

    sqlite = config.db_address == '-'

    if sqlite:
        sql_commands = ['PRAGMA foreign_keys = 0;']
        truncate_command = "DELETE FROM"
    else:
        sql_commands =  ['SET FOREIGN_KEY_CHECKS=0;']
        truncate_command = "TRUNCATE"

    sql_commands += [
                    f'{truncate_command} Session;',
                    f'{truncate_command} RS_Protein_Session;',
                    f'{truncate_command} RS_Position_Session;',
                    f'{truncate_command} RS_Indel_Session;',
                    f'{truncate_command} RS_SNV_Session;',
                    f'{truncate_command} RS_Multi_Mutation_Session;'
                    ]

    if sqlite:
        sql_commands.append('PRAGMA foreign_keys = 1;')
    else:
        sql_commands.append('SET FOREIGN_KEY_CHECKS=1;')

    for sql in sql_commands:
        try:
            # Execute the SQL command
            cursor.execute(sql)
            # Commit your changes in the database
            # db.commit()
        except:
            # Rollback in case there is any error
            # db.rollback()
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc(g)
            print("Error: ", e, f, g)

    db.close()
    return


def reset(config, keep_structures=False):
    db, cursor = config.getDB()
    sqlite = (config.db_address == '-')
    if sqlite:
        sql_commands = ['PRAGMA foreign_keys = 0;']
        truncate_command = "DELETE FROM"
    else:
        sql_commands =  ['SET FOREIGN_KEY_CHECKS=0;']
        truncate_command = "TRUNCATE"

    tables =    [
                'Protein', 'Gene', 'Position', 'Interface', 'SNV', 'Multi_Mutation',
                'RS_Isoform', 'Indel', 'GO_Term', 'Pathway', 'Session', 'Alignment',
                'Position_Position_Interaction', 'Protein_Protein_Interaction', 
                'RS_Protein_Session', 'RS_Protein_GO_Term', 'RS_Position_Session',
                'RS_Position_Interface', 'RS_SNV_Session', 'RS_Multi_Mutation_Session', 
                'RS_Indel_Session', 'RS_Protein_Pathway'
                ]

    if not keep_structures:
        tables += [
            'Ligand', 'Structure', 'Residue', 'Complex', 'RS_Ligand_Structure',
            'RS_Residue_Residue', 'RS_Residue_Interface'
        ]

    for table in tables:
        if sqlite:
            sql_commands.append(f'{truncate_command} {table};')
            sql_commands.append(f"UPDATE SQLITE_SEQUENCE SET SEQ=0 WHERE NAME='{table}';")
        else:
            sql_commands.append(f'{truncate_command} {table};')

    if sqlite:
        sql_commands.append('PRAGMA foreign_keys = 1;')
    else:
        sql_commands.append('SET FOREIGN_KEY_CHECKS=1;')

    for sql in sql_commands:
        try:
            # Execute the SQL command
            cursor.execute(sql)
            # Commit your changes in the database
            db.commit()
        except:
            # Rollback in case there is any error
            # db.rollback()
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            print("Error: ", e, f, g)

    db.close()


def reduceToStructures(config, infile):
    f = open(infile, 'r')
    lines = f.readlines()
    f.close()

    complexes = set()
    for line in lines:
        words = line.replace('\t', ' ').split()
        pdb_tuple = words[0]
        if pdb_tuple.count(':') != 1:
            continue
        pdb_id = pdb_tuple.split(':')[0]
        complexes.add(pdb_id)

    table = 'Structure'
    rows = ['Structure_Id', 'PDB']
    results = select(config, rows, table)

    list_of_doom = []
    for row in results:
        if row[1] in complexes:
            continue
        list_of_doom.append(row[0])

    results = select(config, ['Complex_Id', 'PDB'], 'Complex')

    list_of_doom_c = []
    for row in results:
        if row[1] in complexes:
            continue
        list_of_doom_c.append(row[0])

    statement = 'DELETE FROM Residue WHERE Structure IN (%s)' % (','.join(['%s'] * len(list_of_doom)))
    db, cursor = config.getDB()

    cursor.execute(statement, list_of_doom)
    results = cursor.fetchall()
    db.commit()

    db.close()

    statement = 'DELETE FROM Structure WHERE Structure_Id IN (%s)' % (','.join(['%s'] * len(list_of_doom)))
    db, cursor = config.getDB()

    cursor.execute(statement, list_of_doom)
    results = cursor.fetchall()
    db.commit()

    db.close()

    statement = 'DELETE FROM Complex WHERE Complex_Id IN (%s)' % (','.join(['%s'] * len(list_of_doom_c)))
    db, cursor = config.getDB()

    cursor.execute(statement, list_of_doom_c)
    results = cursor.fetchall()
    db.commit()

    db.close()


def export(config, target_path):
    target_path = '%s/%s.sql.gz' % (resolve_path(target_path), config.db_name)
    cmds = ' '.join(['mysqldump', '-n', '-u', config.db_user_name, '-h', config.db_address, '--password=%s' % config.db_password, '--single-transaction', '--databases', config.db_name,
                     '|', 'pigz', '--best', '-p', '8', '>', '"%s"' % target_path])

    p = subprocess.Popen(cmds, shell=True)
    p.wait()


def empty(config):
    reset(config)

def clear(config):
    reset(config, keep_structures=True)

def destroy(config):
    try:
        if config.db_address == '-':
            os.remove(config.local_db_path)
        else:
            db, cursor = config.getDB()
            sql = 'DROP DATABASE `%s`' % config.db_name
            cursor.execute(sql)
            db.close()
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        errortext = '\n'.join([str(e), str(f), str(g)]) + '\n\n'
        config.errorlog.add_error(errortext)
        raise NameError("Error in destroy: %s" % sql)

def insert_meta_data(config):
    db, cursor = config.getDB()
    sql = f'Insert INTO Database_Metadata (StructMAn_Version, PPI_Feature) VALUES ("{_version.__version__}", "{int(config.compute_ppi)}");'

    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        print('\n'.join([str(e), str(f), str(g)]))
        db.close()
    db.close()

def load(config):
    if config.db_address != '-':
        try:
            db, cursor = config.getDB(silent = True)
            db.close()
            print(f'Database with name {config.db_name} already exists. If you want to recreate it, please call [structman database destroy] first.')
            return
        except:
            if config.verbosity >= 3:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                print('\n'.join([str(e), str(f), str(g)]))
            pass
        db, cursor = config.getDB(server_connection=True)
    else:
        if os.path.exists(config.local_db_path):
            print(f'Database with path {config.local_db_path} already exists. If you want to recreate it, please call [structman database destroy] first.')
            return
        else:
            db = sqlite3.connect(config.local_db_path)
            cursor = db.cursor()


    if config.db_address != '-':
        sql = 'CREATE DATABASE `%s`' % config.db_name

        try:
            cursor.execute(sql)
            #db.close()
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            print(f'SQL command: [{sql}] failed:\n{e}\n{f}\n{g}')
            try:
                db.close()
            except:
                pass

        if config.database_source_path[-3:] == '.gz':
            f = gzip.open(config.database_source_path, 'rb')
            binary = True
        else:
            f = open(config.database_source_path, 'r')
            binary = False
        lines = f.readlines()
        f.close()

        sql_commands = []
        new_command = ''
        for line in lines:
            if binary:
                line = line.decode('ascii')
            line = line.strip()
            if line[:4] == 'USE ':
                sql_commands.append('USE `%s`;' % config.db_name)
            else:
                if line == '':
                    continue
                if line[0] == '-':
                    continue
                if line[:1] == '/*':
                    continue
                new_command += line
                if line[-1] == ';':
                    sql_commands.append(new_command)
                    new_command = ''

        for sql_command in sql_commands:

            try:
                cursor.execute(sql_command)
            except:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                print(f'SQL command: [{sql_command}] failed:\n{e}\n{f}\n{g}')
                try:
                    db.close()
                except:
                    pass
    else:
        if config.sqlite_database_source_path[-3:] == '.gz':
            f = gzip.open(config.sqlite_database_source_path, 'rb')
            binary = True
        else:
            f = open(config.sqlite_database_source_path, 'r')
            binary = False
        lines = f.readlines()
        f.close()

        sql_commands = []
        new_command = ''
        for line in lines:
            if binary:
                line = line.decode('ascii')
            line = line.strip()

            if line == '':
                continue
            if line[0] == '-':
                continue
            if line[:1] == '/*':
                continue
            new_command += line
            if line[-1] == ';':
                sql_commands.append(new_command)
                new_command = ''

        for sql_command in sql_commands:

            try:
                cursor.execute(sql_command)
            except:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                print('\n'.join([str(e), str(f), str(g), sql_command]))
                db.close()


    db.close()

    insert_meta_data(config)


def check_structman_version(config):

    meta = select(config, ['StructMAn_Version', 'PPI_Feature'], 'Database_Metadata')

    if len(meta) == 0:
        insert_meta_data(config)
        meta = select(config, ['StructMAn_Version', 'PPI_Feature'], 'Database_Metadata')

    version_in_db, ppi_feature = meta[0]

    if version_in_db != _version.__version__:
        major_version_number, database_version_number, minor_version_number = _version.__version__.split('.')
        major_version_number_in_db, database_version_number_in_db, minor_version_number_in_db = version_in_db.split('.')

        if major_version_number != major_version_number_in_db:
            print(f'\n\n!!!CRITICAL WARNING!!!\nMajor version number ({major_version_number}) does not match with the major version number stored in the database ({major_version_number_in_db}).\nDatabase will most likely not be compatible!\n\n')
        elif database_version_number != database_version_number_in_db:
            print(f'\nWARNING\nVersion number ({_version.__version__}) does not match with the version number stored in the database ({version_in_db}).\nDatabase might not be compatible!\n')

    if not ppi_feature and config.compute_ppi:
        print('\nWARNING\nThe given database does not support the PPI feature and will not be computed.\n')
    elif ppi_feature and not config.compute_ppi:
        config.compute_ppi = True

# destroys and reinitializes the database
def reinit(config):
    empty(config)
    destroy(config)
    load(config)


def main(config, reclassify_null=False):
    # called with --rcn
    if reclassify_null:
        reclass_null(config)


if __name__ == "__main__":
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "c:", ['rcn'])
    except getopt.GetoptError:
        print("Illegal Input")
        sys.exit(2)

    reclassify_null = False
    config_path = None

    for opt, arg in opts:
        if opt == '-c':
            config_path = arg
        elif opt == '--rcn':
            reclassify_null = True

    config = structman.Config(config_path, None, None, None, None, False)

    main(config, reclassify_null=reclassify_null)
