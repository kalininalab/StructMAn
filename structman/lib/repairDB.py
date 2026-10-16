#!/usr/bin/python3
import getopt
import gzip
import os
import subprocess
import sys
import traceback
import sqlite3
import time
from concurrent.futures import ThreadPoolExecutor
import ray

from structman.lib.database.database_core_functions import binningSelect, select, remove, insert
#from structman.base_utils.config_class import Config

from structman import _version

table_names =    [
            'Protein', 'Gene',
            'Position_0',
            'Position_1',
            'Position_2',
            'Position_3',
            'Position_4',
            'Position_5',
            'Position_6',
            'Position_7',
            'Position_8',
            'Position_9',
            'Position_10',
            'Position_11',
            'Position_12',
            'Position_13',
            'Position_14',
            'Position_15',
            'Position_16',
            'Position_17',
            'Position_18',
            'Position_19',
            'Interface', 'SNV', 'Multi_Mutation',
            'RS_Isoform', 'Indel', 'GO_Term', 'Pathway', 'Session',
            'Alignment_0',
            'Alignment_1',
            'Alignment_2',
            'Alignment_3',
            'Alignment_4',
            'Alignment_5',
            'Alignment_6',
            'Alignment_7',
            'Alignment_8',
            'Alignment_9',
            'Alignment_10',
            'Alignment_11',
            'Alignment_12',
            'Alignment_13',
            'Alignment_14',
            'Alignment_15',
            'Alignment_16',
            'Alignment_17',
            'Alignment_18',
            'Alignment_19',
            'Protein_Protein_Interaction_0', 
            'Protein_Protein_Interaction_1', 
            'Protein_Protein_Interaction_2', 
            'Protein_Protein_Interaction_3', 
            'Protein_Protein_Interaction_4', 
            'Protein_Protein_Interaction_5', 
            'Protein_Protein_Interaction_6', 
            'Protein_Protein_Interaction_7', 
            'Protein_Protein_Interaction_8', 
            'Protein_Protein_Interaction_9',
            'Protein_Protein_Interaction_10', 
            'Protein_Protein_Interaction_11', 
            'Protein_Protein_Interaction_12', 
            'Protein_Protein_Interaction_13', 
            'Protein_Protein_Interaction_14', 
            'Protein_Protein_Interaction_15', 
            'Protein_Protein_Interaction_16', 
            'Protein_Protein_Interaction_17', 
            'Protein_Protein_Interaction_18', 
            'Protein_Protein_Interaction_19', 
            'RS_Protein_Session', 'RS_Protein_GO_Term', 'RS_Position_Session',
            'RS_Position_Interface', 'RS_SNV_Session', 'RS_Multi_Mutation_Session', 
            'RS_Indel_Session', 'RS_Protein_Pathway'
            ]

structure_table_names = [
        'Ligand',
        'Structure_0',
        'Structure_1',
        'Structure_2',
        'Structure_3',
        'Structure_4',
        'Structure_5',
        'Structure_6',
        'Structure_7',
        'Structure_8',
        'Structure_9',
        'Structure_10',
        'Structure_11',
        'Structure_12',
        'Structure_13',
        'Structure_14',
        'Structure_15',
        'Structure_16',
        'Structure_17',
        'Structure_18',
        'Structure_19',
        'Complex', 'RS_Ligand_Structure',
    ]

session_tables = [
    'Session', 
    'RS_Protein_Session', 'RS_Position_Session',
    'RS_SNV_Session', 'RS_Multi_Mutation_Session', 
    'RS_Indel_Session'
]

table_names_in_order = [
            [
                'Gene', 'Complex',
                'Structure_0',
            'Structure_1',
            'Structure_2',
            'Structure_3',
            'Structure_4',
            'Structure_5',
            'Structure_6',
            'Structure_7',
            'Structure_8',
            'Structure_9',
            'Structure_10',
            'Structure_11',
            'Structure_12',
            'Structure_13',
            'Structure_14',
            'Structure_15',
            'Structure_16',
            'Structure_17',
            'Structure_18',
            'Structure_19',
            'Session', 'GO_Term', 'Pathway', 'Ligand', 'Database_Metadata'
            ],
            [
                'Protein', 'Interface', 'RS_Ligand_Structure', 
            ],
            [
                'Alignment_0',
                'Alignment_1',
                'Alignment_2',
                'Alignment_3',
                'Alignment_4',
                'Alignment_5',
                'Alignment_6',
                'Alignment_7',
                'Alignment_8',
                'Alignment_9',
                'Alignment_10',
                'Alignment_11',
                'Alignment_12',
                'Alignment_13',
                'Alignment_14',
                'Alignment_15',
                'Alignment_16',
                'Alignment_17',
                'Alignment_18',
                'Alignment_19',
                'Position_0',
                'Position_1',
                'Position_2',
                'Position_3',
                'Position_4',
                'Position_5',
                'Position_6',
                'Position_7',
                'Position_8',
                'Position_9',
                'Position_10',
                'Position_11',
                'Position_12',
                'Position_13',
                'Position_14',
                'Position_15',
                'Position_16',
                'Position_17',
                'Position_18',
                'Position_19',
                'Multi_Mutation',
                'Protein_Protein_Interaction_0', 
                'Protein_Protein_Interaction_1', 
                'Protein_Protein_Interaction_2', 
                'Protein_Protein_Interaction_3', 
                'Protein_Protein_Interaction_4', 
                'Protein_Protein_Interaction_5', 
                'Protein_Protein_Interaction_6', 
                'Protein_Protein_Interaction_7', 
                'Protein_Protein_Interaction_8', 
                'Protein_Protein_Interaction_9',
                'Protein_Protein_Interaction_10', 
                'Protein_Protein_Interaction_11', 
                'Protein_Protein_Interaction_12', 
                'Protein_Protein_Interaction_13', 
                'Protein_Protein_Interaction_14', 
                'Protein_Protein_Interaction_15', 
                'Protein_Protein_Interaction_16', 
                'Protein_Protein_Interaction_17', 
                'Protein_Protein_Interaction_18', 
                'Protein_Protein_Interaction_19',
                'RS_Protein_Session', 'Indel',
                'RS_Protein_GO_Term', 'RS_Protein_Pathway', 
            ],
            [
                'SNV',
                'RS_Position_Session',
                'RS_Multi_Mutation_Session', 
                'RS_Indel_Session', 'RS_Position_Interface',
                'RS_Isoform'
            ],
            [
                'RS_SNV_Session'
            ]
        ]


def remove_sessions(config):
    db, cursor = config.getDB()

    sqlite = config.db_address == '-'

    if sqlite:
        sql_commands = ['PRAGMA foreign_keys = 0;']
        truncate_command = "DELETE FROM"
    else:
        sql_commands =  ['SET FOREIGN_KEY_CHECKS=0;']
        truncate_command = "TRUNCATE"


    for table_name in session_tables:
        sql_commands.append(f'{truncate_command} {table_name};')

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
            g = traceback.format_exc()
            print("Error: ", e, f, g)

    # TRUNCATE commits on its own in MySQL, but SQLite's DELETE FROM is rolled back on close.
    db.commit()
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

    tables = table_names

    if not keep_structures:
        tables += structure_table_names

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


def get_base_info(cursor):
    cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_allowed_packet'")
    max_package_size = cursor.fetchone()[1]

    cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_connections'")
    max_connections = cursor.fetchone()[1]
    return max_package_size, max_connections


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


def remove_empty_structures(config):
    results = select(config, ['Structure_Id'], 'Structure')

    stored_structures = []
    for row in results:
        stored_structures.append(row[0])

    results = select(config, ['Structure'], 'Residue')
    structures_with_residues = set()

    for row in results:
        structures_with_residues.add(row[0])

    empty_structures = []
    for structure_id in stored_structures:
        if structure_id not in structures_with_residues:
            empty_structures.append(structure_id)

    remove(config, 'Structure', in_rows={'Structure_Id' : empty_structures})


def fix_database(config):
    remove_empty_structures(config)


def reroll(config):
    results = select(config, ['Session_Id'], 'Session')

    last_session = None
    for row in results:
        if last_session is None:
            last_session = row[0]
        elif last_session < row[0]:
            last_session = row[0]

@ray.remote(max_calls=1)
def r_clone_table(config: any, table_name: str, target_db: str, template_db: str) -> None:
    clone_table_v2(config, table_name, target_db, template_db)

def clone_table(config, table_name, target_db):
    t0 = time.time()
    print(f'Cloning {table_name}')

    sql_cmd: str = f"SHOW COLUMNS FROM {table_name};"

    db, cursor = config.getDB()
   
    try:
        cursor.execute(sql_cmd)
        results: tuple[tuple] = cursor.fetchall()
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        print('\n'.join([str(e), str(f), str(g)]))
        db.close()
    db.close()

    column_names: list[str] = []
    for row in results:
        column_names.append(row[0])

    results: tuple[tuple] = select(config, column_names, table_name)

    t1 = time.time()
    print(f'  Time for select: {t1-t0} {len(results)}')

    insert(table_name, column_names, results, config, db_name = target_db)

    t2 = time.time()
    print(f'  Time for insert: {t2-t1}')

def clone_table_v2(config, table_name, target_db, template_db):
    t0 = time.time()
    print(f'Cloning {table_name}')
    db, cursor = config.getDB(server_connection=True)

    sql = f'CREATE Table `{target_db}`.{table_name} LIKE `{template_db}`.{table_name};'

    cursor.execute(sql)
    db.commit()
    t1 = time.time()
    print(f'  Time for creating structure of {table_name}: {t1-t0}')

    sql = f'INSERT INTO `{target_db}`.{table_name} SELECT * FROM `{template_db}`.{table_name};'

    cursor.execute(sql)
    db.commit()
    t2 = time.time()
    print(f'  Time for copying data of {table_name}: {t2-t1}')

    db.close()

def clone_afdb(config, target_db_name: str, template_db: str | None = None):
    if template_db is None:
        template_db = config.afdb

    db, cursor = config.getDB(server_connection=True, db_name=template_db)
    cursor.execute(f"SHOW TABLES FROM {template_db}")

    tables = [[], [], [], []]
    alternate = 0
    for row in cursor.fetchall():
        tables[alternate].append(row[0])
        alternate += 1
        if alternate == 4:
            alternate = 0

    db.close()

    clone_db(config, target_db_name, tables = tables, template_db=template_db)

def clone_db(config, target_db_name: str, tables: None | list[list[str]] = None, template_db: str | None = None) -> None:
    t0 = time.time()
    if tables is None:
        tables: list[list[str]] = table_names_in_order

    if template_db is None:
        template_db = config.db_name

    print('============ Database clone ===============')
    print(f'From {template_db} to {target_db_name}\n')

    # Standalone entry point, so ray is usually down here; an unguarded put would
    # auto-start an unconfigured cluster (see structman/__init__.py). Imported here
    # rather than at module level: config_class imports this module, so a top level
    # ray_utils import would close the circle config_class -> repairDB -> ray_utils.
    from structman.base_utils.ray_utils import ray_init
    ray_init(config, stop_existing = False)
    conf_dump = ray.put(config)

    db, cursor = config.getDB(server_connection=True)
    sql = f'CREATE DATABASE `{target_db_name}`'

    cursor.execute(sql)
    db.commit()

    db.close()

    for table_names in tables:
        clone_proc_ids = []
        for table_name in table_names:
            clone_proc_ids.append(r_clone_table.remote(conf_dump, table_name, target_db_name, template_db))
            #clone_table_v2(config, table_name, target_db_name)

        ray.get(clone_proc_ids)

    t1 = time.time()
    print(f'Total time for cloning: {t1-t0}')

# ---------------------------------------------------------------------------------------------
# database repair: remove the proteins a crashed run left unfinished
# ---------------------------------------------------------------------------------------------
#
# A chunk that dies after protCheck leaves its new proteins stored but unfinished, and a rerun
# cannot repair that by itself: protCheck finds the Protein row, positionCheck then marks the
# protein completely stored, and alignment and classification both skip stored proteins - so it
# stays unclassified for good. Removing exactly those proteins, with every row that hangs off
# them, turns the rerun into a resume: they come back as new and get computed, while everything
# that finished is found stored and skipped.
#
# What tells them apart: classification writes Position_Data for every position of every protein
# it classifies (an uncovered position gets a packed (None, None), not NULL), and all alignments
# are written before structural analysis begins. A protein with alignments and a position without
# Position_Data was therefore cut off in structural analysis or classification. Proteins without
# templates have no alignments and are never flagged. A crash before the alignments were written
# leaves proteins that look exactly like template-less ones; for that, from_protein_id names the
# chunk's first protein id (protCheck logs it), and every protein from there on that did not
# finish goes too.

REPAIR_SHARDS = 20
REPAIR_BATCH_SIZE = 500


def _wildcard(config) -> str:
    return '?' if config.db_address == '-' else '%s'


def _in_size(config) -> int:
    # sqlite caps the parameters of one statement (999 before 3.32).
    return 900 if config.db_address == '-' else 5000


def _connect(config):
    db, cursor = config.getDB()
    if cursor is None:
        if db is not None:
            db.close()
        raise ConnectionError(f'No database connection ({config.db_address}, {config.db_name})')
    return db, cursor


def _select_in(cursor, config, columns: str, table: str, key: str, values, condition: str = '') -> list:
    """SELECT columns FROM table WHERE key IN (values) [condition], in parts of _in_size."""
    values = sorted(set(values))
    wildcard = _wildcard(config)
    size = _in_size(config)
    rows = []
    for start in range(0, len(values), size):
        part = values[start:start + size]
        cursor.execute(
            f'SELECT {columns} FROM {table} WHERE {key} IN ({",".join([wildcard] * len(part))}){condition}', part)
        rows.extend(cursor.fetchall())
    return rows


def _delete_in(cursor, config, table: str, key: str, values) -> int:
    """DELETE FROM table WHERE key IN (values), in parts of _in_size. Returns the rows deleted."""
    values = sorted(set(values))
    wildcard = _wildcard(config)
    size = _in_size(config)
    deleted = 0
    for start in range(0, len(values), size):
        part = values[start:start + size]
        cursor.execute(f'DELETE FROM {table} WHERE {key} IN ({",".join([wildcard] * len(part))})', part)
        deleted += max(cursor.rowcount, 0)
    return deleted


def _existing_tables(config) -> set[str]:
    db, cursor = _connect(config)
    try:
        if config.db_address == '-':
            cursor.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
        else:
            cursor.execute('SHOW TABLES')
        return {row[0] for row in cursor.fetchall()}
    finally:
        db.close()


def _all_sessions(config) -> list[tuple[int, str]]:
    db, cursor = _connect(config)
    try:
        cursor.execute('SELECT Session_Id, Input_File FROM Session')
        return sorted((row[0], row[1]) for row in cursor.fetchall())
    finally:
        db.close()


def _by_shard(prot_ids) -> dict[int, list[int]]:
    shards: dict[int, list[int]] = {}
    for prot_id in prot_ids:
        shards.setdefault(prot_id % REPAIR_SHARDS, []).append(prot_id)
    return shards


def _table_label(table: str) -> str:
    # The 20 shards of a table are reported as one.
    head, _, tail = table.rpartition('_')
    return f'{head}_n' if tail.isdigit() else table


def sessions_of_proteins(config) -> list[tuple[int, int, int, int, str | None]]:
    """(session id, proteins, lowest and highest Protein_Id, input file) per session that first
    stored proteins, from Protein.Original_Session. That column outlives the session tables, so
    this still names the sessions after a truncation; their input file is then None. It is what
    --session_id takes once -i can no longer find a session."""
    db, cursor = _connect(config)
    try:
        cursor.execute('SELECT Original_Session, COUNT(*), MIN(Protein_Id), MAX(Protein_Id) FROM Protein '
                       'WHERE Original_Session IS NOT NULL GROUP BY Original_Session')
        per_session = [tuple(int(value) for value in row) for row in cursor.fetchall()]
        cursor.execute('SELECT Session_Id, Input_File FROM Session')
        input_files = {row[0]: row[1] for row in cursor.fetchall()}
    finally:
        db.close()
    return sorted(row + (input_files.get(row[0]),) for row in per_session)


def print_sessions_of_proteins(config) -> None:
    rows = sessions_of_proteins(config)
    if len(rows) == 0:
        print('No protein in the database names the session that stored it.')
        return
    print('Sessions that stored proteins (Protein.Original_Session):')
    for session, n, lowest, highest, input_file in rows:
        source = input_file if input_file is not None else 'no Session row any more'
        print(f'  session {session}: {n} proteins, Protein_Id {lowest} .. {highest}, {source}')
    print('Session ids start again at 1 once the session tables are truncated, so one id can name an old '
          'and a new run: add --from_protein_id to keep to one of them.')


def _candidate_protein_ids(config, session_id, from_protein_id, original_session=None) -> list[int]:
    wildcard = _wildcard(config)
    db, cursor = _connect(config)
    try:
        if original_session is not None:
            cursor.execute(f'SELECT Protein_Id FROM Protein WHERE Original_Session = {wildcard}', (original_session,))
        elif session_id is not None:
            cursor.execute(f'SELECT DISTINCT Protein FROM RS_Protein_Session WHERE Session = {wildcard}', (session_id,))
        elif from_protein_id is not None:
            cursor.execute(f'SELECT Protein_Id FROM Protein WHERE Protein_Id >= {wildcard}', (from_protein_id,))
        else:
            cursor.execute('SELECT Protein_Id FROM Protein')
        prot_ids = {int(row[0]) for row in cursor.fetchall()}
    finally:
        db.close()
    if from_protein_id is not None:
        prot_ids = {prot_id for prot_id in prot_ids if prot_id >= from_protein_id}
    return sorted(prot_ids)


def _scan_shard(config, shard: int, prot_ids: list[int]) -> tuple[set[int], set[int]]:
    """Of one shard's proteins: those with alignments, and those of them that have a position
    without Position_Data."""
    db, cursor = _connect(config)
    try:
        aligned = {row[0] for row in _select_in(
            cursor, config, 'DISTINCT Protein', f'Alignment_{shard}', 'Protein', prot_ids)}
        unfinished = {row[0] for row in _select_in(
            cursor, config, 'DISTINCT Protein', f'Position_{shard}', 'Protein', aligned,
            condition=' AND Position_Data IS NULL')}
    finally:
        db.close()
    return aligned, unfinished


def find_unfinished_proteins(config, session_id=None, from_protein_id=None,
                             original_session=None) -> tuple[list[int], int]:
    """The proteins a crashed run left unfinished, and how many proteins were checked.

    Without from_protein_id or original_session: every protein with alignments and a position
    without Position_Data. With either: every protein from that id on, or stored first by that
    session, that is not finished - which also takes the ones whose chunk died before their
    alignments were written.
    """
    candidates = _candidate_protein_ids(config, session_id, from_protein_id, original_session)
    jobs = sorted(_by_shard(candidates).items())

    # Reads only, one connection per shard. sqlite gets a single thread: its file lock would
    # serialise them anyway.
    if config.db_address == '-' or len(jobs) <= 1:
        results = [_scan_shard(config, shard, prot_ids) for shard, prot_ids in jobs]
    else:
        with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
            results = list(pool.map(lambda job: _scan_shard(config, *job), jobs))

    aligned: set[int] = set()
    unfinished: set[int] = set()
    for shard_aligned, shard_unfinished in results:
        aligned |= shard_aligned
        unfinished |= shard_unfinished

    if from_protein_id is None and original_session is None:
        victims = unfinished
    else:
        victims = set(candidates) - (aligned - unfinished)
    return sorted(victims), len(candidates)


def describe_unfinished_proteins(config, prot_ids: list[int], tables: set[str]) -> dict:
    """What the repair report shows: example identifiers, the sessions the proteins belong to,
    and how many structures aligned to them have no residue data."""
    info = {'examples': [], 'sessions': {}, 'structures': 0, 'structures_without_residues': 0}
    if len(prot_ids) == 0:
        return info

    db, cursor = _connect(config)
    try:
        info['examples'] = sorted(
            (row[0], row[1]) for row in _select_in(
                cursor, config, 'Protein_Id,Primary_Protein_Id', 'Protein', 'Protein_Id', prot_ids[:10]))

        proteins_per_session: dict[int, int] = {}
        for session, _prot_id in set(_select_in(
                cursor, config, 'Session,Protein', 'RS_Protein_Session', 'Protein', prot_ids)):
            proteins_per_session[session] = proteins_per_session.get(session, 0) + 1
        input_files = {row[0]: row[1] for row in _select_in(
            cursor, config, 'Session_Id,Input_File', 'Session', 'Session_Id', proteins_per_session.keys())}
        info['sessions'] = {
            session: (input_files.get(session), n) for session, n in proteins_per_session.items()}

        # Complex and Structure rows are written during alignment, their residues only during
        # structural analysis. A crash in between leaves structures that a rerun takes as stored
        # and never analyses. They are shared with the alignments of other proteins, so the
        # repair leaves them alone and only counts them. Structure ids are 100 * Structure_Id +
        # shard here, as everywhere else.
        structure_ids: set[int] = set()
        for shard, shard_prot_ids in _by_shard(prot_ids).items():
            structure_ids.update(row[0] for row in _select_in(
                cursor, config, 'DISTINCT Structure', f'Alignment_{shard}', 'Protein', shard_prot_ids))
        info['structures'] = len(structure_ids)

        per_structure_shard: dict[int, list[int]] = {}
        for internal_id in structure_ids:
            per_structure_shard.setdefault(internal_id % 100, []).append(internal_id // 100)
        for shard, ids in per_structure_shard.items():
            table = f'Structure_{shard}'
            if table not in tables:
                continue
            info['structures_without_residues'] += len(_select_in(
                cursor, config, 'Structure_Id', table, 'Structure_Id', ids, condition=' AND Residue_Data IS NULL'))
    finally:
        db.close()
    return info


def delete_proteins(config, prot_ids: list[int], tables: set[str]) -> dict[str, int]:
    """Delete these proteins and every row that hangs off them, in one transaction.

    One connection and one commit: if anything fails, the rollback leaves the batch exactly as it
    was, so the next repair finds it again and nothing is ever left half deleted. Returns the
    rows deleted per table, the shards of a table counted together.

    The session tables are left alone: the repair truncates them as a whole once every batch is
    through, which is far cheaper than deleting their rows protein by protein.
    """
    by_shard = _by_shard(prot_ids)
    deleted: dict[str, int] = {}
    db, cursor = _connect(config)
    try:
        def ids_of(column: str, table: str, key: str, values) -> set:
            if table not in tables or len(values) == 0:
                return set()
            return {row[0] for row in _select_in(cursor, config, column, table, key, values)}

        # Positions are referenced as 100 * Position_Id + shard, as everywhere else.
        position_ids = []
        for shard, shard_prot_ids in by_shard.items():
            position_ids.extend(100 * row[0] + shard for row in _select_in(
                cursor, config, 'Position_Id', f'Position_{shard}', 'Protein', shard_prot_ids))

        interface_ids = ids_of('Interface_Id', 'Interface', 'Protein', prot_ids)
        snv_ids = ids_of('SNV_Id', 'SNV', 'Position', position_ids)
        multi_mutation_ids = (ids_of('Multi_Mutation_Id', 'Multi_Mutation', 'Wildtype_Protein', prot_ids)
                              | ids_of('Multi_Mutation_Id', 'Multi_Mutation', 'Mutant_Protein', prot_ids))
        indel_ids = (ids_of('Indel_Id', 'Indel', 'Wildtype_Protein', prot_ids)
                     | ids_of('Indel_Id', 'Indel', 'Mutant_Protein', prot_ids))

        deletes = [(f'Protein_Protein_Interaction_{shard}', 'Protein_A', shard_prot_ids)
                   for shard, shard_prot_ids in by_shard.items()]
        # A PPI row is written from one side only, so a finished protein's row can name a
        # protein removed here as Protein_B - an id that would no longer exist. The removed
        # protein writes the interaction again from its own side when it is recomputed.
        deletes += [(f'Protein_Protein_Interaction_{shard}', 'Protein_B', prot_ids)
                    for shard in range(REPAIR_SHARDS)]
        deletes += [
            ('RS_Position_Interface', 'Interface', interface_ids),
            ('Interface', 'Interface_Id', interface_ids),
            ('SNV', 'SNV_Id', snv_ids),
            ('RS_Isoform', 'Multi_Mutation', multi_mutation_ids),
            ('RS_Isoform', 'Protein_A', prot_ids),
            ('RS_Isoform', 'Protein_B', prot_ids),
            ('Multi_Mutation', 'Multi_Mutation_Id', multi_mutation_ids),
            ('Indel', 'Indel_Id', indel_ids),
            ('RS_Protein_GO_Term', 'Protein', prot_ids),
            ('RS_Protein_Pathway', 'Protein', prot_ids),
        ]
        for shard, shard_prot_ids in by_shard.items():
            deletes.append((f'Alignment_{shard}', 'Protein', shard_prot_ids))
            deletes.append((f'Position_{shard}', 'Protein', shard_prot_ids))
        deletes.append(('Protein', 'Protein_Id', prot_ids))

        for table, key, values in deletes:
            if table not in tables or len(values) == 0:
                continue
            label = _table_label(table)
            deleted[label] = deleted.get(label, 0) + _delete_in(cursor, config, table, key, values)

        db.commit()
    except BaseException:
        db.rollback()
        raise
    finally:
        db.close()
    return deleted


def repair_crashed_chunks(config, session_id=None, from_protein_id=None, apply: bool = False,
                          batch_size: int = REPAIR_BATCH_SIZE, original_session=None) -> dict:
    """structman database repair: find the proteins a crashed run left unfinished, report them,
    and with apply remove them, so that rerunning the input computes them again.

    With apply it then truncates Session and the RS_*_Session tables (session_tables, as
    'database remove_sessions' does), so that every session is gone and a rerun links its proteins
    exactly once, instead of a second time on top of the links the crashed run left. That comes
    last, after every batch: an interrupted repair still finds its session to be resumed with.

    session_id scopes the search to the proteins linked to a session (-i); original_session to
    the proteins a session stored first (--session_id), which still works once the session
    tables are gone.

    Genes and all structure tables are kept. Returns what was found and deleted.
    """
    t0 = time.time()
    tables = _existing_tables(config)

    if original_session is not None:
        scope = f'the proteins session {original_session} stored first'
    elif session_id is not None:
        scope = f'session {session_id}'
    else:
        scope = 'the whole database'
    if from_protein_id is not None:
        scope = f'{scope}, Protein_Id >= {from_protein_id}'
    print(f'Looking for unfinished proteins in {scope} of {config.db_name} ...', flush=True)

    victims, n_checked = find_unfinished_proteins(config, session_id, from_protein_id, original_session)
    print(f'{n_checked} proteins checked in {time.time() - t0:.1f}s, {len(victims)} of them unfinished', flush=True)
    if original_session is not None and n_checked == 0:
        print_sessions_of_proteins(config)

    summary = {'victims': victims, 'deleted': {}, 'sessions_removed': 0}
    info = describe_unfinished_proteins(config, victims, tables)
    summary.update(info)

    if len(victims) > 0:
        print(f'  Protein_Id {victims[0]} .. {victims[-1]}, for example '
              + ', '.join(f'{primary_id} ({prot_id})' for prot_id, primary_id in info['examples']))
    for session, (input_file, n) in sorted(info['sessions'].items()):
        print(f'  {n} of them belong to session {session}: {input_file}')
    if info['structures_without_residues'] > 0:
        print(f'  Warning: {info["structures_without_residues"]} of the {info["structures"]} structures aligned to '
              'them have no residue data. If the run crashed before its structural analysis finished, a rerun '
              'takes these structures as stored and does not analyse them again. Structures whose analysis '
              'failed are counted here too.')

    all_sessions = _all_sessions(config) if 'Session' in tables else []

    if not apply:
        print(f'\nDry run, nothing was changed. --apply removes these proteins and then truncates '
              f'{", ".join(session_tables)}: all {len(all_sessions)} session(s) go, and every input has to be rerun '
              'for its output. Use it while no StructMAn run uses this database.')
        return summary

    t_apply = time.time()
    n_batches = (len(victims) + batch_size - 1) // batch_size
    report_every = max(1, n_batches // 20)
    deleted: dict[str, int] = {}
    for batch_nr, start in enumerate(range(0, len(victims), batch_size), start=1):
        batch = victims[start:start + batch_size]
        try:
            batch_deleted = delete_proteins(config, batch, tables)
        except BaseException:
            print(f'Batch {batch_nr}/{n_batches} failed and was rolled back. The {start} proteins removed before it '
                  'stay removed; run the repair again to continue.', flush=True)
            raise
        for label, n in batch_deleted.items():
            deleted[label] = deleted.get(label, 0) + n
        if batch_nr % report_every == 0 or batch_nr == n_batches:
            print(f'  {start + len(batch)}/{len(victims)} proteins removed ({time.time() - t_apply:.0f}s)', flush=True)

    summary['deleted'] = deleted
    if len(deleted) > 0:
        print('Rows removed: ' + ', '.join(f'{label} {n}' for label, n in sorted(deleted.items()) if n > 0))

    remove_sessions(config)
    summary['sessions_removed'] = len(all_sessions)
    print(f'Truncated {", ".join(session_tables)}: {len(all_sessions)} session(s) removed.')

    print('Rerun the inputs you need output for, with the command they were run with. Finished proteins are '
          'found stored and skipped, removed ones are computed again' + (':' if len(all_sessions) > 0 else '.'))
    for _session, input_file in all_sessions[:20]:
        print(f'  {input_file}')
    if len(all_sessions) > 20:
        print(f'  ... and {len(all_sessions) - 20} more')
    return summary


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

    #Circular import, no direct calling!
    #config = Config(config_path, None, None, None, None, False)

    main(config, reclassify_null=reclassify_null)
