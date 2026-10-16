import os
import sys
import csv
import shutil
import subprocess
import tarfile
import ray
import time
import gzip
import structman.base_utils.mysql_driver as pymysql
import numpy as np
import random
from typing import Callable
#from compression import zstd
from ftplib import FTP


from structman.base_utils.base_utils import (
    reset_logger_for_remotes, calc_checksum, get_af_model_folder, 
    topfolder_id_to_chunk_id, decode_zstd, expand_af_model_id, fetch_af_msa,
    disect_filename
)
from structman.base_utils.config_class import Config
from structman.lib.database.database_core_functions import insert
from structman.lib.database.afdb_retrieval import afdb_create_mmseqs_index
from structman.lib.rin_utils import calcRIN, get_entry_path, call_in_out_flow

compression_types = ['gzip' , 'zstd']
file_types = [
    'pdb',
    'centralities',
    'rin',
    'msa',
    'residue_list',
    'interaction',
    'network',
    'in_out_flow',
    'unknown']
column_headers = ['ID', 'Structure_data', 'Sd_compr_type', 'Centralities_data', 'Cd_compr_type', 'RIN_data', 'Rd_compr_type', 'In_out_flow_data', 'Iod_compr_type', 'MSA_data', 'Md_compr_type', 'Version']


def check_db(config: Config, just_count=False, fetch_msas=False):
    process_chunks = chunkify_file_db(config)
    processes = []
    for chunk_id in process_chunks:
        topfolder_list = process_chunks[chunk_id]
        processes.append(check_chunk.remote(config, topfolder_list, just_count=just_count, fetch_msas=fetch_msas))

    results = ray.get(processes)

    if just_count:
        complete = 0
        only_struct = 0
        struct_and_rin = 0
        struct_and_msa = 0
        for (sub_complete, sub_only_struct, sub_struct_and_rin, sub_struct_and_msa) in results:
            complete += sub_complete
            only_struct += sub_only_struct
            struct_and_rin += sub_struct_and_rin
            struct_and_msa += sub_struct_and_msa

        print(f'{complete=} {only_struct=} {struct_and_rin=} {struct_and_msa=}')


def init_database_chunks(config: Config, chunk_id: str):
    init_database_chunk(config, f'main_{chunk_id}', config.afdb)
    init_database_chunk(config, f'collab_{chunk_id}', config.afdb)

def init_database_chunk(config: Config, chunk_id: str, dbname: str):

    sql = f"""CREATE TABLE `AFDB_{chunk_id}` (
        `ID` varchar(64) NOT NULL,
        `Structure_data` longblob DEFAULT NULL,
        `Sd_compr_type` tinyint(4) DEFAULT NULL,
        `Centralities_data` longblob DEFAULT NULL,
        `Cd_compr_type` tinyint(4) DEFAULT NULL,
        `RIN_data` longblob DEFAULT NULL,
        `Rd_compr_type` tinyint(4) DEFAULT NULL,
        `In_out_flow_data` longblob DEFAULT NULL,
        `Iod_compr_type` tinyint(4) DEFAULT NULL,
        `MSA_data` longblob DEFAULT NULL,
        `Md_compr_type` tinyint(4) DEFAULT NULL,
        `Version` tinyint(4) NOT NULL
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;
        """
    db, cursor = config.getDB(db_name = dbname)
    try:
        cursor.execute(sql)
        db.commit()

        sql = f"""ALTER TABLE `AFDB_{chunk_id}`
            ADD PRIMARY KEY (`ID`);
            """
        
        cursor.execute(sql)
        db.commit()
    except pymysql.err.OperationalError:
        sql =  f'TRUNCATE AFDB_{chunk_id};'
        cursor.execute(sql)
        db.commit()
    finally:
        db.close()

def call_af_model_rin_calc(config: Config, struct_compr: int, struct_file: str, afdb_id: str, custom_db_path=None):
    if struct_compr == 0:
        f = gzip.open(struct_file, 'rb')
        page = f.read()
        f.close()
    elif struct_compr == 1:
        f = open(struct_file, 'rb')
        page = f.read()
        f.close()
        page = decode_zstd(page)
    elif struct_compr == -1:
        inflow_file = struct_file
        outfile = inflow_file.replace('_intsc.ea.gz','_in_out_flow.csv')
        call_in_out_flow(inflow_file, outfile)
        return
    else:
        try:
            f = open(struct_file, 'rb')
        except FileNotFoundError:
            return
        page = f.read()
        f.close()

    if custom_db_path is not None:
        path_to_model_db = None
        custom_db_path = custom_db_path
    else:
        path_to_model_db = config.path_to_model_db
        custom_db_path = None

    folder_path, path_stem, structure_id = get_entry_path(expand_af_model_id(afdb_id), config.rin_db_path, path_to_model_db=path_to_model_db, custom_db_path=custom_db_path, temp_folder=config.tmp_folder)
    try:
        calcRIN(config, page, folder_path, path_stem, structure_id, config.rinerator_path, True)
    except UnboundLocalError as err:
        raise UnboundLocalError(f'{afdb_id=} {struct_file=} {struct_compr=} {err=}')
    except FileNotFoundError as err:
        raise FileNotFoundError(f'{afdb_id=} {struct_file=} {struct_compr=} {err=}')

def check_subfolder(config: Config, subfolder: str, topfolder: str, fetch_msas: bool = False, inqueue = None, just_count=False, from_scratch=False):
    N = 0
    id_batch_dict = {}
    for fn in os.listdir(f'{config.path_to_model_db}/{topfolder}/{subfolder}'):
        afdb_id, file_type, compr_type, version_nr = disect_filename(fn)

        if afdb_id is None:
            if not just_count:
                print(f'Found faulty file: {fn=}, trying to delete it...')
                if os.path.exists(f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}'):
                    os.remove(f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}')
            continue

        if afdb_id not in id_batch_dict:
            id_batch_dict[afdb_id] = {version_nr: [(file_type, compr_type, fn)]}
        elif version_nr not in id_batch_dict[afdb_id]:
            id_batch_dict[afdb_id][version_nr] = [(file_type, compr_type, fn)]
        else:
            id_batch_dict[afdb_id][version_nr].append((file_type, compr_type, fn))
    
    complete = 0
    only_struct = 0
    struct_and_rin = 0
    struct_and_msa = 0
    print_first = False

    for afdb_id in id_batch_dict:
        
        if -1 in id_batch_dict[afdb_id]:
            for ver in id_batch_dict[afdb_id]:
                if ver == -1:
                    continue
                id_batch_dict[afdb_id][ver] += id_batch_dict[afdb_id][-1]

        #newest_v = max(list(id_batch_dict[afdb_id].keys()))
        versions_in_order = sorted(id_batch_dict[afdb_id].keys(), reverse=True)
        no_rin = True
        no_cent = True
        no_iof = True
        no_msa = True
        struct_file = None
        rin_file = None
        struct_compr = None
        for version_number in versions_in_order:
            for file_triple in id_batch_dict[afdb_id][version_number]:
                ft, compr, fn = file_triple
                match ft:
                    case 0:
                        if struct_file is None:
                            struct_file = f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}'
                            struct_compr = compr
                    case 1:
                        no_cent = False
                    case 2:
                        no_rin = False
                        rin_file = f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}'
                    case 3:
                        no_msa = False
                    case 7:
                        no_iof = False
        if struct_file is None:
            continue

        if no_rin or no_cent or from_scratch:
            if just_count:
                if no_msa:
                    only_struct += 1
                    if not print_first:
                        print(f'Entry with only struct: {afdb_id=} {id_batch_dict[afdb_id]=} {struct_file=} {struct_compr=}')
                        print_first = True
                else:
                    struct_and_msa += 1
            else:
                if inqueue is not None:
                    inqueue.put((struct_compr, struct_file, afdb_id))
                    N += 1
                else:
                    call_af_model_rin_calc(config, struct_compr, struct_file, afdb_id)
        elif just_count:
            if no_msa:
                struct_and_rin += 1
            else:
                complete += 1

        elif no_iof and not no_rin:
            
            if inqueue is not None:
                inqueue.put((-1, rin_file, afdb_id))
                N += 1
            else:
                call_af_model_rin_calc(config, -1, rin_file, afdb_id)

        if no_msa and fetch_msas and afdb_id.count('-') > 0 and not just_count:
            fetch_af_msa(config, afdb_id, f'{config.path_to_model_db}/{topfolder}/{subfolder}', encode_it=True)
    if just_count:
        return (complete, only_struct, struct_and_rin, struct_and_msa)
    return N

@ray.remote
def check_chunk(config, topfolder_list, just_count=False, fetch_msas=False):
    reset_logger_for_remotes(config)
    complete = 0
    only_struct = 0
    struct_and_rin = 0
    struct_and_msa = 0
    for topfolder in topfolder_list:
        subfolder_list = os.listdir(f'{config.path_to_model_db}/{topfolder}')
        
        for subfolder in subfolder_list:
            returns = check_subfolder(config, subfolder, topfolder, fetch_msas= fetch_msas, just_count=just_count)
            if just_count:
                (sub_complete, sub_only_struct, sub_struct_and_rin, sub_struct_and_msa) = returns
                complete += sub_complete
                only_struct += sub_only_struct
                struct_and_rin += sub_struct_and_rin
                struct_and_msa += sub_struct_and_msa

    return (complete, only_struct, struct_and_rin, struct_and_msa)

def extract_values_from_subfolder(
        config: Config,
        topfolder: str,
        subfolder: str,
        main_values: list[tuple],
        contig_shift: int = 0,
        high_confs: np.ndarray[bool] = None,
        collab_values: list[tuple] = None,
        detection_delimiter: str = None
        ):
    #print(f'=====================\nSubfolder {subfolder}\n=====================')
    if len(os.listdir(f'{topfolder}/{subfolder}')) == 0:
        shutil.rmtree(f'{topfolder}/{subfolder}')
        return

    id_batch_dict = {}
    for fn in os.listdir(f'{topfolder}/{subfolder}'):
        afdb_id, file_type, compr_type, version_nr = disect_filename(fn)

        if afdb_id is None:
            print(f'Found faulty file: {fn=}, trying to delete it...')
            try:
                os.remove(f'{topfolder}/{subfolder}/{fn}')
            except FileNotFoundError:
                pass
            continue

        if high_confs is not None:
            id_contig = int(afdb_id)
            shifted_contig = (id_contig-contig_shift)-1
            is_high_conf_entry = high_confs[shifted_contig]

            if not is_high_conf_entry:
                continue

        if afdb_id not in id_batch_dict:
            id_batch_dict[afdb_id] = {version_nr: [(file_type, compr_type, fn)]}
        elif version_nr not in id_batch_dict[afdb_id]:
            id_batch_dict[afdb_id][version_nr] = [(file_type, compr_type, fn)]
        else:
            id_batch_dict[afdb_id][version_nr].append((file_type, compr_type, fn))

    
    for afdb_id in id_batch_dict:
        if -1 in id_batch_dict[afdb_id]:
            for ver in id_batch_dict[afdb_id]:
                if ver == -1:
                    continue
                id_batch_dict[afdb_id][ver] += id_batch_dict[afdb_id][-1]

        #print(f'{afdb_id=} {id_batch_dict[afdb_id]=}')
        newest_v = max(list(id_batch_dict[afdb_id].keys()))
        value_precursor = [afdb_id, None, None, None, None, None, None, None, None, None, None, newest_v]
        for file_triple in id_batch_dict[afdb_id][newest_v]:
            ft, compr, fn = file_triple
            match ft:
                case 0:
                    f = open(f'{topfolder}/{subfolder}/{fn}', 'rb')
                    value_precursor[1] = f.read()
                    f.close()
                    value_precursor[2] = compr

                case 1:
                    f = open(f'{topfolder}/{subfolder}/{fn}', 'rb')
                    value_precursor[3] = f.read()
                    f.close()
                    value_precursor[4] = compr

                case 2:
                    f = open(f'{topfolder}/{subfolder}/{fn}', 'rb')
                    value_precursor[5] = f.read()
                    f.close()
                    value_precursor[6] = compr

                case 7:
                    f = open(f'{topfolder}/{subfolder}/{fn}', 'rb')
                    value_precursor[7] = f.read()
                    f.close()
                    value_precursor[8] = compr

                case 3:
                    if compr is not None:
                        f = open(f'{topfolder}/{subfolder}/{fn}', 'rb')
                        value_precursor[9] = f.read()
                        f.close()
                        value_precursor[10] = compr

        #print(f'{afdb_id=} {value_precursor=}')

        if collab_values is not None:
            if afdb_id.count(detection_delimiter) > 0:
                main_values.append(value_precursor)
            else:
                collab_values.append(value_precursor)
        else:
            main_values.append(value_precursor)


@ray.remote
def process_chunk(config, chunk_id, topfolder_list, lock):
    reset_logger_for_remotes(config)
    t0 = time.time()
    for tnr, topfolder in enumerate(topfolder_list):
        subfolder_list = os.listdir(f'{config.path_to_model_db}/{topfolder}')
        main_values = []
        collab_values = []
        pack = 0
        for subfolder in subfolder_list:
            extract_values_from_subfolder(config, f'{config.path_to_model_db}/{topfolder}', subfolder, main_values, collab_values=collab_values, detection_delimiter='-')
            
            if pack == 10:
                if len(main_values) > 0:
                    insert(f'AFDB_main_{chunk_id}', column_headers, main_values, config, db_name=config.afdb, db_lock=lock)

                if len(collab_values) > 0:
                    insert(f'AFDB_collab_{chunk_id}', column_headers, collab_values, config, db_name=config.afdb, db_lock=lock)
                main_values = []
                collab_values = []
                pack = 0
            else:
                pack += 1

        if len(main_values) > 0:
            insert(f'AFDB_main_{chunk_id}', column_headers, main_values, config, db_name=config.afdb, db_lock=lock)

        if len(collab_values) > 0:
            insert(f'AFDB_collab_{chunk_id}', column_headers, collab_values, config, db_name=config.afdb, db_lock=lock)

        subfolder_list = os.listdir(f'{config.path_to_model_db}/{topfolder}')
        if len(subfolder_list) == 0:
            shutil.rmtree(f'{config.path_to_model_db}/{topfolder}')

    t1 = time.time()

    print(f'Time for {chunk_id=}: {t1-t0}')

@ray.remote
def clone_folder_to_sql(folder_list: list[str], config_dump: tuple[Config, int], dbname: str, high_confs: np.ndarray[bool]):
    config, contig_shift = config_dump

    reset_logger_for_remotes(config)
    
    for tnr, topfolder in enumerate(folder_list):
        t0 = time.time()
        subfolder_list = os.listdir(topfolder)
        main_values = []

        for subfolder in subfolder_list:
            if len(subfolder) != 2:
                continue
            if not os.path.isdir(f'{topfolder}/{subfolder}'):
                continue

            extract_values_from_subfolder(config, topfolder, subfolder, main_values, contig_shift=contig_shift, high_confs=high_confs)
            
            if len(main_values) > 0:
                insert(f'AFDB_main_{subfolder}', column_headers, main_values, config, db_name=dbname)

            main_values = []
                

        if len(main_values) > 0:
            insert(f'AFDB_main_{subfolder}', column_headers, main_values, config, db_name=dbname)

        t1 = time.time()

        print(f'Time for {topfolder=}: {t1-t0}')

def chunkify_file_db(config: Config, iter_func: Callable | None = None):
    topfolder_list: list[str] = os.listdir(config.path_to_model_db)
    process_chunks = {}

    for tnr, topfolder in enumerate(topfolder_list):
        #print(f'\n=====================\nTopfolder {topfolder} ({tnr+1}/{len(topfolder_list)})\n=====================\n')
        if len(topfolder) > 2:
            continue

        chunk_id = topfolder_id_to_chunk_id(topfolder)
        if chunk_id not in process_chunks:
            if iter_func is not None:
                iter_func(config, chunk_id)
            process_chunks[chunk_id] = []

        process_chunks[chunk_id].append(topfolder)
    return process_chunks

def clone_to_sql(config: Config):
    print(f'Cloning AFDB to mySQL: {config.path_to_model_db=}')

    process_chunks = chunkify_file_db(config, iter_func=init_database_chunks)

    nr_of_locks = 20
    current_lock = 0
    processes = []
    for chunk_id in process_chunks:
        topfolder_list = process_chunks[chunk_id]
        processes.append(process_chunk.remote(config, chunk_id, topfolder_list, str(current_lock)))
        current_lock += 1
        if current_lock == nr_of_locks:
            current_lock = 0

    ray.get(processes)
    return

def clone_complex_to_sql(
        config: Config,
        db_name: str,
        folder_prefix: str,
        high_conf_models: np.ndarray[bool],
        contig_shift: int
        ):

    db, cursor = config.getDB(db_name = db_name)
    sql = f'DROP DATABASE `{db_name}`'
    try:
        cursor.execute(sql)
    except pymysql.err.OperationalError:
        pass
    db.close()

    db, cursor = config.getDB(server_connection=True)
    sql = f'CREATE DATABASE `{db_name}`'

    cursor.execute(sql)
    db.commit()

    db.close()

    folder_list: list[str] = os.listdir(config.path_to_af_complex_db)

    chunks: list[list[str]] = []

    current_chunk = 0

    num_of_chunks = 50 #limited by sql server connections

    conf_dump = ray.put((config, contig_shift))

    chunks_inited = set()

    for fol in folder_list:
        if not fol.startswith(folder_prefix):
            continue

        fp = f'{config.path_to_af_complex_db}/{fol}'

        if not os.path.isdir(fp):
            continue

        if len(chunks) == current_chunk:
            chunks.append([])

        subfolder_list = os.listdir(fp)

        for subfolder in subfolder_list:
            if len(subfolder) != 2:
                continue
            if not os.path.isdir(f'{fp}/{subfolder}'):
                continue

            if subfolder not in chunks_inited:
                init_database_chunk(config, f'main_{subfolder}', db_name)
                chunks_inited.add(subfolder)

        chunks[current_chunk].append(fp)

        current_chunk += 1
        if current_chunk == num_of_chunks:
            current_chunk = 0

    stored_high_confs = ray.put(high_conf_models)

    process_ids = []
    for chunk in chunks:
        process_ids.append(clone_folder_to_sql.remote(chunk, conf_dump, db_name, stored_high_confs))

    ray.get(process_ids)

    return


def write_high_conf_models(high_conf_models: np.ndarray[bool], contig_shift: int, outfile: str):
    outlines = []
    for shifted_contig in range(len(high_conf_models)):
        id_contig = shifted_contig + 1 + contig_shift
        is_high_conf_entry = high_conf_models[shifted_contig]

        if is_high_conf_entry:
            outlines.append(f'{id_contig}\n')

    f = open(outfile, 'w')
    f.write(''.join(outlines))
    f.close()

def clone_complex_dbs_to_sql(config: Config):

    hetero_meta_file = f'{config.path_to_af_complex_db}/heterodimer_metadata.csv'
    hetero_high_conf_models, hetero_contig_shift = get_high_confidence_entries(hetero_meta_file)

    clone_complex_to_sql(config, config.afdb_hetero, 'heteromer_shard', hetero_high_conf_models, hetero_contig_shift)

    mmseqs_fasta_path = config.afdb_hetero_mmseqs_path

    afdb_create_mmseqs_index(config, config.afdb_hetero, mmseqs_fasta_path)

    homo_meta_file = f'{config.path_to_af_complex_db}/homodimer_metadata.csv'
    homo_high_conf_models, homo_contig_shift = get_high_confidence_entries(homo_meta_file)

    clone_complex_to_sql(config, config.afdb_homo, 'homomer_chunk', homo_high_conf_models, homo_contig_shift)

    afdb_create_mmseqs_index(config, config.afdb_homo, config.afdb_homo_mmseqs_path)

def pull_meta_file(config: Config):
    print('Call of pull_meta_file')
    cmds = ' '.join(['wget', 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.path_to_model_db)
    p.wait()


def versioncheck(config: Config):
    path_to_meta_file = f'{config.path_to_model_db}/download_metadata.json'
    if not os.path.exists(path_to_meta_file):
        pull_meta_file(config)
        return False
    current_meta_file_checksum = calc_checksum(path_to_meta_file)
    print(f'{current_meta_file_checksum=}')
    os.remove(path_to_meta_file)
    pull_meta_file(config)
    new_meta_file_checksum = calc_checksum(path_to_meta_file)
    print(f'{new_meta_file_checksum=}')

    return current_meta_file_checksum == new_meta_file_checksum


def retrieve_collab_data(config: Config):
    collabs = ['atbc', 'bfvd', 'ntdx']
    for collab in collabs:
        for subfolder in ['models', 'msas']:
            cmds = ' '.join(['wget', f'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/{collab}/{subfolder}/*'])
            p = subprocess.Popen(cmds, shell=True, cwd = config.path_to_model_db)
            p.wait()

@ray.remote
def process_complex_data_folder(folder_list: list[str], conf_dump):

    config, hetero_high_conf_models, hetero_contig_shift, homo_high_conf_models, homo_contig_shift = conf_dump

    for fp in folder_list:
        content = os.listdir(fp)

        fol = fp.rsplit('/',1)[1]

        if fol.startswith('heteromer_shard'):
            high_confs = hetero_high_conf_models
            contig_shift = hetero_contig_shift
        else:
            high_confs = homo_high_conf_models
            contig_shift = homo_contig_shift

        for fn in content:
            combined_path = f'{fp}/{fn}'
            if os.path.isdir(combined_path):

                continue

            if not fn.startswith('AF-'):
                continue

            afdb_id, file_type, compr_type, version_nr = disect_filename(fn)

            id_contig = int(afdb_id)
            shifted_contig = (id_contig-contig_shift)-1
            is_high_conf_entry = high_confs[shifted_contig]

            if not is_high_conf_entry:
                continue

            sf_code = afdb_id[-2:]

            subfolder = f'{fp}/{sf_code}'

            if not os.path.isdir(subfolder):
                os.makedirs(subfolder)

            target_path = f'{subfolder}/{fn}'

            shutil.move(combined_path, target_path)

            call_af_model_rin_calc(config, compr_type, target_path, afdb_id, custom_db_path=subfolder)


def process_complex_data(config: Config):
    folder_list: list[str] = os.listdir(config.path_to_af_complex_db)

    chunks: list[list[str]] = []

    current_chunk = 0

    hetero_meta_file = f'{config.path_to_af_complex_db}/heterodimer_metadata.csv'
    hetero_high_conf_models, hetero_contig_shift = get_high_confidence_entries(hetero_meta_file)

    write_high_conf_models(hetero_high_conf_models, hetero_contig_shift, f'{config.path_to_af_complex_db}/heterodimer_high_conf_ids.txt')

    homo_meta_file = f'{config.path_to_af_complex_db}/homodimer_metadata.csv'
    homo_high_conf_models, homo_contig_shift = get_high_confidence_entries(homo_meta_file)

    write_high_conf_models(homo_high_conf_models, homo_contig_shift, f'{config.path_to_af_complex_db}/homodimer_high_conf_ids.txt')

    sys.exit()

    conf_dump = ray.put((config, hetero_high_conf_models, hetero_contig_shift, homo_high_conf_models, homo_contig_shift))


    random.shuffle(folder_list)
    for fol in folder_list:
        fp = f'{config.path_to_af_complex_db}/{fol}'

        if not os.path.isdir(fp):
            continue

        if not (fol.startswith('heteromer_shard') or fol.startswith('homomer_chunk')):
            continue

        if len(chunks) == current_chunk:
            chunks.append([])

        chunks[current_chunk].append(fp)

        current_chunk += 1
        if current_chunk == config.proc_n:
            current_chunk = 0

    process_ids = []
    for chunk in chunks:
        process_ids.append(process_complex_data_folder.remote(chunk, conf_dump))

    ray.get(process_ids)

def af_id_to_int(af_id: str) -> int:
    number_part = af_id[3:]
    return int(number_part)

def get_high_confidence_entries(meta_file: str) -> tuple[np.ndarray[bool], int]:
    print(f'Parsing metadata file: {meta_file}')
    f = open(meta_file, newline='')
    rows = list(csv.reader(f))
    f.close()

    headers = rows[0]

    for pos, header in enumerate(headers):
        if header == 'modelEntityId':
            id_pos = pos

        if header == 'ipSAE_AB':
            ipsae_ab_pos = pos
        if header == 'ipSAE_BA':
            ipsae_ba_pos = pos

        if header == 'pDockQ':
            pdock_pos = pos

        if header == 'N_clash_backbone':
            clash_pos = pos

    print(f'{headers=}')

    contigs = []
    for words in rows[1:]:
        model_id = words[id_pos]
        model_contig = af_id_to_int(model_id)
        contigs.append(model_contig)

    min_contig = min(contigs)
    max_contig = max(contigs)

    print(f'Contig range: {min_contig} - {max_contig}')

    high_conf_models = np.zeros(max_contig-min_contig, dtype=bool) 
    high_conf_chunks = set()

    for words in rows[1:]:

        model_id = words[id_pos]
        ipsae_ab = float(words[ipsae_ab_pos])
        ipsae_ba = float(words[ipsae_ba_pos])
        pdockq = float(words[pdock_pos])

        chunk_name = words[-1]

        try:
            n_clash = int(words[clash_pos])
        except ValueError:
            try:
                n_clash = int(float(words[clash_pos]))
            except ValueError as err:
                raise ValueError(f'{err=} {words=}')

        ipsae_min = min([ipsae_ab, ipsae_ba])

        if pdockq >= 0.23 and ipsae_min >= 0.6 and n_clash <= 10:
            model_contig = af_id_to_int(model_id)
            high_conf_models[(model_contig-min_contig)-1] = True
            high_conf_chunks.add(chunk_name)

    print(f'{np.sum(high_conf_models)=} {len(high_conf_chunks)=}')

    return high_conf_models, min_contig

def retrieve_complex_data(config: Config):
    if not os.path.isdir(config.path_to_af_complex_db):
        os.makedirs(config.path_to_af_complex_db)

    heteromer_folder = 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/heterodimers/'

    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login("anonymous", '')

    # Change to target directory and grab clean filenames
    ftp.cwd("/pub/databases/alphafold/collaborations/nvda/heterodimers")
    filenames = ftp.nlst()
    ftp.quit()
    
    for name in filenames:
        
        if name.startswith('shard'):
            subfolder = f'{config.path_to_af_complex_db}/heteromer_{name[:-4]}'
            if not os.path.isdir(subfolder):
                os.makedirs(subfolder)

                cmds = ["curl", f"{heteromer_folder}/{name}"]

                p1 = subprocess.Popen(
                    cmds,
                    cwd = subfolder,
                    stdout=subprocess.PIPE,)
                
                cmds = ["tar", "-xvf", "-", "--wildcards", "*.pdb.zst"]

                p2 = subprocess.Popen(
                    cmds,
                    cwd = subfolder,
                    stdin=p1.stdout,)
                p1.stdout.close()
                p2.communicate()

    homomer_folder = 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/homodimers/'
    
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login("anonymous", '')

    # Change to target directory and grab clean filenames
    ftp.cwd("/pub/databases/alphafold/collaborations/nvda/homodimers")
    filenames = ftp.nlst()
    ftp.quit()

    for name in filenames:

        if name.startswith('chunk'):
            subfolder = f'{config.path_to_af_complex_db}/homomer_{name[:-4]}'
            if not os.path.isdir(subfolder):
                os.makedirs(subfolder)

                cmds = ["curl", f"{homomer_folder}/{name}"]

                p1 = subprocess.Popen(
                    cmds,
                    cwd = subfolder,
                    stdout=subprocess.PIPE,)
                
                cmds = ["tar", "-xvf", "-", "--wildcards", "*.pdb.zst"]

                p2 = subprocess.Popen(
                    cmds,
                    cwd = subfolder,
                    stdin=p1.stdout,)
                p1.stdout.close()
                p2.communicate()


def retrieve_raw_data(config: Config):
    cmds = ' '.join(['wget', 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP*'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.path_to_model_db)
    p.wait()
    cmds = ' '.join(['wget', 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest/swissprot_pdb_v2.tar'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.path_to_model_db)
    p.wait()

    retrieve_collab_data(config)

    #https://alphafold.ebi.ac.uk/files/msa/

def alphafold_model_id_to_folder_path(model_id, config, create_mode = False, get_version=False):

    full_path, topfolder_id = get_af_model_folder(model_id, config, get_topfolder_id=True)

    version_number = int(model_id.split('_')[1].split('.')[0][1:])

    if not create_mode:
        if get_version:
            return full_path, version_number
        return full_path
    
    if not os.path.isdir(f'{config.path_to_model_db}/{topfolder_id}'):
        os.mkdir(f'{config.path_to_model_db}/{topfolder_id}')

    if not os.path.isdir(full_path):
        os.mkdir(full_path)

    if get_version:
        return full_path, version_number
    return full_path

def purge_old_versions(folder, new_version_number, new_model_id):
    new_model_stem = new_model_id.split('_')[0]

    for filename in os.listdir(folder):
        if filename.count('_') < 1:
            continue
        stem, tail = filename.split('_')[:2]
        if stem == new_model_stem:
            version = int(tail.split('_')[0].split('.')[0][1:])
            if new_version_number > version:
                os.remove(f'{folder}/{filename}')

def expand_raw_collab_data(config):
    for filename in os.listdir(config.path_to_model_db):
        if filename[:5] != 'chunk':
            continue
        
        tarball = tarfile.open(f'{config.path_to_model_db}/{filename}')
        for member in tarball:
            if len(member.name) < 3:
                continue
            memberfilename = member.name.split('/')[-1]
            if memberfilename.count('.pdb') == 0 and memberfilename.count('.a3m') == 0:
                continue
            folder_path, version_number = alphafold_model_id_to_folder_path(memberfilename, config, create_mode = True, get_version=True)
            purge_old_versions(folder_path, version_number, memberfilename)
            
            if memberfilename[-7:] == 'a3m.zst':
                target = f'{folder_path}/{memberfilename}'

                with open(target, 'wb') as tf:
                    f = tarball.extractfile(member)
                    try:
                        tf.write(f.read())
                    except tarfile.ReadError as e:
                        print(f'ERROR: reading from {filename=} {member.name=} failed: {e=}')
                        break


            else:
                tarball.extract(member, path = f'{folder_path}')
        tarball.close()

        os.remove(f'{config.path_to_model_db}/{filename}')

def expand_raw_data(config):
    for filename in os.listdir(config.path_to_model_db):
        if filename[-4:] != '.tar':
            continue
        tarball = tarfile.open(f'{config.path_to_model_db}/{filename}')
        for member in tarball:
            if member.name.count('.pdb') > 0:
                folder_path, version_number = alphafold_model_id_to_folder_path(member.name, config, create_mode = True, get_version=True)
                purge_old_versions(folder_path, version_number, member.name)
                tarball.extract(member, path = f'{folder_path}')
        tarball.close()

        os.remove(f'{config.path_to_model_db}/{filename}')

    expand_raw_collab_data(config)

def main(config):
    print('Updating local Alphafold model DB')

    #retrieve_complex_data(config)
    #sys.exit()

    #check_db(config, just_count=False)
    #clone_to_sql(config)
    
    process_complex_data(config)
    clone_complex_dbs_to_sql(config)
    sys.exit()

    if versioncheck(config):
        print('Local Alphafold model DB is already up-to-date.')
        return
    retrieve_raw_data(config)
    expand_raw_data(config)

    p = subprocess.Popen(['chmod', '777', '-R', '.'], cwd=config.path_to_model_db)
    p.wait()

    clone_to_sql(config)