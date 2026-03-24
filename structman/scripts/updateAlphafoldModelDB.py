import os
import sys
import shutil
import subprocess
import tarfile
import ray
import time
#from compression import zstd

from structman.base_utils.base_utils import calc_checksum, get_af_model_folder, topfolder_id_to_chunk_id
from structman.base_utils.config_class import Config
from structman.lib.database.database_core_functions import insert

compression_types = ['gzip' , 'zstd']
file_types = ['pdb', 'centralities', 'rin', 'msa', 'residue_list', 'interaction', 'network', 'unknown']
column_headers = ['ID', 'Structure_data', 'Sd_compr_type', 'Centralities_data', 'Cd_compr_type', 'RIN_data', 'Rd_compr_type', 'MSA_data', 'Md_compr_type', 'Version']

def disect_filename(fn: str) -> tuple[str, int, int, int]:
    parts = fn.split('.')
    if parts[-1] == 'zst':
        compr = 1
        ft_pos = -2
    elif parts[-1] == 'gz':
        compr = 0
        ft_pos = -2
    else:
        compr = None
        ft_pos = -1

    filetype_ext = parts[ft_pos]
    name_precursor = parts[0]

    try:
        raw_af_id, fn_add = name_precursor.split('_v')
    except ValueError:
        raw_af_id = name_precursor
        fn_add = ''
        if raw_af_id[:3] != 'AF-':
            return None, None, None, None

    if len(fn_add) > 0:
        version_nr = int(fn_add[0])
        if raw_af_id[-6:] == '-model':
            af_id = raw_af_id[3:-6]
        else:
            af_id = raw_af_id[3:-4]
        if len(fn_add) > 2:
            fn_add = fn_add[2:]
    else:
        version_nr = 0
        af_id = raw_af_id[3:-4]

    ft = 7

    match filetype_ext:
        case 'pdb':
            ft = 0
        
        case 'txt':
            if fn_add == 'btw_cent':
                ft = 1
            else:
                ft = 4
        
        case 'ea':
            if fn_add == 'intsc':
                ft = 2
            else:
                ft = 5

        case 'a3m':
            ft = 3
        
        case 'fasta':
            ft = 3

        case 'sif':
            ft = 6

    return af_id, ft, compr, version_nr


def init_database_chunks(config: Config, chunk_id: str):
    init_database_chunk(config, f'main_{chunk_id}')
    init_database_chunk(config, f'collab_{chunk_id}')

def init_database_chunk(config: Config, chunk_id: str):
    sql = f"""CREATE TABLE `AFDB_{chunk_id}` (
        `ID` varchar(64) NOT NULL,
        `Structure_data` longblob DEFAULT NULL,
        `Sd_compr_type` tinyint(4) DEFAULT NULL,
        `Centralities_data` longblob DEFAULT NULL,
        `Cd_compr_type` tinyint(4) DEFAULT NULL,
        `RIN_data` longblob DEFAULT NULL,
        `Rd_compr_type` tinyint(4) DEFAULT NULL,
        `MSA_data` longblob DEFAULT NULL,
        `Md_compr_type` tinyint(4) DEFAULT NULL,
        `Version` tinyint(4) NOT NULL
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8mb3_unicode_ci;
        """
    db, cursor = config.getDB(db_name = config.afdb)
    cursor.execute(sql)
    db.commit()

    sql = f"""ALTER TABLE `AFDB_{chunk_id}`
        ADD PRIMARY KEY (`ID`);
        """
    
    cursor.execute(sql)
    db.commit()
    db.close()

@ray.remote
def process_chunk(config, chunk_id, topfolder_list, lock):
    t0 = time.time()
    for tnr, topfolder in enumerate(topfolder_list):
        subfolder_list = os.listdir(f'{config.path_to_model_db}/{topfolder}')
        main_values = []
        collab_values = []
        pack = 0
        for subfolder in subfolder_list:
            #print(f'=====================\nSubfolder {subfolder}\n=====================')
            if len(os.listdir(f'{config.path_to_model_db}/{topfolder}/{subfolder}')) == 0:
                shutil.rmtree(f'{config.path_to_model_db}/{topfolder}/{subfolder}')
                continue

            id_batch_dict = {}
            for fn in os.listdir(f'{config.path_to_model_db}/{topfolder}/{subfolder}'):
                afdb_id, file_type, compr_type, version_nr = disect_filename(fn)

                if afdb_id is None:
                    print(f'Found faulty file: {fn=}, trying to delete it...')
                    os.remove(f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}')
                    continue
        
                if afdb_id not in id_batch_dict:
                    id_batch_dict[afdb_id] = {version_nr: [(file_type, compr_type, fn)]}
                elif version_nr not in id_batch_dict[afdb_id]:
                    id_batch_dict[afdb_id][version_nr] = [(file_type, compr_type, fn)]
                else:
                    id_batch_dict[afdb_id][version_nr].append((file_type, compr_type, fn))

            
            for afdb_id in id_batch_dict:
                
                newest_v = max(list(id_batch_dict[afdb_id].keys()))
                value_precursor = [afdb_id, None, None, None, None, None, None, None, None, newest_v]
                for file_triple in id_batch_dict[afdb_id][newest_v]:
                    ft, compr, fn = file_triple
                    match ft:
                        case 0:
                            f = open(f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}', 'rb')
                            value_precursor[1] = f.read()
                            f.close()
                            value_precursor[2] = compr

                        case 1:
                            f = open(f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}', 'rb')
                            value_precursor[3] = f.read()
                            f.close()
                            value_precursor[4] = compr

                        case 2:
                            f = open(f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}', 'rb')
                            value_precursor[5] = f.read()
                            f.close()
                            value_precursor[6] = compr

                        case 3:
                            if compr is not None:
                                f = open(f'{config.path_to_model_db}/{topfolder}/{subfolder}/{fn}', 'rb')
                                value_precursor[7] = f.read()
                                f.close()
                                value_precursor[8] = compr

                if afdb_id.count('-') > 0:
                    main_values.append(value_precursor)
                else:
                    collab_values.append(value_precursor)

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

def clone_to_sql(config: Config):
    print(f'Cloning AFDB to mySQL: {config.path_to_model_db=}')

    topfolder_list: list[str] = os.listdir(config.path_to_model_db)
    process_chunks = {}

    for tnr, topfolder in enumerate(topfolder_list):
        print(f'\n=====================\nTopfolder {topfolder} ({tnr+1}/{len(topfolder_list)})\n=====================\n')
        if len(topfolder) > 2:
            continue

        chunk_id = topfolder_id_to_chunk_id(topfolder)
        if chunk_id not in process_chunks:
            init_database_chunks(config, chunk_id)
            process_chunks[chunk_id] = []

        process_chunks[chunk_id].append(topfolder)

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
    clone_to_sql(config)
    sys.exit()
    if versioncheck(config):
        print('Local Alphafold model DB is already up-to-date.')
        return
    retrieve_raw_data(config)
    expand_raw_data(config)


    p = subprocess.Popen(['chmod', '777', '-R', '.'], cwd=config.path_to_model_db)
    p.wait()