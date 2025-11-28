import os
import subprocess
import tarfile

from structman.base_utils.base_utils import calc_checksum

def pull_meta_file(config):
    print('Call of pull_meta_file')
    cmds = ' '.join(['wget', 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.path_to_model_db)
    p.wait()

def versioncheck(config):
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


def retrieve_raw_data(config):
    cmds = ' '.join(['wget', 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP*'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.path_to_model_db)
    p.wait()
    cmds = ' '.join(['wget', 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest/swissprot_pdb_v2.tar'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.path_to_model_db)
    p.wait()

    #https://alphafold.ebi.ac.uk/files/msa/

def alphafold_model_id_to_folder_path(model_id, config, create_mode = False, get_version=False):
    uniprot_ac = model_id.split('-')[1]
    version_number = int(model_id.split('_')[1].split('.')[0][1:])
    topfolder_id = uniprot_ac[-2:]
    subfolder_id = uniprot_ac[-4:]

    full_path = f'{config.path_to_model_db}/{topfolder_id}/{subfolder_id}'

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

def main(config):
    print('Updating local Alphafold model DB')
    if versioncheck(config):
        print('Local Alphafold model DB is already up-to-date.')
        return
    retrieve_raw_data(config)
    expand_raw_data(config)
