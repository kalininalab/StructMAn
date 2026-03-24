import os
import gzip
import subprocess
import shutil

from structman.base_utils.base_utils import calc_checksum

def retrieve_raw_data(config):
    if not os.path.isdir(config.complex_model_db_path):
        os.makedirs(config.complex_model_db_path)
    checksum_file_path = f'{config.complex_model_db_path}/tarball_checksum.txt'

    cmds = ' '.join(['wget', 'https://conglab.swmed.edu/humanPPI/downloads/best_models.tar.gz', '--no-check-certificate'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.complex_model_db_path)
    p.wait()
    tar_path = f'{config.complex_model_db_path}/best_models.tar.gz'
    current_file_checksum = calc_checksum(tar_path)

    if os.path.isfile(checksum_file_path):
        f = open(checksum_file_path, 'r')
        try:
            old_checksum = int(f.read())
        except ValueError:
            old_checksum = None
        f.close()

        if old_checksum == current_file_checksum:
            return False
        
    f = open(checksum_file_path, 'w')
    f.write(str(current_file_checksum))
    f.close()
    return True

def complex_model_id_to_folder_path(config):
    cmds = ' '.join(['tar', '-xvzf', 'best_models.tar.gz'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.complex_model_db_path)
    p.wait()

    cmds = ' '.join(['rm', 'best_models.tar.gz'])
    p = subprocess.Popen(cmds, shell=True, cwd = config.complex_model_db_path)
    p.wait()

    for foldername in os.listdir(os.path.join(config.complex_model_db_path, 'best_models')):
        if not os.path.isdir(f'{config.complex_model_db_path}/best_models/{foldername}'):
            continue
        splitted = foldername.split("_")
        first_entry = splitted[0]
        first_folder = first_entry[-2:]
        second_folder = first_entry[-4:-2]

        target_dir = f'{config.complex_model_db_path}/{first_folder}/{second_folder}'
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        for file in os.listdir(f'{config.complex_model_db_path}/best_models/{foldername}'):
            if file.endswith('.pdb'):
                src = f'{config.complex_model_db_path}/best_models/{foldername}/{file}'
                dst = f'{target_dir}/{file}.gz'
                with open(src, 'rb') as f_in, gzip.open(dst, 'wb') as f_out:
                    f_out.writelines(f_in)
                
        shutil.rmtree(f'{config.complex_model_db_path}/best_models/{foldername}')
    shutil.rmtree(f'{config.complex_model_db_path}/best_models')

def main(config):
    print('Updating local Complex model DB')
    database_got_updated = retrieve_raw_data(config)
    if database_got_updated:
        complex_model_id_to_folder_path(config)
    else:
        os.remove(f'{config.complex_model_db_path}/best_models.tar.gz')

    p = subprocess.Popen(['chmod', '777', '-R', '.'], cwd=config.complex_model_db_path)
    p.wait()