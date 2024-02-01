#!/usr/bin/python3

import os
import shutil
import sys

from structman import settings
from structman.base_utils.base_utils import calc_checksum

def update_file(source, target):
    if not os.path.exists(target):
        shutil.copy(source, target)
        return
    if calc_checksum(source) == calc_checksum(target):
        return
    shutil.copy(source, target)
    return

def update_folder(source_folder_path, target_folder_path, recursive = False, blacklist = None):
    if blacklist is None:
        blacklist = set([])
    blacklist.add('__pycache__')
    if not os.path.isdir(target_folder_path):
        os.mkdir(target_folder_path)
    for entry in os.listdir(source_folder_path):
        if entry in blacklist:
            continue
        source = f'{source_folder_path}/{entry}'
        target = f'{target_folder_path}/{entry}'
        if os.path.isfile(source):
            update_file(source, target)
        elif os.path.isdir(target_folder_path) and recursive:
            update_folder(source, target, recursive = True, blacklist = blacklist)

if __name__ == "__main__":

    update_script_folder = os.path.abspath(sys.argv[0].rsplit('/',1)[0])

    target_folder = sys.argv[1]

    structman_target_folder = f'{target_folder}/structman'
    lib_target_folder = f'{structman_target_folder}/lib'
    rinerator_target_folder = f'{lib_target_folder}/rinerator'

    update_file(f'{update_script_folder}/../../changelog.txt', f'{target_folder}/changelog.txt')
    update_file(f'{update_script_folder}/../../MANIFEST.in', f'{target_folder}/MANIFEST.in')
    update_file(f'{update_script_folder}/../../public_config_template.txt', f'{target_folder}/config_template.txt')
    update_file(f'{update_script_folder}/../../public_install.sh', f'{target_folder}/install.sh')
    update_file(f'{update_script_folder}/../../setup.py', f'{target_folder}/setup.py')

    update_folder(settings.ROOT_DIR, structman_target_folder)

    utils_lib_path = f'{settings.ROOT_DIR}/base_utils'
    utils_target_path = f'{target_folder}/structman/base_utils'
    update_folder(utils_lib_path, utils_target_path)

    scripts_lib_path = f'{settings.ROOT_DIR}/scripts'
    scripts_target_path = f'{target_folder}/structman/scripts'
    scripts_blacklist = set(['updateDockerSource.py', 'updatePublicVersion.py'])
    update_folder(scripts_lib_path, scripts_target_path)

    resources_target_path = f'{target_folder}/structman/resources'
    update_folder(settings.RESOURCES_DIR, resources_target_path)

    lib_target_path =  f'{target_folder}/structman/lib'
    update_folder(settings.LIB_DIR, lib_target_path, recursive = True)    


