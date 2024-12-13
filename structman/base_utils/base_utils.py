import math
import os
import time
import traceback

from zlib import adler32

import zstd
import msgpack
import subprocess

from structman.base_utils.custom_encoder import custom_encoder, custom_decoder, custom_encoder_complete

class Errorlog:
    def __init__(self, path=None, warn_path=None, print_all_errors=False, print_all_warns=False):
        self.path = path
        self.warn_path = warn_path
        self.error_counter = 0
        self.warning_counter = 0
        self.print_all_errors = print_all_errors
        self.print_all_warns = print_all_warns

    def start(self, nfname, session):
        if self.warn_path is not None:
            date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
            errortext = "###############################################################################\n%s:\n%s - Session: %s\n" % (date, nfname, str(session))
            f = open(self.warn_path, 'a')
            f.write(errortext)
            f.close()

        if self.path is None:
            return
        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
        errortext = "###############################################################################\n%s:\n%s - Session: %s\n" % (date, nfname, str(session))
        f = open(self.path, 'a')
        f.write(errortext)
        f.close()

    def add_error(self, error_text):
        self.error_counter += 1
        g = ''.join(traceback.format_list(traceback.extract_stack()))
        error_text = 'Error %s:\n%s\n%s\n' % (str(self.error_counter), error_text, g)
        if self.path is None or self.print_all_errors or self.print_all_warns:
            print(error_text)
            return
        f = open(self.path, 'a')
        f.write(error_text)
        f.close()

    def stop(self):
        if self.warn_path is not None:
            if self.warning_counter == 0:
                errortext = "Finished without any warnings\n###############################################################################\n"
                f = open(self.warn_path, 'a')
                f.write(errortext)
                f.close()
            else:
                errortext = "###############################################################################\n"
                f = open(self.warn_path, 'a')
                f.write(errortext)
                f.close()
                print(f"\n\nAt least one warning occured, please check the warninglog:\n{self.warn_path}\n\n")
        if self.path is None:
            return
        if self.error_counter == 0:
            errortext = "Finished without any error\n###############################################################################\n"
            f = open(self.path, 'a')
            f.write(errortext)
            f.close()
        else:
            errortext = "###############################################################################\n"
            f = open(self.path, 'a')
            f.write(errortext)
            f.close()
            print(f"\n\nAt least one error occured, please check the errorlog:\n{self.path}\n\n")

    def add_warning(self, warn_text, lock=None):
        self.warning_counter += 1
        g = ''.join(traceback.format_list(traceback.extract_stack()))
        warn_text = 'Warning %s:\n%s\n%s\n' % (str(self.warning_counter), warn_text, g)
        if self.warn_path is None or self.print_all_warns:
            print(warn_text)
            return
        if lock is not None:
            with lock:
                f = open(self.warn_path, 'a')
                f.write(warn_text)
                f.close()
            return
        f = open(self.warn_path, 'a')
        f.write(warn_text)
        f.close()

structure_id_types = [
    'PDB',
    'Asymmetric unit PDB',
    'Alphafold model'
]

def identify_structure_id_type(structure_id):
    structe_id_type_key = identify_structure_id_type_key(structure_id)
    if structe_id_type_key is None:
        return 'Unknown structure identifier'
    return structure_id_types[structe_id_type_key]


def identify_structure_id_type_key(structure_id):
    l = len(structure_id)
    if l == 7:
        if structure_id[-3:] == '_AU':
            return 1

    if len(structure_id) == 4:
        return 0

    if structure_id[:3] == 'AF-':
        return 2

    return None

def is_alphafold_model(structure_id):
    structe_id_type_key = identify_structure_id_type_key(structure_id)
    return (structe_id_type_key == 2)

def alphafold_model_id_to_file_path(model_id, config):
    uniprot_ac = model_id.split('-')[1]
    topfolder_id = uniprot_ac[-2:]
    subfolder_id = uniprot_ac[-4:]

    folder_path = f'{config.path_to_model_db}/{topfolder_id}/{subfolder_id}'

    for fn in os.listdir(folder_path):
        if fn.count(model_id) == 1 and fn.count('.pdb.gz') == 1:
            return f'{config.path_to_model_db}/{topfolder_id}/{subfolder_id}/{fn}'

    config.errorlog.add_warning(f'Did not find file {model_id} in alphafold model db: {config.path_to_model_db}')
    return None

def is_alphafold_db_valid(config):
    if not os.path.exists(config.path_to_model_db):
        return False
    if not os.path.exists(f'{config.path_to_model_db}/00'):
        return False
    return True


def resolve_path(path):
    return os.path.abspath(os.path.realpath(path))


def get_conda_infos():
    p = subprocess.Popen(["conda", "env", "list"], stdout=subprocess.PIPE)
    out, _ = p.communicate()
    for line in out.decode('ascii').split('\n'):
        if line == '':
            continue
        
        if line[0] == '#':
            continue

        words = line.split()
        if len(words) == 3:
            current_env = True
            env_name = words[0]
            env_path = words[2]
            current_env_name = env_name
            current_env_path = env_path
        elif len(words) == 2:
            current_env = False
            env_name = words[0]
            env_path = words[1]
        else:
            env_name = None
        if env_name == 'base':
            base_path = env_path

    return base_path, current_env_name, current_env_path

def set_modeller_key(modeller_key):
    p = subprocess.Popen("conda list | grep -w modeller", stdout=subprocess.PIPE, shell=True)
    out, _ = p.communicate()
    words = out.decode('ascii').split()

    modeller_version = words[1]
    _, _, env_path = get_conda_infos()

    modeller_config_path = f'{env_path}/lib/modeller-{modeller_version}/modlib/modeller/config.py'
    f = open(modeller_config_path, 'r')
    lines = f.readlines()
    f.close()

    newlines = []
    for line in lines:
        words = line.split()
        if words[0] != 'license':
            newlines.append(line)
        else:
            license_line = f"license = r'{modeller_key}'"
            newlines.append(license_line)
    
    f = open(modeller_config_path, 'w')
    f.write(''.join(newlines))
    f.close()

def median(l):
    n = len(l)
    l = sorted(l)
    if n == 1:
        return l[0]
    if n == 0:
        return None
    if n % 2 == 0:
        med = (l[(n // 2) - 1] + l[n // 2]) / 2.0
    else:
        med = l[(n - 1) // 2]
    return med


def distance(coord1, coord2):
    diff = [coord1[0] - coord2[0], coord1[1] - coord2[1], coord1[2] - coord2[2]]
    return math.sqrt(diff[0]**2.0 + diff[1]**2.0 + diff[2]**2.0)


def calc_checksum(filename):
    try:
        with open(filename, 'rb') as file:
            return adler32(file.read())
    except FileNotFoundError:
        return adler32(bytes(filename, 'utf-8'))

def calculate_chunksizes(n_of_chunks, n_of_items):
    if n_of_chunks == 0:
        n_of_chunks = 1
    small_chunksize = n_of_items // n_of_chunks
    big_chunksize = small_chunksize + 1
    n_of_small_chunks = n_of_chunks * big_chunksize - n_of_items
    n_of_big_chunks = n_of_chunks - n_of_small_chunks
    if n_of_big_chunks == 0:
        big_chunksize = 0
    return small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks


def serialize(item):
    """Serialization using MessagePack

    Args:
        item: Any item (int, str, list, dict) to serialize

    Returns: Serialized object

    """
    return msgpack.dumps(item)

def deserialize(item):
    """Deserialization using MessagePack

    Args:
        item: Any MessagePack serialized object

    Returns: Deserialized object

    """
    return msgpack.loads(item)


def compress(value):
    """Serialize and Zstandard compress the value

    Args:
        value: Any value, could be list, dict, string, int

    Returns: Compressed serialized bytes

    """
    return zstd.compress(value, 9, 2)

def decompress(compressed_value):
    """Zstandard decompress and deserialize the compressed value

    Args:
        compressed_value: Any bytes that was compressed using 'compress' function

    Returns: Decompressed and deserialized value if not None, else empty list as default value is returned

    """
    if compressed_value is None:
        return None
    return zstd.decompress(compressed_value)

def pack(some_object, complete = False):
    #packed_object = compress(pickletools.optimize(pickle.dumps(some_object, protocol = pickle.HIGHEST_PROTOCOL)))
    if complete:
        packed_object = compress(msgpack.packb(some_object, default = custom_encoder_complete))
    else:
        packed_object = compress(msgpack.packb(some_object, default = custom_encoder))
    return packed_object

def unpack(packed_object):
    #some_object = pickle.loads(decompress(packed_object))
    dec_obj = decompress(packed_object)
    if dec_obj is None:
        return None
    some_object = msgpack.unpackb(dec_obj, object_hook = custom_decoder, strict_map_key = False, use_list = False)
    return some_object
