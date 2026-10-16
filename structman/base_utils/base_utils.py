import math
import os
import sys
import time
import traceback
import numpy as np
from zlib import adler32, crc32

import socket
import wget
import requests

from collections import Counter
import zstd
import zstandard
import hashlib
import io
import msgpack
import subprocess
import logging
import http
from numba import njit
from typing import TypeVar, Generic

from structman.base_utils.custom_encoder import custom_encoder, custom_decoder, clear_decode_caches

T = TypeVar('T')

def sizeof_fmt(num, suffix='B'):
    ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def deep_get_size_of(obj):
    if isinstance(obj, list) or isinstance(obj, set):
        size_sum = 0
        for s_obj in obj:
            size_sum += deep_get_size_of(s_obj)
        return size_sum
    
    if isinstance(obj, dict):
        size_sum = 0
        for s_obj in obj.values():
            size_sum += deep_get_size_of(s_obj)
        return size_sum
    
    if isinstance(obj, Slotted_obj):
        return obj.__getsize__()

    return sys.getsizeof(obj)


class Slotted_obj:
    def __serialize__(self) -> list:
        serialized_object = []
        try:
            slot_mask: list[bool] = self.slot_mask
            for slot_number, attribute_name in enumerate(self.__slots__):
                try:
                    if slot_mask[slot_number]:
                        try:
                            serialized_object.append(self.__getattribute__(attribute_name))
                        except AttributeError:
                            serialized_object.append(None)

                except IndexError as e:
                    raise IndexError(f'{e}\n{type(self)} {slot_number=} {attribute_name=}')
        except AttributeError:
            for attribute_name in self.__slots__:
                try:
                    serialized_object.append(self.__getattribute__(attribute_name))
                except AttributeError:
                    serialized_object.append(None)
        return serialized_object
    
    def deactivate_slot_mask(self) -> None:
        slot_mask_bu = []
        for pos in range(len(self.slot_mask)):
            slot_mask_bu.append(self.slot_mask[pos])
            self.__class__.slot_mask[pos] = True
        # The decoder caches the serialised layout per class, and this changes it.
        clear_decode_caches()
        return slot_mask_bu

    def reactivate_slot_mask(self, slot_mask_bu) -> None:
        for pos in range(len(slot_mask_bu)):
            self.__class__.slot_mask[pos] = slot_mask_bu[pos]
        clear_decode_caches()

    def __getsize__(self):
        size_sum = 0
        for attr_name in self.__slots__:
            try:
                obj = self.__getattribute__(attr_name)
            except AttributeError:
                continue
            size_sum += deep_get_size_of(obj)
        return size_sum

    def log_attr_sizes(self, logger, label = '', print_values=False):
        tup_list = []
        for attr_name in self.__slots__:
            obj = self.__getattribute__(attr_name)            
            size = deep_get_size_of(obj)
            if print_values:
                tup_list.append((attr_name, size, obj))
            else:
                tup_list.append((attr_name, size))

        tup_list = sorted(tup_list, key= lambda x :x[1], reverse=True)
        if print_values:
            for attr_name, size, obj in tup_list:
                logger.info(f'{label}{attr_name} {sizeof_fmt(size)} {obj}')
        else:
            for attr_name, size in tup_list[:10]:
                logger.info(f'{label}{attr_name} {sizeof_fmt(size)}')

class SparseArray(Slotted_obj):
    __slots__ = ['ranged_lists']

    def __init__(self):
        self.ranged_lists: list[tuple[int, int, list[any]]] = []

    def insert(self, key: int, value: any) -> None:
        for ranged_list_id, (min_key, max_key, value_list) in enumerate(self.ranged_lists):
            if key >= min_key:
                if key <= max_key: #New element inside one of the ranged lists
                    try:
                        value_list[key-min_key] = value
                    except TypeError:
                        value_list = list(value_list)
                        value_list[key-min_key] = value
                    return
            elif key == (min_key-1): #New element just one index to the left of current ranged list
                if ranged_list_id == 0: #And it is the most left ranged list -> fuse with current ranged list
                    new_list: list[any] = [value]
                    new_list.extend(value_list)
                    self.ranged_lists[0] = (min_key-1, max_key, new_list)
                    return
                elif key == self.ranged_lists[ranged_list_id -1][1]+1: #It is also just one index to the right of the previous ranged list -> fuse lists
                    self.ranged_lists[ranged_list_id -1][2].append(value)
                    self.ranged_lists[ranged_list_id -1][2].extend(value_list)
                    self.ranged_lists[ranged_list_id -1] = (self.ranged_lists[ranged_list_id -1][0], max_key, self.ranged_lists[ranged_list_id - 1][2])
                    self.ranged_lists.pop(ranged_list_id)
                    return
                else: #fuse with current ranged list
                    new_list: list[any] = [value]
                    new_list.extend(value_list)
                    self.ranged_lists[ranged_list_id] = (min_key-1, max_key, new_list)
                    return
            elif ranged_list_id == 0: #New element outside boundarys of current ranged lists to the left side -> spawn new ranged list
                new_ranged_list: tuple[int, int, list[any]] = (key, key, [value])
                new_ranged_lists: list[tuple[int, int, list[any]]] = [new_ranged_list]
                new_ranged_lists.extend(self.ranged_lists)
                self.ranged_lists = new_ranged_lists
                return
            elif key == self.ranged_lists[ranged_list_id -1][1]+1:#New element just one index to the right of previous ranged list -> append it there
                self.ranged_lists[ranged_list_id -1][2].append(value)
                self.ranged_lists[ranged_list_id -1] = (self.ranged_lists[ranged_list_id -1][0], self.ranged_lists[ranged_list_id -1][1] + 1, self.ranged_lists[ranged_list_id - 1][2])
                return
            else: #New element in the middle ground between two ranged lists -> spawn new ranged list
                new_ranged_list: tuple[int, int, list[any]] = (key, key, [value])
                new_ranged_lists: list[tuple[int, int, list[any]]] = self.ranged_lists[:ranged_list_id]
                new_ranged_lists.append(new_ranged_list)
                new_ranged_lists.extend(self.ranged_lists[ranged_list_id:])
                self.ranged_lists = new_ranged_lists
                return
        #If not returned yet, the new element lies to right of current ranged lists
        if len(self.ranged_lists) == 0:
            self.ranged_lists.append((key, key, [value]))
        elif key == self.ranged_lists[-1][1]+1:
            l = self.ranged_lists[-1][2]
            try:
                l.append(value)
            except AttributeError:
                l = list(self.ranged_lists[-1][2])
                l.append(value)
            try:
                self.ranged_lists[-1] = (self.ranged_lists[-1][0], self.ranged_lists[-1][1] + 1, l)
            except TypeError:
                self.ranged_lists = list(self.ranged_lists)
                self.ranged_lists[-1] = (self.ranged_lists[-1][0], self.ranged_lists[-1][1] + 1, l)
        else:
            new_ranged_list: tuple[int, int, list[any]] = (key, key, [value])
            self.ranged_lists.append(new_ranged_list)

    def get(self, key: int) -> any:
        for min_key, max_key, value_list in self.ranged_lists:
            if key >= min_key:
                if key <= max_key:
                    return value_list[key-min_key]
        return None

    def get_keys(self) -> list[int]:
        keys: list[int] = []
        for min_key, max_key, _ in self.ranged_lists:
            keys.extend(range(min_key, max_key+1))
        return keys

    def get_first_item(self) -> any:
        return self.ranged_lists[0][2][0]

    def __len__(self) -> int:
        n: int = 0
        for min_key, max_key, _ in self.ranged_lists:
            n += ((max_key + 1) - min_key)
        return n

    def to_numpy(self) -> np.ndarray:
        conc_array = []
        for _, _, value_list in self.ranged_lists:
            conc_array.append(np.array(value_list))
        conc_array = np.concatenate(conc_array)
        
        return conc_array



class Residue_Map(Slotted_obj, Generic[T]):
    __slots__ = ['value_list', 'value_dict']

    def __init__(self):
        self.value_list = SparseArray()
        self.value_dict = {}

    #def __setitem__(self, key: int | str, value: any) -> None:
    #    return self.add_item(key, value)

    def add_item(self, key: int | str, value: any) -> None:
        if isinstance(key, int):
            self.add_intkey_item(key, value)
        else:
            try:
                int_key = int(key)
                self.add_intkey_item(int_key, value)
            except:
                self.add_strkey_item(key, value)

    def add_intkey_item(self, intkey: int, value: any) -> None:
        self.value_list.insert(intkey, value)

    def add_strkey_item(self, strkey, value):
        self.value_dict[strkey] = value

    #def __getitem__(self, key: int | str) -> any:
    #    return self.get_item(key)

    def get_item(self, key: int | str) -> any:
        if isinstance(key, int):
            return self.value_list.get(key)
        try:
            intkey = int(key)
            value = self.value_list.get(intkey)
            return value
        except:
            try:
                return self.value_dict[key]
            except:
                return None

    def get_keys(self) -> list[int | str]:
        return self.value_list.get_keys() + list(self.value_dict.keys())
    
    def get_ordered_keys(self) -> list[int | str]:
        return self.value_list.get_keys() + list(sorted(self.value_dict.keys()))

    def get_int_keys(self):
        return self.value_list.get_keys()

    def contains(self, key: int | str) -> bool:
        return self.get_item(key) is not None

    def __len__(self) -> int:
        return len(self.value_list) + len(self.value_dict)
    
    def to_numpy(self) -> np.ndarray:
        list_array: np.ndarray = self.value_list.to_numpy()
        dict_keys = sorted(self.value_dict.keys())
        dict_array: np.array = np.empty(len(self.value_dict))
        for key_nr, key in enumerate(dict_keys):
            dict_array[key_nr] = self.value_dict[key]

        full_array = np.concatenate((list_array, dict_array))
        return full_array


# Per-protein rows - sequences and their position db ids - used to be handed to
# the remote workers as dense matrices padded to the longest protein. One titin
# in the input made every other row 34k entries wide, so the padding, not the
# data, dominated: a 20k protein run spent gigabytes storing zeros. These build
# the same rows as a flat concatenation plus offsets, which costs exactly the
# real data and keeps row lookup an O(1) slice.

def build_ragged_int_store(rows: list[np.ndarray], dtype = np.int32) -> tuple[np.ndarray, np.ndarray]:
    offsets: np.ndarray = np.zeros(len(rows) + 1, dtype=np.int64)
    for row_index, row in enumerate(rows):
        offsets[row_index + 1] = offsets[row_index] + len(row)

    if len(rows) > 0:
        values: np.ndarray = np.concatenate(rows, dtype=dtype)
    else:
        values = np.empty(0, dtype=dtype)

    return values, offsets


def ragged_row(store: tuple[np.ndarray, np.ndarray], row_index: int) -> np.ndarray:
    values, offsets = store
    return values[offsets[row_index]:offsets[row_index + 1]]


def build_sequence_store(sequences: list[str]) -> tuple[np.ndarray, np.ndarray]:
    encoded: list[bytes] = [sequence.encode('utf-8') for sequence in sequences]

    offsets: np.ndarray = np.zeros(len(encoded) + 1, dtype=np.int64)
    for seq_index, sequence in enumerate(encoded):
        offsets[seq_index + 1] = offsets[seq_index] + len(sequence)

    blob: np.ndarray = np.frombuffer(b''.join(encoded), dtype=np.uint8)

    return blob, offsets


def sequence_from_store(store: tuple[np.ndarray, np.ndarray], seq_index: int) -> str:
    blob, offsets = store
    return blob[offsets[seq_index]:offsets[seq_index + 1]].tobytes().decode('utf-8')


class Errorlog:
    def __init__(self, path=None, warn_path=None, print_all_errors=False, print_all_warns=False):
        self.path = path
        self.warn_path = warn_path
        self.error_counter = 0
        self.warning_counter = 0
        self.print_all_errors = print_all_errors
        self.print_all_warns = print_all_warns

    def start(self, nfname, session):
        # Which code wrote this log. The errorlog is what gets read first after a crash, and three
        # crashed runs were once blamed on code they never ran - a stale install of another branch.
        try:
            from structman._build_stamp import describe_running
            build_line = describe_running() + '\n'
        except Exception:
            build_line = ''

        if self.warn_path is not None:
            date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
            errortext = "###############################################################################\n%s:\n%s - Session: %s\n%s" % (date, nfname, str(session), build_line)
            f = open(self.warn_path, 'a')
            f.write(errortext)
            f.close()

        if self.path is None:
            return
        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
        errortext = "###############################################################################\n%s:\n%s - Session: %s\n%s" % (date, nfname, str(session), build_line)
        f = open(self.path, 'a')
        f.write(errortext)
        f.close()

    def add_error(self, error_text, logger = None):
        self.error_counter += 1
        g = ''.join(traceback.format_list(traceback.extract_stack()))
        error_text = 'Error %s:\n%s\n%s\n' % (str(self.error_counter), error_text, g)
        if self.path is None or self.print_all_errors or self.print_all_warns:
            print(error_text)
            return
        f = open(self.path, 'a')
        f.write(error_text)
        f.close()
        if logger is not None:
            logger.error(error_text)

    def stop(self) -> int:
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
            return 0
        if self.error_counter == 0:
            errortext = "Finished without any error\n###############################################################################\n"
            f = open(self.path, 'a')
            f.write(errortext)
            f.close()
            return 0
        else:
            errortext = "###############################################################################\n"
            f = open(self.path, 'a')
            f.write(errortext)
            f.close()
            print(f"\n\nAt least one error occured, please check the errorlog:\n{self.path}\n\n")
            return 1

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

def add_to_times(times: list[tuple[float, str]], prev_timepoint: float, label: str = '') -> float:
    now = time.time()
    times.append([now - prev_timepoint, label])
    return now

def aggregate_times(total_times: list[float | list[float]], times: list[float | list[float]], scope = 0) -> list[float | list[float]]:

    if times is None:
        return total_times
    try:
        for pos, entry in enumerate(times):
            if isinstance(entry, tuple):
                entry = list(entry)

            if isinstance(entry, float):
                entry = [entry, '']
                new_scope = False
            elif not isinstance(entry, list):
                raise TypeError(f'{total_times=} {times=} {scope=}')
            elif len(entry) == 2:
                if isinstance(entry[1], str):
                    new_scope = False
                else:
                    new_scope = True
            else:
                new_scope = True

            if pos >= len(total_times):
                d = (1+pos) - len(total_times)
                for _ in range(d):
                    total_times.append([0., ''])

            if not new_scope:
                try:
                    total_times[pos][0] += entry[0]
                    total_times[pos][1] = entry[1]
                except (IndexError, TypeError):
                    total_times[pos] = entry
            else:
                if len(total_times[pos]) == 2:
                    if isinstance(total_times[pos][1], str):
                        total_times[pos] = []
                total_times[pos] = aggregate_times(total_times[pos], times[pos], scope=scope + 1)
    except IndexError as err:
        raise IndexError(f'{err=} {total_times=} {times=} {scope=}')
    except TypeError as err:
        raise TypeError(f'{err=} {total_times=} {times=} {scope=}')

    return total_times

def print_times(times, label = '', logger = None, scope = 0, scope_nrs = None):
    if len(times) == 0:
        if logger is None:
            print(f'empty times: {label}')
        else:
            logger.info(f'empty times: {label}')
        return
    
    if scope_nrs is None:
        scope_nrs = []

    try:
        for pos, entry in enumerate(times):
            if isinstance(entry, float):
                entry = [entry, '']
                new_scope = False
            elif not isinstance(entry, list):
                raise TypeError(f'{times=} {scope=}')
            elif len(entry) == 2:
                if isinstance(entry[1], str):
                    new_scope = False
                else:
                    new_scope = True
            else:
                new_scope = True

            if not new_scope:
                if scope == 0:
                    intend = ''
                    part = f'{pos+1}'
                else:
                    intend = '  '*scope
                    part = '.'.join(scope_nrs) + f'.{pos+1}'
                message = f'{intend}{label} time part {part}: {entry[0]:_} {entry[1]}'
                if logger is None:
                    print(message)
                else:
                    logger.info(message)

            else:
                next_scope_nrs = scope_nrs.copy() + [f'{pos+1}']
                print_times(times[pos], label = label, scope=scope + 1, scope_nrs=next_scope_nrs, logger=logger)
                

    except IndexError as err:
        raise IndexError(f'{err=} {times=} {scope=}')
    except TypeError as err:
        raise TypeError(f'{err=} {times=} {scope=}')

structure_id_types = [
    'PDB',
    'Asymmetric unit PDB',
    'Alphafold model',
    'Complex model'
]

def identify_structure_id_type(structure_id):
    structe_id_type_key = identify_structure_id_type_key(structure_id)
    if structe_id_type_key is None:
        return 'Unknown structure identifier'
    return structure_id_types[structe_id_type_key]


def identify_structure_id_type_key(structure_id: str) -> int | None:
    le = len(structure_id)
    if le == 7:
        if structure_id[-3:] == '_AU':
            return 1

    if le == 4:
        return 0

    if structure_id[:3] == 'AF-':
        return 2
    
    if structure_id.startswith('afdb_'):
        return 2

    if structure_id[-6:] == 'pdb.gz':
        return 3

    if structure_id.count('__') > 0:
        return 3

    return None

def is_alphafold_model(structure_id):
    structe_id_type_key = identify_structure_id_type_key(structure_id)
    return (structe_id_type_key == 2)

def is_model(structure_id):
    structe_id_type_key = identify_structure_id_type_key(structure_id)
    try:
        return (structe_id_type_key >= 2)
    except TypeError:
        return f'Unknown {structure_id=}'

def topfolder_id_to_chunk_id(topfolder_id: str):
    a = topfolder_id[0]
    b = topfolder_id[1]

    if a.isdigit():
        if b.isdigit():
            bi = int(b)
            sec = (bi//2)*2
            return f'{a}{sec}'
        else:
            return f'{a}_'
    else:
        if b.isdigit():
            return f'_{b}'
        else:    
            return '__'

def get_af_model_folder(model_id, config, get_topfolder_id = False):
    try:
        uniprot_ac = model_id.split('-')[1]
    except IndexError:
        uniprot_ac = model_id

    topfolder_id = uniprot_ac[-2:]
    subfolder_id = uniprot_ac[-4:]

    folder_path = f'{config.path_to_model_db}/{topfolder_id}/{subfolder_id}'
    if get_topfolder_id:
        return folder_path, topfolder_id
    return folder_path

def alphafold_model_id_to_db_key(model_id: str):
    id_parts = model_id.split('-')
    main_db: bool = len(id_parts) > 2
    if len(id_parts) == 1:
        afdb_id_core = id_parts[0]
        main_db = True
        chunk_id = afdb_id_core[-2:]
        return afdb_id_core, main_db, chunk_id
    else:
        try:
            afdb_id_core = id_parts[1]
        except IndexError as err:
            raise IndexError(f'alphafold_model_id_to_db_key failed: {model_id=} {err=}')
    if not main_db and len(afdb_id_core) < 11:
        main_db = True
        if not afdb_id_core.endswith('-F1'):
            afdb_id = f'{afdb_id_core}-F1'
        else:
            afdb_id = afdb_id_core
    elif main_db and len(id_parts) > 1:
        afdb_id = f'{afdb_id_core}-{id_parts[2]}'
    else:
        afdb_id = afdb_id_core
    topfolder_id = afdb_id_core[-2:]
    chunk_id = topfolder_id_to_chunk_id(topfolder_id)
    return afdb_id, main_db, chunk_id

def expand_af_model_id(model_id: str):
    if model_id[:3] != 'AF-':
        model_id = f'AF-{model_id}'

    if len(model_id) < 14:
        if not model_id.endswith('-F1'):
            model_id = f'{model_id}-F1'

    return model_id

def alphafold_model_id_to_file_path(model_id: str, config):
    if model_id.startswith('afdb_'):
        if model_id.startswith('afdb_homo_'):
            return 'complex_homo'
        else:
            return 'complex_hetero'
    
    folder_path = get_af_model_folder(model_id, config)

    if not os.path.isdir(folder_path):
        return None
    
    for fn in os.listdir(folder_path):
        if fn.count(model_id) == 1 and (fn.count('.pdb.gz') == 1 or fn.count('.pdb.zst')):
            return f'{folder_path}/{fn}'

    config.errorlog.add_warning(f'Did not find file {model_id} in alphafold model db: {config.path_to_model_db}')
    return None

def is_alphafold_db_valid(config):
    if not os.path.exists(config.path_to_model_db):
        return False
    if not os.path.exists(f'{config.path_to_model_db}/00'):
        return False
    return True

def complex_id_to_file_path(structure_id, config):
    first_uac = structure_id.split('_')[0]
    topfolder_id = first_uac[-2:]
    subfolder_id = first_uac[-4:-2]
    if structure_id[-7:] != '.pdb.gz':
        structure_id = f'{structure_id}.pdb.gz'
    return f'{config.complex_model_db_path}/{topfolder_id}/{subfolder_id}/{structure_id}'


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

@njit
def distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    squared_diff = (coord1 - coord2)**2
    return squared_diff.sum()**0.5

def calc_checksum(filename):
    try:
        with open(filename, 'rb') as file:
            return adler32(file.read())
    except FileNotFoundError:
        return adler32(bytes(filename, 'utf-8'))

def sequence_hash(sequence: str | None) -> int | None:
    """CRC32 of an amino acid sequence, or None when there is no sequence.

    Used to find duplicate sequences without comparing every pair: proteins are
    bucketed by this value and only the members of a bucket get compared for real.
    32 bits means collisions between different sequences are expected, so every
    caller must confirm a match by comparing the sequences themselves.

    None rather than crc32(b'') for the empty case, because 0 is a legitimate hash
    value and folding every sequence-less protein onto it would make them all look
    like duplicates of each other. Anything that is not a string counts as no sequence
    too: a failed fetch leaves the int 0 in Protein.sequence, and hashing that ended a
    whole-transcriptome run in addProtInfos.

    Returns a plain int, never a numpy scalar - see the interface_hash bug.
    """
    if not isinstance(sequence, str) or len(sequence) == 0:
        return None
    return int(crc32(sequence.encode('utf-8')))


def calculate_chunksizes(n_of_chunks: int, n_of_items: int) -> tuple[int, int, int, int]:
    if n_of_chunks == 0:
        n_of_chunks = 1
    small_chunksize: int = n_of_items // n_of_chunks
    big_chunksize: int = small_chunksize + 1
    n_of_small_chunks: int = n_of_chunks * big_chunksize - n_of_items
    n_of_big_chunks: int = n_of_chunks - n_of_small_chunks
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
    return zstd.compress(value, 3, 2)

def decompress(compressed_value):
    """Zstandard decompress and deserialize the compressed value

    Args:
        compressed_value: Any bytes that was compressed using 'compress' function

    Returns: Decompressed and deserialized value if not None, else empty list as default value is returned

    """
    if compressed_value is None:
        return None
    return zstd.decompress(compressed_value)

def pack(some_object):
    #packed_object = compress(pickletools.optimize(pickle.dumps(some_object, protocol = pickle.HIGHEST_PROTOCOL)))
    packed_object = compress(msgpack.packb(some_object, default = custom_encoder))
    return packed_object


# Every zstd frame opens with this. unpack tells its two input forms apart by it, so
# a payload that skipped compression is still read correctly - and so is every blob
# written before pack_frame existed.
_ZSTD_MAGIC: bytes = bytes.fromhex('28b52ffd')


def pack_frame(some_object):
    """pack without the compression pass, for payloads that are already compressed.

    Position_Data is a pair of separately packed blobs, and running zstd over two
    zstd frames spent 17us a position to make the result 10 bytes longer. The framing
    is all that is wanted there, so this leaves the compression out.
    """
    return msgpack.packb(some_object, default = custom_encoder)


def unpack(packed_object):
    #some_object = pickle.loads(decompress(packed_object))
    if packed_object is None:
        return None
    if packed_object[:4] != _ZSTD_MAGIC:
        # Written by pack_frame. A well formed msgpack payload cannot begin with the
        # zstd magic: 0x28 is a positive fixint, so anything following it is trailing
        # data that msgpack itself would reject.
        return msgpack.unpackb(packed_object, object_hook = custom_decoder, strict_map_key = False, use_list = False)
    dec_obj = decompress(packed_object)
    if dec_obj is None:
        return None
    some_object = msgpack.unpackb(dec_obj, object_hook = custom_decoder, strict_map_key = False, use_list = False)
    return some_object

def reset_logger_for_remotes(config):
    if config.verbosity < 8:
        numba_logger = logging.getLogger('numba')
        numba_logger.setLevel(logging.WARNING)
    fl_logger = logging.getLogger('filelock')
    fl_logger.setLevel(logging.WARNING)
    main_logger = logging.getLogger(__name__)
    logging.basicConfig(filename=config.logfile, encoding='utf-8', level=logging.DEBUG)

    # Only rebind when config.logger has nothing of its own, which is the ray worker
    # case: the Config arrives unpickled and its logger has no handlers there, so the
    # basicConfig above is what gives the worker its file.
    #
    # In the driver config.logger is the logger structman_main gave the run's own
    # FileHandler to, and the root handler still points at the default_log.log next to
    # the config file (Config.__init__ ran basicConfig first, so structman_main's is a
    # no-op). Replacing it with a handler-less logger therefore does not lose the
    # messages - it silently moves the rest of the run's log into that default file.
    # This function used to be reachable only from workers; the serial database paths
    # call the same bodies in the driver now.
    if not config.logger.handlers:
        config.logger = main_logger

def is_seq_valid(seq: str):
    
    freq = dict(Counter(seq))
    
    if len(freq) < 5:
        return False, f'Too few different amino acids ({len(freq)})'
    
    return True, None

def is_connected(url):
    if url.count('//') == 0:
        return False, 'Invalid url'
    
    domain = url.split('//')[1].split('/')[0]

    if domain.count(':') == 1:
        domain, port = domain.split(':')
    else:
        if url[:5] == 'https':
            port = 443
        else:
            port = 80
    
    try:
        # connect to the host -- tells us if the host is actually
        # reachable
        socket.create_connection((domain, port))
        return True, None
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        return False, f'{e}\n{f}\n{g}'

def fetch_af_msa(config, model_id, folder_path, current_afdb_version = 'v6', tmp_target=None, encode_it = False, retry = 0):
    full_id = None
    if os.path.isdir(folder_path):
        for fn in os.listdir(folder_path):
            if fn.count(model_id) == 1 and fn.count('-model') > 0:
                full_id = fn.split('-model')[0]
                break

        if full_id is None:
            config.errorlog.add_warning(f'Did not find file {model_id} in alphafold model db: {config.path_to_model_db}')
            return None

    else:
        full_id = expand_af_model_id(model_id)

    if config.verbosity >= 6:
        config.logger.info(f'Call of fetch_af_msa: {model_id=} {folder_path=} {full_id=} {tmp_target=} {encode_it=} {retry=}')


    af_msa_filename = f'{full_id}-msa_{current_afdb_version}.a3m'

    if config.verbosity >= 6:
        print(f'Fetching MSA for {model_id} from AFDB')

 
    url = f'https://alphafold.ebi.ac.uk/files/msa/{af_msa_filename}'
    connected, error = is_connected(url)
    if not connected:
        config.errorlog.add_warning(f'Did not find file {model_id} in online alphafold model db: {url} {error=}')
        return None
    
    try:
        r = requests.head(url)
    except (requests.exceptions.SSLError, requests.exceptions.ConnectionError) as err:
        config.errorlog.add_warning(f'requests failed for {model_id=} {full_id=} msa in online alphafold model db: {url} {err=}')
        return None

    if tmp_target is None:
        target_folder = folder_path
    else:
        target_folder = tmp_target.rsplit('/',1)[0]

    if r.status_code < 300:
        try:
            wget.download(url, out=target_folder)
        except (OSError, ValueError, http.client.IncompleteRead) as err:
            config.errorlog.add_warning(f'wget failed for {model_id=} {full_id=} msa in online alphafold model db: {url} {err=}')
            return None
    elif r.status_code == 429:
        if retry == 0:
            return None
        retry += 1
        time.sleep(retry**2)
        return fetch_af_msa(config, model_id, folder_path, tmp_target=tmp_target, encode_it = encode_it, retry = retry)
    else:
        config.errorlog.add_warning(f'Invalid http code {r.status_code} trying to fetch {model_id} msa in online alphafold model db: {url}')
        return None

    if tmp_target is None:
        if not encode_it:
            os.rename(f'{folder_path}/{af_msa_filename}', f'{folder_path}/{model_id}_msa.fasta')
        else:
            cctx = zstandard.ZstdCompressor()
            with open(f'{folder_path}/{af_msa_filename}', 'rb') as ifh, open(f'{folder_path}/{af_msa_filename}.zst', 'wb') as ofh:
                cctx.copy_stream(ifh, ofh)
    else:
        os.rename(f'{target_folder}/{af_msa_filename}', tmp_target)

def decode_zstd(encoded):
    dctx = zstandard.ZstdDecompressor()
    with io.BytesIO(encoded) as compr_stream:
        with dctx.stream_reader(compr_stream) as reader:
            decoded = reader.read()

    return decoded

def struct_id_to_bin_nr(struct_id: str, nr_of_bins: int = 20):
    # 1. Create a deterministic SHA-256 hash of the string
    sha256_hash = hashlib.sha256(struct_id.encode()).hexdigest()
    
    # 2. Convert the hexadecimal hash into a base-16 integer
    hash_as_int = int(sha256_hash, 16)
    
    return hash_as_int % nr_of_bins


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
        version_nr = -1
        splits = raw_af_id[3:].split('_', 1)
        af_id = splits[0]
        if len(splits) > 1:
            fn_add = splits[1]

        #print(f'No version number: {raw_af_id=} {af_id=}')

    if len(af_id) < 11:
        if not af_id.endswith('-F1'):
            af_id = f'{af_id}-F1'

    ft = 78

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
        
        case 'csv':
            ft = 7


    if compr is None and ft != 7:
        version_nr = -1
    elif ft == 3 or ft == 7:
        version_nr = -1

    return af_id, ft, compr, version_nr

def create_mmseqs_index(infile, outfolder, mmseqs2_tmp):
    p = subprocess.Popen(['mmseqs', 'createdb', infile, f'{infile}_search_db'], cwd=outfolder)
    p.wait()

    p = subprocess.Popen(['rm', '-R', mmseqs2_tmp], cwd=outfolder)
    p.wait()

    if not os.path.isdir(mmseqs2_tmp):
        os.mkdir(mmseqs2_tmp)

    p = subprocess.Popen(['chmod', '777', '-R', mmseqs2_tmp], cwd=outfolder)
    p.wait()

    p = subprocess.Popen(['mmseqs', 'createindex', f'{infile}_search_db', mmseqs2_tmp, '-s', '7.5'], cwd=outfolder)
    p.wait()