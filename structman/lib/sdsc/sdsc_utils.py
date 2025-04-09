import sys
import os
import traceback
import numpy as np
import socket
import time
from numba.experimental import jitclass
from numba.typed import List as NumbaList
from numba.types import ListType, int64
#from numba import typeof


try:
    from reprlib import repr
except ImportError:
    pass

from structman.lib.sdsc.consts import ligands, residues
from structman.base_utils.base_utils import identify_structure_id_type_key

class Slotted_obj:
    def __serialize__(self) -> dict[str, any]:
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

class SparseArray(Slotted_obj):
    __slots__ = ['ranged_lists']

    def __init__(self):
        self.ranged_lists: list[tuple[int, int, list[any]]] = []

    def insert(self, key: int, value: any) -> None:
        for ranged_list_id, (min_key, max_key, value_list) in enumerate(self.ranged_lists):
            if key >= min_key:
                if key <= max_key: #New element inside one of the ranged lists
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
                    self.ranged_lists[0] = (min_key-1, max_key, new_list)
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

"""
spec = [
    ('ranged_lists', ListType((int64, int64, ListType((str, int64)))))
]
@jitclass()
class NBSparseArray(object):

    def __init__(self):
        a = NumbaList.empty_list((str, int64))

        self.ranged_lists = NumbaList([(0,5,a)])
        #self.ranged_lists = NumbaList.empty_list((int32, int32, ListType((str, int32))))

    def insert(self, key: int, value: any) -> None:
        for ranged_list_id, (min_key, max_key, value_list) in enumerate(self.ranged_lists):
            if key >= min_key:
                if key <= max_key: #New element inside one of the ranged lists
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
                    self.ranged_lists[0] = (min_key-1, max_key, new_list)
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
            self.ranged_lists[-1][2].append(value)
            self.ranged_lists[-1] = (self.ranged_lists[-1][0], self.ranged_lists[-1][1] + 1, self.ranged_lists[-1][2])
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
"""
        
def invert_deletion_flanks(del_flanks):
    wt_inserts = []
    for lf, rf in del_flanks:
        insert = []
        for i in range(lf, (rf+1)):
            insert.append(i)
        wt_inserts.append(insert)
    return wt_inserts

def invert_insertion_positions(insertion_positions):
    wt_del_flanks = []
    for insert in insertion_positions:
        lf = min(insert)
        rf = max(insert)
        wt_del_flanks.append((lf, rf))
    return wt_del_flanks

def translate(nuc_seq):
    codon_counter = 0
    aa_seq = ''
    while True:
        if codon_counter >= len(nuc_seq):
            break
        codon = nuc_seq[codon_counter:codon_counter + 3]
        codon_counter += 3
        aa = codons.CODONS[codon]
        if aa == '_':
            break
        aa_seq += aa
    return aa_seq


def doomsday_protocol(obj, immunity = None):
    for attribute_name in obj.__slots__:
        try:
            if immunity is not None:
                if attribute_name in immunity:
                    continue
            delattr(obj, attribute_name)
        except:
            pass

            
def boring(abr: str) -> bool:
    if abr in residues.METAL_ATOMS or abr in residues.ION_ATOMS or abr in ligands.NON_BORING_SHORT_LIGANDS:
        return False
    if abr in ligands.BORING_LIGANDS:
        return True
    if len(abr) < 3:
        return True
    return False


def bin_boring(abr: bytes) -> bool:
    if abr in residues.BIN_METAL_ATOMS or abr in residues.BIN_ION_ATOMS or abr in ligands.BIN_NON_BORING_SHORT_LIGANDS:
        return False
    if abr in ligands.BIN_BORING_LIGANDS:
        return True
    if len(abr) < 3:
        return True
    return False


def is_mutant_ac(ac):
    if not ac.count('del') == 1 and not ac.count('ins') == 1:
        return False
    else:
        return True


def locate(rsa, config, binary_decision=False):
    if rsa is None:
        return None
    else:
        if rsa > config.surface_threshold:
            return "Surface"
        elif rsa > config.buried_threshold:
            if binary_decision:
                return "Core"
            return "Buried"
        else:
            return "Core"


def triple_locate(rsa, mc_rsa, sc_rsa, config, binary_decision=False):
    loc = locate(rsa, config, binary_decision=binary_decision)
    mc_loc = locate(mc_rsa, config, binary_decision=binary_decision)
    sc_loc = locate(sc_rsa, config, binary_decision=binary_decision)
    return (loc, mc_loc, sc_loc)


def parseFasta(path=None, new_file=None, lines=None, page=None, left_split=None, right_split=' '):
    if lines is None and page is None:
        f = open(path, 'r')
        lines = f.read().split('\n')
        f.close()
    elif lines is None:
        lines = page.split('\n')

    seq_map = {}
    n = 0

    if new_file is not None:
        new_lines = []

    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '>':
            entry_id = line[1:]
            if left_split is not None:
                entry_id = entry_id.split(left_split, 1)[1]
            if right_split is not None:
                entry_id = entry_id.split(right_split, 1)[0]
            seq_map[entry_id] = ''
            n += 1
            if new_file is not None:
                new_lines.append(line)
        else:
            seq_map[entry_id] += line
            if new_file is not None:
                new_lines.append(line)

    if new_file is not None:
        f = open(new_file, 'w')
        f.write('\n'.join(new_lines))
        f.close()

    return seq_map


def process_recommend_structure_str(recommended_structure_str):
    if recommended_structure_str is not None and recommended_structure_str != '-':
        words = recommended_structure_str.split(';')
        if len(words) < 4:
            resolution = '-'
            cov = '-'
            seq_id = '-'
            recommended_structure = '-'
        else:
            recommended_structure, seq_id, cov, resolution = words
    else:
        resolution = '-'
        cov = '-'
        seq_id = '-'
        recommended_structure = '-'
    return recommended_structure, seq_id, cov, resolution


def classify(interaction_profile, location):
    if interaction_profile is None:
        return None, None
    raw_rin_class, raw_rin_simple_class = interaction_profile.getClass()
    if raw_rin_class == 'No interaction':
        rin_class = location
    else:
        rin_class = raw_rin_class
    if raw_rin_simple_class == 'No interaction':
        rin_simple_class = location
    else:
        rin_simple_class = raw_rin_simple_class
    return rin_class, rin_simple_class


def process_alignment_data(alignment):
    if alignment is None:
        return 'Cannot process None'
    lines = alignment.split("\n")
    nextline = False
    target_start = False
    template_start = False
    target_lines = []
    template_lines = []
    target_name = ""
    for line in lines:
        if len(line) > 0:
            if line[0] == ">":
                ids = line.split(";")
                if target_name == "":
                    target_name = ids[1]
                    target_name = target_name.replace(" ", "")
                    target_name = target_name.replace("\n", "")
                nextline = True
            elif nextline:
                if not target_start:
                    target_start = True
                else:
                    target_start = False
                    template_start = True
                    words = line.split(":")
                    startres = words[2]
                    endres = words[4]
                    chain = words[3]
                nextline = False
            elif line[0] == "\n":
                template_start = False
            elif target_start:
                target_lines.append(line)
            elif template_start:
                template_lines.append(line)

    target_seq = "".join(target_lines)
    target_seq = target_seq.replace("*", "")
    template_seq = "".join(template_lines)
    template_seq = template_seq.replace("*", "")
    return target_seq, template_seq


def get_shortest_distances(chains, lig_dists, chain_distances, homomer_distances):
    min_ld = None
    min_md = None
    min_id = None
    min_cd = None
    min_dd = None
    min_rd = None
    min_hd = None

    ldists = {}
    mdists = {}
    idists = {}
    if lig_dists is not None:
        for lig_id in lig_dists:
            lig_name, res, chain = lig_id.split('_')
            (dist, atom_pair) = lig_dists[lig_id]
            if lig_name in residues.METAL_ATOMS:
                mdists[dist] = lig_name, res, chain
            elif lig_name in residues.ION_ATOMS:
                idists[dist] = lig_name, res, chain
            else:
                ldists[dist] = lig_name, res, chain

    min_lig = None
    min_metal = None
    min_ion = None

    if len(ldists) > 0:
        min_ld = min(ldists.keys())
        min_lig = ldists[min_ld]

    if len(mdists) > 0:
        min_md = min(mdists.keys())
        min_metal = mdists[min_md]

    if len(idists) > 0:
        min_id = min(idists.keys())
        min_ion = idists[min_id]

    cdists = {}
    ddists = {}
    rdists = {}

    iacs = {}

    if chain_distances is not None:
        for chain_id in chain_distances:
            (dist, atom_pair, min_resi) = chain_distances[chain_id]
            if dist is None:
                continue
            if chain_id not in chains:
                chaintype = 'Protein'
            else:
                chaintype = chains[chain_id]

            if chaintype == "Protein" or chaintype == 'Peptide':
                cdists[dist] = chain_id
            elif chaintype == "RNA":
                rdists[dist] = chain_id
            elif chaintype == "DNA":
                ddists[dist] = chain_id

    if len(cdists) > 0:
        min_cd = min(cdists.keys())
        iacs['Protein'] = cdists[min_cd]
    if len(rdists) > 0:
        min_rd = min(rdists.keys())
        iacs['RNA'] = rdists[min_rd]
    if len(ddists) > 0:
        min_dd = min(ddists.keys())
        iacs['DNA'] = ddists[min_dd]

    homo_dists = []
    if homomer_distances is not None:
        for homo_chain in homomer_distances:
            dist = homomer_distances[homo_chain]
            homo_dists.append(dist)
    if len(homo_dists) > 0:
        min_hd = min(homo_dists)

    minimal_distances = []
    if min_cd is not None:
        minimal_distances.append(min_cd)
    if min_dd is not None:
        minimal_distances.append(min_dd)
    if min_rd is not None:
        minimal_distances.append(min_rd)
    if min_ld is not None:
        minimal_distances.append(min_ld)
    if min_md is not None:
        minimal_distances.append(min_md)
    if min_id is not None:
        minimal_distances.append(min_id)

    if len(minimal_distances) == 0:
        min_minimal_distances = 2.0
    else:
        min_minimal_distances = min(minimal_distances)

    if min_minimal_distances < 1.2:
        return None

    return min_hd, min_ld, min_md, min_id, min_cd, min_rd, min_dd, min_lig, min_metal, min_ion, iacs


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


def connection_sleep_cycle(verbosity, url):
    connected, error = is_connected(url)
    if connected:
        return
    while not connected:
        connected, error = is_connected(url)
        if verbosity >= 1:
            print(f'No connection: {error}\n Sleeping a bit and then try again')
        time.sleep(30)

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