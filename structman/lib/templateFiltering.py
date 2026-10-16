import math
import os
import subprocess
import sys
import time
import gc
import psutil
from itertools import cycle
import numpy as np
from numba import njit
from numba.core import types as numba_types
from numba.typed import List as NumbaList
from numba.typed import Dict as NumbaDict

import ray

from structman.lib import pdbParser
from structman.lib.rin import (
    lookup, Interaction_profile, Inflow_values,
    profile_degree_index, profile_matrix_from_profiles, total_interface_scores_from_matrix,
    )
from structman.lib.database.insertion_lib import remote_insertResidues, insertComplexes
from structman.lib.database.retrieval import retrieve_structure_info_from_struct_ids
from structman.base_utils.base_utils import (
    median, distance,
    is_alphafold_model,
    is_model, complex_id_to_file_path,
    aggregate_times, print_times, add_to_times,
    reset_logger_for_remotes,
    Residue_Map, struct_id_to_bin_nr, pack, Slotted_obj
)
from structman.base_utils.config_class import Config
from structman.base_utils.ray_utils import parallel_db_access, ray_init
from structman.base_utils.chain_consts import (
    is_usable_chain_id, chain_id_number, MAX_NUMBER_OF_CHAINS, CHAIN_ID_ORDER)
from structman.lib.sdsc.consts.residues import CORRECT_COUNT, BIN_THREE_TO_ONE, THREE_TO_ONE, BLOSUM62, ONE_TO_THREE, RESIDUE_MAX_ACC, HYDROPATHY, METAL_ATOMS, ION_ATOMS
from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc import residue as residue_package
from structman.lib.sdsc import interface as interface_package
from structman.lib.aggregation import classification


# called by database
def candidateScore(lig_sub_dist, chain_sub_dist, lig_wf=1.0, chain_wf=1.0, useBlosum=False, aac="", blosum_wf=0.5):
    if useBlosum:
        if aac == "":
            raise NameError("If useBlosum is True, an aac is needed")

        try:
            blosum_value = 0.6 - float(BLOSUM62[(aac[0], aac[-1])]) / 10
        except KeyError:
            blosum_value = 0.6 - float(BLOSUM62[(aac[-1], aac[0])]) / 10
        blosum_value = max(blosum_value, 0.0)
    else:
        blosum_value = 1.0
        blosum_wf = 0.0

    # project the criteria to [0,1] via logistic regression
    lig_value = (1 + math.exp(lig_sub_dist - 10))**(-1)
    chain_value = (1 + math.exp(chain_sub_dist - 10))**(-1)

    if lig_sub_dist == -1:
        lig_value = 0.0
    if chain_sub_dist == -1:
        chain_value = 0.0

    ws = sum((lig_wf, chain_wf, blosum_wf))

    candidate = sum((lig_value * lig_wf, chain_value * chain_wf, blosum_value * blosum_wf)) / ws
    return candidate


def calcDSSP(path, DSSP, angles=False):

    dssp_dict: dict[str, Residue_Map[tuple[float, float, float, str]]] = {}
    errorlist = []

    try:
        p = subprocess.Popen([DSSP, path], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
    except OSError:
        out = ''
        err = f'{sys.exc_info()[0]!s} {sys.exc_info()[1]!s}'
    # Alert user for errors
    if err.strip() and not out.strip():
        
        errorlist.append(f'DSSP failed to produce an output\n{path}\n{err}\n')
        return dssp_dict, errorlist

    lines = out.split('\n')

    i = 1

    for line in lines:
        if line == '':
            continue
        words = line.split()
        if words[0] == '#':
            break
        i += 1

    for line in lines[i:]:
        if line == '':
            continue
        if len(line) < 46:
            continue
        if line[9] == " ":
            continue

        res = (line[5:11]).replace(" ", "")  # this includes the insertion code
        
        chain = line[11]
        aa_type_one_letter = line[13]
        ssa = line[16]
        racc = None
        relative_main_chain_acc = None
        relative_side_chain_acc = None

        try:
            acc = float(line[34:38])
            main_chain_acc = float(line[38:42])
            side_chain_acc = float(line[42:46])
        except ValueError:
            errorlist.append(f'DSSP float conversion failed for: \n{path}\nChain: {chain}, Residue: {res},{aa_type_one_letter}\n{line}\n\n')
            acc = None

        if acc is not None:
            try:
                if aa_type_one_letter == 'X':
                    macc = 220.0
                elif aa_type_one_letter.islower():  # This happens for SS-bond cysteines
                    aa = ONE_TO_THREE['C']
                    macc = RESIDUE_MAX_ACC['Sander'][aa]
                else:
                    aa = ONE_TO_THREE[aa_type_one_letter]
                    macc = RESIDUE_MAX_ACC['Sander'][aa]
                racc = acc / macc
                relative_main_chain_acc = main_chain_acc / macc
                relative_side_chain_acc = side_chain_acc / macc
            except (KeyError, ZeroDivisionError):
                errorlist.append(f'DSSP ({DSSP}) failed for: \n{path}\nChain: {chain}, Residue: {res},{aa_type_one_letter}\n{line}\n\n')
                racc = None
                relative_main_chain_acc = None
                relative_side_chain_acc = None

        phi = None
        psi = None
        tco = None
        kappa = None
        alpha = None

        try:
            if len(line) > 123:
                tco = float(line[93:99].replace(" ", ""))
                kappa = float(line[99:105].replace(" ", ""))
                alpha = float(line[105:111].replace(" ", ""))
                phi = float(line[111:117].replace(" ", ""))
                psi = float(line[117:123].replace(" ", ""))
                
        except ValueError:
            errorlist.append(f'DSSP failed for (angles): \n{path}\nChain: {chain}, Residue: {res},{aa_type_one_letter}\n{line}\n\n')

        if chain not in dssp_dict:
            dssp_dict[chain] = Residue_Map()
        dssp_dict[chain].add_item(res, (racc, relative_main_chain_acc, relative_side_chain_acc, ssa))
        if angles:
            dssp_dict[chain].add_item(res, (racc, relative_main_chain_acc, relative_side_chain_acc, ssa, phi, psi, tco, kappa, alpha))
    return dssp_dict, errorlist



# --- the parsed structure ----------------------------------------------------
#
# parsePDB used to hand back ten dicts keyed by chain id. Everything downstream
# then re-looked-up the same chain in each of them, and nothing keyed that way can
# be reached from a kernel. They are columns indexed by the chain's position in
# CHAIN_ID_ORDER now - the same 62 slot alphabet the rest of the pipeline already
# uses - with chain_numbers listing the chains that are actually present, in the
# order the parser met them, because several consumers depend on that order.
#
# The maps whose keys are residue numbers rather than chain ids (modres, ssbond,
# link, the two cis maps, res_contig) stay dicts; only the chain level changes.

_PS_CHAIN_TYPE_NAMES = ('Protein', 'RNA', 'DNA', 'Peptide', 'MISC')
_PS_PROTEIN, _PS_RNA, _PS_DNA, _PS_PEPTIDE, _PS_MISC = range(5)
_PS_CHAIN_TYPE_CODES = {name: code for code, name in enumerate(_PS_CHAIN_TYPE_NAMES)}


class Parsed_structure(Slotted_obj):
    """Everything parsePDB extracts, as chain-number-indexed columns."""

    __slots__ = [
        'chain_numbers', 'residue_names', 'residue_atoms', 'hetero_names',
        'hetero_atoms', 'box', 'chain_type', 'b_factors', 'res_contig', 'modres',
        'ssbond', 'link', 'cis_conformation', 'cis_follower', 'ligands', 'metals',
        'ions']

    def __init__(self):
        slots = MAX_NUMBER_OF_CHAINS
        self.chain_numbers = np.empty(0, dtype=np.int64)
        self.residue_names = [None] * slots
        self.residue_atoms = [None] * slots
        self.hetero_names = [None] * slots
        self.hetero_atoms = [None] * slots
        self.box = np.full((slots, 6), np.nan, dtype=np.float32)
        self.chain_type = np.full(slots, -1, dtype=np.int8)
        self.b_factors = [None] * slots
        self.res_contig = [None] * slots
        self.modres = [None] * slots
        self.ssbond = [None] * slots
        self.link = [None] * slots
        self.cis_conformation = [None] * slots
        self.cis_follower = [None] * slots
        self.ligands = set()
        self.metals = set()
        self.ions = set()

    def chain_id(self, chain_number: int) -> str:
        return CHAIN_ID_ORDER[chain_number]

    def chain_id_number_of(self, chain_id: str) -> int:
        return chain_id_number(chain_id)

    def chain_ids(self) -> list[str]:
        return [CHAIN_ID_ORDER[number] for number in self.chain_numbers.tolist()]

    def is_protein(self, chain_number: int) -> bool:
        return self.chain_type[chain_number] == _PS_PROTEIN

    def type_name(self, chain_number: int) -> str | None:
        code = self.chain_type[chain_number]
        if code == -1:
            return None
        return _PS_CHAIN_TYPE_NAMES[code]

    # --- the dict form, for the consumers that still want it -----------------
    #
    # rin.lookup reads RIN files and is keyed by chain id throughout; it is called
    # once per structure and does file i/o, so it keeps its interface and gets the
    # dicts rebuilt here rather than being dragged into the column form.

    def chain_type_map(self) -> dict[str, str]:
        # Every slot, not just chain_numbers: a chain that only ever appeared in a
        # MODRES record has a type but no coordinates, and the dict this replaces
        # held it.
        return {CHAIN_ID_ORDER[number]: _PS_CHAIN_TYPE_NAMES[self.chain_type[number]]
                for number in range(MAX_NUMBER_OF_CHAINS)
                if self.chain_type[number] != -1}

    def protein_chain_numbers(self) -> list[int]:
        return [number for number in range(MAX_NUMBER_OF_CHAINS)
                if self.chain_type[number] == _PS_PROTEIN]

    def res_contig_map(self) -> dict:
        return {CHAIN_ID_ORDER[number]: self.res_contig[number]
                for number in self.chain_numbers.tolist()
                if self.res_contig[number] is not None}

    def chain_keyed(self, values: set) -> set:
        """A {(chain number, residue)} set back as {(chain id, residue)}."""
        return {(CHAIN_ID_ORDER[number], residue) for number, residue in values}


def contig_map_to_list_dict(
        res_contig_map: dict[str, Residue_Map[tuple[int, str]]],
        chain_type_map: dict[str, str]
        ) -> dict[str, list[None | str]]:
    contig_lists: dict[str, list[None | str] | None] = {}

    for chain in res_contig_map:
        if chain_type_map[chain] != 'Protein':
            contig_lists[chain] = None
            continue
        try:
            keys: list[int | str] = res_contig_map[chain].get_keys()
            contig_list: list[None | str] = [None] * (len(keys) + 1)
            for key in keys:
                res_contig, _ = res_contig_map[chain].get_item(key)
                contig_list[res_contig] = key
            contig_lists[chain] = contig_list
        except IndexError as err:
            contigs = []
            for key in res_contig_map[chain].get_keys():
                res_contig, _ = res_contig_map[chain].get_item(key)
                contigs.append(res_contig)
            raise IndexError(f'contig_map_to_list_dict failed for {chain=} {chain_type_map=} {err=} {len(keys)=} {len(contigs)=} {sorted(contigs)=}')

    return contig_lists

def safe_int(value: str) -> int | str:
    try:
        return int(value)
    except ValueError:
        return value

def reindex_by_contig(res_nr_map: dict, res_contig_map_chain: Residue_Map) -> dict:
    reindexed_map = {}
    for res_nr in res_nr_map:
        try:
            res_index = res_contig_map_chain.get_item(res_nr)[0] - 1
        except TypeError:
            continue
        reindexed_map[res_index] = res_nr_map[res_nr]
    return reindexed_map

def parsePDB(input_page: str) -> "Parsed_structure":
    """
    Parses a PDB-file and takes all atomic coordinates.

    Input:
    input_page: String ; content of a pdb file

    Output:
    coordinate_map: ,{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
    """

    lines = input_page.split('\n')

    # Chain-indexed columns rather than dicts keyed by chain id: chain ids are
    # single characters from the fixed CHAIN_ID_ORDER alphabet, so 62 slots cover
    # any structure, and a chain's slot is a table lookup rather than a string
    # hash on every atom line. chain_order keeps the order chains were first met
    # in, which several consumers depend on.
    slots = MAX_NUMBER_OF_CHAINS
    chain_order: list[int] = []
    coordinate_map: list = [None] * slots
    res_contig_map: list = [None] * slots
    contig_help_map: list[int] = [0] * slots

    modres_map: list = [None] * slots

    ssbond_map: list = [None] * slots
    link_map: list = [None] * slots
    cis_conformation_map: list = [None] * slots
    cis_follower_map: list = [None] * slots

    chain_type_map: list = [None] * slots
    chain_type = '-'

    peptide_count: list = [None] * slots

    all_coords: list = [None] * slots

    # Keyed by chain number, like everything else; Parsed_structure.chain_keyed
    # turns them back into chain ids for the RIN lookup.
    ligands: set[tuple[int, int | str]] = set()
    metals: set[tuple[int, int | str]] = set()
    ions: set[tuple[int, int | str]] = set()
    rare_residues = set()
    b_factors: list = [None] * slots

    firstAltLoc = None

    for line in lines:
        if len(line) > 5:
            record_name: str = line[0:6].replace(" ", "")
            if record_name == "ENDMDL":
                break
            elif record_name == 'SEQRES':
                for tlc in line[19:].split():
                    tlc = tlc.strip()
                    if len(tlc) != 3:
                        continue
                    if tlc not in THREE_TO_ONE:
                        rare_residues.add(tlc)
        # ignore short lines
        if len(line) > 20:

            if record_name.count('ATOM') > 0 and record_name != 'ATOM':  # 100k atom bug fix
                atom_nr: int = int(f'{record_name[4:]}{line[6:11].replace(" ", "")}')
                record_name = 'ATOM'
            else:
                try:
                    atom_nr: int | None = int(line[6:11].replace(" ", ""))
                except ValueError:
                    atom_nr = None

            atom_name = line[12:16].replace(" ", "")
            res_name = line[17:20].replace(" ", "")

            if len(line) > 21:
                if record_name == 'MODRES':
                    chain_id = line[16]
                    res_nr = safe_int(line[18:23].replace(" ", ""))
                    res_name = line[24:27].replace(" ", "")
                    if modres_map[chain_number] is None:
                        modres_map[chain_number] = {}
                    if res_nr not in modres_map[chain_number]:
                        modres_map[chain_number][res_nr] = res_name

                if record_name == 'SSBOND':
                    chain_1 = line[15]
                    res_nr_1 = safe_int(line[17:22].replace(" ", ""))
                    chain_2 = line[29]
                    res_nr_2 = safe_int(line[31:36].replace(" ", ""))

                    ssbond_len = float(line[73:78].replace(' ', ''))
                    # The record names two chains; either may be outside the alphabet.
                    number_1 = chain_id_number(chain_1)
                    number_2 = chain_id_number(chain_2)
                    if number_1 != -1:
                        if ssbond_map[number_1] is None:
                            ssbond_map[number_1] = {}
                        ssbond_map[number_1][res_nr_1] = (chain_2, res_nr_2, ssbond_len)
                    if number_2 != -1:
                        if ssbond_map[number_2] is None:
                            ssbond_map[number_2] = {}
                        ssbond_map[number_2][res_nr_2] = (chain_1, res_nr_1, ssbond_len)

                if record_name == 'LINK  ':
                    atom_1 = line[12:16].replace(" ", "")
                    res_name_1 = line[17:20].replace(" ", "")
                    res_nr_1 = safe_int(line[22:27].replace(" ", ""))
                    chain_1 = line[21]

                    atom_2 = line[42:46].replace(" ", "")
                    res_name_2 = line[47:50].replace(" ", "")
                    res_nr_2 = safe_int(line[52:57].replace(" ", ""))
                    chain_2 = line[51]

                    altLoc_1 = line[16]
                    altLoc_2 = line[46]
                    if firstAltLoc is None and altLoc_1 != ' ':
                        firstAltLoc = altLoc_1  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                    if altLoc_1 != ' ' and altLoc_1 != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue
                    if altLoc_2 != ' ' and altLoc_2 != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue

                    link_dist = float(line[73:78].replace(" ", ""))

                    number_1 = chain_id_number(chain_1)
                    number_2 = chain_id_number(chain_2)
                    if number_1 != -1:
                        if link_map[number_1] is None:
                            link_map[number_1] = {}
                        link_map[number_1][res_nr_1] = (atom_1, res_name_1, atom_2, res_name_2, res_nr_2, chain_2, link_dist)
                    if number_2 != -1:
                        if link_map[number_2] is None:
                            link_map[number_2] = {}
                        link_map[number_2][res_nr_2] = (atom_2, res_name_2, atom_1, res_name_1, res_nr_1, chain_1, link_dist)

                if record_name == 'CISPEP':
                    res_name_1 = line[11:14].replace(" ", "")
                    chain_1 = line[15]
                    res_nr_1 = safe_int(line[17:22].replace(" ", ""))
                    res_name_2 = line[25:28].replace(" ", "")
                    chain_2 = line[29]
                    res_nr_2 = safe_int(line[31:36].replace(" ", ""))
                    angle = float(line[53:59].replace(" ", ""))

                    number_1 = chain_id_number(chain_1)
                    number_2 = chain_id_number(chain_2)
                    if number_1 != -1:
                        if cis_conformation_map[number_1] is None:
                            cis_conformation_map[number_1] = {}
                        cis_conformation_map[number_1][res_nr_1] = (res_name_2, chain_2, res_nr_2, angle)
                    if number_2 != -1:
                        if cis_follower_map[number_2] is None:
                            cis_follower_map[number_2] = {}
                        cis_follower_map[number_2][res_nr_2] = (res_name_1, chain_1, res_nr_1, angle)

                chain_id = line[21]

                # Chain ids are single characters from the fixed CHAIN_ID_ORDER alphabet;
                # anything else (a blank chain id above all) has no chain number and is
                # dropped rather than carried through as a key nothing can index.
                chain_number = chain_id_number(chain_id)
                if chain_number == -1:
                    continue

                if modres_map[chain_number] is None:
                    modres_map[chain_number] = {}

                res_nr = safe_int(line[22:27].replace(" ", ""))  # [22:27] includes the insertion_code

                if record_name == "ATOM" or record_name == "HETATM":
                    altLoc = line[16]
                    if firstAltLoc is None and altLoc != ' ':
                        firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                    if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue

            if len(line) > 60:
                try:
                    b_factor = float(line[61:67].replace(" ", ""))
                except ValueError:
                    b_factor = 0.0
            else:
                b_factor = 0.0  # reset per line so a short line can't reuse the previous line's value
            if record_name == "ATOM" or record_name == 'MODRES' or record_name == "HETATM":
                if chain_type_map[chain_number] is None:
                    if record_name == "ATOM" or record_name == 'MODRES':
                        if len(res_name) == 1:
                            chain_type = "RNA"
                        elif len(res_name) == 2:
                            chain_type = "DNA"
                        elif len(res_name) == 3:
                            chain_type = "Protein"
                        chain_type_map[chain_number] = chain_type
                        peptide_count[chain_number] = [0, 0]
                    elif record_name == 'HETATM':
                        if res_name in THREE_TO_ONE or not sdsc_utils.boring(res_name) or (res_name in rare_residues):
                            chain_type_map[chain_number] = "Protein"
                            peptide_count[chain_number] = [0, 0]
                        elif len(res_name) == 1:
                            chain_type_map[chain_number] = "RNA"
                            peptide_count[chain_number] = [0, 0]
                        elif len(res_name) == 2:
                            chain_type_map[chain_number] = "DNA"
                            peptide_count[chain_number] = [0, 0]
                        else:
                            chain_type_map[chain_number] = 'MISC'


                if record_name == 'HETATM':
                    if res_name in THREE_TO_ONE or not sdsc_utils.boring(res_name) or (res_name in rare_residues):  # For a hetero peptide 'boring' hetero amino acids are allowed as well as other non boring molecules not in THREE_TO_ONE, which are hopefully some kind of anormal amino acids
                        peptide_count[chain_number][1] += 1
                    elif len(res_name) < 3:
                        peptide_count[chain_number][0] += 1

                elif record_name == "ATOM" or record_name == 'MODRES':
                    if res_name != 'UNK':
                        peptide_count[chain_number][0] += 1
                    else:
                        peptide_count[chain_number][1] += 1

                    if len(res_name) == 1:
                        chain_type_map[chain_number] = "RNA"
                    elif len(res_name) == 2:
                        chain_type_map[chain_number] = "DNA"

            if record_name == "ATOM" and len(line) > 50 and atom_name[0] not in ('H', 'D'):
                x = float(line[30:38].replace(" ", ""))
                y = float(line[38:46].replace(" ", ""))
                z = float(line[46:54].replace(" ", ""))
                if coordinate_map[chain_number] is None:
                    coordinate_map[chain_number] = [Residue_Map(), Residue_Map()]
                    chain_order.append(chain_number)
 
                    all_coords[chain_number] = []

                if not coordinate_map[chain_number][0].contains(res_nr):
                    coordinate_map[chain_number][0].add_item(res_nr, (res_name, []))

                coordinate_map[chain_number][0].get_item(res_nr)[1].append((atom_nr, x, y, z))

                if b_factors[chain_number] is None:
                    b_factors[chain_number] = {}
                if res_nr not in b_factors[chain_number]:
                    b_factors[chain_number][res_nr] = []
                b_factors[chain_number][res_nr].append(b_factor)

                all_coords[chain_number].append((x, z, y))

                if res_contig_map[chain_number] is None:
                    res_contig_map[chain_number] = Residue_Map()
                    res_contig_map[chain_number].add_item(res_nr, (1, res_name))
                    contig_help_map[chain_number] = 1
                elif not res_contig_map[chain_number].contains(res_nr):
                    contig_help_map[chain_number] += 1
                    res_contig_map[chain_number].add_item(res_nr, (contig_help_map[chain_number], res_name))

            if record_name == "HETATM" and len(line) > 50:
                x = float(line[30:38].replace(" ", ""))
                y = float(line[38:46].replace(" ", ""))
                z = float(line[46:54].replace(" ", ""))
                if coordinate_map[chain_number] is None:
                    coordinate_map[chain_number] = [Residue_Map(), Residue_Map()]
                    chain_order.append(chain_number)

                if res_nr in modres_map[chain_number] or (res_name in THREE_TO_ONE) or (res_name in rare_residues):  # If it is a modified residue, than add it to the normal residues...
                    if atom_name[0] in ('H', 'D'):
                        continue
                    if res_nr not in modres_map[chain_number]:
                        modres_map[chain_number][res_nr] = res_name
                    if not coordinate_map[chain_number][0].contains(res_nr):
                        coordinate_map[chain_number][0].add_item(res_nr, (res_name, []))
                    coordinate_map[chain_number][0].get_item(res_nr)[1].append((atom_nr, x, y, z))

                    if b_factors[chain_number] is None:
                        b_factors[chain_number] = {}
                    if res_nr not in b_factors[chain_number]:
                        b_factors[chain_number][res_nr] = []
                    b_factors[chain_number][res_nr].append(b_factor)

                    if all_coords[chain_number] is None:
                        all_coords[chain_number] = []

                    all_coords[chain_number].append((x, z, y))

                    if res_contig_map[chain_number] is None:
                        res_contig_map[chain_number] = Residue_Map()
                        res_contig_map[chain_number].add_item(res_nr, (1, res_name))
                        contig_help_map[chain_number] = 1
                    elif not res_contig_map[chain_number].contains(res_nr):
                        contig_help_map[chain_number] += 1
                        res_contig_map[chain_number].add_item(res_nr, (contig_help_map[chain_number], res_name))

                else:
                    if not coordinate_map[chain_number][1].contains(res_nr):  # If not, then add it to the ligands
                        coordinate_map[chain_number][1].add_item(res_nr, (res_name, []))
                    coordinate_map[chain_number][1].get_item(res_nr)[1].append((atom_nr, x, y, z))

                    if all_coords[chain_number] is None:
                        all_coords[chain_number] = []

                    all_coords[chain_number].append((x, z, y))

                    if res_name in METAL_ATOMS:
                        metals.add((chain_number, res_nr))
                    elif res_name in ION_ATOMS:
                        ions.add((chain_number, res_nr))
                    else:
                        ligands.add((chain_number, res_nr))


    box_map: list = [None] * slots
    for chain_number in chain_order:
        if all_coords[chain_number] is None:
            continue
        # PDB coordinates carry three decimals below 1e4, i.e. seven significant
        # digits, which is exactly what float32 holds. The distance kernels these
        # feed already work in float32.
        all_coords[chain_number] = np.array(all_coords[chain_number], dtype=np.float32)
        box_map[chain_number] = np.empty(6, dtype=np.float32)

        box_map[chain_number][0] = min(all_coords[chain_number][...,0])
        box_map[chain_number][1] = max(all_coords[chain_number][...,0])

        box_map[chain_number][2] = min(all_coords[chain_number][...,1])
        box_map[chain_number][3] = max(all_coords[chain_number][...,1])

        box_map[chain_number][4] = min(all_coords[chain_number][...,2])
        box_map[chain_number][5] = max(all_coords[chain_number][...,2])

    # New Peptide detection counts irregular amino acids
    for chain_number in range(slots):
        if peptide_count[chain_number] is None:
            continue
        if chain_type_map[chain_number] == 'DNA' or chain_type_map[chain_number] == 'RNA':
            continue

        if peptide_count[chain_number][0] <= peptide_count[chain_number][1] * 2:
            chain_type_map[chain_number] = 'Peptide'
        elif (peptide_count[chain_number][0] + peptide_count[chain_number][1]) < 300:  # Total number of atoms
            chain_type_map[chain_number] = 'Peptide'

    ordered_b_factors: list = [None] * slots
    for chain_number in chain_order:
        ordered_keys = coordinate_map[chain_number][0].get_ordered_keys()
        lens = []
        for key in ordered_keys:
            lens.append(len(coordinate_map[chain_number][0].get_item(key)[1]))
        if len(lens) > 0:
            max_atoms = max(lens)

            # float32 holds the coordinates exactly (see the box_map comment
            # above) and the atom number alongside them: the 100k-atom fix in
            # parsePDB caps atom numbers at 7 digits, well inside float32's exact
            # integer range of 2**24.
            residue_matrix: np.array = np.full((len(ordered_keys), max_atoms, 4), np.nan, dtype = np.float32)
            res_contigs = np.empty(len(ordered_keys), dtype='|S3')

            if chain_type_map[chain_number] == 'Protein':
                ordered_b_factors[chain_number] = np.full((len(ordered_keys), max_atoms), np.nan, dtype=np.float32)

                for res_index, key in enumerate(ordered_keys):
                    res_name, atom_quads = coordinate_map[chain_number][0].get_item(key)
                    res_contigs[res_index] = res_name
                    for index, atom_quad in enumerate(atom_quads):
                        try:
                            residue_matrix[res_index, index] = atom_quad
                        except IndexError as err:
                            raise IndexError(f'{err=} {chain_number=} {res_name=} {res_index=} {input_page.split('\n')[:100]=}')
                    ordered_b_factors[chain_number][res_index][:len(b_factors[chain_number][key])] = np.array(b_factors[chain_number][key], dtype=np.float32)

            coordinate_map[chain_number][0] = res_contigs, residue_matrix

            if modres_map[chain_number] is not None:
                modres_map[chain_number] = reindex_by_contig(modres_map[chain_number], res_contig_map[chain_number])

            if ssbond_map[chain_number] is not None:
                ssbond_map[chain_number] = reindex_by_contig(ssbond_map[chain_number], res_contig_map[chain_number])

            if link_map[chain_number] is not None:
                link_map[chain_number] = reindex_by_contig(link_map[chain_number], res_contig_map[chain_number])

            if cis_conformation_map[chain_number] is not None:
                cis_conformation_map[chain_number] = reindex_by_contig(cis_conformation_map[chain_number], res_contig_map[chain_number])

            if cis_follower_map[chain_number] is not None:
                cis_follower_map[chain_number] = reindex_by_contig(cis_follower_map[chain_number], res_contig_map[chain_number])
        else:
            coordinate_map[chain_number][0] = [[],[]]

        ordered_keys = coordinate_map[chain_number][1].get_ordered_keys()
        lens = []
        for key in ordered_keys:
            lens.append(len(coordinate_map[chain_number][1].get_item(key)[1]))
        if len(lens) > 0:
            max_atoms = max(lens)

            residue_matrix: np.array = np.full((len(ordered_keys), max_atoms, 4), np.nan, dtype = np.float32)
            res_contigs = np.empty(len(ordered_keys), dtype='|S3')

            for res_index, key in enumerate(ordered_keys):
                res_name, atom_quads = coordinate_map[chain_number][1].get_item(key)
                res_contigs[res_index] = res_name
                for index, atom_quad in enumerate(atom_quads):
                    residue_matrix[res_index, index] = atom_quad

            coordinate_map[chain_number][1] = res_contigs, residue_matrix
        else:
            coordinate_map[chain_number][1] = [], []

    # --- into the returned structure ------------------------------------------
    # Everything above is already indexed by chain number, so this only moves the
    # columns onto the container and records which chains were present.
    parsed = Parsed_structure()
    parsed.chain_numbers = np.array(chain_order, dtype=np.int64)

    for chain_number in chain_order:
        entry = coordinate_map[chain_number]
        parsed.residue_names[chain_number], parsed.residue_atoms[chain_number] = entry[0]
        parsed.hetero_names[chain_number], parsed.hetero_atoms[chain_number] = entry[1]

        if box_map[chain_number] is not None:
            parsed.box[chain_number] = box_map[chain_number]
        if chain_type_map[chain_number] is not None:
            parsed.chain_type[chain_number] = _PS_CHAIN_TYPE_CODES[chain_type_map[chain_number]]
        parsed.b_factors[chain_number] = ordered_b_factors[chain_number]
        parsed.res_contig[chain_number] = res_contig_map[chain_number]
        parsed.modres[chain_number] = modres_map[chain_number]
        parsed.ssbond[chain_number] = ssbond_map[chain_number]
        parsed.link[chain_number] = link_map[chain_number]
        parsed.cis_conformation[chain_number] = cis_conformation_map[chain_number]
        parsed.cis_follower[chain_number] = cis_follower_map[chain_number]

    parsed.ligands = ligands
    parsed.metals = metals
    parsed.ions = ions

    return parsed

def get_atomic_coverage(parsed: "Parsed_structure", chain_number: int) -> float:
    expected_number_of_atoms = 0
    resolved_number_of_atoms = 0
    for res_nr, res_name in enumerate(parsed.residue_names[chain_number]):
        res_name: str
        atoms = parsed.residue_atoms[chain_number][res_nr]
        try:
            one_letter = BIN_THREE_TO_ONE[res_name]
        except KeyError:
            one_letter = b'X'

        if one_letter == b'X':
            expected_number_of_atoms += 1
            continue

        if one_letter.decode('ascii') not in CORRECT_COUNT:
            expected_number_of_atoms += 1
            continue

        expected_number_of_atoms += CORRECT_COUNT[one_letter.decode('ascii')]
        resolved_number_of_atoms += len(atoms)
    if expected_number_of_atoms == 0:
        return 0
    return (resolved_number_of_atoms / expected_number_of_atoms)

@njit
def exact_top20_search(
        target_res_id: np.int32,
        top20_res: np.ndarray[np.int32],
        top20_dists: np.ndarray[np.float32],
        target_coords: np.ndarray,
        inter_coords: np.ndarray,
        distance_threshold: np.float32
        ):
    
    target_atom_quads: np.ndarray = target_coords[target_res_id]

    inside_sphere_r: list[int] = []
    inside_sphere_d: list[float] = []

    sub_ds: list[float] = []
    ress: list[int] = []
    target_atoms: list[int] = []
    chain_atoms: list[int] = []

    for index in range(len(top20_res)):
        res: int = top20_res[index]
        d: float = top20_dists[index]
        inside_ress: list[int] = []
        inside_dists: list[float] = []
        for atom_quad in inter_coords[res]:
            coord1 = atom_quad[1:4]
            if np.isnan(coord1[0]):
                continue
            for target_atom_quad in target_atom_quads:
                coord2 = target_atom_quad[1:4]
                if np.isnan(coord2[0]):
                    continue
                d = distance(coord1, coord2)
                if d < 1.2:
                    continue
                if d < distance_threshold:
                    inside_ress.append(res)
                    inside_dists.append(d)

                sub_ds.append(d)
                ress.append(res)
                target_atoms.append(target_atom_quad[0])
                chain_atoms.append(atom_quad[0])

        if len(inside_dists) > 0:
            min_index: int = np.argmin(np.array(inside_dists))
            min_res: int = inside_ress[min_index]
            min_d: float = inside_dists[min_index]
            inside_sphere_r.append(min_res)
            inside_sphere_d.append(min_d)

    if len(sub_ds) > 0:
        sub_ds: np.ndarray = np.array(sub_ds)
        min_index: int = np.argmin(sub_ds)
        min_sub_d: float = sub_ds[min_index]
        min_res: int = ress[min_index]
        min_atom_sub: int = target_atoms[min_index]
        min_atom_chain: int = chain_atoms[min_index]
    else:
        min_sub_d: float = np.nan
        min_res: int = -1
        min_atom_sub: int = -1
        min_atom_chain: int = -1

    inside_sphere_r: np.ndarray[np.int32] = np.array(inside_sphere_r, dtype=np.int32)
    inside_sphere_d: np.ndarray[np.float32] = np.array(inside_sphere_d, dtype = np.float32)
    return min_sub_d, min_res, min_atom_sub, min_atom_chain, inside_sphere_r, inside_sphere_d

def is_overlapping(x1: float, x2: float, y1: float, y2: float) -> bool:
    return max(x1,y1) <= min(x2,y2)

def box_check(
        box_1: tuple[int, int, int, int, int, int],
        box_2: tuple[int, int, int, int, int, int],
        distance_threshold: float = 5.0) -> bool:

    [min_x_1, max_x_1, min_y_1, max_y_1, min_z_1, max_z_1] = box_1
    [min_x_2, max_x_2, min_y_2, max_y_2, min_z_2, max_z_2] = box_2

    touching: bool = False
    if is_overlapping(min_x_1 - distance_threshold, max_x_1 + distance_threshold, min_x_2, max_x_2):
        touching = True
    if is_overlapping(min_y_1 - distance_threshold, max_y_1 + distance_threshold, min_y_2, max_y_2):
        touching = True
    if is_overlapping(min_z_1 - distance_threshold, max_z_1 + distance_threshold, min_z_2, max_z_2):
        touching = True

    return touching

@njit
def distance_pre_check(coord_1: np.ndarray, coord_2: np.ndarray, distance: np.float32) -> bool:
    if abs(coord_1[0] - coord_2[0]) > distance:
        return True
    if abs(coord_1[1] - coord_2[1]) > distance:
        return True
    if abs(coord_1[2] - coord_2[2]) > distance:
        return True    
    return False

@njit
def internal_build_fuzzy_dm_loop(
    residue_matrix_chain_1: np.ndarray,
    residue_matrix_chain_2: np.ndarray,
    chain_intern: bool,
    distance_threshold: np.float32,
    atom_list_1: np.ndarray,
    atom_list_2: np.ndarray
    ):
    fuzzy_dm_1 = np.full((len(residue_matrix_chain_1), len(residue_matrix_chain_2)), np.nan, dtype = np.float32)
    fuzzy_dm_2 = np.full((len(residue_matrix_chain_2), len(residue_matrix_chain_1)), np.nan, dtype = np.float32)

    for res_index in range(len(residue_matrix_chain_1)):
        test_coord: np.ndarray = residue_matrix_chain_1[res_index, 0][1:4]
        
        for res_index_2 in range(len(residue_matrix_chain_2)):
            if chain_intern and res_index == res_index_2:
                continue
            test_coord_2: np.ndarray = residue_matrix_chain_2[res_index_2, 0][1:4]

            if distance_pre_check(test_coord, test_coord_2, distance_threshold):
                continue

            d = distance(test_coord, test_coord_2)
            if d > distance_threshold:
                continue
            
            d = getMinDist(atom_list_1[res_index], atom_list_2[res_index_2])

            fuzzy_dm_1[res_index, res_index_2] = d
            if chain_intern:
                fuzzy_dm_1[res_index_2, res_index] = d
            else:
                fuzzy_dm_2[res_index_2, res_index] = d

    return fuzzy_dm_1, fuzzy_dm_2

# Key of the fuzzy distance matrix store: the two chain numbers.
_FUZZY_KEY = numba_types.UniTuple(numba_types.int64, 2)


@njit(cache=True)
def _fuzzy_box_check(box_1, box_2, distance_threshold):
    """box_check, in the kernel. Same three-axis test, same or-semantics."""
    touching = False
    if max(box_1[0] - distance_threshold, box_2[0]) <= min(box_1[1] + distance_threshold, box_2[1]):
        touching = True
    if max(box_1[2] - distance_threshold, box_2[2]) <= min(box_1[3] + distance_threshold, box_2[3]):
        touching = True
    if max(box_1[4] - distance_threshold, box_2[4]) <= min(box_1[5] + distance_threshold, box_2[5]):
        touching = True
    return touching


@njit(cache=True)
def _calc_fuzzy_dm(chain_numbers, residue_atoms, box, chain_type, protein_code,
                   distance_threshold, box_threshold, fuzzy_dm):
    """Every chain pair of a complex, into a (chain, chain) keyed store.

    The pair walk is what used to sit in python around
    internal_build_fuzzy_dm_loop, reading three chain-keyed dicts per pair. The
    chains are columns now, so the whole thing runs here; the dict of dicts becomes
    one store keyed by the two chain numbers.
    """
    n_chains = chain_numbers.shape[0]
    for index_1 in range(n_chains):
        chain = chain_numbers[index_1]
        atom_list_1 = residue_atoms[chain]
        if atom_list_1.shape[0] == 0:
            continue

        # `processed_chains` in the python version: a pair is only built once, from
        # the first chain of the two that the outer loop reaches.
        for index_2 in range(index_1, n_chains):
            chain_2 = chain_numbers[index_2]
            if not (chain_type[chain] == protein_code or chain_type[chain_2] == protein_code):
                continue
            if not _fuzzy_box_check(box[chain], box[chain_2], box_threshold):
                continue

            atom_list_2 = residue_atoms[chain_2]
            if atom_list_2.shape[0] == 0:
                continue

            chain_intern = chain == chain_2
            fuzzy_dm_1, fuzzy_dm_2 = internal_build_fuzzy_dm_loop(
                atom_list_1, atom_list_2, chain_intern, distance_threshold,
                atom_list_1, atom_list_2)

            fuzzy_dm[(chain, chain_2)] = fuzzy_dm_1
            if not chain_intern:
                fuzzy_dm[(chain_2, chain)] = fuzzy_dm_2


def calcFuzzyDM(
        parsed: "Parsed_structure",
        config: Config,
        distance_threshold: float):
    """Residue-level distance matrices for every touching pair of chains.

    Returns a numba typed dict keyed by (chain number, chain number); use
    fuzzy_row(store, a, b) to read one out.
    """
    residue_atoms = NumbaList()
    for number in range(MAX_NUMBER_OF_CHAINS):
        atoms = parsed.residue_atoms[number]
        if atoms is None or len(atoms) == 0:
            residue_atoms.append(_EMPTY_ATOMS)
        else:
            residue_atoms.append(np.ascontiguousarray(atoms, dtype=np.float32))

    fuzzy_dm = NumbaDict.empty(_FUZZY_KEY, numba_types.float32[:, ::1])
    _calc_fuzzy_dm(parsed.chain_numbers, residue_atoms, parsed.box, parsed.chain_type,
                   np.int8(_PS_PROTEIN), np.float32(distance_threshold),
                   np.float32(config.short_distance_threshold), fuzzy_dm)
    return fuzzy_dm


def fuzzy_row(fuzzy_dm, chain_number: int, other_chain_number: int):
    """The (chain, other) matrix, or None when the pair was never built."""
    key = (chain_number, other_chain_number)
    if key in fuzzy_dm:
        return fuzzy_dm[key]
    return None

@njit
def calc_min_dist_between_atomic_coordinates_w(
    atomic_coordinates_res_1: np.ndarray,
    atomic_coordinates_res_2: np.ndarray
    ) -> tuple[float, int, int]:

    min_d: float = float('inf')
    for p1, coord_1 in enumerate(atomic_coordinates_res_1):
        if np.isnan(coord_1[0]):
            continue
        for p2, coord_2 in enumerate(atomic_coordinates_res_2):
            if np.isnan(coord_2[0]):
                continue
            d: float = distance(coord_1, coord_2)
            if d < min_d:
                min_d = d
                min_atom_p_1: int = p1
                min_atom_p_2: int = p2

    return min_d, min_atom_p_1, min_atom_p_2

@njit
def calc_min_dist_between_atomic_coordinates(
    atomic_coordinates_res_1: np.ndarray,
    atomic_coordinates_res_2: np.ndarray
    ) -> float:

    min_d: float = float('inf')
    for coord_1 in atomic_coordinates_res_1:
        if np.isnan(coord_1[0]):
            continue
        for coord_2 in atomic_coordinates_res_2:
            if np.isnan(coord_2[0]):
                continue
            d: float = distance(coord_1, coord_2)
            if d < min_d:
                min_d = d

    return min_d

@njit
def getMinDist(
        atom_list_res_1: np.ndarray,
        atom_list_res_2: np.ndarray,
        ) -> float:
    
    atomic_coordinates_res_1: np.ndarray = atom_list_res_1[...,1:4]
    atomic_coordinates_res_2: np.ndarray = atom_list_res_2[...,1:4]

    min_d = calc_min_dist_between_atomic_coordinates(atomic_coordinates_res_1, atomic_coordinates_res_2)
    return min_d

@njit
def getMinDist_w_atom_names(
        atom_list_res_1: np.ndarray,
        atom_list_res_2: np.ndarray,
        
        ) -> tuple[float, int, int]:
    
    atom_nrs: np.ndarray = atom_list_res_1[...,0]
    atomic_coordinates_res_1: np.ndarray = atom_list_res_1[...,1:4]

    atom_nrs_2: np.ndarray = atom_list_res_2[...,0]
    atomic_coordinates_res_2: np.ndarray = atom_list_res_2[...,1:4]

    min_d, min_atom_p_1, min_atom_p_2 = calc_min_dist_between_atomic_coordinates_w(atomic_coordinates_res_1, atomic_coordinates_res_2)

    min_atom: int = int(atom_nrs[min_atom_p_1])
    min_atom2: int = int(atom_nrs_2[min_atom_p_2])

    return min_d, min_atom, min_atom2

def calculate_interfaces(
        IAmap: list[tuple[str, int, str, int, str, str , float]],
        parsed: "Parsed_structure"
        ) -> dict[str, list[tuple[str, np.ndarray[float]]]]:

    # IAmap comes out of the RIN files keyed by chain id, and the result goes on to
    # proteins.complexes keyed the same way, so the chain ids stay. Only the type
    # test moves onto the parsed columns: the two conditions the dict served were
    # "has a type at all" and "is a Protein", which are one set membership each.
    protein_chains = {CHAIN_ID_ORDER[number] for number in parsed.protein_chain_numbers()}

    filtered_list = [(x[0], x[1], x[2], x[3], x[6]) for x in IAmap
                     if x[4] == 'combi' and x[5] == 'all'
                     and x[0] in protein_chains and x[2] in protein_chains]

    filtered_list = [x for x in filtered_list if x[0] != x[2]]

    interfaces: dict[str, dict[str, tuple[int, int, float]]] = {}

    for chain, res, chain_b, res_b, score in filtered_list:
        if chain not in interfaces:
            interfaces[chain] = {}
        if chain_b not in interfaces[chain]:
            interfaces[chain][chain_b] = []
        
        interfaces[chain][chain_b].append((res, res_b, score))

    ret_interfaces: dict[str, list[tuple[str, np.ndarray[np.int64], np.ndarray[np.float32]]]] = {}
    for chain in interfaces:
        ret_interfaces[chain] = []
        for chain_b in interfaces[chain]:
            triple_list = interfaces[chain][chain_b]

            rows = np.array([x[0] for x in triple_list], dtype=np.int32)
            cols = np.array([x[1] for x in triple_list], dtype=np.int32)
            # Interaction scores, so float32 is ample. These are both stored and
            # held concatenated in the ray object store during aggregation.
            vals = np.array([x[2] for x in triple_list], dtype=np.float32)

            # 1. Encode coordinates into a single 64-bit integer
            # We shift the row index by 32 bits to make space for the column index
            combined_keys = (rows.astype(np.int64) << 32) | cols.astype(np.int64)

            # 2. Sort everything by these keys once (before putting in Ray)
            sort_idx = np.argsort(combined_keys)
            sorted_keys = combined_keys[sort_idx]
            sorted_vals = vals[sort_idx]

            ret_interfaces[chain].append((chain_b, sorted_keys, sorted_vals))

    return ret_interfaces

@ray.remote(max_calls=1)
def annotate_wrapper(conf_dump, chunk: list[str]):
    config: Config = conf_dump
    reset_logger_for_remotes(config)
    return annotate(config, chunk, locked = True)

def annotate(config: Config, chunk: list[str], locked = False):

    total_times = []
    t0 = time.time()

    structural_analysis_list: list[tuple[str, dict[str, list[np.ndarray[float]]]]] = []
    interface_map: dict[str, dict[str, dict[str, interface_package.Interface]]] = {}

    struct_bins = [[] for _ in range(20)]
    for structure_id in chunk:
        bin_nr = struct_id_to_bin_nr(structure_id)
        struct_bins[bin_nr].append(structure_id)

    structure_ids: dict[str, dict[str, int]] = {}

    cols = ['PDB', 'Structure_Id', 'Chain', 'Complex', 'Resolution', 'Homooligomer']

    compl_id_pdb_id_dict = {}
    homomer_dict: dict[str, set[str]] = {}
    for db_nr, binned_structure_ids in enumerate(struct_bins):
        if len(binned_structure_ids) > 0:
            results = retrieve_structure_info_from_struct_ids(config, binned_structure_ids, db_nr, cols = cols)
            for row in results:
                struct_id, struct_db_id, chain, compl_db_id, resolution, homomers = row

                homomer_set = set(homomers)
                homomer_dict[struct_id] = homomer_set

                if resolution is None: #Emulate a resolution of 2.0A for AF structure and other models
                    resolution = 2.0

                if resolution > config.option_res_thresh:
                    continue

                if struct_id not in structure_ids:
                    structure_ids[struct_id] = {}
                    compl_id_pdb_id_dict[struct_id] = compl_db_id
                internal_struct_db_id = 100*struct_db_id + db_nr
                structure_ids[struct_id][chain] = internal_struct_db_id

    t0 = add_to_times(total_times, t0, label='retrieve structure info')

    comp_times = []

    total_loop_times = []
    for structure_id in chunk:
        if structure_id not in structure_ids:
            continue

        loop_times = []
        tl = time.time()
        tl0 = time.time()

        model_path = None
        # set model path to pdb file path in custom db path 
        if config.custom_db_path:
            custom_path_list = os.listdir(config.custom_db_path)
            for model in custom_path_list:
                if model.startswith(structure_id):
                    model_path = config.custom_db_path + "/%s" %model

        elif config.model_db_active:
            if is_model(structure_id):
                if is_alphafold_model(structure_id):
                    model_path = None
                else:
                    model_path = complex_id_to_file_path(structure_id, config)
        interfaces: dict[str, dict[str, interface_package.Interface]]
        structural_analysis_dict: dict[str, list[np.ndarray[float]]]
        (structural_analysis_dict, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles,
         interfaces, times) = structuralAnalysis(structure_id, config, homomer_dict[structure_id], model_path = model_path)

        interface_map[structure_id] = interfaces

        loop_times.append(times)
        tl = add_to_times(loop_times, tl, 'structuralAnalysis')

        # The four profile dicts used to be packed here and written to Complex.Profiles.
        # Nothing ever selected that column, so the pack and the UPDATE it fed are both
        # gone. The column is left in the schema so existing databases are unaffected.

        structural_analysis_list.append((structure_id, structural_analysis_dict))

        tl = add_to_times(loop_times, tl , 'process structural analysis results')
    
        total_loop_times = aggregate_times(total_loop_times, loop_times)

        tl1 = time.time()

        comp_times.append((structure_id,tl1-tl0))

        if config.verbosity >= 6:
            config.logger.info(f'Time for analyzing {structure_id} {tl1-tl0}')

    total_times.append(total_loop_times)

    t0 = add_to_times(total_times, t0, label='chain analysis loop')

    remote_insertResidues(structural_analysis_list, structure_ids, interface_map, config, locked=locked)

    t0 = add_to_times(total_times, t0, label='insertResidues')

    if len(comp_times) > 0:
        max_comp_time = sorted(comp_times, key=lambda x:x[1], reverse=True)[0]
    else:
        max_comp_time = (None, 0)

    return total_times, max_comp_time

def _degree_column(profile_matrix: np.ndarray, degree_index: int) -> np.ndarray:
    # Same semantics as Interaction_type.getDegree over a whole chain at once:
    # set_data_from_list truncates the degree with int(), and a missing (NaN) slot counts as 0.
    column = profile_matrix[:, degree_index]
    return np.where(np.isnan(column), 0.0, np.trunc(column))


def build_chain_profile_lookups(
        profiles: dict[str, list[Interaction_profile]]
        ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    # Everything the sphere loops in analysis_chain need from a neighbouring residue's
    # interaction profile is a per-residue constant, but the loops used to re-derive it
    # for every target residue whose sphere that residue falls into. A residue sits in
    # tens of such spheres, and each lookup walked the four-level Interaction_profile
    # dispatch. Deriving these once per complex turns the inner loops into array reads.
    ppi_cnt_index = profile_degree_index('sc', 'cnt', 'Protein')
    ppi_hbond_index = profile_degree_index('sc', 'hbond', 'Protein')
    hom_cnt_index = profile_degree_index('sc', 'cnt', 'Homomer')
    hom_hbond_index = profile_degree_index('sc', 'hbond', 'Homomer')

    lookups: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}

    for chain, chain_profiles in profiles.items():
        profile_matrix, valid = profile_matrix_from_profiles(chain_profiles)

        has_ppi = (_degree_column(profile_matrix, ppi_cnt_index)
                   + _degree_column(profile_matrix, ppi_hbond_index)) > 0
        has_hom = (_degree_column(profile_matrix, hom_cnt_index)
                   + _degree_column(profile_matrix, hom_hbond_index)) > 0

        interface_scores = total_interface_scores_from_matrix(profile_matrix, valid)

        # `valid` is carried out as well: the sphere loops test
        # `chain_profiles[res] is not None`, which is exactly this mask, and a kernel
        # cannot look at a list of objects.
        lookups[chain] = (has_ppi & valid, has_hom & valid, interface_scores, valid,
                          profile_matrix)

    return lookups



# --- the per-chain distance loop ---------------------------------------------
#
# Layout of the kernel's two result buffers. Distances are NaN where the python
# version left a None, and the degrees are -1 there.

_R_LIG_DIST = 0
_R_METAL_DIST = 1
_R_ION_DIST = 2
_R_CHAIN_DIST = 3
_R_RNA_DIST = 4
_R_DNA_DIST = 5
_R_INTER_MEDIAN_KD = 6
_R_INTER_DW_KD = 7
_R_INTER_MEDIAN_RSA = 8
_R_INTER_DW_RSA = 9
_R_INTER_INTERACTIONS_MEDIAN = 10
_R_INTER_INTERACTIONS_DW = 11
_R_INTER_WEIGHTS = 12
_R_INTRA_MEDIAN_KD = 13
_R_INTRA_DW_KD = 14
_R_INTRA_MEDIAN_RSA = 15
_R_INTRA_DW_RSA = 16
_R_INTRA_WEIGHTS = 17
_R_INTRA_INTERACTIONS_MEDIAN = 18
_R_INTRA_INTERACTIONS_DW = 19
_R_INTRA_INTERACTIONS_WEIGHTS = 20
_R_LEN = 21

_I_INTER_PPI = 0
_I_INTER_HOM = 1
_I_INTRA_PPI = 2
_I_INTRA_HOM = 3
_I_BROKEN_CHAIN = 4
_I_LEN = 5

_CHAIN_LOOP_OK = 0
_CHAIN_LOOP_NO_INTER_PROFILES = 1
_CHAIN_LOOP_NO_INTRA_PROFILES = 2


@njit(cache=True)
def _sum_float32(values):
    """sum() over a list of np.float32, which is a naive left-to-right float32 add.

    CPython's compensated summation only kicks in for exact python floats, and a
    float32 is not one, so the fallback `0 + x + y + ...` is what ran - in float32,
    because a python int and a python float are both weak next to it.
    """
    total = np.float32(0.0)
    for index in range(len(values)):
        total += values[index]
    return total


@njit(cache=True)
def _sum_float64(values):
    """Same, for the float64 lists: interface scores are float64, so score / dist is."""
    total = np.float64(0.0)
    for index in range(len(values)):
        total += values[index]
    return total


@njit(cache=True)
def _median(values):
    """base_utils.median over a numeric sequence, same halves and same rounding."""
    n = len(values)
    if n == 0:
        return np.nan
    ordered = np.sort(np.array(values))
    if n == 1:
        return ordered[0]
    if n % 2 == 0:
        return (ordered[(n // 2) - 1] + ordered[n // 2]) / 2.0
    return ordered[(n - 1) // 2]


@njit(cache=True)
def _top20(n_of_res, chain_intern, target_res_id, dist_matrix, top):
    """get_top20 in the kernel; the caller used to run this loop in python."""
    unsorted_res = []
    unsorted_dists = []
    for res in range(n_of_res):
        if chain_intern and res == target_res_id:
            continue
        d = dist_matrix[target_res_id, res]
        if np.isnan(d):
            continue
        unsorted_res.append(res)
        unsorted_dists.append(d)

    res_array = np.array(unsorted_res, dtype=np.int32)
    dist_array = np.array(unsorted_dists, dtype=np.float32)
    order = np.argsort(dist_array)
    return res_array[order][:top], dist_array[order][:top]


@njit(cache=True)
def _min_sub_dist(target_atoms, chain_atoms, fuzzy_row, has_fuzzy, target_res_id,
                  chain_intern, distance_threshold):
    """getMinSubDist, with its get_top20 and exact_top20_search folded in."""
    if has_fuzzy == 0:
        empty_r = np.empty(0, dtype=np.int32)
        empty_d = np.empty(0, dtype=np.float32)
        return np.nan, -1, -1, -1, empty_r, empty_d

    top20_res, top20_dists = _top20(len(chain_atoms), chain_intern, target_res_id,
                                    fuzzy_row, 20)
    if len(top20_res) == 0:
        empty_r = np.empty(0, dtype=np.int32)
        empty_d = np.empty(0, dtype=np.float32)
        return np.nan, -1, -1, -1, empty_r, empty_d

    return exact_top20_search(target_res_id, top20_res, top20_dists,
                              target_atoms, chain_atoms, distance_threshold)


@njit(cache=True)
def _chain_distance_loop(
        target_res_id, target_atoms, target_number,
        chain_order, branch, is_protein, has_dssp, has_profiles, has_fuzzy, type_code,
        atoms, fuzzy, kd, rsa, rsa_present, ppi, hom, interface_scores, profile_valid,
        het_atoms, het_kind,
        milieu_threshold, intra_milieu_threshold, ligand_sphere,
        result, result_int):
    """Every chain of the complex, for one target residue.

    Replaces the python loop that walked the chain records: the two sphere passes,
    the neighbour accumulations and the hetero-residue distances. The ligand and
    chain distances are reduced to their per-category minima here rather than
    collected into dicts, because that is all get_shortest_distances did with them.

    The accumulation widths follow the python ones exactly: a distance comes out of
    exact_top20_search as float32, so `kd / dist` and the running weight totals are
    float32 sums, while `interface_score / dist` is float64 because the score is.
    """
    for index in range(result.shape[0]):
        result[index] = np.nan
    for index in range(result_int.shape[0]):
        result_int[index] = -1

    inter_kds = [np.float64(x) for x in range(0)]
    inter_dw_kds = [np.float32(x) for x in range(0)]
    inter_rsas = [np.float64(x) for x in range(0)]
    inter_dw_rsas = [np.float32(x) for x in range(0)]
    inter_interactions = [np.float64(x) for x in range(0)]
    inter_interactions_dw = [np.float64(x) for x in range(0)]
    inter_weights = np.float32(0.0)
    inter_ppi_degree = -1
    inter_hom_degree = -1

    min_ligand = np.nan
    min_metal = np.nan
    min_ion = np.nan
    min_chain = np.nan
    min_rna = np.nan
    min_dna = np.nan

    for order_index in range(chain_order.shape[0]):
        chain = chain_order[order_index]
        chain_branch = branch[chain]
        chain_atoms = atoms[chain]

        if chain_branch == 1:
            min_chain_dist, _min_res, _atom_sub, _atom_chain, sphere_r, sphere_d = _min_sub_dist(
                target_atoms, chain_atoms, fuzzy[chain], has_fuzzy[chain],
                target_res_id, False, milieu_threshold)

            # get_shortest_distances bucketed these by the chain's type, counting a
            # chain it had no type for as Protein.
            if not np.isnan(min_chain_dist):
                chain_type = type_code[chain]
                if chain_type == _CHAIN_TYPE_PROTEIN:
                    if np.isnan(min_chain) or min_chain_dist < min_chain:
                        min_chain = min_chain_dist
                elif chain_type == _CHAIN_TYPE_RNA:
                    if np.isnan(min_rna) or min_chain_dist < min_rna:
                        min_rna = min_chain_dist
                elif chain_type == _CHAIN_TYPE_DNA:
                    if np.isnan(min_dna) or min_chain_dist < min_dna:
                        min_dna = min_chain_dist

            if is_protein[chain] == 1 and has_dssp[chain] == 1 and len(sphere_r) > 0:
                chain_kd = kd[chain]
                chain_rsa = rsa[chain]
                chain_rsa_present = rsa_present[chain]
                chain_valid = profile_valid[chain]
                chain_ppi = ppi[chain]
                chain_hom = hom[chain]
                chain_scores = interface_scores[chain]
                for sphere_index in range(len(sphere_r)):
                    neighbour = sphere_r[sphere_index]
                    dist = sphere_d[sphere_index]
                    inter_weights += np.float32(1.0) / dist

                    inter_kds.append(chain_kd[neighbour])
                    inter_dw_kds.append(np.float32(chain_kd[neighbour]) / dist)
                    if chain_rsa_present[neighbour] == 1:
                        inter_rsas.append(chain_rsa[neighbour])
                        inter_dw_rsas.append(np.float32(chain_rsa[neighbour]) / dist)

                    if has_profiles[chain] == 0:
                        result_int[_I_BROKEN_CHAIN] = chain
                        return _CHAIN_LOOP_NO_INTER_PROFILES

                    if chain_valid[neighbour] == 1:
                        score = chain_scores[neighbour]
                        inter_interactions.append(score)
                        inter_interactions_dw.append(score / np.float64(dist))

                        if chain_ppi[neighbour] == 1:
                            if inter_ppi_degree == -1:
                                inter_ppi_degree = 1
                            else:
                                inter_ppi_degree += 1
                        if chain_hom[neighbour] == 1:
                            if inter_hom_degree == -1:
                                inter_hom_degree = 1
                            else:
                                inter_hom_degree += 1

        elif chain_branch == 2:
            if is_protein[chain] == 1:
                if has_dssp[chain] == 0:
                    continue

                _min_d, _min_res, _atom_sub, _atom_chain, sphere_r, sphere_d = _min_sub_dist(
                    target_atoms, chain_atoms, fuzzy[chain], has_fuzzy[chain],
                    target_res_id, True, milieu_threshold)

                intra_kds = [np.float64(x) for x in range(0)]
                intra_dw_kds = [np.float32(x) for x in range(0)]
                intra_rsas = [np.float64(x) for x in range(0)]
                intra_dw_rsas = [np.float32(x) for x in range(0)]
                intra_ppi_degree = -1
                intra_hom_degree = -1
                intra_weights = np.float32(0.0)

                chain_kd = kd[chain]
                chain_rsa = rsa[chain]
                chain_rsa_present = rsa_present[chain]
                chain_valid = profile_valid[chain]
                chain_ppi = ppi[chain]
                chain_hom = hom[chain]
                chain_scores = interface_scores[chain]

                for sphere_index in range(len(sphere_r)):
                    neighbour = sphere_r[sphere_index]
                    dist = sphere_d[sphere_index]
                    intra_weights += np.float32(1.0) / dist

                    intra_kds.append(chain_kd[neighbour])
                    intra_dw_kds.append(np.float32(chain_kd[neighbour]) / dist)
                    if chain_rsa_present[neighbour] == 1:
                        intra_rsas.append(chain_rsa[neighbour])
                        intra_dw_rsas.append(np.float32(chain_rsa[neighbour]) / dist)

                    # An absent profile chain used to surface as a TypeError rather
                    # than the KeyError the old handler caught, so it is simply a
                    # missing profile here.
                    if has_profiles[chain] == 1 and chain_valid[neighbour] == 1:
                        if chain_ppi[neighbour] == 1:
                            if intra_ppi_degree == -1:
                                intra_ppi_degree = 1
                            else:
                                intra_ppi_degree += 1
                        if chain_hom[neighbour] == 1:
                            if intra_hom_degree == -1:
                                intra_hom_degree = 1
                            else:
                                intra_hom_degree += 1

                result_int[_I_INTRA_PPI] = intra_ppi_degree
                result_int[_I_INTRA_HOM] = intra_hom_degree

                if intra_weights > 0:
                    result[_R_INTRA_MEDIAN_KD] = _median(intra_kds)
                    result[_R_INTRA_DW_KD] = _sum_float32(intra_dw_kds) / intra_weights
                    result[_R_INTRA_MEDIAN_RSA] = _median(intra_rsas)
                    result[_R_INTRA_DW_RSA] = _sum_float32(intra_dw_rsas) / intra_weights
                    result[_R_INTRA_WEIGHTS] = intra_weights

                _min_d, _min_res, _atom_sub, _atom_chain, sphere_r, sphere_d = _min_sub_dist(
                    target_atoms, chain_atoms, fuzzy[chain], has_fuzzy[chain],
                    target_res_id, True, intra_milieu_threshold)

                intra_interactions = [np.float64(x) for x in range(0)]
                intra_interactions_dw = [np.float64(x) for x in range(0)]
                intra_weights = np.float32(0.0)
                if len(sphere_r) > 0:
                    if has_profiles[chain] == 0:
                        result_int[_I_BROKEN_CHAIN] = chain
                        return _CHAIN_LOOP_NO_INTRA_PROFILES

                    for sphere_index in range(len(sphere_r)):
                        neighbour = sphere_r[sphere_index]
                        dist = sphere_d[sphere_index]
                        intra_weights += np.float32(1.0) / dist

                        if chain_valid[neighbour] == 1:
                            score = chain_scores[neighbour]
                            intra_interactions.append(score)
                            intra_interactions_dw.append(score / np.float64(dist))

                if intra_weights > 0:
                    result[_R_INTRA_INTERACTIONS_MEDIAN] = _median(intra_interactions)
                    result[_R_INTRA_INTERACTIONS_DW] = (
                        _sum_float64(intra_interactions_dw) / np.float64(intra_weights))
                    result[_R_INTRA_INTERACTIONS_WEIGHTS] = intra_weights

        # --- hetero residues, for every chain whatever branch it took
        chain_het_atoms = het_atoms[chain]
        chain_het_kind = het_kind[chain]
        for het_index in range(chain_het_kind.shape[0]):
            min_d, _atom, _atom2 = getMinDist_w_atom_names(
                target_atoms[target_res_id], chain_het_atoms[het_index])
            if min_d > ligand_sphere or np.isnan(min_d):
                continue
            kind = chain_het_kind[het_index]
            if kind == _HET_KIND_METAL:
                if np.isnan(min_metal) or min_d < min_metal:
                    min_metal = min_d
            elif kind == _HET_KIND_ION:
                if np.isnan(min_ion) or min_d < min_ion:
                    min_ion = min_d
            else:
                if np.isnan(min_ligand) or min_d < min_ligand:
                    min_ligand = min_d

    result[_R_LIG_DIST] = min_ligand
    result[_R_METAL_DIST] = min_metal
    result[_R_ION_DIST] = min_ion
    result[_R_CHAIN_DIST] = min_chain
    result[_R_RNA_DIST] = min_rna
    result[_R_DNA_DIST] = min_dna

    result_int[_I_INTER_PPI] = inter_ppi_degree
    result_int[_I_INTER_HOM] = inter_hom_degree

    if inter_weights > 0:
        result[_R_INTER_MEDIAN_KD] = _median(inter_kds)
        result[_R_INTER_DW_KD] = _sum_float32(inter_dw_kds) / inter_weights
        result[_R_INTER_MEDIAN_RSA] = _median(inter_rsas)
        result[_R_INTER_DW_RSA] = _sum_float32(inter_dw_rsas) / inter_weights
        result[_R_INTER_INTERACTIONS_MEDIAN] = _median(inter_interactions)
        result[_R_INTER_INTERACTIONS_DW] = (
            _sum_float64(inter_interactions_dw) / np.float64(inter_weights))
        result[_R_INTER_WEIGHTS] = inter_weights

    return _CHAIN_LOOP_OK


# --- residue columns and the residue loop ------------------------------------
#
# The per-residue loop used to build a Residue object, set forty-odd attributes on
# it and call encode_to_numpy, which walked the Interaction_profile object tree to
# rebuild the same 132 floats every time. Everything it wrote is a fixed offset in
# the encoded row, so the loop fills an (n_residues x 195) matrix directly and the
# object never exists. The offsets come from Residue.values_per_slot, so they
# cannot drift from the class.

def _encoded_slot_offsets() -> dict[str, int]:
    offsets = {}
    position = 0
    for slot, width in zip(residue_package.Residue.data_slots,
                           residue_package.Residue.values_per_slot):
        offsets[slot] = position
        position += width
    return offsets


_SLOT_OFFSET: dict[str, int] = _encoded_slot_offsets()
_ENCODED_WIDTH: int = sum(residue_package.Residue.values_per_slot)

_O_RSA = _SLOT_OFFSET['RSA']
_O_MC_ACC = _SLOT_OFFSET['relative_main_chain_acc']
_O_SC_ACC = _SLOT_OFFSET['relative_side_chain_acc']
_O_SSA = _SLOT_OFFSET['SSA']
_O_PROFILE = _SLOT_OFFSET['interaction_profile']
_O_CENTRALITIES = _SLOT_OFFSET['centralities']
_O_INFLOW = _SLOT_OFFSET['inflow_map']
_O_PHI = _SLOT_OFFSET['phi']
_O_INTRA_SSBOND = _SLOT_OFFSET['intra_ssbond']
_O_INTER_SSBOND = _SLOT_OFFSET['inter_ssbond']
_O_SSBOND_LENGTH = _SLOT_OFFSET['ssbond_length']
_O_INTRA_LINK = _SLOT_OFFSET['intra_link']
_O_INTER_LINK = _SLOT_OFFSET['inter_link']
_O_LINK_LENGTH = _SLOT_OFFSET['link_length']
_O_CIS_CONFORMATION = _SLOT_OFFSET['cis_conformation']
_O_CIS_FOLLOWER = _SLOT_OFFSET['cis_follower']
_O_INTER_MEDIAN_KD = _SLOT_OFFSET['inter_chain_median_kd']
_O_INTRA_MEDIAN_KD = _SLOT_OFFSET['intra_chain_median_kd']
_O_INTER_INTERACTIONS_MEDIAN = _SLOT_OFFSET['inter_chain_interactions_median']
_O_INTRA_INTERACTIONS_MEDIAN = _SLOT_OFFSET['intra_chain_interactions_median']
_O_B_FACTOR = _SLOT_OFFSET['b_factor']
_O_PLDDT = _SLOT_OFFSET['pLDDT']
_O_MODRES = _SLOT_OFFSET['modres']
_O_AA = _SLOT_OFFSET['aa']
_O_LIG_DIST = _SLOT_OFFSET['lig_dist']
_O_INTRA_PPI_DEGREE = _SLOT_OFFSET['intra_chain_ppi_degree']
_O_INTRA_HOM_DEGREE = _SLOT_OFFSET['intra_chain_hom_degree']
_O_INTER_PPI_DEGREE = _SLOT_OFFSET['inter_chain_ppi_degree']
_O_INTER_HOM_DEGREE = _SLOT_OFFSET['inter_chain_hom_degree']

_N_CENTRALITY_VALUES = residue_package.Residue.values_per_slot[
    residue_package.Residue.data_slots.index('centralities')]
_N_INFLOW_VALUES = residue_package.Residue.values_per_slot[
    residue_package.Residue.data_slots.index('inflow_map')]
_N_PROFILE_VALUES = residue_package.Residue.values_per_slot[
    residue_package.Residue.data_slots.index('interaction_profile')]

_AA_ENCODES = residue_package.string_encodes[
    residue_package.Residue.data_slots.index('aa')]
_SSA_ENCODES = residue_package.string_encodes[
    residue_package.Residue.data_slots.index('SSA')]
_MODRES_ENCODES = residue_package.string_encodes[
    residue_package.Residue.data_slots.index('modres')]


def _aa_code_column(residue_names: np.ndarray) -> np.ndarray:
    """The encoded one-letter code per residue."""
    codes = np.empty(len(residue_names), dtype=np.float64)
    for index in range(len(residue_names)):
        three_letter = residue_names[index]
        one_letter = 'X'
        if three_letter in BIN_THREE_TO_ONE:
            decoded = BIN_THREE_TO_ONE[three_letter].decode('ascii')[0]
            if decoded in ONE_TO_THREE:
                one_letter = decoded
        codes[index] = _AA_ENCODES[one_letter]
    return codes


def _dssp_columns(target_dssp, has_target_dssp, dssp_enabled, n_residues):
    """The dssp rows as columns, plus how many leading residues actually get them.

    The python loop turned dssp off for the rest of the chain the first time a
    residue had no dssp row, so the rows apply to a prefix; that boundary is found
    once here instead of being rediscovered per residue.
    """
    values = np.full((n_residues, 9), np.nan, dtype=np.float64)
    if not dssp_enabled:
        return values, 0
    if not has_target_dssp:
        return values, 0

    limit = n_residues
    for index in range(n_residues):
        row = target_dssp[index]
        if row is None:
            limit = index
            break
        values[index, 0] = row[0]
        values[index, 1] = row[1]
        values[index, 2] = row[2]
        try:
            values[index, 3] = _SSA_ENCODES[row[3]]
        except KeyError as err:
            raise IndexError(f'slot_index=3 SSA entry={row[3]!r} {err=}')
        for offset in range(4, 9):
            values[index, offset] = row[offset]
    return values, limit


def _object_matrix(entries, width, n_residues):
    """A list of to_list()-able objects (or None) as a matrix plus a validity mask."""
    matrix = np.full((n_residues, width), np.nan, dtype=np.float64)
    valid = np.zeros(n_residues, dtype=np.uint8)
    if entries is None:
        return matrix, valid
    for index in range(n_residues):
        entry = entries[index]
        if entry is None:
            continue
        matrix[index] = entry.to_list()
        valid[index] = 1
    return matrix, valid


def _build_residue_columns(
        target_residue_names, target_dssp, has_target_dssp, dssp_enabled,
        target_profile_matrix, target_profile_valid, target_profiles, has_target_profiles,
        target_centralities, has_target_centralities,
        target_inflow, has_target_inflow,
        target_b_factors, target_modres, has_target_modres,
        ssbond_map, link_map, cis_conformation_map, cis_follower_map, target_chain):
    """Everything the residue loop writes that does not depend on the chain loop."""
    n_residues = len(target_residue_names)

    aa_code = _aa_code_column(target_residue_names)
    dssp_values, dssp_limit = _dssp_columns(
        target_dssp, has_target_dssp, dssp_enabled, n_residues)

    if has_target_profiles and target_profile_matrix is not None:
        profile_matrix = np.ascontiguousarray(target_profile_matrix, dtype=np.float64)
        profile_valid = np.ascontiguousarray(target_profile_valid).astype(np.uint8)
    elif has_target_profiles:
        profile_matrix, profile_valid = _object_matrix(
            target_profiles, _N_PROFILE_VALUES, n_residues)
    else:
        profile_matrix = np.full((n_residues, _N_PROFILE_VALUES), np.nan, dtype=np.float64)
        profile_valid = np.zeros(n_residues, dtype=np.uint8)

    centrality_matrix, centrality_valid = _object_matrix(
        target_centralities if has_target_centralities else None,
        _N_CENTRALITY_VALUES, n_residues)

    # set_data_from_list followed by to_list is an identity round trip over the
    # eight slots, so the stored list is the encoded row and no object is built.
    inflow_matrix = np.full((n_residues, _N_INFLOW_VALUES), np.nan, dtype=np.float64)
    if has_target_inflow:
        for index in range(n_residues):
            inflow_matrix[index] = target_inflow[index]

    # nansum over the atoms of a residue, divided by how many were not NaN. The
    # stored b factors are float32 and the divisor was a python float, so both the
    # sum and the division happened in float32; widening them here would shift the
    # last bit.
    b_factor_matrix = np.asarray(target_b_factors)
    counts = (~np.isnan(b_factor_matrix)).sum(axis=1).astype(b_factor_matrix.dtype)
    with np.errstate(invalid='ignore', divide='ignore'):
        avg_b_factor = np.nansum(b_factor_matrix, axis=1) / counts

    modres_code = np.full(n_residues, np.nan, dtype=np.float64)
    if has_target_modres:
        for index in range(n_residues):
            if index not in target_modres:
                continue
            try:
                modres_code[index] = _MODRES_ENCODES[target_modres[index]]
            except KeyError:
                modres_code[index] = _MODRES_ENCODES['UNK']

    ssbond_length = np.full(n_residues, np.nan, dtype=np.float64)
    ssbond_intra = np.zeros(n_residues, dtype=np.uint8)
    link_intra = np.zeros(n_residues, dtype=np.uint8)
    link_present = np.zeros(n_residues, dtype=np.uint8)
    link_length = np.full(n_residues, np.nan, dtype=np.float64)
    cis_conformation = np.full(n_residues, np.nan, dtype=np.float64)
    cis_follower = np.full(n_residues, np.nan, dtype=np.float64)
    # The caller passes this chain's own maps, keyed by residue. The loop used to
    # ask a chain-keyed map for the tuple key (target_chain, target_res_id), which
    # it can never hold - so ssbond, link and both cis slots were never written.
    chain_ssbond = ssbond_map or {}
    chain_link = link_map or {}
    chain_cis_conformation = cis_conformation_map or {}
    chain_cis_follower = cis_follower_map or {}
    for index in range(n_residues):
        if index in chain_ssbond:
            chain_2, _res_nr_2, length = chain_ssbond[index]
            ssbond_length[index] = length
            ssbond_intra[index] = 1 if target_chain == chain_2 else 0
        if index in chain_link:
            link_present[index] = 1
            link_intra[index] = 1 if target_chain == chain_link[index][5] else 0
            link_length[index] = chain_link[index][6]
        if index in chain_cis_conformation:
            cis_conformation[index] = chain_cis_conformation[index][3]
        if index in chain_cis_follower:
            cis_follower[index] = chain_cis_follower[index][3]

    ssbond_present = ~np.isnan(ssbond_length)

    return (aa_code, dssp_values, dssp_limit, profile_matrix, profile_valid,
            centrality_matrix, centrality_valid, inflow_matrix,
            np.uint8(1 if has_target_inflow else 0), avg_b_factor, modres_code,
            ssbond_length, ssbond_present.astype(np.uint8), ssbond_intra,
            link_present, link_intra, link_length, cis_conformation, cis_follower)


@njit(cache=True)
def _residue_loop(
        n_residues, target_atom_matrix, target_number,
        chain_order, branch, is_protein, has_dssp, has_profiles, has_fuzzy, type_code,
        atoms, fuzzy, kd, rsa, rsa_present, ppi, hom, interface_scores, profile_valid_col,
        het_atoms, het_kind,
        milieu_threshold, intra_milieu_threshold, ligand_sphere,
        aa_code, dssp_values, dssp_limit, profile_matrix, profile_valid,
        centrality_matrix, centrality_valid, inflow_matrix, inflow_present,
        avg_b_factor, modres_code, ssbond_length, ssbond_present, ssbond_intra,
        link_present, link_intra, link_length, cis_conformation, cis_follower,
        is_alphafold, encoded, chain_result, chain_result_int):
    """Every residue of the target chain, straight into the encoded matrix.

    Returns (status, residue) - the residue being the one the distance loop gave
    up on, or -1. The row of a residue it gave up on is still written, because the
    python loop only abandoned the chain on the *next* iteration and so had already
    appended that row.
    """
    for index in range(encoded.shape[0]):
        for column in range(encoded.shape[1]):
            encoded[index, column] = np.nan

    for res_id in range(n_residues):
        status = _chain_distance_loop(
            res_id, target_atom_matrix, target_number,
            chain_order, branch, is_protein, has_dssp, has_profiles, has_fuzzy,
            type_code, atoms, fuzzy, kd, rsa, rsa_present, ppi, hom,
            interface_scores, profile_valid_col, het_atoms, het_kind,
            milieu_threshold, intra_milieu_threshold, ligand_sphere,
            chain_result, chain_result_int)

        row = encoded[res_id]

        # --- distances, once get_shortest_distances' own cutoff has been applied
        if status == _CHAIN_LOOP_OK:
            closest = np.inf
            for offset in range(_R_LIG_DIST, _R_DNA_DIST + 1):
                value = chain_result[offset]
                if not np.isnan(value) and value < closest:
                    closest = value
            if closest == np.inf:
                closest = 2.0
            if closest >= 1.2:
                for offset in range(6):
                    row[_O_LIG_DIST + offset] = chain_result[_R_LIG_DIST + offset]

        # --- what the distance loop aggregated
        if chain_result_int[_I_INTRA_PPI] != -1:
            row[_O_INTRA_PPI_DEGREE] = chain_result_int[_I_INTRA_PPI]
        if chain_result_int[_I_INTRA_HOM] != -1:
            row[_O_INTRA_HOM_DEGREE] = chain_result_int[_I_INTRA_HOM]
        if chain_result_int[_I_INTER_PPI] != -1:
            row[_O_INTER_PPI_DEGREE] = chain_result_int[_I_INTER_PPI]
        if chain_result_int[_I_INTER_HOM] != -1:
            row[_O_INTER_HOM_DEGREE] = chain_result_int[_I_INTER_HOM]

        if not np.isnan(chain_result[_R_INTRA_WEIGHTS]):
            row[_O_INTRA_MEDIAN_KD] = chain_result[_R_INTRA_MEDIAN_KD]
            row[_O_INTRA_MEDIAN_KD + 1] = chain_result[_R_INTRA_DW_KD]
            row[_O_INTRA_MEDIAN_KD + 2] = chain_result[_R_INTRA_MEDIAN_RSA]
            row[_O_INTRA_MEDIAN_KD + 3] = chain_result[_R_INTRA_DW_RSA]
        if not np.isnan(chain_result[_R_INTRA_INTERACTIONS_WEIGHTS]):
            row[_O_INTRA_INTERACTIONS_MEDIAN] = chain_result[_R_INTRA_INTERACTIONS_MEDIAN]
            row[_O_INTRA_INTERACTIONS_MEDIAN + 1] = chain_result[_R_INTRA_INTERACTIONS_DW]
        if not np.isnan(chain_result[_R_INTER_WEIGHTS]):
            row[_O_INTER_MEDIAN_KD] = chain_result[_R_INTER_MEDIAN_KD]
            row[_O_INTER_MEDIAN_KD + 1] = chain_result[_R_INTER_DW_KD]
            row[_O_INTER_MEDIAN_KD + 2] = chain_result[_R_INTER_MEDIAN_RSA]
            row[_O_INTER_MEDIAN_KD + 3] = chain_result[_R_INTER_DW_RSA]
            row[_O_INTER_INTERACTIONS_MEDIAN] = chain_result[_R_INTER_INTERACTIONS_MEDIAN]
            row[_O_INTER_INTERACTIONS_MEDIAN + 1] = chain_result[_R_INTER_INTERACTIONS_DW]

        # --- dssp, for the prefix of residues that had a row
        if res_id < dssp_limit:
            row[_O_RSA] = dssp_values[res_id, 0]
            row[_O_MC_ACC] = dssp_values[res_id, 1]
            row[_O_SC_ACC] = dssp_values[res_id, 2]
            row[_O_SSA] = dssp_values[res_id, 3]
            for offset in range(5):
                row[_O_PHI + offset] = dssp_values[res_id, 4 + offset]

        if profile_valid[res_id] == 1:
            for offset in range(profile_matrix.shape[1]):
                row[_O_PROFILE + offset] = profile_matrix[res_id, offset]

        if centrality_valid[res_id] == 1:
            for offset in range(centrality_matrix.shape[1]):
                row[_O_CENTRALITIES + offset] = centrality_matrix[res_id, offset]

        if inflow_present == 1:
            for offset in range(inflow_matrix.shape[1]):
                row[_O_INFLOW + offset] = inflow_matrix[res_id, offset]

        # --- bonds. The four flags are written whatever happens, because the python
        # version assigned True or False to them on every residue; link_length is
        # left alone because nothing ever set it (see the note at the call site).
        if ssbond_present[res_id] == 1:
            row[_O_SSBOND_LENGTH] = ssbond_length[res_id]
            row[_O_INTRA_SSBOND] = 1.0 if ssbond_intra[res_id] == 1 else 0.0
            row[_O_INTER_SSBOND] = 0.0 if ssbond_intra[res_id] == 1 else 1.0
        else:
            row[_O_INTRA_SSBOND] = 0.0
            row[_O_INTER_SSBOND] = 0.0

        if link_present[res_id] == 1:
            row[_O_INTRA_LINK] = 1.0 if link_intra[res_id] == 1 else 0.0
            row[_O_INTER_LINK] = 0.0 if link_intra[res_id] == 1 else 1.0
            row[_O_LINK_LENGTH] = link_length[res_id]
        else:
            row[_O_INTRA_LINK] = 0.0
            row[_O_INTER_LINK] = 0.0

        row[_O_CIS_CONFORMATION] = cis_conformation[res_id]
        row[_O_CIS_FOLLOWER] = cis_follower[res_id]

        row[_O_MODRES] = modres_code[res_id]
        row[_O_AA] = aa_code[res_id]

        if is_alphafold == 1:
            row[_O_PLDDT] = avg_b_factor[res_id]
        else:
            row[_O_B_FACTOR] = avg_b_factor[res_id]

        if status != _CHAIN_LOOP_OK:
            return status, res_id

    return _CHAIN_LOOP_OK, -1


# --- chain columns -----------------------------------------------------------
#
# The distance loop over a complex's chains used to walk a list of python tuples
# holding, per chain, a dssp list, a list of Interaction_profile objects, and
# chain-keyed dict lookups. None of that is addressable from a kernel, so every
# per-chain input is flattened here into a column indexed by the chain's position
# in CHAIN_ID_ORDER. Chain ids are single characters from that fixed 62 symbol
# alphabet, so the columns are 62 slots wide whatever the structure holds, and an
# absent chain is simply an empty row nothing indexes.

_CHAIN_TYPE_PROTEIN = 0
_CHAIN_TYPE_RNA = 1
_CHAIN_TYPE_DNA = 2
_CHAIN_TYPE_OTHER = 3

_HET_KIND_LIGAND = 0
_HET_KIND_METAL = 1
_HET_KIND_ION = 2

_EMPTY_ATOMS = np.zeros((0, 1, 4), dtype=np.float32)
_EMPTY_FUZZY = np.zeros((0, 0), dtype=np.float32)
_EMPTY_F8 = np.zeros(0, dtype=np.float64)
_EMPTY_U1 = np.zeros(0, dtype=np.uint8)
_EMPTY_I1 = np.zeros(0, dtype=np.int8)


def _chain_type_code(chain_type):
    # get_shortest_distances counts a chain it has no type for as Protein, and puts
    # Protein and Peptide in the same bucket.
    if chain_type is None or chain_type == 'Protein' or chain_type == 'Peptide':
        return _CHAIN_TYPE_PROTEIN
    if chain_type == 'RNA':
        return _CHAIN_TYPE_RNA
    if chain_type == 'DNA':
        return _CHAIN_TYPE_DNA
    return _CHAIN_TYPE_OTHER


def _hydropathy_column(residue_names: np.ndarray) -> np.ndarray:
    """Kyte-Doolittle value per residue, as accumulate_kd_rsa derived it."""
    column = np.empty(len(residue_names), dtype=np.float64)
    for index in range(len(residue_names)):
        try:
            one_letter = BIN_THREE_TO_ONE[residue_names[index]].decode('ascii')[0]
        except KeyError:
            one_letter = 'X'
        column[index] = HYDROPATHY[one_letter]
    return column


def _rsa_column(chain_dssp):
    """RSA per residue plus the mask of rows dssp actually produced."""
    if chain_dssp is None:
        return _EMPTY_F8, _EMPTY_U1
    rsa = np.zeros(len(chain_dssp), dtype=np.float64)
    present = np.zeros(len(chain_dssp), dtype=np.uint8)
    for index, row in enumerate(chain_dssp):
        if row is None:
            continue
        rsa[index] = row[0]
        present[index] = 1
    return rsa, present


def _het_kind_column(het_names: np.ndarray) -> np.ndarray:
    """Which distance bucket each hetero residue falls into.

    get_shortest_distances split the ligand distances by name, off the identifier
    the loop formatted per hetero residue; the split only depends on the name, so
    it is decided once here and the identifier never has to be built.
    """
    kinds = np.empty(len(het_names), dtype=np.int8)
    for index in range(len(het_names)):
        raw = het_names[index]
        name = raw.decode('ascii') if isinstance(raw, bytes) else str(raw)
        if name in METAL_ATOMS:
            kinds[index] = _HET_KIND_METAL
        elif name in ION_ATOMS:
            kinds[index] = _HET_KIND_ION
        else:
            kinds[index] = _HET_KIND_LIGAND
    return kinds


def build_structure_columns(parsed, dssp_ordered, profiles, chain_profile_lookups):
    """The per-chain inputs of the distance loop that do not depend on the target.

    Only two of the sixteen columns vary with which chain is being analysed: the
    branch a chain takes (is it the target?) and the target's row of the fuzzy
    distance matrix. The other fourteen are properties of the chain alone, and
    rebuilding them for every analysed chain made this quadratic in the number of
    protein chains - decoding every residue name into a hydropathy value once per
    chain per target, and re-deriving each chain's dssp and profile columns the
    same way. They are built once per structure here instead.
    """
    n_slots = MAX_NUMBER_OF_CHAINS

    is_protein = np.zeros(n_slots, dtype=np.uint8)
    has_dssp = np.zeros(n_slots, dtype=np.uint8)
    has_profiles = np.zeros(n_slots, dtype=np.uint8)
    type_code = np.full(n_slots, _CHAIN_TYPE_PROTEIN, dtype=np.int8)
    n_residues = np.zeros(n_slots, dtype=np.int64)

    atoms = NumbaList()
    kd = NumbaList()
    rsa = NumbaList()
    rsa_present = NumbaList()
    ppi = NumbaList()
    hom = NumbaList()
    interface_scores = NumbaList()
    profile_valid = NumbaList()
    het_atoms = NumbaList()
    het_kind = NumbaList()
    for _ in range(n_slots):
        atoms.append(_EMPTY_ATOMS)
        kd.append(_EMPTY_F8)
        rsa.append(_EMPTY_F8)
        rsa_present.append(_EMPTY_U1)
        ppi.append(_EMPTY_U1)
        hom.append(_EMPTY_U1)
        interface_scores.append(_EMPTY_F8)
        profile_valid.append(_EMPTY_U1)
        het_atoms.append(_EMPTY_ATOMS)
        het_kind.append(_EMPTY_I1)

    for number in parsed.chain_numbers.tolist():
        chain_id = parsed.chain_id(number)
        residue_names = parsed.residue_names[number]
        n_residues[number] = 0 if residue_names is None else len(residue_names)

        type_code[number] = _chain_type_code(parsed.type_name(number))
        is_protein[number] = 1 if parsed.is_protein(number) else 0

        if n_residues[number] > 0:
            atoms[number] = np.ascontiguousarray(parsed.residue_atoms[number], dtype=np.float32)
            kd[number] = _hydropathy_column(residue_names)

        if chain_id in dssp_ordered:
            has_dssp[number] = 1
            rsa[number], rsa_present[number] = _rsa_column(dssp_ordered[chain_id])

        if chain_id in profiles:
            has_profiles[number] = 1
        lookup = chain_profile_lookups.get(chain_id)
        if lookup is not None:
            chain_ppi, chain_hom, chain_scores, chain_valid, _matrix = lookup
            ppi[number] = np.ascontiguousarray(chain_ppi).astype(np.uint8)
            hom[number] = np.ascontiguousarray(chain_hom).astype(np.uint8)
            interface_scores[number] = np.ascontiguousarray(chain_scores, dtype=np.float64)
            profile_valid[number] = np.ascontiguousarray(chain_valid).astype(np.uint8)

        het_names = parsed.hetero_names[number]
        if het_names is not None and len(het_names) > 0:
            het_atoms[number] = np.ascontiguousarray(parsed.hetero_atoms[number], dtype=np.float32)
            het_kind[number] = _het_kind_column(het_names)

    static_flags = (is_protein, has_dssp, has_profiles, type_code, n_residues)
    static_columns = (atoms, kd, rsa, rsa_present, ppi, hom, interface_scores,
                      profile_valid, het_atoms, het_kind)
    return static_flags, static_columns


def build_chain_columns(parsed, structure_columns, fuzzy_dist_matrix, target_number):
    """The two target-dependent columns, on top of the structure's own.

    branch says which of the distance loop's two paths a chain takes, and fuzzy
    carries the target's row of the residue distance matrix; both change with the
    chain being analysed. Everything else comes straight from build_structure_columns.
    """
    (is_protein, has_dssp, has_profiles, type_code, n_residues) = structure_columns[0]
    (atoms, kd, rsa, rsa_present, ppi, hom, interface_scores, profile_valid,
     het_atoms, het_kind) = structure_columns[1]

    n_slots = MAX_NUMBER_OF_CHAINS
    branch = np.zeros(n_slots, dtype=np.int8)
    has_fuzzy = np.zeros(n_slots, dtype=np.uint8)
    fuzzy = NumbaList()
    for _ in range(n_slots):
        fuzzy.append(_EMPTY_FUZZY)

    for number in parsed.chain_numbers.tolist():
        if number != target_number and n_residues[number] > 0:
            branch[number] = 1
        elif number == target_number and parsed.residue_atoms[number] is not None:
            branch[number] = 2

        chain_fuzzy = fuzzy_row(fuzzy_dist_matrix, target_number, number)
        if chain_fuzzy is not None:
            fuzzy[number] = np.ascontiguousarray(chain_fuzzy, dtype=np.float32)
            has_fuzzy[number] = 1

    flags = (branch, is_protein, has_dssp, has_profiles, has_fuzzy, type_code)
    columns = (atoms, fuzzy, kd, rsa, rsa_present, ppi, hom, interface_scores,
               profile_valid, het_atoms, het_kind)
    return np.asarray(parsed.chain_numbers, dtype=np.int64), flags, columns


def analysis_chain(
        target_chain: str,
        config: Config,
        profiles: dict[str, list[Interaction_profile]],
        centralities,
        inflow_map: dict[str, np.ndarray], 
        analysis_dump,
        chain_profile_lookups: dict[str, tuple[np.ndarray, np.ndarray]] | None = None,
        structure_columns=None):
    dssp_ordered: dict[str, list[tuple[float, float, float, str]]]

    (pdb_id, parsed, fuzzy_dist_matrix, dssp, dssp_ordered, page_altered) = analysis_dump
    parsed: Parsed_structure
    target_number = chain_id_number(target_chain)
    errorlist = []

    if config.verbosity >= 5:
        if len(parsed.residue_names[target_number]) <= 1 or config.verbosity >= 7:
            config.logger.info(f'Call of analysis_chain: {pdb_id} {target_chain} {parsed.residue_names[target_number]=}')
        else:
            config.logger.info(f'Call of analysis_chain: {pdb_id} {target_chain} {len(parsed.residue_names[target_number])=}')

    structural_analysis_list: list[np.ndarray[float]] | None = []

    sub_loop_broken = False

    target_res_id: int

    total_loop_times = []

    if chain_profile_lookups is None:
        chain_profile_lookups = build_chain_profile_lookups(profiles)

    # Target-independent, so the caller builds it once per structure; rebuilding it
    # per analysed chain is what made this quadratic in the chain count.
    if structure_columns is None:
        structure_columns = build_structure_columns(
            parsed, dssp_ordered, profiles, chain_profile_lookups)

    # Everything getMinSubDist needs about a chain is constant across target residues:
    # res_chain is always target_chain, so its fuzzy_dist_matrix row is fixed, and both
    # coordinate_map lookups depend only on the chain. Resolved once here instead of on
    # every (residue, chain) pair.
    target_atom_lists = parsed.residue_atoms[target_number]
    target_residue_names = parsed.residue_names[target_number]

    # target_chain is fixed for the whole call, so every map keyed by it resolves once
    # here instead of on every residue. Presence is kept as an explicit flag rather than
    # a None check, so a stored None cannot change the branch that used to test `in`.
    has_target_dssp = target_chain in dssp_ordered
    target_dssp = dssp_ordered[target_chain] if has_target_dssp else None
    has_target_profiles = target_chain in profiles
    target_profiles = profiles[target_chain] if has_target_profiles else None
    # The encoded rows of those same profiles, built once per chain by
    # build_chain_profile_lookups. encode_to_numpy used to rebuild each row from the
    # object tree on every residue, which is the same 132 floats either way.
    target_lookup = chain_profile_lookups.get(target_chain)
    target_profile_matrix = target_lookup[4] if target_lookup is not None else None
    target_profile_valid = target_lookup[3] if target_lookup is not None else None
    has_target_centralities = target_chain in centralities
    target_centralities = centralities[target_chain] if has_target_centralities else None
    has_target_inflow = target_chain in inflow_map
    target_inflow = inflow_map[target_chain] if has_target_inflow else None
    target_modres = parsed.modres[target_number]
    has_target_modres = target_modres is not None
    # Unguarded in the original, so its KeyError must still only fire when the residue
    # loop would actually run.
    target_b_factors = parsed.b_factors[target_number] if len(target_residue_names) > 0 else None

    # Same idea for the other chain-keyed maps: none of these depend on the residue, so
    # they are resolved once per chain rather than on every (residue, chain) pair.
    # Presence flags again, so a stored None cannot flip a branch that tested `in`.
    chain_loop_iterations = 0

    # Every per-chain input the distance loop needs, as columns indexed by the
    # chain's position in CHAIN_ID_ORDER, so the loop can run as one kernel.
    chain_order, chain_flags, chain_columns = build_chain_columns(
        parsed, structure_columns, fuzzy_dist_matrix, target_number)
    (col_branch, col_is_protein, col_has_dssp, col_has_profiles, col_has_fuzzy,
     col_type_code) = chain_flags
    (col_atoms, col_fuzzy, col_kd, col_rsa, col_rsa_present, col_ppi, col_hom,
     col_interface_scores, col_profile_valid, col_het_atoms, col_het_kind) = chain_columns
    # The kernel indexes this one directly rather than through the column list.
    target_atom_matrix = (col_atoms[target_number] if target_number != -1
                          else _EMPTY_ATOMS)

    # Refilled by every residue rather than reallocated.
    chain_result = np.empty(_R_LEN, dtype=np.float64)
    chain_result_int = np.empty(_I_LEN, dtype=np.int64)

    n_residues = len(target_residue_names)

    # Everything the residue loop writes that does not come out of the distance
    # loop, gathered once for the chain. link_length is among them now: the python
    # version assigned the parsed distance to `residue.link_dist`, which is not one
    # of Residue's slots - Slotted_obj declares no __slots__, so that silently
    # created a shadow attribute and the encoded link_length was never written.
    residue_columns = _build_residue_columns(
        target_residue_names, target_dssp, has_target_dssp, dssp,
        target_profile_matrix, target_profile_valid, target_profiles, has_target_profiles,
        target_centralities, has_target_centralities,
        target_inflow, has_target_inflow,
        target_b_factors, target_modres, has_target_modres,
        parsed.ssbond[target_number], parsed.link[target_number],
        parsed.cis_conformation[target_number], parsed.cis_follower[target_number],
        target_chain)
    (aa_code, dssp_values, dssp_limit, profile_matrix, profile_valid,
     centrality_matrix, centrality_valid, inflow_matrix, inflow_present,
     avg_b_factor, modres_code, ssbond_length, ssbond_present, ssbond_intra,
     link_present, link_intra, link_length, cis_conformation, cis_follower) = residue_columns

    encoded = np.empty((n_residues, _ENCODED_WIDTH), dtype=np.float32)
    chain_result = np.empty(_R_LEN, dtype=np.float64)
    chain_result_int = np.empty(_I_LEN, dtype=np.int64)

    t_loop = time.time()
    loop_times = []

    chain_loop_iterations = n_residues * chain_order.shape[0]
    loop_status, broken_res_id = _residue_loop(
        n_residues, target_atom_matrix, np.int64(target_number),
        chain_order, col_branch, col_is_protein, col_has_dssp, col_has_profiles,
        col_has_fuzzy, col_type_code,
        col_atoms, col_fuzzy, col_kd, col_rsa, col_rsa_present, col_ppi, col_hom,
        col_interface_scores, col_profile_valid, col_het_atoms, col_het_kind,
        np.float32(config.milieu_threshold), np.float32(config.intra_milieu_threshold),
        np.float32(config.ligand_interest_sphere),
        aa_code, dssp_values, np.int64(dssp_limit), profile_matrix, profile_valid,
        centrality_matrix, centrality_valid, inflow_matrix, inflow_present,
        avg_b_factor, modres_code, ssbond_length, ssbond_present, ssbond_intra,
        link_present, link_intra, link_length, cis_conformation, cis_follower,
        np.uint8(1 if is_alphafold_model(pdb_id) else 0),
        encoded, chain_result, chain_result_int)

    t_loop = add_to_times(loop_times, t_loop, label = 'residue loop')

    # Both reports were emitted by a particular residue, so which of them happened -
    # and in which order - follows from the residue each belongs to. The dssp report
    # belongs to residue dssp_limit, the give-up report to broken_res_id, and within
    # one residue the distance loop ran before the dssp block.
    last_residue = broken_res_id if loop_status != _CHAIN_LOOP_OK else n_residues - 1
    dssp_reported = dssp and dssp_limit < n_residues and dssp_limit <= last_residue

    def _report_dssp():
        atomic_coverage = get_atomic_coverage(parsed, target_number)
        if isinstance(atomic_coverage, str):
            errorlist.append(f'Error in get_atomic_coverage: {pdb_id}, {target_chain}\n{atomic_coverage}')
        elif atomic_coverage > 0.5 and not page_altered: #DSSP does not work for chains with missing atomic coordinates or chains over the 10K atoms mark, no need for a warning here
            errorlist.append(f"dssp error: chain not in dssp_dict; {pdb_id}; {target_chain}")

    if dssp_reported and dssp_limit < last_residue:
        _report_dssp()

    if loop_status == _CHAIN_LOOP_OK:
        structural_analysis_list = list(encoded)
    else:
        broken_chain = CHAIN_ID_ORDER[chain_result_int[_I_BROKEN_CHAIN]]
        if loop_status == _CHAIN_LOOP_NO_INTER_PROFILES:
            errorlist.append(f'(1) interacting_chain not in profiles: {pdb_id} {broken_chain} {list(profiles.keys())}')
        else:
            errorlist.append(f'(2) target_chain not in profiles: {pdb_id} {broken_chain} {list(profiles.keys())}')
        # The python loop noticed only on its next turn, so it kept the rows it had
        # already appended - including the one it gave up on - and only discarded
        # the chain if there was a next residue at all.
        if broken_res_id < n_residues - 1:
            structural_analysis_list = None
        else:
            structural_analysis_list = list(encoded)

    if dssp_reported and dssp_limit >= last_residue:
        _report_dssp()

    total_loop_times = aggregate_times(total_loop_times, loop_times)

    if config.verbosity >= 5:
        if structural_analysis_list is None or config.verbosity >= 7:
            config.logger.info(f'Results of analysis_chain: {pdb_id} {target_chain} {structural_analysis_list}')
        else:
            config.logger.info(f'Results of analysis_chain: {pdb_id} {target_chain} {len(structural_analysis_list)=}')

    if config.verbosity >= 4:
        # The per-(residue, chain) timing in this loop costs ~3 us an iteration and is
        # unconditional, so this is what sizes it against the phase total.
        config.logger.info(
            f'[analysis_chain] {pdb_id} {target_chain}: chain_loop_iterations={chain_loop_iterations:_} '
            f'residues={len(target_residue_names)} chains={chain_order.shape[0]}')

    return target_chain, structural_analysis_list, errorlist, total_loop_times


def structuralAnalysis(
        structure_id: str,
        config: Config,
        homomer_set: set[str],
        model_path: str = None, keep_rin_files: bool =True,
        ):

    dssp: bool = True
    dssp_path: str = config.dssp_path
    rin_db_path: str = config.rin_db_path
    verbosity: int = config.verbosity

    times = []
    ta = time.time()

    if verbosity >= 4:  
        config.logger.info(f'Start structuralAnalysis of {structure_id=} {model_path=} {config.pdb_path=} {config.afdb=}')

    page, page_altered, path, _ = pdbParser.standardParsePDB(structure_id, config, return_10k_bool=True, get_is_local=True, model_path=model_path)

    ta = add_to_times(times, ta, label = 'get structure file')

    if verbosity >= 8:
        config.logger.info(f'\n{structure_id}\n\n{page}')

    parsed: Parsed_structure = parsePDB(page)

    if config.verbosity >= 5:
        for chain_number in parsed.chain_numbers.tolist():
            config.logger.info(f'{structure_id=} chain={parsed.chain_id(chain_number)} '
                               f'{len(parsed.residue_names[chain_number])=}')

    ta = add_to_times(times, ta, label = 'parse structure file')

    #distance_threshold = max([2. * config.short_distance_threshold, 2. * config.milieu_threshold])
    distance_threshold = 2. * config.short_distance_threshold

    fuzzy_dist_matrix = calcFuzzyDM(parsed, config, distance_threshold)

    ta = add_to_times(times, ta, label = 'calc fuzzy DM')

    if dssp:
        # write temp pdb file only if we had to fix the 10k atom bug or the file is not stored locally
        if page_altered or path is None:
            tmp_path = f'{config.temp_folder}/tmp_{structure_id}.pdb'
            f = open(tmp_path, 'w')
            f.write(page)
            f.close()
        else:
            tmp_path = path

        # call DSSP  -structure of dssp_dict: {Chain:{Res_id:(racc,ssa)}}
        dssp_dict, _ = calcDSSP(tmp_path, dssp_path, angles=True)

        dssp_ordered: dict[str, list[tuple[float, float, float, str]]] = {}
        for chain in dssp_dict:
            chain_number = chain_id_number(chain)
            if chain_number == -1 or parsed.res_contig[chain_number] is None:
                continue
            dssp_ordered[chain] = []
            for key in parsed.res_contig[chain_number].get_ordered_keys():
                dssp_ordered[chain].append(dssp_dict[chain].get_item(key))

        # remove tmp_file
        if page_altered:
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        if dssp_dict == {}:
            dssp = False

    ta = add_to_times(times, ta, label = 'calc DSSP')

    profiles: dict[str, Residue_Map[Interaction_profile]] = {}
    centralities = {}
    inflow_map = {}
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}
    IAmap = {}
    lookup_times = []

    # lookup reads the RIN files and is keyed by chain id from end to end. It runs
    # once per structure and is i/o bound, so it keeps its dict interface and gets
    # the maps rebuilt here rather than being pulled into the column form.
    chain_type_map = parsed.chain_type_map()
    res_contig_map = parsed.res_contig_map()
    try:
        lookup_out = lookup(structure_id, page, config, parsed.chain_keyed(parsed.ligands),
                            parsed.chain_keyed(parsed.metals), parsed.chain_keyed(parsed.ions),
                            res_contig_map,
                            rin_db_path, chain_type_map, homomer_set,
                            model_path=model_path, keep_tmp_files=keep_rin_files)
    except KeyError as err:
        raise KeyError(f'lookup failed for: {structure_id=} {parsed.chain_ids()=} {chain_type_map=} {err=}')

    if not isinstance(lookup_out, str):
        (profiles, centralities, inflow_map, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, IAmap, lookup_times) = lookup_out
        if verbosity >= 5:
            config.logger.info(f'Chains in profiles: {list(profiles.keys())}')

    times.append(lookup_times)

    ta = add_to_times(times, ta, label = 'RIN lookup')
                      
    if config.compute_ppi:
        interfaces: dict[str, dict[str, interface_package.Interface]] = calculate_interfaces(IAmap, parsed)
    else:
        interfaces: dict[str, dict[str, interface_package.Interface]] = {}

    ta = add_to_times(times, ta, label = 'calc interfaces')

    number_protein_chains = 0
    for chain_number in parsed.chain_numbers.tolist():
        if parsed.res_contig[chain_number] is None:
            continue
        if parsed.chain_type[chain_number] == -1:
            config.errorlog.add_warning(f'{parsed.chain_id(chain_number)} not in chain_type_map for: {structure_id}')
            continue
        if not parsed.is_protein(chain_number):
            continue
        number_protein_chains += 1

    analysis_dump = (structure_id, parsed, fuzzy_dist_matrix, dssp, dssp_ordered, page_altered)
    
    structural_analysis_dict: dict[str, list[np.ndarray[float]]] = {}
    total_analysis_times = []

    # Depends only on profiles, so it is derived once per complex rather than once per
    # analysed chain.
    chain_profile_lookups = build_chain_profile_lookups(profiles)
    # Same idea one level up: these columns describe the chains, not the chain
    # being analysed, so they are built once for the whole complex instead of
    # once per analysed chain.
    structure_columns = build_structure_columns(
        parsed, dssp_ordered, profiles, chain_profile_lookups)

    for chain_number in parsed.chain_numbers.tolist():
        if parsed.res_contig[chain_number] is None or not parsed.is_protein(chain_number):
            continue
        target_chain = parsed.chain_id(chain_number)

        target_chain, chain_structural_analysis_list, chain_errorlist, chain_analysis_times = analysis_chain(target_chain, config, profiles, centralities, inflow_map, analysis_dump, chain_profile_lookups=chain_profile_lookups, structure_columns=structure_columns)
        aggregate_times(total_analysis_times, chain_analysis_times)
        if chain_structural_analysis_list is None:
            config.errorlog.add_warning(f'Structural analysis failed for {structure_id=} {target_chain=} chain_type={parsed.type_name(chain_number)}')
            continue
        structural_analysis_dict[target_chain] = chain_structural_analysis_list
        if config.verbosity >= 5:
            config.logger.info(f'Received structural analysis results (1) from: {structure_id} {target_chain} {len(chain_structural_analysis_list)=} {chain_errorlist=}')

    times.append(total_analysis_times)

    ta = add_to_times(times, ta, label = 'analyze chains')

    return structural_analysis_dict, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, interfaces, times

def trim_down_backmaps(config: Config, proteins):
    removed_backmaps = 0
    for structure_id in proteins.complexes:
        chains_with_interfaces = set()
        for chain_a in proteins.complexes[structure_id].interfaces:
            chains_with_interfaces.add(chain_a)
            if chain_a not in proteins.structures[structure_id]:
                continue
            
            for chain_b in proteins.complexes[structure_id].interfaces[chain_a]:
                new_back_maps = {}
                for res in proteins.complexes[structure_id].interfaces[chain_a][chain_b].residues:
                    for prot_id in proteins.structures[structure_id][chain_a].backmaps:
                        if prot_id not in new_back_maps:
                            new_back_maps[prot_id] = Residue_Map()
                        try:
                            new_back_maps[prot_id].add_item(res, proteins.structures[structure_id][chain_a].backmaps[prot_id].get_item(res))
                        except (KeyError, IndexError):
                            pass
            proteins.structures[structure_id][chain_a].backmaps = new_back_maps

        for chain in proteins.complexes[structure_id].chains:
            if chain in chains_with_interfaces:
                continue
            if structure_id in proteins.structures:
                if chain in proteins.structures[structure_id]:
                    removed_backmaps += len(proteins.structures[structure_id][chain].backmaps)
                    proteins.structures[structure_id][chain].backmaps = {}
                    
    if config.verbosity >= 4:
        config.logger.info(f'Removed backmaps: {removed_backmaps}')
    return

# How long a chunk whose worker failed to start waits before it goes out again, by attempt.
# The chunk waits; the driver does not. These used to be a time.sleep inside the dispatch loop,
# once per failed start and one after another, and on a whole-PDB chunk 512 failed starts
# arriving in bursts stalled the driver for up to 4.3 h: every task that finished meanwhile left
# its core idle, because nothing was there to hand it the next chunk.
ANNOTATION_RETRY_BACKOFF_SECONDS = (30, 60, 120)

# How often paraAnnotate logs its occupancy at verbosity >= 2.
ANNOTATION_HEARTBEAT_SECONDS = 300


def annotation_retry_backoff(attempts: int) -> float:
    """Seconds a chunk waits before its attempt number 'attempts' (1-based) goes out."""
    index = min(max(attempts, 1), len(ANNOTATION_RETRY_BACKOFF_SECONDS)) - 1
    return ANNOTATION_RETRY_BACKOFF_SECONDS[index]


def take_due_retries(deferred: list, now: float) -> tuple[list, list]:
    """Split deferred (not_before, chunk, attempts) entries into (due, still_waiting)."""
    due = [entry for entry in deferred if entry[0] <= now]
    waiting = [entry for entry in deferred if entry[0] > now]
    return due, waiting


def annotation_chunksize(remaining_cost: int, config: Config) -> int:
    """Atoms one annotation chunk may carry.

    The cost-proportional part splits the remaining work over the cores with enough
    chunks left (the factor 8) that a slow one cannot leave the others idle at the end.
    On its own it scales the wrong way: the bigger the input, the bigger every chunk,
    because total_cost is the numerator. A whole-PDB run landed on 2.08M atoms per
    chunk, and annotate holds a chunk's whole structural_analysis_list before inserting
    it, so proc_n workers each sat on a chunk that large. That is what left the cgroup
    without the memory to start the next worker process and killed the run with ray's
    "Failed to startup worker after retrying 5 times".

    The cap breaks that dependency: per-worker peak memory is the same on a whole-PDB
    input as on a small one, paid for with more, cheaper chunks.

    Never 0. A remainder that floor-divides away would otherwise dispatch a chunk per
    complex, since chunk_cost >= 0 holds immediately.
    """
    chunksize = remaining_cost // (8 * config.proc_n)
    cap = getattr(config, 'annotation_chunk_atom_cap', 0)
    if cap > 0:
        chunksize = min(chunksize, cap)
    return max(1, chunksize)


def resolve_annotation_refs(ready, chunk_of_ref, attempts_of_ref, config, get=None):
    """Resolve the refs ray.wait handed back, keeping the run alive through a failed worker start.

    Returns (outs, retry_chunks): the results of the chunks that ran, and the
    (chunk, attempts) pairs whose worker ray never managed to start and which may
    therefore go out again.

    Resolved one ref at a time rather than as a batch. A RaySystemError names one lease,
    and the batched ray.get this replaces threw away the results of every other ready
    chunk along with it.

    Only RaySystemError is retried, and the distinction matters: it means ray could not
    bring a worker process up, so the task never ran and nothing of it was inserted. A
    RayTaskError or a WorkerCrashedError means the task did run and may have committed
    part of its residues through remote_insertResidues, so handing the chunk out again
    would insert them twice - those keep propagating and still fail the run.

    The usual cause is a start that does not register with the raylet in time: every chunk
    starts a fresh process (max_calls=1) from an environment on NFS, and many at once make
    some of them too slow. Asking again later gets the chunk through, instead of losing the
    hours that came before the stage.

    This function never waits: the caller defers the returned chunks (see
    annotation_retry_backoff) and keeps dispatching meanwhile. 'get' is injectable so the retry
    policy can be exercised without ray.
    """
    if get is None:
        get = ray.get

    outs = []
    retry_chunks: list[tuple[list[str], int]] = []
    retry_limit = getattr(config, 'annotation_retry_limit', 3)

    for ref in ready:
        chunk = chunk_of_ref.pop(ref, None)
        attempts = attempts_of_ref.pop(ref, 0)
        try:
            outs.append(get(ref))
        except ray.exceptions.RaySystemError as err:
            if chunk is None or attempts >= retry_limit:
                config.errorlog.add_error(
                    f'Ray failed to start a worker for an annotation chunk and it will '
                    f'not be retried again (attempts: {attempts}, limit: {retry_limit}). '
                    f'Check the raylet logs for workers that did not register in time; '
                    f'ray_worker_register_timeout_seconds and ray_worker_startup_concurrency '
                    f'govern that. {err}'
                )
                raise
            if config.verbosity >= 1:
                config.logger.warning(
                    f'Ray failed to start a worker for an annotation chunk of '
                    f'{len(chunk)} complexes (attempt {attempts + 1} of {retry_limit}); '
                    f'retrying in {annotation_retry_backoff(attempts + 1):.0f} s. {err}'
                )
            retry_chunks.append((chunk, attempts + 1))

    return outs, retry_chunks


def paraAnnotate(config: Config, proteins, indel_analysis_follow_up=False):

    ta = time.time()
    times = []

    # new structure for paralell annotation: give each annotation process a pdb and dict of chain and structure information,
    # all chains without structure information shall go into the database as fully annotated structures, which are not aligned (yet)
    # these are unmapped chains, or protein interaction partners, their presence in the database is necesary for the computation of the structural neighborhood
    # a new structure-structure entry has to be created and has to be inserted into the database before the residue insertion
    # homooligomer information is not present for such cases, this information has to be updated, if the chain is mapped at a later iteration or run

    size_list: list[tuple[str, int]] = []

    structure_list = proteins.get_structure_list()

    complex_list: list[str] = proteins.get_complex_list()

    n_of_stored_complexes = 0
    n_of_comps = 0

    total_cost: int = 0

    for pdb_id in complex_list:
        if proteins.is_complex_stored(pdb_id):
            n_of_stored_complexes += 1
            continue
        #s: int = len(proteins.get_complex_chains(pdb_id, only_protein=True))
        ac: int = proteins.get_atom_count(pdb_id)
        if ac is None:
            config.errorlog.add_warning('Complex with None as atom count: %s' % pdb_id)
            ac = 0
        
        size_list.append((pdb_id,ac))
        total_cost += ac

        if config.verbosity >= 6:
            config.logger.info(f'PDB entry {pdb_id} added to size_list, {ac=}')
        n_of_comps += 1

    chunksize: int = annotation_chunksize(total_cost, config)

    if config.verbosity >= 2:
        config.logger.info(f'Starting structural analysis with {n_of_comps} complexes. {n_of_stored_complexes} complexes are already stored. Chunksize: {chunksize}, n_of_chain_thresh: {config.n_of_chain_thresh}')
        config.logger.info(f'Total amount in structure_list: {len(structure_list)}. Total cost: {total_cost:_}, Threads: {config.proc_n}')

    size_list = sorted(size_list, key=lambda x:x[1] ,reverse=True)

    anno_result_ids = []

    if parallel_db_access(config):
        ray_init(config, stop_existing = False)
        conf_dump = ray.put(config)
    else:
        conf_dump = None

    n_started = 0
    started_cost = 0
    current_chunk: list[str] = []
    chunk_cost = 0

    # The chunk behind every pending ObjectRef, so a chunk whose worker ray never
    # managed to start can be handed out again instead of ending the run. Entries are
    # dropped as soon as the ref resolves, so this holds at most proc_n chunks.
    chunk_of_ref: dict = {}
    # Resubmissions already spent on that chunk, carried across from the failed ref.
    attempts_of_ref: dict = {}
    # When each pending ref went out, for the stage's utilisation summary.
    dispatch_time_of_ref: dict = {}
    first_dispatch_at = None

    def dispatch(chunk: list[str], target_list: list, cost_of_chunk: int, attempts: int = 0) -> None:
        """Hand one chunk to a worker, or run it here when ray is not in play.

        cost_of_chunk is the chunk's whole cost, not the last complex's: started_cost
        is what the dynamic chunksize below subtracts from total_cost, and adding one
        complex per chunk left it barely moving over a whole run. A resubmission passes
        0 - that cost was already counted when the chunk first went out.
        """
        nonlocal n_started, started_cost, first_dispatch_at
        if parallel_db_access(config):
            ref = annotate_wrapper.remote(conf_dump, chunk)
            chunk_of_ref[ref] = chunk
            attempts_of_ref[ref] = attempts
            dispatch_time_of_ref[ref] = time.time()
            if first_dispatch_at is None:
                first_dispatch_at = dispatch_time_of_ref[ref]
            target_list.append(ref)
            n_started += 1
            started_cost += cost_of_chunk
        else:
            target_list.append(annotate(config, chunk))

    input_procession_counter = 0
    for pdb_id, ac in size_list:
        if config.verbosity >= 5:
            config.logger.info(f'Scheduling remote structural analysis: {pdb_id=} {ac:_=}')
        input_procession_counter += 1
        cost = ac
        if config.proc_n < 2:
            cost = 0

        current_chunk.append(pdb_id)

        chunk_cost += cost
        if chunk_cost >= chunksize:
            dispatch(current_chunk, anno_result_ids, chunk_cost)
            current_chunk = []
            chunk_cost = 0

        if n_started >= config.proc_n:
            break

    if config.proc_n < 2 and len(current_chunk) > 0:
        anno_result_ids.append(annotate(config, current_chunk))
        current_chunk = []
        chunk_cost = 0

    if config.verbosity >= 2:
        ta = add_to_times(times, ta, label = f"{n_started=} {input_procession_counter=}")

    max_comp_time_structure = None
    struct_comp_times = []

    total_computed = 0

    loop_counter = 0

    #gc.disable()

    total_times = []

    # Chunks whose worker failed to start, as (not_before, chunk, attempts). They wait out their
    # backoff here while the loop keeps dispatching and collecting; see annotation_retry_backoff.
    deferred: list[tuple[float, list[str], int]] = []
    start_failures = 0
    # Occupancy of the stage, for the heartbeat and the closing summary. Before these the stage
    # logged nothing between its first and its last line at verbosity 3, so a CPU drop in it
    # could only be seen on the machine, never in the log.
    stage_t0 = time.time()
    next_heartbeat = stage_t0 + ANNOTATION_HEARTBEAT_SECONDS
    busy_seconds = 0.0

    while True:
        loop_counter += 1
        chunksize: int = annotation_chunksize(total_cost - started_cost, config)  # dynamically adjust chunksize
        finished = 0

        not_ready = []
        if len(anno_result_ids) > 0:
            if config.verbosity >= 5:
                config.logger.info(f'Before wait, loop counter: {loop_counter} {n_started=} {input_procession_counter=} ({len(size_list)=})')

            if parallel_db_access(config):
                # Woken in time for the next deferred chunk even if no task returns meanwhile.
                timeout = None
                if deferred:
                    timeout = max(0.0, min(entry[0] for entry in deferred) - time.time())
                ready, not_ready = ray.wait(anno_result_ids, timeout = timeout)

                if config.verbosity >= 5:
                    config.logger.info(f'Wait happened: {len(ready)=} {len(not_ready)=}')

                returned_at = time.time()
                for ref in ready:
                    dispatched_at = dispatch_time_of_ref.pop(ref, None)
                    if dispatched_at is not None:
                        busy_seconds += returned_at - dispatched_at

                chunk_struct_outs, retry_chunks = resolve_annotation_refs(
                    ready, chunk_of_ref, attempts_of_ref, config)
                # Every chunk going back out had its slot freed when its worker failed, and the
                # refill below hands that slot to new work while the chunk waits.
                n_started -= len(retry_chunks)
                start_failures += len(retry_chunks)
                for chunk, attempts in retry_chunks:
                    deferred.append((returned_at + annotation_retry_backoff(attempts), chunk, attempts))
            else:
                chunk_struct_outs = anno_result_ids

            for anno_times, max_comp_time in chunk_struct_outs:
                n_started -= 1  # finished chunks free 1 cost
                finished += 1
                struct_comp_times.append(max_comp_time)
                total_times = aggregate_times(total_times, anno_times)
        elif deferred:
            # Nothing is running and everything left is waiting out a backoff: the one case in
            # which the driver has nothing better to do than wait.
            time.sleep(max(0.0, min(entry[0] for entry in deferred) - time.time()))

        new_anno_result_ids = []

        # Deferred chunks whose backoff is over go out before new work.
        due, deferred = take_due_retries(deferred, time.time())
        for _not_before, chunk, attempts in due:
            dispatch(chunk, new_anno_result_ids, 0, attempts = attempts)

        for pdb_id, ac in size_list[input_procession_counter:]:
            input_procession_counter += 1
            cost = ac
            if config.proc_n < 2:
                cost = 0
            current_chunk.append(pdb_id)
            chunk_cost += cost
            if chunk_cost >= chunksize:
                dispatch(current_chunk, new_anno_result_ids, chunk_cost)
                current_chunk = []
                chunk_cost = 0
            if n_started >= config.proc_n:
                break

        if (input_procession_counter >= len(size_list)) and (len(current_chunk) > 0):
            dispatch(current_chunk, new_anno_result_ids, chunk_cost)
            current_chunk = []
            chunk_cost = 0

        anno_result_ids = new_anno_result_ids + not_ready

        if finished > 0:
            total_computed += finished
            if config.verbosity >= 5:
                config.logger.info(f'Newly finished structural analysis of {finished}. Total: {total_computed} Current package: {len(anno_result_ids)}, Amount of new started processes: {len(new_anno_result_ids)} Amount of pending processes: {len(not_ready)}')
                config.logger.info(f'RAM memory % used: {psutil.virtual_memory()[2]}')

        if config.verbosity >= 2 and time.time() >= next_heartbeat:
            next_heartbeat = time.time() + ANNOTATION_HEARTBEAT_SECONDS
            config.logger.info(
                f'[paraAnnotate] t={time.time() - stage_t0:.0f}s in_flight={len(anno_result_ids)}/{config.proc_n} '
                f'finished={total_computed} complexes_queued={max(0, len(size_list) - input_procession_counter)} '
                f'deferred={len(deferred)} start_failures={start_failures}')

        if (len(anno_result_ids) == 0 and input_procession_counter >= len(size_list)
                and not deferred):
            break

    if config.verbosity >= 2 and first_dispatch_at is not None:
        # Share of the cores the stage kept busy between its first dispatch and its last
        # return - the figure a CPU drop shows up in.
        wall = max(time.time() - first_dispatch_at, 1e-9)
        config.logger.info(
            f'[paraAnnotate] {total_computed} chunks over {wall:.0f}s, core utilisation '
            f'{100.0 * busy_seconds / (wall * max(1, config.proc_n)):.1f}% of {config.proc_n}, '
            f'{start_failures} failed worker starts retried')

    if config.verbosity >= 2:
        times.append(total_times)
        if len(struct_comp_times) > 0:
            max_comp_time_structure = sorted(struct_comp_times, key=lambda x:x[1], reverse=True)[0]
        ta = add_to_times(times, ta, label = f'Longest computation for: {max_comp_time_structure=}.')

    del conf_dump

    if config.verbosity >= 5:
        try:
            times_stat_file = f'{config.temp_folder}/../struct_wise_times.tsv'
            config.logger.info(f'Writing struct wise times to: {times_stat_file}')
            with open(times_stat_file, 'w') as f:
                f.write('\n'.join(f'{structure_id}\t{comp_time}' for structure_id, comp_time in struct_comp_times))
        except OSError:
            pass

    #insertComplexes(proteins, config)

    if config.verbosity >= 2:
        ta = add_to_times(times, ta, label = 'insertComplexes')
        print_times(times, label = 'paraAnnotate', logger=config.logger)

    classification(proteins, config, None)

    if not indel_analysis_follow_up:
        if config.verbosity >= 3:
            config.logger.info('Proteins object semi deconstruction')
        proteins.semi_deconstruct()

    gc.collect()
        
    return 0


# Return True, if the template is not good
def weakCriteria(seq, res, rel_aln_len, seq_thresh, res_thresh, rel_aln_len_thresh):

    if float(seq) <= 1.0:
        seq = float(seq) * 100.0
    if float(seq) < seq_thresh:
        return True
    if float(res) > res_thresh:
        return True
    if float(rel_aln_len) < rel_aln_len_thresh:
        return True
    return False


# template-structure: [pdb_id,seq_id,chain,aln_length,resolution,ligands,r-value,templateSelectionScore]


def good(structure, seq_threshold, resolution_threshold, cov_threshold):
    seq_id = structure['Seq_Id']
    resolution = structure['Resolution']
    cov = structure['Coverage']
    if weakCriteria(seq_id, resolution, cov, seq_threshold, resolution_threshold, cov_threshold):
        return False
    return True

