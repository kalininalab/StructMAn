import math
import os
import subprocess
import sys
import time
import traceback
import gc
import psutil
import numpy as np
from numba import njit

import ray

from structman.lib import pdbParser, rin
from structman.lib.database.retrieval import getStoredResidues
from structman.lib.database.insertion_lib import insertResidues, remote_insertResidues,insertClassifications, insertComplexes, insert_interface_residues, insert_interfaces
from structman.base_utils.base_utils import median, distance, pack, unpack, is_alphafold_model, alphafold_model_id_to_file_path, aggregate_times, print_times
from structman.lib.sdsc.consts.residues import CORRECT_COUNT, THREE_TO_ONE, BLOSUM62, ONE_TO_THREE, RESIDUE_MAX_ACC, HYDROPATHY, METAL_ATOMS, ION_ATOMS
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
        except:
            blosum_value = 0.6 - float(BLOSUM62[(aac[-1], aac[0])]) / 10
        if blosum_value < 0.0:
            blosum_value = 0.0
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


def calcDSSP(path, DSSP, angles=False, verbosity_level=0):

    dssp_dict = {}
    errorlist = []

    try:
        p = subprocess.Popen([DSSP, path], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
    except:
        out = ''
        err = '%s %s' % (str(sys.exc_info()[0]), str(sys.exc_info()[1]))
    # Alert user for errors
    if err.strip():
        if not out.strip():
            if verbosity_level >= 4:
                print('DSSP failed to produce an output\n%s\n%s\n' % (path, err))
            errorlist.append('DSSP failed to produce an output\n%s\n%s\n' % (path, err))
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
        try:
            acc = float(line[34:38])

            main_chain_acc = float(line[38:42])

            side_chain_acc = float(line[42:46])
        except:
            errorlist.append('DSSP float conversion failed for: \n%s\nChain: %s, Residue: %s,%s\n%s\n\n' % (path, chain, res, aa_type_one_letter, line))
            racc = None
            relative_main_chain_acc = None
            relative_side_chain_acc = None

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
        except:
            errorlist.append('DSSP (%s) failed for: \n%s\nChain: %s, Residue: %s,%s\n%s\n\n' % (DSSP, path, chain, res, aa_type_one_letter, line))
            racc = None
            relative_main_chain_acc = None
            relative_side_chain_acc = None

        try:
            if len(line) > 123:
                phi = float(line[111:117].replace(" ", ""))
                psi = float(line[117:123].replace(" ", ""))
        except:
            errorlist.append('DSSP failed for (angles): \n%s\nChain: %s, Residue: %s,%s\n%s\n\n' % (path, chain, res, aa_type_one_letter, line))
            phi = None
            psi = None

        if chain not in dssp_dict:
            dssp_dict[chain] = residue_package.Residue_Map()
        dssp_dict[chain].add_item(res, (racc, relative_main_chain_acc, relative_side_chain_acc, ssa))
        if angles:
            dssp_dict[chain].add_item(res, (racc, relative_main_chain_acc, relative_side_chain_acc, ssa, phi, psi))
    return dssp_dict, errorlist


def parsePDB(input_page: str) -> tuple[
    dict[str, tuple[residue_package.Residue_Map, residue_package.Residue_Map]],
    dict[str, residue_package.Residue_Map[tuple[int, str]]],
    set[tuple[str, int | str]],
    set[tuple[str, int | str]],
    set[tuple[str, int | str]],
    dict,
    dict[str, str],
    list,
    dict,
    dict,
    dict,
    dict,
    dict,
    dict
    ]:
    """
    Parses a PDB-file and takes all atomic coordinates.

    Input:
    input_page: String ; content of a pdb file

    Output:
    coordinate_map: ,{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
    """
    
    lines = input_page.split('\n')

    coordinate_map: dict[str, tuple[residue_package.Residue_Map, residue_package.Residue_Map]] = {}

    modres_map = {}
    res_contig_map: dict[str, residue_package.Residue_Map[tuple[int, str]]] = {}
    contig_help_map: dict[str, int] = {}

    ssbond_map = {}
    link_map = {}
    cis_conformation_map = {}
    cis_follower_map = {}

    chain_type_map: dict[str, str] = {}
    chain_type = '-'
    chainlist = []

    peptide_count = {}

    box_map = {}

    ligands: set[tuple[str, int | str]] = set()
    metals: set[tuple[str, int | str]] = set()
    ions: set[tuple[str, int | str]] = set()
    rare_residues = set()
    b_factors = {}

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
                    res_nr = line[18:23].replace(" ", "")
                    try:
                        res_nr = int(res_nr)
                    except:
                        pass
                    res_name = line[24:27].replace(" ", "")
                    if not chain_id in modres_map:
                        modres_map[chain_id] = {}
                    if not res_nr in modres_map[chain_id]:
                        modres_map[chain_id][res_nr] = res_name

                if record_name == 'SSBOND':
                    chain_1 = line[15]
                    res_nr_1 = line[17:22].replace(" ", "")
                    try:
                        res_nr_1 = int(res_nr_1)
                    except:
                        pass
                    chain_2 = line[29]
                    res_nr_2 = line[31:36].replace(" ", "")
                    try:
                        res_nr_2 = int(res_nr_2)
                    except:
                        pass

                    ssbond_len = float(line[73:78].replace(' ', ''))
                    ssbond_map[(chain_1, res_nr_1)] = (chain_2, res_nr_2, ssbond_len)
                    ssbond_map[(chain_2, res_nr_2)] = (chain_1, res_nr_1, ssbond_len)

                if record_name == 'LINK  ':
                    atom_1 = line[12:16].replace(" ", "")
                    res_name_1 = line[17:20].replace(" ", "")
                    res_nr_1 = line[22:27].replace(" ", "")
                    try:
                        res_nr_1 = int(res_nr_1)
                    except:
                        pass                    
                    chain_1 = line[21]

                    atom_2 = line[42:46].replace(" ", "")
                    res_name_2 = line[47:50].replace(" ", "")
                    res_nr_2 = line[52:57].replace(" ", "")
                    try:
                        res_nr_2 = int(res_nr_2)
                    except:
                        pass
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

                    link_map[(chain_1, res_nr_1)] = (atom_1, res_name_1, atom_2, res_name_2, res_nr_2, chain_2, link_dist)
                    link_map[(chain_2, res_nr_2)] = (atom_2, res_name_2, atom_1, res_name_1, res_nr_1, chain_1, link_dist)

                if record_name == 'CISPEP':
                    res_name_1 = line[11:14].replace(" ", "")
                    chain_1 = line[15]
                    res_nr_1 = line[17:22].replace(" ", "")
                    try:
                        res_nr_1 = int(res_nr_1)
                    except:
                        pass
                    res_name_2 = line[25:28].replace(" ", "")
                    chain_2 = line[29]
                    res_nr_2 = line[31:36].replace(" ", "")
                    try:
                        res_nr_2 = int(res_nr_2)
                    except:
                        pass
                    angle = float(line[53:59].replace(" ", ""))

                    cis_conformation_map[(chain_1, res_nr_1)] = (res_name_2, chain_2, res_nr_2, angle)
                    cis_follower_map[(chain_2, res_nr_2)] = (res_name_1, chain_1, res_nr_1, angle)

                chain_id = line[21]

                if not chain_id in modres_map:
                    modres_map[chain_id] = {}

                try:
                    res_nr = int(line[22:27].replace(" ", ""))  # [22:27] includes the insertion_code
                except:
                    res_nr = line[22:27].replace(" ", "")  # [22:27] includes the insertion_code

                if record_name == "ATOM" or record_name == "HETATM":
                    altLoc = line[16]
                    if firstAltLoc is None and altLoc != ' ':
                        firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                    if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue

            if len(line) > 60:
                try:
                    b_factor = float(line[61:67].replace(" ", ""))
                except:
                    b_factor = 0.0
            if record_name == "ATOM" or record_name == 'MODRES' or record_name == "HETATM":
                if chain_id not in chain_type_map:
                    chainlist.append(chain_id)
                    if record_name == "ATOM" or record_name == 'MODRES':
                        if len(res_name) == 1:
                            chain_type = "RNA"
                        elif len(res_name) == 2:
                            chain_type = "DNA"
                        elif len(res_name) == 3:
                            chain_type = "Protein"
                        chain_type_map[chain_id] = chain_type
                        peptide_count[chain_id] = [0, 0]
                    elif record_name == 'HETATM':
                        if res_name in THREE_TO_ONE or not sdsc_utils.boring(res_name) or (res_name in rare_residues):
                            chain_type_map[chain_id] = "Protein"
                            peptide_count[chain_id] = [0, 0]
                        elif len(res_name) == 1:
                            chain_type_map[chain_id] = "RNA"
                            peptide_count[chain_id] = [0, 0]
                        elif len(res_name) == 2:
                            chain_type_map[chain_id] = "DNA"
                            peptide_count[chain_id] = [0, 0]
                        else:
                            chain_type_map[chain_id] = 'MISC'


                if record_name == 'HETATM':
                    if res_name in THREE_TO_ONE or not sdsc_utils.boring(res_name) or (res_name in rare_residues):  # For a hetero peptide 'boring' hetero amino acids are allowed as well as other non boring molecules not in THREE_TO_ONE, which are hopefully some kind of anormal amino acids
                        peptide_count[chain_id][1] += 1
                    elif len(res_name) < 3:
                        peptide_count[chain_id][0] += 1

                elif record_name == "ATOM" or record_name == 'MODRES':
                    peptide_count[chain_id][0] += 1
                    if len(res_name) == 1:
                        chain_type_map[chain_id] = "RNA"
                    elif len(res_name) == 2:
                        chain_type_map[chain_id] = "DNA"

            if record_name == "ATOM":
                if len(line) > 50 and not atom_name[0] in ('H', 'D'):
                    x = float(line[30:38].replace(" ", ""))
                    y = float(line[38:46].replace(" ", ""))
                    z = float(line[46:54].replace(" ", ""))
                    if chain_id not in coordinate_map:
                        coordinate_map[chain_id] = [residue_package.Residue_Map(), residue_package.Residue_Map()]
 
                        box_map[chain_id] = [x, x, y, y, z, z]

                    if not coordinate_map[chain_id][0].contains(res_nr):
                        coordinate_map[chain_id][0].add_item(res_nr, (res_name, sdsc_utils.SparseArray()))

                    coordinate_map[chain_id][0].get_item(res_nr)[1].insert(atom_nr, (atom_name, np.array([x, y, z])))


                    if chain_id not in b_factors:
                        b_factors[chain_id] = {}
                    if res_nr not in b_factors[chain_id]:
                        b_factors[chain_id][res_nr] = []
                    b_factors[chain_id][res_nr].append(b_factor)

                    if x < box_map[chain_id][0]:
                        box_map[chain_id][0] = x
                    if x > box_map[chain_id][1]:
                        box_map[chain_id][1] = x
                    if y < box_map[chain_id][2]:
                        box_map[chain_id][2] = y
                    if y > box_map[chain_id][3]:
                        box_map[chain_id][3] = y
                    if z < box_map[chain_id][4]:
                        box_map[chain_id][4] = z
                    if z > box_map[chain_id][1]:
                        box_map[chain_id][4] = z

                    if chain_id not in res_contig_map:
                        res_contig_map[chain_id] = residue_package.Residue_Map()
                        res_contig_map[chain_id].add_item(res_nr, (1, res_name))
                        contig_help_map[chain_id] = 1
                    elif not res_contig_map[chain_id].contains(res_nr):
                        contig_help_map[chain_id] += 1
                        res_contig_map[chain_id].add_item(res_nr, (contig_help_map[chain_id], res_name))

            if record_name == "HETATM":
                if len(line) > 50:
                    x = float(line[30:38].replace(" ", ""))
                    y = float(line[38:46].replace(" ", ""))
                    z = float(line[46:54].replace(" ", ""))
                    if chain_id not in coordinate_map:
                        coordinate_map[chain_id] = [residue_package.Residue_Map(), residue_package.Residue_Map()]
                        box_map[chain_id] = [x, x, y, y, z, z]

                    if res_nr in modres_map[chain_id] or (res_name in THREE_TO_ONE) or (res_name in rare_residues):  # If it is a modified residue, than add it to the normal residues...
                        if atom_name[0] in ('H', 'D'):
                            continue
                        if not res_nr in modres_map[chain_id]:
                            modres_map[chain_id][res_nr] = res_name
                        if not coordinate_map[chain_id][0].contains(res_nr):
                            coordinate_map[chain_id][0].add_item(res_nr, (res_name, sdsc_utils.SparseArray()))
                        coordinate_map[chain_id][0].get_item(res_nr)[1].insert(atom_nr, (atom_name, np.array([x, y, z])))

                        if chain_id not in b_factors:
                            b_factors[chain_id] = {}
                        if res_nr not in b_factors[chain_id]:
                            b_factors[chain_id][res_nr] = []
                        b_factors[chain_id][res_nr].append(b_factor)

                        if x < box_map[chain_id][0]:
                            box_map[chain_id][0] = x
                        if x > box_map[chain_id][1]:
                            box_map[chain_id][1] = x
                        if y < box_map[chain_id][2]:
                            box_map[chain_id][2] = y
                        if y > box_map[chain_id][3]:
                            box_map[chain_id][3] = y
                        if z < box_map[chain_id][4]:
                            box_map[chain_id][4] = z
                        if z > box_map[chain_id][1]:
                            box_map[chain_id][4] = z

                        if chain_id not in res_contig_map:
                            res_contig_map[chain_id] = residue_package.Residue_Map()
                            res_contig_map[chain_id].add_item(res_nr, (1, res_name))
                            contig_help_map[chain_id] = 1
                        elif not res_contig_map[chain_id].contains(res_nr):
                            contig_help_map[chain_id] += 1
                            res_contig_map[chain_id].add_item(res_nr, (contig_help_map[chain_id], res_name))
                    else:
                        if not coordinate_map[chain_id][1].contains(res_nr):  # If not, then add it to the ligands
                            coordinate_map[chain_id][1].add_item(res_nr, (res_name, sdsc_utils.SparseArray()))
                        coordinate_map[chain_id][1].get_item(res_nr)[1].insert(atom_nr, (atom_name, np.array([x, y, z])))

                        if res_name in METAL_ATOMS:
                            metals.add((chain_id, res_nr))
                        elif res_name in ION_ATOMS:
                            ions.add((chain_id, res_nr))
                        else:
                            ligands.add((chain_id, res_nr))

    # New Peptide detection counts irregular amino acids
    for chain_id in peptide_count:
        if chain_type_map[chain_id] == 'DNA' or chain_type_map[chain_id] == 'RNA':
            continue
        if peptide_count[chain_id][0] <= peptide_count[chain_id][1] * 2:
            chain_type_map[chain_id] = 'Peptide'
        elif (peptide_count[chain_id][0] + peptide_count[chain_id][1]) < 500:  # Total number of atoms
            chain_type_map[chain_id] = 'Peptide'

    return (coordinate_map, res_contig_map, ligands, metals, ions, box_map,
            chain_type_map, chainlist, b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map)

def get_atomic_coverage(coordinate_map: dict[str, tuple[residue_package.Residue_Map, residue_package.Residue_Map]], chain: str) -> float:
    expected_number_of_atoms = 0
    resolved_number_of_atoms = 0
    for res_nr in coordinate_map[chain][0].get_keys():
        res_name: str
        atoms: sdsc_utils.SparseArray
        res_name, atoms = coordinate_map[chain][0].get_item(res_nr)
        try:
            one_letter = THREE_TO_ONE[res_name]
        except:
            return f'{res_name} not in THREE_TO_ONE'
        
        if one_letter == 'X':
            expected_number_of_atoms += 1
            continue

        if one_letter not in CORRECT_COUNT:
            expected_number_of_atoms += 1
            continue

        expected_number_of_atoms += CORRECT_COUNT[one_letter]
        resolved_number_of_atoms += len(atoms)
    if expected_number_of_atoms == 0:
        return 0
    return (resolved_number_of_atoms / expected_number_of_atoms)

def getMinSubDist(
        coordinate_map: dict[str, tuple[residue_package.Residue_Map, residue_package.Residue_Map]],
        fuzzy_dist_matrix: dict[str, dict[str, residue_package.Residue_Map]],
        target_res_id, res_chain, chain, distance_threshold=10.0):
    top20: list[tuple[int | str, float]] = []
    top: int = 20
    for res in coordinate_map[chain][0].get_keys():
        if chain == res_chain and res == target_res_id:
            continue
        try:
            d: float = fuzzy_dist_matrix[res_chain][chain].get_item(target_res_id).get_item(res)
            if d is None:
                continue
        except:
            continue

        if len(top20) < top:
            top20.append((res, d))
        elif len(top20) == top:
            top20.append((res, d))
            top20 = sorted(top20, key=lambda x: x[1])
        elif d < top20[top][1]:
            top20[top] = (res, d)
            top20 = sorted(top20, key=lambda x: x[1])

    min_sub_d = None
    min_res = None
    min_atom_sub = None
    min_atom_chain = None
    atom_tuple = coordinate_map[res_chain][0].get_item(target_res_id)
    if atom_tuple is None:
        return min_sub_d, min_res, min_atom_sub, min_atom_chain, None
    atomlist: sdsc_utils.SparseArray = atom_tuple[1]
    inside_sphere = {}
    for (res, d) in top20:
        for atomnr in coordinate_map[chain][0].get_item(res)[1].get_keys():
            # exclude Hydrogens
            if coordinate_map[chain][0].get_item(res)[1].get(atomnr)[0][0] != 'H':
                coord1 = coordinate_map[chain][0].get_item(res)[1].get(atomnr)[1]
                for atomnr2 in atomlist.get_keys():
                    if atomlist.get(atomnr2)[0][0] != 'H':
                        coord2 = atomlist.get(atomnr2)[1]
                        d = distance(coord1, coord2)
                        if d < 1.2:
                            continue
                        if d < distance_threshold:
                            if res not in inside_sphere:
                                inside_sphere[res] = d
                            elif d < inside_sphere[res]:
                                inside_sphere[res] = d
                        if min_sub_d is None or d < min_sub_d:
                            min_sub_d = d
                            min_res = res
                            min_atom_sub = atomnr2
                            min_atom_chain = atomnr

    return min_sub_d, min_res, min_atom_sub, min_atom_chain, inside_sphere

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
def distance_pre_check(coord_1: np.ndarray, coord_2: np.ndarray, distance: float) -> bool:
    if abs(coord_1[0] - coord_2[0]) > distance:
        return True
    if abs(coord_1[1] - coord_2[1]) > distance:
        return True
    if abs(coord_1[2] - coord_2[2]) > distance:
        return True    
    return False

def calcFuzzyDM(
        coordinate_map: dict[str, tuple[residue_package.Residue_Map, residue_package.Residue_Map]],
        box_map, config,
        distance_threshold: float,
        target_chains=None
        ) -> dict[str, dict[str, residue_package.Residue_Map]]:
    
    fuzzy_dm: dict[str, dict[str, residue_package.Residue_Map]] = {}
    processed_chains = set()
    for chain in coordinate_map:
        chain_1_atom_residue_map: residue_package.Residue_Map = coordinate_map[chain][0]
        residue_list_chain_1: list[int | str] = chain_1_atom_residue_map.get_keys()

        if chain not in fuzzy_dm:
            fuzzy_dm[chain] = {}

        for chain_2 in coordinate_map:
            if chain_2 in processed_chains:
                continue
            if target_chains is not None:
                if not (chain in target_chains or chain_2 in target_chains):
                    continue
            neighbors = box_check(box_map[chain], box_map[chain_2], distance_threshold=config.short_distance_threshold)
            if not neighbors:
                continue

            if chain_2 not in fuzzy_dm[chain]:
                fuzzy_dm[chain][chain_2] = residue_package.Residue_Map()

            if chain_2 not in fuzzy_dm:
                fuzzy_dm[chain_2] = {}
            if chain not in fuzzy_dm[chain_2]:
                fuzzy_dm[chain_2][chain] = residue_package.Residue_Map()

            chain_2_atom_residue_map: residue_package.Residue_Map = coordinate_map[chain_2][0]
            residue_list_chain_2: list[int | str] = chain_2_atom_residue_map.get_keys()

            for res in residue_list_chain_1:
                test_coord: np.ndarray = chain_1_atom_residue_map.get_item(res)[1].get_first_item()[1]
                

                for res_2 in residue_list_chain_2:
                    if chain == chain_2 and res == res_2:
                        continue
                    test_coord_2: np.ndarray = chain_2_atom_residue_map.get_item(res_2)[1].get_first_item()[1]

                    if distance_pre_check(test_coord, test_coord_2, distance_threshold):
                        continue

                    d = distance(test_coord, test_coord_2)
                    if d > distance_threshold:
                        continue
                    
                    d = getMinDist(coordinate_map, res, chain, res_2, chain_2)

                    if not fuzzy_dm[chain][chain_2].contains(res):
                        fuzzy_dm[chain][chain_2].add_item(res, residue_package.Residue_Map())
                    fuzzy_dm[chain][chain_2].get_item(res).add_item(res_2, d)

                    if not fuzzy_dm[chain_2][chain].contains(res_2):
                        fuzzy_dm[chain_2][chain].add_item(res_2, residue_package.Residue_Map())

                    fuzzy_dm[chain_2][chain].get_item(res_2).add_item(res, d)

        processed_chains.add(chain)
    return fuzzy_dm

@njit
def calc_min_dist_between_atomic_coordinates_w(
    atomic_coordinates_res_1: np.ndarray,
    atomic_coordinates_res_2: np.ndarray
    ) -> tuple[float, int, int]:

    min_d: float = float('inf')
    for p1, coord_1 in enumerate(atomic_coordinates_res_1):
        for p2, coord_2 in enumerate(atomic_coordinates_res_2):
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
        for coord_2 in atomic_coordinates_res_2:
            d: float = distance(coord_1, coord_2)
            if d < min_d:
                min_d = d

    return min_d


def getMinDist(
        coordinate_map: dict[str, tuple[residue_package.Residue_Map, residue_package.Residue_Map]],
        res: int | str, chain: str,
        res2: int | str, chain2: str,
        het: int = 0, het2: int = 0,
        return_atom_names: bool = False,
        ) -> None | tuple[float, str, str] | float:
    
    
    if not coordinate_map[chain][het].contains(res):
        return None
    if not coordinate_map[chain2][het2].contains(res2):
        return None

    atom_list_res_1: sdsc_utils.SparseArray = coordinate_map[chain][het].get_item(res)[1]
    atom_nrs: list[int] =  atom_list_res_1.get_keys()
    atomic_coordinates_res_1: np.ndarray = np.zeros((len(atom_nrs), 3))
    for pos, atomnr in enumerate(atom_nrs):
        atomic_coordinates_res_1[pos] = atom_list_res_1.get(atomnr)[1]

    atom_list_res_2: sdsc_utils.SparseArray = coordinate_map[chain2][het2].get_item(res2)[1]
    atom_nrs_2: list[int] =  atom_list_res_2.get_keys()
    atomic_coordinates_res_2: np.ndarray = np.zeros((len(atom_nrs_2), 3))
    for pos, atomnr in enumerate(atom_nrs_2):
        atomic_coordinates_res_2[pos] = atom_list_res_2.get(atomnr)[1]    

    if return_atom_names:
        min_d, min_atom_p_1, min_atom_p_2 = calc_min_dist_between_atomic_coordinates_w(atomic_coordinates_res_1, atomic_coordinates_res_2)

        min_atom: str = atom_list_res_1.get(atom_nrs[min_atom_p_1])[0]
        min_atom2: str = atom_list_res_2.get(atom_nrs_2[min_atom_p_2])[0]

        return min_d, min_atom, min_atom2

    else:
        min_d = calc_min_dist_between_atomic_coordinates(atomic_coordinates_res_1, atomic_coordinates_res_2)
        return min_d

def calculate_interfaces(IAmap, dssp_dict, chain_type_map, config):
    interfaces = {}

    for chain, res, chain_b, res_b, btype, ctype, score in IAmap:
        if btype != 'combi':
            continue
        if ctype != 'all':
            continue
        if chain not in chain_type_map:
            continue
        if chain_type_map[chain] != 'Protein':
            continue
        
        if chain_b not in chain_type_map:
            continue
        if chain_type_map[chain_b] != 'Protein': #At the moment, we are only interested in PPIs
            continue
        if chain == chain_b:
            continue
        if not chain in interfaces:
            interfaces[chain] = {}
        if not chain_b in interfaces[chain]:
            interfaces[chain][chain_b] = interface_package.Interface(chain, chain_b)
        
        interfaces[chain][chain_b].add_interaction(res, res_b, score)
        #print('Add interaction:', chain, chain_b, res, res_b)

    """
    #find interface edge triangles
    for chain in interfaces:
        for chain_b in interfaces[chain]:
            for res in interfaces[chain][chain_b].interactions:
                for chain_c in IAmap[chain].get_item(res):
                    if chain_c != chain:
                        continue
                    if not chain_c in dssp_dict:
                        continue
                    for res_c in IAmap[chain].get_item(res)[chain_c].get_keys():
                        if res_c in interfaces[chain][chain_b].interactions:
                            continue
                        if not dssp_dict[chain_c].contains(res_c):
                            continue
                        sc_rsa = dssp_dict[chain_c].get_item(res_c)[2]
                        if sdsc_utils.locate(sc_rsa, config, binary_decision=True) != 'Surface':
                            continue
                        edge_count = 0
                        if not IAmap[chain_c].contains(res_c):
                            continue
                        for chain_d in IAmap[chain_c].get_item(res_c):
                            if chain_d != chain:
                                continue
                            for res_d in IAmap[chain_c].get_item(res_c)[chain_d].get_keys():
                                if res_d in interfaces[chain][chain_b].interactions:
                                    edge_count += 1
                        if edge_count >= 2:
                            interfaces[chain][chain_b].add_support(res_c)
    """
    return interfaces

@ray.remote(max_calls = 1)
def analysis_chain_remote_wrapper(target_chain, analysis_dump):

    config, profiles, centralities, packed_analysis_dump = analysis_dump

    result = [analysis_chain(target_chain, config, profiles, centralities, packed_analysis_dump)]

    #print('Finished nested chain:', target_chain)

    result = pack(result)

    #print('Finished packaging of nested chain:', target_chain)

    return result

@ray.remote(max_calls = 1)
def analysis_chain_package_wrapper(package, analysis_dump):
    results = []
    config, profiles, centralities, packed_analysis_dump = analysis_dump
    analysis_dump = packed_analysis_dump
    for target_chain in package:
        results.append(analysis_chain(target_chain, config, profiles, centralities, analysis_dump))
    return pack(results)

@ray.remote(max_calls = 1)
def annotate_wrapper(conf_dump, chunk, structure_ids):
    config = conf_dump
    return annotate(config, chunk, structure_ids, locked = True)

def annotate(config, chunk, structure_ids, locked = False):

    outputs = []
    nested_outputs = []
    nested_processes_dump = []
    total_times = []
    total_nested_times = []
    total_serialized_times = []
    structural_analysis_list: list[tuple[str, dict[str, residue_package.Residue_Map]]] = []
    for pdb_id, target_dict, nested_cores_limiter, nested_call in chunk:
        t0 = time.time()
        model_path = None
        # set model path to pdb file path in custom db path 
        if config.custom_db_path:
            custom_path_list = os.listdir(config.custom_db_path)
            for model in custom_path_list:
                if model.startswith(pdb_id):
                    model_path = config.custom_db_path + "/%s" %model

        elif config.model_db_active:
            if is_alphafold_model(pdb_id):
                model_path = alphafold_model_id_to_file_path(pdb_id, config)
        (structural_analysis_dict, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles,
         chain_type_map, chainlist, nested_processes, interfaces, times, nested_times, serialized_times) = structuralAnalysis(pdb_id, config, target_dict=target_dict, nested_cores_limiter = nested_cores_limiter, nested_call = nested_call, model_path = model_path, chain_ids = structure_ids[pdb_id])
        t1 = time.time()

        total_times = aggregate_times(total_times, times)
        total_nested_times = aggregate_times(total_nested_times, nested_times)
        total_serialized_times = aggregate_times(total_serialized_times, serialized_times)

        if config.verbosity >= 5:
            print(f'Time for structural analysis of {pdb_id}: {t1 - t0}, nested call: {nested_call}')

        if nested_call:
            nested_outputs.append((pdb_id, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, interfaces, t1 - t0))
            nested_processes_dump.append(nested_processes)
        else:
            outputs.append((pdb_id, structural_analysis_dict, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, interfaces, t1 - t0))
            structural_analysis_list.append((pdb_id, structural_analysis_dict))

    outputs = pack(outputs)
    nested_outputs = pack(nested_outputs)

    if not nested_call:
         background_process = remote_insertResidues(structural_analysis_list, structure_ids, config, locked=locked)

    return (outputs, nested_outputs, nested_processes_dump, total_times, total_nested_times, total_serialized_times, background_process)

def analysis_chain(
        target_chain: str,
        config,
        profiles: dict[str, residue_package.Residue_Map[rin.Interaction_profile]],
        centralities, analysis_dump):
    coordinate_map: dict[str, tuple[residue_package.Residue_Map, residue_package.Residue_Map]]
    res_contig_map: dict[str, residue_package.Residue_Map[tuple[int, str]]]
    target_residues: dict[str, residue_package.Residue_Map[tuple[int, str]]]
    chain_type_map: dict[str, str]

    (pdb_id, target_residues,
     res_contig_map, coordinate_map, fuzzy_dist_matrix, chain_type_map,
     b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map,
     dssp, dssp_dict, page_altered) = analysis_dump
    errorlist = []

    if config.verbosity >= 5:
        if len(target_residues[target_chain]) <= 1:
            print(f'Call of analysis_chain: {pdb_id} {target_chain} {target_residues[target_chain]}')
        else:
            print(f'Call of analysis_chain: {pdb_id} {target_chain} {len(target_residues[target_chain])}')

    structural_analysis_dict = residue_package.Residue_Map()

    sub_loop_broken = False

    for target_res_id in target_residues[target_chain].get_keys():
        if not res_contig_map[target_chain].contains(target_res_id):
            #print(f'{pdb_id} {target_chain} {target_res_id}')
            continue
        if sub_loop_broken:
            break
        three_letter = res_contig_map[target_chain].get_item(target_res_id)[1]
        if three_letter in THREE_TO_ONE:
            if THREE_TO_ONE[three_letter][0] in ONE_TO_THREE:
                one_letter = THREE_TO_ONE[three_letter][0]
            else:
                one_letter = 'X'
        else:
            one_letter = 'X'

        residue: residue_package.Residue = residue_package.Residue()
        residue.res_num = target_res_id
        residue.aa = one_letter

        lig_dists = {}
        min_chain_dists = {}

        inter_chain_kds = []
        inter_chain_dist_weighted_kds = []
        inter_chain_rsas = []
        inter_chain_dist_weighted_rsas = []
        total_inter_dist_weights = 0.

        inter_chain_interactions = []
        inter_chain_interactions_dw = []

        # Distance Calculations
        for chain in coordinate_map:

            #Sub - Chain - Calculations
            if chain != target_chain and len(coordinate_map[chain][0]) > 0:
                min_chain_dist, min_res, atom_sub, atom_chain, inside_sphere = getMinSubDist(coordinate_map, fuzzy_dist_matrix, target_res_id,
                                                                                             target_chain, chain, distance_threshold=config.milieu_threshold)
                min_chain_dists[chain] = (min_chain_dist, (atom_chain, atom_sub), min_res)
                if chain_type_map[chain] == 'Protein':  # Filter DNA and RNA chains
                    if inside_sphere is not None:
                        for inter_chain_res in inside_sphere:
                            if chain not in dssp_dict:
                                continue
                            if not dssp_dict[chain].contains(inter_chain_res):
                                continue
                            dist = inside_sphere[inter_chain_res]
                            total_inter_dist_weights += 1 / dist
                            if not res_contig_map[chain].contains(inter_chain_res):
                                continue
                            if not res_contig_map[chain].get_item(inter_chain_res)[1] in THREE_TO_ONE:
                                continue
                            inter_chain_res_one_letter = THREE_TO_ONE[res_contig_map[chain].get_item(inter_chain_res)[1]][0]
                            kd = HYDROPATHY[inter_chain_res_one_letter]
                            inter_chain_kds.append(kd)
                            inter_chain_dist_weighted_kds.append(kd / dist)
                            rsa = dssp_dict[chain].get_item(inter_chain_res)[0]
                            if rsa is not None:
                                inter_chain_rsas.append(rsa)
                                inter_chain_dist_weighted_rsas.append(rsa / dist)
                            if chain not in profiles:
                                errorlist.append('(1) target_chain not in profiles: %s %s %s' % (pdb_id, chain, list(profiles.keys())))
                                sub_loop_broken = True
                                break
                            elif profiles[chain].contains(inter_chain_res):
                                interface_score = profiles[chain].get_item(inter_chain_res).getTotalInterfaceInteractionScore()
                                inter_chain_interactions.append(interface_score)
                                inter_chain_interactions_dw.append(interface_score / dist)

            # Residue-Residue Calculations
            elif chain == target_chain and len(coordinate_map[chain][0]) > 0:
                if chain_type_map[chain] == 'Protein':  # Filter DNA and RNA chains
                    if chain not in dssp_dict:
                        continue
                    min_chain_dist, min_res, atom_sub, atom_chain, inside_sphere = getMinSubDist(coordinate_map, fuzzy_dist_matrix,
                                                                                                 target_res_id, target_chain, chain,
                                                                                                 distance_threshold=config.milieu_threshold)
                    intra_chain_kds = []
                    intra_chain_dist_weighted_kds = []
                    intra_chain_rsas = []
                    intra_chain_dist_weighted_rsas = []
                    total_dist_weights = 0.


                    if inside_sphere is not None:
                        if config.verbosity >= 8:
                            print('Residue-Residue calc:', pdb_id, chain, target_res_id, 'inside the sphere:', len(inside_sphere))

                        for intra_chain_res in inside_sphere:
                            if not dssp_dict[chain].contains(intra_chain_res):
                                continue
                            if not res_contig_map[chain].get_item(intra_chain_res)[1] in THREE_TO_ONE:
                                continue
                            dist = inside_sphere[intra_chain_res]
                            total_dist_weights += 1 / dist
                            intra_chain_res_one_letter = THREE_TO_ONE[res_contig_map[chain].get_item(intra_chain_res)[1]][0]
                            kd = HYDROPATHY[intra_chain_res_one_letter]
                            intra_chain_kds.append(kd)
                            intra_chain_dist_weighted_kds.append(kd / dist)
                            rsa = dssp_dict[chain].get_item(intra_chain_res)[0]
                            if rsa is not None:
                                intra_chain_rsas.append(rsa)
                                intra_chain_dist_weighted_rsas.append(rsa / dist)

                    if total_dist_weights > 0:
                        residue.intra_chain_median_kd = median(intra_chain_kds)
                        residue.intra_chain_dist_weighted_kd = sum(intra_chain_dist_weighted_kds) / total_dist_weights
                        residue.intra_chain_median_rsa = median(intra_chain_rsas)
                        residue.intra_chain_dist_weighted_rsa = sum(intra_chain_dist_weighted_rsas) / total_dist_weights

                    min_chain_dist, min_res, atom_sub, atom_chain, inside_sphere = getMinSubDist(coordinate_map, fuzzy_dist_matrix,
                                                                                                 target_res_id, target_chain, chain,
                                                                                                 distance_threshold=config.intra_milieu_threshold)

                    intra_chain_interactions = []
                    intra_chain_interactions_dw = []
                    total_dist_weights = 0.
                    if inside_sphere is not None:
                        for intra_chain_res in inside_sphere:
                            dist = inside_sphere[intra_chain_res]
                            total_dist_weights += 1 / dist
                            if not chain in profiles:
                                errorlist.append('(2) target_chain not in profiles: %s %s %s' % (pdb_id, chain, list(profiles.keys())))
                                sub_loop_broken = True
                                break
                            elif profiles[chain].contains(intra_chain_res):
                                interface_score = profiles[chain].get_item(intra_chain_res).getTotalInterfaceInteractionScore()
                                intra_chain_interactions.append(interface_score)
                                intra_chain_interactions_dw.append(interface_score / dist)

                    if total_dist_weights > 0:
                        residue.intra_chain_interactions_median = median(intra_chain_interactions)
                        residue.intra_chain_interactions_dist_weighted = sum(intra_chain_interactions_dw) / total_dist_weights

            #Sub - Lig - Calculations
            for hetres in coordinate_map[chain][1].get_keys():
                gmd_out = getMinDist(coordinate_map, target_res_id, target_chain, hetres, chain, het2=1, return_atom_names=True)
                if gmd_out is None:
                    continue
                min_d, atom, atom2 = gmd_out
                # only return ligand distances that are inside a given threshold
                if min_d > config.ligand_interest_sphere:
                    continue
                abr = coordinate_map[chain][1].get_item(hetres)[0]
                ligand_identifier = "%s_%s_%s" % (abr, hetres, chain)
                lig_dists[ligand_identifier] = (min_d, (atom, atom2))

        residue.lig_dists = lig_dists
        residue.chain_distances = min_chain_dists

        if total_inter_dist_weights > 0:
            residue.inter_chain_median_kd = median(inter_chain_kds)
            residue.inter_chain_dist_weighted_kd = sum(inter_chain_dist_weighted_kds) / total_inter_dist_weights
            residue.inter_chain_median_rsa = median(inter_chain_rsas)
            residue.inter_chain_dist_weighted_rsa = sum(inter_chain_dist_weighted_rsas) / total_inter_dist_weights
            residue.inter_chain_interactions_median = median(inter_chain_interactions)
            residue.inter_chain_interactions_dist_weighted = sum(inter_chain_interactions_dw) / total_inter_dist_weights

        if dssp:
            if target_chain in dssp_dict:
                if dssp_dict[target_chain].contains(target_res_id):
                    (residue.RSA, residue.relative_main_chain_acc, residue.relative_side_chain_acc, residue.SSA, residue.phi, residue.psi) = dssp_dict[target_chain].get_item(target_res_id)

            else:
                atomic_coverage = get_atomic_coverage(coordinate_map, chain)
                if isinstance(atomic_coverage, str):
                    errorlist.append(f'Error in get_atomic_coverage: {pdb_id}, {target_chain}\n{atomic_coverage}')
                elif atomic_coverage > 0.5 and not page_altered: #DSSP does not work for chains with missing atomic coordinates or chains over the 10K atoms mark, no need for a warning here
                    errorlist.append(("dssp error: chain not in dssp_dict; %s; %s" % (pdb_id, target_chain)))
                dssp = False


        #print(f'Before Adding profiles: {target_chain} {target_res_id} {list(profiles.keys())}')
        if target_chain in profiles:
            #print(profiles[target_chain].get_keys())
            if profiles[target_chain].contains(target_res_id):
                #print(f'Setting {profiles[target_chain].get_item(target_res_id)=} to {target_chain} {target_res_id}')
                residue.interaction_profile = profiles[target_chain].get_item(target_res_id)

        if target_chain in centralities:
            if centralities[target_chain].contains(target_res_id):
                residue.centralities = centralities[target_chain].get_item(target_res_id)

        avg_b_factor = sum(b_factors[target_chain][target_res_id]) / float(len(b_factors[target_chain][target_res_id]))

        if target_chain in modres_map:
            if target_res_id in modres_map[target_chain]:
                residue.modres = True
            else:
                residue.modres = False
        else:
            residue.modres = False

        if (target_chain, target_res_id) in ssbond_map:
            (chain_2, res_nr_2, residue.ssbond_length) = ssbond_map[(target_chain, target_res_id)]
            if target_chain == chain_2:
                residue.intra_ssbond = True
                residue.inter_ssbond = False
            else:
                residue.intra_ssbond = False
                residue.inter_ssbond = True
        else:
            residue.intra_ssbond = False
            residue.inter_ssbond = False

        if (target_chain, target_res_id) in link_map:
            (atom_1, res_name_1, atom_2, res_name_2, res_nr_2, chain_2, residue.link_dist) = link_map[(target_chain, target_res_id)]
            if target_chain == chain_2:
                residue.intra_link = True
                residue.inter_link = False
            else:
                residue.intra_link = False
                residue.inter_link = True
        else:
            residue.intra_link = False
            residue.inter_link = False

        if (target_chain, target_res_id) in cis_conformation_map:
            (res_name_2, chain_2, res_nr_2, angle) = cis_conformation_map[(target_chain, target_res_id)]
            residue.cis_conformation = angle

        if (target_chain, target_res_id) in cis_follower_map:
            (res_name_2, chain_2, res_nr_2, angle) = cis_follower_map[(target_chain, target_res_id)]
            residue.cis_follower = angle

        #print(f'{pdb_id} {is_alphafold_model(pdb_id)} {avg_b_factor}')
        if is_alphafold_model(pdb_id):
            residue.pLDDT = avg_b_factor
        else:
            residue.b_factor = avg_b_factor

        structural_analysis_dict.add_item(target_res_id, residue)

    if config.verbosity >= 5:
        if len(structural_analysis_dict) <= 10:
            print(f'Results of analysis_chain: {pdb_id} {target_chain} {structural_analysis_dict}')
        else:
            print(f'Results of analysis_chain: {pdb_id} {target_chain} {len(structural_analysis_dict)}')

    return target_chain, structural_analysis_dict, errorlist


def structuralAnalysis(
        pdb_id: str, config: any,
        target_dict: list[str] | None = None,
        model_path: str = None, keep_rin_files: bool =True,
        nested_cores_limiter: int | None = None, nested_call: bool = False,
        chain_ids: dict[str, int] = None):

    dssp: bool = True
    dssp_path: str = config.dssp_path
    pdb_path: str = config.pdb_path
    rin_db_path: str = config.rin_db_path
    verbosity: int = config.verbosity

    times = []
    t0 = time.time()

    if verbosity >= 6:  
        print('Start structuralAnalysis of:', pdb_id)

    page, page_altered, path, _ = pdbParser.standardParsePDB(pdb_id, pdb_path, return_10k_bool=True, get_is_local=True, model_path=model_path)

    t1 = time.time()
    times.append(t1-t0)

    if verbosity >= 8:
        print('\n',pdb_id,'\n\n',page)

    errorlist = []

    if page == '':
        errorlist.append("Error while parsing: %s" % pdb_id)
        return {}, errorlist, {}, {}, {}, {}, {}, {}, [], {}, None

    (coordinate_map, res_contig_map, ligands, metals, ions, box_map, chain_type_map, chainlist,
        b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map) = parsePDB(page)

    t2 = time.time()
    times.append(t2-t1)

    if target_dict is None:
        target_residues: dict[str, residue_package.Residue_Map[tuple[int, str]]] = {}
        for chain in res_contig_map:
            if chain_type_map[chain] != 'Protein':
                continue
            if chain not in chain_ids:
                continue
            target_residues[chain] = res_contig_map[chain]
    elif isinstance(target_dict, list):
        target_residues = {}
        for chain in target_dict:
            target_residues[chain] = res_contig_map[chain]

    t3 = time.time()
    times.append(t3-t2)

    #distance_threshold = max([2. * config.short_distance_threshold, 2. * config.milieu_threshold])
    distance_threshold = 2. * config.short_distance_threshold

    fuzzy_dist_matrix = calcFuzzyDM(coordinate_map, box_map, config, distance_threshold, target_chains=set(target_residues.keys()))

    t4 = time.time()
    times.append(t4-t3)

    structural_analysis_dict = {}

    if dssp:
        # write temp pdb file only if we had to fix the 10k atom bug or the file is not stored locally
        if page_altered or path is None:
            tmp_path = f'{config.temp_folder}/tmp_{pdb_id}.pdb'
            f = open(tmp_path, 'w')
            f.write(page)
            f.close()
        else:
            tmp_path = path

        # call DSSP  -structure of dssp_dict: {Chain:{Res_id:(racc,ssa)}}
        dssp_dict, errlist = calcDSSP(tmp_path, dssp_path, angles=True, verbosity_level=verbosity)
        errorlist += errlist
        # remove tmp_file
        if page_altered:
            try:
                os.remove(tmp_path)
            except:
                pass
        if dssp_dict == {}:
            dssp = False

    t5 = time.time()
    times.append(t5-t4)

    profiles: dict[str, residue_package.Residue_Map[rin.Interaction_profile]] = {}
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}
    IAmap = {}
    lookup_times = []

    lookup_out = rin.lookup(pdb_id, page, config, ligands,
                            metals, ions, res_contig_map,
                            rin_db_path, chain_type_map, model_path=model_path, keep_tmp_files=keep_rin_files)
    if isinstance(lookup_out, str):
        errorlist.append(lookup_out)
    else:
        (profiles, centralities, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, IAmap, lookup_times) = lookup_out
        if len(profiles) == 0:
            errorlist.append('Empty profiles %s %s %s' % (pdb_id, str(len(page)), str(len(res_contig_map))))
        elif verbosity >= 5:
            print('Chains in profiles:', list(profiles.keys()))


    times.append(lookup_times)

    t6 = time.time()
    times.append(t6-t5)

    if config.compute_ppi:
        interfaces = calculate_interfaces(IAmap, dssp_dict, chain_type_map, config)
    else:
        interfaces = {}

    t7 = time.time()
    times.append(t7-t6)

    number_protein_chains = 0
    for target_chain in target_residues:
        if not target_chain in chain_type_map:
            config.errorlog.add_warning(f'{target_chain} not in chain_type_map for: {pdb_id}')
            continue
        if chain_type_map[target_chain] != 'Protein':
            continue
        number_protein_chains += 1

    if nested_cores_limiter is None:
        nested_cores_limiter = config.proc_n - 1

    if nested_call and not config.low_mem_system:
        nested_times = []


        analysis_dump = ray.put((config, profiles, (pdb_id, target_residues,
                                 res_contig_map, coordinate_map, fuzzy_dist_matrix, chain_type_map,
                                 b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map,
                                 dssp, dssp_dict, page_altered)))

        t8 = time.time()
        nested_times.append(t8-t7)

        core_results = []

        nested_cores = min([config.proc_n - 1, nested_cores_limiter])

        if number_protein_chains < nested_cores:
            packaged = False
            for target_chain in target_residues:
                if target_chain not in chain_type_map:
                    continue
                if chain_type_map[target_chain] != 'Protein':
                    continue

                core_results.append(analysis_chain_remote_wrapper.remote(target_chain, analysis_dump))

        else:
            if verbosity >= 5:
                print('Packaging a big structure:',pdb_id,'with',number_protein_chains,'chains')
            packaged = True
            max_package_size = number_protein_chains // nested_cores
            if (number_protein_chains % nested_cores) != 0:
                max_package_size += 1

            package = []
            for target_chain in target_residues:
                if target_chain not in chain_type_map:
                    continue
                if chain_type_map[target_chain] != 'Protein':
                    continue
                package.append(target_chain)
                if len(package) == max_package_size:
                    core_results.append(analysis_chain_package_wrapper.remote(package, analysis_dump))
                    package = []

            if len(package) != 0:
                core_results.append(analysis_chain_package_wrapper.remote(package, analysis_dump))

        if verbosity >= 5:
            print('Reached nested remotes:', pdb_id, number_protein_chains)

        t9 = time.time()
        nested_times.append(t9-t8)

        del target_residues
        del res_contig_map
        del coordinate_map
        del fuzzy_dist_matrix
        del b_factors
        del modres_map
        del ssbond_map
        del link_map
        del cis_conformation_map
        del cis_follower_map
        del profiles
        del dssp
        del dssp_dict
        del page_altered

        return {}, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, core_results, interfaces, times, nested_times, None

    else:
        serialized_times = []

        analysis_dump = (pdb_id, target_residues,
                         res_contig_map, coordinate_map, fuzzy_dist_matrix, chain_type_map,
                         b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map,
                         dssp, dssp_dict, page_altered)
        for target_chain in target_residues:
            if target_chain not in chain_type_map:
                continue
            if chain_type_map[target_chain] != 'Protein':
                continue

            target_chain, chain_structural_analysis_dict, chain_errorlist = analysis_chain(target_chain, config, profiles, centralities, analysis_dump)
            structural_analysis_dict[target_chain] = chain_structural_analysis_dict
            errorlist += ['%s - %s' % (x, pdb_id) for x in chain_errorlist[:5]]
            if config.verbosity >= 5:
                print(f'Received structural analysis results (1) from: {pdb_id} {target_chain} {len(chain_structural_analysis_dict)} {chain_errorlist}')

        t8 = time.time()
        serialized_times.append(t8-t7)

    return structural_analysis_dict, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, [], interfaces, times, None, serialized_times

def trim_down_backmaps(config, proteins):
    removed_backmaps = 0
    for structure_id in proteins.complexes:
        chains_with_interfaces = set()
        for chain_a in proteins.complexes[structure_id].interfaces:
            if chain_a not in proteins.structures[structure_id]:
                continue
            chains_with_interfaces.add(chain_a)
            for chain_b in proteins.complexes[structure_id].interfaces[chain_a]:
                new_back_maps = {}
                for res in proteins.complexes[structure_id].interfaces[chain_a][chain_b].residues:
                    for prot_id in proteins.structures[structure_id][chain_a].backmaps:
                        if prot_id not in new_back_maps:
                            new_back_maps[prot_id] = residue_package.Residue_Map()
                        try:
                            new_back_maps[prot_id].add_item(res, proteins.structures[structure_id][chain_a].backmaps[prot_id].get_item(res))
                        except:
                            pass
            proteins.structures[structure_id][chain_a].backmaps = new_back_maps

        #print(structure_id, chains_with_interfaces)

        
        for chain in proteins.complexes[structure_id].chains:
            if chain in chains_with_interfaces:
                continue
            if structure_id in proteins.structures:
                if chain in proteins.structures[structure_id]:
                    removed_backmaps += len(proteins.structures[structure_id][chain].backmaps)
                    proteins.structures[structure_id][chain].backmaps = {}
                    

    if config.verbosity >= 4:
        print(f'Removed backmaps: {removed_backmaps}')
    return

def paraAnnotate(config, proteins, indel_analysis_follow_up=False):

    t0 = time.time()

    # new structure for paralell annotation: give each annotation process a pdb and dict of chain and structure information,
    # all chains without structure information shall go into the database as fully annotated structures, which are not aligned (yet)
    # these are unmapped chains, or protein interaction partners, their presence in the database is necesary for the computation of the structural neighborhood
    # a new structure-structure entry has to be created and has to be inserted into the database before the residue insertion
    # homooligomer information is not present for such cases, this information has to be updated, if the chain is mapped at a later iteration or run

    size_list: list[tuple[str, int]] = []

    total_structural_analysis = {}
    background_insert_residues_process = None
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
            print(f'PDB entry {pdb_id} added to size_list, {ac=}')
        n_of_comps += 1

    chunksize: int = total_cost // (8 * config.proc_n)

    if config.verbosity >= 2:
        print(f'Starting structural analysis with {n_of_comps} complexes. {n_of_stored_complexes} complexes are already stored. Chunksize: {chunksize}, n_of_chain_thresh: {config.n_of_chain_thresh}')
        print(f'Total amount in structure_list: {len(structure_list)}. Total cost: {total_cost}, Threads: {config.proc_n}')

    size_list = sorted(size_list, key=lambda x:x[1] ,reverse=True)

    anno_result_ids = []

    assigned_costs = {}
    temporarily_freed = {}
    if config.proc_n > 1:
        conf_dump = ray.put(config)
    else:
        conf_dump = None

    n_started = 0
    started_cost = 0
    current_chunk = []
    structure_ids: dict[str, dict[str, int]] = {}
    chunk_cost = 0

    nested_packages = set()

    input_procession_counter = 0
    for pdb_id, ac in size_list:
        input_procession_counter += 1
        cost = ac
        if config.proc_n < 2:
            cost = 0

        current_chunk.append((pdb_id, None, None, False))
        structure_ids[pdb_id] = {}
        for chain in proteins.structures[pdb_id]:
            structure_ids[pdb_id][chain] = proteins.structures[pdb_id][chain].database_id
        chunk_cost += cost
        if chunk_cost >= chunksize:
            if config.proc_n > 1:
                anno_result_ids.append(annotate_wrapper.remote(conf_dump, current_chunk, structure_ids))
                n_started += 1
                started_cost += cost
            else:
                anno_result_ids.append(annotate(config, current_chunk, structure_ids))
            current_chunk = []
            structure_ids = {}
            chunk_cost = 0
        
            """ old routine to started nested proceses
            if config.proc_n > 1:
                anno_result_ids.append(annotate_wrapper.remote(conf_dump, [(pdb_id, None, cost, True)]))
            else:
                anno_result_ids.append(annotate(config, [(pdb_id, None, None, False)]))
            package_cost += cost  # big structures with nested remotes cost the amount of nested calls
            nested_packages.add(pdb_id)
            if config.verbosity >= 4:
                print('started nested package:', pdb_id, cost, s)
            """
        if n_started >= config.proc_n:
            break

    if config.proc_n < 2 and len(current_chunk) > 0:
        anno_result_ids.append(annotate(config, current_chunk, structure_ids))
        current_chunk = []
        structure_ids = {}
        chunk_cost = 0

    if config.verbosity >= 2:
        t11 = time.time()
        print(f"Annotation Part 1.1: {t11 - t0} {n_started=} {input_procession_counter=}")

    #getStoredResidues(proteins, config, exclude_interacting_chains = not config.compute_ppi)  # has to be called after insertStructures
    getStoredResidues(proteins, config, exclude_interacting_chains=True)

    if config.verbosity >= 2:
        t1 = time.time()
        print("Annotation Part 1.2: %s" % (str(t1 - t11)))
        print("Annotation Part 1: %s" % (str(t1 - t0)))

    max_comp_time = 0.
    total_comp_time = 0.
    max_comp_time_structure = None
    struct_comp_times = []
    amount_of_structures = 0
    amount_of_chains_in_analysis_dict = 0

    total_computed = 0

    no_result_run = 1

    core_not_ready_dict = {}

    loop_counter = 0

    #gc.disable()

    returned_pdbs = []
    total_times = []
    total_nested_times = []
    total_serialized_times = []

    while True:
        loop_counter += 1
        chunksize: int = (total_cost-started_cost) // (8 * config.proc_n)  # dynamically adjust chunksize
        finished = 0
        del_list = []
        # this loop collects the results from the nested remotes
        wait_counter = 0
        for ret_pdb_id in core_not_ready_dict:
            wait_counter += 1
            wait_counter_2 = 0
            if ret_pdb_id not in total_structural_analysis:
                total_structural_analysis[ret_pdb_id] = {}
            while True:
                wait_counter_2 += 1
                core_ready, core_not_ready = ray.wait(core_not_ready_dict[ret_pdb_id], timeout=0.1)
                if config.verbosity >= 5:
                    print(f'{wait_counter} {wait_counter_2} Nested wait:', ret_pdb_id, len(core_ready), len(core_not_ready))
                    print(f'RAM memory % used: {psutil.virtual_memory()[2]}')
                if len(core_ready) > 0:
                    errorlist = []
                    structural_analysis_dict = {}
                    core_outs = ray.get(core_ready)

                    for core_package in core_outs:
                        core_package = unpack(core_package)
                        for target_chain, chain_structural_analysis_dict, chain_errorlist in core_package:
                            structural_analysis_dict[target_chain] = chain_structural_analysis_dict
                            errorlist += chain_errorlist
                            package_cost -= 1  # relieve cost for each finished nested remote (individual chains)
                            temporarily_freed[ret_pdb_id] += 1
                            if config.verbosity >= 5:
                                print(f'Received structural analysis results (2) from: {ret_pdb_id} {target_chain} {len(chain_structural_analysis_dict)} {chain_errorlist}')
                    if len(errorlist) > 0:
                        for error_text in errorlist:
                            config.errorlog.add_warning(error_text)
                    for chain in structural_analysis_dict:
                        
                        total_structural_analysis[ret_pdb_id][chain] = structural_analysis_dict[chain]
                        amount_of_chains_in_analysis_dict += 1

                    #for core_ray_id in core_ready:
                    #    ray.cancel(core_ray_id, force=True)
                if len(core_not_ready) == 0:
                    del_list.append(ret_pdb_id)
                    finished += 1
                    break
                else:
                    core_not_ready_dict[ret_pdb_id] = core_not_ready
                if len(core_ready) == 0:
                    break


        for ret_pdb_id in del_list:
            del core_not_ready_dict[ret_pdb_id]
            nested_packages.remove(ret_pdb_id)
            package_cost = (package_cost - assigned_costs[ret_pdb_id]) + temporarily_freed[ret_pdb_id] # correct package cost, when the number of protein chains got estimated incorrectly
        
        not_ready = []
        if len(anno_result_ids) > 0:
            if config.verbosity >= 5:
                print('Before wait, loop counter:', loop_counter)

            '''
            ready, not_ready = ray.wait(anno_result_ids, timeout = 0.00001)
            if len(ready) == 0:
                t_gc_0 = time.time()
                gc.collect()
                t_gc_1 = time.time()
                if config.verbosity >= 5:
                    print('Manual gc happend:', (t_gc_1 - t_gc_0))
            '''
            if config.proc_n > 1:
                ready, not_ready = ray.wait(anno_result_ids)

                if config.verbosity >= 5:
                    print('Wait happened:', no_result_run)

                t_0 = time.time()

                if len(ready) > 0:
                    chunk_struct_outs = ray.get(ready)
                    no_result_run = 1
                else:
                    chunk_struct_outs = []
                    no_result_run = min([10, no_result_run + 1])

                t_1 = time.time()
                if config.verbosity >= 5:
                    print('Time for loop part 1:', (t_1-t_0))
            else:
                t_1 = time.time()
                chunk_struct_outs = anno_result_ids

            t_integrate_0 = 0.
            t_integrate_1 = 0.
            t_integrate_2 = 0.
            t_integrate_3 = 0.
            t_integrate_4 = 0.
            t_pickle_0 = 0.
            t_pickle_1 = 0.
            t_pickle_2 = 0.

            for chunk_struct_out, nested_struct_outs, nested_processes_list, times, nested_times, serialized_times, background_insert_residues_process in chunk_struct_outs:
                n_started -= 1  # finished chunks free 1 cost
                t_pickle_0 += time.time()
                chunk_struct_out = unpack(chunk_struct_out)
                nested_struct_outs = unpack(nested_struct_outs)
                t_pickle_1 += time.time()
                for struct_out in chunk_struct_out:
                    t_integrate_0 += time.time()
                    (ret_pdb_id, structural_analysis_dict, errorlist, ligand_profiles, metal_profiles,
                        ion_profiles, chain_chain_profiles, chain_type_map, chainlist, interfaces, comp_time) = struct_out

                    returned_pdbs.append(ret_pdb_id)

                    finished += 1

                    if comp_time > max_comp_time:
                        max_comp_time = comp_time
                        max_comp_time_structure = ret_pdb_id
                    if config.verbosity >= 3:
                        s = len(proteins.get_complex_chains(ret_pdb_id, only_protein=True))
                        ac = proteins.get_atom_count(ret_pdb_id)
                        struct_comp_times.append(f'{ret_pdb_id}\t{s}\t{ac}\t{comp_time}\n')
                    amount_of_structures += 1
                    total_comp_time += comp_time

                    total_times = aggregate_times(total_times, times)
                    total_nested_times = aggregate_times(total_nested_times, nested_times)
                    total_serialized_times = aggregate_times(total_serialized_times, serialized_times)

                    #This part adds the results from the structural analysis into the central datastructure
                    t_integrate_1 += time.time()
                    if len(errorlist) > 0:
                        for error_text in errorlist:
                            config.errorlog.add_warning(error_text)

                    if ret_pdb_id not in total_structural_analysis:
                        total_structural_analysis[ret_pdb_id] = {}

                    for chain in structural_analysis_dict:
                        if config.verbosity >= 5:
                            print(f'Received structural analysis results (3) from: {ret_pdb_id} {chain} {len(structural_analysis_dict[chain])}')
                        proteins.structures[ret_pdb_id][chain].residues = structural_analysis_dict[chain]
                        total_structural_analysis[ret_pdb_id][chain] = structural_analysis_dict[chain]
                        amount_of_chains_in_analysis_dict += 1

                    t_integrate_2 += time.time()

                    proteins.complexes[ret_pdb_id].lig_profile = ligand_profiles
                    proteins.complexes[ret_pdb_id].ion_profile = ion_profiles
                    proteins.complexes[ret_pdb_id].metal_profile = metal_profiles
                    proteins.complexes[ret_pdb_id].chain_chain_profile = chain_chain_profiles
                    #proteins.set_IAmap(ret_pdb_id, IAmap)
                    proteins.set_interfaces(ret_pdb_id, interfaces)

                    t_integrate_3 += time.time()

                    proteins.set_chain_type_map(ret_pdb_id, chain_type_map, chainlist)
                    del chain_type_map

                    if config.low_mem_system:
                        if amount_of_chains_in_analysis_dict > 2 * (config.chunksize):
                            if background_insert_residues_process is not None:
                                background_insert_residues_process.join()
                            
                            background_insert_residues_process = insertResidues(total_structural_analysis, proteins, config)
                            total_structural_analysis = {}
                            amount_of_chains_in_analysis_dict = 0
                            gc.collect()

                    t_integrate_4 += time.time()

                if len(nested_struct_outs) > 0:
                    config.errorlog.add_warning('Currently there should be no nested processes in structural analysis')
                for nested_process_number, struct_out in enumerate(nested_struct_outs):
                    package_cost += 1 #repay one cost, since this a nested out
                    t_integrate_0 += time.time()
                    (ret_pdb_id, errorlist, ligand_profiles, metal_profiles,
                        ion_profiles, chain_chain_profiles, chain_type_map, chainlist, interfaces, comp_time) = struct_out
                    
                    nested_processes = nested_processes_list[nested_process_number]
                    returned_pdbs.append(ret_pdb_id)

                    core_ready, core_not_ready = ray.wait(nested_processes, timeout=0.01)
                    structural_analysis_dict = {}
                    if len(core_ready) > 0:
                        core_outs = ray.get(core_ready)

                        for core_package in core_outs:
                            core_package = unpack(core_package)
                            for target_chain, chain_structural_analysis_dict, chain_errorlist in core_package:
                                structural_analysis_dict[target_chain] = chain_structural_analysis_dict
                                errorlist += chain_errorlist
                                package_cost -= 1  # relieve cost for each finished nested remote (individual chains)
                                temporarily_freed[ret_pdb_id] += 1
                                if config.verbosity >= 5:
                                    print(f'Received structural analysis results (4) from: {ret_pdb_id} {target_chain} {len(chain_structural_analysis_dict)} {chain_errorlist}')

                    if len(core_not_ready) > 0:
                        core_not_ready_dict[ret_pdb_id] = core_not_ready
                    else:
                        finished += 1
                        if ret_pdb_id in core_not_ready_dict:
                            del core_not_ready_dict[ret_pdb_id]
                        nested_packages.remove(ret_pdb_id)
                        package_cost = (package_cost - assigned_costs[ret_pdb_id]) + temporarily_freed[ret_pdb_id] # correct package cost, when the number of protein chains got estimated incorrectly

                    if comp_time > max_comp_time:
                        max_comp_time = comp_time
                        max_comp_time_structure = ret_pdb_id

                    if config.verbosity >= 3:
                        s = len(proteins.get_complex_chains(ret_pdb_id, only_protein=True))
                        ac = proteins.get_atom_count(ret_pdb_id)
                        struct_comp_times.append(f'{ret_pdb_id}\t{s}\t{ac}\t{comp_time}\n')

                    amount_of_structures += 1
                    total_comp_time += comp_time

                    total_times = aggregate_times(total_times, times)
                    total_nested_times = aggregate_times(total_nested_times, nested_times)
                    total_serialized_times = aggregate_times(total_serialized_times, serialized_times)

                    #This part adds the results from the structural analysis into the central datastructure
                    t_integrate_1 += time.time()
                    if len(errorlist) > 0:
                        for error_text in errorlist:
                            config.errorlog.add_warning(error_text)

                    if ret_pdb_id not in total_structural_analysis:
                        total_structural_analysis[ret_pdb_id] = {}
                    for chain in structural_analysis_dict:
                        proteins.structures[ret_pdb_id][chain].residues = structural_analysis_dict[chain]
                        total_structural_analysis[ret_pdb_id][chain] = structural_analysis_dict[chain]
                        amount_of_chains_in_analysis_dict += 1

                    t_integrate_2 += time.time()

                    proteins.set_lig_profile(ret_pdb_id, ligand_profiles)
                    proteins.set_ion_profile(ret_pdb_id, ion_profiles)
                    proteins.set_metal_profile(ret_pdb_id, metal_profiles)
                    proteins.set_chain_chain_profile(ret_pdb_id, chain_chain_profiles)
                    #proteins.set_IAmap(ret_pdb_id, IAmap)
                    proteins.set_interfaces(ret_pdb_id, interfaces)

                    

                    t_integrate_3 += time.time()

                    proteins.set_chain_type_map(ret_pdb_id, chain_type_map, chainlist)
                    del chain_type_map

                    if config.low_mem_system:
                        if amount_of_chains_in_analysis_dict > 2 * (config.chunksize):
                            if background_insert_residues_process is not None:
                                background_insert_residues_process.join()

                            background_insert_residues_process = insertResidues(total_structural_analysis, proteins, config)
                            total_structural_analysis = {}
                            amount_of_chains_in_analysis_dict = 0
                            gc.collect()
                    t_integrate_4 += time.time()
            t_2 = time.time()

            
            if config.verbosity >= 4:
                print('Time for result integration:', (t_integrate_4 - t_integrate_0), 'for amount of analyzed chunks:', len(chunk_struct_out),
                      'individual times:', (t_integrate_1 - t_integrate_0),
                                           (t_integrate_2 - t_integrate_1),
                                           (t_integrate_3 - t_integrate_2),
                                           (t_integrate_4 - t_integrate_3))
                print('Time for unpacking:', (t_pickle_1 - t_pickle_0))

                if (t_pickle_2 - t_pickle_1) > 60:
                    print('Long depickling:', returned_pdbs)
                print(f'RAM memory % used: {psutil.virtual_memory()[2]}')

            if config.verbosity >= 5:
                print('Time for loop part 2:', (t_2-t_1))

        new_anno_result_ids = []

        for pdb_id, ac in size_list[input_procession_counter:]:
            input_procession_counter += 1
            cost = ac
            if config.proc_n < 2:
                cost = 0
            current_chunk.append((pdb_id, None, None, False))
            structure_ids[pdb_id] = {}
            for chain in proteins.structures[pdb_id]:
                structure_ids[pdb_id][chain] = proteins.structures[pdb_id][chain].database_id
            chunk_cost += cost
            if chunk_cost >= chunksize:
                if config.proc_n > 1:
                    new_anno_result_ids.append(annotate_wrapper.remote(conf_dump, current_chunk, structure_ids))
                    n_started += 1
                    started_cost += cost
                else:
                    new_anno_result_ids.append(annotate(config, current_chunk, structure_ids))
                current_chunk = []
                structure_ids = {}
                chunk_cost = 0
            if n_started >= config.proc_n:
                break

        if (input_procession_counter >= len(size_list)) and (len(current_chunk) > 0):
            if config.proc_n > 1:
                new_anno_result_ids.append(annotate_wrapper.remote(conf_dump, current_chunk, structure_ids))
                n_started += 1
                started_cost += cost
            else:
                new_anno_result_ids.append(annotate(config, current_chunk, structure_ids))
            current_chunk = []
            structure_ids = {}
            chunk_cost = 0
        

        if len(new_anno_result_ids) > 0 and config.verbosity >= 5:
            t_3 = time.time()
            print('Time for loop part 3:', (t_3-t_2))

        #for ray_id in ready:
        #    ray.cancel(ray_id, force=True)

        anno_result_ids = new_anno_result_ids + not_ready

        if finished > 0:
            total_computed += finished
            if config.verbosity >= 5:
                print('Newly finished structural analysis of', finished, '. Total:', total_computed, 'Current package:', len(anno_result_ids), 'and nested packages:', str(nested_packages), ', number of pending nested packaged chains:', len(core_not_ready_dict), 'Amount of new started processes:', len(new_anno_result_ids), 'Amount of pending processes:', len(not_ready))
                print(f'RAM memory % used: {psutil.virtual_memory()[2]}')
        if len(anno_result_ids) == 0 and input_procession_counter >= len(size_list) and len(core_not_ready_dict) == 0:
            break


    del conf_dump
    t_1_5 = time.time()
    #gc.enable()

    if config.verbosity >= 3:
        print_times(total_times)
        print_times(total_nested_times)
        print_times(total_serialized_times)

    if config.verbosity >= 2:
        t2 = time.time()
        print('Time for reenabling gc:', (t2- t_1_5))
        print("Annotation Part 2: %s" % (str(t2 - t1)))
        print('Longest computation for:', max_comp_time_structure, 'with:', max_comp_time, 'In total', amount_of_structures, 'structures', 'Accumulated time:', total_comp_time)

    if config.verbosity >= 3:
        try:
            times_stat_file = f'{config.temp_folder}/../struct_wise_times.tsv'
            print(f'Writing struct wise times to: {times_stat_file}')
            with open(times_stat_file, 'w') as f:
                f.write(''.join(struct_comp_times))
        except:
            pass

    if config.verbosity >= 2:
        t32 = time.time()
        print(f'Annotation Part 3.1: {t32 - t2}')

    insertComplexes(proteins, config)

    if config.verbosity >= 2:
        t33 = time.time()
        print('Annotation Part 3.2:', t33 - t32)

    #if len(total_structural_analysis) > 0:
    #    background_insert_residues_process = insertResidues(total_structural_analysis, proteins, config)

    if config.verbosity >= 2:
        t34 = time.time()
        print('Annotation Part 3.3:', t34 - t33)

    trim_down_backmaps(config, proteins)

    if config.verbosity >= 2:
        t35 = time.time()
        print('Annotation Part 3.4:', t35 - t34)

    if config.verbosity >= 2:
        t3 = time.time()
        print("Annotation Part 3: %s" % (str(t3 - t2)))

    if not config.fast_pdb_annotation:
        classification(proteins, config, background_insert_residues_process)
    elif background_insert_residues_process is not None:
        if config.verbosity >= 2:
            t41 = time.time()
        background_insert_residues_process.join()
        background_insert_residues_process.close()
        background_insert_residues_process = None
        if config.verbosity >= 2:
            t42 = time.time()
            print(f'Annotation Part 4.1, pdb fast annotation, bg inserst residues: {t42-t41}')
        insert_interfaces(proteins, config)
        if config.verbosity >= 2:
            t43 = time.time()
            print(f'Annotation Part 4.2, pdb fast annotation, insert interfaces: {t43-t42}')

    if config.verbosity >= 2:
        t4 = time.time()
        print(f'Annotation Part 4: {t4-t3}')

    insert_interface_residues(proteins, config)

    if config.verbosity >= 2:
        t5 = time.time()
        print(f"Annotation Part 5: {t5-t4}")

    if not indel_analysis_follow_up:
        if config.verbosity >= 3:
            print('Proteins object semi deconstruction')
        #for protein_id in protein_ids:
        #    proteins.remove_protein_annotations(protein_id)
        proteins.semi_deconstruct()

    if not config.fast_pdb_annotation:
        insertClassifications(proteins, config)

    if config.verbosity >= 2:
        t6 = time.time()
        print(f"Annotation Part 6: {t6-t5}")

    gc.collect()

    if config.verbosity >= 2:
        t7 = time.time()
        print(f"Time for garbage collection: {t7-t6}")

    return amount_of_structures


# Return True, if the template is not good
def weakCriteria(seq, res, rel_aln_len, seq_thresh, res_thresh, rel_aln_len_thresh):

    if float(seq) <= 1.0:
        seq = float(seq) * 100.0
    if float(seq) < seq_thresh:
        return True
    if float(res) > res_thresh:
        return True
    else:
        return False


# template-structure: [pdb_id,seq_id,chain,aln_length,resolution,ligands,r-value,templateSelectionScore]


def good(structure, seq_threshold, resolution_threshold, cov_threshold):
    seq_id = structure['Seq_Id']
    resolution = structure['Resolution']
    cov = structure['Coverage']
    if weakCriteria(seq_id, resolution, cov, seq_threshold, resolution_threshold, cov_threshold):
        return False
    return True

