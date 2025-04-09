import gzip
import os

import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request

from structman.lib.sdsc.consts import residues as residue_consts
from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc import residue as residue_package

chain_order = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz'

def getCurrentPDBID(pdb_id, pdb_path, debug=False):
    obsolete_path = '%s/data/structures/obsolete/pdb/%s/pdb%s.ent.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())

    if os.path.isfile(obsolete_path):
        buf = gzip.open(obsolete_path, 'rb')
        # The replacement PDB ID is always given in the second line
        next(buf)
        try:
            pdb_id = buf.readline().decode('ascii').split()[3]
        except:
            print('Error in obsolete_check:', pdb_id)

        buf.close()
    elif debug:
        print('Is not a file:', obsolete_path)

    return pdb_id

def test_for_AU(pdb_id, pdb_path):
    path = '%s/data/biounit/PDB/divided/%s/%s.pdb1.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
    au_path = '%s/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
    AU = os.path.isfile(au_path) and not os.path.isfile(path)
    return AU

def getStructureBuffer(pdb_id, pdb_path, AU=False, obsolete_check=False, get_is_local=False, verbosity=0, model_path = None):
    if model_path is None:
        buf = getPDBBuffer(pdb_id, pdb_path, AU=AU, obsolete_check=obsolete_check, verbosity=verbosity)
    elif model_path[-3:] == '.gz':
        buf = gzip.open(model_path, 'rb')
    else:
        buf = open(model_path, 'rb')
    if get_is_local:
        return buf, model_path
    return buf


def getPDBBuffer(pdb_id: str, pdb_path: str, AU: bool =False, obsolete_check: bool =False, get_is_local: bool =False, verbosity: int=0) -> None | gzip.GzipFile:
    if obsolete_check:
        pdb_id = getCurrentPDBID(pdb_id, pdb_path)

    if AU:
        path = '%s/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
        if not os.path.isfile(path):
            pdb_id = getCurrentPDBID(pdb_id, pdb_path)
            path = '%s/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
            if not os.path.isfile(path):
                url = 'https://files.rcsb.org/view/%s.pdb' % pdb_id
                sdsc_utils.connection_sleep_cycle(verbosity, url)
                if pdb_path != '' and verbosity >= 2:
                    print("Did not find asymetric unit entry in local pdb, searching online: ", url)
                try:
                    request = urllib.request.Request(url)
                    if get_is_local:
                        return urllib.request.urlopen(request, timeout=60), None
                    return urllib.request.urlopen(request, timeout=60)
                except:
                    if verbosity >= 2:
                        print("Did not find the PDB-file (AU is true): %s" % pdb_id)
                    if get_is_local:
                        return None, None
                    return None
            else:
                if get_is_local:
                    return gzip.open(path, 'rb'), path
                return gzip.open(path, 'rb')

        else:
            if get_is_local:
                return gzip.open(path, 'rb'), path
            return gzip.open(path, 'rb')
    else:
        path = '%s/data/biounit/PDB/divided/%s/%s.pdb1.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
        if not os.path.isfile(path):
            pdb_id = getCurrentPDBID(pdb_id, pdb_path)
            path = '%s/data/biounit/PDB/divided/%s/%s.pdb1.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
            if not os.path.isfile(path):
                if test_for_AU(pdb_id, pdb_path):
                    return getPDBBuffer(pdb_id, pdb_path, AU=True, obsolete_check=obsolete_check, get_is_local=get_is_local, verbosity=verbosity)

                url = 'https://files.rcsb.org/view/%s.pdb1' % pdb_id
                sdsc_utils.connection_sleep_cycle(verbosity, url)
                if pdb_path != '' and verbosity >= 2:
                    print("Did not find entry in local pdb, searching online: ", url)
                try:
                    request = urllib.request.Request(url)
                    if get_is_local:
                        return urllib.request.urlopen(request, timeout=60), None
                    return urllib.request.urlopen(request, timeout=60)
                except:
                    if verbosity >= 3:
                        print("Did not find the PDB-file: %s, trying AU" % pdb_id)
                    return getPDBBuffer(pdb_id, pdb_path, AU=True, obsolete_check=obsolete_check, get_is_local=get_is_local, verbosity=verbosity)
            else:
                if get_is_local:
                    return gzip.open(path, 'rb'), path
                return gzip.open(path, 'rb')

        else:
            if get_is_local:
                return gzip.open(path, 'rb'), path
            return gzip.open(path, 'rb')

# Newly Written
def getMMCIFBuffer(pdb_id, pdb_path, AU=False, obsolete_check=False):
    if obsolete_check:
        pdb_id = getCurrentPDBID(pdb_id, pdb_path)

    if AU:
        path = '%s/data/structures/divided/pdb/%s/pdb%s.cif.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
        if not os.path.isfile(path):
            url = 'https://files.rcsb.org/view/%s.cif' % pdb_id

            if path != '':
                print("Did not find asymetric unit entry in local mmcif, searching online: ", url)
            try:
                request = urllib.request.Request(url)
                return urllib.request.urlopen(request)
            except:
                print("Did not find the MMCIF-file: %s" % pdb_id)
                return None
        else:
            return gzip.open(path, 'rb')
    else:
        path = '%s/data/biounit/PDB/divided/%s/%s.cif.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
        if not os.path.isfile(path):
            url = 'https://files.rcsb.org/view/%s.cif' % pdb_id

            if pdb_path != '':
                print("Did not find entry in local mmcif, searching online: ", url)
            try:
                request = urllib.request.Request(url)
                print('this is request', request)
                print('this is url open', urllib.request.urlopen(request))
                return urllib.request.urlopen(request)
            except:
                print("Did not find the MMCIF-file: %s" % pdb_id)
                return None
        else:
            return gzip.open(path, 'rb')


# nucleic acid parser called by mmcif parser with only buf if chain is nuc acid
def nuc_acid_parser(buf):
    i = 1
    res_pos_map = {}
    seq = ""
    used_res = ()
    temp = 0

    x = []
    for line in open(buf):  # What is this?
        splitted = line.split()
        if len(line.strip()) == 0:
            continue
        if splitted[0] == 'ATOM' or splitted[0] == 'HETATM':
            chain_id = splitted[6]
            if chain_id != chain:  # BUG: undefined variable
                continue
            if temp == splitted[8]:
                continue
            res_name = splitted[5]
            nuc = res_name
            seq = seq + nuc
            temp = splitted[8]
            res_pos_map[splitted[8]] = i
            i += 1
    buf.close()
    return seq, res_pos_map


# mmcif parser with only takes buf as an argument
def parseMMCIFSequence(buf):
    i = 1
    res_pos_map = {}
    seq = ""
    used_res = ()
    temp = 0

    x = []
    for line in open(buf):
        splitted = line.split()
        if len(line.strip()) == 0:
            continue
        if splitted[0] == 'ATOM' or splitted[0] == 'HETATM':
            chain_id = splitted[6]
            nuc_acid = splitted[17]
            if nuc_acid == 'T' or nuc_acid == 'C' or nuc_acid == 'G' or nuc_acid == 'U' or nuc_acid == 'A':
                seq, res_pos_map = nuc_acid_parser(buf)
                break
            else:
                res_name = splitted[5]
                if res_name not in residue_consts.THREE_TO_ONE:
                    continue
                if temp == splitted[8]:
                    continue
                aa = residue_consts.THREE_TO_ONE[res_name]
                seq = seq + aa
                temp = splitted[8]
                res_pos_map[splitted[8]] = i
                i += 1
    buf.close()
    return seq, res_pos_map


def getSequence(pdb_id: str, chains: set[str], pdb_path: str) -> dict[str, tuple[str, dict[str, int]]]:
    AU: bool = test_for_AU(pdb_id, pdb_path)
    obsolete_check: bool = not AU
    buf: gzip.GzipFile = getPDBBuffer(pdb_id, pdb_path, AU=AU, obsolete_check=obsolete_check)
    
    if buf is None:
        return None

    chain_seq_dict: dict[str, tuple[str, dict[str, int]]] = parsePDBSequence(buf, chains)
    buf.close()

    return chain_seq_dict

# called by serializedPipeline

def getSequences(pdb_id_tuples: list[str], pdb_path: str) -> list[tuple[str, str, dict]]:
    pdb_sequences: list[tuple[str, str, dict]] = []

    pdb_chains_tuples: dict[str, set[str]] = {}
    for pdb_chain_tuple in pdb_id_tuples:
        [pdb, chain] = pdb_chain_tuple.split(':')
        if pdb not in pdb_chains_tuples:
            pdb_chains_tuples[pdb] = set([chain])
        else:
            pdb_chains_tuples[pdb].add(chain)
            
    for pdb in pdb_chains_tuples:
        chains: set[str] = pdb_chains_tuples[pdb]
        chain_seq_dict: dict[str, tuple[str, dict[str, int]]] = getSequence(pdb, chains, pdb_path)
        for chain in chain_seq_dict:
            seq, res_pos_map = chain_seq_dict[chain]
            pdb_chain_tuple = f'{pdb}:{chain}'
            pdb_sequences.append((pdb_chain_tuple, seq, res_pos_map))

    return pdb_sequences

def get_chains_from_pdb_id(pdb_id, pdb_path):
    AU = test_for_AU(pdb_id, pdb_path)
    buf = getPDBBuffer(pdb_id, pdb_path, AU=AU, obsolete_check=True)
    return get_all_chains(buf)


def parsePDBSequence(buf: gzip.GzipFile, chains: set[str]) -> dict[str, tuple[str, dict[str, int]]]:
    chain_seq_dict: dict[str, tuple[str | list[str], dict[str, int]]] = {}
    used_res: dict[str, set[str]] = {}
    
    firstAltLoc = None
    try:
        for line in buf:
            line = line.decode('ascii')

            if len(line) > 26:
                record_name = line[0:6].strip()
                if not record_name == 'ATOM' \
                        and not record_name == 'HETATM':
                    continue

                altLoc = line[16]
                if firstAltLoc is None and altLoc != ' ':
                    firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                    continue

                res_name = line[17:20].strip()
                if record_name == 'HETATM':
                    if res_name not in residue_consts.THREE_TO_ONE:
                        continue

                chain_id = line[21]
                if not chain_id in chains:
                    continue
                if chain_id not in chain_seq_dict:
                    chain_seq_dict[chain_id] = ([], {})
                    used_res[chain_id] = set()

                #atom_nr = line[6:11].strip()
                res_nr = line[22:27].strip()  # this includes the insertion code

                if res_nr not in used_res[chain_id]:
                    if res_name not in residue_consts.THREE_TO_ONE:
                        continue
                    aa = residue_consts.THREE_TO_ONE[res_name][0]
                    chain_seq_dict[chain_id][0].append(aa)
                    used_res[chain_id].add(res_nr)
                    chain_seq_dict[chain_id][1][res_nr] = len(chain_seq_dict[chain_id][1]) + 1
    except OSError:
        pass

    for chain_id in chain_seq_dict:
        chain_seq_dict[chain_id] = (''.join(chain_seq_dict[chain_id][0]), chain_seq_dict[chain_id][1])

    return chain_seq_dict



def parse_chaintype_map(input_page):
    """
    Parses a PDB-file and takes all atomic coordinates.

    Input:
    input_page: String ; content of a pdb file

    Output:
    chaintype_map
    """

    lines = input_page.split('\n')

    chain_type_map = {}
    chain_type = '-'

    peptide_count = {}

    rare_residues = set()

    firstAltLoc = None

    for line in lines:
        if len(line) > 5:
            record_name = line[0:6].replace(" ", "")
            if record_name == "ENDMDL":
                break
            elif record_name == 'SEQRES':
                for tlc in line[19:].split():
                    tlc = tlc.strip()
                    if len(tlc) != 3:
                        continue
                    if tlc not in residue_consts.THREE_TO_ONE:
                        rare_residues.add(tlc)

        # ignore short lines
        if len(line) > 20:
            atom_nr = line[6:11].replace(" ", "")
            if record_name.count('ATOM') > 0 and record_name != 'ATOM':  # 100k atom bug fix
                atom_nr = '%s%s' % (record_name[4:], atom_nr)
                record_name = 'ATOM'

            res_name = line[17:20].replace(" ", "")

            if len(line) > 21:

                chain_id = line[21]

                if record_name == "ATOM" or record_name == "HETATM":
                    altLoc = line[16]
                    if firstAltLoc is None and altLoc != ' ':
                        firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                    if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue


            if record_name == "ATOM" or record_name == 'MODRES' or record_name == "HETATM":
                if chain_id not in chain_type_map:
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
                        if res_name in residue_consts.THREE_TO_ONE or not sdsc_utils.boring(res_name) or (res_name in rare_residues):
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
                    if res_name in residue_consts.THREE_TO_ONE or not sdsc_utils.boring(res_name) or (res_name in rare_residues):  # For a hetero peptide 'boring' hetero amino acids are allowed as well as other non boring molecules not in THREE_TO_ONE, which are hopefully some kind of anormal amino acids
                        peptide_count[chain_id][1] += 1
                    elif len(res_name) < 3:
                        peptide_count[chain_id][0] += 1

                elif record_name == "ATOM" or record_name == 'MODRES':
                    peptide_count[chain_id][0] += 1
                    if len(res_name) == 1:
                        chain_type_map[chain_id] = "RNA"
                    elif len(res_name) == 2:
                        chain_type_map[chain_id] = "DNA"


    # New Peptide detection counts irregular amino acids
    for chain_id in peptide_count:
        if chain_type_map[chain_id] == 'DNA' or chain_type_map[chain_id] == 'RNA':
            continue
        if peptide_count[chain_id][0] <= peptide_count[chain_id][1] * 2:
            chain_type_map[chain_id] = 'Peptide'
        elif (peptide_count[chain_id][0] + peptide_count[chain_id][1]) < 150:  # Total number of atoms
            chain_type_map[chain_id] = 'Peptide'

    return chain_type_map


def get_all_chains(buf):
    chains = set()
    for line in buf:
        line = line.decode('ascii')

        if len(line) > 26:
            chain_id = line[21]
            chains.add(chain_id)
    return chains

# NewMMCIFParser
def nuc_acid_parser(buf, chain):
    i = 1
    res_pos_map = {}
    seq = ""
    used_res = ()
    temp = 0

    x = []
    for line in buf:
        line = line.decode('ascii')
        splitted = line.split()
        if len(line.strip()) == 0:
            continue
        if splitted[0] == 'ATOM' or splitted[0] == 'HETATM':
            chain_id = splitted[6]
            if chain_id != chain:
                continue
            if temp == splitted[8]:
                continue
            res_name = splitted[5]
            nuc = res_name
            seq = seq + nuc
            temp = splitted[8]
            res_pos_map[splitted[8]] = i
            i += 1

    return seq, res_pos_map


def parseMMCIFSequence(buf, chain):
    i = 1
    res_pos_map = {}
    seq = ""
    used_res = ()
    temp = 0

    x = []
    for line in buf:
        line = line.decode('ascii')
        splitted = line.split()
        if len(line.strip()) == 0:
            continue
        if splitted[0] == 'ATOM' or splitted[0] == 'HETATM':
            chain_id = splitted[6]
            if chain_id != chain:
                continue
            nuc_acid = splitted[17]
            if nuc_acid == 'T' or nuc_acid == 'C' or nuc_acid == 'G' or nuc_acid == 'U' or nuc_acid == 'A':
                seq, res_pos_map = nuc_acid_parser(buf, chain)
                break
            else:
                res_name = splitted[5]
                if res_name not in residue_consts.THREE_TO_ONE:
                    continue
                if temp == splitted[8]:
                    continue
                aa = residue_consts.THREE_TO_ONE[res_name]
                seq = seq + aa
                temp = splitted[8]
                res_pos_map[splitted[8]] = i
                i += 1

    return seq, res_pos_map

# called by serializedPipeline
# called by templateFiltering
# called by sdsc


def standardParsePDB(pdb_id, pdb_path, obsolete_check=False, return_10k_bool=False, get_is_local=False, verbosity=0, model_path=None, only_first_model = False):
    AU = False
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True
        only_first_model = True

    buf = getStructureBuffer(pdb_id, pdb_path, AU=AU, obsolete_check=obsolete_check, verbosity=verbosity, get_is_local=get_is_local, model_path = model_path)

    if get_is_local:
        buf, path = buf

    if buf is None:
        if return_10k_bool and get_is_local:
            return '', None, None, 0
        elif return_10k_bool:
            return '', None, 0
        elif get_is_local:
            return '', None, 0
        return '', 0
    # Added automatization here, if pdb file not exist, automatically calls MMCIF version
    #    print('PDB file does not exist, trying MMCIF file')
    #    return standardParseMMCIF(pdb_id,pdb_path,obsolete_check=False)

    chain_ids = set()
    rare_residues = set()
    newlines = []

    page_altered = False
    not_crystal = False
    multi_model_mode = False
    multi_model_chain_dict = {}
    chain_order_pos_marker = 0
    asymetric_chain_number = None
    current_model_id = 0
    multi_model_chain_dict = {current_model_id: {}}

    current_serial = 1
    firstAltLoc = None
    for line in buf:
        if len(line) >= 27:
            record_name = line[0:6].rstrip()
            line = line.rstrip(b'\n')
        else:
            continue

        if record_name == b'EXPDTA':
            words = line.split()
            if len(words) < 3:
                continue
            if words[2] == b'NMR':
                not_crystal = True

            if words[1] == b'SOLUTION' and (words[2] == b'SCATTERING' or words[2] == b'SCATTERING;' or words[2] == b'NMR;' or words[2] == b'NMR'):
                not_crystal = True

            if line.count(b'NMR') > 0:
                not_crystal = True
            newlines.append(line)
            continue

        elif record_name == b'MODEL':
            current_model_id = line[10:14]
            if not current_model_id in multi_model_chain_dict:
                multi_model_chain_dict[current_model_id] = {}
            if not multi_model_mode:
                newlines.append(line)
            elif len(multi_model_chain_dict) * asymetric_chain_number >= len(chain_order):
                break
            continue

        elif record_name == b'ENDMDL':
            if not_crystal or only_first_model:
                newlines.append(line)
                break
            elif not multi_model_mode:
                multi_model_mode = True
                page_altered = True
                asymetric_chain_number = len(chain_ids)
            continue

        elif record_name == b'SEQRES':
            for tlc in line[19:].split():
                tlc = tlc.strip().decode('ascii')
                if len(tlc) != 3:
                    continue
                if tlc not in residue_consts.THREE_TO_ONE:
                    rare_residues.add(tlc)
        elif record_name == b'SSBOND':
            newlines.append(line)
            continue
        elif record_name == b'LINK  ':
            newlines.append(line)
            continue
        elif record_name == b'CISPEP':
            newlines.append(line)
            continue
        elif record_name == b'SEQRES':
            newlines.append(line)
            continue
        elif record_name == b'TER':
            # This can happen, if the whole chain consists of boring ligands,
            # which is very boring
            if chain_id not in chain_ids:
                continue
        elif record_name == b'MODRES':
            pass
        elif record_name == b'ATOM' or record_name == b'HETATM':
            altLoc = chr(line[16])
            if firstAltLoc is None and altLoc != ' ':
                firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
            if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                continue
            res_name = line[17:20].strip().decode('ascii')
            chain_id = line[21:22].decode('ascii')

            if record_name == b'HETATM':
                if res_name not in residue_consts.THREE_TO_ONE and sdsc_utils.boring(res_name) and (len(res_name) == 3) and (res_name not in rare_residues):
                    continue

            if chain_id not in chain_ids:
                chain_ids.add(chain_id)
                multi_model_chain_dict[current_model_id][chain_id] = chain_id
            elif multi_model_mode:
                if not chain_id in multi_model_chain_dict[current_model_id]:
                    new_chain_id = None
                    while new_chain_id is None:
                        if not chain_order[chain_order_pos_marker] in chain_ids:
                            new_chain_id = chain_order[chain_order_pos_marker]
                        chain_order_pos_marker += 1
                    multi_model_chain_dict[current_model_id][chain_id] = new_chain_id
                    chain_id = new_chain_id
                else:
                    chain_id = multi_model_chain_dict[current_model_id][chain_id]

        else:
            # We are not interested in other record types
            continue

        # Adjust serials after skipping records
        if current_serial > 99999:
            page_altered = True
            newline = record_name[0:5].ljust(5, b' ') + str(current_serial).encode('ascii').rjust(6, b' ')
        else:
            newline = record_name.ljust(6, b' ') + str(current_serial).encode('ascii').rjust(5, b' ')

        newline += line[11:]

        if multi_model_mode:
            newline = newline[:21] + chain_id.encode('ascii') + newline[22:]

        newlines.append(newline)
        current_serial += 1

    newlines.append(b'ENDMDL                                                                          \n')

    template_page = b'\n'.join(newlines).decode('ascii')

    buf.close()

    if return_10k_bool and get_is_local:
        return template_page, page_altered, path, current_serial
    elif return_10k_bool:
        return template_page, page_altered, current_serial
    elif get_is_local:
        return template_page, path, current_serial
    return template_page, current_serial

# called by modelling


def relocate_hetatm(page, filter_chains=None, filter_het=None):
    chain_lines = {}
    chain_order = []
    removed_ligands = {}
    no_double_remove = {}
    rare_residues = set()
    for line in page.split('\n'):
        if len(line) >= 27:
            record_name = line[0:6].strip()
        else:
            continue

        if record_name == 'ENDMDL':
            break
        elif record_name == 'SEQRES':
            for tlc in line[19:].split():
                tlc = tlc.strip()
                if len(tlc) != 3:
                    continue
                if tlc not in residue_consts.THREE_TO_ONE:
                    rare_residues.add(tlc)

        elif record_name == 'ATOM' or record_name == 'HETATM':
            chain_id = line[21:22]

            if filter_chains is not None:
                if chain_id not in filter_chains:
                    continue

            if filter_het is not None:
                res_name = line[17:20].strip()
                res_nr = line[22:27].strip()
                if res_name in filter_het:
                    if chain_id not in removed_ligands:
                        removed_ligands[chain_id] = 0
                        no_double_remove[chain_id] = set()
                    if res_nr in no_double_remove[chain_id]:
                        continue
                    removed_ligands[chain_id] += 1
                    no_double_remove[chain_id].add(res_nr)
                    continue
            #atom = line[77]
            # if atom == 'H':
            #    continue
            res_name = line[17:20].strip()

            if record_name == 'HETATM':
                if res_name not in residue_consts.THREE_TO_ONE and sdsc_utils.boring(res_name) and (res_name not in rare_residues):
                    continue

            if chain_id not in chain_lines:
                chain_lines[chain_id] = []
                chain_order.append(chain_id)
            chain_lines[chain_id].append(line)

    new_page = ''
    for chain in chain_order:
        new_page += '\n'.join(chain_lines[chain]) + '\nTER\n'
    new_page += 'ENDMDL\n'
    return new_page, removed_ligands


# Called by standardParsePDB if PDB file does not exist
def standardParseMMCIF(pdb_id, pdb_path, obsolete_check=False):
    AU = False
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True

    buf = getMMCIFBuffer(pdb_id, pdb_path, AU=AU, obsolete_check=obsolete_check)

    if buf is None:
        return None

    chain_ids = set()
    newlines = []

    current_serial = 1
    firstAltLoc = None
    for line in buf:
        if len(line) >= 27:
            record_name = line[0:6].rstrip()
            line = line.rstrip(b'\n')
            splitted = line.split()
        else:
            continue

        if record_name == b'ENDMDL':
            newlines.append(line)
            break
        elif record_name == b'TER':
            # This can happen, if the whole chain consists of boring ligands,
            # which is very boring
            if chain_id not in chain_ids:
                continue
        elif record_name == b'MODRES':
            pass
        elif record_name == b'ATOM' or record_name == b'HETATM':
            altLoc = chr(line[17])
            if firstAltLoc is None and altLoc != ' ':
                firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
            if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                continue
            res_name = splitted[5]
            chain_id = splitted[6]
            res_nr = splitted[8]

            if record_name == b'HETATM':
                if res_name not in residue_consts.THREE_TO_ONE and sdsc_utils.boring(res_name):
                    continue

            if chain_id not in chain_ids:
                chain_ids.add(chain_id)

        else:
            # We are not interested in other record types
            continue

        # Adjust serials after skipping records
        if current_serial > 99999:
            newline = record_name[0:5].ljust(5, b' ') + str(current_serial).encode('ascii').rjust(6, b' ')
        else:
            newline = record_name.ljust(6, b' ') + str(current_serial).encode('ascii').rjust(5, b' ')
        newline += line[11:]
        newlines.append(newline)
        current_serial += 1

    template_page = b'\n'.join(newlines).decode('ascii')

    buf.close()

    return template_page


# called by serializedPipeline
def getStandardizedPdbFile(
        pdb_id: str,
        pdb_path: str,
        oligo: set[str] = set(),
        verbosity: int = 0,
        model_path: str | None = None,
        obsolete_check: bool = False) -> tuple[list[bytes], list, dict[str, str], set[str], int, list[str], set[str]]:

    AU: bool = False

    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True

    buf = getStructureBuffer(pdb_id, pdb_path, AU=AU, obsolete_check=obsolete_check, verbosity=verbosity, model_path = model_path)

    if buf is None:
        # Added automatization here, if pdb file not exist, automatically calls MMCIF version
        #print('PDB file does not found, trying MMCIF file')
        # return getStandardizedMMCIFFile(pdb_id=pdb_id,pdb_path=pdb_path)
        return 'Buffer was None'

    chain_type_map: dict[str, str] = {}
    chain_list: list[str] = []
    modres_map: dict[str, residue_package.Residue_Map] = {}

    lig_set: dict[str, set[int | bytes]] = {}

    newlines: list[bytes] = []

    interaction_partners: list[tuple] = []

    rare_residues: set[bytes] = set()
    chain_ids: set[str] = set()

    not_crystal: bool = False
    multi_model_mode: bool = False
    chain_order_pos_marker: int = 0
    asymetric_chain_number: int | None = None
    current_model_id: int = 0
    multi_model_chain_dict: dict[int, dict[str, str]] = {current_model_id:{}}

    current_serial: int = 1
    firstAltLoc: bytes | None = None

    for line in buf:
        
        if len(line) >= 27:
            try:
                record_name: bytes = line[0:6].rstrip()
                line: bytes = line.rstrip(b'\n')
            except:
                return f'Error in parsing PDB file {pdb_id}, line:\n{line}'
        else:
            continue

        if record_name == b'EXPDTA':
            words: list[bytes] = line.split()
            if len(words) < 3:
                continue
            if words[2] == b'NMR':
                not_crystal = True

            if words[1] == b'SOLUTION' and (words[2] == b'SCATTERING' or words[2] == b'SCATTERING;' or words[2] == b'NMR;' or words[2] == b'NMR'):
                not_crystal = True

            if line.count(b'NMR') > 0:
                not_crystal = True
            newlines.append(line)
            continue

        elif record_name == b'MODEL':
            current_model_id = int(line[10:14])
            if not current_model_id in multi_model_chain_dict:
                multi_model_chain_dict[current_model_id] = {}
            if not multi_model_mode:
                newlines.append(line)
            elif len(multi_model_chain_dict) * asymetric_chain_number >= len(chain_order):
                break
            continue

        elif record_name == b'ENDMDL':
            if not_crystal or AU:
                newlines.append(line)
                break
            elif not multi_model_mode:
                multi_model_mode = True
                asymetric_chain_number = len(chain_ids)
            continue

        elif record_name == b'SEQRES':
            for tlc in line[19:].split():
                tlc = tlc.strip()
                if len(tlc) != 3:
                    continue
                if tlc not in residue_consts.BIN_THREE_TO_ONE:
                    rare_residues.add(tlc)

        elif record_name == b'TER':
            # This can happen, if the whole chain consists of boring ligands,
            # which is very boring
            if chain_id not in chain_type_map:
                continue

        elif record_name == b'MODRES':
            chain_id: str = line[16:17].decode('ascii')
            res_nr: str = line[18:23].strip().decode('ascii')
            res_name: str = line[24:27].strip().decode('ascii')
            if chain_id not in modres_map:
                modres_map[chain_id] = residue_package.Residue_Map()
            modres_map[chain_id].add_item(res_nr, res_name)

        elif record_name == b'ATOM' or record_name == b'HETATM':
            altLoc: bytes = line[16:17]
            if firstAltLoc is None and altLoc != b' ':
                firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
            if altLoc != b' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                continue
            res_name: str = line[17:20].strip().decode('ascii')
            chain_id: str = line[21:22].decode('ascii')

            if chain_id not in chain_ids:
                chain_ids.add(chain_id)
                multi_model_chain_dict[current_model_id][chain_id] = chain_id
            elif multi_model_mode:
                if not chain_id in multi_model_chain_dict[current_model_id]:
                    new_chain_id = None
                    while new_chain_id is None:
                        if not chain_order[chain_order_pos_marker] in chain_ids:
                            new_chain_id = chain_order[chain_order_pos_marker]
                        chain_order_pos_marker += 1
                    multi_model_chain_dict[current_model_id][chain_id] = new_chain_id
                    chain_id = new_chain_id
                else:
                    chain_id = multi_model_chain_dict[current_model_id][chain_id]

            res_nr: str = line[22:27].strip().decode('ascii')

            #occupancy = float(line[54:60].replace(" ",""))

            if chain_id not in chain_type_map:
                if record_name == b'ATOM':
                    if len(res_name) == 1:
                        chain_type = 'RNA'
                    elif len(res_name) == 2:
                        chain_type = 'DNA'
                    elif len(res_name) == 3:
                        chain_type = 'Protein'
                    chain_type_map[chain_id] = chain_type
                    chain_list.append(chain_id)
                elif record_name == b'HETATM':
                    # For a hetero peptide 'boring' hetero amino acids are
                    # allowed [as well as other non boring molecules not in
                    # threeToOne, which are hopefully some kind of abnormal amino
                    # acids] did not work
                    if res_name in residue_consts.BIN_THREE_TO_ONE or res_name in rare_residues:  # or not sdsc_utils.boring(res_name):
                        chain_type = 'Peptide'
                        chain_type_map[chain_id] = chain_type
                        chain_list.append(chain_id)

            elif record_name == b'ATOM' and chain_type_map[chain_id] == 'Peptide':
                if len(res_name) == 1:
                    chain_type = "RNA"
                elif len(res_name) == 2:
                    chain_type = "DNA"
                elif len(res_name) == 3:
                    chain_type = "Protein"
                chain_type_map[chain_id] = chain_type

            if record_name == b'HETATM':
                # modified residue are not parsed as ligands
                if (res_name not in residue_consts.BIN_THREE_TO_ONE) and (res_name not in rare_residues):
                    if not sdsc_utils.bin_boring(res_name):
                        if chain_id not in lig_set:
                            lig_set[chain_id] = set()
                        if res_nr not in lig_set[chain_id]:
                            lig_set[chain_id].add(res_nr)

                            if chain_id in modres_map:
                                if modres_map[chain_id].contains(res_nr):
                                    continue

                            interaction_partners.append(["Ligand", res_name, res_nr, chain_id])
                    elif len(res_name) == 3:
                        continue

        else:
            # We are not interested in other record types
            continue

        # Adjust serials after skipping records
        if current_serial > 99999:
            newline = record_name[0:5].ljust(5, b' ') + str(current_serial).encode('ascii').rjust(6, b' ')
        else:
            newline = record_name.ljust(6, b' ') + str(current_serial).encode('ascii').rjust(5, b' ')
        newline += line[11:]

        if multi_model_mode:
            newline = newline[:21] + chain_id.encode('ascii') + newline[22:]

        newlines.append(newline)
        current_serial += 1

    buf.close()

    newlines.append(b'ENDMDL                                                                          \n')

    pep_not_lig = []
    for chain_id in chain_type_map:
        chain_type = chain_type_map[chain_id]
        interaction_partners.append([chain_type, chain_id])

    for chain_id in chain_type_map:
        chain_type = chain_type_map[chain_id]
        if chain_type == 'Peptide':
            for i, iap in enumerate(interaction_partners):
                if iap[0] == 'Ligand':
                    if iap[3] == chain_id:
                        pep_not_lig.append(i)

    pep_not_lig = sorted(pep_not_lig, reverse=True)

    for i in pep_not_lig:
        del interaction_partners[i]

    # Happens for repeated chains in asymetric units, which do not occur in the
    # biological assembly
    irregular_homo_chains = []
    for homo_chain in oligo:
        if homo_chain not in chain_type_map:
            irregular_homo_chains.append(homo_chain)

    for chain in irregular_homo_chains:
        oligo.remove(chain)

    return newlines, interaction_partners, chain_type_map, oligo, current_serial, chain_list, rare_residues

# Newly written MMCIF version
# Called by getStandardizedPDBFile if PDB file does not exist


def getStandardizedMMCIFFile(pdb_id, chain, structure, pdb_path):
    times = [0.0, 0.0, 0.0]
    t0 = time.time()
    AU = False

    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True

    buf = getMMCIFBuffer(pdb_id, pdb_path, AU)

    if buf is None:
        return "", '', "Did not find the MMCIF-file: %s" % (pdb_id), times

    structure['IAP'] = []
    t1 = time.time()

    chain_type_map = {}
    modres_map = {}

    lig_set = {}

    newlines = []

    current_serial = 1
    firstAltLoc = None
    for line in buf:
        line = line.decode('ascii')
        if len(line) >= 27:
            record_name = line[0:6].rstrip()
            line = line.rstrip('\n')
            splitted = line.split()
        else:
            continue

        if record_name == 'ENDMDL':
            newlines.append(line)
            break
        elif record_name == 'TER':
            # This can happen, if the whole chain consists of boring ligands,
            # which is very boring
            if chain_id not in chain_type_map:
                continue
        # No idea id there is modres in mmcif, should check
        elif record_name == 'MODRES':
            chain_id = line[16]
            res_nr = line[18:23].strip()
            res_name = line[24:27].strip()
            if (chain_id, res_nr) not in modres_map:
                modres_map[(chain_id, res_nr)] = res_name
        elif record_name == 'ATOM' or record_name == 'HETATM':
            altLoc = line[17]  # Check
            if firstAltLoc is None and altLoc != ' ':
                firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
            if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                continue
            res_name = splitted[5]
            chain_id = splitted[6]
            res_nr = splitted[8]

            #occupancy = float(line[54:60].replace(" ",""))

            if chain_id not in chain_type_map:
                if record_name == 'ATOM':
                    if len(res_name) == 1:
                        chain_type = 'RNA'
                    elif len(res_name) == 2:
                        chain_type = 'DNA'
                    elif len(res_name) == 3:
                        chain_type = 'Protein'
                    chain_type_map[chain_id] = chain_type
                elif record_name == 'HETATM':
                    # For a hetero peptide 'boring' hetero amino acids are
                    # allowed as well as other non boring molecules not in
                    # residue_consts.THREE_TO_ONE, which are hopefully some kind of anormal amino
                    # acids
                    if res_name in residue_consts.THREE_TO_ONE or not sdsc_utils.boring(res_name):
                        chain_type = 'Peptide'
                        chain_type_map[chain_id] = chain_type
            elif record_name == 'ATOM' and chain_type_map[chain_id] == 'Peptide':
                if len(res_name) == 1:
                    chain_type = "RNA"
                elif len(res_name) == 2:
                    chain_type = "DNA"
                elif len(res_name) == 3:
                    chain_type = "Protein"
                chain_type_map[chain_id] = chain_type

            if record_name == 'HETATM':
                # modified residue are not parsed as ligands
                if res_name not in residue_consts.THREE_TO_ONE and not sdsc_utils.boring(res_name):
                    if chain_id not in lig_set:
                        lig_set[chain_id] = set()
                    if res_nr not in lig_set[chain_id]:
                        lig_set[chain_id].add(res_nr)

                        if (chain_id, res_nr) in modres_map:
                            continue

                        structure['IAP'].append(["Ligand", res_name, res_nr, chain_id])

        else:
            # We are not interested in other record types
            continue

        # Adjust serials after skipping records
        if current_serial > 99999:
            newline = record_name[0:5].ljust(5, ' ') + str(current_serial).rjust(6, ' ')
        else:
            newline = record_name.ljust(6, ' ') + str(current_serial).rjust(5, ' ')
        newline += line[11:]
        newlines.append(newline)
        current_serial += 1

    pep_not_lig = []
    for chain_id in chain_type_map:
        chain_type = chain_type_map[chain_id]
        structure['IAP'].append([chain_type, chain_id])

    for chain_id in chain_type_map:
        chain_type = chain_type_map[chain_id]
        if chain_type == 'Peptide':
            for i, iap in enumerate(structure['IAP']):
                if iap[0] == 'Ligand':
                    if iap[3] == chain_id:
                        pep_not_lig.append(i)

    pep_not_lig = sorted(pep_not_lig, reverse=True)

    for i in pep_not_lig:
        del structure['IAP'][i]

    # Happens for repeated chains in asymetric units, which do not occur in the
    # biological assembly
    irregular_homo_chains = []
    for i, homo_chain in enumerate(structure['Oligo']):
        if homo_chain not in chain_type_map:
            irregular_homo_chains.append(i)

    irregular_homo_chains = sorted(irregular_homo_chains, reverse=True)
    for i in irregular_homo_chains:
        del structure['Oligo'][i]

    t2 = time.time()

    template_page = '\n'.join(newlines)

    t3 = time.time()
    times = [t1 - t0, t2 - t1, t3 - t2]

    return template_page, None, times, interaction_partners, chain_type_map, oligo  # BUG: undefined variable


# called by database
def parseLigandDB(smiles_path, inchi_path):
    ligand_db = {}

    f = open(smiles_path, 'rt')
    for line in f:
        words = line.rstrip('\n').split('\t')
        if len(words) < 2:
            continue
        smiles = words[0]
        name = words[1]
        ligand_db[name] = smiles
    f.close()

    f = open(inchi_path, 'rt')
    for line in f:
        words = line.rstrip('\n').split('\t')
        if len(words) < 2:
            continue
        inchi = words[0]
        name = words[1]
        smiles = ligand_db[name]
        if not isinstance(smiles, str):
            continue
        ligand_db[name] = (smiles, inchi)
    f.close()

    return ligand_db


# called by database
def updateLigandDB(new_ligands, smiles_path, inchi_path):
    smiles_lines = []
    inchi_lines = []
    for (name, smiles, inchi) in new_ligands:
        smiles_lines.append('%s\t%s' % (smiles, name))
        inchi_lines.append('%s\t%s' % (inchi, name))

    f = open(smiles_path, 'a')
    f.write('\n'.join(smiles_lines) + '\n')
    f.close()

    f = open(inchi_path, 'a')
    f.write('\n'.join(inchi_lines) + '\n')
    f.close()


# get called by database
def getSI(pdb_id, name, res, chain, pdb_path, config):

    AU = False
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]
        AU = True

    buf = getPDBBuffer(pdb_id, pdb_path, AU=AU, obsolete_check=False, verbosity=config.verbosity)

    newlines = []

    for line in buf:
        if len(line) >= 27:
            record_name = line[0:6].rstrip()
        else:
            continue
        if record_name == b'ENDMDL':
            break

        if record_name == b'HETATM':
            res_name = line[17:20].strip().decode('ascii')
            chain_id = line[21:22].decode('ascii')
            res_nr = line[22:27].strip().decode('ascii')
            if chain == chain_id and name == res_name and res_nr == res:
                newlines.append(line)

    buf.close()

    page = (b''.join(newlines)).decode('ascii')

    FNULL = open(os.devnull, 'w')
    try:
        i_err = None
        s_err = None
        inchiproc = subprocess.Popen(["babel", "-i", "pdb", "-o", "inchi"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        inchi, i_err = inchiproc.communicate(page, timeout=600)
        smilesproc = subprocess.Popen(["babel", "-i", "pdb", "-o", "smi"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        smiles, s_err = smilesproc.communicate(page, timeout=600)
    except:
        config.errorlog.add_warning(f"Babel Package seems not to be installed: {pdb_id} {name} {res} {chain}\n{i_err}\n{s_err}")
        return ("", "")

    smiles = smiles.split()
    if len(smiles) > 0:
        smiles = smiles[0]
    else:
        config.errorlog.add_warning("Ligand not found in database and Babel failed: %s %s %s %s\n%s\n%s" % (name, res, chain, pdb_id, i_err, s_err))
        return ("", "")

    inchi = inchi.replace('\n', '')

    return (smiles, inchi)


def getPDBHeaderBuffer(pdb_id: bytes | str, pdb_path: str, tries: int =0, verbosity: int =0):
    if isinstance(pdb_id, bytes):
        pdb_id = pdb_id.decode('ascii')
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]

    path = '%s/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
    if not os.path.isfile(path):
        url = 'https://files.rcsb.org/header/%s.pdb' % pdb_id
        # if pdb_path != '':
        #    print("Did not find asymetric unit entry in local pdb, searching online: ",url)
        sdsc_utils.connection_sleep_cycle(verbosity, url)
        try:
            request = urllib.request.Request(url)
            return urllib.request.urlopen(request, timeout=(tries + 1) * 10)
        except:
            if tries < 2:
                #print("Unable to connect to PDB for the Header: %s\n%s\nThis might due to bad connection, let's try again ...'" % (pdb_id,url))
                return getPDBHeaderBuffer(pdb_id, pdb_path, tries=tries + 1)
            else:
                #print("Unable to connect to PDB for the Header: %s\n%s'" % (pdb_id,url))
                return None
    else:
        return gzip.open(path, 'rb')


# Same thing with mmcif
def getMMCIFHeaderBuffer(pdb_id, pdb_path, tries=0):
    if pdb_id.count('_AU') == 1:
        pdb_id = pdb_id[0:4]

    path = '%s/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
    if not os.path.isfile(path):
        url = 'https://files.rcsb.org/header/%s.cif' % pdb_id
        if pdb_path != '':
            print("Did not find asymetric unit entry in local pdb, searching online: ", url)
        try:
            request = urllib.request.Request(url)
            return urllib.request.urlopen(request, timeout=(tries + 1) * 10)
        except:
            if tries < 2:
                print("Unable to connect ot PDB for the Header: %s\n%s\nThis might due to bad connection, let's try again ...'" % (pdb_id, url))
                return getMMCIFHeaderBuffer(pdb_id, pdb_path, tries=tries + 1)
            else:
                print("Unable to connect ot PDB for the Header: %s\n%s'" % (pdb_id, url))
                return None
    else:
        return gzip.open(path, 'rb')


# called by templateSelection
def getInfo(pdb_id, pdb_path) -> tuple[float | None, dict[str, set[str]]]:
    buf = getPDBHeaderBuffer(pdb_id, pdb_path)

    # Added automatization here
    if buf is None:
        #buf = getMMCIFHeaderBuffer(pdb_id,pdb_path)
        # if buf is None:
        #    return None,{}
        return None, {}
    abort = False
    resolution = None

    not_crystal = False

    homomer_dict: dict[str, set[str]] = {}

    multiple_chain_line = False

    for line in buf:
        if abort:
            break

        line = line.rstrip(b'\n')
        words = line.split()

        if len(words) > 0:
            if words[0] == b'COMPND':
                if line[11:16] == b'CHAIN' or multiple_chain_line:
                    if line[11:17] == b'CHAIN;':
                        continue
                    if line.count(b';') == 0:
                        multiple_chain_line = True
                    else:
                        multiple_chain_line = False

                    chains = set(line[10:].replace(b'CHAIN:', b'').replace(b' ', b'').replace(b';', b'').decode('ascii').split(','))

                    for chain in chains:
                        homomer_dict[chain] = chains

            # Filter out NMR-structures, if you decide to use them in future, remember Recoord
            elif words[0] == b'EXPDTA':
                if words[2] == b'NMR':
                    resolution = 2.5
                    not_crystal = True

                if words[1] == b'POWDER':
                    resolution = 100.0
                    abort = True

                if words[1] == b'SOLUTION' and (words[2] == b'SCATTERING' or words[2] == b'SCATTERING;' or words[2] == b'NMR;' or words[2] == b'NMR'):
                    resolution = 2.5
                    not_crystal = True

                # if line.count(b'ELECTRON MICROSCOPY'):
                #    resolution = 2.5
                #    not_crystal = True

                if line.count(b'NMR') > 0:
                    resolution = 2.5
                    not_crystal = True

            elif words[0] == b'REMARK':
                if not_crystal:
                    continue
                if len(words) > 2:
                    if words[1] == b'2':
                        if words[2] == b'RESOLUTION.':
                            if words[3] == b'NULL' or words[3] == b'NOT':
                                resolution = None
                            else:
                                resolution = float(words[3])
                            abort = True

            # Filter out C-alpha only structures
            elif words[0] == b'MDLTYP':
                if len(words) > 2:
                    if b' '.join(words[1:3]) == b'CA ATOMS':
                        resolution = 100.0
                        abort = True

            elif words[0] == b'SEQRES':
                break

            # the authors of 1L9U were not able to make a proper MDLTYP entry
            if words[0] == b'REMARK':
                if b'COORDINATES' in words \
                        and b'CONTAIN' in words \
                        and b'ONLY' in words \
                        and b'CA' in words:
                    resolution = 100.0
                    abort = True
    buf.close()

    if resolution is None:
        resolution = 100.0
    try:
        resolution = float(resolution)
    except:
        raise NameError("Resolution buggy for: %s" % pdb_id)
    return resolution, homomer_dict


chain_id_list = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]
