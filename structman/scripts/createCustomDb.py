import gzip
import os
import sys

import psutil
import ray

# if running as script, add local structman package to path
if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__))))))

import structman

from structman.lib import serializedPipeline
from structman.lib.sdsc.consts import residues as residue_consts
from structman.base_utils.ray_utils import ray_init, ray_hack

"""

New python script to create custom db, functions has been obtained from 
createPdbBaDb.py

"""

def sanity_check_for_filename(filename):

    not_accepted_characters = [",",";","!","?","/","%"]

    for char in not_accepted_characters:
        if char in filename:
            return False
    
    return True


@ray.remote(num_cpus=1)
def parseByAtom(parse_dump, model_db = False):
    ray_hack()
    chain_order = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz'
    divide_folder = parse_dump

    files = os.listdir(divide_folder)
    great_seq_map = {}
    bio_entries = set()

    error_files = []

    for filename in files:

        print("parseByAtom started for file: %s" %filename)

        if os.path.isdir("%s/%s" % (divide_folder, filename)):
            continue
        else:
            if not sanity_check_for_filename:
                print("Error: filename %s" % filename)
                print("Filename cannot contain specific characters: , ; ! ? / %")
                sys.exit()
            if not model_db:
                if not len(filename) > 0:
                    continue
                pdb_id = filename[:4].upper()

                bio_entries.add(pdb_id)
                path = "%s/%s" % (divide_folder, filename)

                if filename.endswith(".gz"):
                    f = gzip.open(path, 'rb')
                else:
                    f = open(path,'rb')
                
            else:
                pdb_id = '-'.join(filename.split('-')[:2])
                path = f'{divide_folder}/{filename}'

                if filename.endswith(".gz"):
                    f = gzip.open(path, 'rb')
                else:
                    f = open(path,'rb')
            page = f.read()
            f.close()

            try:
                page = page.decode('ascii')
            except:
                print("Couldn't decode page:", path)
                error_files.append(path)
                continue

            lines = page.split('\n')

            seqs = {}
            chain_ids = set()

            not_crystal = False
            multi_model_mode = False
            multi_model_chain_dict = {}
            chain_order_pos_marker = 0
            asymetric_chain_number = None
            current_model_id = 0
            multi_model_chain_dict[current_model_id] = {}

            used_res = {}
            firstAltLoc = None
            for line in lines:
                record_name = line[0:6].replace(" ", "")
                atom_nr = line[6:11].replace(" ", "")
                atom_name = line[12:16].replace(" ", "")
                res_name = line[17:20].replace(" ", "")

                if record_name == 'EXPDTA':
                    words = line.split()
                    if len(words) < 3:
                        continue
                    if words[2] == 'NMR':
                        not_crystal = True

                    if words[1] == 'SOLUTION' and (words[2] == 'SCATTERING' or words[2] == 'SCATTERING;' or words[2] == 'NMR;' or words[2] == 'NMR'):
                        not_crystal = True

                    if line.count('NMR') > 0:
                        not_crystal = True
                    continue

                if record_name == 'MODEL':
                    current_model_id = line[10:14]
                    if not current_model_id in multi_model_chain_dict:
                        multi_model_chain_dict[current_model_id] = {}
                    if multi_model_mode:
                        if len(multi_model_chain_dict) * asymetric_chain_number >= len(chain_order):
                            break
                    continue
                elif record_name == 'ENDMDL':
                    if not_crystal:
                        break
                    elif not multi_model_mode:
                        multi_model_mode = True
                        asymetric_chain_number = len(chain_ids)
                    continue

                if len(res_name) < 3:
                    continue
                if len(line) > 21:
                    if record_name == "ATOM" or record_name == 'HETATM':
                        chain_id = line[21]

                        if chain_id not in chain_ids:
                            chain_ids.add(chain_id)
                            multi_model_chain_dict[current_model_id][chain_id] = chain_id
                        elif multi_model_mode:
                            if not chain_id in multi_model_chain_dict[current_model_id]:
                                new_chain_id = None
                                while new_chain_id is None:
                                    if chain_order_pos_marker == len(chain_order):
                                        print(pdb_id, filename)
                                    if not chain_order[chain_order_pos_marker] in chain_ids:
                                        new_chain_id = chain_order[chain_order_pos_marker]
                                    chain_order_pos_marker += 1
                                multi_model_chain_dict[current_model_id][chain_id] = new_chain_id
                                chain_id = new_chain_id
                            else:
                                chain_id = multi_model_chain_dict[current_model_id][chain_id]

                        res_nr = line[22:27].replace(" ", "")

                        altLoc = line[16]

                        if firstAltLoc is None and altLoc != ' ':
                            firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                        if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                            continue
                        if chain_id not in seqs:
                            seqs[chain_id] = ""
                            used_res[chain_id] = set()
                        if record_name == 'HETATM' and res_name not in residue_consts.THREE_TO_ONE:
                            continue
                        if res_nr not in used_res[chain_id]:
                            if res_name in residue_consts.THREE_TO_ONE:
                                aa = residue_consts.THREE_TO_ONE[res_name][0]
                            else:
                                aa = 'X'
                            seqs[chain_id] += aa
                            used_res[chain_id].add(res_nr)

            for chain in seqs:
                seq = seqs[chain]
                if len(seq) > 10:
                    if seq not in great_seq_map:
                        great_seq_map[seq] = [(pdb_id, chain)]
                    else:
                        great_seq_map[seq].append((pdb_id, chain))
    return great_seq_map, bio_entries, error_files


def parseBySeqres(path, pdb_id, great_seq_map):
    
    if path.endswith(".gz"):
        f = gzip.open(path, 'rb')
    else:
        f = open(path,'rb')

    page = f.read()
    f.close()

    page = page.decode('ascii')

    lines = page.split('\n')

    modres = False

    seqs = {}
    #modres_map = {}
    for line in lines:
        record_name = line[0:6].replace(" ", "")

        if record_name == "MODRES":
            print(line)

        if record_name != 'SEQRES':
            if record_name == 'ATOM':
                break
            elif modres and record_name == "HETNAM":
                #mod = line[12:15]
                #tlc = line[24:27]
                # print(mod,tlc)
                # modres_map[mod] = residue_consts.THREE_TO_ONE[tlc]"""
                print(line)
            continue

        chain = line[11]
        if chain not in seqs:
            seqs[chain] = ""
        seq = line[19:].split()
        for aa in seq:
            if len(aa) == 3:
                if aa in residue_consts.THREE_TO_ONE:
                    if residue_consts.THREE_TO_ONE[aa][0] in residue_consts.ONE_TO_THREE:
                        seqs[chain] += residue_consts.THREE_TO_ONE[aa][0]
                    else:
                        seqs[chain] += 'X'
                else:

                    modres = True
                    seqs[chain] += 'X'
            else:
                break

    """
    if modres:
        for chain in seqs:
            seq = seqs[chain]
            for mod in modres_map:
                seq = seq.replace("X"+mod,modres_map[mod])
            seqs[chain] = seq
    """

    for chain in seqs:
        seq = seqs[chain]
        if len(seq) > 10:
            if seq not in great_seq_map:
                great_seq_map[seq] = [(pdb_id, chain)]
            else:
                great_seq_map[seq].append((pdb_id, chain))


def create_custom_db_fasta(config):
    great_seq_map = {}
    parse_results = []

    parse_dump = ray.put(config.custom_db_path)

    parse_results.append(parseByAtom.remote(parse_dump, config.custom_db_path, model_db = True))

    parse_out = ray.get(parse_results)
    for seq_map, bio_entries_part, error_files in parse_out:
        if len(error_files) > 0:
            for error_file in error_files:
                config.errorlog.add_warning('Error with file: %s' % error_file)
        for seq in seq_map:
            if seq not in great_seq_map:
                great_seq_map[seq] = seq_map[seq]
            else:
                great_seq_map[seq] += seq_map[seq]

    entries = []

    id_list = []

    for seq in great_seq_map:
        if seq == '':
            continue

        for pdb, chain in great_seq_map[seq]:
            id_list.append('%s-%s' % (pdb, chain))
            #seed_id_map_lines.append('%s\t%s' % (seed_id,','.join(id_list)))

        block_seqs = []
        while len(seq) > 0:
            m = min((len(seq), 80))
            block_seqs.append(seq[:m])
            seq = seq[m:]

        if id_list == []:
            continue

        entry = ">%s\n%s" % (','.join(id_list), '\n'.join(block_seqs))
        entries.append(entry)
        id_list = []

    page = '\n'.join(entries)

    f = open(config.custom_db_fasta_path, 'w') # get set in scripts/update.py
    f.write(page)
    f.close()

    return len(great_seq_map) > 0


def main(config):
    seq_file = config.custom_db_fasta_path # get set in scripts/update.py

    ray_init(config)

    custom_path = config.custom_db_path

    divide_folder = "%s/" % custom_path

    great_seq_map = {}

    bio_entries = set()

    parse_results = []
    parse_dump = ray.put(divide_folder)

    parse_results.append(parseByAtom.remote(parse_dump, divide_folder))

    parse_out = ray.get(parse_results)
    for seq_map, bio_entries_part, error_files in parse_out:
        if len(error_files) > 0:
            for error_file in error_files:
                config.errorlog.add_warning('Error with file: %s' % error_file)
        for seq in seq_map:
            if seq not in great_seq_map:
                great_seq_map[seq] = seq_map[seq]
            else:
                great_seq_map[seq] += seq_map[seq]
        bio_entries = bio_entries | bio_entries_part

    entries = []

    id_list = []

    for seq in great_seq_map:
        if seq == '':
            continue

        for pdb, chain in great_seq_map[seq]:
            id_list.append('%s-%s' % (pdb, chain))

        block_seqs = []
        while len(seq) > 0:
            m = min((len(seq), 80))
            block_seqs.append(seq[:m])
            seq = seq[m:]

        if id_list == []:
            continue

        entry = ">%s\n%s" % (','.join(id_list), '\n'.join(block_seqs))
        entries.append(entry)
        id_list = []

    page = '\n'.join(entries)

    f = open(seq_file, 'w')
    f.write(page)
    f.close()

    ray.shutdown()


if __name__ == "__main__":

    main(seq_file, set(), fromScratch=True)  # BUG: undefined variable
