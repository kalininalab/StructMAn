import gzip
import os
import sys
import subprocess
import ray
import traceback

# if running as script, add local structman package to path
if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__))))))

from structman.lib.sdsc.consts import residues as residue_consts
from structman.base_utils.ray_utils import ray_init, ray_hack

def parseSequencesByAtom(full_path = None, divide_folder = None, sub_folder = None, filename = None, bio_entries = None, model_db = False, error_files = None):
    if bio_entries is None:
        bio_entries = set()

    if error_files is None:
        error_files = []

    if full_path is None:
        path = f'{divide_folder}/{sub_folder}/{filename}'
    else:
        filename = full_path.rsplit('/', 1)[-1]
        path = full_path

    if not model_db:
        if not filename.count('.pdb1.gz') > 0:
            return {}, None, bio_entries, error_files
        pdb_id = filename[:4].upper()

        bio_entries.add(pdb_id)

        f = gzip.open(path, 'rb')
    else:
        if len(filename) < 7:
            return {}, None, bio_entries, error_files
        if filename[-6:] != 'pdb.gz':
            return {}, None, bio_entries, error_files
        pdb_id = '-'.join(filename.split('-')[:2])
        path = f'{divide_folder}/{sub_folder}/{filename}'

        f = gzip.open(path,'rb')
    page = f.read()
    f.close()

    try:
        page = page.decode('ascii')
    except:
        try:
            decoded_lines = []
            for line in page.split(b'\n'):
                try:
                    decoded_lines.append(f"{line.decode('ascii')}\n")
                except:
                    pass
            page = ''.join(decoded_lines)
        except:
            [e, f, g] = sys.exc_info()
            try:
                g = traceback.format_exc(g)
            except:
                pass
            error = f"Couldn't decode page: {path}\n{e}\n{f}\n{g}"
            error_files.append(error)
            return {}, pdb_id, bio_entries, error_files

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
    chain_order = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz'

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
    return seqs, pdb_id, bio_entries, error_files

def parseSequencesByAtomAU(full_path = None, AU_folder = None, sub_folder = None, filename = None):

    if full_path is None:
        path = f'{AU_folder}/{sub_folder}/{filename}'
    else:
        filename = full_path.rsplit('/', 1)[-1]
        path = full_path

    if not filename.count('.ent.gz') > 0:
        return {}, None
    
    pdb_id = filename[3:7].upper()

    pdb_id = '%s_AU' % pdb_id

    f = gzip.open(path, 'rb')
    page = f.read()
    f.close()

    try:
        page = page.decode('ascii')
    except:
        return {}, None

    lines = page.split('\n')

    seqs = {}
    used_res = {}
    firstAltLoc = None
    for line in lines:
        record_name = line[0:6].replace(" ", "")
        atom_nr = line[6:11].replace(" ", "")
        atom_name = line[12:16].replace(" ", "")
        res_name = line[17:20].replace(" ", "")
        if record_name == 'ENDMDL':
            break
        if len(res_name) < 3:
            continue
        if len(line) > 21:

            chain_id = line[21]
            res_nr = line[22:27].replace(" ", "")

            if record_name == "ATOM" or record_name == 'HETATM':
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
    return seqs, pdb_id

@ray.remote(num_cpus=1)
def parseByAtom(parse_dump, sub_folder, model_db = False):
    ray_hack()
    
    divide_folder = parse_dump

    files = os.listdir("%s/%s" % (divide_folder, sub_folder))
    great_seq_map = {}
    bio_entries = set()

    error_files = []

    for filename in files:
        seqs, pdb_id, bio_entries, error_files = parseSequencesByAtom(divide_folder = divide_folder, sub_folder = sub_folder, filename = filename, bio_entries = bio_entries, model_db = model_db, error_files = error_files)

        if len(seqs) == 0:
            continue

        for chain in seqs:
            seq = seqs[chain]
            if len(seq) > 10:
                if seq not in great_seq_map:
                    great_seq_map[seq] = [(pdb_id, chain)]
                else:
                    great_seq_map[seq].append((pdb_id, chain))
    return great_seq_map, bio_entries, error_files


@ray.remote(num_cpus=1)
def parseByAtomAU(sub_folder, AU_dump):
    ray_hack()

    AU_folder, bio_entries = AU_dump
    great_seq_map = {}
    files = os.listdir("%s/%s" % (AU_folder, sub_folder))
    for filename in files:
        seqs, pdb_id = parseSequencesByAtomAU(AU_folder = AU_folder, sub_folder = sub_folder, filename = filename)
        
        if len(seqs) == 0:
            continue

        for chain in seqs:
            seq = seqs[chain]
            if len(seq) > 10:
                if seq not in great_seq_map:
                    great_seq_map[seq] = [(pdb_id, chain)]
                else:
                    great_seq_map[seq].append((pdb_id, chain))
    return great_seq_map


def parseBySeqres(path, pdb_id, great_seq_map):
    f = gzip.open(path, 'rb')
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
                    print(pdb_id)
                    print(aa)
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


def parse_seq_file(seq_file, recently_modified_structures):
    f = open(seq_file, 'rb')
    lines = f.read().decode('ascii').split('\n')
    f.close()

    pdb_ids = set()

    great_seq_map = {}

    seq = ""

    chains = []

    for line in lines:
        if line[0] != '>':
            seq += line
        else:
            if not seq == '':
                great_seq_map[seq] = chains
                seq = ''
                chains = []
            parts = line[1:].split(',')
            for part in parts:
                if part.count('-') < 1:
                    continue
                pdb_id, chain = part.split('-')
                if pdb_id.lower() in recently_modified_structures:
                    continue
                pdb_ids.add(pdb_id)
                if not os.path.isfile('%s/data/structures/obsolete/pdb/%s/pdb%s.ent.gz' % (pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())):  # BUG: undefined variable
                    chains.append((pdb_id, chain))

    great_seq_map[seq] = chains

    return pdb_ids, great_seq_map

def create_model_db_fasta(config):
    great_seq_map = {}
    parse_results = []

    parse_dump = ray.put(config.path_to_model_db)

    for topfolder in os.listdir(config.path_to_model_db):
        if not os.path.isdir(f'{config.path_to_model_db}/{topfolder}'):
            continue
        for sub_folder in os.listdir(f'{config.path_to_model_db}/{topfolder}'):
            parse_results.append(parseByAtom.remote(parse_dump, f'{topfolder}/{sub_folder}', model_db = True))

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

    f = open(config.model_db_fasta_path, 'w') # get set in scripts/update.py
    f.write(page)
    f.close()

    return len(great_seq_map) > 0

def main(config, update_source = None):
    seq_file = config.pdb_fasta_path # get set in scripts/update.py

    ray_init(config)

    print('After ray init')
    pdb_path = config.pdb_path

    divide_folder = "%s/data/biounit/PDB/divided/" % pdb_path

    sub_folders = os.listdir(divide_folder)

    AU_folder = "%s/data/structures/divided/pdb/" % pdb_path

    AU_sub_folders = os.listdir(AU_folder)

    # if not fromScratch:
    #    safe_list,great_seq_map = parse_seq_file(seq_file,recently_modified_structures)
    # else:
    #safe_list = set()
    great_seq_map = {}

    sub_folders = sorted(sub_folders)
    bio_entries = set()

    parse_results = []
    parse_dump = ray.put(divide_folder)

    for sub_folder in sub_folders:
        parse_results.append(parseByAtom.remote(parse_dump, sub_folder))

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

    AU_dump = ray.put((AU_folder, bio_entries))
    parse_results = []

    sub_folders = sorted(AU_sub_folders)
    for sub_folder in sub_folders:
        parse_results.append(parseByAtomAU.remote(sub_folder, AU_dump))

    parse_out = ray.get(parse_results)
    for seq_map in parse_out:
        for seq in seq_map:
            if seq not in great_seq_map:
                great_seq_map[seq] = seq_map[seq]
            else:
                great_seq_map[seq] += seq_map[seq]

    entries = []

    #seed_id_map_lines = []

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

    f = open(seq_file, 'w')
    f.write(page)
    f.close()

    if (config.path_to_model_db is not None) and (config.path_to_model_db != ''):
        config.model_db_fasta_created = create_model_db_fasta(config)
    else:
        config.model_db_fasta_created = False

    if update_source is not None:
        target_folder = f'{update_source}/structman/resources'
        with open(f'{target_folder}/pdbba.fasta.gz', 'w') as outfile:
            p = subprocess.run(['gzip', '-9', '-c', seq_file], stdout=outfile)
            
        if config.model_db_fasta_created:
            with open(f'{target_folder}/modeldb.fasta.gz', 'w') as outfile:
                p = subprocess.run(['gzip', '-9', '-c', config.model_db_fasta_path], stdout=outfile)             

    ray.shutdown()

    #f = open(id_map_file,'w')
    # f.write('\n'.join(seed_id_map_lines))
    # f.close()
