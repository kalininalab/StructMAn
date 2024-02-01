import gzip
import os
import sys

# if running as script, add local structman package to path
if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__))))))
from structman import settings


def parseByAtom(path, pdb_id, great_seq_map):
    f = gzip.open(path, 'rb')
    page = f.read()
    f.close()
    lines = page.split('\n')

    seqs = {}
    used_res = {}
    for line in lines:
        record_name = line[0:6].replace(" ", "")
        atom_nr = line[6:11].replace(" ", "")
        atom_name = line[12:16].replace(" ", "")
        res_name = line[17:20].replace(" ", "")
        if record_name == 'ENDMDL':
            break
        if len(res_name) != 1:
            continue
        if len(line) > 21:
            chain_id = line[21]
            res_nr = line[22:27].replace(" ", "")

            if record_name == "ATOM" or record_name == 'HETATM':
                if chain_id not in seqs:
                    seqs[chain_id] = ""
                    used_res[chain_id] = set()
                if record_name == 'HETATM':
                    continue
                if res_nr not in used_res[chain_id]:

                    seqs[chain_id] += res_name
                    used_res[chain_id].add(res_nr)

    for chain in seqs:
        seq = seqs[chain]
        if len(seq) > 10:
            if seq not in great_seq_map:
                great_seq_map[seq] = [(pdb_id, chain)]
            else:
                great_seq_map[seq].append((pdb_id, chain))


def parse_seq_file(seq_file):
    f = open(seq_file, 'r')
    lines = f.read().split('\n')
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
                pdb_ids.add(pdb_id)
                if not os.path.isfile('%s/data/structures/obsolete/pdb/%s/pdb%s.ent.gz' % (settings.PDB_DIR, pdb_id[1:-1].lower(), pdb_id.lower())):
                    chains.append((pdb_id, chain))

    great_seq_map[seq] = chains

    return pdb_ids, great_seq_map


def main(seq_file, fromScratch=False):
    global sub_folders
    global AU_sub_folders
    global divide_folder

    if not fromScratch:
        safe_list, great_seq_map = parse_seq_file(seq_file)  # update this
    else:
        safe_list = set()
        great_seq_map = {}

    sub_folders = sorted(sub_folders)
    bio_entries = set()

    for sub_folder in sub_folders:
        files = os.listdir("%s%s" % (divide_folder, sub_folder))

        print(sub_folder)
        for filename in files:
            if filename.count('.pdb1.gz') > 0:
                pdb_id = filename[:4].upper()
                if pdb_id in safe_list:
                    continue
                bio_entries.add(pdb_id)
                path = "%s%s/%s" % (divide_folder, sub_folder, filename)
                parseByAtom(path, pdb_id, great_seq_map)
                # parseBySeqres(path,pdb_id,great_seq_map)

    sub_folders = sorted(AU_sub_folders)
    for sub_folder in sub_folders:
        files = os.listdir("%s%s" % (AU_folder, sub_folder))
        print(sub_folder, ' AU')
        for filename in files:
            if filename.count('.ent.gz') > 0:
                pdb_id = filename[3:7].upper()
                if pdb_id in bio_entries or pdb_id in safe_list:
                    continue
                pdb_id = '%s_AU' % pdb_id
                if pdb_id in safe_list:
                    continue
                path = "%s%s/%s" % (AU_folder, sub_folder, filename)
                parseByAtom(path, pdb_id, great_seq_map)

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

    #f = open(id_map_file,'w')
    # f.write('\n'.join(seed_id_map_lines))
    # f.close()


if __name__ == "__main__":
    seq_file = '/TL/sin/work/agress/pdbba_nuc'

    divide_folder = "%s/data/biounit/PDB/divided/" % settings.PDB_DIR

    sub_folders = os.listdir(divide_folder)

    AU_folder = "%s/data/structures/divided/pdb/" % settings.PDB_DIR

    AU_sub_folders = os.listdir(AU_folder)

    main(seq_file, fromScratch=True)
