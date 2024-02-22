import gzip
import os
import subprocess
import time
import ray

from structman.base_utils.base_utils import pack, unpack, is_alphafold_model, alphafold_model_id_to_file_path
from structman.base_utils.ray_utils import ray_init
#from structman.lib.MMseqs2 import apply_mmseqs
from structman.scripts.createPdbBaDb import parseSequencesByAtom, parseSequencesByAtomAU
from structman.lib.createRINdb import get_recently_modified_structures
from structman.lib.pdbParser import getCurrentPDBID

site_thresh_1 = 0.99
site_thresh_2 = 0.7
site_thresh_1_str = str(int(site_thresh_1*100))
site_thresh_2_str = str(int(site_thresh_2*100))

full_thresh_1 = 0.9
full_thresh_2 = 0.5
full_thresh_1_str = str(int(full_thresh_1*100))
full_thresh_2_str = str(int(full_thresh_2*100))

feature_names = [
    'Backbone_RMSD_seq_id_high_weight', 'All_atom_RMSD_seq_id_high_weight', 'nof_site_residue_seq_id_high_weight', 'Site_LDDT_seq_id_high_weight',
    'Backbone_RMSD_seq_id_low_weight', 'All_atom_RMSD_seq_id_low_weight', 'nof_site_residue_seq_id_low_weight', 'Site_LDDT_seq_id_low_weight',
    'Backbone_RMSD_seq_id_greater_90', 'All_atom_RMSD_seq_id_greater_90', 'nof_site_residue_seq_id_greater_90', 'Site_LDDT_seq_id_greater_90',
    'Backbone_RMSD_seq_id_between_50_and_90', 'All_atom_RMSD_seq_id_between_50_and_90', 'nof_site_residue_seq_id_between_50_and_90', 'Site_LDDT_seq_id_between_50_and_90',
    'Backbone_RMSD_seq_id_lower_50', 'All_atom_RMSD_seq_id_lower_50', 'nof_site_residue_seq_id_lower_50', 'Site_LDDT_seq_id_lower_50',
    'Backbone_RMSD_site_id_greater_99', 'All_atom_RMSD_site_id_greater_99', 'nof_site_residue_site_id_greater_99', 'Site_LDDT_site_id_greater_99',
    'Backbone_RMSD_site_id_between_70_and_99', 'All_atom_RMSD_site_id_between_70_and_99', 'nof_site_residue_site_id_between_70_and_99', 'Site_LDDT_site_id_between_70_and_99',
    'Backbone_RMSD_site_id_lower_70', 'All_atom_RMSD_site_id_lower_70', 'nof_site_residue_site_id_lower_70', 'Site_LDDT_site_id_lower_70',
    'Backbone_RMSD_seq_id_greater_90_site_id_greater_99', 'All_atom_RMSD_seq_id_greater_90_site_id_greater_99', 'nof_site_residue_seq_id_greater_90_site_id_greater_99', 'Site_LDDT_seq_id_greater_90_site_id_greater_99',
    'Backbone_RMSD_seq_id_between_50_and_90_site_id_greater_99', 'All_atom_RMSD_seq_id_between_50_and_90_site_id_greater_99', 'nof_site_residue_seq_id_between_50_and_90_site_id_greater_99', 'Site_LDDT_seq_id_between_50_and_90_site_id_greater_99',
    'Backbone_RMSD_seq_id_lower_50_site_id_greater_99', 'All_atom_RMSD_seq_id_lower_50_site_id_greater_99', 'nof_site_residue_seq_id_lower_50_site_id_greater_99', 'Site_LDDT_seq_id_lower_50_site_id_greater_99',
    'Backbone_RMSD_seq_id_greater_90_site_id_between_70_and_99', 'All_atom_RMSD_seq_id_greater_90_site_id_between_70_and_99', 'nof_site_residue_seq_id_greater_90_site_id_between_70_and_99', 'Site_LDDT_seq_id_greater_90_site_id_between_70_and_99',
    'Backbone_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99', 'All_atom_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99', 'nof_site_residue_seq_id_between_50_and_90_site_id_between_70_and_99', 'Site_LDDT_seq_id_between_50_and_90_site_id_between_70_and_99',
    'Backbone_RMSD_seq_id_lower_50_site_id_between_70_and_99', 'All_atom_RMSD_seq_id_lower_50_site_id_between_70_and_99', 'nof_site_residue_seq_id_lower_50_site_id_between_70_and_99', 'Site_LDDT_seq_id_lower_50_site_id_between_70_and_99',
    'Backbone_RMSD_seq_id_greater_90_site_id_lower_70', 'All_atom_RMSD_seq_id_greater_90_site_id_lower_70', 'nof_site_residue_seq_id_greater_90_site_id_lower_70', 'Site_LDDT_seq_id_greater_90_site_id_lower_70',
    'Backbone_RMSD_seq_id_between_50_and_90_site_id_lower_70', 'All_atom_RMSD_seq_id_between_50_and_90_site_id_lower_70', 'nof_site_residue_seq_id_between_50_and_90_site_id_lower_70', 'Site_LDDT_seq_id_between_50_and_90_site_id_lower_70',
    'Backbone_RMSD_seq_id_lower_50_site_id_lower_70', 'All_atom_RMSD_seq_id_lower_50_site_id_lower_70', 'nof_site_residue_seq_id_lower_50_site_id_lower_70', 'Site_LDDT_seq_id_lower_50_site_id_lower_70'
]

def get_mm_search_db_path(micro_miner_db_path, k):
    microminer_search_db = f'{micro_miner_db_path}/index_k_{k}/index'
    return microminer_search_db

def get_pdb_file_path(pdb_id, config):
    if is_alphafold_model(pdb_id):
        return alphafold_model_id_to_file_path(pdb_id, config)
    pdb_id = getCurrentPDBID(pdb_id, config.pdb_path)
    path = '%s/data/biounit/PDB/divided/%s/%s.pdb' % (config.pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
    if not os.path.isfile(path):
        path = '%s/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (config.pdb_path, pdb_id[1:-1].lower(), pdb_id.lower())
    if not os.path.isfile(path):
        return None
    return path

def mm_lookup(pdb_id, config, n_sub_procs = 1):
    if pdb_id[-3:] == '_AU':
        pdb_id = pdb_id[0:4]

    mmdb_raw_entry, mmdb_aggregated_entry = get_mmdb_entry_path(pdb_id, config.microminer_db, identity = config.microminer_identity, kmer_matching_rate = config.microminer_kmer_matching_rate, fragment_length = config.microminer_fragment_length, k = config.microminer_k)
    
    if not os.path.isfile(mmdb_aggregated_entry) and config.update_microminer_active:
        if config.verbosity >= 4:
            print(f'Running microminer in mm_lookup: {pdb_id} {mmdb_raw_entry}')
        outfolder = mmdb_raw_entry.rsplit('/', 1)[0]
        microminer_search_db = get_mm_search_db_path(config.microminer_db, config.microminer_k)
        pdb_file_path = get_pdb_file_path(pdb_id, config)
        apply_microminer(config, microminer_search_db, pdb_file_path, outfolder, config.microminer_identity, kmer_matching_rate = config.microminer_kmer_matching_rate, fragment_length = config.microminer_fragment_length, n_procs = n_sub_procs)
        mm_post_procession(outfolder, mmdb_raw_entry, mmdb_aggregated_entry)

    if os.path.isfile(mmdb_aggregated_entry):
        feature_dict, feature_names = parse_processed_mmdb_entry(mmdb_aggregated_entry)
    else:
        feature_dict = {}
        feature_names = []
        
    return feature_dict, feature_names

def apply_microminer(config, microminer_search_db, infile, outfolder, identity = 0.6, kmer_matching_rate = 0.4, fragment_length = 8, n_procs = 1):
    cmds =  [
                config.microminer_exe, 'search', 
                '--identity', str(identity), 
                '--kmer_matching_rate', str(kmer_matching_rate),
                '--fragment_length', str(fragment_length),
                '-q', infile,
                '-s', microminer_search_db,
                '-o', outfolder,
                '-c', str(n_procs)
            ]
    if config.verbosity >= 3:
        print(f'Calling MicroMiner with cmds: {cmds}')

    p = subprocess.run(cmds, capture_output = True, text = True)

    if len(p.stdout) > 0:
        outlog = f'{outfolder}/microminer_outlog.log'
        f = open(outlog, 'w')
        f.write(p.stdout)
        f.close()
    if len(p.stderr) > 0:
        errlog = f'{outfolder}/microminer_errlog.log'
        f = open(errlog, 'w')
        f.write(p.stderr)
        f.close()


def get_mmdb_entry_path(pdb_id, mmdb_path, identity = 0.6, kmer_matching_rate = 0.4, fragment_length = 8, k = 7):
    if is_alphafold_model(pdb_id):
        uniprot_ac = pdb_id.split('-')[1]
        topfolder_id = uniprot_ac[-2:]
        subfolder_id = uniprot_ac[-4:]
        topfolder = f'{mmdb_path}/{topfolder_id}'
        if not os.path.isdir(topfolder):
            os.mkdir(topfolder)
        subfolder = f'{topfolder}/{subfolder_id}'
        if not os.path.isdir(subfolder):
            os.mkdir(subfolder)
        entry_folder = f'{subfolder}/{pdb_id}'
        if not os.path.isdir(entry_folder):
            os.mkdir(entry_folder)
        mmdb_entry_path_stem = f'{entry_folder}/{pdb_id}_i{int(identity*100)}_kmr_{int(kmer_matching_rate*100)}_f_{fragment_length}_k_{k}'

    else:
        pdb_id = pdb_id.lower()
        mmdb_entry_path_stem = f'{mmdb_path}/{pdb_id[1:-1]}/{pdb_id}/{pdb_id}_i{int(identity*100)}_kmr_{int(kmer_matching_rate*100)}_f_{fragment_length}_k_{k}'
    mmdb_raw_entry = f'{mmdb_entry_path_stem}.csv'
    mmdb_aggregated_entry = f'{mmdb_entry_path_stem}_processed.tsv'
    return mmdb_raw_entry, mmdb_aggregated_entry

def unpack_bioassembly(pdb_id, path, outfile_path):
    pdb_id = pdb_id.upper()

    f = gzip.open(path, 'rt')
    lines = f.readlines()
    f.close()

    lines[0] = lines[0].replace('XXXX', pdb_id)

    f = open(outfile_path, 'w')
    f.write(''.join(lines))
    f.close()
    return

def write_great_fasta(prot_seq_map, seq_file):
    lines = []
    for prot in prot_seq_map:
        for chain in prot_seq_map[prot]:
            seq = prot_seq_map[prot][chain]
            if seq == '':
                continue

            block_seqs = []
            while len(seq) > 0:
                m = min((len(seq), 80))
                block_seqs.append(seq[:m])
                seq = seq[m:]

            lines.append(f'>{prot}:{chain}\n')
            for block_seq in block_seqs:
                lines.append(f'{block_seq}\n')

    page = ''.join(lines)

    f = open(seq_file, 'w')
    f.write(page)
    f.close()

def parse_MM_results_file(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()

    results_dict = {}
    query_pdb_id = None

    for line in lines[1:]:
        words = line.split()
        if len(words) == 14:
            query_pdb_id, query_aa, query_chain, query_pos, hit_pdb_id, hit_aa, hit_chain, hit_pos, site_identity, bb_rmsd, all_rmsd, n_site_res, site_lddt, full_seq_id = words
        else:
            query_pdb_id, query_aa, query_chain, query_pos, hit_pdb_id, hit_aa, site_identity, bb_rmsd, all_rmsd, n_site_res, site_lddt, full_seq_id = words
            hit_chain = 'gap'
            hit_pos = 'gap'
        if (query_chain, query_pos) not in results_dict:
            results_dict[(query_chain, query_pos)] = [query_aa, {}]
        if hit_pdb_id not in results_dict[(query_chain, query_pos)][1]:
            results_dict[(query_chain, query_pos)][1][hit_pdb_id] = {}
        if (hit_chain, hit_pos) not in results_dict[(query_chain, query_pos)][1][hit_pdb_id]:
            results_dict[(query_chain, query_pos)][1][hit_pdb_id][(hit_chain, hit_pos)] = (hit_aa, float(site_identity), float(bb_rmsd), float(all_rmsd), int(n_site_res), float(site_lddt), float(full_seq_id))

    return query_pdb_id, results_dict

def aggregate_multi_mers(mm_dict):
    for qchain, qpos in mm_dict:
        for hit_pdb_id in list(mm_dict[(qchain, qpos)][1].keys()):
            if len(mm_dict[(qchain, qpos)][1][hit_pdb_id]) < 2:
                continue
            max_site_identity = 0
            max_p = None
            for (hit_chain, hit_pos) in mm_dict[(qchain, qpos)][1][hit_pdb_id]:
                site_identity = mm_dict[(qchain, qpos)][1][hit_pdb_id][(hit_chain, hit_pos)][1]
                if site_identity > max_site_identity:
                    max_site_identity = site_identity
                    max_p = (hit_chain, hit_pos)

            mm_dict[(qchain, qpos)][1][hit_pdb_id] = {max_p : mm_dict[(qchain, qpos)][1][hit_pdb_id][max_p]}

def add_to_sums(sums, weight, bb_rmsd, all_rmsd, n_site_res, site_lddt):
    sums[0] += weight
    sums[1] += bb_rmsd * weight
    sums[2] += all_rmsd * weight
    sums[3] += n_site_res * weight
    sums[4] += site_lddt * weight

def weight_sums(sums):
    weighted_sums = []
    if sums[0] == 0:
        return [0. , 0., 0., 0.]
    weighted_sums.append(sums[1]/sums[0])
    weighted_sums.append(sums[2]/sums[0])
    weighted_sums.append(sums[3]/sums[0])
    weighted_sums.append(sums[4]/sums[0])
    return weighted_sums

def aggregate_mm_dict(mm_dict, query_pdb_id):
    aggregated_dict = {}
    for qchain, qpos in mm_dict:
        seq_id_high_weight_sums = [0., 0. , 0., 0., 0.]
        seq_id_low_weight_sums = [0., 0. , 0., 0., 0.]

        seq_id_greater_90_sums = [0., 0. , 0., 0., 0.]
        seq_id_between_50_and_90_sums = [0., 0. , 0., 0., 0.]
        seq_id_lower_50_sums = [0., 0. , 0., 0., 0.]

        site_id_greater_90_sums = [0., 0. , 0., 0., 0.]
        site_id_between_50_and_90_sums = [0., 0. , 0., 0., 0.]
        site_id_lower_50_sums = [0., 0. , 0., 0., 0.]

        seq_id_greater_90_site_id_greater_90_sums = [0., 0. , 0., 0., 0.]
        seq_id_between_50_and_90_site_id_greater_90_sums = [0., 0. , 0., 0., 0.]
        seq_id_lower_50_site_id_greater_90_sums = [0., 0. , 0., 0., 0.]

        seq_id_greater_90_site_id_between_50_and_90_sums = [0., 0. , 0., 0., 0.]
        seq_id_between_50_and_90_site_id_between_50_and_90_sums = [0., 0. , 0., 0., 0.]
        seq_id_lower_50_site_id_between_50_and_90_sums = [0., 0. , 0., 0., 0.]

        seq_id_greater_90_site_id_lower_50_sums = [0., 0. , 0., 0., 0.]
        seq_id_between_50_and_90_site_id_lower_50_sums = [0., 0. , 0., 0., 0.]
        seq_id_lower_50_site_id_lower_50_sums = [0., 0. , 0., 0., 0.]
        for hit_pdb_id in mm_dict[(qchain, qpos)][1]:
            if query_pdb_id == hit_pdb_id:
                continue
            for (hit_chain, hit_pos) in mm_dict[(qchain, qpos)][1][hit_pdb_id]:
                (hit_aa, site_identity, bb_rmsd, all_rmsd, n_site_res, site_lddt, full_seq_id) = mm_dict[(qchain, qpos)][1][hit_pdb_id][(hit_chain, hit_pos)]
                seq_id_high_weight = site_identity * full_seq_id
                add_to_sums(seq_id_high_weight_sums, seq_id_high_weight, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                seq_id_low_weight = (1.5-site_identity) * (1.5-full_seq_id)
                add_to_sums(seq_id_low_weight_sums, seq_id_low_weight, bb_rmsd, all_rmsd, n_site_res, site_lddt)

                if full_seq_id >= full_thresh_1:
                    add_to_sums(seq_id_greater_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    if site_identity >= site_thresh_1:
                        add_to_sums(seq_id_greater_90_site_id_greater_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    elif site_identity >= site_thresh_2:
                        add_to_sums(seq_id_greater_90_site_id_between_50_and_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    else:
                        add_to_sums(seq_id_greater_90_site_id_lower_50_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                elif full_seq_id >= full_thresh_2:
                    add_to_sums(seq_id_between_50_and_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    if site_identity >= site_thresh_1:
                        add_to_sums(seq_id_between_50_and_90_site_id_greater_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    elif site_identity >= site_thresh_2:
                        add_to_sums(seq_id_between_50_and_90_site_id_between_50_and_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    else:
                        add_to_sums(seq_id_between_50_and_90_site_id_lower_50_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                else:
                    add_to_sums(seq_id_lower_50_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    if site_identity >= site_thresh_1:
                        add_to_sums(seq_id_lower_50_site_id_greater_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    elif site_identity >= site_thresh_2:
                        add_to_sums(seq_id_lower_50_site_id_between_50_and_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                    else:
                        add_to_sums(seq_id_lower_50_site_id_lower_50_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                if site_identity >= site_thresh_1:
                    add_to_sums(site_id_greater_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                elif site_identity >= site_thresh_2:
                    add_to_sums(site_id_between_50_and_90_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)
                else:
                    add_to_sums(site_id_lower_50_sums, 1, bb_rmsd, all_rmsd, n_site_res, site_lddt)

        weighted_seq_id_high_weight = weight_sums(seq_id_high_weight_sums)
        weighted_seq_id_low_weight = weight_sums(seq_id_low_weight_sums)

        weighted_seq_id_greater_90 = weight_sums(seq_id_greater_90_sums)
        weighted_seq_id_between_50_and_90 = weight_sums(seq_id_between_50_and_90_sums)
        weighted_seq_id_lower_50 = weight_sums(seq_id_lower_50_sums)

        weighted_site_id_greater_90 = weight_sums(site_id_greater_90_sums)
        weighted_site_id_between_50_and_90 = weight_sums(site_id_between_50_and_90_sums)
        weighted_site_id_lower_50 = weight_sums(site_id_lower_50_sums)

        weighted_seq_id_greater_90_site_id_greater_90 = weight_sums(seq_id_greater_90_site_id_greater_90_sums)
        weighted_seq_id_between_50_and_90_site_id_greater_90 = weight_sums(seq_id_between_50_and_90_site_id_greater_90_sums)
        weighted_seq_id_lower_50_site_id_greater_90 = weight_sums(seq_id_lower_50_site_id_greater_90_sums)

        weighted_seq_id_greater_90_site_id_between_50_and_90 = weight_sums(seq_id_greater_90_site_id_between_50_and_90_sums)
        weighted_seq_id_between_50_and_90_site_id_between_50_and_90 = weight_sums(seq_id_between_50_and_90_site_id_between_50_and_90_sums)
        weighted_seq_id_lower_50_site_id_between_50_and_90 = weight_sums(seq_id_lower_50_site_id_between_50_and_90_sums)

        weighted_seq_id_greater_90_site_id_lower_50 = weight_sums(seq_id_greater_90_site_id_lower_50_sums)
        weighted_seq_id_between_50_and_90_site_id_lower_50 = weight_sums(seq_id_between_50_and_90_site_id_lower_50_sums)
        weighted_seq_id_lower_50_site_id_lower_50 = weight_sums(seq_id_lower_50_site_id_lower_50_sums)

        aggregated_dict[(qchain, qpos)] = [
            weighted_seq_id_high_weight,
            weighted_seq_id_low_weight,
            weighted_seq_id_greater_90,
            weighted_seq_id_between_50_and_90,
            weighted_seq_id_lower_50,
            weighted_site_id_greater_90,
            weighted_site_id_between_50_and_90,
            weighted_site_id_lower_50,
            weighted_seq_id_greater_90_site_id_greater_90,
            weighted_seq_id_between_50_and_90_site_id_greater_90,
            weighted_seq_id_lower_50_site_id_greater_90,
            weighted_seq_id_greater_90_site_id_between_50_and_90,
            weighted_seq_id_between_50_and_90_site_id_between_50_and_90,
            weighted_seq_id_lower_50_site_id_between_50_and_90,
            weighted_seq_id_greater_90_site_id_lower_50,
            weighted_seq_id_between_50_and_90_site_id_lower_50,
            weighted_seq_id_lower_50_site_id_lower_50
        ]
    return aggregated_dict

def write_aggregated_dict(aggregated_dict, outfile):
    headers = ['Chain', 'Pos']
    base_strings = [
        'seq_id_high_weight',
        'seq_id_low_weight',
        f'seq_id_greater_{full_thresh_1_str}',
        f'seq_id_between_{full_thresh_2_str}_and_{full_thresh_1_str}',
        f'seq_id_lower_{full_thresh_2_str}',
        f'site_id_greater_{site_thresh_1_str}',
        f'site_id_between_{site_thresh_2_str}_and_{site_thresh_1_str}',
        f'site_id_lower_{site_thresh_2_str}',
        f'seq_id_greater_{full_thresh_1_str}_site_id_greater_{site_thresh_1_str}',
        f'seq_id_between_{full_thresh_2_str}_and_{full_thresh_1_str}_site_id_greater_{site_thresh_1_str}',
        f'seq_id_lower_{full_thresh_2_str}_site_id_greater_{site_thresh_1_str}',
        f'seq_id_greater_{full_thresh_1_str}_site_id_between_{site_thresh_2_str}_and_{site_thresh_1_str}',
        f'seq_id_between_{full_thresh_2_str}_and_{full_thresh_1_str}_site_id_between_{site_thresh_2_str}_and_{site_thresh_1_str}',
        f'seq_id_lower_{full_thresh_2_str}_site_id_between_{site_thresh_2_str}_and_{site_thresh_1_str}',
        f'seq_id_greater_{full_thresh_1_str}_site_id_lower_{site_thresh_2_str}',
        f'seq_id_between_{full_thresh_2_str}_and_{full_thresh_1_str}_site_id_lower_{site_thresh_2_str}',
        f'seq_id_lower_{full_thresh_2_str}_site_id_lower_{site_thresh_2_str}'
    ]
    feat_strs = [
        'Backbone RMSD',
        'All atom RMSD',
        'nof_site_residue',
        'Site LDDT'
    ]
    for bs in base_strings:
        for fs in feat_strs:
            headers.append(f'{fs} {bs}')
    header = '\t'.join(headers) + '\n'

    lines = [header]

    for (qchain, qpos) in aggregated_dict:
        words = [qchain, qpos]
        for feature_set in aggregated_dict[(qchain, qpos)]:
            for feat in feature_set:
                words.append(str(feat))
        lines.append('\t'.join(words) + '\n')

    f = open(outfile, 'w')
    f.write(''.join(lines))
    f.close()

def parse_processed_mmdb_entry(mmdb_aggregated_entry):

    f = open(mmdb_aggregated_entry, 'r')
    lines = f.readlines()
    f.close()

    feature_dict = {}

    headers = lines[0][:-1].split('\t')

    feature_names = headers[2:]
    feature_type_code_list = []
    for featname in feature_names:
        if featname[:4] == 'nof_':
            feature_type_code_list.append(1)
        else:
            feature_type_code_list.append(0)

    for line in lines[1:]:
        words = line[:-1].split('\t')
        chain_id = words[0]
        res_id = words[1]
        feature_values = words[2:]

        feature_dict[(chain_id, res_id)] = []
        for pos, feat_val_str in enumerate(feature_values):
            if feature_type_code_list[pos] == 0:
                try:
                    feat_val = float(feat_val_str)
                except:
                    feat_val = None
                if not feat_val == feat_val:
                    feat_val = None
                elif feat_val == 0.:
                    feat_val = None
            else:
                try:
                    feat_val = int(feat_val_str)
                except:
                    try:
                        feat_val = round(float(feat_val_str))
                    except:
                        feat_val = None
                if feat_val == 0:
                    feat_val = None
            feature_dict[(chain_id, res_id)].append(feat_val)

    return feature_dict, feature_names
    


@ray.remote(max_calls = 1)
def calc_microminer_package(package, store):
    package = unpack(package)

    config, microminer_search_db, identity, kmer_matching_rate, fragment_length, k, n_sub_procs, fromScratch = store

    for infile, pdb_id in package:
        mmdb_raw_entry, mmdb_aggregated_entry = get_mmdb_entry_path(pdb_id, config.microminer_db, identity = identity, kmer_matching_rate = kmer_matching_rate, fragment_length = fragment_length, k = k)
        if not fromScratch and os.path.isfile(mmdb_aggregated_entry):
            continue
        outfolder = mmdb_raw_entry.rsplit('/', 1)[0]
        apply_microminer(config, microminer_search_db, infile, outfolder, identity = identity, kmer_matching_rate = kmer_matching_rate, fragment_length = fragment_length, n_procs = n_sub_procs)
        mm_post_procession(outfolder, mmdb_raw_entry, mmdb_aggregated_entry)

def mm_post_procession(outfolder, mmdb_raw_entry, mmdb_aggregated_entry):
    try:
        os.rename(f'{outfolder}/resultStatistic.csv', mmdb_raw_entry)
    except:
        pass
    query_pdb_id, mm_dict = parse_MM_results_file(mmdb_raw_entry)
    aggregate_multi_mers(mm_dict)
    aggregated_dict = aggregate_mm_dict(mm_dict, query_pdb_id)
    write_aggregated_dict(aggregated_dict, mmdb_aggregated_entry)

@ray.remote(max_calls = 1)
def parse_subfolders(package, store):
    bio_assembly_path, config, recently_modified_structures, fromScratch, identity, kmer_matching_rate, fragment_length, k = store
    results = []
    bio_pdbs = set()
    for subfolder in package:
        sub_path = "%s/%s" % (bio_assembly_path, subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (config.microminer_db, subfolder)):
            os.mkdir("%s/%s" % (config.microminer_db, subfolder))

        for fn in files:
            if fn.count('.pdb1.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path, fn)
                if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                    continue

                pdb_id = fn.replace('.pdb1.gz', '')
                bio_pdbs.add(pdb_id)

                if recently_modified_structures is not None:
                    if not pdb_id in recently_modified_structures:
                        continue

                _, mmdb_aggregated_entry = get_mmdb_entry_path(pdb_id, config.microminer_db, identity = identity, kmer_matching_rate = kmer_matching_rate, fragment_length = fragment_length, k = k)

                if not fromScratch and os.path.exists(mmdb_aggregated_entry):
                    continue

                if not os.path.exists(f"{config.microminer_db}/{subfolder}/{pdb_id}"):
                    os.mkdir(f"{config.microminer_db}/{subfolder}/{pdb_id}")

                unpacked_path = pdbgz_path[:-4]
                if not os.path.exists(unpacked_path):
                    unpack_bioassembly(pdb_id, pdbgz_path, unpacked_path)

                seqs, _, bio_pdbs, _ = parseSequencesByAtom(full_path = pdbgz_path, bio_entries = bio_pdbs)

                seq_len_map = {}
                for chain in seqs:
                    seq_len_map[chain] = len(seqs[chain])

                results.append((unpacked_path, pdb_id, seq_len_map))
    return pack((results, bio_pdbs))

def update_microminer_db(config, fromScratch = False):
    t0 = time.time()
    if config.verbosity >= 2:
        print("Starting update Microminer DB routine")

    ray_init(config)

    bio_assembly_path = f'{config.pdb_path}/data/biounit/PDB/divided'
    AU_path = f'{config.pdb_path}/data/structures/divided/pdb'
    status_path = f'{config.pdb_path}/data/status'

    meta_file_path = f'{config.microminer_db}/meta.txt'

    identity = config.microminer_identity
    kmer_matching_rate = config.microminer_kmer_matching_rate
    fragment_length = config.microminer_fragment_length
    k = config.microminer_k

    microminer_search_db = get_mm_search_db_path(config.microminer_db, k)

    if os.path.isfile(meta_file_path):
        f = open(meta_file_path, 'r')
        time_of_last_update = f.read()
        f.close()
    else:
        time_of_last_update = None

    f = open(meta_file_path,'w')
    f.write(time.asctime(time.gmtime()))
    f.close()

    if not fromScratch:
        #recently_modified_structures = get_recently_modified_structures(time_of_last_update, status_path)
        recently_modified_structures = None
    else:
        recently_modified_structures = None

    bio_pdbs = set()
    subfolders = os.listdir(bio_assembly_path)

    N = 0
    to_process = []
    seq_len_map = {}

    #prot_seq_map = {}
    #small_proteins = set()

    packages = []
    current_package = 0

    parse_process_ids = []

    for subfolder in subfolders:

        if current_package == len(packages):
            packages.append([])

        packages[current_package].append(subfolder)


        sub_path = "%s/%s" % (bio_assembly_path, subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (config.microminer_db, subfolder)):
            os.mkdir("%s/%s" % (config.microminer_db, subfolder))

        for fn in files:
            if fn.count('.pdb1.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path, fn)
                if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                    continue

                pdb_id = fn.replace('.pdb1.gz', '')
                bio_pdbs.add(pdb_id)

                if recently_modified_structures is not None:
                    if not pdb_id in recently_modified_structures:
                        continue

                _, mmdb_aggregated_entry = get_mmdb_entry_path(pdb_id, config.microminer_db, identity = identity, kmer_matching_rate = kmer_matching_rate, fragment_length = fragment_length, k = k)

                if not fromScratch and os.path.exists(mmdb_aggregated_entry):
                    continue

                unpacked_path = pdbgz_path[:-4]
                if not os.path.exists(unpacked_path):
                    unpack_bioassembly(pdb_id, pdbgz_path, unpacked_path)

                seqs, _, bio_pdbs, _ = parseSequencesByAtom(full_path = pdbgz_path, bio_entries = bio_pdbs)

                seq_len_map[pdb_id] = {}
                for chain in seqs:
                    #if len(seqs[chain]) < 100:
                    #    small_proteins.add(f'{pdb_id}:{chain}')
                    seq_len_map[pdb_id][chain] = len(seqs[chain])

                to_process.append((unpacked_path, pdb_id))
                N += 1


        current_package += 1
        if current_package >= config.proc_n:
            current_package = 0

    t1 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 1: {t1-t0}')

    store = ray.put((bio_assembly_path, config, recently_modified_structures, fromScratch, identity, kmer_matching_rate, fragment_length, k))

    for package in packages:
        parse_process_ids.append(parse_subfolders.remote(package, store))

    t2 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 2: {t2-t1}')

    result_list = ray.get(parse_process_ids)

    t3 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 3: {t3-t2}')

    for packaged_result in result_list:
        results, bio_pdbs_subset = unpack(packaged_result)
        bio_pdbs = bio_pdbs.union(bio_pdbs_subset)
        for unpacked_path, pdb_id, seq_len_sub_map in results:
            to_process.append((unpacked_path, pdb_id))
            seq_len_map[pdb_id] = seq_len_sub_map

    t4 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 4: {t4-t3}')

    subfolders = os.listdir(AU_path)
    for subfolder in subfolders:

        sub_path = "%s/%s" % (AU_path, subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists(f"{config.microminer_db}/{subfolder}"):
            os.mkdir(f"{config.microminer_db}/{subfolder}")

        for fn in files:
            if fn.count('.ent.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path, fn)
                if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                    continue
                pdb_id = fn[3:7]

                if recently_modified_structures is not None:
                    if not pdb_id in recently_modified_structures:
                        continue

                if pdb_id in bio_pdbs:
                    continue

                _, mmdb_aggregated_entry = get_mmdb_entry_path(pdb_id, config.microminer_db, identity = identity, kmer_matching_rate = kmer_matching_rate, fragment_length = fragment_length, k = k)

                if not fromScratch and os.path.exists(mmdb_aggregated_entry):
                    continue

                if not os.path.exists(f"{config.microminer_db}/{subfolder}/{pdb_id}"):
                    os.mkdir(f"{config.microminer_db}/{subfolder}/{pdb_id}")

                seqs, _ = parseSequencesByAtomAU(full_path = pdbgz_path)

                #prot_seq_map[pdb_id] = seqs
                seq_len_map[pdb_id] = {}
                for chain in seqs:
                    #if len(seqs[chain]) < 100:
                    #    small_proteins.add(f'{pdb_id}:{chain}')
                    seq_len_map[pdb_id][chain] = len(seqs[chain])

                to_process.append((pdbgz_path, pdb_id))

    t5 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 5: {t5-t4}')

    #seq_file = f'{config.microminer_db}/search_sequences.fasta'

    #write_great_fasta(prot_seq_map, seq_file)

    t6 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 6: {t6-t5}')

    """
    mmseqs_tmp_folder = '%s/mmseqs_tmp' % config.mmseqs_tmp_folder
    if not os.path.exists(mmseqs_tmp_folder):
        os.mkdir(mmseqs_tmp_folder)

    hits, pdb_ids, _ = apply_mmseqs(mmseqs_tmp_folder, config.mmseqs2_path, seq_file, config.mmseqs2_db_path, config.gigs_of_ram, config.errorlog, config.option_seq_thresh, small_proteins = small_proteins, verbosity = config.verbosity)
    """
    t7 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 7: {t7-t6}')

        
    cost_list = []

    for path, pdb_id in to_process:
        total_cost = 0
        """
        for chain in prot_seq_map[pdb_id]:
            seq_len = len(prot_seq_map[pdb_id][chain])
            pdb_tuple = f'{pdb_id}:{chain}'
            if pdb_tuple in hits:
                n_of_hits = len(hits[pdb_tuple])
            else:
                n_of_hits = 1
            total_cost += seq_len * n_of_hits
        """
        for chain in seq_len_map[pdb_id]:
            total_cost += seq_len_map[pdb_id][chain]
        if total_cost == 0:
            continue
        cost_list.append((total_cost, path, pdb_id))

    t8 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 8: {t8-t7}')

    cost_list = sorted(cost_list, key=lambda x: x[0], reverse = True)

    t9 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 9: {t9-t8}')

    #cost_list = cost_list[-10000:-9950]
    #print(cost_list)

    packages = []
    current_package = 0
    n_sub_procs = 4
    n_of_packages = config.proc_n // n_sub_procs
    n = 0

    for cost, path, pdb_id in cost_list:
        if current_package == len(packages):
            packages.append([])
        
        packages[current_package].append((path, pdb_id))

        current_package += 1
        if current_package >= n_of_packages:
            current_package = 0

        if n % 1000 == 0:
            print(f'Structure number: {n} -> Size: {cost}')
        n += 1
    
    t10 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 10: {t10-t9}')

    microminer_processes = []

    store = ray.put((config, microminer_search_db, identity, kmer_matching_rate, fragment_length, k, n_sub_procs, fromScratch))

    t11 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 11: {t11-t10}')

    for package in packages:
        microminer_processes.append(calc_microminer_package.remote(pack(package), store))

    t12 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 12: {t12-t11}')

    ray.get(microminer_processes)

    t13 = time.time()
    if config.verbosity >= 2:
        print(f'Time for updating MMDB part 13: {t13-t12}')