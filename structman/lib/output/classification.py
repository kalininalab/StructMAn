import os
import time
import sys
import traceback
from typing import IO

import ray
#from numba import jit, njit
#import numba as nb

from structman import settings
from structman.lib import rin
from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc import mappings as mappings_package
from structman.lib.sdsc.consts import residues as residue_consts
from structman.lib.database import database
from structman.lib.database.database_core_functions import select, binningSelect
from structman.lib.output.feature import init_feature_table
from structman.lib.output.interfaces import init_aggregated_interface_table, init_protein_protein_table, init_position_position_table
from structman.lib.output import out_generator
from filelock import FileLock
from structman.base_utils import base_utils
from structman.base_utils import ray_utils
from structman.lib.output import out_utils

def writeStatFile(out_file, mutation_dict, class_dict, tag_map, stat_dict=None):
    class Count_helper:
        def __init__(self, tag_name):
            self.tag_name = tag_name
            self.prot_type_map = {}
            self.total_positions = 0
            self.mapped_positions = 0
            self.unmapped_positions = 0
            self.positions_mapped_to_corresponding_structure = 0
            self.positions_mapped_to_homolog = 0
            self.positions_mapped_to_model = 0

    seq_id_threshold = 0.99
    startline = 'Tag\tTotal proteins\tTotal positions\tUnmapped proteins\tProteins mapped to at least one corresponding structure (seq-id > %s%%\tProteins mapped only to structure of homologs (seq-id <= %s%%)\tProteins mapped only to a model\tMapped positions\tMapped into at least one corresponding structure (seq-id > %s%%)\tMapped only in homologs (seq-id <= %s%%)\tMapped only in a model\tUnmapped' % (str(seq_id_threshold), str(seq_id_threshold), str(seq_id_threshold), str(seq_id_threshold))
    outmap = {'All': Count_helper('All')}

    if stat_dict is not None:
        class_dict = stat_dict
        max_seq_key = 1
        rec_struct_key = 2
    else:
        max_seq_key = 21
        rec_struct_key = 22
    for m in tag_map:
        if tag_map[m] is not None and tag_map[m] != '':
            raw_tags = tag_map[m].split(',')
        else:
            raw_tags = []
        tags = set(['All'])
        for tag in raw_tags:
            if tag[0] == '#':
                tag = tag[1:].split(':')[0]
            tags.add(tag)

        for tag in tags:
            if tag not in outmap:
                outmap[tag] = Count_helper(tag)
            ct = outmap[tag]
            prot_id = mutation_dict[m][1]
            if prot_id not in ct.prot_type_map:
                ct.prot_type_map[prot_id] = 1
            ct.total_positions += 1

            if m in class_dict:
                clas = class_dict[m][0]
                max_seq_id = class_dict[m][max_seq_key]
                recommended_structure = class_dict[m][rec_struct_key]
                is_model = recommended_structure[:3] == 'AF-'
                if clas is not None:
                    ct.mapped_positions += 1
                    if max_seq_id > seq_id_threshold:
                        ct.positions_mapped_to_corresponding_structure += 1
                        ct.prot_type_map[prot_id] = 1
                    elif not is_model:
                        ct.positions_mapped_to_homolog += 1
                        if not ct.prot_type_map[prot_id] == 1:
                            ct.prot_type_map[prot_id] = 2
                    else:
                        ct.positions_mapped_to_model += 1
                        if not (ct.prot_type_map[prot_id] == 1 or ct.prot_type_map[prot_id] == 2):
                            ct.prot_type_map[prot_id] = 3

                else:
                    ct.unmapped_positions += 1
                    if not ct.prot_type_map[prot_id] > 1:
                        ct.prot_type_map[prot_id] = 0
            else:
                ct.unmapped_positions += 1
                if not ct.prot_type_map[prot_id] > 1:
                    ct.prot_type_map[prot_id] = 0
    if None in outmap:
        del outmap[None]

    lines = [startline]
    for tag in outmap:
        ct = outmap[tag]
        tot_prot = len(ct.prot_type_map)
        tot_pos = ct.total_positions
        mapped = ct.mapped_positions
        unmapped = ct.unmapped_positions
        mapped_to_corr = ct.positions_mapped_to_corresponding_structure
        mapped_to_homolog = ct.positions_mapped_to_homolog
        mapped_to_model = ct.positions_mapped_to_model

        prot_numbers = [0, 0, 0, 0]
        for prot_id in ct.prot_type_map:
            prot_numbers[ct.prot_type_map[prot_id]] += 1

        if float(tot_pos) == 0.0:
            continue
        line = '%s\t%s\t%s\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)' % (
            tag,
            str(tot_prot),
            str(tot_pos),
            str(prot_numbers[0]), str(100. * float(prot_numbers[0]) / float(tot_prot)),
            str(prot_numbers[1]), str(100. * float(prot_numbers[1]) / float(tot_prot)),
            str(prot_numbers[2]), str(100. * float(prot_numbers[2]) / float(tot_prot)),
            str(prot_numbers[3]), str(100. * float(prot_numbers[3]) / float(tot_prot)),
            str(mapped), str(100. * float(mapped) / float(tot_pos)),
            str(mapped_to_corr), str(100. * float(mapped_to_corr) / float(tot_pos)),
            str(mapped_to_homolog), str(100. * float(mapped_to_homolog) / float(tot_pos)),
            str(mapped_to_model), str(100. * float(mapped_to_model) / float(tot_pos)),
            str(unmapped), str(100. * float(unmapped) / float(tot_pos)))
        lines.append(line)
    f = open(out_file, 'w')
    f.write("\n".join(lines))
    f.close()


def create_structguy_project_file(config, outfolder, session_id, structguy_project_file, feature_file, sequences_file, most_common_snv_value_tag, snv_labels, multi_savs_file):
    project_file_template_path = f'{settings.SCRIPTS_DIR}/structguy_project_template.conf'
    f = open(project_file_template_path,'r')
    lines = f.readlines()
    f.close()

    abs_path_feature_file = os.path.abspath(feature_file)
    infile = database.getSessionFile(session_id, config)
    if infile.count('/') > 0:
        dataset_name = infile.rsplit("/",1)[1].rsplit('.',1)[0]
    else:
        dataset_name = infile.rsplit('.',1)[0]

    msa_db_path = f'{config.tmp_folder}/MSA/'

    if lines[-1][-1] != '\n':
        lines[-1] += '\n'

    lines.append(f'Path_to_dataset_file = {infile}\n')
    lines.append(f'Path_to_multi_savs_table = {multi_savs_file}\n')
    lines.append(f'dataset_name = {dataset_name}\n')
    lines.append(f'Path_to_structman_features_file = {abs_path_feature_file}\n')
    lines.append(f'Path_to_sequences_file = {sequences_file}\n')
    lines.append(f'outfolder = {outfolder}\n')

    lines.append(f'mmseqs_tmp_folder = {config.tmp_folder}\n')
    lines.append(f'mmseqs_path = {config.mmseqs2_path}\n')
    lines.append(f'pdb_path = {config.pdb_path}\n')


    lines.append(f'msa_db = {msa_db_path}\n')
    if most_common_snv_value_tag is not None:
        lines.append(f'target_values = #{most_common_snv_value_tag}\n')
    if len(snv_labels) > 0:
        lines.append(f'class_labels = {",".join(snv_labels)}')

    f = open(structguy_project_file, 'w')
    f.write(''.join(lines))
    f.close()
    return



def init_classification_table(class_file, obj_only = False):
    classification_output = out_generator.OutputGenerator()

    if obj_only:
        return classification_output

    f = open(class_file, 'a')

    return classification_output, f

def generate_multi_sav_table(config, outfile, session_id, snv_position_dict, position_dict, protein_dict):
    table = 'RS_Multi_Mutation_Session'
    columns = ['Multi_Mutation', 'Tags']
    eq_cols = {'Session': session_id}

    results = select(config, columns, table, equals_rows=eq_cols)

    mm_tag_dict = {}
    for row in results:
        mm_tag_dict[row[0]] = row[1]

    table = 'Multi_Mutation'
    columns = ['Multi_Mutation_Id', 'SNVs', 'Indels', 'Wildtype_Protein', 'Mutant_Protein']

    results = binningSelect(mm_tag_dict.keys(), columns, table, config)

    multi_sav_list = []
    for row in results:
        mm_db_id = row[0]
        if mm_db_id not in mm_tag_dict:
            continue

        snv_db_ids_str = row[1]
        if snv_db_ids_str is None:
            continue
        if snv_db_ids_str.count(',') == 0:
            continue
        snv_db_ids = snv_db_ids_str.split(',')

        indel_db_ids_str = row[2]
        if indel_db_ids_str != None:
            if indel_db_ids_str.count(',') > 0:
                continue #Filter all Multi Mutations that contain non-SAV mutations (Indels)

        multi_sav_list.append((mm_db_id, snv_db_ids, row[3], row[4]))

    if config.verbosity >= 3:
        print(f'Size of Multi Sav List: {len(multi_sav_list)}')

    lines = ['Protein ID\tSAVs\tTags\n']

    for mm_db_id, snv_db_ids, wt_prot_db_id, mut_prot_db_id in multi_sav_list:
        aacs = []
        for snv_db_id in snv_db_ids:
            pos_db_id, aa2 = snv_position_dict[int(snv_db_id)]
            pos_nr, aa1, prot_db_id = position_dict[pos_db_id]
            prot_id = protein_dict[prot_db_id][0]
            aacs.append(f'{aa1}{pos_nr}{aa2}')
        line = f'{prot_id}\t{":".join(aacs)}\t{mm_tag_dict[mm_db_id]}\n'
        lines.append(line)
        
    f = open(outfile, 'w')
    f.write(''.join(lines))
    f.close()

def class_and_feature_table_main_loop(
        all_results: list[list],
        snv_map: dict,
        mutation_dict: dict,
        protein_dict: dict,
        tag_map: dict[int, str],
        pos_info_map: dict,
        position_interface_map: dict,
        interface_dict: dict,
        snv_tag_map: dict,
        config: any,
        classification_output: any,
        f: IO,
        interface_output: any,
        interface_f: IO,
        feature_output: any,
        feat_f: IO,
        max_rows_at_a_time: int = 5000000,
) -> tuple[dict, dict, dict, dict, dict]:

    main_loop_counter = 0
    stat_dict = {}
    interface_numbers = {}

    feature_header_written = False
    classification_header_written = False

    position_dict = {}

    while((main_loop_counter)*max_rows_at_a_time <= len(all_results)):
        if config.verbosity >= 2:
                t4 = time.time()
        results = all_results[(main_loop_counter*max_rows_at_a_time):((main_loop_counter+1)*max_rows_at_a_time)]

        if config.verbosity >= 3:
            print('Processing row', main_loop_counter*max_rows_at_a_time, 'till', (main_loop_counter+1)*max_rows_at_a_time)

        main_loop_counter += 1

        
        if config.verbosity >= 2:
            t5 = time.time()
            print("Time for classificationOutput part 5: ", t5 - t4, main_loop_counter)

        if config.verbosity >= 7:
            if len(snv_map) > 0:
                print(f'SNV map going into unpacked rows part:\n{snv_map}')

        warn_count = 0

        for row in results:
            m = row[0]
            position_number = row[1]
            wt_aa = row[2]
            prot_db_id = row[3]

            position_dict[m] = (position_number, wt_aa, prot_db_id)

            if config.verbosity >= 6:
                print(f'Processing in unpacked rows part: {m} {position_number} {wt_aa} {prot_db_id}')

            if row[4] is not None:
                (recommended_structure_str, max_seq_structure_str) = base_utils.unpack(row[4])
            else:
                recommended_structure_str = None
                max_seq_structure_str = None
            row = list(row)
            if row[5] is None:
                mappings = mappings_package.Mappings(init_none = True)
                microminer_features = mappings_package.Microminer_features(init_none=True)
                if config.verbosity >= 6:
                    print(f'Position Data is None for Mutation db id {m} (Protein db id {prot_db_id})')
            else:
                try:
                    packed_mappings, packed_mm_features = base_utils.unpack(row[5])
                    mappings = base_utils.unpack(packed_mappings)
                    microminer_features = base_utils.unpack(packed_mm_features)
                    if microminer_features is None:
                        microminer_features = mappings_package.Microminer_features(init_none=True)
                    if mappings is None:
                        mappings = mappings_package.Mappings(init_none=True)
                except:
                    warn_count += 1
                    if warn_count <= 10:
                        print(f'Error in retrieving position data, couldnt unpack: {m} {position_number} {prot_db_id=}')
                    continue

            rin_based_features = mappings.rin_based_features

            structural_features = mappings.structural_features

            integrated_features = mappings.integrated_features

            if config.verbosity >= 6:
                print(f'Position data unpacked rows part: {m} {position_number} {wt_aa} {prot_db_id}:\n{mappings}')

            recommended_structure, seq_id, cov, resolution = sdsc_utils.process_recommend_structure_str(recommended_structure_str)
            max_seq_structure, max_seq_seq_id, max_seq_cov, max_seq_resolution = sdsc_utils.process_recommend_structure_str(max_seq_structure_str)

            input_res_id = mutation_dict[m][2]
            if input_res_id is None:
                input_res_id = ''


            if max_seq_seq_id == '-':
                stat_dict[m] = (integrated_features.structural_classification, 0., recommended_structure)
            else:
                stat_dict[m] = (integrated_features.structural_classification, float(max_seq_seq_id), recommended_structure)

            (prot_id, u_ac, refseq, u_id, _,_, input_id, _, _) = protein_dict[prot_db_id]

            if input_id is None:
                input_id = prot_id

            classification_output.add_value('Input Protein ID', input_id)
            classification_output.add_value('Primary Protein ID', prot_id)
            classification_output.add_value('Uniprot-Ac', u_ac)
            classification_output.add_value('Uniprot ID', u_id)
            classification_output.add_value('Refseq', refseq)
            classification_output.add_value('Residue ID', input_res_id)
            classification_output.add_value('Amino Acid', wt_aa)
            classification_output.add_value('Position', position_number)
            classification_output.add_value('Tags', tag_map[m])
            classification_output.add_value('Sidechain Location', integrated_features.sidechain_location)
            classification_output.add_value('Classification', integrated_features.structural_classification)
            classification_output.add_value('Simple Class', integrated_features.simple_class)
            classification_output.add_value('Secondary Structure', structural_features.ssa)
            classification_output.add_value('Recommended Structure', recommended_structure)
            classification_output.add_value('Sequence-ID', seq_id)
            classification_output.add_value('Coverage', cov)
            classification_output.add_value('Resolution', resolution)
            classification_output.add_value('Max Seq Id Structure', max_seq_structure)
            classification_output.add_value('Max Sequence-ID', max_seq_seq_id)
            classification_output.add_value('Max Seq Id Coverage', max_seq_cov)
            classification_output.add_value('Max Seq Id Resolution', max_seq_resolution)
            classification_output.add_value('Amount of mapped structures', mappings.amount_of_structures)

            if not classification_header_written:
                f.write(classification_output.get_header())
                classification_header_written = True
            f.write(classification_output.pop_line())

            if config.compute_ppi:
                pos_info_map[m] = wt_aa, position_number

                if m in position_interface_map:
                    interface_db_id, recommended_interface_residue = position_interface_map[m]
                    if not prot_db_id in interface_numbers:
                        interface_numbers[prot_db_id] = {}
                    if not interface_db_id in interface_numbers[prot_db_id]:
                        interface_numbers[prot_db_id][interface_db_id] = len(interface_numbers[prot_db_id]) + 1
                    interface_number = interface_numbers[prot_db_id][interface_db_id]
                    interface_structure_recommendation = interface_dict[interface_db_id][1]

                    interface_output.add_value('Input Protein ID', input_id)
                    interface_output.add_value('Primary Protein ID', prot_id)
                    interface_output.add_value('Uniprot-Ac', u_ac)
                    interface_output.add_value('WT Amino Acid', wt_aa)
                    interface_output.add_value('Position', position_number)
                    interface_output.add_value('Interface Number', interface_number)
                    interface_output.add_value('Interface Structure Recommendation', interface_structure_recommendation)
                    interface_output.add_value('Position Structure Recommendation', recommended_interface_residue)

                    interface_f.write(interface_output.pop_line())

            if m not in snv_map:
                snv_map[m] = {0: wt_aa}
                snv_tag_map[0] = tag_map[m]
            else:
                if config.verbosity >= 6:
                    print(f'Processing SNVs in unpacked rows part: {m} {wt_aa} {position_number} {prot_id}')

            for snv_database_id in snv_map[m]:
                new_aa = snv_map[m][snv_database_id]

                tags = snv_tag_map[snv_database_id]
                aac = f"{wt_aa}{position_number}{new_aa}"
                feature_output.add_value('Input Protein ID', input_id)
                feature_output.add_value('Primary Protein ID', prot_id)
                feature_output.add_value('Uniprot-Ac', u_ac)
                feature_output.add_value('WT Amino Acid', wt_aa)
                feature_output.add_value('Position', position_number)
                feature_output.add_value('Mut Amino Acid', new_aa)
                feature_output.add_value('AA change', f'{wt_aa}{new_aa}')
                feature_output.add_value('Tags', tags)

                integrated_features.add_to_output_object(feature_output)
                structural_features.add_to_output_object(feature_output)
                rin_based_features.add_to_output_object(feature_output)
                microminer_features.add_to_output_object(feature_output)

                try:
                    KDmean = abs(residue_consts.HYDROPATHY[wt_aa] - residue_consts.HYDROPATHY[new_aa])
                except:
                    KDmean = 0
                    warn_count += 1
                    if warn_count < 10:
                        print(f'Calculating KDmean failed: Protein: {prot_id} Position: {position_number} WTAA: {wt_aa} MUTAA: {new_aa} SNV DB Id: {snv_database_id}')
                feature_output.add_value('KD mean', KDmean)

                try:
                    d_vol = abs(residue_consts.VOLUME[wt_aa] - residue_consts.VOLUME[new_aa])
                except:
                    d_vol = 0
                feature_output.add_value('Volume mean', d_vol)

                chemical_distance = database.getChemicalDistance(aac)
                feature_output.add_value('Chemical distance', chemical_distance)

                blosum_value = database.getBlosumValue(aac)
                feature_output.add_value('Blosum62', blosum_value)

                aliphatic_change = int((wt_aa in residue_consts.AA_MAP_ALIPHATIC) != (new_aa in residue_consts.AA_MAP_ALIPHATIC))
                hydrophobic_change = int((wt_aa in residue_consts.AA_MAP_HYDROPHOBIC) != (new_aa in residue_consts.AA_MAP_HYDROPHOBIC))
                aromatic_change = int((wt_aa in residue_consts.AA_MAP_AROMATIC) != (new_aa in residue_consts.AA_MAP_AROMATIC))
                positive_change = int((wt_aa in residue_consts.AA_MAP_POSITIVE) != (new_aa in residue_consts.AA_MAP_POSITIVE))
                polar_change = int((wt_aa in residue_consts.AA_MAP_POLAR) != (new_aa in residue_consts.AA_MAP_POLAR))
                negative_change = int((wt_aa in residue_consts.AA_MAP_NEGATIVE) != (new_aa in residue_consts.AA_MAP_NEGATIVE))
                charged_change = int((wt_aa in residue_consts.AA_MAP_CHARGED) != (new_aa in residue_consts.AA_MAP_NEGATIVE))
                small_change = int((wt_aa in residue_consts.AA_MAP_SMALL) != (new_aa in residue_consts.AA_MAP_SMALL))
                tiny_change = int((wt_aa in residue_consts.AA_MAP_TINY) != (new_aa in residue_consts.AA_MAP_TINY))
                total_change = aliphatic_change + hydrophobic_change + aromatic_change + positive_change + polar_change + negative_change + charged_change + small_change + tiny_change
                feature_output.add_value('Aliphatic change', aliphatic_change)
                feature_output.add_value('Hydrophobic change', hydrophobic_change)
                feature_output.add_value('Aromatic change', aromatic_change)
                feature_output.add_value('Positive charged change', positive_change)
                feature_output.add_value('Polar change', polar_change)
                feature_output.add_value('Negative charge change', negative_change)
                feature_output.add_value('Charged change', charged_change)
                feature_output.add_value('Small change', small_change)
                feature_output.add_value('Tiny change', tiny_change)
                feature_output.add_value('Total change', total_change)

                if not feature_header_written:
                    feat_f.write(feature_output.get_header())
                    feature_header_written = True

                feat_f.write(feature_output.pop_line())

        if warn_count >= 10:
            print(f'In total {warn_count} warnings happened')

        if config.verbosity >= 2:
            t6 = time.time()
            print("Time for classificationOutput part 5.2: ", t6 - t5)

    f.close()
    feat_f.close()
    if config.compute_ppi:
        interface_f.close()

    return pos_info_map, snv_tag_map, interface_numbers, position_dict, stat_dict

#@njit
def process_pos_pos_table_entries(
        results: list[tuple[int, int, int, int, float]],
        tag_map: set[int]
    )-> tuple[dict[int, dict[int, tuple[int, int, float]]], set[int], set[int]]:
    #)-> list[tuple[int, int, int, int, float]]:

    pos_pos_map: dict[int, dict[int, tuple[int, int, float]]] = {}
 
    #TupleType = tuple((int, int, float))
    #pos_pos_map = nb.typed.typeddict.Dict.empty(int, nb.typed.typeddict.Dict.empty(int, TupleType))
 
    #pos_pos_map = nb.typed.List() #list[tuple[int, int, int, int, float]]
    #external_ids = nb.typed.Dict.empty(key_type = int, value_type = set)

    residue_ids: set[int] = set()
    session_less_position_ids: set[int] = set()

    #print(results)

    for row in results:
        pos_a: int = row[0]

        pos_b: int = row[1]
        #if pos_a not in external_ids:
        #    external_ids[pos_a] = set()
        #elif pos_b in external_ids[pos_a]:
        #    continue

        #external_ids[pos_a].add(pos_b)
        #pos_pos_map.append(row)

        
        if pos_a not in pos_pos_map:
            pos_pos_map[pos_a] = {pos_b: (row[2], row[3], row[4])}
        else:
            if pos_b in pos_pos_map[pos_a]:
                print(f'Duplicate in pos pos table: {pos_a} {pos_b}')
            pos_pos_map[pos_a][pos_b] =  (row[2], row[3], row[4])
        
            
        residue_ids.add(row[2])
        residue_ids.add(row[3])
        if pos_a not in tag_map:
            session_less_position_ids.add(pos_a)
        if pos_b not in tag_map:
            session_less_position_ids.add(pos_b)

    return pos_pos_map, residue_ids, session_less_position_ids

def classificationOutput(
        config: any, 
        outfolder: str, 
        session_name: str, 
        session_id: int | str,
        enhance_network: bool = True) -> None:
    
    outfile: str = f'{outfolder}/{session_name}'

    if config.verbosity >= 2:
        t0 = time.time()

    table = 'RS_Position_Session'
    rows = ['Position', 'Tag']
    eq_rows = {'Session': session_id}

    results: list[list] = select(config, rows, table, equals_rows=eq_rows)

    if config.verbosity >= 2:
        t1 = time.time()
        print("Time for classificationOutput part 1: ", t1 - t0)

    tag_map: dict[int, str] = {}
    for row in results:
        mut_id: int = row[0]
        tags: str = row[1]
        tag_map[mut_id] = tags

    if config.verbosity >= 3:
        print(f'Size of tag_map: {len(tag_map)}')

    mutation_dict: dict[int, tuple[int | None, int, str | None, str]] = database.getMutationDict(list(tag_map.keys()), config)

    if config.verbosity >= 2:
        t2 = time.time()
        print("Time for classificationOutput part 2: ", t2 - t1)

    table = 'RS_SNV_Session'
    columns = ['SNV', 'Tag']

    results = select(config, columns, table, equals_rows=eq_rows)

    snv_tag_map = {}
    snv_ids = []

    most_common_snv_value_tag = None
    snv_value_tags_count = {}
    snv_labels = set()

    for row in results:
        snv_ids.append(row[0])
        snv_tags = row[1]

        for snv_tag in snv_tags.split(','):
            snv_tag_map[row[0]] = snv_tag
            if len(snv_tag) == 0:
                continue
            if snv_tag[0] == '#':
                try:
                    snv_value_tag, value_str = snv_tag[1:].split(':')
                except:
                    try:
                        snv_value_tag, value_str = snv_tag[1:].split('=')
                    except:
                        print(f'Splitting of SNV value Tag failed for SNV: {row[0]}, value tag: {snv_tag}')
                if not snv_value_tag in snv_value_tags_count:
                    snv_value_tags_count[snv_value_tag] = 1
                else:
                    snv_value_tags_count[snv_value_tag] += 1
            elif snv_tag[:6] == 'label:':
                snv_label = snv_tag[6:]
                snv_labels.add(snv_label)

    max_count = 0
    for snv_value_tag in snv_value_tags_count:
        count = snv_value_tags_count[snv_value_tag]
        if count > max_count:
            max_count = count
            most_common_snv_value_tag = snv_value_tag

    if config.verbosity >= 2:
        t21 = time.time()
        print("Time for classificationOutput part 2.1: ", t21 - t2)


    table = 'SNV'
    columns = ['SNV_Id', 'Position', 'New_AA']

    results = binningSelect(snv_ids, columns, table, config)

    snv_map = {}
    snv_position_dict = {}
    for row in results:
        if not row[1] in snv_map:
            snv_map[row[1]] = {}
        snv_map[row[1]][row[0]] = row[2]
        snv_position_dict[row[0]] = (row[1], row[2])

    if config.verbosity >= 7:
        print(f'SNV map in classificationOutput:\n{snv_map}')

    prot_id_list = set()
    input_res_id_dict = {}
    for m in mutation_dict:
        prot_id_list.add(mutation_dict[m][1])
        input_res_id_dict[m] = mutation_dict[m][2]

    if config.verbosity >= 2:
        t22 = time.time()
        print("Time for classificationOutput part 2.2: ", t22 - t21)

    protein_dict, gene_db_ids = database.getProteinDict(prot_id_list, session_id, config, includeSequence = True)

    if config.verbosity >= 2:
        t23 = time.time()
        print("Time for classificationOutput part 2.3: ", t23 - t22)

    gene_id_map = database.retrieve_gene_id_map(gene_db_ids, config)

    if config.verbosity >= 2:
        t24 = time.time()
        print("Time for classificationOutput part 2.4: ", t24 - t23)

    if config.compute_ppi:
        table = 'Interface'
        columns = ['Protein', 'Interface_Id', 'Structure_Recommendation']

        results = binningSelect(prot_id_list, columns, table, config)

        interface_dict = {}
        for row in results:
            interface_dict[row[1]] = (row[0], row[2])

        if config.verbosity >= 2:
            t251 = time.time()
            print("Time for classificationOutput part 2.5.1: ", t251 - t24)

        table = 'Protein_Protein_Interaction'
        columns = ['Interface_A', 'Interface_B', 'Complex', 'Chain_A', 'Chain_B', 'Interaction_Score']

        results = binningSelect(list(interface_dict.keys()), columns, table, config)

        prot_prot_file = out_utils.generate_ppi_filename(outfolder, session_name)

        if os.path.isfile(prot_prot_file):
            os.remove(prot_prot_file)
        prot_prot_output, prot_prot_f = init_protein_protein_table(prot_prot_file)

        ppi_map = {}
        complex_ids = set()
        interacting_interfaces = set()
        for row in results:
            if row[0] not in ppi_map:
                ppi_map[row[0]] = [row[1:]]
                interacting_interfaces.add(row[1])
            elif len(ppi_map[row[0]]) < 5:
                ppi_map[row[0]].append(row[1:])
                interacting_interfaces.add(row[1])
            elif row[5] > ppi_map[row[0]][0][4]:
                ppi_map[row[0]] = [row[1:]]
                interacting_interfaces.add(row[1])
            complex_ids.add(row[2])
            

        for int_a in ppi_map:
            interacting_interfaces.add(ppi_map[int_a][0][0])

        if config.verbosity >= 2:
            t252 = time.time()
            print("Time for classificationOutput part 2.5.2: ", t252 - t251)

        if enhance_network:
            table = 'Interface'
            columns = ['Interface_Id', 'Protein', 'Structure_Recommendation']

            results = binningSelect(interacting_interfaces, columns, table, config)

            expanded_protein_ids = set()

            for row in results:
                if row[0] not in interacting_interfaces:
                    continue
                interface_dict[row[0]] = (row[1], row[2])
                expanded_protein_ids.add(row[1])


            if config.verbosity >= 2:
                t253 = time.time()
                print("Time for classificationOutput part 2.5.3: ", t253 - t252)

            expanded_protein_dict, expanded_gene_db_ids = database.getProteinDict(expanded_protein_ids, None, config)


            if config.verbosity >= 2:
                t254 = time.time()
                print("Time for classificationOutput part 2.5.4: ", t254 - t253)

            gene_id_map.update(database.retrieve_gene_id_map(expanded_gene_db_ids, config))

            """
            if len(interacting_interfaces) < 100: #This adds edges between the enhanced interaction nodes
                table = 'Protein_Protein_Interaction'
                columns = ['Interface_A', 'Interface_B', 'Complex', 'Chain_A', 'Chain_B', 'Interaction_Score']

                results = binningSelect(interacting_interfaces, columns, table, config)
                for row in results:
                    if row[1] not in interacting_interfaces:
                        continue
                    if row[0] not in ppi_map:
                        ppi_map[row[0]] = [row[1:]]
                    else:
                        ppi_map[row[0]].append(row[1:])
                    complex_ids.add(row[2])
            """

        else:
            expanded_protein_dict = protein_dict
            if config.verbosity >= 2:
                t254 = time.time()

        table = 'Complex'
        columns = ['Complex_Id', 'PDB']

        results = binningSelect(list(complex_ids), columns, table, config)

        complex_id_map = {}
        for row in results:
            complex_id_map[row[0]] = row[1]


        if config.verbosity >= 2:
            t255 = time.time()
            print("Time for classificationOutput part 2.5.5: ", t255 - t254)

        table = 'RS_Position_Interface'
        columns = ['Interface', 'Position', 'Recommended_Residue']

        results = binningSelect(list(interface_dict.keys()), columns, table, config)

        position_interface_map = {}
        for row in results:
            position_interface_map[row[1]] = (row[0], row[2])

        if config.verbosity >= 2:
            t256 = time.time()
            print("Time for classificationOutput part 2.5.6: ", t256 - t255)

        table = 'Position_Position_Interaction'
        columns = ['Position_A', 'Position_B', 'Residue_A', 'Residue_B', 'Interaction_Score']

        results = binningSelect(list(tag_map.keys()), columns, table, config)

        pos_pos_map, residue_ids, session_less_position_ids = process_pos_pos_table_entries(results, set(tag_map.keys()))

        if config.verbosity >= 2:
            t257 = time.time()
            print(f'Time for classificationOutput part 2.5.7: {t257 - t256}, {len(pos_pos_map)}')

        table = 'Residue'
        columns = ['Residue_Id', 'Number', 'Structure']

        results = binningSelect(list(residue_ids), columns, table, config)

        res_nr_map = {}
        structure_ids = set()
        for row in results:
            res_nr_map[row[0]] = (row[1], row[2])
            structure_ids.add(row[2])

        if config.verbosity >= 2:
            t258 = time.time()
            print("Time for classificationOutput part 2.5.8: ", t258 - t257)

        table = 'Position'
        columns = ['Position_Id', 'Position_Number', 'Wildtype_Residue']

        results = binningSelect(list(session_less_position_ids), columns, table, config)

        pos_info_map = {}

        for row in results:
            pos_info_map[row[0]] = row[2], row[1]

        table = 'Structure'
        columns = ['Structure_Id', 'PDB', 'Chain']

        results = binningSelect(structure_ids, columns, table, config)

        structure_info_dict = {}
        for row in results:
            structure_info_dict[row[0]] = (row[1], row[2])

        if config.verbosity >= 2:
            t259 = time.time()
            print("Time for classificationOutput part 2.5.9: ", t259 - t258)
    else:
        interface_dict = None
        position_interface_map = None

    if config.verbosity >= 2:
        t25 = time.time()
        print("Time for classificationOutput part 2.5: ", t25 - t24)

    class_files = []

    class_file = f"{outfile}.classification.tsv"
    class_files.append(class_file)
    if os.path.isfile(class_file):
        os.remove(class_file)

    feature_file = f"{outfile}.features.tsv"
    if os.path.isfile(feature_file):
        os.remove(feature_file)

    multi_sav_table_file = f'{outfile}.multi_savs.tsv'
    if os.path.isfile(multi_sav_table_file):
        os.remove(multi_sav_table_file)

    structguy_project_file = f'{outfile}.structguy_project.conf'
    if os.path.isfile(structguy_project_file):
        os.remove(structguy_project_file)

    sequences_file = f'{outfile}.protein_sequences.fasta'
    if os.path.isfile(sequences_file):
        os.remove(sequences_file)

    fasta_lines = []
    for prot_db_id in protein_dict:
        prot_id = protein_dict[prot_db_id][0]
        seq = protein_dict[prot_db_id][-1]
        fasta_lines.append(f'>{prot_id}\n')
        fasta_lines.append(f'{seq}\n')

    f = open(sequences_file,'w')
    f.write(''.join(fasta_lines))
    f.close()

    stat_file = f"{outfile}.statistics.tsv"
    if os.path.isfile(stat_file):
        os.remove(stat_file)

    if config.compute_ppi:
        interface_file = f'{outfile}.interfaces.tsv'
        if os.path.isfile(interface_file):
            os.remove(interface_file)

        pos_pos_file = f'{outfile}.position_positions_interactions.tsv'
        if os.path.isfile(pos_pos_file):
            os.remove(pos_pos_file)
        pos_pos_output, pos_pos_f = init_position_position_table(pos_pos_file)
    else:
        interface_file = None

    if config.verbosity >= 2:
        t26 = time.time()
        print("Time for classificationOutput part 2.6: ", t26 - t25)

    if config.verbosity >= 2:
        t3 = time.time()
        print("Time for classificationOutput part 3: ", t3 - t2)

    columns = ['Position_Id', 'Position_Number', 'Wildtype_Residue', 'Protein', 'Recommended_Structure_Data', 'Position_Data']

    table = 'Position'
    all_results = binningSelect(mutation_dict.keys(), columns, table, config)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Time for classificationOutput part 4: ", t4 - t3)
        print('Total number of positions:', len(all_results))

    classification_output, f = init_classification_table(class_file)
    feature_output, feat_f = init_feature_table(feature_file)
    create_structguy_project_file(config, outfolder, session_id, structguy_project_file, feature_file, sequences_file, most_common_snv_value_tag, snv_labels, multi_sav_table_file)

    if config.compute_ppi:
        interface_output, interface_f = init_aggregated_interface_table(interface_file)
    else:
        interface_output = None
        interface_f = None

    max_rows_at_a_time = 5000000

    #use_ray = len(all_results) > config.proc_n
    use_ray = False #Currently no working with ppi output

    if use_ray:
        ray_utils.ray_init(config)
        data_store = ray.put((config, input_res_id_dict, class_file, feature_file, interface_file, interface_dict, position_interface_map))

    if config.verbosity >= 2:
        t5 = time.time()
        print("Time for classificationOutput part 4.1: ", t5 - t4)


    pos_info_map, snv_tag_map, interface_numbers, position_dict, stat_dict = class_and_feature_table_main_loop(
        all_results, snv_map, mutation_dict, protein_dict, tag_map, pos_info_map, position_interface_map, interface_dict, snv_tag_map,
        config,
        classification_output, f,
        interface_output, interface_f,
        feature_output, feat_f,
        max_rows_at_a_time = max_rows_at_a_time)

    if config.compute_ppi:
        for interface_a_db_id in ppi_map:

            prot_a_db_id, structure_rec_prot_a = interface_dict[interface_a_db_id]
            try:
                (prot_a_id, u_ac_a, refseq_a, u_id_a, error_code_a, error_a, input_id_a, gene_db_id_a, seq_a) = protein_dict[prot_a_db_id]
            except:
                (prot_a_id, u_ac_a, refseq_a, u_id_a, error_code_a, error_a, input_id_a, gene_db_id_a, seq_a) = expanded_protein_dict[prot_a_db_id]

            if input_id_a is None:
                input_id_a = prot_a_id

            if gene_db_id_a is not None:
                gene_name_a, _ = gene_id_map[gene_db_id_a]
            else:
                gene_name_a = None

            for interface_b_db_id, complex_db_id, chain_a, chain_b, ppi_interaction_score in ppi_map[interface_a_db_id]:

                try:
                    prot_b_db_id, structure_rec_prot_b = interface_dict[interface_b_db_id]
                    if not prot_b_db_id in interface_numbers:
                        interface_numbers[prot_b_db_id] = {}
                    if not interface_b_db_id in interface_numbers[prot_b_db_id]:
                        interface_numbers[prot_b_db_id][interface_b_db_id] = len(interface_numbers[prot_b_db_id]) + 1
                    interface_b_number = interface_numbers[prot_b_db_id][interface_b_db_id]
                except KeyError:
                    if config.verbosity >= 6:
                        [e, f, g] = sys.exc_info()
                        g = traceback.format_exc()
                        print(f'PPI skipped: {prot_a_id} {interface_a_db_id=} {interface_b_db_id=}\n{e}\n{f}\n{g}')
                    continue

                prot_prot_output.add_value('Gene A', gene_name_a)
                prot_prot_output.add_value('Input Protein ID A', input_id_a)
                prot_prot_output.add_value('Primary Protein ID A', prot_a_id)
                prot_prot_output.add_value('Uniprot-Ac A', u_ac_a)

                if not prot_a_db_id in interface_numbers:
                    interface_numbers[prot_a_db_id] = {}
                if not interface_a_db_id in interface_numbers[prot_a_db_id]:
                    interface_numbers[prot_a_db_id][interface_a_db_id] = len(interface_numbers[prot_a_db_id]) + 1

                interface_a_number = interface_numbers[prot_a_db_id][interface_a_db_id]
                prot_prot_output.add_value('Interface Number A', interface_a_number)

                (prot_b_id, u_ac_b, refseq_b, u_id_b, error_code_b, error_b, input_id_b, gene_db_id_b, _) = expanded_protein_dict[prot_b_db_id]

                if input_id_b is None:
                    input_id_b = prot_b_id

                if gene_db_id_b is not None:
                    gene_name_b, _ = gene_id_map[gene_db_id_b]
                else:
                    gene_name_b = None
                prot_prot_output.add_value('Gene B', gene_name_b)
                prot_prot_output.add_value('Input Protein ID B', input_id_b)
                prot_prot_output.add_value('Primary Protein ID B', prot_b_id)
                prot_prot_output.add_value('Uniprot-Ac B', u_ac_b)

                if not prot_b_db_id in interface_numbers:
                    interface_numbers[prot_b_db_id] = {}
                if not interface_b_db_id in interface_numbers[prot_b_db_id]:
                    interface_numbers[prot_b_db_id][interface_b_db_id] = len(interface_numbers[prot_b_db_id]) + 1

                interface_b_number = interface_numbers[prot_b_db_id][interface_b_db_id]
                prot_prot_output.add_value('Interface Number B', interface_b_number)

                rec_struct_id = complex_id_map[complex_db_id]

                prot_prot_output.add_value('Structure Recommendation', f'{rec_struct_id}:{chain_a}:{chain_b}')
                prot_prot_output.add_value('Interaction Score', ppi_interaction_score)

                prot_prot_f.write(prot_prot_output.pop_line())

        prot_prot_f.close()

        for pos_a_db_id in pos_pos_map:
            for pos_b_db_id in pos_pos_map[pos_a_db_id]:
                res_a_db_id, res_b_db_id, pos_pos_interaction_score = pos_pos_map[pos_a_db_id][pos_b_db_id]
                if pos_a_db_id not in position_interface_map:
                    #This shouldnt actually happen, investigate deeper ...
                    continue

                interface_db_id, recommended_interface_residue = position_interface_map[pos_a_db_id]

                if not pos_b_db_id in position_interface_map:
                    continue

                try:
                    interface_b_db_id, complex_db_id, chain_a, chain_b, ppi_interaction_score = ppi_map[interface_db_id][0]
                    rec_struct_id = complex_id_map[complex_db_id]
                except:

                    if not pos_b_db_id in position_interface_map:
                        continue
                    interface_b_db_id, recommended_interface_residue_b = position_interface_map[pos_b_db_id]
                    rec_struct_id, chain_a = recommended_interface_residue.split()[0].split(':')
                    chain_b = recommended_interface_residue_b.split()[0].split(':')[1]

                prot_a_db_id, structure_rec_prot_a = interface_dict[interface_db_id]
                try:
                    prot_b_db_id, structure_rec_prot_b = interface_dict[interface_b_db_id]
                    if not prot_db_id in interface_numbers:
                        interface_numbers[prot_b_db_id] = {}
                    if not interface_b_db_id in interface_numbers[prot_b_db_id]:
                        interface_numbers[prot_b_db_id][interface_b_db_id] = len(interface_numbers[prot_b_db_id]) + 1
                    interface_b_number = interface_numbers[prot_b_db_id][interface_b_db_id]
                except:
                    continue

                (prot_a_id, u_ac_a, refseq_a, u_id_a, error_code_a, error_a, input_id_a, gene_db_id_a, seq_a) = protein_dict[prot_a_db_id]

                try:
                    (prot_b_id, u_ac_b, refseq_b, u_id_b, error_code_b, error_b, input_id_b, gene_db_id_b, _) = expanded_protein_dict[prot_b_db_id]
                except:
                    (prot_b_id, u_ac_b, refseq_b, u_id_b, error_code_b, error_b, input_id_b, gene_db_id_b, _) = protein_dict[prot_a_db_id]

                if input_id_a is None:
                    input_id_a = prot_a_id

                pos_pos_output.add_value('Input Protein ID A', input_id_a)
                pos_pos_output.add_value('Primary Protein ID A', prot_a_id)
                pos_pos_output.add_value('Uniprot-Ac A', u_ac_a)

                wt_aa_a, position_number_a = pos_info_map[pos_a_db_id]
                pos_pos_output.add_value('WT Amino Acid A', wt_aa_a)
                pos_pos_output.add_value('Position A', position_number_a)

                if not prot_a_db_id in interface_numbers:
                    interface_numbers[prot_a_db_id] = {}
                if not interface_db_id in interface_numbers[prot_a_db_id]:
                    interface_numbers[prot_a_db_id][interface_db_id] = len(interface_numbers[prot_a_db_id]) + 1
                interface_a_number = interface_numbers[prot_a_db_id][interface_db_id]

                pos_pos_output.add_value('Interface Number A', interface_a_number)

                
                if input_id_b is None:
                    input_id_b = prot_b_id

                pos_pos_output.add_value('Input Protein ID B', input_id_b)
                pos_pos_output.add_value('Primary Protein ID B', prot_b_id)
                pos_pos_output.add_value('Uniprot-Ac B', u_ac_b)

                wt_aa_b, position_number_b = pos_info_map[pos_b_db_id]
                pos_pos_output.add_value('WT Amino Acid B', wt_aa_b)
                pos_pos_output.add_value('Position B', position_number_b)

                pos_pos_output.add_value('Interface Number B', interface_b_number)

                res_nr_a, struct_a_db_id = res_nr_map[res_a_db_id]
                res_nr_b, struct_b_db_id = res_nr_map[res_b_db_id]
                rec_pdb_id, chain_a = structure_info_dict[struct_a_db_id]
                rec_pdb_id, chain_b = structure_info_dict[struct_b_db_id]

                pos_pos_output.add_value('Structure Recommendation', f'{rec_pdb_id}:{chain_a}-{res_nr_a}:{chain_b}-{res_nr_b}')
                pos_pos_output.add_value('Interaction Score', pos_pos_interaction_score)
                pos_pos_f.write(pos_pos_output.pop_line())

        pos_pos_f.close()

    generate_multi_sav_table(config, multi_sav_table_file, session_id, snv_position_dict, position_dict, protein_dict)

    if os.path.isfile(f'{feature_file}.lock'):
        os.remove(f'{feature_file}.lock')

    if os.path.isfile(f'{class_file}.lock'):
        os.remove(f'{class_file}.lock')

    if config.compute_ppi:
        if os.path.isfile(f'{interface_file}.lock'):
            os.remove(f'{interface_file}.lock')

    if config.verbosity >= 2:
        t6 = time.time()
        print("Time for classificationOutput part 6: ", t6 - t5)

    writeStatFile(stat_file, mutation_dict, {}, tag_map, stat_dict=stat_dict)

    if config.verbosity >= 2:
        t7 = time.time()
        print("Time for classificationOutput part 7: ", t7 - t6)
        print(f'=== Total time for classificationOutput: {t7-t0} ===')

