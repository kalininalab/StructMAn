import os
import time
from typing import IO

import ray
#from numba import jit, njit
#import numba as nb

from structman import settings
from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc import mappings as mappings_package
from structman.lib.sdsc.interface import Position_Position_Interaction
from structman.lib.sdsc.residue import string_decodes
from structman.base_utils import residue_consts
from structman.lib.database import database
from structman.lib.database.database_core_functions import select, binningSelect
from structman.lib.output.feature import init_feature_table
from structman.lib.output.interfaces import init_aggregated_interface_table, init_protein_protein_table, init_position_position_table
from structman.lib.output import out_generator
from structman.base_utils import base_utils
from structman.base_utils import ray_utils
from structman.base_utils.config_class import Config
from structman.lib.output import out_utils
from structman.lib.output import alignments as alignments_package

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
    lines.append(f'msa_folder_path = {config.outpath}/msas\n')

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

def generate_multi_sav_table(config: Config, outfile, session_id, snv_position_dict, position_dict, protein_dict):
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
        if indel_db_ids_str is not None:
            if indel_db_ids_str.count(',') > 0:
                continue #Filter all Multi Mutations that contain non-SAV mutations (Indels)

        multi_sav_list.append((mm_db_id, snv_db_ids, row[3], row[4]))

    if config.verbosity >= 3:
        config.logger.info(f'Size of Multi Sav List: {len(multi_sav_list)}')

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
        snv_map: dict,
        position_list: list,
        protein_dict: dict,
        tag_map: dict[int, str],
        structure_data_dict,
        pos_info_map: dict,
        position_interface_map: dict,
        interface_dict: dict,
        snv_tag_map: dict,
        config: Config,
        classification_output: any,
        f: IO,
        interface_output: any,
        interface_f: IO,
        feature_output: any,
        feat_f: IO,
        max_rows_at_a_time: int = 5000000,
        alignment_annotations: dict | None = None,
) -> tuple[dict, dict, dict, dict, dict]:
    # alignment_annotations, when given, is filled in place with
    # {prot_db_id: {position number: alignments.position_annotation(...)}}, the
    # per-position values the --fetch_alignments matrices take from this table.

    main_loop_counter = 0
    stat_dict = {}
    interface_numbers = {}

    feature_header_written = False
    classification_header_written = False

    position_dict: dict[int, tuple[int, str, int]] = {}

    while((main_loop_counter)*max_rows_at_a_time <= len(position_list)):
        if config.verbosity >= 2:
                t4 = time.time()
  

        if config.verbosity >= 3:
            config.logger.info(f'Processing row {main_loop_counter*max_rows_at_a_time} till {(main_loop_counter+1)*max_rows_at_a_time}')
        
        if config.verbosity >= 2:
            t5 = time.time()
            config.logger.info(f"Time for classificationOutput part 5: {t5 - t4} {main_loop_counter=}")

        if config.verbosity >= 7:
            if len(snv_map) > 0:
                config.logger.info(f'SNV map going into unpacked rows part:\n{snv_map}')

        warn_count = 0

        for row in position_list[(main_loop_counter*max_rows_at_a_time):((main_loop_counter+1)*max_rows_at_a_time)]:
            (
                internal_pos_db_id,
                position_number,
                prot_db_id,
                input_res_id,
                wt_aa,
                recommended_structure,
                seq_id,
                cov,
                resolution,
                max_seq_structure,
                max_seq_seq_id,
                max_seq_cov,
                max_seq_resolution,
                packed_pos_data
            ) = row

            position_dict[internal_pos_db_id] = (position_number, wt_aa, prot_db_id)

            if config.verbosity >= 6:
                config.logger.info(f'Processing in unpacked rows part: {internal_pos_db_id} {position_number} {wt_aa} {prot_db_id}')

            
            if packed_pos_data is None:
                mappings = mappings_package.Mappings(init_none = True)
                microminer_features = mappings_package.Microminer_features(init_none=True)
                if config.verbosity >= 6:
                    config.logger.info(f'Position Data is None for Mutation db id {internal_pos_db_id} (Protein db id {prot_db_id})')
            else:
                try:
                    packed_mappings, packed_mm_features = base_utils.unpack(packed_pos_data)
                    mappings = base_utils.unpack(packed_mappings)
                    microminer_features = base_utils.unpack(packed_mm_features)
                    if microminer_features is None:
                        microminer_features = mappings_package.Microminer_features(init_none=True)
                        if config.verbosity >= 6:
                            config.logger.info(f'Microminer features are None for {internal_pos_db_id=} {prot_db_id=}')
                    if mappings is None:
                        mappings = mappings_package.Mappings(init_none=True)
                except:
                    warn_count += 1
                    if warn_count <= 10:
                        config.logger.info(f'Error in retrieving position data, couldnt unpack: {internal_pos_db_id} {position_number} {prot_db_id=}')
                    continue

            rin_based_features = mappings.rin_based_features

            structural_features = mappings.structural_features

            integrated_features = mappings.integrated_features

            if config.verbosity >= 6:
                config.logger.info(f'Position data unpacked rows part: {internal_pos_db_id=} {position_number=} {wt_aa=} {prot_db_id=} {recommended_structure=}:\n{mappings}')

            if input_res_id is None:
                input_res_id = ''

            if max_seq_seq_id == '-':
                stat_dict[internal_pos_db_id] = (integrated_features.structural_classification, 0., recommended_structure)
            else:
                stat_dict[internal_pos_db_id] = (integrated_features.structural_classification, float(max_seq_seq_id), recommended_structure)

            (prot_id, u_ac, refseq, u_id, _,_, input_id, prot_specific_tags, _, _) = protein_dict[prot_db_id]

            function_type = out_utils.get_function_type_from_tags(prot_specific_tags)

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
            classification_output.add_value('Tags', tag_map[internal_pos_db_id])
            classification_output.add_value('Sidechain Location', integrated_features.sidechain_location)
            classification_output.add_value('Classification', integrated_features.structural_classification)
            classification_output.add_value('Simple Class', integrated_features.simple_class)
            try:
                decoded_ssa_string = string_decodes[3][int(structural_features.ssa)]
            except:
                decoded_ssa_string = '-'
            classification_output.add_value('Secondary Structure', decoded_ssa_string)

            # Before recommended_structure is rewritten below: the matrices need its
            # residue contig index, not the residue number the table shows.
            if alignment_annotations is not None:
                if prot_db_id not in alignment_annotations:
                    alignment_annotations[prot_db_id] = {}
                alignment_annotations[prot_db_id][position_number] = alignments_package.position_annotation(
                    integrated_features.simple_class, recommended_structure)

            """
            try:
                recommended_structure_words = recommended_structure.split(':')
                if len(recommended_structure_words) == 3:
                    recommended_structure = ':'.join(recommended_structure_words) + f':{wt_aa}'
                elif len(recommended_structure_words) == 4:
                    struct_aa = string_decodes[31][int(recommended_structure_words[3])]
                    recommended_structure = ':'.join(recommended_structure_words[:3]) + f':{struct_aa}'
            except:
                ...
            """
            if recommended_structure != '-' and recommended_structure != '':
                try:
                    struct_id_chain_id, res_contig_aa = recommended_structure.split()
                except ValueError as err:
                    raise ValueError(f'{err=} {recommended_structure=}')

                structure_id, chain = struct_id_chain_id.split(':')
                res_contig, aa = res_contig_aa.split(':')

                if structure_id != '':
                    try:
                        res_nr = structure_data_dict[structure_id][chain][int(res_contig)-1]
                    except IndexError as err:
                        raise IndexError(f'{err=} {structure_id=} {chain=} {res_contig=}')
                    except KeyError as err:
                        raise KeyError(f'{err=} {internal_pos_db_id=} {recommended_structure=}')
                else:
                    config.errorlog.add_warning(f'{structure_id=} {internal_pos_db_id=} {recommended_structure=}')
                    res_nr = 'NA'

                resolution = out_utils.hide_emulated_resolution(structure_id, resolution)

                recommended_structure = f'{struct_id_chain_id} {res_nr}:{aa}'

            classification_output.add_value('Recommended Structure', recommended_structure)
            classification_output.add_value('Sequence-ID', seq_id)
            classification_output.add_value('Coverage', cov)
            classification_output.add_value('Resolution', resolution)

            """
            try:
                max_seq_structure_words = max_seq_structure.split(':')
                if len(max_seq_structure_words) == 3:
                    max_seq_structure = ':'.join(max_seq_structure_words) + f':{wt_aa}'
                elif len(max_seq_structure_words) == 4:
                    struct_aa = string_decodes[31][int(max_seq_structure_words[3])]
                    max_seq_structure = ':'.join(max_seq_structure_words[:3]) + f':{struct_aa}'
            except:
                ...
            """
            
            if max_seq_structure != '-' and max_seq_structure != '':
                struct_id_chain_id, res_contig_aa = max_seq_structure.split()

                structure_id, chain = struct_id_chain_id.split(':')
                res_contig, aa = res_contig_aa.split(':')

                if structure_id != '':
                    try:
                        res_nr = structure_data_dict[structure_id][chain][int(res_contig)-1]
                    except KeyError as err:
                        raise KeyError(f'{err=} {internal_pos_db_id=} {recommended_structure=} {max_seq_structure=}')
                else:
                    res_nr = 'NA'

                max_seq_resolution = out_utils.hide_emulated_resolution(structure_id, max_seq_resolution)

                max_seq_structure = f'{struct_id_chain_id} {res_nr}:{aa}'

            classification_output.add_value('Max Seq Id Structure', max_seq_structure)
            classification_output.add_value('Max Sequence-ID', max_seq_seq_id)
            classification_output.add_value('Max Seq Id Coverage', max_seq_cov)
            classification_output.add_value('Max Seq Id Resolution', max_seq_resolution)
            classification_output.add_value('Amount of mapped structures', mappings.amount_of_structures)

            if not classification_header_written:
                f.write(classification_output.get_header())
                classification_header_written = True
            f.write(classification_output.pop_line())

            if config.ppi_data_available:
                pos_info_map[internal_pos_db_id] = wt_aa, position_number

                if internal_pos_db_id in position_interface_map:
                    interface_db_id, recommended_interface_residue = position_interface_map[internal_pos_db_id]
                    if prot_db_id not in interface_numbers:
                        interface_numbers[prot_db_id] = {}
                    if interface_db_id not in interface_numbers[prot_db_id]:
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

            if internal_pos_db_id not in snv_map:
                snv_map[internal_pos_db_id] = {0: wt_aa}
                snv_tag_map[0] = tag_map[internal_pos_db_id]
            else:
                if config.verbosity >= 6:
                    config.logger.info(f'Processing SNVs in unpacked rows part: {internal_pos_db_id} {wt_aa} {position_number} {prot_id}')

            for snv_database_id in snv_map[internal_pos_db_id]:
                new_aa = snv_map[internal_pos_db_id][snv_database_id]

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

                feature_output.add_value('Function Type', function_type)

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
                        config.logger.info(f'Calculating KDmean failed: Protein: {prot_id} Position: {position_number} WTAA: {wt_aa} MUTAA: {new_aa} SNV DB Id: {snv_database_id}')
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

        
        main_loop_counter += 1

        if warn_count >= 10:
            config.logger.info(f'In total {warn_count} warnings happened')

        if config.verbosity >= 2:
            t6 = time.time()
            config.logger.info(f"Time for classificationOutput part 5.2: {t6 - t5}")

    f.close()
    feat_f.close()
    if config.ppi_data_available:
        interface_f.close()

    return pos_info_map, snv_tag_map, position_dict, stat_dict

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
            #if pos_b in pos_pos_map[pos_a]:
            #    print(f'Duplicate in pos pos table: {pos_a} {pos_b}')
            pos_pos_map[pos_a][pos_b] =  (row[2], row[3], row[4])
        
            
        residue_ids.add(row[2])
        residue_ids.add(row[3])
        if pos_a not in tag_map:
            session_less_position_ids.add(pos_a)
        if pos_b not in tag_map:
            session_less_position_ids.add(pos_b)

    return pos_pos_map, residue_ids, session_less_position_ids

def transpose_pos_pos_data(
        pos_pos_dict: dict[int, dict[int, Position_Position_Interaction]]
        ) -> dict[int, dict[int, Position_Position_Interaction]]:
    """Pos_Pos_Data seen from the other protein of the pair.

    The blob is nested {position of the stored Protein_A: {position of the stored
    Protein_B: interaction}}, and everything downstream keys it by the session
    protein's own positions - a position belonging to the partner is dropped by the
    `pos_a_db_id not in position_dict` check in the output loop. A row picked up by
    the Protein_B pass is therefore turned around here.

    recommended_interaction has to turn with it: it is (chain_a, res_a, chain_b,
    res_b) in the orientation the row was stored with, and the output attributes its
    first half to protein A and its second half to protein B. Nothing else in the
    object is orientation dependent - protein_a/position_a/protein_b/position_b are
    outside Position_Position_Interaction's slot mask and never survive the pack, so
    the nesting is the only record of which side is which.
    """
    transposed: dict[int, dict[int, Position_Position_Interaction]] = {}

    for pos_a_db_id, partners in pos_pos_dict.items():
        for pos_b_db_id, interaction in partners.items():
            try:
                chain_a, res_a, chain_b, res_b = interaction.recommended_interaction
                interaction.recommended_interaction = (chain_b, res_b, chain_a, res_a)
            except (AttributeError, TypeError, ValueError):
                # An interaction that never got a recommendation carries None here, and
                # the output skips it anyway; leave it rather than failing the row.
                pass

            if pos_b_db_id not in transposed:
                transposed[pos_b_db_id] = {}
            transposed[pos_b_db_id][pos_a_db_id] = interaction

    return transposed


def collect_protein_protein_interactions(
        prot_id_set: set[int],
        tag_map: dict,
        output_pos_pos: bool,
        config: Config
        ) -> tuple[dict, set[int], set[int], dict, set[int]]:
    """Every stored interaction that involves a protein of this session.

    Returns (ppi_map, complex_ids, expanded_protein_ids, pos_pos_map,
    session_less_position_ids). ppi_map is keyed by the session protein and holds
    up to five (partner, complex, chain_a, chain_b, score) entries seen from that
    protein's side; expanded_protein_ids collects partners that are not themselves
    in the session, so their names can be looked up separately.
    """
    ppi_map: dict[int, list[tuple[int, int, str, str, float]]] = {}
    complex_ids: set[int] = set()
    expanded_protein_ids: set[int] = set()
    pos_pos_map: dict[int, dict[int, Position_Position_Interaction]] = {}
    session_less_position_ids: set[int] = set()

    # A Protein_Protein_Interaction row is written from one side only: Protein_A is
    # whichever protein had its interfaces computed in the run that stored the row
    # (insertion_lib.remote_insert_interfaces), and no reciprocal row is written.
    # binningSelect keys on columns[0], so selecting on Protein_A alone found an
    # interaction only from the side that happened to be annotated - a session
    # holding just the partner saw nothing, though the row was right there in the
    # database. The second pass keys on Protein_B (indexed on every one of the 20
    # tables) and returns the same tuple shape with the two sides already swapped.
    #
    # Pos_Pos_Data is a large packed blob that is only ever consumed inside the
    # output_pos_pos branch below. Selecting and unpacking it regardless costs
    # both the DB transfer and the deserialization for nothing.
    forward_columns = ['Protein_A', 'Protein_B', 'Complex', 'Chain_A', 'Chain_B', 'Interaction_Score']
    reverse_columns = ['Protein_B', 'Protein_A', 'Complex', 'Chain_B', 'Chain_A', 'Interaction_Score']
    if output_pos_pos:
        forward_columns = forward_columns + ['Pos_Pos_Data']
        reverse_columns = reverse_columns + ['Pos_Pos_Data']

    # Both rows exist whenever both proteins have been annotated at some point, so
    # the same interaction reaches the same perspective twice - once per pass.
    seen_ppi: set[tuple[int, int, int, str, str]] = set()

    for i in range(20):
        table = f'Protein_Protein_Interaction_{i}'

        for columns, swapped in ((forward_columns, False), (reverse_columns, True)):

            results = binningSelect(prot_id_set, columns, table, config)

            for row in results:
                prot_a_db_id, prot_b_db_id, complex_db_id, chain_a, chain_b, ppi_score = row[:6]
                transpose = swapped

                # A self interaction is its own mirror image: it matches both passes,
                # and a homodimer also stores one row per interface, so the same
                # contact arrives as (chain_x, chain_y) and as (chain_y, chain_x).
                # Ordering the chain pair when both sides are the same protein folds
                # those onto a single entry, and picking the order rather than taking
                # whichever the database happened to return first keeps the output
                # stable between runs. The pos_pos nesting follows the chains.
                if (prot_a_db_id == prot_b_db_id and chain_a is not None
                        and chain_b is not None and chain_b < chain_a):
                    chain_a, chain_b = chain_b, chain_a
                    transpose = not transpose

                ppi_key = (prot_a_db_id, prot_b_db_id, complex_db_id, chain_a, chain_b)
                if ppi_key in seen_ppi:
                    continue
                seen_ppi.add(ppi_key)

                ppi_entry = (prot_b_db_id, complex_db_id, chain_a, chain_b, ppi_score)

                if prot_a_db_id not in ppi_map:
                    ppi_map[prot_a_db_id] = [ppi_entry]
                    #interacting_interfaces.add(prot_b_db_id)
                elif len(ppi_map[prot_a_db_id]) < 5:
                    ppi_map[prot_a_db_id].append(ppi_entry)
                    #interacting_interfaces.add(prot_b_db_id)
                elif ppi_score > ppi_map[prot_a_db_id][0][4]:
                    ppi_map[prot_a_db_id] = [ppi_entry]
                    #interacting_interfaces.add(prot_b_db_id)
                complex_ids.add(complex_db_id)
                if prot_b_db_id not in prot_id_set:
                    expanded_protein_ids.add(prot_b_db_id)

                #config.logger.info(f'{row[:5]=}')
                if output_pos_pos:
                    pos_pos_dict = base_utils.unpack(row[6])

                    if transpose:
                        pos_pos_dict = transpose_pos_pos_data(pos_pos_dict)

                    for pos_a_db_id in pos_pos_dict:
                        if pos_a_db_id == 0:
                            continue
                        if pos_a_db_id not in tag_map:
                            session_less_position_ids.add(pos_a_db_id)
                        if pos_a_db_id not in pos_pos_map:
                            pos_pos_map[pos_a_db_id] = {}
                        for pos_b_db_id in pos_pos_dict[pos_a_db_id]:
                            if pos_b_db_id == 0:
                                continue
                            if pos_b_db_id not in tag_map:
                                session_less_position_ids.add(pos_b_db_id)
                            pos_pos_map[pos_a_db_id][pos_b_db_id] = pos_pos_dict[pos_a_db_id][pos_b_db_id]
                            complex_ids.add(pos_pos_map[pos_a_db_id][pos_b_db_id].recommended_complex)

    return ppi_map, complex_ids, expanded_protein_ids, pos_pos_map, session_less_position_ids


def classificationOutput(
        config: any, 
        outfolder: str, 
        session_name: str, 
        session_id: int | str,
        enhance_network: bool = True) -> dict | None:
    
    outfile: str = f'{outfolder}/{session_name}'

    if config.verbosity >= 2:
        t0 = time.time()

    table = 'RS_Position_Session'
    rows = ['Position', 'Tag']
    eq_rows = {'Session': session_id}

    results: list[list] = select(config, rows, table, equals_rows=eq_rows)

    if config.verbosity >= 2:
        t1 = time.time()
        config.logger.info(f"Time for classificationOutput part 1: {t1 - t0}")

    tag_map: dict[int, str] = {}
    for row in results:
        mut_id: int = row[0]
        tags: str = row[1]
        tag_map[mut_id] = tags

    if config.verbosity >= 3:
        config.logger.info(f'Size of tag_map: {len(tag_map)}')

    position_list, struct_ids = database.getPositions(list(tag_map.keys()), config)

    if config.verbosity >= 2:
        t2 = time.time()
        config.logger.info(f"Time for classificationOutput part 2: {t2 - t1}")

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
                        config.logger.info(f'Splitting of SNV value Tag failed for SNV: {row[0]}, value tag: {snv_tag}')
                if snv_value_tag not in snv_value_tags_count:
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
        config.logger.info(f"Time for classificationOutput part 2.1: {t21 - t2}")


    table = 'SNV'
    columns = ['SNV_Id', 'Position', 'New_AA']

    results = binningSelect(snv_ids, columns, table, config)

    snv_map = {}
    snv_position_dict = {}
    for row in results:
        if row[1] not in snv_map:
            snv_map[row[1]] = {}
        snv_map[row[1]][row[0]] = row[2]
        snv_position_dict[row[0]] = (row[1], row[2])

    if config.verbosity >= 7:
        config.logger.info(f'SNV map in classificationOutput:\n{snv_map}')

    prot_id_set: set[int] = set([x[2] for x in position_list])

    if len(prot_id_set) > 50:
        output_pos_pos = False
    else:
        output_pos_pos = True

    if config.verbosity >= 2:
        t22 = time.time()
        config.logger.info(f"Time for classificationOutput part 2.2: {t22 - t21} {len(prot_id_set)=} {output_pos_pos=}")

    protein_dict, gene_db_ids = database.getProteinDict(prot_id_set, session_id, config, includeSequence = True)

    if config.verbosity >= 2:
        t23 = time.time()
        config.logger.info(f"Time for classificationOutput part 2.3: {t23 - t22}")

    gene_id_map = database.retrieve_gene_id_map(gene_db_ids, config)

    if config.verbosity >= 2:
        t24 = time.time()
        config.logger.info(f"Time for classificationOutput part 2.4: {t24 - t23}")

    if config.ppi_data_available:
        table = 'Interface'
        columns = ['Protein', 'Interface_Id', 'Structure_Recommendation']

        results = binningSelect(prot_id_set, columns, table, config)

        interface_dict = {}
        for row in results:
            interface_dict[row[1]] = (row[0], row[2])


        #table = 'Position_Position_Interaction'
        #columns = ['Position_A', 'Position_B', 'Residue_A', 'Residue_B', 'Interaction_Score']

        #results = binningSelect(list(tag_map.keys()), columns, table, config)

        #pos_pos_map, residue_ids, session_less_position_ids = process_pos_pos_table_entries(results, set(tag_map.keys()))

        table = 'RS_Position_Interface'
        columns = ['Position', 'Interface', 'Recommended_Residue']

        results = binningSelect(list(tag_map.keys()), columns, table, config)

        interacting_interfaces = set()
        position_interface_map = {}
        for row in results:
            position_interface_map[row[0]] = (row[1], row[2])
            interacting_interfaces.add(row[1])

        if config.verbosity >= 2:
            t251 = time.time()
            config.logger.info(f"Time for classificationOutput part 2.5.1: {t251 - t24}")

        prot_prot_file = out_utils.generate_ppi_filename(outfolder, session_name)
        network_file = f'{outfolder}/{session_name}_ppi.sif'
        sif_lines = []

        if os.path.isfile(prot_prot_file):
            os.remove(prot_prot_file)
        prot_prot_output, prot_prot_f = init_protein_protein_table(prot_prot_file)

        (ppi_map, complex_ids, expanded_protein_ids,
         pos_pos_map, session_less_position_ids) = collect_protein_protein_interactions(
            prot_id_set, tag_map, output_pos_pos, config)

        if config.verbosity >= 2:
            t252 = time.time()
            config.logger.info(f"Time for classificationOutput part 2.5.2: {t252 - t251} {len(complex_ids)=} {len(ppi_map)=} {len(pos_pos_map)=}")

        if enhance_network:
            table = 'Interface'
            columns = ['Interface_Id', 'Protein', 'Structure_Recommendation']

            results = binningSelect(interacting_interfaces, columns, table, config)

            for row in results:
                if row[0] not in interacting_interfaces:
                    continue
                interface_dict[row[0]] = (row[1], row[2])

            if config.verbosity >= 2:
                t253 = time.time()
                config.logger.info(f"Time for classificationOutput part 2.5.3: {t253 - t252}")

            expanded_protein_dict, expanded_gene_db_ids = database.getProteinDict(expanded_protein_ids, None, config)


            if config.verbosity >= 2:
                t254 = time.time()
                config.logger.info(f"Time for classificationOutput part 2.5.4: {t254 - t253}")

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
            struct_ids.add(row[1])


        if config.verbosity >= 2:
            t255 = time.time()
            config.logger.info(f"Time for classificationOutput part 2.5.5: {t255 - t254} {len(complex_id_map)=}")


        if config.verbosity >= 2:
            t256 = t257 = time.time()
            config.logger.info(f'Time for classificationOutput part 2.5.6: {t256 - t255}, {len(pos_pos_map)=}')

        #table = 'Residue'
        #columns = ['Residue_Id', 'Number', 'Structure']

        #results = binningSelect(list(residue_ids), columns, table, config)

        #res_nr_map = {}
        #structure_ids = set()
        #for row in results:
        #    res_nr_map[row[0]] = (row[1], row[2])
        #    structure_ids.add(row[2])

        if config.verbosity >= 2:
            t258 = time.time()
            config.logger.info(f"Time for classificationOutput part 2.5.8: {t258 - t257}")

        
        # Protein is selected alongside the position itself: the partner side of a
        # pos-pos interaction is almost never a position of this session, so it is
        # missing from position_dict (which only ever sees position_list), and the
        # output loop needs its protein to name the B side. Without the column,
        # every partner position resolved to nothing and the whole pos-pos table
        # came out empty. Same idea as expanded_protein_dict one level down.
        columns = ['Position_Id', 'Position_Number', 'Wildtype_Residue', 'Protein']

        pos_db_bins = [[] for _ in range(20)]
        for internal_pos_db in session_less_position_ids:
            db_nr = internal_pos_db % 100
            pos_db_id = internal_pos_db // 100
            pos_db_bins[db_nr].append(pos_db_id)

        pos_info_map = {}
        expanded_position_dict: dict[int, tuple[int, str, int]] = {}

        for db_nr, pos_ids in enumerate(pos_db_bins):
            table = f'Position_{db_nr}'
            results = binningSelect(pos_ids, columns, table, config)

            for row in results:
                internal_pos_db_id = 100*row[0] + db_nr
                pos_info_map[internal_pos_db_id] = row[2], row[1]
                # Same (number, wt_aa, protein) shape as position_dict, so the
                # output loop can fall through from one to the other.
                expanded_position_dict[internal_pos_db_id] = (row[1], row[2], row[3])

        #table = 'Structure'
        #columns = ['Structure_Id', 'PDB', 'Chain']

        #results = binningSelect(structure_ids, columns, table, config)

        #structure_info_dict = {}
        #for row in results:
        #    structure_info_dict[row[0]] = (row[1], row[2])

        if config.verbosity >= 2:
            t259 = time.time()
            config.logger.info(f"Time for classificationOutput part 2.5.9: {t259 - t258}")
    else:
        interface_dict = None
        position_interface_map = None
        expanded_position_dict = {}
        # Only ever filled by the PPI branch (for positions of *other* sessions that
        # a pos-pos partner points at); the main loop fills its own entries.
        pos_info_map = {}

    # Hoisted out of the PPI branch: the recommended-structure and max-seq-structure
    # columns of the classification table resolve their residue numbers through this
    # map, so it is needed whether or not PPI data exists. Built after the branch
    # because the PPI branch extends struct_ids with the complexes it selected.
    structure_data_list = database.get_parse_data_from_struct_ids(list(struct_ids), config)

    structure_data_dict = {}
    for struct_db_id, struct_id, chain, packed_parse_data in structure_data_list:
        if struct_id not in structure_data_dict:
            structure_data_dict[struct_id] = {}
        (seq_res_map, template_seq) = base_utils.unpack(packed_parse_data)
        structure_data_dict[struct_id][chain] = seq_res_map

    if config.verbosity >= 2:
        t25 = time.time()
        config.logger.info(f"Time for classificationOutput part 2.5: {t25 - t24}")

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

    if config.ppi_data_available:
        interface_file = f'{outfile}.interfaces.tsv'
        if os.path.isfile(interface_file):
            os.remove(interface_file)

        if output_pos_pos:
            pos_pos_file = f'{outfile}.position_positions_interactions.tsv'
            if os.path.isfile(pos_pos_file):
                os.remove(pos_pos_file)
            pos_pos_output, pos_pos_f = init_position_position_table(pos_pos_file)
    else:
        interface_file = None

    if config.verbosity >= 2:
        t26 = time.time()
        config.logger.info(f"Time for classificationOutput part 2.6: {t26 - t25}")

    if config.verbosity >= 2:
        t3 = time.time()
        config.logger.info(f"Time for classificationOutput part 3: {t3 - t2}")

    
    if config.verbosity >= 2:
        t4 = time.time()
        config.logger.info(f"Time for classificationOutput part 4: {t4 - t3}")

    classification_output, f = init_classification_table(class_file)
    feature_output, feat_f = init_feature_table(feature_file)
    create_structguy_project_file(config, outfolder, session_id, structguy_project_file, feature_file, sequences_file, most_common_snv_value_tag, snv_labels, multi_sav_table_file)

    if config.ppi_data_available:
        interface_output, interface_f = init_aggregated_interface_table(interface_file)
    else:
        interface_output = None
        interface_f = None

    max_rows_at_a_time = 5000000

    if config.verbosity >= 2:
        t5 = time.time()
        config.logger.info(f"Time for classificationOutput part 4.1: {t5 - t4}")


    # Only kept when the alignment matrices are written, which read it back.
    alignment_annotations = {} if config.fetch_alignments else None

    position_dict: dict[int, tuple[int, str, int]]
    pos_info_map, snv_tag_map, position_dict, stat_dict = class_and_feature_table_main_loop(
        snv_map, position_list, protein_dict, tag_map, structure_data_dict,
        pos_info_map, position_interface_map, interface_dict, snv_tag_map,
        config,
        classification_output, f,
        interface_output, interface_f,
        feature_output, feat_f,
        max_rows_at_a_time = max_rows_at_a_time,
        alignment_annotations = alignment_annotations)

    if config.ppi_data_available:
        for prot_a_db_id in ppi_map:

            try:
                (prot_a_id, u_ac_a, _, _, _, _, input_id_a, _, gene_db_id_a, _) = protein_dict[prot_a_db_id]
            except KeyError:
                (prot_a_id, u_ac_a, _, _, _, _, input_id_a, _, gene_db_id_a, _) = expanded_protein_dict[prot_a_db_id]

            if input_id_a is None:
                input_id_a = prot_a_id

            if gene_db_id_a is not None:
                gene_name_a, _ = gene_id_map[gene_db_id_a]
            else:
                gene_name_a = None

            for prot_b_db_id, complex_db_id, chain_a, chain_b, ppi_interaction_score in ppi_map[prot_a_db_id]:

                prot_prot_output.add_value('Gene A', gene_name_a)
                prot_prot_output.add_value('Input Protein ID A', input_id_a)
                prot_prot_output.add_value('Primary Protein ID A', prot_a_id)
                prot_prot_output.add_value('Uniprot-Ac A', u_ac_a)

                try:
                    (prot_b_id, u_ac_b, _, _, _, _, input_id_b, _, gene_db_id_b, _) = protein_dict[prot_b_db_id]
                except KeyError:
                    (prot_b_id, u_ac_b, _, _, _, _, input_id_b, _, gene_db_id_b, _) = expanded_protein_dict[prot_b_db_id]

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

                rec_struct_id = complex_id_map[complex_db_id]

                prot_prot_output.add_value('Structure Recommendation', f'{rec_struct_id}:{chain_a}:{chain_b}')
                prot_prot_output.add_value('Interaction Score', ppi_interaction_score)

                prot_prot_f.write(prot_prot_output.pop_line())

                sif_lines.append(f'{input_id_a} pp {input_id_b}\n')

        f = open(network_file, 'w')
        f.write(''.join(sif_lines))
        f.close()

        prot_prot_f.close()

        if output_pos_pos:
            count_1 = 0
            count_2 = 0
            for pos_a_db_id in pos_pos_map:
                if pos_a_db_id not in position_dict:
                    count_1 += 1
                    continue
                prot_a_db_id = position_dict[pos_a_db_id][2]

                if pos_a_db_id not in position_interface_map:
                    config.logger.info('This shouldnt actually happen, investigate deeper ...')
                    continue

                (prot_a_id, u_ac_a, _, _, _, _, input_id_a, _, gene_db_id_a, _) = protein_dict[prot_a_db_id]

                inter_count = 0
                for pos_b_db_id in pos_pos_map[pos_a_db_id]:
                    if pos_b_db_id in position_dict:
                        prot_b_db_id = position_dict[pos_b_db_id][2]
                    else:
                        # The partner position belongs to a protein that is not in
                        # this session, so it was fetched separately above rather
                        # than coming through position_list.
                        prot_b_db_id = expanded_position_dict.get(pos_b_db_id, (None, None, None))[2]

                    if prot_b_db_id is None:
                        inter_count += 1
                        continue

                    _, complex_db_id, chain_a, chain_b, ppi_interaction_score = ppi_map[prot_a_db_id][0]

                    pos_pos_i_obj: Position_Position_Interaction = pos_pos_map[pos_a_db_id][pos_b_db_id]
                    
                    rec_struct_id = complex_id_map[complex_db_id]
                    
                    try:
                        (prot_b_id, u_ac_b, _, _, _, _, input_id_b, _, gene_db_id_b, _) = expanded_protein_dict[prot_b_db_id]
                    except:
                        (prot_b_id, u_ac_b, _, _, _, _, input_id_b, _, gene_db_id_b, _) = protein_dict[prot_b_db_id]

                    if input_id_a is None:
                        input_id_a = prot_a_id

                    pos_pos_output.add_value('Input Protein ID A', input_id_a)
                    pos_pos_output.add_value('Primary Protein ID A', prot_a_id)
                    pos_pos_output.add_value('Uniprot-Ac A', u_ac_a)

                    wt_aa_a, position_number_a = pos_info_map[pos_a_db_id]
                    pos_pos_output.add_value('WT Amino Acid A', wt_aa_a)
                    pos_pos_output.add_value('Position A', position_number_a)

                    if input_id_b is None:
                        input_id_b = prot_b_id

                    pos_pos_output.add_value('Input Protein ID B', input_id_b)
                    pos_pos_output.add_value('Primary Protein ID B', prot_b_id)
                    pos_pos_output.add_value('Uniprot-Ac B', u_ac_b)

                    wt_aa_b, position_number_b = pos_info_map[pos_b_db_id]

                    

                    pos_pos_output.add_value('WT Amino Acid B', wt_aa_b)
                    pos_pos_output.add_value('Position B', position_number_b)

                    
                    rec_compl_db_id = pos_pos_i_obj.recommended_complex
                    rec_pdb_id = complex_id_map[rec_compl_db_id]
                    (chain_a, res_contig_a, chain_b, res_contig_b) = pos_pos_i_obj.recommended_interaction
                    pos_pos_interaction_score = pos_pos_i_obj.interaction_score

                    #config.logger.info(f'{prot_a_db_id=} {prot_a_id=} {pos_a_db_id=} {wt_aa_a=} {prot_b_db_id=} {prot_b_id=} {pos_b_db_id=} {wt_aa_b=} {position_number_b=} {pos_pos_i_obj.recommended_interaction=}')

                    try:
                        # -1 because res_contig is a 1-based contig index (contig_help_map
                        # starts at 1 in templateFiltering) while seq_res_map is a 0-based
                        # list of residue numbers in sequence order. This lookup was the
                        # only one of the four in this file missing it, so chain A's
                        # residue came out one position too far along - and the last
                        # residue of a chain raised IndexError and was silently dropped by
                        # the `continue` below.
                        res_nr_a = structure_data_dict[rec_pdb_id][chain_a][int(res_contig_a)-1]
                    except IndexError as err:
                        if config.verbosity >=3:
                            config.logger.info(f'{err=} {rec_pdb_id=} {chain_a=} {res_contig_a=} {structure_data_dict[rec_pdb_id][chain_a]=} {prot_a_id=} {pos_a_db_id=} {wt_aa_a=} {position_number_a=}')
                        #raise IndexError(f'{err=} {rec_pdb_id=} {chain_a=} {res_contig_a=} {structure_data_dict[rec_pdb_id][chain_a]=} {prot_a_id=} {pos_a_db_id=} {wt_aa_a=} {position_number_a=}')
                        continue
                    try:
                        res_nr_b = structure_data_dict[rec_pdb_id][chain_b][int(res_contig_b)-1]
                    except IndexError as err:
                        if config.verbosity >=3:
                            config.logger.info(f'{err=} {rec_pdb_id=} {chain_b=} {res_contig_b=} {structure_data_dict[rec_pdb_id][chain_b]=} {prot_b_id=} {pos_b_db_id=} {wt_aa_b=} {position_number_b=}')
                        #raise IndexError(f'{err=} {rec_pdb_id=} {chain_b=} {res_contig_b=} {structure_data_dict[rec_pdb_id][chain_b]=} {prot_b_id=} {pos_b_db_id=} {wt_aa_b=} {position_number_b=}')
                        continue

                    pos_pos_output.add_value('Structure Recommendation', f'{rec_pdb_id}:{chain_a}-{res_nr_a}:{chain_b}-{res_nr_b}')
                    pos_pos_output.add_value('Interaction Score', pos_pos_interaction_score)
                    pos_pos_f.write(pos_pos_output.pop_line())
                if inter_count == len(pos_pos_map[pos_a_db_id]):
                    count_2 += 1
            pos_pos_f.close()

            if config.verbosity >= 2:
                config.logger.info(f"After pos pos table write {count_1=} {count_2=}")

    generate_multi_sav_table(config, multi_sav_table_file, session_id, snv_position_dict, position_dict, protein_dict)

    if os.path.isfile(f'{feature_file}.lock'):
        os.remove(f'{feature_file}.lock')

    if os.path.isfile(f'{class_file}.lock'):
        os.remove(f'{class_file}.lock')

    if config.ppi_data_available:
        if os.path.isfile(f'{interface_file}.lock'):
            os.remove(f'{interface_file}.lock')

    if config.verbosity >= 2:
        t6 = time.time()
        config.logger.info(f"Time for classificationOutput part 6: {t6 - t5}")

    #writeStatFile(stat_file, position_list, {}, tag_map, stat_dict=stat_dict)

    if config.verbosity >= 2:
        t7 = time.time()
        config.logger.info(f"Time for classificationOutput part 7: {t7 - t6}")
        config.logger.info(f'=== Total time for classificationOutput: {t7-t0} ===')

    return alignment_annotations

