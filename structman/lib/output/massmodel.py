import os
import shutil

import ray

from structman.lib.database import database
from structman.lib.globalAlignment import mutate_aligned_sequence
from structman.lib.sdsc.sdsc_utils import invert_deletion_flanks, invert_insertion_positions
try:
    from structman.lib import modelling
except:
    pass

from structman.lib.output import out_generator

def get_optimal_templates(prot_id, proteins, config):

    best_interface_structures = {}
    for interface_number, (positions, interface_structure_map) in enumerate(proteins[prot_id].aggregated_interface_map):
        best_score = 0.
        best_interface_structure = None
        for structure_id in interface_structure_map:
            chain, i_chain, interface_size = interface_structure_map[structure_id]
            seq_id = proteins[prot_id].structure_annotations[(structure_id, chain)].sequence_identity
            cov = proteins[prot_id].structure_annotations[(structure_id, chain)].coverage
            score = interface_size * seq_id * cov
            if score > best_score:
                best_score = score
                best_interface_structure = (structure_id, chain)
        if not best_interface_structure in best_interface_structures:
            best_interface_structures[best_interface_structure] = [interface_number]
        else:
            best_interface_structures[best_interface_structure].append(interface_number)

    for position in proteins[prot_id].positions:
        for (structure_id, tchain) in proteins[prot_id].structure_annotations:
            sub_info = protein_obj.structure_annotations[(pdb_id, chain)].get_sub_info(pos)
            res_nr = sub_info[0]
    #TODO

# called by structman_main
def mass_model(session_id, config, outfolder, include_snvs=False, template_selection_scheme = "functional_annotation", with_multi_mutations = True,
                prot_db_ids = None, model_wts = False):
    current_chunk = 1
    proteins, prot_db_id_dict = database.proteinsFromDb(session_id, config, with_residues=True, with_mappings = True,
                                       with_snvs=True, mutate_snvs=include_snvs, with_alignments=True, with_multi_mutations = with_multi_mutations,
                                       with_complexes=True, keep_ray_alive=True, current_chunk = current_chunk, prot_db_id_dict = prot_db_ids)

    while proteins is not None:

        print(f'Mass Modelling Pipeline\n####################### Processing Chunk {current_chunk} ########################################\n')

        prot_ids = proteins.get_protein_ids()

        if config.proc_n <= 1:
            para = False
        else:
            para = True
        if para:
            modelling_result_ids = []
            conf_dump = ray.put(config)
        else:
            models = []

        model_jobs = {}

        for prot_id in prot_ids:

            print(f'Processing: {prot_id}')

            #rec_struct = proteins[prot_id].get_major_recommend_structure(prefer_mutation_positions = True)
            wildtype_protein = proteins[prot_id].wildtype_protein

            #We are more interested in the wildtype structure, especially parts that get deleted should be present in the model
            proteins[wildtype_protein].invert_mutation_positions(proteins[prot_id])
            rec_struct = proteins[wildtype_protein].get_recommended_structure(scheme = template_selection_scheme, prefer_mutation_positions = True)

            if rec_struct is None:
                print('No mutation rec struct:', prot_id, wildtype_protein)
                if wildtype_protein is None:
                    print('Error wt prot is None:', prot_id)
                rec_struct = proteins[prot_id].get_recommended_structure(scheme = template_selection_scheme, prefer_mutation_positions = True)
                #proteins[wildtype_protein].invert_mutation_positions(proteins[prot_id])
                #rec_struct = proteins[wildtype_protein].get_major_recommend_structure(prefer_mutation_positions = True)
                if rec_struct is None:
                    print('No rec struct:', prot_id, wildtype_protein)
                    continue
            pdb_id, tchain = rec_struct.split(':')
            compl_obj = proteins.complexes[pdb_id]
            structures = proteins.get_complex_structures(pdb_id)

            print(f'Recommended structure: {pdb_id}:{tchain}, Mutant Type: {proteins[prot_id].mutant_type}')

            """
            del_flanks = proteins[prot_id].deletion_flanks

            if proteins[prot_id].mutant_type == 'SAV':
                aligned_seq, aligned_temp_seq = proteins.get_alignment(wildtype_protein, pdb_id, tchain)
                pos = proteins[prot_id].sav_positions[0]
                mut_aa = proteins[prot_id].positions[pos].wt_aa
                aligned_seq = mutate_aligned_sequence(aligned_seq, pos, mut_aa)
                alignment_tuple = aligned_seq, aligned_temp_seq
                seq_id = proteins[wildtype_protein].structure_annotations[(pdb_id, tchain)].sequence_identity
                cov = proteins[wildtype_protein].structure_annotations[(pdb_id, tchain)].coverage
                wt_sav_positions = proteins[prot_id].sav_positions

                if wildtype_protein not in model_jobs:
                    model_jobs[wildtype_protein] = {}

                    print(f'Adding an additional wildtype model {wildtype_protein} for {prot_id} on template {pdb_id}:{tchain}')
                    model_jobs[wildtype_protein][pdb_id] = [compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, wildtype_protein, None, None, None]
                elif pdb_id not in model_jobs[wildtype_protein]:
                    print(f'Adding an additional wildtype model {wildtype_protein} for {prot_id} on template {pdb_id}:{tchain}')
                    model_jobs[wildtype_protein][pdb_id] = [compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, wildtype_protein, None, None, None]

                for wt_target_struct in model_jobs[wildtype_protein]:
                    if wt_sav_positions is not None and model_jobs[wildtype_protein][wt_target_struct][8] is not None:
                        model_jobs[wildtype_protein][wt_target_struct][8] += wt_sav_positions
                    elif wt_sav_positions is not None:
                        model_jobs[wildtype_protein][wt_target_struct][8] = wt_sav_positions

            else:
                
                try:
                    alignment_tuple = proteins.get_alignment(prot_id, pdb_id, tchain)
                    seq_id = proteins[prot_id].structure_annotations[(pdb_id, tchain)].sequence_identity
                    cov = proteins[prot_id].structure_annotations[(pdb_id, tchain)].coverage
                except:
                    continue

                if proteins[prot_id].mutant_type == 'indel':
                    #if rec_struct != proteins[wildtype_protein].get_major_recommend_structure():

                    wt_alignment_tuple = proteins.get_alignment(wildtype_protein, pdb_id, tchain)
                    if isinstance(wt_alignment_tuple, str):
                        print(f'Wildtype {wildtype_protein} was not mapped to {pdb_id} {tchain}:\n{wt_alignment_tuple}')
                        continue
                    if (pdb_id, tchain) in proteins[wildtype_protein].structure_annotations:
                        wt_seq_id = proteins[wildtype_protein].structure_annotations[(pdb_id, tchain)].sequence_identity
                        wt_cov = proteins[wildtype_protein].structure_annotations[(pdb_id, tchain)].coverage
                    else:
                        wt_seq_id = None
                        wt_cov = None

                    wt_inserts, wt_del_flanks = proteins[prot_id].get_mutation_position_inverts()

                    #print(f'wt_inserts of {prot_id}: {wt_inserts}')

                    if wildtype_protein not in model_jobs:
                        model_jobs[wildtype_protein] = {}

                        print(f'Adding an additional wildtype model {wildtype_protein} for {prot_id} on template {pdb_id}:{tchain}')
                        model_jobs[wildtype_protein][pdb_id] = [compl_obj, structures, wt_alignment_tuple, wt_seq_id, wt_cov, pdb_id, tchain, wildtype_protein, None, None, None]
                    elif pdb_id not in model_jobs[wildtype_protein]:
                        print(f'Adding an additional wildtype model {wildtype_protein} for {prot_id} on template {pdb_id}:{tchain}')
                        model_jobs[wildtype_protein][pdb_id] = [compl_obj, structures, wt_alignment_tuple, wt_seq_id, wt_cov, pdb_id, tchain, wildtype_protein, None, None, None]

                    for wt_target_struct in model_jobs[wildtype_protein]:

                        if wt_inserts is not None and model_jobs[wildtype_protein][wt_target_struct][9] is not None:
                            model_jobs[wildtype_protein][wt_target_struct][9] += wt_inserts
                        elif wt_inserts is not None:
                            model_jobs[wildtype_protein][wt_target_struct][9] = wt_inserts

                        if wt_del_flanks is not None and model_jobs[wildtype_protein][wt_target_struct][10] is not None:
                            model_jobs[wildtype_protein][wt_target_struct][10] += wt_del_flanks
                        elif wt_inserts is not None:
                            model_jobs[wildtype_protein][wt_target_struct][model_nr][10] = wt_del_flanks
            """


            if with_multi_mutations:
                mm_position_dict = proteins[prot_id].get_mutation_positions_by_mm()

                if config.verbosity >= 3:
                    print(f'mm_position_dict for {prot_id}\n{mm_position_dict}\n')

                for mut_prot_id in mm_position_dict:
                    sav_pos, inserts, del_flanks = mm_position_dict[mut_prot_id]

                    inv_inserts = invert_deletion_flanks(del_flanks)
                    inv_del_flanks = invert_insertion_positions(inserts)

                    #print('inserts', inserts)
                    #print('inv_inserts', inv_inserts)
                    #print('del_flanks', del_flanks)
                    #print('inv_del_flanks', inv_del_flanks)

                    mm_rec_struct = proteins[wildtype_protein].get_recommended_structure(scheme = template_selection_scheme, prefer_mutation_positions = True, custom_mutations = (sav_pos, inv_inserts, inv_del_flanks), functionally_weighted = True)

                    if mm_rec_struct is not None:
                        rec_struct = mm_rec_struct

                        pdb_id, tchain = rec_struct.split(':')
                        compl_obj = proteins.complexes[pdb_id]
                        structures = proteins.get_complex_structures(pdb_id)
                        try:
                            alignment_tuple = proteins.get_alignment(prot_id, pdb_id, tchain)
                            seq_id = proteins[prot_id].structure_annotations[(pdb_id, tchain)].sequence_identity
                            cov = proteins[prot_id].structure_annotations[(pdb_id, tchain)].coverage
                        except:
                            print(f'Error: getting the alignments failed in the modelling process of {prot_id} with highlighting {mut_prot_id}')
                            continue

                    if prot_id not in model_jobs:
                        model_jobs[prot_id] = {}
                    if not pdb_id in model_jobs[prot_id]:
                        model_jobs[prot_id][pdb_id] = []

                    model_jobs[prot_id][pdb_id].append([compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, sav_pos, inv_inserts, inv_del_flanks, mut_prot_id])


            if prot_id not in model_jobs or model_wts:
                model_jobs[prot_id] = {}
                if pdb_id not in model_jobs[prot_id]:
                    model_jobs[prot_id][pdb_id] = []
                model_jobs[prot_id][pdb_id].append([compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, [], [], [], ''])


        if not os.path.exists(outfolder):
            os.mkdir(outfolder)

        model_database = '%s/models' % outfolder
        summary_file = '%s/model_summary.tsv' % outfolder

        if not os.path.exists(model_database):
            os.mkdir(model_database)

        for prot_id in model_jobs:
            for pdb_id in model_jobs[prot_id]:
                for compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, sav_positions, insertion_positions, deletion_flanks, label_add in model_jobs[prot_id][pdb_id]:

                    print('Modelling of:', prot_id, pdb_id, tchain, label_add)

                    if para:
                        modelling_result_ids.append(modelling.model_remote_wrapper.remote(conf_dump, compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, label_add = label_add, skip_analysis=True, force_modelling = config.force_modelling, mutation_positions = [sav_positions, insertion_positions, deletion_flanks], target_path = model_database))
                    else:
                        models.append(modelling.model(config, compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, label_add = label_add, skip_analysis=True, force_modelling = config.force_modelling, mutation_positions = [sav_positions, insertion_positions, deletion_flanks], target_path = model_database))
        if para:
            models = ray.get(modelling_result_ids)

        if len(models) == 0:
            print('No models got produced :(')
            current_chunk += 1
            proteins, prot_db_id_dict = database.proteinsFromDb(session_id, config, with_residues=True,
                                       with_snvs=True, mutate_snvs=include_snvs, with_alignments=True,
                                       with_complexes=True, keep_ray_alive=True, current_chunk = current_chunk, prot_db_id_dict = prot_db_id_dict)
            continue

        if os.path.exists(summary_file) and current_chunk == 1:
            os.remove(summary_file)

        if current_chunk == 1:
            model_output = out_generator.OutputGenerator()
            headers = ['Input ID', 'Wildtype protein', 'Protein ID', 'Template PDB ID', 'Template chain',
                       'Model chain', 'Template resolution', 'Sequence identity', 'Coverage', 'File location']
            model_output.add_headers(headers)
            f = open(summary_file, 'a')
            f.write(model_output.get_header())
        else:
            f = open(summary_file, 'a')

        for model in models:
            if isinstance(model, str):
                config.errorlog.add_warning(model)
            else:

                # model.clear_tmp()
                pdb_id, tchain = model.template_structure
                model_chain = model.chain_id_map[tchain]
                label_add = model.label_add

                target_file_name = f'{model.truncated_prot_id}-{label_add}_{pdb_id}:{tchain}.pdb'

                target_path = f'{model_database}/{target_file_name}'
                shutil.copy(model.path, target_path)

                input_id = proteins[model.target_protein].input_id
                if input_id is None:
                    input_id = proteins[proteins[model.target_protein].wildtype_protein].input_id

                model_output.add_value('Input ID', input_id)
                model_output.add_value('Wildtype protein', proteins[model.target_protein].wildtype_protein)
                model_output.add_value('Protein ID', model.target_protein)
                model_output.add_value('Template PDB ID', pdb_id)
                model_output.add_value('Template chain', tchain)
                model_output.add_value('Model chain', model_chain)
                model_output.add_value('Template resolution', model.template_resolution)
                model_output.add_value('Sequence identity', model.sequence_identity)
                model_output.add_value('Coverage', model.coverage)
                model_output.add_value('File location', target_path)

                f.write(model_output.pop_line())

        f.close()

        current_chunk += 1
        proteins, prot_db_id_dict = database.proteinsFromDb(session_id, config, with_residues=True,
                                       with_snvs=True, mutate_snvs=include_snvs, with_alignments=True,
                                       with_complexes=True, keep_ray_alive=True, current_chunk = current_chunk, prot_db_id_dict = prot_db_id_dict)

    ray.shutdown()

