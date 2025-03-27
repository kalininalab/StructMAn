import os
import psutil
import ray
import time

from structman.lib import serializedPipeline
from structman.lib.output.out_generator import OutputGenerator
from structman.lib.database.insertion_lib import insert_indel_results
try:
    from structman.lib import modelling
except:
    pass


def separate_structure_annotations(wt_structure_annotations, mut_structure_annotations, config):
    sep_wt = {}
    sep_mut = {}
    for struct_tuple in wt_structure_annotations:
        if struct_tuple in mut_structure_annotations:
            if wt_structure_annotations[struct_tuple].sequence_identity is None:
                config.errorlog.add_error('sequence identity was None for: %s' % str(struct_tuple))
                continue
            # This is especially important to separate structures for substitution-type indels
            if wt_structure_annotations[struct_tuple].sequence_identity >= mut_structure_annotations[struct_tuple].sequence_identity:
                sep_wt[struct_tuple] = wt_structure_annotations[struct_tuple]
            else:
                sep_mut[struct_tuple] = mut_structure_annotations[struct_tuple]
            del mut_structure_annotations[struct_tuple]
        else:
            sep_wt[struct_tuple] = wt_structure_annotations[struct_tuple]

    for struct_tuple in mut_structure_annotations:
        sep_mut[struct_tuple] = mut_structure_annotations[struct_tuple]

    return sep_wt, sep_mut


def get_region_class(classification_results, l, r):
    if l not in classification_results:
        l = 1
    if r not in classification_results:
        r = len(classification_results) - 1
    if r < l:
        return None
    if r == l:
        return classification_results[l].rin_simple_class
    class_counts = {}
    for i in range(l, r + 1):
        rcl = classification_results[i].rin_simple_class
        if rcl is None:
            continue
        if rcl not in class_counts:
            class_counts[rcl] = 1
        else:
            class_counts[rcl] += 1

    max_count = 0
    max_rcl = None
    for rcl in class_counts:
        if class_counts[rcl] >= max_count:
            max_count = class_counts[rcl]
            max_rcl = rcl
    return max_rcl


neutral_classes = set(['Disorder', 'Surface'])


def flanking_region_analysis(classification_results, left_flank_pos, right_flank_pos, region_size=4):
    left_flank_region_class = get_region_class(classification_results, (left_flank_pos - region_size) + 1, left_flank_pos)
    right_flank_region_class = get_region_class(classification_results, right_flank_pos, (right_flank_pos + region_size) - 1)

    if left_flank_region_class is None:
        return None
    if right_flank_region_class is None:
        return None

    if left_flank_region_class == right_flank_region_class:
        if left_flank_region_class in neutral_classes:
            return 'Identical neutral'
        else:
            return 'Identical potent'
    else:
        if left_flank_region_class in neutral_classes:
            if right_flank_region_class in neutral_classes:
                return 'Transition neutral,neutral'
            else:
                return 'Transition neutral,potent'
        else:
            if right_flank_region_class in neutral_classes:
                return 'Transition neutral,potent'
            else:
                return 'Transition potent,potent'


@ray.remote(num_cpus=1)
def indel_analysis(config, indel_obj, wt_structure_annotations, mut_structure_annotations):
    serializedPipeline.ray_hack()

    # print('1:',indel_obj.wt_prot,indel_obj.get_notation(),len(wt_structure_annotations),len(mut_structure_annotations))

    # The structures has to be checked that they include the flanks
    wt_structure_annotations, mut_structure_annotations = indel_obj.inclusion_check(wt_structure_annotations, mut_structure_annotations)

    # print('2:',indel_obj.wt_prot,indel_obj.get_notation(),len(wt_structure_annotations),len(mut_structure_annotations))

    wt_structure_annotations, mut_structure_annotations = separate_structure_annotations(wt_structure_annotations, mut_structure_annotations, config)

    # print('3:',indel_obj.wt_prot,indel_obj.get_notation(),len(wt_structure_annotations),len(mut_structure_annotations))

    analysis_list = indel_obj.get_scenario(len(wt_structure_annotations), len(mut_structure_annotations))

    # print('4:',indel_obj.wt_prot,indel_obj.get_notation(),analysis_list,len(wt_structure_annotations),len(mut_structure_annotations))

    indel_type, terminal, terminal_end = indel_obj.get_type()

    wt_seq_ids = []
    mut_seq_ids = []
    for struct_tuple in wt_structure_annotations:
        wt_seq_ids.append(wt_structure_annotations[struct_tuple].get_sequence_id())
    for struct_tuple in mut_structure_annotations:
        mut_seq_ids.append(mut_structure_annotations[struct_tuple].get_sequence_id())

    if len(wt_seq_ids) > 0:
        avg_wt_seq_id = sum(wt_seq_ids) / len(wt_seq_ids)
    else:
        avg_wt_seq_id = None
    if len(mut_seq_ids) > 0:
        avg_mut_seq_id = sum(mut_seq_ids) / len(mut_seq_ids)
    else:
        avg_mut_seq_id = None

    annotations = {indel_obj.wt_prot: list(wt_structure_annotations.keys()), indel_obj.mut_prot: list(mut_structure_annotations.keys())}

    flanking_region_left_pos, flanking_region_right_pos = indel_obj.get_flanking_region_pos()

    ret_tuple = (
        indel_type, terminal, terminal_end, len(wt_structure_annotations), len(mut_structure_annotations),
        avg_wt_seq_id, avg_mut_seq_id, annotations, analysis_list, indel_obj.tags, flanking_region_left_pos, flanking_region_right_pos,
        indel_obj.wt_prot
    )

    return ret_tuple


def compare_classifications(proteins, config, prot_id, models, sub_infos, write_output=None):
    if write_output is not None:
        dC_output = OutputGenerator()
        headers = ['Position', 'Original classification', 'Recommended structure', 'Model ID', 'Model chain', 'Model residue', 'Model classification']
        dC_output.add_headers(headers)
        outfile = '%s/%s_dC_output_%s.tsv' % (config.outfolder, prot_id, write_output)
        if os.path.exists(outfile):
            os.remove(outfile)
        f = open(outfile, 'a')
        f.write(dC_output.get_header())

    delta_classification = 0
    rec_structs = proteins[prot_id].getRecommendedStructures()
    for rec_struct in rec_structs:
        pdb_id, tchain = rec_struct.split(':')

        if not (pdb_id, tchain) in models:
            continue
        for pos in rec_structs[rec_struct]:
            original_class = proteins.get_classification(prot_id, pos)
            model = models[(pdb_id, tchain)]
            model_chain = model.chain_id_map[tchain]
            sub_info = sub_infos[model.model_id]
            res_id = sub_info[pos][0]
            if write_output is not None:
                dC_output.add_value('Position', pos)
                dC_output.add_value('Original classification', original_class)
                dC_output.add_value('Recommended structure', rec_struct)
                dC_output.add_value('Model ID', model.model_id)
                dC_output.add_value('Model chain', model_chain)
                dC_output.add_value('Model residue', res_id)
            if res_id is None:
                delta_classification += 1
                if config.verbosity >= 2:
                    print('Position not annotated in model:', prot_id, pos, pdb_id, tchain, model.model_id, original_class)
                if write_output is not None:
                    dC_output.add_value('Model classification', 'Unmapped')
            else:
                if model_chain not in model.structural_analysis_dict:
                    config.errorlog.add_warning('Chain not in structural analysis of model: %s %s %s' % (model_chain, model.model_id, model.path))
                    model_class = None
                elif res_id not in model.structural_analysis_dict[model_chain]:
                    config.errorlog.add_warning('Res_id not in structural analysis of model: %s %s %s %s' % (res_id, model_chain, model.model_id, model.path))
                    model_class = None
                else:
                    model_class = model.structural_analysis_dict[model_chain][res_id].get_classification(config)[1]
                if write_output is not None:
                    dC_output.add_value('Model classification', model_class)
                if original_class != model_class:
                    delta_classification += 1
            if write_output is not None:
                f.write(dC_output.pop_line())

    f.close()

    return delta_classification


def para_indel_analysis(proteins, config):
    indel_result_ids = []
    modelling_result_ids = []

    prot_specific_mapping_dumps = {}

    mut_backmap = {}

    model_stack = set()

    if config.verbosity >= 2:
        print('Starting indel analysis with', len(proteins.indels), 'indels')
        t0 = time.time()

    complex_store_references = []
    structures_store_references = []

    conf_dump = None

    if config.model_indel_structures:
        conf_dump = ray.put(config)

        for prot_id in proteins.indels:
            seq = proteins.get_sequence(prot_id)
            aaclist = proteins.getAACList(prot_id)
            prot_specific_mapping_dumps[prot_id] = ray.put((prot_id, seq, aaclist))
            del seq
            del aaclist
            wt_structure_annotations = proteins.get_protein_structure_annotations(prot_id)
            if len(wt_structure_annotations) == 0:
                if config.verbosity >= 3:
                    print('Indel analysis not possible, no structures annotated', prot_id)
                continue

            rec_structs = proteins[prot_id].getRecommendedStructures()

            if len(rec_structs) == 0:
                if config.verbosity >= 3:
                    print('Indel analysis not possible, no recommended structures for', prot_id)
                continue

            for indel_notation in proteins.indels[prot_id]:
                indel_obj = proteins.indels[prot_id][indel_notation]

                mut_structure_annotations = proteins.get_protein_structure_annotations(indel_obj.mut_prot)
                #indel_result_ids.append(indel_analysis.remote(conf_dump, indel_obj, wt_structure_annotations, mut_structure_annotations))
                mut_backmap[indel_obj.mut_prot] = indel_obj.wt_prot

                for rec_struct in rec_structs:
                    pdb_id, tchain = rec_struct.split(':')
                    if not (pdb_id, tchain) in mut_structure_annotations:
                        if config.verbosity >= 3:
                            print('Indel analysis not possible, structure not annotated to mut protein', prot_id, indel_obj.wt_prot, pdb_id, tchain)
                        continue
                    # modelling the WT for all recommended structures
                    compl_obj = ray.put(proteins.complexes[pdb_id])
                    complex_store_references.append(compl_obj)
                    structures = ray.put(proteins.get_complex_structures(pdb_id))
                    structures_store_references.append(structures)

                    if not (pdb_id, tchain, indel_obj.wt_prot) in model_stack:
                        alignment_tuple = proteins.get_alignment(indel_obj.wt_prot, pdb_id, tchain)
                        seq_id = proteins[indel_obj.wt_prot].structure_annotations[(pdb_id, tchain)].sequence_identity
                        cov = proteins[indel_obj.wt_prot].structure_annotations[(pdb_id, tchain)].coverage
                        modelling_result_ids.append(modelling.model.remote(conf_dump, compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, indel_obj.wt_prot))

                        model_stack.add((pdb_id, tchain, indel_obj.wt_prot))

                    # modelling the MUT for all recommend structures
                    if not (pdb_id, tchain, indel_obj.mut_prot) in model_stack:
                        alignment_tuple = proteins.get_alignment(indel_obj.mut_prot, pdb_id, tchain)
                        seq_id = proteins[indel_obj.mut_prot].structure_annotations[(pdb_id, tchain)].sequence_identity
                        cov = proteins[indel_obj.mut_prot].structure_annotations[(pdb_id, tchain)].coverage
                        modelling_result_ids.append(modelling.model.remote(conf_dump, compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, indel_obj.mut_prot))

                        model_stack.add((pdb_id, tchain, indel_obj.mut_prot))

    if config.verbosity >= 2:
        t1 = time.time()
        print('Indel analysis, part 1:', t1 - t0)

    #indel_analysis_outs = ray.get(indel_result_ids)
    if config.model_indel_structures:
        models = ray.get(modelling_result_ids)

    del conf_dump
    for i in reversed(range(len(structures_store_references))):
        del structures_store_references[i]
    del structures_store_references
    for i in reversed(range(len(complex_store_references))):
        del complex_store_references[i]
    del complex_store_references

    if config.verbosity >= 2:
        t2 = time.time()
        print('Indel analysis, part 2:', t2 - t1)

    alignment_results = []
    model_map = {}

    if config.model_indel_structures:
        for model in models:
            if isinstance(model, str):
                config.errorlog.add_warning(model)
            else:
                if model.target_protein not in model_map:
                    model_map[model.target_protein] = {}
                model_map[model.target_protein][model.template_structure] = model
                pdb_id, tchain = model.template_structure
                if model.target_protein in mut_backmap:
                    align_prot = mut_backmap[model.target_protein]
                else:
                    align_prot = model.target_protein

                model_chain = model.chain_id_map[tchain]

                structure_infos = [(model.model_id, model_chain, [])]

                chunk =  [(prot_specific_mapping_dumps[align_prot], structure_infos)]

                alignment_results.append(serializedPipeline.align_remote_wrapper.remote(conf_dump, chunk, model_path=model.path))

    if config.verbosity >= 2:
        t3 = time.time()
        print('Indel analysis, part 3:', t3 - t2)

    if config.model_indel_structures:
        result_chunks = ray.get(alignment_results)

    for prot_id in list(prot_specific_mapping_dumps.keys()):
        del prot_specific_mapping_dumps[prot_id]
    del prot_specific_mapping_dumps

    if config.verbosity >= 2:
        t4 = time.time()
        print('Indel analysis, part 4:', t4 - t3)

    warn_map = set()
    sub_info_map = {}

    if config.model_indel_structures:

        for align_outs in result_chunks:
            for out in align_outs[0]:

                if len(out) == 3:
                    config.errorlog.add_error('Illegal alignment output: %s' % (str(out)))
                    continue

                if len(out) == 4:
                    (prot_id, pdb_id, chain, warn_text) = out
                    if prot_id not in warn_map:
                        config.errorlog.add_warning(warn_text)
                        warn_map.add(prot_id)
                    continue

                if len(out) == 5:
                    (prot_id, pdb_id, chain, sub_infos, align_info) = out
                    if config.verbosity >= 4:
                        seq_id, alignment_text = align_info
                        print('Alignment got filtered:', prot_id, pdb_id, chain, len(sub_infos), seq_id, alignment_text)
                    continue

                (prot_id, model_id, model_chain, alignment, seq_id, coverage, interaction_partners, chain_type_map,
                oligo, sub_infos, atom_count, last_residue, first_residue, chainlist, rare_residues) = out

                if prot_id not in sub_info_map:
                     sub_info_map[prot_id] = {}
                sub_info_map[prot_id][model_id] = sub_infos

    if config.verbosity >= 2:
        t5 = time.time()
        print('Indel analysis, part 5:', t5 - t4)

    n_filtered = 0
    n_success = 0

    for prot_id in proteins.indels:
        for indel_notation in proteins.indels[prot_id]:
            indel_obj = proteins.indels[prot_id][indel_notation]

            itype, terminal, terminal_end = indel_obj.get_type()
            indel_obj.aggregate(proteins, config)
            indel_obj.flank_aggregates(proteins, config)


            if itype == 'Substitution':
                pass #TODO
            elif itype == 'Insertion':
                pass #TODO
            else:
                pass

            if len(proteins.get_protein_structure_annotations(indel_obj.wt_prot)) == 0:
                n_filtered += 1
                continue

            if config.model_indel_structures:

                if indel_obj.wt_prot not in model_map:
                    config.errorlog.add_warning('Id not in model map: %s %s %s' % (prot_id, indel_obj.wt_prot, indel_notation))

                    continue
                if indel_obj.wt_prot not in sub_info_map:
                    config.errorlog.add_warning('Id not in sub info map: %s %s %s' % (prot_id, indel_obj.wt_prot, indel_notation))
                    continue
                wt_model_delta_classification = compare_classifications(proteins, config, indel_obj.wt_prot, model_map[indel_obj.wt_prot], sub_info_map[indel_obj.wt_prot], write_output='WT')
                mut_model_delta_classification = compare_classifications(proteins, config, indel_obj.wt_prot, model_map[indel_obj.mut_prot], sub_info_map[indel_obj.wt_prot], write_output='MUT')

                #print('ddC analysis:',prot_id,indel_notation,wt_model_delta_classification,mut_model_delta_classification)
                indel_obj.delta_delta_classification = mut_model_delta_classification - wt_model_delta_classification
            n_success += 1


    if config.verbosity >= 2:
        t6 = time.time()
        print('Indel analysis, part 6:', t6 - t5, n_filtered, n_success)

    insert_indel_results(proteins, config)

    if config.verbosity >= 2:
        t7 = time.time()
        print('Indel analysis, part 7:', t7 - t6)
