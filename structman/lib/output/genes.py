import os
import sys
import traceback
import math
import igraph as ig
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import ray

try:
    import biotite.sequence as seq
    import biotite.sequence.align as align
    import biotite.sequence.graphics as graphics
except:
    print('Warning: biotite imports failed, gene report functionalities will be impaired.')

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import math
import colorsys

import powerlaw
import numpy

import markdown
import pdfkit

from structman.lib.database import database
from structman.lib.globalAlignment import init_bp_aligner_class
from structman.lib.database.database_core_functions import select
from structman.lib.lib_utils import fuse_multi_mutations, generate_multi_mutation
from structman.lib.output import out_generator, out_utils, massmodel
from structman.base_utils.base_utils import resolve_path, calculate_chunksizes
from structman.lib.sdsc.mutations import MultiMutation

import structman.scripts.spectrumany as pymol_spectrumany

class Gene:

    def __init__(self, gene_id, name, isoform_pair_multimutations, isoform_pair_scores, isoform_expressions, isoform_prot_objs):
        self.gene_id = gene_id
        self.name = name
        isoforms = {}
        self.all_expression_values_none = True
        for iso in isoform_prot_objs:
            if iso in isoform_expressions:
                expression_wt_condition, expression_d_condition = isoform_expressions[iso]
                self.all_expression_values_none = False
            else:
                expression_wt_condition = None
                expression_d_condition = None
            isoforms[iso] = expression_wt_condition, expression_d_condition, isoform_prot_objs[iso]
        isoform_pairs = {}
        for iso_a, iso_b, multi_mutation in isoform_pair_multimutations:
            if (iso_a.primary_protein_id, iso_b.primary_protein_id) in isoform_pair_scores:
                SDS = isoform_pair_scores[(iso_a.primary_protein_id, iso_b.primary_protein_id)]
            else:
                SDS = None
            isoform_pairs[(iso_a.primary_protein_id, iso_b.primary_protein_id)] = (SDS, multi_mutation)
        self.isoforms = isoforms
        self.isoform_pairs = isoform_pairs

    def generate_input_for_sequence_plot(self, isoform_id, insignificant_isoforms):
        prot = self.isoforms[isoform_id][2]
        seq = prot.sequence
        classifications = prot.get_ordered_classifications()
        isoform_differences = {}
        isoform_sequences = {}
        isoform_classifications = {}
        for iso_a, iso_b in self.isoform_pairs:
            if iso_a != isoform_id:
                continue
            if iso_a in insignificant_isoforms:
                continue
            if iso_b in insignificant_isoforms:
                continue
            insertions, dels, subs = self.isoform_pairs[(iso_a, iso_b)][1].get_affected_positions()
            #print(iso_a, iso_b)
            #self.isoform_pairs[(iso_a, iso_b)][1].print_indels()
            isoform_differences[iso_b] = insertions, dels, subs
            isoform_sequences[iso_b] = self.isoforms[iso_b][2].sequence
            isoform_classifications[iso_b] = self.isoforms[iso_b][2].get_ordered_classifications()
        return seq, isoform_id, classifications, isoform_differences, isoform_sequences, isoform_classifications

def init_gene_isoform_table(out_file, obj_only = False):
    genes_output = out_generator.OutputGenerator()
    headers = [
        'Gene ID', 'Gene Isoform Score', 'Canonical Isoform', 'Weighted GI Score',
        'Highest Contribution Isoform', 'Contribution Score'
    ]


    genes_output.add_headers(headers)
    if obj_only:
        return genes_output

    if os.path.exists(out_file):
        os.remove(out_file)

    f = open(out_file, 'a')
    f.write(genes_output.get_header())

    return genes_output, f

def decideOnWildtype(protein_list, wt_condition_tag = None):
    wt_prot = None
    longest_seq_prot = None
    max_tag_value = None
    max_tag_prot = None
    for prot in protein_list:
        
        if longest_seq_prot is None:
            longest_seq_prot = prot
            wt_prot = prot
        elif len(prot.sequence) > len(longest_seq_prot.sequence):
            longest_seq_prot = prot

        if wt_condition_tag is not None:
            tag_value = prot.retrieve_tag_value(wt_condition_tag)
            if tag_value is not None:
                if max_tag_value is None:
                    max_tag_value = tag_value
                    max_tag_prot = prot
                elif tag_value > max_tag_value:
                    max_tag_value = tag_value
                    max_tag_prot = prot

    if wt_condition_tag is None:
        wt_prot = longest_seq_prot
    elif max_tag_prot is not None:
        wt_prot = max_tag_prot
    return wt_prot

def calcWeightedGI(isoforms, proteins, protein_list, wt_condition_tag, disease_condition_tag):
    gi_score_sum = 0.
    max_contribution_isoform = None
    max_contribution = 0.
    isoform_expressions = {}
    isoform_pair_scores = {}
    for (prot_a, prot_b, multi_mutation) in isoforms:
        expression_wt_condition_prot_a = prot_a.retrieve_tag_value(wt_condition_tag)
        expression_d_condition_prot_a = prot_a.retrieve_tag_value(disease_condition_tag)
        if prot_a.primary_protein_id not in isoform_expressions:
            isoform_expressions[prot_a.primary_protein_id] = expression_wt_condition_prot_a, expression_d_condition_prot_a
        expression_wt_condition_prot_b = prot_b.retrieve_tag_value(wt_condition_tag)
        expression_d_condition_prot_b = prot_b.retrieve_tag_value(disease_condition_tag)

        if expression_wt_condition_prot_a is None:
            isoform_pair_scores[(prot_a.primary_protein_id, prot_b.primary_protein_id)] = 0.
            continue
        if expression_d_condition_prot_a is None:
            isoform_pair_scores[(prot_a.primary_protein_id, prot_b.primary_protein_id)] = 0.
            continue
        if expression_wt_condition_prot_b is None:
            isoform_pair_scores[(prot_a.primary_protein_id, prot_b.primary_protein_id)] = 0.
            continue
        if expression_d_condition_prot_b is None:
            isoform_pair_scores[(prot_a.primary_protein_id, prot_b.primary_protein_id)] = 0.
            continue

        try:
            unweighted_score = multi_mutation.get_score(proteins)
        except:
            unweighted_score = 0.
        isoform_pair_scores[(prot_a.primary_protein_id, prot_b.primary_protein_id)] = unweighted_score
        if expression_wt_condition_prot_a <= expression_d_condition_prot_a: #prot_a needs to decrease in expression
            continue        
        if expression_wt_condition_prot_b >= expression_d_condition_prot_b: #prot_b needs to increase in expression
            continue
        weighted_score = unweighted_score * (expression_wt_condition_prot_a - expression_d_condition_prot_a) * (expression_d_condition_prot_b - expression_wt_condition_prot_b)
        gi_score_sum += weighted_score
        if weighted_score > max_contribution:
            max_contribution = weighted_score
            max_contribution_isoform = prot_b.primary_protein_id

    total_expression = 0.
    for prot in protein_list:
        expression_wt_condition = prot.retrieve_tag_value(wt_condition_tag)
        expression_d_condition = prot.retrieve_tag_value(disease_condition_tag)
        if expression_wt_condition is None:
            continue
        if expression_d_condition is None:
            continue
        total_expression += expression_wt_condition + expression_d_condition

    if total_expression == 0:
        weighted_gi_score = None
        normalized_max_contribution = None
        max_contribution_isoform = None
    else:
        weighted_gi_score = gi_score_sum/total_expression
        normalized_max_contribution = max_contribution/total_expression    

    return weighted_gi_score, normalized_max_contribution, max_contribution_isoform, isoform_pair_scores, isoform_expressions

    
def calculateGeneIsoformScore(isoforms, proteins, wt = None, wt_condition_tag = None, disease_condition_tag = None):
    if wt is None or disease_condition_tag is not None:
        protein_set = set()
        for (prot_a, prot_b, multi_mutation) in isoforms:
            protein_set.add(prot_a)
            protein_set.add(prot_b)
        protein_list = list(protein_set)
    if wt is None:
        wt = decideOnWildtype(protein_list, wt_condition_tag = wt_condition_tag)

    #print(f'Wildtype in calculateGeneIsoformScore: {wt.primary_protein_id}')

    if disease_condition_tag is not None:
        weighted_gi_score, normalized_max_contribution, max_contribution_isoform, isoform_pair_scores, isoform_expressions = calcWeightedGI(isoforms, proteins, protein_list, wt_condition_tag, disease_condition_tag)
        #weighted_gi_score, normalized_max_contribution, max_contribution_isoform = calcWeightedGI(isoforms, proteins, wt, disease_condition_tag)
        
    else:
        weighted_gi_score = None
        normalized_max_contribution = None
        max_contribution_isoform = None
        isoform_pair_scores = {}
        isoform_expressions = {}

    wt_reduced_isoforms = []
    wt_reduced_multi_mutations = []
    for (prot_a, prot_b, multi_mutation) in isoforms:
        if prot_a.primary_protein_id != wt.primary_protein_id:
            continue
        wt_reduced_isoforms.append((prot_a, prot_b, multi_mutation))
        wt_reduced_multi_mutations.append(multi_mutation)
        #print(f'Isoform: {prot_b.primary_protein_id} with indels:')
        #multi_mutation.print_indels()

    if len(wt_reduced_multi_mutations) == 0:
        GI_score = None
    else:
        fused_multi_mutation = fuse_multi_mutations(wt_reduced_multi_mutations)

        #print(f'Fused Isoform:')
        #fused_multi_mutation.print_indels()

        try:
            GI_score = fused_multi_mutation.get_score(proteins)
        except:
            GI_score = None

    return GI_score, wt.primary_protein_id, weighted_gi_score, max_contribution_isoform, normalized_max_contribution, isoform_pair_scores, isoform_expressions

def init_gene_maps(session_id, config, proteins):
    table = 'RS_Isoform'
    rows = ['Protein_A', 'Protein_B', 'Gene', 'Multi_Mutation']
    eq_rows = {'Session': session_id}

    results = select(config, rows, table, equals_rows=eq_rows)

    gene_db_ids = set()
    gene_isoform_map = {}
    completeness = {}
    for row in results:
        prot_a_db_id = row[0]
        prot_b_db_id = row[1]
        gene_db_id = row[2]
        multi_mutation_db_id = row[3]

        gene_db_ids.add(gene_db_id)
        if not gene_db_id in gene_isoform_map:
            gene_isoform_map[gene_db_id] = []
            completeness[gene_db_id] = {}
        gene_isoform_map[gene_db_id].append((prot_a_db_id, prot_b_db_id, multi_mutation_db_id, None))
        completeness[gene_db_id][prot_a_db_id] = [prot_b_db_id]

    gene_id_map = database.retrieve_gene_id_map(gene_isoform_map.keys(), config)
    aligner_class = init_bp_aligner_class()
    processed = set()
    gene_id_map_updates = []
    del_list = []
    for gene_db_id_a in gene_id_map:
        if gene_db_id_a not in gene_isoform_map:
            continue
        gene_id_a, gene_name_a = gene_id_map[gene_db_id_a]

        isoforms = set()
        for (prot_a_db_id, _, _, _) in gene_isoform_map[gene_db_id_a]:
            isoforms.add(prot_a_db_id)

        for iso_a in isoforms:
            for iso_b in isoforms:
                if iso_a == iso_b:
                    continue
                if iso_b not in completeness[gene_db_id_a][iso_a]:
                    prot_a = proteins.getByDbId(iso_a)
                    prot_b = proteins.getByDbId(iso_b)
                    multi_mutation_obj = generate_new_multi_mutation(proteins, config, aligner_class, prot_a, prot_b)
                    gene_isoform_map[gene_db_id_a].append((iso_a, iso_b, None, multi_mutation_obj))

        if config.verbosity >= 3:
            print(f'Size of gene_isoform_map of {gene_id_a}: {len(gene_isoform_map[gene_db_id_a])}')
        for gene_db_id_b in gene_id_map:
            if gene_db_id_b not in gene_isoform_map:
                continue
            if gene_db_id_b in processed:
                continue
            gene_id_b, gene_name_b = gene_id_map[gene_db_id_b]
            if gene_id_a == gene_id_b:
                continue
            if gene_name_a == gene_name_b:
                fused_isoform_pairs = fuse_genes(gene_isoform_map[gene_db_id_a], gene_isoform_map[gene_db_id_b], config, proteins, aligner_class)
                gene_isoform_map[gene_db_id_a] = fused_isoform_pairs
                del_list.append(gene_db_id_b)

                fused_gene_id = f'{gene_id_a}|{gene_id_b}'
                gene_id_map_updates.append((gene_db_id_a, fused_gene_id, gene_name_a))

        processed.add(gene_db_id_a)

    for gene_db_id, fused_gene_id, gene_name in gene_id_map_updates:
        gene_id_map[gene_db_id] = fused_gene_id, gene_name

    for gene_db_id in del_list:
        del gene_isoform_map[gene_db_id]

    return gene_isoform_map, gene_id_map

def fuse_genes(isoform_pairs_a, isoform_pairs_b, config, proteins, aligner_class):
    isoforms_a = set()
    for (prot_a_db_id, _, _, _) in isoform_pairs_a:
        isoforms_a.add(prot_a_db_id)

    isoforms_b = set()
    for (prot_a_db_id, _, _, _) in isoform_pairs_b:
        isoforms_b.add(prot_a_db_id)

    fused_isoform_pairs = isoform_pairs_a + isoform_pairs_b

    for prot_a_db_id in isoforms_a:
        for prot_b_db_id in isoforms_b:
            prot_a = proteins.getByDbId(prot_a_db_id)
            prot_b = proteins.getByDbId(prot_b_db_id)

            multi_mutation_obj = generate_new_multi_mutation(proteins, config, aligner_class, prot_a, prot_b)
            fused_isoform_pairs.append((prot_a_db_id, prot_b_db_id, None, multi_mutation_obj))

            multi_mutation_obj = generate_new_multi_mutation(proteins, config, aligner_class, prot_b, prot_a)
            fused_isoform_pairs.append((prot_b_db_id, prot_a_db_id, None, multi_mutation_obj))

    return fused_isoform_pairs

def generate_new_multi_mutation(proteins, config, aligner_class, prot_a, prot_b):
    seq_a = prot_a.sequence
    seq_b = prot_b.sequence
    _, multi_mutation_list = generate_multi_mutation(seq_a, seq_b, config, aligner_class = aligner_class)
    multi_mutation_obj = MultiMutation(prot_a.primary_protein_id, f'{prot_a.primary_protein_id}_{prot_b.primary_protein_id}', multi_mutation_list)
    
    prot_a.multi_mutations.append((multi_mutation_list, prot_b.primary_protein_id, None))
    if not prot_a.primary_protein_id in proteins.multi_mutations:
        proteins.multi_mutations[prot_a.primary_protein_id] = []
    proteins.multi_mutations[prot_a.primary_protein_id].append(multi_mutation_obj)
    return multi_mutation_obj

def calculate_gene_isoform_score_dict(session_id, config, proteins = None, create_gene_class = False, target_genes = None, return_prot_db_ids = False):
    if proteins is None:
        proteins = database.proteinsFromDb(session_id, config, with_multi_mutations = True, with_mappings = True, with_alignments=True)

    gene_isoform_map, gene_id_map = init_gene_maps(session_id, config, proteins)

    gene_isoform_score_dict = {}

    prot_db_ids = set()

    overwrite_targets = []

    for gene_db_id in gene_isoform_map:
        gene_id, gene_name = gene_id_map[gene_db_id]

        if target_genes is not None:
            if gene_id not in target_genes:
                continue
        else:
            overwrite_targets.append(gene_id)

        isoforms = []
        isoform_prot_objs = {}
        for (prot_a_db_id, prot_b_db_id, multi_mutation_db_id, multi_mutation) in gene_isoform_map[gene_db_id]:
            prot_a = proteins.getByDbId(prot_a_db_id)
            prot_b = proteins.getByDbId(prot_b_db_id)
            if multi_mutation_db_id is not None:
                multi_mutation = proteins.multi_mutation_back_map[multi_mutation_db_id]
            
            isoforms.append((prot_a, prot_b, multi_mutation))
            if prot_a.primary_protein_id not in isoform_prot_objs:
                isoform_prot_objs[prot_a.primary_protein_id] = prot_a
                if return_prot_db_ids:
                    prot_db_ids.add(prot_a_db_id)
            

        if config.verbosity >= 3:
            print(f'Call of calculateGeneIsoformScore for {gene_id} {len(isoforms)}\n')

        gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score, isoform_pair_scores, isoform_expressions = calculateGeneIsoformScore(isoforms, proteins, wt_condition_tag = config.condition_1_tag, disease_condition_tag = config.condition_2_tag)

        if config.verbosity >= 3:
            print(f'GI Score for {gene_id}: {gi_score} {len(isoform_pair_scores)}')

        if create_gene_class:
            gene_class = Gene(gene_id, gene_name, isoforms, isoform_pair_scores, isoform_expressions, isoform_prot_objs)
            gene_isoform_score_dict[gene_id] = gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score, gene_class
        else:
            gene_isoform_score_dict[gene_id] = gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score

    if len(overwrite_targets) > 0:
        target_genes = overwrite_targets

    if return_prot_db_ids:
        return gene_isoform_score_dict, target_genes, prot_db_ids, proteins

    return gene_isoform_score_dict, target_genes

def generate_gene_isoform_scores(session_id, config, outfolder, session_name):
    outfile = out_utils.generate_gis_filename(outfolder, session_name)

    gene_isoform_score_dict, _ = calculate_gene_isoform_score_dict(session_id, config)

    genes_output, f = init_gene_isoform_table(outfile)

    for gene_name in gene_isoform_score_dict:
        gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score = gene_isoform_score_dict[gene_name]

        genes_output.add_value('Gene ID', gene_name)
        genes_output.add_value('Gene Isoform Score', gi_score)
        genes_output.add_value('Canonical Isoform', wt)
        genes_output.add_value('Weighted GI Score', weighted_gi_score)
        genes_output.add_value('Highest Contribution Isoform', max_contribution_isoform)
        genes_output.add_value('Contribution Score', contribution_score)        
        f.write(genes_output.pop_line())

    f.close()

def make_barplot(title, label_a, label_b, value_a, value_b, outfile, color = 'grey', measure_label = 'Transcript Usage'):
    x_pos = np.arange(2)
    plt.bar(x_pos, [value_a, value_b], color = color)
    plt.title(title)
    plt.ylabel(measure_label)
    plt.xticks(x_pos, [label_a, label_b])
    plt.savefig(outfile)
    plt.clf()
    plt.cla()
    plt.close()


def create_isoform_relations_plot(config, gene_name, gene_isoform_score_dict, cmap, outfile_stem, condition_1_label = 'ExC1', condition_2_label = 'ExC2'):
    if config.verbosity >= 3:
        print(f'Creating isoform relations plot for: {gene_name}')

    gene_class = gene_isoform_score_dict[gene_name][5]

    if gene_class.all_expression_values_none:
        print(f'Creating isoform relations invalid for {gene_class.gene_id}, expression values are not defined. Please set the condition tags with --condition_1 "tag name" and --condition_2 "tag name"')
        return None

    barplot_folder = f'{outfile_stem}_barplots'
    if not os.path.isdir(barplot_folder):
        try:
            os.mkdir(barplot_folder)
        except:
            pass

    relation_plot_folder = f'{outfile_stem}_relationplots'
    if not os.path.isdir(relation_plot_folder):
        os.mkdir(relation_plot_folder)    

    node_id_list = []
    node_ids = {}
    node_labels = []
    max_ex_change = 0
    min_ex_change = 0
    raw_ex_changes = {}
    max_seq_len = 0

    total_expression_condition_1 = 0
    total_expression_condition_2 = 0
    for isoform in gene_class.isoforms:
        expression_wt_condition, expression_d_condition, prot_object = gene_class.isoforms[isoform]
        if (expression_d_condition == 0 and expression_wt_condition == 0) or expression_wt_condition is None or expression_d_condition is None:
            continue
        total_expression_condition_1 += expression_wt_condition
        total_expression_condition_2 += expression_d_condition

    path_to_gene_expression_plot = f'{barplot_folder}/{gene_name}_total_expression_plot.png'
    make_barplot(gene_name, 'Condition 1', 'Condition 2', total_expression_condition_1, total_expression_condition_2, path_to_gene_expression_plot, color = 'grey', measure_label = 'Total Expression')

    identicals = set()
    insignificant_isoforms = set()
    top_4_isoforms = []
    longest_isoform = None
    usage_plots = {}
    for isoform in gene_class.isoforms:

        expression_wt_condition, expression_d_condition, prot_object = gene_class.isoforms[isoform]
        if (expression_d_condition == 0 and expression_wt_condition == 0) or expression_wt_condition is None or expression_d_condition is None:
            insignificant_isoforms.add(isoform)
            continue
        if isoform in identicals:
            continue

        fused_label = isoform
        fused_label_short = 0
        for isoform_b in gene_class.isoforms:
            if isoform_b == isoform:
                continue

            try:
                _, multi_mutation = gene_class.isoform_pairs[(isoform, isoform_b)]
            except:
                try:
                    _, multi_mutation = gene_class.isoform_pairs[(isoform_b, isoform)]
                except:
                    continue


            l_deleted, l_inserted, l_substituted = multi_mutation.give_count()

            if l_deleted == 0 and l_inserted == 0 and l_substituted == 0:
                identicals.add(isoform_b)
                fused_label = f'{fused_label}\n{isoform_b}'
                fused_label_short += 1
                expression_wt_condition_b, expression_d_condition_b, _ = gene_class.isoforms[isoform_b]
                expression_wt_condition += expression_wt_condition_b
                expression_d_condition += expression_d_condition_b

        ex_change = expression_wt_condition - expression_d_condition
        if ex_change > max_ex_change:
            max_ex_change = ex_change
        if ex_change < min_ex_change:
            min_ex_change = ex_change

        raw_ex_changes[isoform] = ex_change

        if expression_wt_condition > 0 and expression_d_condition > 0:
            log_fold_change = round(math.log2(expression_d_condition/expression_wt_condition), 2)
        else:
            log_fold_change = 'inf'
        usage_condition_1 = round(expression_wt_condition/total_expression_condition_1, 2)
        usage_condition_2 = round(expression_d_condition/total_expression_condition_2, 2)

        combined_usage = usage_condition_1 + usage_condition_2
        #Filter out insignificant isoforms
        if combined_usage < 0.01:
            insignificant_isoforms.add(isoform)
            continue

        top_4_isoforms.append((isoform, combined_usage))

        if fused_label_short == 0:
            fused_label_short = isoform
        else:
            fused_label_short = f'{isoform}_and_{fused_label_short}_others'

        outfile = f'{barplot_folder}/{gene_name}_{fused_label_short}_usage_plot.png'
        make_barplot(fused_label, 'Condition 1', 'Condition 2', usage_condition_1, usage_condition_2, outfile, color = 'grey', measure_label = 'Transcript Usage')
        usage_plots[isoform] = outfile

        usage_change = round(usage_condition_2 - usage_condition_1, 2)

        node_id_list.append(isoform)
        node_ids[isoform] = len(node_ids)
        seq_len = len(prot_object.sequence)
        if seq_len > max_seq_len:
            max_seq_len = seq_len
            longest_isoform = isoform
        label = f'{fused_label}, Size: {seq_len}\n{condition_1_label}: {round(expression_wt_condition, 5)}\n{condition_2_label}: {round(expression_d_condition, 5)}\nl2fc: {log_fold_change}\nUsage {condition_1_label}: {usage_condition_1}\n Usage {condition_2_label}: {usage_condition_2}\nChange in usage: {usage_change}'
        #print(label)
        node_labels.append(label)

    #print(list(gene_class.isoforms.keys()))

    top_4_isoforms = sorted(top_4_isoforms, key=lambda x:x[1], reverse=True)[:4]
    top_4_isoforms = set([x[0] for x in top_4_isoforms])
    #print(top_4_isoforms)
    for isoform in gene_class.isoforms:
        if isoform in identicals:
            insignificant_isoforms.add(isoform)
        if isoform not in top_4_isoforms:
            insignificant_isoforms.add(isoform)

    #print(insignificant_isoforms)

    node_colors = []
    for isoform in node_id_list:
        ex_change = raw_ex_changes[isoform]
        if ex_change < 0:
            scaled_ex_change = ex_change/abs(min_ex_change)
        else:
            scaled_ex_change = ex_change/max_ex_change

        scaled_ex_change = (scaled_ex_change + 1)/2
        node_colors.append(cmap(scaled_ex_change))

    edge_list = []
    edge_labels = []
    edge_widths = []
    edge_colors = []
    identicals = set()
    max_sds = 0
    uninterresting_isoform_pairs = []
    for (iso_a, iso_b) in gene_class.isoform_pairs:
        #print(iso_a, iso_b)
        if iso_a not in node_ids:
            uninterresting_isoform_pairs.append((iso_a, iso_b))
            continue
        if iso_b not in node_ids:
            uninterresting_isoform_pairs.append((iso_a, iso_b))
            continue
        if iso_a in insignificant_isoforms:
            uninterresting_isoform_pairs.append((iso_a, iso_b))
            continue
        if iso_b in insignificant_isoforms:
            uninterresting_isoform_pairs.append((iso_a, iso_b))
            continue

        expression_wt_condition_a, expression_d_condition_a, _ = gene_class.isoforms[iso_a]
        expression_wt_condition_b, expression_d_condition_b, _ = gene_class.isoforms[iso_b]

        #Show only edges from higher expressed isoforms to lower expressed isoforms
        if expression_wt_condition_b + expression_d_condition_b > expression_wt_condition_a + expression_d_condition_a:
            continue

        SDS, multi_mutation = gene_class.isoform_pairs[(iso_a, iso_b)]

        l_deleted, l_inserted, l_substituted = multi_mutation.give_count()

        if l_deleted == 0 and l_inserted == 0 and l_substituted == 0:
            if (iso_b, iso_a) in identicals:
                continue
            edge_list.append((node_ids[iso_a], node_ids[iso_b]))
            edge_label = 'Identical protein'
            edge_labels.append(edge_label)
            edge_widths.append('max')
            edge_colors.append('blue')
            identicals.add((iso_a, iso_b))
            uninterresting_isoform_pairs.append((iso_a, iso_b))
            continue

        #if raw_ex_changes[iso_a] < 0 or raw_ex_changes[iso_b] > 0:
        #    continue

        if SDS == 0 and len(node_labels) > 6:
            uninterresting_isoform_pairs.append((iso_a, iso_b))
            continue

        edge_list.append((node_ids[iso_a], node_ids[iso_b]))
        edge_label = f'SDS:{round(SDS, 2)}\nDel:{l_deleted}\nIns:{l_inserted}\nSubs:{l_substituted}'
        edge_labels.append(edge_label)
        edge_widths.append(SDS)
        edge_colors.append('black')
        if SDS > max_sds:
            max_sds = SDS

    if max_sds == 0:
        max_sds = 1
    for pos, ew in enumerate(edge_widths):
        if ew == 'max':
            edge_widths[pos] = max_sds
        else:
            edge_widths[pos] = ew

    g = ig.Graph(n = len(node_ids), edges = edge_list, directed = True)
    g.vs['name'] = node_id_list
    g.vs['label'] = node_labels

    g.es['label'] = edge_labels

    node_size = 100
    margin = node_size
    graph_size = ((2*math.sqrt(len(node_id_list))*node_size) // 1) + 2*margin

    layout = g.layout('circle')
    visual_style = {}
    visual_style['layout'] = layout
    visual_style['margin'] = margin
    visual_style['bbox'] = (graph_size, graph_size)

    visual_style['vertex_label'] = g.vs['label']
    visual_style['vertex_label_size'] = node_size // 15
    visual_style['vertex_size'] = node_size
    visual_style['vertex_color'] = node_colors

    visual_style['edge_label'] = g.es['label']
    visual_style['edge_width'] = [((0.9*x)/max_sds) + 0.1 for x in edge_widths]
    visual_style['edge_label_size'] = node_size // 15
    visual_style['edge_color'] = edge_colors
    
    svg_file = f'{relation_plot_folder}/{gene_name}_isoform_relations.svg'

    ig.plot(g, svg_file, **visual_style)

    return svg_file, longest_isoform, uninterresting_isoform_pairs, path_to_gene_expression_plot, usage_plots, insignificant_isoforms


def generate_isoform_relations(session_id, config, outfolder, session_name):
    outfile_stem = os.path.join(outfolder, f'{session_name}')

    gene_isoform_score_dict, _ = calculate_gene_isoform_score_dict(session_id, config, create_gene_class = True)

    cmap=mpl.colormaps['RdYlGn']

    for gene_name in gene_isoform_score_dict:
        create_isoform_relations_plot(config, gene_name, gene_isoform_score_dict, cmap, outfile_stem)
       


def parse_gis_file(gis_file):
    f = open(gis_file, 'r')
    lines = f.readlines()
    f.close()

    gene_isoform_score_dict = {}

    for line in lines[1:]:
        words = line[:-1].split('\t')
        gene_name = words[0]
        gi_score = float(words[1])
        wt = words[2]
        try:
            weighted_gi_score = float(words[3])
        except:
            weighted_gi_score = None
        max_contribution_isoform = words[4]
        try:
            contribution_score = float(words[5])
        except:
            contribution_score = None

        gene_isoform_score_dict[gene_name] = gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score

    return gene_isoform_score_dict

def parse_ppi_file(ppi_file):
    f = open(ppi_file, 'r')
    lines = f.readlines()
    f.close()

    gene_gene_network = {}

    for line in lines[1:]:
        words = line.split('\t')
        (
            gene_a, prot_a_input_id, prot_a_primary_id, prot_a_u_ac, interface_number_a,
            gene_b, prot_b_input_id, prot_b_primary_id, prot_b_u_ac, interface_number_b,
            recommended_structure_str, interaction_score_str
        ) = words
        if (gene_a, gene_b) not in gene_gene_network:
            gene_gene_network[(gene_a, gene_b)] = {}
        gene_gene_network[(gene_a, gene_b)][(prot_a_primary_id, prot_b_primary_id)] = float(interaction_score_str)
    return gene_gene_network

def get_weight(prot, tag_label):
    if tag_label is None:
        return 1.0
    weight = prot.retrieve_tag_value(tag_label)
    if weight is None:
        weight = 1.0
    return weight

def generate_gene_interaction_network(session_id, config, outfolder, session_name):
    outfile_stem = os.path.join(outfolder, f'{session_name}_gene_gene_network')
    ppi_file = out_utils.generate_ppi_filename(outfolder, session_name)
    gis_file = out_utils.generate_gis_filename(outfolder, session_name)
    gml_file = f'{outfile_stem}.gml'

    gene_gene_network = parse_ppi_file(ppi_file)

    if os.path.isfile(gis_file):
        proteins = database.proteinsFromDb(session_id, config, with_multi_mutations = True)
        gene_isoform_score_dict = parse_gis_file(gis_file)
    else:
        proteins = database.proteinsFromDb(session_id, config, with_multi_mutations = True, with_mappings = True)
        gene_isoform_score_dict, _ = calculate_gene_isoform_score_dict(session_id, config, proteins = proteins)

    node_ids = {}
    node_id_list = []
    node_score_list = []
    edge_list = []
    edge_weights = []
    edge_weights_d_condition = []
    edge_weight_changes = []

    max_gi_score = 0.
    max_edge_weight = 0.
    max_edge_weight_change = 0.
    min_edge_weight_change = 0.
    degree_map = {}

    undirected_check = set()

    for (gene_a, gene_b) in gene_gene_network:
        if (gene_b, gene_a) in undirected_check:
            continue
        undirected_check.add((gene_a, gene_b))
        total_score_wt_condition = 0.
        wt_weight_sum = 0.
        total_score_d_condition = 0.
        d_weight_sum = 0.
        for (prot_a_primary_id, prot_b_primary_id) in gene_gene_network[(gene_a, gene_b)]:
            raw_interaction_score = gene_gene_network[(gene_a, gene_b)][(prot_a_primary_id, prot_b_primary_id)]
            prot_a = proteins[prot_a_primary_id]
            prot_b = proteins[prot_b_primary_id]
            
            wt_weight_a = get_weight(prot_a, config.condition_1_tag)
            wt_weight_b = get_weight(prot_b, config.condition_1_tag)
            wt_weight = wt_weight_a * wt_weight_b
            wt_weight_sum += wt_weight
            d_weight_a = get_weight(prot_a, config.condition_2_tag)
            d_weight_b = get_weight(prot_b, config.condition_2_tag)
            d_weight = d_weight_a * d_weight_b
            d_weight_sum += d_weight

            total_score_wt_condition += wt_weight * raw_interaction_score
            total_score_d_condition += d_weight * raw_interaction_score

        if wt_weight_sum > 0:
            weighted_wt_score = total_score_wt_condition / wt_weight_sum
        else:
            weighted_wt_score = 0.

        if d_weight_sum > 0:
            weighted_d_score = total_score_d_condition / d_weight_sum
        else:
            weighted_d_score = 0.

        weighted_wt_score = max([0, weighted_wt_score])
        weighted_d_score = max([0, weighted_d_score])

        if gene_a not in node_ids:
            node_ids[gene_a] = len(node_ids)
            node_id_list.append(gene_a)
            if gene_a in gene_isoform_score_dict:
                gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score = gene_isoform_score_dict[gene_a]
            else:
                weighted_gi_score = 0.

            if weighted_gi_score is None:
                weighted_gi_score = 0.

            node_score_list.append(weighted_gi_score)

            if weighted_gi_score > max_gi_score:
                max_gi_score = weighted_gi_score

        if gene_b not in node_ids:
            node_ids[gene_b] = len(node_ids)
            node_id_list.append(gene_b)
            if gene_b in gene_isoform_score_dict:
                gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score = gene_isoform_score_dict[gene_b]
            else:
                weighted_gi_score = 0.

            if weighted_gi_score is None:
                weighted_gi_score = 0.

            node_score_list.append(weighted_gi_score)

            if weighted_gi_score > max_gi_score:
                max_gi_score = weighted_gi_score
        
        node_id_a = node_ids[gene_a]
        node_id_b = node_ids[gene_b]

        edge = (node_id_a, node_id_b)
        edge_list.append(edge)

        if not node_id_a in degree_map:
            degree_map[node_id_a] = 0
        if not node_id_b in degree_map:
            degree_map[node_id_b] = 0
        degree_map[node_id_a] += 1
        degree_map[node_id_b] += 1

        if weighted_wt_score > 0:
            if weighted_wt_score >= weighted_d_score:
                perc_change = (100 * (weighted_d_score - weighted_wt_score)) / weighted_wt_score
            else:
                perc_change = (100 * (weighted_d_score - weighted_wt_score)) / weighted_d_score

            edge_weights.append(weighted_wt_score)
            edge_weights_d_condition.append(weighted_d_score)
            if weighted_wt_score > max_edge_weight:
                max_edge_weight = weighted_wt_score

            edge_weight_changes.append(perc_change)
            if perc_change > max_edge_weight_change:
                max_edge_weight_change = perc_change
            if perc_change < min_edge_weight_change:
                min_edge_weight_change = perc_change

        else:
            perc_change = None

            edge_weights.append(weighted_wt_score)
            edge_weights_d_condition.append(weighted_d_score)
            if weighted_wt_score > max_edge_weight:
                max_edge_weight = weighted_wt_score

            edge_weight_changes.append(perc_change)

    cmap=mpl.colormaps['RdYlGn']
    bound = max([max_edge_weight_change, abs(min_edge_weight_change)])

    edge_colors = []
    edge_color_codes = []
    scaled_ewgs = []
    for pos, ewg in enumerate(edge_weight_changes):
        if ewg is None: #new edges
            if edge_weights_d_condition[pos] == 0: #Bad edges (probably overlaps)
                edge_colors.append('white')
                edge_color_codes.append(-1)
                scaled_ewgs.append(None)
                continue
            edge_colors.append('black')
            edge_color_codes.append(-2)
            scaled_ewgs.append(None)
            continue
        if bound == 0:
            scaled_ewg = 0.
        else:
            scaled_ewg = (ewg+bound)/(2.*bound)
        scaled_ewgs.append(scaled_ewg)
        if ewg == -100: #lost edges:
            edge_colors.append('purple')
            edge_color_codes.append(-3)
        else:
            edge_colors.append(cmap(scaled_ewg))
            edge_color_codes.append(scaled_ewg)

    if len(node_score_list) > 1000:
        sorted_scores = sorted(node_score_list)
        perc_number = len(node_score_list) // 100
        nintynine_perc_value = sorted_scores[-perc_number]

        scaled_weighted_gi_scores = [gi_score/nintynine_perc_value for gi_score in node_score_list]
    elif len(node_score_list) > 0:
        max_score = max(node_score_list)
        if max_score > 0:
            scaled_weighted_gi_scores = [gi_score/max_score for gi_score in node_score_list]
        else:
            scaled_weighted_gi_scores = node_score_list
    else:
        scaled_weighted_gi_scores = []

    if len(edge_weights) > 1000:
        sorted_edge_weights = sorted(edge_weights)
        perc_number = len(edge_weights) // 100
        nintynine_perc_edge_weight_value = sorted_edge_weights[-perc_number]
    elif len(edge_weights) > 0:
        nintynine_perc_edge_weight_value = max(edge_weights)
    else:
        nintynine_perc_edge_weight_value = 0.

    scaled_edge_weights = []
    for pos, w in enumerate(edge_weights):
        if nintynine_perc_edge_weight_value > 0:
            scaled_wt_value = min([1.0, w/nintynine_perc_edge_weight_value])
            scaled_d_value = min([1.0, edge_weights_d_condition[pos]/nintynine_perc_edge_weight_value])
        else:
            scaled_wt_value = min([1.0, w])
            scaled_d_value = min([1.0, edge_weights_d_condition[pos]])
        if edge_colors[pos] == 'white': #improve if statement
            scaled_edge_weights.append(0)
            continue
        scaled_edge_weights.append(max([scaled_wt_value, scaled_d_value]))

    degree_distribution = {}
    for node_id in degree_map:
        deg = degree_map[node_id]
        if deg not in degree_distribution:
            degree_distribution[deg] = 0
        degree_distribution[deg] += 1

    print(f'Number of nodes: {len(node_id_list)}, Number of edges: {len(edge_weights)}')

    node_size = 10
    margin = 2*node_size
    graph_size = ((20*math.sqrt(len(node_id_list))*node_size) // 1) + 2*margin

    g = ig.Graph(n = len(node_ids), edges = edge_list)
    g.vs['name'] = node_id_list
    g.vs['weighted_gi_score'] = node_score_list
    g.vs['scaled_weighted_gi_score'] = scaled_weighted_gi_scores

    g.es['interaction_strength'] = edge_weights
    g.es['interaction_strength_d_condition'] = edge_weights_d_condition
    g.es['interaction_change'] = edge_weight_changes
    g.es['scaled_interaction_change'] = scaled_ewgs
    g.es['edge_color_codes'] = edge_color_codes


    layout = g.layout('drl')
    #layout = g.layout('circle')
    visual_style = {}
    visual_style['layout'] = layout
    visual_style['margin'] = margin
    visual_style['bbox'] = (graph_size, graph_size)

    visual_style['vertex_label'] = g.vs['name']
    visual_style['vertex_label_size'] = node_size // 2
    visual_style['vertex_size'] = node_size
    visual_style['vertex_color'] = [cmap(1.0 - score) for score in scaled_weighted_gi_scores]
    #g.vs['vertex_color'] = visual_style['vertex_color']

    visual_style['edge_width'] = [0.5*node_size*weight for weight in scaled_edge_weights]
    g.es['edge_width'] = visual_style['edge_width']
    visual_style['edge_color'] = edge_colors
    #g.es['edge_color'] = visual_style['edge_color']

    svg_file = f'{outfile_stem}.svg'

    ig.plot(g, svg_file, **visual_style)

    g.write_gml(gml_file)
    
    degree_distribution_file = f'{outfile_stem}_degree_distribution.tsv'
    f = open(degree_distribution_file, 'w')
    header = 'Degree\tCount\n'
    f.write(header)
    for deg in sorted(degree_distribution.keys()):
        line = f'{deg}\t{degree_distribution[deg]}\n'
        f.write(line)
    f.close()

    if len(degree_map) == 0:
        return
    
    try:
        pl_fit = powerlaw.Fit(numpy.array(list(degree_map.values())),xmin=1, xmax=1500)
        figure = pl_fit.plot_pdf(color = 'b')
        plt.xlabel('Node degree')
        plt.ylabel('PDF (blue), CCDF (red)')
        pl_fit.power_law.plot_pdf(color = 'b', linestyle = 'dashed', ax = figure)
        pl_fit.plot_ccdf(color = 'r', ax = figure)
        pl_fit.power_law.plot_ccdf(color = 'r', linestyle = 'dashed', ax = figure)    

        powerlaw_file = f'{outfile_stem}_powerlaw_plot.png'
        plt.savefig(powerlaw_file)
        plt.clf()
        plt.cla()
        plt.close()
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        print(f"Creating the power law plot failed due to:\n{e}\n{f}\n{g}")

def parse_target_list_file(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()

    targets = []
    for line in lines:
        target = line.replace('\n','')
        if target == '':
            continue
        targets.append(target)

    return targets

@ray.remote(max_calls = 1)
def gene_report_preprocess(chunk, store):
    (config, gene_isoform_score_dict, cmap, outfile_stem, out_f, session_name) = store
    whole_graph = None
    uninterresting_isoform_pairs_dict = {}
    plot_dict = {}
    errors = []
    for target in chunk:
        gene_class = gene_isoform_score_dict[target][5]
        if config.verbosity >= 3:
            print(f'Creating Gene Report for {gene_class.gene_id} - {gene_class.name}')

        path_to_isoform_relations_plot, longest_isoform, uninterresting_isoform_pairs, path_to_gene_expression_plot, usage_plots, insignificant_isoforms = create_isoform_relations_plot(config, target, gene_isoform_score_dict, cmap, outfile_stem, condition_1_label = config.condition_1_tag, condition_2_label = config.condition_2_tag)

        uninterresting_isoform_pairs_dict[target] = (insignificant_isoforms, uninterresting_isoform_pairs)

        if longest_isoform is None:
            errors.append(f'Creating isoform relations plot failed for {target}')
            continue
        
        input_for_sequence_plot = gene_class.generate_input_for_sequence_plot(longest_isoform, insignificant_isoforms)
        if config.verbosity >= 5:
            print(f'Input for sequence plot for {target}:\n{input_for_sequence_plot}')

        try:
            paths_to_isoform_plot = create_isoform_plot(input_for_sequence_plot, out_f, session_name, target)
        except:
            paths_to_isoform_plot = None
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            errors.append(f"Creating the isoform plot failed due to:\n{e}\n{f}\n{g}")

        plot_dict[target] = (path_to_isoform_relations_plot, paths_to_isoform_plot, path_to_gene_expression_plot, usage_plots)

        #whole_graph = retrieve_sub_graph(gml_file, config, out_f, session_name, target, depth = 1, filtered = False, whole_graph = whole_graph, return_whole_graph=True)
    return plot_dict, uninterresting_isoform_pairs_dict, errors
    

def create_gene_report(session_id, target, config, out_f, session_name):

    if config.verbosity >= 2:
        print(f'Calling create_gene_report: session_id {session_id}, target {target}, out folder {out_f}, session_name {session_name}')
    if os.path.isfile(resolve_path(target)):
        targets = parse_target_list_file(resolve_path(target))
    else:
        targets = [target]

    if target == 'all':
        targets = None

    #gene_gene_stem = os.path.join(out_f, f'{session_name}_gene_gene_network')
    #gml_file = f'{gene_gene_stem}.gml'
    #if not os.path.isfile(gml_file):
    #    generate_gene_interaction_network(session_id, config, out_f, session_name)

    gene_isoform_score_dict, targets, prot_db_ids, proteins = calculate_gene_isoform_score_dict(session_id, config, create_gene_class = True, target_genes = targets, return_prot_db_ids = True)

    outfile_stem = os.path.join(out_f, f'{session_name}')
    cmap=mpl.colormaps['RdYlGn']

    store = ray.put((config, gene_isoform_score_dict, cmap, outfile_stem, out_f, session_name))

    small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(config.proc_n, len(targets))

    chunks = []
    chunk = []
    for target in targets:
        chunk.append(target)
        if n_of_big_chunks > 0 and len(chunks) < n_of_big_chunks:
            if len(chunk) == big_chunksize:
                chunks.append(chunk)
                chunk = []
        elif n_of_small_chunks > 0:
            if len(chunk) == small_chunksize:
                chunks.append(chunk)
                chunk = []

    print(f'Amount of chunks: {len(chunks)}')

    gene_report_preprocess_ids = []
    for chunk in chunks:
        gene_report_preprocess_ids.append(gene_report_preprocess.remote(chunk, store))

    """
    whole_graph = None
    uninterresting_isoform_pairs_dict = {}
    plot_dict = {}
    for target in targets:
        gene_class = gene_isoform_score_dict[target][5]
        if config.verbosity >= 3:
         print(f'Creating Gene Report for {gene_class.gene_id} - {gene_class.name}')

        path_to_isoform_relations_plot, longest_isoform, uninterresting_isoform_pairs, path_to_gene_expression_plot, usage_plots, insignificant_isoforms = create_isoform_relations_plot(config, target, gene_isoform_score_dict, cmap, outfile_stem, condition_1_label = config.condition_1_tag, condition_2_label = config.condition_2_tag)

        uninterresting_isoform_pairs_dict[target] = (insignificant_isoforms, uninterresting_isoform_pairs)


        if longest_isoform is None:
            print(f'Creating isoform relations plot failed for {target}')
            continue

        
        input_for_sequence_plot = gene_class.generate_input_for_sequence_plot(longest_isoform)
        if config.verbosity >= 5:
            print(f'Input for sequence plot for {target}:\n{input_for_sequence_plot}')

        try:
            paths_to_isoform_plot = create_isoform_plot(input_for_sequence_plot, out_f, session_name, target)
        except:
            paths_to_isoform_plot = None
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            print(f"Creating the isoform plot failed due to:\n{e}\n{f}\n{g}")

        plot_dict[target] = (path_to_isoform_relations_plot, paths_to_isoform_plot, path_to_gene_expression_plot, usage_plots)

        #whole_graph = retrieve_sub_graph(gml_file, config, out_f, session_name, target, depth = 1, filtered = False, whole_graph = whole_graph, return_whole_graph=True)
    
    """

    uninterresting_isoform_pairs_dict = {}
    plot_dict = {}
    for chunk_result in ray.get(gene_report_preprocess_ids):
        sub_plot_dict, sub_uninterresting_isoform_pairs_dict, errors = chunk_result
        for error in errors:
            print(error)
        plot_dict.update(sub_plot_dict)
        uninterresting_isoform_pairs_dict.update(sub_uninterresting_isoform_pairs_dict)

    blacklist = set()
    for target in uninterresting_isoform_pairs_dict:
        (insignificant_isoforms, uninterresting_isoform_pairs) = uninterresting_isoform_pairs_dict[target]
        blacklist = blacklist | insignificant_isoforms

    if config.verbosity >= 3:
        print(f'Blacklist: {blacklist}')

    massmodel.mass_model(session_id, config, out_f, with_multi_mutations = True, prot_db_ids = prot_db_ids, skip_individual_indel_mutants = True, blacklist = blacklist, update_proteins = proteins.protein_map)
    md_lines = [f"# **--- Gene Report for {session_name} ---**\n\n\n"]
    for target in targets:
        gene_class = gene_isoform_score_dict[target][5]
        model_plot_dicts = {}
        for isoform in gene_class.isoforms:
           model_plot_dicts[isoform] = create_pymol_plot(config, out_f, isoform, uninterresting_isoform_pairs_dict[target][1])
        if target in plot_dict:
            md_lines += create_markdown(config, out_f, model_plot_dicts, gene_class.name, target, plot_dict[target])
        else:
            print(f'\n\nERROR: Creating Gene Report FAILED for {target} !!!\n\n')
    
    fuse_markdown(md_lines, session_name, out_f)

def retrieve_sub_graph(ggin_file, config, out_f, session_name, target, depth = 2, filtered = False, whole_graph = None, return_whole_graph = False):
    if config.verbosity >= 2:
        print(f'Calling retrieve_sub_graph for {target} and {ggin_file}, whole_graph: {whole_graph}')

    if whole_graph is None:
        try:
            g = ig.Graph.Read_GML(ggin_file)
        except:
            [e, f, g] = sys.exc_info()
            print(f"Couldn't load graph file in retrieve_sub_graph: {ggin_file}\n{e}")
            return None
    else:
        g = whole_graph

    try:
        target_node = g.vs.find(name = target)
    except:
        [e, f, g] = sys.exc_info()
        print(f"Couldn't locate target: {target} in graph: {ggin_file}\n{e}")
        return g
    
    new_nodes = [int(target_node['id'])]

    print(target_node.attributes())
    print(target_node['id'])

    nodes_to_process = [target_node]
    seen_nodes = set()
    current_depth = 0
    while current_depth < depth:
        new_nodes_to_process = []
        for node in nodes_to_process:
            node_id = int(node['id'])
            seen_nodes.add(node_id)
            for neigh in node.neighbors():
                neigh_id = int(neigh['id'])
                if neigh_id in seen_nodes:
                    continue
                if filtered:
                    try:
                        edge_id = g.get_eid(node_id, neigh_id)
                        edge = g.es[edge_id]
                        #edge = g.es.find(_between = ((node_id,), (neigh_id,)))
                    except:
                        edge_id = g.get_eid(neigh_id, node_id)
                        edge = g.es[edge_id]
                        #edge = g.es.find(_between = ((neigh_id,), (node_id,)))
                    if edge['scaledinteractionchange'] == 0.5:
                        continue
                new_nodes_to_process.append(neigh)
                new_nodes.append(neigh_id)
        nodes_to_process = new_nodes_to_process
        current_depth += 1

    g_sub = g.subgraph(new_nodes)

    max_edge_weight_change = 0.
    min_edge_weight_change = 0.
    for pos,scaled_change in enumerate(g_sub.es['interactionchange']):
        if scaled_change == float('inf'):
            continue
        if g_sub.es['edgecolorcodes'][pos] < 0:
            continue
        if scaled_change > max_edge_weight_change:
            max_edge_weight_change = scaled_change
        if scaled_change < min_edge_weight_change:
            min_edge_weight_change = scaled_change

    bound = max([max_edge_weight_change, abs(min_edge_weight_change)])

    rescaled_changes = []

    for scaled_change in g_sub.es['interactionchange']:
        if bound == 0:
            rescaled_change = 0.
        else:
            rescaled_change = (scaled_change+bound)/(2.*bound)

        rescaled_changes.append(rescaled_change)

    g_sub.es['scaledinteractionchange'] = rescaled_changes

    cmap=mpl.colormaps['RdYlGn']

    edge_colors = []
    for pos, edge_color_code in enumerate(g_sub.es['edgecolorcodes']):
        if edge_color_code == -1:
            edge_colors.append('white')
        elif edge_color_code == -2:
            edge_colors.append('black')
        elif edge_color_code == -3:
            edge_colors.append('purple')
        else:
            edge_colors.append(cmap(rescaled_changes[pos]))

    filtered_str = ''
    if filtered:
        filtered_str = '_filtered'

    svg_file = f'{out_f}/{target}_D{depth}_subgraph{filtered_str}.svg'

    print(svg_file)

    node_size = 20
    margin = 2*node_size
    graph_size = ((4*math.sqrt(len(new_nodes))*node_size) // 1) + 2*margin
    layout = g_sub.layout('drl')

    visual_style = {}
    visual_style['layout'] = layout
    visual_style['margin'] = margin
    visual_style['bbox'] = (graph_size, graph_size)
    visual_style['vertex_label'] = g_sub.vs['name']
    visual_style['vertex_label_size'] = node_size // 2
    visual_style['vertex_size'] = node_size
    visual_style['vertex_color'] = [cmap(1.0 - score) for score in g_sub.vs["scaledweightedgiscore"]]
    visual_style['edge_width'] = g_sub.es["edgewidth"]
    visual_style['edge_color'] = edge_colors

    gml_file = f'{out_f}/{target}_D{depth}_subgraph{filtered_str}.gml'
    g_sub.write_gml(gml_file)
    if return_whole_graph:
        return g


def create_pymol_plot(config, out_f, target, uninterresting_isoform_pairs, dpi=200):
    try:
        import pymol
        import pymol.cmd as pymol_cmd
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        print(f"importing pymol failed due to:\n{e}\n{f}\n{g}")
        return None
    
    #parse models from summary
    model_summary = pd.read_csv(f"{out_f}/model_summary.tsv", sep='\t')

    pymol_out_f = f"{out_f}/models_pymol_plots"

    if config.verbosity >= 2:
        print(f'Calling of create_pymol_plot: target protein: {target}, outfolder: {pymol_out_f}')

    list_of_paths = model_summary['File location']
    prot_ids = model_summary['Input ID']
    mut_prot_ids = model_summary['Additional Label']
    plot_dict = {}
    
    if not os.path.exists(pymol_out_f):
        os.makedirs(pymol_out_f)

    #produce one png with all structures
    #pymol_cmd.set("grid_mode", 1)

    for pos, prot_id in enumerate(prot_ids):
        if prot_id != target:
            continue
        file = list_of_paths[pos]
        mut_prot_id = mut_prot_ids[pos]
        if (prot_id, mut_prot_id) in uninterresting_isoform_pairs:
            continue
        pymol_cmd.delete('all')
        pymol_cmd.load(file)

        # indel plot
        # colors given here are not the final colors for the plot but similar
        # substitution, interacting chain, insertion, deletion, no chamge
        # if selection is set to all, change also in 1393, colors all interacting molecules by assigned b factor
        #pymol_cmd.spectrum("b", "blue_white_red", "n. ca", minimum=0.00, maximum=100.00)
        #pymol_spectrumany.spectrumany("b", "blue white cyan brightorange red", "n. ca", minimum=0.00, maximum=100.00)
        pymol_spectrumany.spectrumany("b", "black white cyan blue red", "n. ca", minimum=0.00, maximum=100.00)
        
        pymol_cmd.orient()

        outfile = f"{pymol_out_f}/{target}_{mut_prot_id}_pymol_indel.png"
        pymol_cmd.png(outfile, width=1200, dpi=dpi, ray=1)

        #produce legend for indel plot
        pymol.color_list = []
        rgb_list = []
        pymol_cmd.iterate('n. ca', 'pymol.color_list.append(color)')
        for color_id in pymol.color_list:
            rgb_list.append(pymol_cmd.get_color_tuple(color_id))

        # needs adaption if more colors are added or changed
        rgb_dict = {}
        for color in rgb_list:
            if color in rgb_dict:
                continue
            else:
                if color[0] > 0.9 and color[1] > 0.9 and color[2] > 0.9:
                    rgb_dict['white'] = color
                elif color[0] > 0.9 and color[1] < 0.2 and color[2] < 0.2:
                    rgb_dict['red'] = color
                elif color[0] < 0.1 and color[1] < 0.1 and color[2] > 0.9: 
                    rgb_dict['blue'] = color
                elif color[0] < 0.1 and color[1] < 0.1 and color[2] < 0.1:
                    rgb_dict['black'] = color
                elif color[0] < 0.1 and color[1] > 0.9 and color[2] > 0.9:
                    rgb_dict['cyan'] = color
                
                """
                if color[0] > 0.9 and color[1] > 0.9 and color[2] > 0.9:
                    rgb_dict['white'] = color
                elif color[0] > 0.9 and color[2] < 0.2:
                    rgb_dict['red'] = color
                elif color[2] == 1.0 and color[1] == 0: 
                    rgb_dict['blue'] = color
                elif color[0] > 0.9 and color[2] > 0:
                    rgb_dict['orange'] = color
                else:
                    rgb_dict['cyan'] = color
                """

        image = plt.imread(outfile)
        fig = plt.figure()
        ax = plt.subplot()
        patches = []

        if 'red' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['red'], label='no change'))
        if 'cyan' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['cyan'], label='insertion'))
        if 'black' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['black'], label='substitution'))
        if 'blue' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['blue'], label='deletion'))
        if 'white' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['white'], label='interacting chain'))
        
        """
        if 'red' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['red'], label='no change'))
        if 'cyan' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['cyan'], label='insertion'))
        if 'blue' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['blue'], label='substitution'))
        if 'orange' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['orange'], label='deletion'))
        if 'white' in rgb_dict.keys():
            patches.append(mpatches.Patch(color=rgb_dict['white'], label='interacting chain'))
        """

        ax.legend(handles=patches, prop={'size': 4}, loc='upper left')
        ax.set_axis_off()
        ax.imshow(image)
        fig.savefig(outfile, dpi=dpi, bbox_inches='tight', pad_inches=0)
        plt.close()

        # plot colored by chain
        # pymol_cmd.util.cbc(selection="n. ca")
        # outfile = f"{pymol_out_f}/{target}_{mut_prot_id}_pymol_by_chain.png"
        # pymol_cmd.png(outfile, width=1200, dpi=300, ray=1)
        # pymol_cmd.delete('all')

        plot_dict[mut_prot_id] = outfile

    return plot_dict

def add_figure_to_md(md_lines, figure_title, path_to_plot):
    md_lines.append(f'![{figure_title}]({path_to_plot})\n')

def create_markdown(config, out_f, model_plot_dicts, gene_name, target, plots):
    if config.verbosity >= 2:
        print(f'Calling of create_markdown: target {target}')

    path_to_isoform_relations_plot, paths_to_isoform_plot, path_to_gene_expression_plot, usage_plots = plots
    md_lines = []
    md_lines.append(f"\n\n# Report for {target} - {gene_name}\n\n")
    if paths_to_isoform_plot is not None:
        md_lines.append(f"## Alignment plot: \n")
        for path_to_isoform_plot in paths_to_isoform_plot:
            add_figure_to_md(md_lines, 'isoforms', path_to_isoform_plot)
    
    md_lines.append('<div style="page-break-after: always;"></div>')

    md_lines.append(f'## Gene expression plot:\n')
    add_figure_to_md(md_lines, 'gene_expressions', path_to_gene_expression_plot)
    md_lines.append(f'## Isoform usages:\n')
    for isoform in usage_plots:
        add_figure_to_md(md_lines, f'isoform_usage_{isoform}', usage_plots[isoform])

    md_lines.append(f"## Isoform relation graph: \n")
    add_figure_to_md(md_lines, 'isoforms', path_to_isoform_relations_plot)
    
    for isoform in model_plot_dicts:
        md_lines.append(f'## Structure models for: {isoform}:\n')
        for mutant_isoform in model_plot_dicts[isoform]:
            md_lines.append(f'### Difference to {mutant_isoform}:\n')
            add_figure_to_md(md_lines, f'structure_{isoform}_{mutant_isoform}', model_plot_dicts[isoform][mutant_isoform])
    md_lines.append('<div style="page-break-after: always;"></div>')
    return md_lines

    

def fuse_markdown(md_lines, session_name, out_f, path_add = ''):
    path_to_md = f"{out_f}/{session_name}{path_add}_gene_report.md"
    path_to_html = f"{out_f}/{session_name}{path_add}_gene_report.html"
    path_to_pdf = f"{out_f}/{session_name}{path_add}_gene_report.pdf"

    f = open(path_to_md, 'w')
    f.write(''.join(md_lines))
    f.close()
        
    html_string = ""
    with open(path_to_md, 'r') as r:
        list = r.readlines()
        for obj in list:
            html_string += markdown.markdown(obj)

    with open(path_to_html, 'w') as w:
        w.write(html_string)

    pdfkit.from_file(path_to_html, path_to_pdf)

class MyPlotter(graphics.LetterPlotter):

    #returns color at trace[column_i, seq_i]
    def set_classifications(self, classification_list: list, classification_color_dict: dict, trace):
        self.class_list = classification_list
        self.class_color_dict = classification_color_dict
        self.trace = trace

    def get_color(self, alignment, column_i, seq_i):
        if seq_i == 0:
            return 'w'
        if self.trace[column_i][seq_i] == -1:
            #case gap, return white
            return 'w'
        color = self.class_color_dict[self.class_list[seq_i - 1][column_i]]
        return color

def align_by_mutation_list(major_sequence_raw, isoform_sequences, isoform_dict, major_isoform_id):
    raw_major_gaps = []
    transformed_isoform_dict = {}
    for isoform_id in isoform_dict:
        insertions = isoform_dict[isoform_id][0]
        deletions = isoform_dict[isoform_id][1]
        transformed_insertions = {}
        transformed_deletions = set(deletions)
        for insertion_pos, inserted_sequence in insertions:
            raw_major_gaps.append((insertion_pos, len(inserted_sequence)))
            transformed_insertions[insertion_pos] = len(inserted_sequence)
        transformed_isoform_dict[isoform_id] = (transformed_insertions, transformed_deletions)


    fused_major_gaps = {}
    for insertion_pos, gap_len in raw_major_gaps:
        if insertion_pos not in fused_major_gaps:
            fused_major_gaps[insertion_pos] = gap_len
        elif gap_len > fused_major_gaps[insertion_pos]:
            fused_major_gaps[insertion_pos] = gap_len
    
    aligned_major_sequence = ''

    if 0 in fused_major_gaps:
        aligned_major_sequence = '-'*fused_major_gaps[0]

    for seq_pos_zero, aa in enumerate(major_sequence_raw):
        seq_pos = seq_pos_zero + 1

        aligned_major_sequence += aa
        if seq_pos in fused_major_gaps:
            aligned_major_sequence += '-'*fused_major_gaps[seq_pos]

    aligned_isoforms = {major_isoform_id : aligned_major_sequence}
    for isoform_id in transformed_isoform_dict:
        insertions = transformed_isoform_dict[isoform_id][0]
        deletions = transformed_isoform_dict[isoform_id][1]
        raw_seq = isoform_sequences[isoform_id]
        aligned_seq = ''
        seq_pos = 0
        iso_pos = 0
        current_insert_len = 0

        for aa in aligned_major_sequence:
            if aa == '-':
                if seq_pos in insertions:
                    insert_len = insertions[seq_pos]
                    if insert_len > current_insert_len:
                        aligned_seq += raw_seq[iso_pos]
                        iso_pos += 1
                        current_insert_len += 1
                    else:
                        aligned_seq += '-'
                else:
                    aligned_seq += '-'
            else:
                current_insert_len = 0
                seq_pos += 1

                if seq_pos in deletions:
                    aligned_seq += '-'
                else:
                    aligned_seq += raw_seq[iso_pos]
                    iso_pos += 1
        aligned_isoforms[isoform_id] = aligned_seq

    return aligned_isoforms

def plot_alignment_chunk(alignment, classification_color_dict, order_classifications, labels, trace, outfile, outlegend, dpi=100, font_size = 12):
    handles = []

    for key in classification_color_dict:
        if key == 'gap':
            continue
        else:
            if key is None:
                patch = mpatches.Patch(color=classification_color_dict[key], label="None")
            else:    
                patch = mpatches.Patch(color=classification_color_dict[key], label=key)
            handles.append(patch)

    trace_array = alignment.trace
    plt.rcParams.update({'font.size':font_size})

    # 50 is an appropriate size of bases or positions per line in the plot
    # legend is automatically plotted at the end of every chunk
    # legend can additionally be accessed by the return value of plot_alignment_chunk if necessary
    if len(alignment) > 50:
        total_height = 0.8*math.ceil(len(alignment)/50)*(len(trace_array[0]))
        fig, ax = plt.subplots(figsize=(font_size*1.2, total_height))
        plotter = MyPlotter(ax)
        plotter.set_classifications(order_classifications, classification_color_dict, trace)

        graphics.plot_alignment(
            ax, alignment, plotter, symbols_per_line=50, labels=labels)
        legend = ax.legend(loc='upper left', handles=handles, bbox_to_anchor=(0, -0.05), ncol=math.ceil(len(handles)/len(trace_array[0])))
        plot_legend(legend, outlegend)
        #ax.legend().remove()
    else:   
        fig, ax = plt.subplots(figsize=(len(alignment), len(trace_array[0])*0.5))
        plotter = MyPlotter(ax)
        plotter.set_classifications(order_classifications, classification_color_dict, trace)
        graphics.plot_alignment(
            ax, alignment, plotter, symbols_per_line=len(alignment), labels=labels)
        legend = ax.legend(loc='upper left', handles=handles, bbox_to_anchor=(0, -0.05), ncol=math.ceil(len(handles)/len(trace_array[0])))
        plot_legend(legend, outlegend)
        #ax.legend().remove()

    plt.savefig(outfile, bbox_inches='tight', dpi = dpi, pad_inches=0)
    plt.clf()
    plt.cla()
    plt.close()
    return legend
    


def create_isoform_plot(obj, out_f, session_name, target):
    major_sequence_raw = obj[0]
    major_isoform_id = obj[1]
    protein_sequence = seq.ProteinSequence(major_sequence_raw)
    classifications = obj[2]
    isoforms = obj[3]
    isoform_sequences = obj[4]
    isoform_classifications = obj[5]

    aligned_sequences = align_by_mutation_list(major_sequence_raw, isoform_sequences, isoforms, major_isoform_id)
    isoform_classifications[major_isoform_id] = classifications

    all_classifications = {}
    #insert gaap classifications
   
    isoform_plot_folder = f'{out_f}/isoform_plots'
    if not os.path.isdir(isoform_plot_folder):
        os.mkdir(isoform_plot_folder)

    for isoform_id in aligned_sequences:
        all_classifications[isoform_id] = []
        i = 0
        for char in aligned_sequences[isoform_id]:

            if char == '-':
                all_classifications[isoform_id].append('gap')
            else:
                original_classification = isoform_classifications[isoform_id][i]
                if original_classification == "None":
                    all_classifications[isoform_id].append("NA")

                else:
                    all_classifications[isoform_id].append(original_classification)
                i += 1

    classification_color_dict = parse_classifications(all_classifications)
    
    labels = ['position', major_isoform_id]
    for key in isoforms.keys():
        labels.append(key)

    ordered_isoforms = [major_isoform_id]
    for isoform_id in isoform_sequences:
        ordered_isoforms.append(isoform_id)

    order_classifications = []
    for isoform in ordered_isoforms:
        order_classifications.append(all_classifications[isoform])

    positions = list(range(1,len(aligned_sequences[major_isoform_id])+1))

    paths = []
    #before 500
    alignment_chunksize = 34*50//(len(labels)+2)

    aligned_sequences_chunks = []
    for isoform_id in ordered_isoforms:
        aligned_sequences_chunks.append(aligned_sequences[isoform_id][:alignment_chunksize])

    chunked_classifications = []
    for classifications in order_classifications:
        chunked_classifications.append(classifications[:alignment_chunksize])

    #chunked_trace = trace_final[:alignment_chunksize]
    raw_chunked_positions = positions[:alignment_chunksize]
    chunked_positions = []
    count = 1
    for pos in raw_chunked_positions:
        if pos < 101:
            chunked_positions.append(pos)
        elif count == 5:
            chunked_positions.append(pos)
            count = 1
        else:
            chunked_positions.append('')
            count += 1

    plotting_done = False
    i = 0
    dpi = 100
    while not plotting_done:
        ungapped_sequences = []
        for aligned_seq in aligned_sequences_chunks:  
            ungapped_sequences.append(aligned_seq.replace('-',''))

        protein_sequences = [chunked_positions]
        for sequence in ungapped_sequences:
            protein_sequence = seq.ProteinSequence(sequence)
            protein_sequences.append(protein_sequence)

        trace = produce_trace(aligned_sequences_chunks)

        j = 0
        trace_final = []

        for entry in trace:
            entry = [j] + entry[0:]
            trace_final.append(entry)
            j += 1

        trace_array = np.asarray(trace_final, dtype=np.int64)
    
        alignment = align.Alignment(protein_sequences, trace_array)

        outfile = f'{isoform_plot_folder}/{target}_isoform_plot_{i}.png'
        outlegend = f'{isoform_plot_folder}/{target}_isoform_plot_legend.png'
        plot_alignment_chunk(alignment, classification_color_dict, chunked_classifications, labels, trace_final, outfile, outlegend, dpi=dpi)
        paths.append(outfile)

        if len(positions) < (i+1)*alignment_chunksize:
            plotting_done = True
            #paths.append(outlegend)
            
        else:
            aligned_sequences_chunks = []
            for isoform_id in ordered_isoforms:
                aligned_sequences_chunks.append(aligned_sequences[isoform_id][(i+1)*alignment_chunksize:(i+2)*alignment_chunksize])

            chunked_classifications = []
            for classifications in order_classifications:
                chunked_classifications.append(classifications[(i+1)*alignment_chunksize:(i+2)*alignment_chunksize])

            #chunked_trace = trace_final[(i+1)*alignment_chunksize:(i+2)*alignment_chunksize]
            raw_chunked_positions = positions[(i+1)*alignment_chunksize:(i+2)*alignment_chunksize]
            chunked_positions = []
            
            count = 1
            for pos in raw_chunked_positions:
                if pos < 101:
                    chunked_positions.append(pos)
                elif count == 5:
                    chunked_positions.append(pos)
                    count = 1
                else:
                    chunked_positions.append('')
                    count += 1
            i += 1
    return paths

def plot_legend(legend, outfile):
    fig = legend.figure
    fig.canvas.draw()
    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(outfile, bbox_inches=bbox)

def divide_in_chunks(alignment, n):
    for i in range(0, len(alignment), n):
        yield alignment[i:i + n]

def produce_trace(aligned_sequences):
    n_rows = 0
    n_cols = len(aligned_sequences)
    for entry in aligned_sequences:
        if len(entry) > n_rows:
            n_rows = len(entry)

    #initialize trace with -2
    trace = [[-2 for _ in range(n_cols)] for _ in range(n_rows)]

    column_i = 0
    seq_j = 0
    for seq_j, sequence in enumerate(aligned_sequences):
        id_counter = 0
        for column_i, symbol in enumerate(sequence):
            if symbol == "-":
                trace_elem = -1
            else:
                trace_elem = id_counter
                id_counter+=1
            trace[column_i][seq_j] = trace_elem
    return trace
    

    
def parse_classifications(data):
    classification_color_dict = {}
    colors = list(mcolors.TABLEAU_COLORS.values())
    
    i = 0
    for isoform_id in data:
        classification_sequence = data[isoform_id]

        classification_dict = {}
        for classification in classification_sequence:
            if classification in classification_color_dict:
                continue
            if classification == 'gap':
                continue
 
            classification_color_dict[classification] = colors[i]
            i += 1
            if i == 10:
                colors = list(mcolors.XKCD_COLORS.values())
    return classification_color_dict