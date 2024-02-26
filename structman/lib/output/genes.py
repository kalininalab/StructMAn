import os
import sys
import traceback
import math
import igraph as ig
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


import biotite.sequence as seq
import biotite.sequence.align as align
import matplotlib.pyplot as plt
import biotite.sequence.graphics as graphics
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import math
import colorsys

import powerlaw
import numpy

import markdown


from structman.lib.database import database
from structman.lib.database.database_core_functions import select
from structman.lib.lib_utils import fuse_multi_mutations
from structman.lib.output import out_generator, out_utils, massmodel
from structman.base_utils.base_utils import resolve_path

class Gene:

    def __init__(self, name, isoform_pair_multimutations, isoform_pair_scores, isoform_expressions, isoform_prot_objs):
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

    def generate_input_for_sequence_plot(self, isoform_id):
        prot = self.isoforms[isoform_id][2]
        seq = prot.sequence
        classifications = prot.get_ordered_classifications()
        isoform_differences = {}
        for iso_a, iso_b in self.isoform_pairs:
            if iso_a != isoform_id:
                continue
            insertions, dels, subs = self.isoform_pairs[(iso_a, iso_b)][1].get_affected_positions()
            #print(iso_a, iso_b)
            #self.isoform_pairs[(iso_a, iso_b)][1].print_indels()
            isoform_differences[iso_b] = insertions, dels, subs
        return seq, isoform_id, classifications, isoform_differences

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

        unweighted_score = multi_mutation.get_score(proteins)
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

    fused_multi_mutation = fuse_multi_mutations(wt_reduced_multi_mutations)

    #print(f'Fused Isoform:')
    #fused_multi_mutation.print_indels()

    try:
        GI_score = fused_multi_mutation.get_score(proteins)
    except:
        GI_score = None

    return GI_score, wt.primary_protein_id, weighted_gi_score, max_contribution_isoform, normalized_max_contribution, isoform_pair_scores, isoform_expressions

def init_gene_maps(session_id, config):
    table = 'RS_Isoform'
    rows = ['Protein_A', 'Protein_B', 'Gene', 'Multi_Mutation']
    eq_rows = {'Session': session_id}

    results = select(config, rows, table, equals_rows=eq_rows)

    gene_db_ids = set()
    gene_isoform_map = {}
    for row in results:
        prot_a_db_id = row[0]
        prot_b_db_id = row[1]
        gene_db_id = row[2]
        multi_mutation_db_id = row[3]

        gene_db_ids.add(gene_db_id)
        if not gene_db_id in gene_isoform_map:
            gene_isoform_map[gene_db_id] = []
        gene_isoform_map[gene_db_id].append((prot_a_db_id, prot_b_db_id, multi_mutation_db_id))
    
    gene_id_map = database.retrieve_gene_id_map(gene_isoform_map.keys(), config)

    return gene_isoform_map, gene_id_map


def calculate_gene_isoform_score_dict(session_id, config, proteins = None, create_gene_class = False, target_genes = None, return_prot_db_ids = False):
    gene_isoform_map, gene_id_map = init_gene_maps(session_id, config)

    if proteins is None:
        proteins = database.proteinsFromDb(session_id, config, with_multi_mutations = True, with_mappings = True)

    gene_isoform_score_dict = {}

    prot_db_ids = set()

    for gene_db_id in gene_isoform_map:
        gene_name = gene_id_map[gene_db_id]

        if target_genes is not None:
            if gene_name not in target_genes:
                continue

        isoforms = []
        isoform_prot_objs = {}
        for (prot_a_db_id, prot_b_db_id, multi_mutation_db_id) in gene_isoform_map[gene_db_id]:
            prot_a = proteins.getByDbId(prot_a_db_id)
            prot_b = proteins.getByDbId(prot_b_db_id)
            multi_mutation = proteins.multi_mutation_back_map[multi_mutation_db_id]
            isoforms.append((prot_a, prot_b, multi_mutation))
            if prot_a.primary_protein_id not in isoform_prot_objs:
                isoform_prot_objs[prot_a.primary_protein_id] = prot_a
                if return_prot_db_ids:
                    prot_db_ids.add(prot_a_db_id)
            

        if config.verbosity >= 3:
            print(f'Call of calculateGeneIsoformScore for {gene_name}\n')

        gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score, isoform_pair_scores, isoform_expressions = calculateGeneIsoformScore(isoforms, proteins, wt_condition_tag = config.condition_1_tag, disease_condition_tag = config.condition_2_tag)

        if config.verbosity >= 3:
            print(f'GI Score for {gene_name}: {gi_score}')

        if create_gene_class:
            gene_class = Gene(gene_name, isoforms, isoform_pair_scores, isoform_expressions, isoform_prot_objs)
            gene_isoform_score_dict[gene_name] = gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score, gene_class
        else:
            gene_isoform_score_dict[gene_name] = gi_score, wt, weighted_gi_score, max_contribution_isoform, contribution_score

    if return_prot_db_ids:
        return gene_isoform_score_dict, prot_db_ids

    return gene_isoform_score_dict

def generate_gene_isoform_scores(session_id, config, outfolder, session_name):
    outfile = out_utils.generate_gis_filename(outfolder, session_name)

    gene_isoform_score_dict = calculate_gene_isoform_score_dict(session_id, config)

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

def create_isoform_relations_plot(gene_name, gene_isoform_score_dict, cmap, outfile_stem, condition_1_label = 'ExC1', condition_2_label = 'ExC2'):
    print(f'Creating isoform relations plot for: {gene_name}')

    gene_class = gene_isoform_score_dict[gene_name][5]

    if gene_class.all_expression_values_none:
        print(f'Creating isoform relations invalid for {gene_class.name}, expression values are not defined. Please set the condition tags with --condition_1 "tag name" and --condition_2 "tag name"')
        return None

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
        if expression_d_condition == 0 and expression_wt_condition == 0:
            continue
        total_expression_condition_1 += expression_wt_condition
        total_expression_condition_2 += expression_d_condition

    for isoform in gene_class.isoforms:
        expression_wt_condition, expression_d_condition, prot_object = gene_class.isoforms[isoform]
        if expression_d_condition == 0 and expression_wt_condition == 0:
            continue

        ex_change = expression_wt_condition - expression_d_condition
        if ex_change > max_ex_change:
            max_ex_change = ex_change
        if ex_change < min_ex_change:
            min_ex_change = ex_change

        raw_ex_changes[isoform] = ex_change

        if expression_wt_condition > 0:
            log_fold_change = round(math.log2(expression_d_condition/expression_wt_condition), 2)
        else:
            log_fold_change = 'inf'
        usage_condition_1 = round(expression_wt_condition/total_expression_condition_1, 2)
        usage_condition_2 = round(expression_d_condition/total_expression_condition_2, 2)
        usage_change = round(usage_condition_2 - usage_condition_1, 2)

        node_id_list.append(isoform)
        node_ids[isoform] = len(node_ids)
        seq_len = len(prot_object.sequence)
        if seq_len > max_seq_len:
            max_seq_len = seq_len
            longest_isoform = isoform
        label = f'{isoform}, Size: {seq_len}\n{condition_1_label}: {round(expression_wt_condition, 5)}\n{condition_2_label}: {round(expression_d_condition, 5)}\nl2fc: {log_fold_change}\nUsage {condition_1_label}: {usage_condition_1}\n Usage {condition_2_label}: {usage_condition_2}\nChange in usage: {usage_change}'
        #print(label)
        node_labels.append(label)

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
    for (iso_a, iso_b) in gene_class.isoform_pairs:
        #print(iso_a, iso_b)
        if iso_a not in node_ids:
            continue
        if iso_b not in node_ids:
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
            continue

        #if raw_ex_changes[iso_a] < 0 or raw_ex_changes[iso_b] > 0:
        #    continue

        if SDS == 0 and len(node_labels) > 6:
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
    
    svg_file = f'{outfile_stem}_{gene_name}_isoform_relations.svg'

    ig.plot(g, svg_file, **visual_style)

    return longest_isoform


def generate_isoform_relations(session_id, config, outfolder, session_name):
    outfile_stem = os.path.join(outfolder, f'{session_name}')

    gene_isoform_score_dict = calculate_gene_isoform_score_dict(session_id, config, create_gene_class = True)

    cmap=mpl.colormaps['RdYlGn']

    for gene_name in gene_isoform_score_dict:
        create_isoform_relations_plot(gene_name, gene_isoform_score_dict, cmap, outfile_stem)
       


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
        gene_isoform_score_dict = calculate_gene_isoform_score_dict(session_id, config, proteins = proteins)

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
    
    pl_fit = powerlaw.Fit(numpy.array(list(degree_map.values())),xmin=1, xmax=1500)
    figure = pl_fit.plot_pdf(color = 'b')
    plt.xlabel('Node degree')
    plt.ylabel('PDF (blue), CCDF (red)')
    pl_fit.power_law.plot_pdf(color = 'b', linestyle = 'dashed', ax = figure)
    pl_fit.plot_ccdf(color = 'r', ax = figure)
    pl_fit.power_law.plot_ccdf(color = 'r', linestyle = 'dashed', ax = figure)    

    powerlaw_file = f'{outfile_stem}_powerlaw_plot.png'
    plt.savefig(powerlaw_file)

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

def create_gene_report(session_id, target, config, out_f, session_name):
    if os.path.isfile(resolve_path(target)):
        targets = parse_target_list_file(resolve_path(target))
    else:
        targets = [target]

    gene_gene_stem = os.path.join(out_f, f'{session_name}_gene_gene_network')
    gml_file = f'{gene_gene_stem}.gml'
    if not os.path.isfile(gml_file):
        generate_gene_interaction_network(session_id, config, out_f, session_name)

    gene_isoform_score_dict, prot_db_ids = calculate_gene_isoform_score_dict(session_id, config, create_gene_class = True, target_genes = targets, return_prot_db_ids = True)

    outfile_stem = os.path.join(out_f, f'{session_name}')
    cmap=mpl.colormaps['RdYlGn']

    whole_graph = None
    for target in targets:
        longest_isoform = create_isoform_relations_plot(target, gene_isoform_score_dict, cmap, outfile_stem, condition_1_label = config.condition_1_tag, condition_2_label = config.condition_2_tag)

        gene_class = gene_isoform_score_dict[target][5]
        input_for_sequence_plot = gene_class.generate_input_for_sequence_plot(longest_isoform)
        if config.verbosity >= 4:
            print(f'Input for sequence plot for {target}:\n{input_for_sequence_plot}')

        path_to_isoform_plot = create_isoform_plot(input_for_sequence_plot, out_f, session_name, target)
        whole_graph = retrieve_sub_graph(gml_file, config, out_f, session_name, target, depth = 1, filtered = False, whole_graph = whole_graph, return_whole_graph=True)
        
    massmodel.mass_model(session_id, config, out_f, with_multi_mutations = True, prot_db_ids = prot_db_ids)
    paths_to_model_plots = create_pymol_plot(out_f)
    create_markdown(out_f, paths_to_model_plots)

def retrieve_sub_graph(ggin_file, config, out_f, session_name, target, depth = 2, filtered = False, whole_graph = None, return_whole_graph = False):
    if config.verbosity >= 2:
        print(f'Calling retrieve_sub_graph for {target} and {ggin_file}')

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


def create_pymol_plot(out_f):
    try:
        import pymol.cmd as pymol_cmd
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        print(f"importing pymol failed due to:\n{e}\n{f}\n{g}")
        return None
    #parse models from summary
    model_summary = pd.read_csv(f"{out_f}/model_summary.tsv", sep='\t')
    list_of_paths = model_summary['File location']
    paths_to_plots = []
    for path_to_pdb in list_of_paths:
        try:
            name_list = path_to_pdb.split("/models/",1)
            name = name_list[1]
            paths_to_plots.append(f'{out_f}/models-plots/{name}.png')
            pymol_cmd.load(path_to_pdb)
            pymol_cmd.spectrum("b", "rainbow", "n. and CA")
            pymol_cmd.orient()
            pymol_cmd.png(f'{out_f}/models-plots/{name}.png', dpi=150, ray=1)
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            print(f"couldn't create pymol plot for {path_to_pdb}:\n{e}\n{f}\n{g}")
    return paths_to_plots

def create_markdown(out_f, paths_to_model_plots):
    path_to_md = f"{out_f}/gene_report.md"
    path_to_html = f"{out_f}/gene_report.html"
    with open(path_to_md, 'w') as f:
        f.write("# Gene Report")
        #ToDo

    html_string = ""
    with open(path_to_md, 'r') as r:
        list = r.readlines()
        for obj in list:
            html_string = html_string = markdown.markdown(obj)

    with open(path_to_html, 'w') as w:
        w.write(html_string)

class MyPlotter(graphics.LetterPlotter):
    class_list = []
    class_color_dict = dict([])
    trace = []
    #returns color at trace[column_i, seq_i]
    def set_classifications(self, classification_list: list, classification_color_dict: dict, trace: list):
        self.class_list = classification_list
        self.class_color_dict = classification_color_dict
        self.trace = trace
    def get_color(self, alignment, column_i, seq_i):
        color = self.class_color_dict[self.class_list[column_i]]
        if seq_i == 0:
            return 'w'
        if self.trace[column_i][seq_i]==-1:
            amount = 0.5
            c = colorsys.rgb_to_hls(*mcolors.to_rgb(color))
            return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
        return color
    
def create_isoform_plot(obj, out_f, session_name, target):
    sequence = obj[0]
    id = obj[1]
    protein_sequence = seq.ProteinSequence(sequence)
    classifications = obj[2]
    isoforms = obj[3]
    
    isoform_triples = isoforms_as_triples(sequence, isoforms)
    classification_color_dict = parse_classifications(classifications)
    
    labels = ['position', id]
    for key in isoforms.keys():
        labels.append(key)
        
       # produce isoform strings, transform them to isoforms as ProteinSequences and produce trace
    isoform_sequences, trace, classification_list = produce_isoform_sequences(sequence, isoform_triples, classifications)
          
    i = 0
    trace_final = []
    positions = []
    for entry in trace:
        positions.append(i+1)
        entry = [i] + entry[0:]
        trace_final.append(entry)
        i += 1
        
    protein_seq_isoforms = [positions, protein_sequence]
    for iso in isoform_sequences:
        p_seq_iso = seq.ProteinSequence(iso)
        protein_seq_isoforms.append(p_seq_iso)
    
    trace_array = np.asarray(trace_final, dtype=np.int64)
    
    alignment = align.Alignment(protein_seq_isoforms, trace_array)


    handles = []
    for key in classification_color_dict:
        patch = mpatches.Patch(color=classification_color_dict[key], label=key)
        handles.append(patch)

    
    if len(alignment) > 50:
        total_height = math.ceil(len(alignment)/50)*(len(trace_array[0]))
        row_height = len(trace_array[0])+2
        plt.rcParams.update({'font.size':20+row_height})
        fig, ax = plt.subplots(figsize=(50+int(total_height/10), total_height))
        plotter = MyPlotter(ax)
        plotter.set_classifications(classification_list, classification_color_dict, trace_final)
        
        graphics.plot_alignment(
            ax, alignment, plotter, symbols_per_line=50, labels=labels)
        ax.legend(loc='upper left', handles=handles, bbox_to_anchor=(0, -0.05), ncol=math.ceil(len(handles)/len(trace_array[0])))
    else:   
         fig, ax = plt.subplots(figsize=(len(alignment), len(trace_array[0])*0.5))
         plotter = MyPlotter(ax)
         plotter.set_classifications(classification_list, classification_color_dict, trace_final)
         graphics.plot_alignment(
            ax, alignment, plotter, symbols_per_line=len(alignment), labels=labels)
         ax.legend(loc='upper left', handles=handles, bbox_to_anchor=(0, -0.05))

    path = f'{out_f}/{target}_isoform_plot.png'
    plt.savefig(path,bbox_inches='tight')
    return path
    
def parse_classifications(data: list) -> list:
    classification_dict = dict([])
    colors = list(mcolors.TABLEAU_COLORS.values())
    i = 0
    for elem in data:
        if i == 10:
            colors = list(mcolors.XKCD_COLORS.values())
        if elem in classification_dict:
            continue
        else:
            classification_dict[elem] = colors[i]
            i += 1
    return classification_dict
    
def isoforms_as_triples(sequence: str, isoforms: dict) -> list:
    triple_list = []
    for key in isoforms.keys():
        iso = isoforms.get(key)
        triple_list.append(iso)
        
    return triple_list


def produce_isoform_sequences(sequence: str, triples: list, classifications: list) -> list:
    isoform_list = []
    trace = []
    inserted_isoforms = []
    i = 0
    for char in sequence:
        trace.append([i])
        # append inserted_isoforms, where the length of inserted isoforms is stored at the index which is the same as the trace index
        inserted_isoforms.append(0)
        i += 1
    isoform_count_i = 1
    
    for isoform in triples:
        insertions = isoform[0]
        deletions = isoform[1]
        substitutions = isoform[2]
        inserted_deletions = []
        
        # insert indices as if isoforms are equal to sequence, insert -1 if entry in trace was created by an insertion
        i = 0
        index  = 0
        j = 0
        for entry in trace:
            if j > 0:
                entry.append(-1)
                j -= 1
            elif (inserted_isoforms[index] != 0):
                j = inserted_isoforms[index]-1
                entry.append(-1)
            else:
                entry.append(i)
                i +=1
            index += 1
            inserted_deletions.append(0)

        seq = sequence
        # substitutions
        # substitutions don't change the trace, they just change the seq
        for entry in substitutions:
            index = entry[0]-1
            new_char = entry[1]
            seq = seq[:index] + new_char + seq[index+1:]

        # deletions
        # deletions insert a -1 in the trace and remove the letter in seq
        # n counts deletions in current isoform
        n = 0
        for entry in deletions: 
            index = entry-1
            # count the isoforms already inserted in the alignment to adapt the trace index
            sum = 0
            i = 0
            while i <= index:
                sum += inserted_isoforms[i]
                i += 1
            # calculate index for trace
            trace_index = index + sum
            if inserted_isoforms[trace_index] != 0:
                trace_index += inserted_isoforms[trace_index]
            # sequence index depends on the deletions in current isoform
            index_seq = index - n

            seq = seq[:index_seq] + seq[index_seq+1:]

            # insert a -1 at position which has to be deleted
            trace[trace_index][isoform_count_i] = -2
            # decrease index in all following entries
            i = 1
            while i < (len(trace)-trace_index):
                if trace[len(trace)-i][isoform_count_i] == -2:
                    i+=1
                else:
                    trace[len(trace)-i][isoform_count_i] -=1
                    i +=1
                    
            inserted_deletions[trace_index] = 1
            n += 1
            
        # insertions
        # insertions add a new entry in the trace, they add chars to the seq
        m = 0
        for entry in insertions:
            # calculate index and modify sequence
            # m tracks insertion lengths in the current isoform
            # sum checks how often stuff has been inserted before index in other isoforms and how long thse insertions are
            index = entry[0]-1
            
            sum_i = 0
            sum_d = 0
            i = 0
            while i < index:
                sum_i += inserted_isoforms[i]
                sum_d += inserted_deletions[i]
                i += 1
            # calculate index for trace
            trace_index = index + sum_i 
            
            if len(inserted_deletions) > trace_index:
                i = index
                while i < trace_index:
                    sum_d += inserted_deletions[i]
                    i += 1
            # calculate index for sequence
            seq_index = index - sum_d + m

            new_char = entry[1]
            seq = seq[:seq_index] + new_char + seq[seq_index:]
            length = len(entry[1])
            m += length
            
            
            # insert an entry in trace filled with -1, insert an entry for every char
            n = length - inserted_isoforms[trace_index] 
            while n > 0:
                insert = []
                for entry in trace[0]:
                    insert.append(-1)
                # insert -1 entry in trace
                left = trace[:trace_index]
                right = trace[trace_index:]
                trace = left + [insert] + right
                
                insert_c = classifications[trace_index]
                left_c = classifications[:trace_index]
                right_c = classifications[trace_index:]
                classifications = left_c + [insert_c] + right_c
                
                left_d = inserted_deletions[:trace_index]
                right_d = inserted_deletions[trace_index:]
                deletions = left_d + [0] + right_d
                n -= 1
            # insert the length of insertion in the insertion tracking list
            if inserted_isoforms[trace_index] != 0:
                # an insertion already exists at this position
                # -> overwrite and insert (length-1) times a 0
                if length > inserted_isoforms[trace_index]:
                    inserted_isoforms[trace_index] = length
                
            else:
                # insertion does not exist yet
                # insert position and (length-1) times a 0
                left_i = inserted_isoforms[:trace_index]
                right_i = inserted_isoforms[trace_index:]
                inserted_isoforms = left_i + [length] + right_i
                
            i = 1
            while i < length:
                left_i = inserted_isoforms[:trace_index+i]
                right_i = inserted_isoforms[trace_index+i:]
                inserted_isoforms = left_i + [0] + right_i
                i += 1

            # insert an non -1 index at the insertion position
            if trace_index == 0:
                start = -1
            else:
                start = trace[trace_index-1][isoform_count_i]
                a = 1
                while start==-2 or start==-1:
                    start = trace[trace_index-1-a][isoform_count_i]
                    a+=1
            trace[trace_index][isoform_count_i] = start + 1
            # adapt following entrys
            j = length - 1
            i = 1
            while trace_index+i < len(trace):
                if j > 0:
                    trace[trace_index+i][isoform_count_i] = start + 2
                    i += 1
                    start += 1
                    j -= 1
                else:    
                    if trace[trace_index+i][isoform_count_i] == -2 or trace[trace_index+i][isoform_count_i] == -1:
                        i += 1
                    else:
                        trace[trace_index+i][isoform_count_i] = start + 2
                        i += 1
                        start += 1

        # finish isoform by adding the sequence to the isoform_list
        isoform_list.append(seq)
        isoform_count_i += 1
        
        # finish trace, replace deletions -2 with -1
        i = 0
        while i < len(trace):
            j = 0
            while j < len(trace[i]):
                if trace[i][j] == -2:
                    trace[i][j] = -1
                j += 1
            i += 1

    return isoform_list, trace, classifications
