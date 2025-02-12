import os

from structman.lib import rin
from structman.lib.database import database
from structman.lib.output import out_generator
from structman.base_utils import base_utils

def add_aggregate_results(aggregates, indel_output, aggregate_type, raw_aggregate_header_base_names):
    if aggregates is None:
        return

    aggregate, disorder_score, disorder_region = aggregates

    if aggregate is None:
        indel_output.add_value('%s IUPred value' % (aggregate_type), disorder_score)
        indel_output.add_value('%s Region structure type' % (aggregate_type), disorder_region)
        return

    (surface_value, mainchain_surface_value, sidechain_surface_value, profile_str, centrality_str, b_factor, modres, ssa, phi, psi, intra_ssbond, ssbond_length, intra_link,
     link_length, cis_conformation, cis_follower, inter_chain_median_kd, inter_chain_dist_weighted_kd, inter_chain_median_rsa, inter_chain_dist_weighted_rsa,
     intra_chain_median_kd, intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa, inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
     intra_chain_interactions_median, intra_chain_interactions_dist_weighted, rin_class, rin_simple_class) = aggregate

    centrality_scores = rin.Centrality_scores(code_str=centrality_str)
    try:
        rin_profile = rin.Interaction_profile(profile_str=profile_str)
    except:
        rin_profile = rin.Interaction_profile()
        print('Warning: problem happened rin profile str decode')

    aggregate_values = [surface_value, mainchain_surface_value, sidechain_surface_value, b_factor, modres, ssa, phi, psi, intra_ssbond, intra_link,
     cis_conformation, cis_follower, inter_chain_median_kd, inter_chain_dist_weighted_kd, inter_chain_median_rsa, inter_chain_dist_weighted_rsa,
     intra_chain_median_kd, intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa, inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
     intra_chain_interactions_median, intra_chain_interactions_dist_weighted, rin_class, rin_simple_class, disorder_score, disorder_region]

    for pos, base_name in enumerate(raw_aggregate_header_base_names):
        value = aggregate_values[pos]
        indel_output.add_value('%s %s' % (aggregate_type, base_name), value)

    indel_output.add_value('%s AbsoluteCentrality' % aggregate_type, centrality_scores.AbsoluteCentrality)
    indel_output.add_value('%s LengthNormalizedCentrality' % aggregate_type, centrality_scores.LengthNormalizedCentrality)
    indel_output.add_value('%s MinMaxNormalizedCentrality' % aggregate_type, centrality_scores.MinMaxNormalizedCentrality)
    indel_output.add_value('%s AbsoluteCentralityWithNegative' % aggregate_type, centrality_scores.AbsoluteCentralityWithNegative)
    indel_output.add_value('%s LengthNormalizedCentralityWithNegative' % aggregate_type, centrality_scores.LengthNormalizedCentralityWithNegative)
    indel_output.add_value('%s MinMaxNormalizedCentralityWithNegative' % aggregate_type, centrality_scores.MinMaxNormalizedCentralityWithNegative)
    indel_output.add_value('%s AbsoluteComplexCentrality' % aggregate_type, centrality_scores.AbsoluteComplexCentrality)
    indel_output.add_value('%s LengthNormalizedComplexCentrality' % aggregate_type, centrality_scores.LengthNormalizedComplexCentrality)
    indel_output.add_value('%s MinMaxNormalizedComplexCentrality' % aggregate_type, centrality_scores.MinMaxNormalizedComplexCentrality)
    indel_output.add_value('%s AbsoluteComplexCentralityWithNegative' % aggregate_type, centrality_scores.AbsoluteComplexCentralityWithNegative)
    indel_output.add_value('%s LengthNormalizedComplexCentralityWithNegative' % aggregate_type, centrality_scores.LengthNormalizedComplexCentralityWithNegative)
    indel_output.add_value('%s MinMaxNormalizedComplexCentralityWithNegative' % aggregate_type, centrality_scores.MinMaxNormalizedComplexCentralityWithNegative)

    for chaintype in ['mc', 'sc']:
        for interaction_type in ['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']:
            feature_name = '%s %s %s score' % (aggregate_type, chaintype, interaction_type)
            value = rin_profile.getChainSpecificCombiScore(chaintype, interaction_type)
            indel_output.add_value(feature_name, value)

            feature_name = '%s %s %s degree' % (aggregate_type, chaintype, interaction_type)
            value = rin_profile.getChainSpecificCombiDegree(chaintype, interaction_type)
            indel_output.add_value(feature_name, value)

            feature_name = '%s %s %s H-bond score' % (aggregate_type, chaintype, interaction_type)
            value = rin_profile.getScore(chaintype, 'hbond', interaction_type)
            indel_output.add_value(feature_name, value)

def create_indel_results_table(config, output_path, session_name, session_id):
    table = 'RS_Indel_Session'
    rows = ['Indel', 'Tags']
    eq_rows = {'Session': session_id}

    results = database.select(config, rows, table, equals_rows=eq_rows)

    tag_map = {}
    ids = []
    for row in results:
        indel_id = row[0]
        tags = row[1]
        tag_map[indel_id] = tags
        ids.append(indel_id)

    cols = ['Indel_Id', 'Indel_Notation', 'Analysis_Results', 'Wildtype_Protein', 'Mutant_Protein']
    results = database.binningSelect(ids, cols, 'Indel', config)

    indel_output = out_generator.OutputGenerator()

    headers = ['Indel', 'Protein', 'Tags', 'Size', 'Delta delta classification']

    raw_aggregate_header_base_names = ['RSA', 'Mainchain RSA', 'Sidechain RSA', 'B Factor', 'Modres score', 'Secondary structure assignment', 'Phi', 'Psi',
        'Intra_SSBOND_Propensity', 'Intra_Link_Propensity', 'CIS_Conformation_Propensity', 'CIS_Follower_Propensity', 'Inter Chain Median KD', 'Inter Chain Distance Weighted KD', 'Inter Chain Median RSA', 'Inter Chain Distance Weighted RSA',
        'Intra Chain Median KD', 'Intra Chain Distance Weighted KD', 'Intra Chain Median RSA', 'Intra Chain Distance Weighted RSA',
        'Inter Chain Interactions Median', 'Inter Chain Interactions Distance Weighted',
        'Intra Chain Interactions Median', 'Intra Chain Interactions Distance Weighted', 'RIN-based classification', 'RIN-based simple classification', 'IUPred value', 'Region structure type']

    centrality_header_names = ['AbsoluteCentrality', 'LengthNormalizedCentrality', 'MinMaxNormalizedCentrality',
        'AbsoluteCentralityWithNegative', 'LengthNormalizedCentralityWithNegative', 'MinMaxNormalizedCentralityWithNegative',
        'AbsoluteComplexCentrality', 'LengthNormalizedComplexCentrality', 'MinMaxNormalizedComplexCentrality',
        'AbsoluteComplexCentralityWithNegative', 'LengthNormalizedComplexCentralityWithNegative', 'MinMaxNormalizedComplexCentralityWithNegative', ]

    aggregate_header_base_names = raw_aggregate_header_base_names + centrality_header_names

    for chaintype in ['mc', 'sc']:
        for interaction_type in ['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']:
            feature_name = '%s %s score' % (chaintype, interaction_type)
            aggregate_header_base_names.append(feature_name)
            feature_name = '%s %s degree' % (chaintype, interaction_type)
            aggregate_header_base_names.append(feature_name)
            feature_name = '%s %s H-bond score' % (chaintype, interaction_type)
            aggregate_header_base_names.append(feature_name)

    aggregate_types = ['Wildtype', 'Mutant', 'WT left flank', 'Mut left flank', 'WT right flank', 'Mut right flank']
    for aggregate_type in aggregate_types:
        for header_base in aggregate_header_base_names:
            headers.append('%s %s' % (aggregate_type, header_base))

    indel_output.add_headers(headers)

    indel_file = '%s/%s_indel_analysis.tsv' % (output_path, session_name)

    if os.path.exists(indel_file):
        os.remove(indel_file)

    prot_id_list = set()
    for row in results:
        prot_id_list.add(row[3])

    if len(prot_id_list) == 0:
        return

    protein_dict, gene_db_ids = database.getProteinDict(prot_id_list, session_id, config)

    f = open(indel_file, 'a')
    f.write(indel_output.get_header())

    for row in results:
        indel_id = row[0]
        indel_output.add_value('Indel', row[1])
        indel_output.add_value('Tags', tag_map[indel_id])

        if row[2] is not None:
            (size, ddC, wt_aggregates, mut_aggregates, left_flank_wt_aggregates, left_flank_mut_aggregates, right_flank_wt_aggregates, right_flank_mut_aggregates) = base_utils.unpack(row[2])
        else:
            (size, ddC, wt_aggregates, mut_aggregates, left_flank_wt_aggregates, left_flank_mut_aggregates, right_flank_wt_aggregates, right_flank_mut_aggregates) = [None] * 8

        indel_output.add_value('Size', size)
        indel_output.add_value('Delta delta classification', ddC)

        if wt_aggregates is None:
            print('WT aggregates is None for:', row[1])

        for pos, aggregate in enumerate([wt_aggregates, mut_aggregates, left_flank_wt_aggregates, left_flank_mut_aggregates, right_flank_wt_aggregates, right_flank_mut_aggregates]):
            add_aggregate_results(aggregate, indel_output, aggregate_types[pos], raw_aggregate_header_base_names)

        indel_notation = row[2]
        wt_prot_db_id = row[3]
        (prot_id, u_ac, refseq, u_id, error_code, error, input_id, gene_db_id, _) = protein_dict[wt_prot_db_id]
        indel_output.add_value('Protein', prot_id)
        mut_prot_db_id = row[4]

        f.write(indel_output.pop_line())

    f.close()
    return
