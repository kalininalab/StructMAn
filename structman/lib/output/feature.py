from structman.lib.output import out_generator

def init_feature_table(feature_file, obj_only = False):
    feature_output = out_generator.OutputGenerator()
    headers = [
        'Input Protein ID', 'Primary Protein ID', 'Uniprot-Ac', 'WT Amino Acid', 'Position', 'Mut Amino Acid', 'AA change', 'Tags',
        'Distance-based classification', 'Distance-based simple classification',
        'RIN-based classification', 'RIN-based simple classification',
        'Classification confidence', 'Structure Location', 'Mainchain Location', 'Sidechain Location',
        'RSA', 'Mainchain RSA', 'Sidechain RSA', 'RSA change score', 'Mainchain RSA change score', 'Sidechain RSA change score', 'Amount of mapped structures',
        'Secondary structure assignment', 'IUPred value', 'Region structure type', 'Modres score',
        'Phi', 'Psi', 'KD mean',
        'Volume mean', 'Chemical distance', 'Blosum62',
        'Aliphatic change', 'Hydrophobic change', 'Aromatic change', 'Positive charged change',
        'Polar change', 'Negative charge change', 'Charged change', 'Small change', 'Tiny change', 'Total change',
        'B Factor',
        'AbsoluteCentrality', 'LengthNormalizedCentrality', 'MinMaxNormalizedCentrality',
        'AbsoluteCentralityWithNegative', 'LengthNormalizedCentralityWithNegative', 'MinMaxNormalizedCentralityWithNegative',
        'AbsoluteComplexCentrality', 'LengthNormalizedComplexCentrality', 'MinMaxNormalizedComplexCentrality',
        'AbsoluteComplexCentralityWithNegative', 'LengthNormalizedComplexCentralityWithNegative', 'MinMaxNormalizedComplexCentralityWithNegative',
        'Intra_SSBOND_Propensity', 'Inter_SSBOND_Propensity', 'Intra_Link_Propensity', 'Inter_Link_Propensity',
        'CIS_Conformation_Propensity', 'CIS_Follower_Propensity',
        'Inter Chain Median KD', 'Inter Chain Distance Weighted KD', 'Inter Chain Median RSA', 'Inter Chain Distance Weighted RSA',
        'Intra Chain Median KD', 'Intra Chain Distance Weighted KD', 'Intra Chain Median RSA', 'Intra Chain Distance Weighted RSA',
        'Inter Chain Interactions Median', 'Inter Chain Interactions Distance Weighted',
        'Intra Chain Interactions Median', 'Intra Chain Interactions Distance Weighted',
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

    for chaintype in ['mc', 'sc']:
        for interaction_type in ['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']:
            feature_name = '%s %s score' % (chaintype, interaction_type)
            headers.append(feature_name)
            feature_name = '%s %s degree' % (chaintype, interaction_type)
            headers.append(feature_name)
            feature_name = '%s %s H-bond score' % (chaintype, interaction_type)
            headers.append(feature_name)

    feature_output.add_headers(headers)
    if obj_only:
        return feature_output

    feat_f = open(feature_file, 'a')
    feat_f.write(feature_output.get_header())
    return feature_output, feat_f
