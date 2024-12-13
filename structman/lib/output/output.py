# Contains any classes/functions that are used externally
import os
import time

from structman.lib import babel
from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc.consts import residues as residue_consts
from structman.lib.database import database

try:
    from structman.lib import modelling
except:
    pass
from structman.lib.output import classification as classification_package
from structman.lib.output import feature as feature_package
from structman.lib.output import indel as indel_package
from structman.lib.output import out_utils


def appendOutput(proteins, outfolder, session_name, out_objects=None):
    outfile = '%s/%s' % (outfolder, session_name)
    class_file = "%s.classification.tsv" % (outfile)

    if out_objects is None:

        if os.path.isfile(class_file):
            os.remove(class_file)

        feature_file = "%s.features.tsv" % (outfile)
        if os.path.isfile(feature_file):
            os.remove(feature_file)

        classification_output, f = classification_package.init_classification_table(class_file)
        feature_output, feat_f = feature_package.init_feature_table(feature_file)
    else:

        feature_file = "%s.features.tsv" % (outfile)
        f = open(class_file, 'a')
        feat_f = open(feature_file, 'a')
        classification_output, feature_output = out_objects

    prot_ids = proteins.get_protein_ids()

    for prot_id in prot_ids:
        u_id = proteins.get_u_id(prot_id)
        refseq = proteins.get_ref_id(prot_id)

        u_ac = proteins[prot_id].u_ac
        input_id = proteins[prot_id].input_id

        for pos in proteins.get_position_ids(prot_id):

            m = (prot_id, pos)

            tags = proteins.get_pos_tags(prot_id, pos)

            position = proteins.get_position(prot_id, pos)

            aa1 = position.wt_aa

            if position is None:
                continue

            mappings = position.mappings

            Class = mappings.Class
            conf = mappings.classification_conf
            weighted_sc = mappings.weighted_location
            recommended_structure_str = mappings.get_recommended_res_str()
            recommended_structure, seq_id, cov, resolution = sdsc_utils.process_recommend_structure_str(recommended_structure_str)
            max_seq_structure_str = mappings.get_max_seq_structure_res_str()
            max_seq_structure, max_seq_seq_id, max_seq_cov, max_seq_resolution = sdsc_utils.process_recommend_structure_str(max_seq_structure_str)

            amount_of_structures = len(mappings.qualities)
            mv_sec_ass = mappings.weighted_ssa
            simple_class = mappings.simple_class
            interaction_str = str(mappings.interaction_recommendations)

            input_res_id = proteins.get_res_id(prot_id, pos)

            classification_output.add_value('Input Protein ID', input_id)
            classification_output.add_value('Primary Protein ID', prot_id)
            classification_output.add_value('Uniprot-Ac', u_ac)
            classification_output.add_value('Uniprot ID', u_id)
            classification_output.add_value('Refseq', refseq)
            classification_output.add_value('Residue ID', input_res_id)
            classification_output.add_value('Amino Acid', aa1)
            classification_output.add_value('Position', pos)
            classification_output.add_value('Tags', tags)
            classification_output.add_value('Weighted Location', weighted_sc)
            classification_output.add_value('Weighted Mainchain Location', mappings.weighted_mainchain_location)
            classification_output.add_value('Weighted Sidechain Location', mappings.weighted_sidechain_location)
            classification_output.add_value('Class', Class)
            classification_output.add_value('Simple Class', simple_class)
            classification_output.add_value('RIN Class', mappings.rin_class)
            classification_output.add_value('RIN Simple Class', mappings.rin_simple_class)
            classification_output.add_value('Individual Interactions', interaction_str)
            classification_output.add_value('Confidence Value', conf)
            classification_output.add_value('Secondary Structure', mv_sec_ass)
            classification_output.add_value('Recommended Structure', recommended_structure)
            classification_output.add_value('Sequence-ID', seq_id)
            classification_output.add_value('Coverage', cov)
            classification_output.add_value('Resolution', resolution)
            classification_output.add_value('Max Seq Id Structure', max_seq_structure)
            classification_output.add_value('Max Sequence-ID', max_seq_seq_id)
            classification_output.add_value('Max Seq Id Coverage', max_seq_cov)
            classification_output.add_value('Max Seq Id Resolution', max_seq_resolution)
            classification_output.add_value('Amount of mapped structures', amount_of_structures)

            f.write(classification_output.pop_line())

            iupred = position.get_disorder_score()
            disorder_state = position.get_disorder_region()
            centrality_scores = mappings.rin_based_features.get_centralities()
            rin_profile = mappings.rin_based_features.get_profile()

            mut_aas = position.get_mut_aas()
            snv_tag_map = {}
            for aa2 in mut_aas:
                snv_tag_map[aa2] = position.get_mut_tags(aa2)
            if aa1 not in snv_tag_map:
                snv_tag_map[aa1] = tags

            for new_aa in snv_tag_map:
                mut_tags = snv_tag_map[new_aa]
                aac = "%s%s%s" % (aa1, str(pos), new_aa)
                feature_output.add_value('Input Protein ID', input_id)
                feature_output.add_value('Primary Protein ID', prot_id)
                feature_output.add_value('Uniprot-Ac', u_ac)
                feature_output.add_value('WT Amino Acid', aa1)
                feature_output.add_value('Position', pos)
                feature_output.add_value('Mut Amino Acid', new_aa)
                feature_output.add_value('AA change', '%s%s' % (aa1, new_aa))
                feature_output.add_value('Tags', mut_tags)
                feature_output.add_value('Distance-based classification', Class)
                feature_output.add_value('Distance-based simple classification', simple_class)
                feature_output.add_value('RIN-based classification', mappings.rin_class)
                feature_output.add_value('RIN-based simple classification', mappings.rin_simple_class)
                feature_output.add_value('Classification confidence', conf)
                feature_output.add_value('Structure Location', weighted_sc)
                feature_output.add_value('Mainchain Location', mappings.weighted_mainchain_location)
                feature_output.add_value('Sidechain Location', mappings.weighted_sidechain_location)
                feature_output.add_value('RSA', mappings.weighted_surface_value)
                feature_output.add_value('Mainchain RSA', mappings.weighted_mainchain_surface_value)
                feature_output.add_value('Sidechain RSA', mappings.weighted_sidechain_surface_value)
                feature_output.add_value('Amount of mapped structures', amount_of_structures)
                feature_output.add_value('Secondary structure assignment', mv_sec_ass)
                feature_output.add_value('IUPred value', iupred)
                feature_output.add_value('Region structure type', disorder_state)
                feature_output.add_value('Modres score', mappings.weighted_modres)
                feature_output.add_value('Phi', mappings.weighted_phi)
                feature_output.add_value('Psi', mappings.weighted_psi)

                KDmean = abs(residue_consts.HYDROPATHY[aa1] - residue_consts.HYDROPATHY[new_aa])
                feature_output.add_value('KD mean', KDmean)

                d_vol = abs(residue_consts.VOLUME[aa1] - residue_consts.VOLUME[new_aa])
                feature_output.add_value('Volume mean', d_vol)

                chemical_distance = database.getChemicalDistance(aac)
                feature_output.add_value('Chemical distance', chemical_distance)

                blosum_value = database.getBlosumValue(aac)
                feature_output.add_value('Blosum62', aac)

                aliphatic_change = int((aa1 in residue_consts.AA_MAP_ALIPHATIC) != (new_aa in residue_consts.AA_MAP_ALIPHATIC))
                hydrophobic_change = int((aa1 in residue_consts.AA_MAP_HYDROPHOBIC) != (new_aa in residue_consts.AA_MAP_HYDROPHOBIC))
                aromatic_change = int((aa1 in residue_consts.AA_MAP_AROMATIC) != (new_aa in residue_consts.AA_MAP_AROMATIC))
                positive_change = int((aa1 in residue_consts.AA_MAP_POSITIVE) != (new_aa in residue_consts.AA_MAP_POSITIVE))
                polar_change = int((aa1 in residue_consts.AA_MAP_POLAR) != (new_aa in residue_consts.AA_MAP_POLAR))
                negative_change = int((aa1 in residue_consts.AA_MAP_NEGATIVE) != (new_aa in residue_consts.AA_MAP_NEGATIVE))
                charged_change = int((aa1 in residue_consts.AA_MAP_CHARGED) != (new_aa in residue_consts.AA_MAP_CHARGED))
                small_change = int((aa1 in residue_consts.AA_MAP_SMALL) != (new_aa in residue_consts.AA_MAP_SMALL))
                tiny_change = int((aa1 in residue_consts.AA_MAP_TINY) != (new_aa in residue_consts.AA_MAP_TINY))
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
                feature_output.add_value('B Factor', mappings.weighted_b_factor)

                if centrality_scores is not None:
                    feature_output.add_value('AbsoluteCentrality', centrality_scores.AbsoluteCentrality)
                    feature_output.add_value('LengthNormalizedCentrality', centrality_scores.LengthNormalizedCentrality)
                    feature_output.add_value('MinMaxNormalizedCentrality', centrality_scores.MinMaxNormalizedCentrality)
                    feature_output.add_value('AbsoluteCentralityWithNegative', centrality_scores.AbsoluteCentralityWithNegative)
                    feature_output.add_value('LengthNormalizedCentralityWithNegative', centrality_scores.LengthNormalizedCentralityWithNegative)
                    feature_output.add_value('MinMaxNormalizedCentralityWithNegative', centrality_scores.MinMaxNormalizedCentralityWithNegative)
                    feature_output.add_value('AbsoluteComplexCentrality', centrality_scores.AbsoluteComplexCentrality)
                    feature_output.add_value('LengthNormalizedComplexCentrality', centrality_scores.LengthNormalizedComplexCentrality)
                    feature_output.add_value('MinMaxNormalizedComplexCentrality', centrality_scores.MinMaxNormalizedComplexCentrality)
                    feature_output.add_value('AbsoluteComplexCentralityWithNegative', centrality_scores.AbsoluteComplexCentralityWithNegative)
                    feature_output.add_value('LengthNormalizedComplexCentralityWithNegative', centrality_scores.LengthNormalizedComplexCentralityWithNegative)
                    feature_output.add_value('MinMaxNormalizedComplexCentralityWithNegative', centrality_scores.MinMaxNormalizedComplexCentralityWithNegative)

                feature_output.add_value('Intra_SSBOND_Propensity', mappings.weighted_intra_ssbond)
                feature_output.add_value('Inter_SSBOND_Propensity', mappings.weighted_inter_ssbond)
                feature_output.add_value('Intra_Link_Propensity', mappings.weighted_intra_link)
                feature_output.add_value('Inter_Link_Propensity', mappings.weighted_inter_link)
                feature_output.add_value('CIS_Conformation_Propensity', mappings.weighted_cis_conformation)
                feature_output.add_value('CIS_Follower_Propensity', mappings.weighted_cis_follower)
                feature_output.add_value('Inter Chain Median KD', mappings.weighted_inter_chain_median_kd)
                feature_output.add_value('Inter Chain Distance Weighted KD', mappings.weighted_inter_chain_dist_weighted_kd)
                feature_output.add_value('Inter Chain Median RSA', mappings.weighted_inter_chain_median_rsa)
                feature_output.add_value('Inter Chain Distance Weighted RSA', mappings.weighted_inter_chain_dist_weighted_rsa)
                feature_output.add_value('Intra Chain Median KD', mappings.weighted_intra_chain_median_kd)
                feature_output.add_value('Intra Chain Distance Weighted KD', mappings.weighted_intra_chain_dist_weighted_kd)
                feature_output.add_value('Intra Chain Median RSA', mappings.weighted_intra_chain_median_rsa)
                feature_output.add_value('Intra Chain Distance Weighted RSA', mappings.weighted_intra_chain_dist_weighted_rsa)
                feature_output.add_value('Inter Chain Interactions Median', mappings.weighted_inter_chain_interactions_median)
                feature_output.add_value('Inter Chain Interactions Distance Weighted', mappings.weighted_inter_chain_interactions_dist_weighted)
                feature_output.add_value('Intra Chain Interactions Median', mappings.weighted_intra_chain_interactions_median)
                feature_output.add_value('Intra Chain Interactions Distance Weighted', mappings.weighted_intra_chain_interactions_dist_weighted)

                for chaintype in ['mc', 'sc']:
                    for interaction_type in ['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']:
                        feature_name = '%s %s score' % (chaintype, interaction_type)
                        value = rin_profile.getChainSpecificCombiScore(chaintype, interaction_type)
                        feature_output.add_value(feature_name, value)

                        feature_name = '%s %s degree' % (chaintype, interaction_type)
                        value = rin_profile.getChainSpecificCombiDegree(chaintype, interaction_type)
                        feature_output.add_value(feature_name, value)

                        feature_name = '%s %s H-bond score' % (chaintype, interaction_type)
                        value = rin_profile.getScore(chaintype, 'hbond', interaction_type)
                        feature_output.add_value(feature_name, value)

                feat_f.write(feature_output.pop_line())

    f.close()
    feat_f.close()

    return classification_output, feature_output


def makeDiffDict(f):
    f = open(f, "r")
    lines = f.readlines()[1:]
    f.close()
    diff_dict = {}
    for line in lines:
        words = line.split("\t")
        diff_dict[words[1]] = [words[0], words[5]]
    return diff_dict


def makeGeneDict(f):
    f = open(f, "r")
    lines = f.readlines()[1:]
    f.close()
    diff_dict = {}
    for line in lines:
        words = line.split("\t")
        diff_dict[words[0]] = words[2]
    return diff_dict


def godiffAna(fileA, fileB, papermode=False):
    go_dict_a = makeDiffDict(fileA)
    go_dict_b = makeDiffDict(fileB)
    result_list = []
    for Id in go_dict_a:
        if Id in go_dict_b:
            diff = float(go_dict_a[Id][1]) - float(go_dict_b[Id][1])
            result_list.append([Id, go_dict_a[Id][0], diff])
    result_list = sorted(result_list, key=lambda x: x[2], reverse=True)

    lines = ["GO-Term\tGO-ID\tScore-Difference"]
    for result in result_list:
        if not papermode:
            lines.append("%s\t%s\t%s" % (result[1], result[0], str(result[2])))
        elif result[1][:2] == ' P':
            lines.append("%s%s\t%s\t%s" % (result[1][3].upper(), result[1][4:], result[0], str(result[2])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    name_b = (fileB.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    outfile = "%s/GO_Diff_%s_%s.tsv" % (fileA.rsplit("/", 1)[0], name_a, name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()


def pathdiffAna(fileA, fileB, papermode=False):
    path_dict_a = makeDiffDict(fileA)
    path_dict_b = makeDiffDict(fileB)
    result_list = []
    for Id in path_dict_a:
        if Id in path_dict_b:
            diff = float(path_dict_a[Id][1]) - float(path_dict_b[Id][1])
            result_list.append([Id, path_dict_a[Id][0], diff])
    result_list = sorted(result_list, key=lambda x: x[2], reverse=True)

    lines = ["Pathway\tReactome-ID\tScore-Difference"]
    for result in result_list:
        lines.append("%s\t%s\t%s" % (result[1], result[0], str(result[2])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    name_b = (fileB.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    outfile = "%s/Path_Diff_%s_%s.tsv" % (fileA.rsplit("/", 1)[0], name_a, name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()


def genediffAna(fileA, fileB):
    gene_dict_a = makeGeneDict(fileA)
    gene_dict_b = makeGeneDict(fileB)
    result_list = []
    for Id in gene_dict_a:
        if Id in gene_dict_b:
            diff = float(gene_dict_a[Id]) - float(gene_dict_b[Id])
            result_list.append([Id, diff])
    result_list = sorted(result_list, key=lambda x: x[1], reverse=True)

    lines = ["Gene\tScore-Difference"]
    for result in result_list:
        lines.append("%s\t%s" % (result[0], str(result[1])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    name_b = (fileB.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    outfile = "%s/Gene_Diff_%s_%s.tsv" % (fileA.rsplit("/", 1)[0], name_a, name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()


# called by structman
def main(sess_id, output_path, config, intertable=False):
    db_name = config.db_name
    db_address = config.db_address
    db_password = config.db_password
    db_user_name = config.db_user_name
    go = config.go
    godiff = config.godiff
    classification = config.classification
    path = config.path
    pathdiff = config.pathdiff
    do_modelling = config.do_modelling
    ligand_file = config.ligand_file
    multi_modelling = config.multi_modelling

    mod_per_gene = config.mod_per_gene
    infile = ''
    tanimoto_cutoff = config.tanimoto_cutoff
    distance_threshold = config.milieu_threshold

    intertable_conf = config.intertable_conf
    intertable = intertable or intertable_conf

    if sess_id == 0 or sess_id is None:
        session_id = database.getSessionId(infile, config)
    else:
        session_id = sess_id

    db, cursor = config.getDB()

    if infile == '':
        infile = database.getSessionFile(session_id, config)

    db.close()

    session_name = (infile.rsplit("/", 1)[1]).rsplit(".", 1)[0]

    t0 = time.time()

    if classification and not config.skip_main_output_generation:
        t00 = time.time()
        classfiles, interfiles = classification_package.classificationOutput(config, output_path, session_name, session_id)
        t01 = time.time()
        if config.verbosity >= 2:
            print("Time for classificationOutput: ", t01 - t00)
        #for classfile in classfiles:
        #    out_utils.classDistributionFromFile(classfile, output_path, session_name, config)
        #    out_utils.classDistributionFromFile(classfile, output_path, session_name, config, rin_classes=True)
        t02 = time.time()
        if config.verbosity >= 2:
            print("Time for producing classification distributions: ", t02 - t01)

        if intertable:
            for interfile in interfiles:
                out_utils.InteractionScoreAveragesFromFile(interfile, output_path, session_name, by_tag=True)
            t03 = time.time()
            if config.verbosity >= 2:
                print("Time for producing Interaction files: ", t03 - t02)

    t1 = time.time()
    if config.verbosity >= 2:
        print("Time for producing classification file: ", t1 - t0)

    #if config.indels_given_by_input:
        #if config.verbosity >= 2:
        #    print('Starting indel analysis output generation')
        #indel_package.create_indel_results_table(config, output_path, session_name, session_id)

    db, cursor = config.getDB()

    if go:
        database.goTermAnalysis(session_id, "%s/%s.goterm.tsv" % (output_path, session_name), db, cursor)
        if godiff:
            files = os.listdir(output_path)
            go_files = []
            for f in files:
                if '.goterm.tsv' in f:
                    go_files.append(f)
            print(go_files)
            print(files)
            if len(go_files) == 2:
                fileA = "%s/%s" % (output_path, go_files[0])
                fileB = "%s/%s" % (output_path, go_files[1])
                godiffAna(fileA, fileB)
    if path:
        database.pathwayAnalysis(session_id, "%s/%s.pathway.tsv" % (output_path, session_name), db, cursor)
        if pathdiff:
            files = os.listdir(output_path)
            path_files = []
            for f in files:
                if '.pathway.tsv' in f:
                    path_files.append(f)
            if len(path_files) == 2:
                fileA = "%s/%s" % (output_path, path_files[0])
                fileB = "%s/%s" % (output_path, path_files[1])
                pathdiffAna(fileA, fileB)
    if do_modelling:
        modelling.massModel(session_id, db, cursor, output_path, total_models=0, model_per_gene=int(mod_per_gene), multiple_mutations=multi_modelling)

    cursor.close()
    db.close()
    if ligand_file is None:
        try:
            ligand_file_names = os.listdir("%s/ligands" % infile.rsplit("/", 1)[0])
            ligand_files = []
            for ligand_file_name in ligand_file_names:
                ligand_files.append("%s/ligands/%s" % (infile.rsplit("/", 1)[0], ligand_file_name))
        except:
            ligand_files = []
    else:
        ligand_files = [ligand_file]
    for ligand_file in ligand_files:
        t0 = time.time()
        anno_dict = babel.ligandAnalyzer(ligand_file, session_id, db_name, db_address, db_user_name, db_password, cutoff=tanimoto_cutoff, distance_threshold=distance_threshold)
        t1 = time.time()
        babel.writeReport(anno_dict, "%s/Ligand_Report_%s_%s.tsv" % (output_path, ligand_file.rsplit('/', 1)[1].rsplit(".", 1)[0], session_name), db_name, db_address, db_user_name, db_password)
        t2 = time.time()
        if config.verbosity >= 2:
            print("Time for ligandAnalyzer: ", t1 - t0)
            print("Time for writeReport: ", t2 - t1)
