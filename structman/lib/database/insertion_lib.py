import time
import multiprocessing

from structman.lib.pdbParser import parseLigandDB, getSI, updateLigandDB
from structman.lib.database.retrieval import getStoredResidues
from structman.lib.sdsc import structure as structure_package
from structman.base_utils.base_utils import pack
from structman.lib.database.database_core_functions import binningSelect, insert, select, update

def insert_interfaces(proteins, config):
    interface_values = []
    prot_prot_values = []
    prot_prot_map = {}
    pos_pos_values = []
    prot_db_ids = set()
    structure_ids_for_residue_retrieval = {}
    pos_pos_pre_values = []
    for protein_id in proteins.protein_map:
        if config.verbosity >= 5:
            print('In insert_interfaces:', protein_id)
        prot_a_db_id = proteins.protein_map[protein_id].database_id
        for aggregated_interface in proteins.protein_map[protein_id].aggregated_interface_map:
            prot_db_ids.add(prot_a_db_id)
            structure_recommendation = f'{aggregated_interface.recommended_complex},{aggregated_interface.chain},{aggregated_interface.interacting_chain}'
            if config.verbosity >= 5:
                print('In insert_interfaces:', protein_id, structure_recommendation)
                aggregated_interface.print_interface()

            #structure_ids_for_residue_retrieval[proteins.structures[(aggregated_interface.recommended_complex, aggregated_interface.chain)].database_id] = (aggregated_interface.recommended_complex, aggregated_interface.chain)

            interface_hash_tuple = aggregated_interface.get_interface_hash()

            if aggregated_interface.database_id is None:
                mi, ma, me, le = interface_hash_tuple
                interface_values.append((prot_a_db_id, structure_recommendation, mi, ma, me, le))
            for pos_a in aggregated_interface.pos_pos_interactions:
                for (i_prot, i_pos) in aggregated_interface.pos_pos_interactions[pos_a]:
                    pos_pos_interaction = aggregated_interface.pos_pos_interactions[pos_a][(i_prot, i_pos)]

                    try:
                        pos_a_db_id = proteins.protein_map[protein_id].positions[pos_a].database_id
                    except:
                        config.errorlog.add_warning(f'Position database id not found: {protein_id} {pos_a}')
                        continue
                    if pos_pos_interaction.protein_b is None or pos_pos_interaction.position_b is None:
                        continue
                    pos_b_db_id = proteins.protein_map[pos_pos_interaction.protein_b].positions[pos_pos_interaction.position_b].database_id
                    
                    (chain_a, res_a, chain_b, res_b) = pos_pos_interaction.recommended_interaction


                    #structure_ids_for_residue_retrieval[proteins.structures[(pos_pos_interaction.recommended_complex, chain_a)].database_id] = (pos_pos_interaction.recommended_complex, chain_a)
                    #structure_ids_for_residue_retrieval[proteins.structures[(pos_pos_interaction.recommended_complex, chain_b)].database_id] = (pos_pos_interaction.recommended_complex, chain_b)

                    pos_pos_pre_values.append((pos_a_db_id, pos_b_db_id, pos_pos_interaction.recommended_complex, chain_a, chain_b, res_a, res_b, pos_pos_interaction.interaction_score))

                    if not pos_pos_interaction.protein_b in prot_prot_map:
                        pos_pos_pre_values.append((pos_b_db_id, pos_a_db_id, pos_pos_interaction.recommended_complex, chain_b, chain_a, res_b, res_a, pos_pos_interaction.interaction_score))

                    if not prot_a_db_id in prot_prot_map:
                        prot_prot_map[prot_a_db_id] = {}

                    prot_b_db_id = proteins.protein_map[pos_pos_interaction.protein_b].database_id

                    prot_db_ids.add(prot_b_db_id)

                    if not prot_b_db_id in prot_prot_map[prot_a_db_id]:
                        prot_prot_map[prot_a_db_id][prot_b_db_id] = [interface_hash_tuple, None, {}, 0]
                    prot_prot_map[prot_a_db_id][prot_b_db_id][3] += pos_pos_interaction.interaction_score

                    if not (pos_pos_interaction.recommended_complex, chain_a, chain_b) in prot_prot_map[prot_a_db_id][prot_b_db_id][2]:
                        prot_prot_map[prot_a_db_id][prot_b_db_id][2][(pos_pos_interaction.recommended_complex, chain_a, chain_b)] = 0

                    prot_prot_map[prot_a_db_id][prot_b_db_id][2][(pos_pos_interaction.recommended_complex, chain_a, chain_b)] += 1

                    if prot_b_db_id in prot_prot_map:
                        if prot_a_db_id in prot_prot_map[prot_b_db_id]:
                            prot_prot_map[prot_b_db_id][prot_a_db_id][1] = interface_hash_tuple

    for structure_id in proteins.complexes:
        for chain in proteins.complexes[structure_id].chains:
            if proteins.complexes[structure_id].chains[chain] == 'Protein' and (structure_id, chain) in proteins.structures:
                if len(proteins.structures[(structure_id, chain)].residues) > 0: #If there are residues inside the structure datastructure ...
                    break_it = False
                    for res in proteins.structures[(structure_id, chain)].residues:
                        if proteins.structures[(structure_id, chain)].residues[res].database_id != None: #look if the first residue has a database id ...
                            break_it = True
                        break
                    if break_it: #If yes, don't queue it for retrieval
                        continue
                structure_ids_for_residue_retrieval[proteins.structures[(structure_id, chain)].database_id] = (structure_id, chain)

    getStoredResidues(proteins, config, custom_ids = structure_ids_for_residue_retrieval, retrieve_only_db_ids = True)

    for (pos_a_db_id, pos_b_db_id, recommended_complex, chain_a, chain_b, res_a, res_b, interaction_score) in pos_pos_pre_values:
        if res_a not in proteins.structures[(recommended_complex, chain_a)].residues:
            continue
        res_a_db_id = proteins.structures[(recommended_complex, chain_a)].residues[res_a].database_id
        if not (recommended_complex, chain_b) in proteins.structures:
            continue
        if not res_b in proteins.structures[(recommended_complex, chain_b)].residues:
            continue
        res_b_db_id = proteins.structures[(recommended_complex, chain_b)].residues[res_b].database_id
        pos_pos_values.append((pos_a_db_id, pos_b_db_id, res_a_db_id, res_b_db_id, interaction_score))

    insert('Interface', ['Protein', 'Structure_Recommendation', 'First_Position', 'Last_Position', 'Mean_Position', 'Interface_Size'], interface_values, config)

    results = binningSelect(prot_db_ids, ['Protein', 'Interface_Id', 'First_Position', 'Last_Position', 'Mean_Position', 'Interface_Size'], 'Interface', config)

    blank_interfaces = {}
    for row in results:
        prot_db_id = row[0]
        interface_db_id = row[1]

        interface_hash_tuple = row[2:6]
        if not prot_db_id in blank_interfaces:
            blank_interfaces[prot_db_id] = {}

        blank_interfaces[prot_db_id][interface_hash_tuple] = interface_db_id

    if config.verbosity >= 5:
        print(f'In insert_interfaces: contents of prot_prot_map before transformation to prot_prot_values:\n{prot_prot_map}')

    for prot_a_db_id in prot_prot_map:
        for prot_b_db_id in prot_prot_map[prot_a_db_id]:

            interface_hash_a, interface_hash_b, pos_pos_recommendations, ppi_score = prot_prot_map[prot_a_db_id][prot_b_db_id]

            max_count = 0
            for recommended_complex, chain_a, chain_b in pos_pos_recommendations:
                count = pos_pos_recommendations[(recommended_complex, chain_a, chain_b)]
                if count > max_count:
                    best_recommended_complex, best_chain_a, best_chain_b = recommended_complex, chain_a, chain_b
                    max_count = count

            recommended_complex, chain_a, chain_b = best_recommended_complex, best_chain_a, best_chain_b


            interface_a_db_id = blank_interfaces[prot_a_db_id][interface_hash_a]
            if interface_hash_b is not None:
                interface_b_db_id = blank_interfaces[prot_b_db_id][interface_hash_b]
            elif len(blank_interfaces[prot_b_db_id]) == 1:
                interface_b_db_id = list(blank_interfaces[prot_b_db_id].values())[0]
            else:
                #TODO retrace correct interface
                interface_b_db_id = list(blank_interfaces[prot_b_db_id].values())[0]

            prot_prot_values.append((interface_a_db_id, interface_b_db_id, proteins.complexes[recommended_complex].database_id, chain_a, chain_b, ppi_score))

            if prot_b_db_id not in prot_prot_map:
                prot_prot_values.append((interface_b_db_id, interface_a_db_id, proteins.complexes[recommended_complex].database_id, chain_b, chain_a, ppi_score))

    if config.verbosity >= 5:
        print(f'In insert_interfaces: contents of prot_prot_values:\n{prot_prot_values}')

    position_interface_values = []

    for protein_id in proteins.protein_map:
        prot_a_db_id = proteins.protein_map[protein_id].database_id
        if prot_a_db_id not in blank_interfaces:
            continue
        for aggregated_interface in proteins.protein_map[protein_id].aggregated_interface_map:
            if aggregated_interface.database_id is not None:
                continue

            interface_db_id = blank_interfaces[prot_a_db_id][aggregated_interface.get_interface_hash()]

            for position in aggregated_interface.positions:
                position_db_id = proteins.protein_map[protein_id].positions[position].database_id
                structure_id, chain, res = aggregated_interface.positions[position]
                position_interface_values.append((interface_db_id, position_db_id, f'{structure_id}:{chain} {res}'))

    insert('RS_Position_Interface', ['Interface', 'Position', 'Recommended_Residue'], position_interface_values, config)

    insert('Position_Position_Interaction', ['Position_A', 'Position_B', 'Residue_A', 'Residue_B', 'Interaction_Score'], pos_pos_values, config)

    insert('Protein_Protein_Interaction', ['Interface_A', 'Interface_B', 'Complex', 'Chain_A', 'Chain_B', 'Interaction_Score'], prot_prot_values, config)

    residue_interface_values = []

    print_count = 0
 
    for structure_id in proteins.complexes:
        #print(structure_id, proteins.complexes[structure_id].interfaces)
        for (chain, i_chain) in proteins.complexes[structure_id].interfaces:
            if proteins.complexes[structure_id].interfaces[(chain, i_chain)].stored:
                continue
            try:
                structure_db_id = proteins.structures[(structure_id, chain)].database_id
            except:
                config.errorlog.add_warning(f'Structure missing: {structure_id} {chain} {i_chain}')
                continue
            for res in proteins.complexes[structure_id].interfaces[(chain, i_chain)].residues:

                if not res in proteins.structures[(structure_id, chain)].residues:
                    if config.verbosity >= 4:
                        if print_count < 100:
                            print(f'Residue {res} filtered from {structure_id} {chain} {i_chain} in insert_interfaces: not in structures object')
                            print_count += 1
                    continue #Happens for ligands declared as part of the chain
                res_db_id = proteins.structures[(structure_id, chain)].residues[res].database_id
                if res_db_id is None: 
                    if config.verbosity >= 4:
                        print(f'Residue {res} filtered from {structure_id} {chain} {i_chain} in insert_interfaces: db_id is None')
                    continue
                if res in proteins.complexes[structure_id].interfaces[(chain, i_chain)].interactions:
                    i_res, score = proteins.complexes[structure_id].interfaces[(chain, i_chain)].interactions[res]
                    i_res_string = f'{i_chain},{i_res},{score}'
                else:
                    i_res_string = f'{i_chain},,'
                residue_interface_values.append((structure_db_id, res_db_id, i_res_string))

    insert('RS_Residue_Interface', ['Structure', 'Residue', 'Interacting_Residue'], residue_interface_values, config)


def insertInteractingChains(interaction_structures, proteins, config):
    if len(interaction_structures) == 0:

        return {}

    interacting_structure_ids = {}

    results = select(config, ['Structure_Id', 'PDB', 'Chain'], 'Structure')

    for row in results:
        pdb_id = row[1]
        chain = row[2]
        if not (pdb_id, chain) in interaction_structures:
            continue
        #s_id = row[0]
        #interacting_structure_ids[(pdb_id,chain)] = s_id
        interaction_structures.remove((pdb_id, chain))

    values = []

    for (pdb_id, chain) in interaction_structures:

        homomers = proteins.get_complex_homomers(pdb_id, chain)

        oligos = ''.join(homomers)

        values.append((pdb_id, chain, oligos))

    if len(values) > 0:
        insert('Structure', ['PDB', 'Chain', 'Homooligomer'], values, config)

    if len(interaction_structures) > 0:
        results = select(config, ['Structure_Id', 'PDB', 'Chain'], 'Structure')

        for row in results:
            pdb_id = row[1]
            chain = row[2]
            if not (pdb_id, chain) in interaction_structures:
                continue
            s_id = row[0]

            interacting_structure_ids[(pdb_id, chain)] = s_id

            proteins.structures[(pdb_id, chain)] = structure_package.Structure(pdb_id, chain, database_id=s_id, new_interacting_chain = True)

    return interacting_structure_ids


def background_insert_residues(values, config):
    """
    columns = ['Structure', 'Number', 'Amino_Acid', 'Sub_Lig_Dist', 'Sub_Chain_Distances', 'Relative_Surface_Access',
               'Relative_Surface_Access_Main_Chain', 'Relative_Surface_Access_Side_Chain', 'Secondary_Structure_Assignment',
               'Homomer_Distances', 'Interaction_Profile', 'Centralities', 'B_Factor', 'Modres', 'PHI', 'PSI', 'Intra_SSBOND', 'SSBOND_Length',
               'Intra_Link', 'Link_Length', 'CIS_Conformation', 'CIS_Follower', 'Inter_Chain_Median_KD', 'Inter_Chain_Dist_Weighted_KD', 'Inter_Chain_Median_RSA',
               'Inter_Chain_Dist_Weighted_RSA', 'Intra_Chain_Median_KD', 'Intra_Chain_Dist_Weighted_KD', 'Intra_Chain_Median_RSA', 'Intra_Chain_Dist_Weighted_RSA',
               'Inter_Chain_Interactions_Median', 'Inter_Chain_Interactions_Dist_Weighted',
               'Intra_Chain_Interactions_Median', 'Intra_Chain_Interactions_Dist_Weighted',
               'Interacting_Chains', 'Interacting_Ligands'
               ]
    """
    columns = ['Structure', 'Number', 'Residue_Data']
    insert('Residue', columns, values, config)


def insertResidues(structural_analysis, interacting_structure_ids, proteins, config):

    if config.verbosity >= 2:
        t0 = time.time()

    values = []

    if len(structural_analysis) == 0:
        if config.verbosity >= 4:
            print('insertResidues was called with an empty structural_analysis dict')
        return

    structure_ids = {}

    t_0 = 0.
    t_1 = 0.
    t_2 = 0.
    t_3 = 0.
    t_4 = 0.
    t_5 = 0.
    t_6 = 0.
    t_7 = 0.
    t_8 = 0.
    t_9 = 0.
    t_10 = 0.
    t_11 = 0.
    t_12 = 0.
    t_13 = 0.

    for (pdb_id, chain) in structural_analysis:
        if (pdb_id, chain) in interacting_structure_ids:
            proteins.structures[(pdb_id, chain)].residues = structural_analysis[(pdb_id, chain)]
            if config.verbosity >= 4:
                print(f'structural_analysis ({len(structural_analysis[(pdb_id, chain)])}) added to {pdb_id} {chain} in insertResidues due to part of interacting_structure_ids')
        if not proteins.contains_structure(pdb_id, chain):
            if config.verbosity >= 5:
                if not (pdb_id, chain) in interacting_structure_ids:
                    print(f'{pdb_id} {chain} from structural_analysis not going into the database. It is not in interacting_structure_ids.')
                else:

                    print(f'{pdb_id} {chain} from structural_analysis not going into the database. contains_structure was False.')
            continue
        if proteins.contains_structure(pdb_id, chain) and proteins.is_structure_stored(pdb_id, chain) and not proteins.structures[(pdb_id, chain)].new_interacting_chain:  # all residues belonging to stored structures must not inserted twice
            if config.verbosity >= 5:
                print(f'{pdb_id} {chain} from structural_analysis not going into the database. Is already stored.')
            continue

        analysis_map = structural_analysis[(pdb_id, chain)]
        if (pdb_id, chain) in interacting_structure_ids:
            s_id = interacting_structure_ids[(pdb_id, chain)]
        else:
            s_id = proteins.get_structure_db_id(pdb_id, chain)
            if s_id is None:
                if config.verbosity >= 5:
                    print(f'{pdb_id} {chain} from structural_analysis not going into the database. s_id was None.')
                continue

            structure_ids[s_id] = (pdb_id, chain)

        if config.verbosity >= 5:
            print(f'{pdb_id} {chain} from structural_analysis going into the database: {len(analysis_map)}')

        for res_id in analysis_map:
            t_0 += time.time()
            residue = analysis_map[res_id]

            one_letter = residue.get_aa()
            lig_dist_str = residue.get_lig_dist_str()
            chain_dist_str = residue.get_chain_dist_str()

            t_1 += time.time()

            (rsa, relative_main_chain_acc, relative_side_chain_acc) = residue.get_rsa(splitted=True)
            ssa = residue.get_ssa()

            t_2 += time.time()

            homo_str = residue.get_homo_dist_str()

            t_3 += time.time()

            profile_str = residue.get_interaction_profile_str()

            t_4 += time.time()

            interacting_chains_str = residue.get_interacting_chains_str()
            interacting_ligands_str = residue.get_interacting_ligands_str()

            t_5 += time.time()

            #if interacting_chains is None:
            #    interacting_chains_str = None
            #else:
            #    interacting_chains_str = ','.join(interacting_chains)

            #if interacting_ligands is None:
            #    interacting_ligands_str = None
            #else:
            #    interacting_ligands_str = ','.join(interacting_ligands)

            t_6 += time.time()

            centrality_score_str = residue.get_centrality_str()

            t_7 += time.time()

            b_factor = residue.get_b_factor()
            modres = residue.get_modres()
            phi, psi = residue.get_angles()

            t_8 += time.time()

            intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower = residue.get_residue_link_information()

            t_9 += time.time()

            (inter_chain_median_kd, inter_chain_dist_weighted_kd, inter_chain_median_rsa, inter_chain_dist_weighted_rsa,
                intra_chain_median_kd, intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa) = residue.get_milieu()

            t_10 += time.time()

            (inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                intra_chain_interactions_median, intra_chain_interactions_dist_weighted) = residue.get_interface_milieu()

            t_11 += time.time()

            packed_res_info = pack((one_letter, lig_dist_str, chain_dist_str, rsa, relative_main_chain_acc, relative_side_chain_acc,
                           ssa, homo_str, profile_str,
                           centrality_score_str, b_factor, modres, phi, psi, intra_ssbond, ssbond_length, intra_link, link_length,
                           cis_conformation, cis_follower, inter_chain_median_kd, inter_chain_dist_weighted_kd,
                           inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                           intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                           inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                           intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                           interacting_chains_str, interacting_ligands_str))
            
            values.append([s_id, res_id, packed_res_info])
            t_12 += time.time()

            if not (pdb_id, chain) in interacting_structure_ids:
                if config.verbosity >= 5:
                    print(f'Calling proteins.add_residue in insertResidues: {pdb_id} {chain} {res_id}')
                proteins.add_residue(pdb_id, chain, res_id, residue)

            t_13 += time.time()

    if config.verbosity >= 2:
        t1 = time.time()
        print('Time for insertResidues part 1:', t1 - t0)
        print('Individual parts: 1 -', (t_1 - t_0), '2 -', (t_2 - t_1), '3 -', (t_3 - t_2), '4 -', (t_4 - t_3), '5 -', (t_5 - t_4), '6 -', (t_6 - t_5),
                '7 -', (t_7 - t_6), '8 -', (t_8 - t_7), '9 -', (t_9 - t_8), '10 -', (t_10 - t_9), '11 -', (t_11 - t_10), '12 -', (t_12 - t_11), '13 -', (t_13 - t_12))

    process = multiprocessing.Process(target=background_insert_residues, args=(values, config))
    process.start()
    #process = None
    #background_insert_residues(values, config)
    return process

def insertClassifications(proteins, config):

    if config.verbosity >= 2:
        t1 = time.time()

    table = 'Position'
    columns = ['Position_Id', 'Recommended_Structure_Data', 'Position_Data']

    values = createClassValues(proteins)

    if config.verbosity >= 2:
        t2 = time.time()
        print('insertClassifications part 1: %s' % (t2 - t1), 'Update Classifications of', len(values), 'positions')

    update(config, table, columns, values)

    if config.verbosity >= 2:
        t3 = time.time()
        print('insertClassifications part 2: %s' % (t3 - t2))


def createClassValues(proteins):
    values = []

    for u_ac in proteins.get_protein_ids():
        if proteins.is_completely_stored(u_ac):
            continue
        positions = proteins.get_position_ids(u_ac)
        for pos in positions:
            aachange = proteins.get_aac_base(u_ac, pos)
            pos = int(aachange[1:])

            if proteins.is_position_stored(u_ac, pos):
                continue

            m = proteins.get_position_database_id(u_ac, pos)

            position = proteins.get_position(u_ac, pos)
            mappings = position.mappings

            values.append((
                            m, pack((mappings.recommended_res, mappings.max_seq_res)),
                            pack((
                                mappings.weighted_location, mappings.weighted_mainchain_location,
                                mappings.weighted_sidechain_location, mappings.weighted_surface_value,
                                mappings.weighted_mainchain_surface_value, mappings.weighted_sidechain_surface_value,
                                mappings.Class, mappings.rin_class, mappings.simple_class, mappings.rin_simple_class,
                                str(mappings.interaction_recommendations), mappings.classification_conf, mappings.weighted_ssa, mappings.amount_of_structures, mappings.get_weighted_profile_str(),
                                mappings.weighted_modres, mappings.weighted_b_factor, mappings.get_weighted_centralities_str(),
                                mappings.weighted_phi, mappings.weighted_psi, mappings.weighted_intra_ssbond, mappings.weighted_inter_ssbond,
                                mappings.weighted_intra_link, mappings.weighted_inter_link, mappings.weighted_cis_conformation,
                                mappings.weighted_cis_follower,
                                mappings.weighted_inter_chain_median_kd, mappings.weighted_inter_chain_dist_weighted_kd,
                                mappings.weighted_inter_chain_median_rsa, mappings.weighted_inter_chain_dist_weighted_rsa,
                                mappings.weighted_intra_chain_median_kd, mappings.weighted_intra_chain_dist_weighted_kd,
                                mappings.weighted_intra_chain_median_rsa, mappings.weighted_intra_chain_dist_weighted_rsa,
                                mappings.weighted_inter_chain_interactions_median, mappings.weighted_inter_chain_interactions_dist_weighted,
                                mappings.weighted_intra_chain_interactions_median, mappings.weighted_intra_chain_interactions_dist_weighted,
                                mappings.rsa_change_score, mappings.mc_rsa_change_score, mappings.sc_rsa_change_score,
                                mappings.Backbone_RMSD_seq_id_high_weight, mappings.All_atom_RMSD_seq_id_high_weight, mappings.nof_site_residue_seq_id_high_weight, mappings.Site_LDDT_seq_id_high_weight,
                                mappings.Backbone_RMSD_seq_id_low_weight, mappings.All_atom_RMSD_seq_id_low_weight, mappings.nof_site_residue_seq_id_low_weight, mappings.Site_LDDT_seq_id_low_weight,
                                mappings.Backbone_RMSD_seq_id_greater_90, mappings.All_atom_RMSD_seq_id_greater_90, mappings.nof_site_residue_seq_id_greater_90, mappings.Site_LDDT_seq_id_greater_90,
                                mappings.Backbone_RMSD_seq_id_between_50_and_90, mappings.All_atom_RMSD_seq_id_between_50_and_90, mappings.nof_site_residue_seq_id_between_50_and_90, mappings.Site_LDDT_seq_id_between_50_and_90,
                                mappings.Backbone_RMSD_seq_id_lower_50, mappings.All_atom_RMSD_seq_id_lower_50, mappings.nof_site_residue_seq_id_lower_50, mappings.Site_LDDT_seq_id_lower_50,
                                mappings.Backbone_RMSD_site_id_greater_99, mappings.All_atom_RMSD_site_id_greater_99, mappings.nof_site_residue_site_id_greater_99, mappings.Site_LDDT_site_id_greater_99,
                                mappings.Backbone_RMSD_site_id_between_70_and_99, mappings.All_atom_RMSD_site_id_between_70_and_99, mappings.nof_site_residue_site_id_between_70_and_99, mappings.Site_LDDT_site_id_between_70_and_99,
                                mappings.Backbone_RMSD_site_id_lower_70, mappings.All_atom_RMSD_site_id_lower_70, mappings.nof_site_residue_site_id_lower_70, mappings.Site_LDDT_site_id_lower_70,
                                mappings.Backbone_RMSD_seq_id_greater_90_site_id_greater_99, mappings.All_atom_RMSD_seq_id_greater_90_site_id_greater_99, mappings.nof_site_residue_seq_id_greater_90_site_id_greater_99, mappings.Site_LDDT_seq_id_greater_90_site_id_greater_99,
                                mappings.Backbone_RMSD_seq_id_between_50_and_90_site_id_greater_99, mappings.All_atom_RMSD_seq_id_between_50_and_90_site_id_greater_99, mappings.nof_site_residue_seq_id_between_50_and_90_site_id_greater_99, mappings.Site_LDDT_seq_id_between_50_and_90_site_id_greater_99,
                                mappings.Backbone_RMSD_seq_id_lower_50_site_id_greater_99, mappings.All_atom_RMSD_seq_id_lower_50_site_id_greater_99, mappings.nof_site_residue_seq_id_lower_50_site_id_greater_99, mappings.Site_LDDT_seq_id_lower_50_site_id_greater_99,
                                mappings.Backbone_RMSD_seq_id_greater_90_site_id_between_70_and_99, mappings.All_atom_RMSD_seq_id_greater_90_site_id_between_70_and_99, mappings.nof_site_residue_seq_id_greater_90_site_id_between_70_and_99, mappings.Site_LDDT_seq_id_greater_90_site_id_between_70_and_99,
                                mappings.Backbone_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99, mappings.All_atom_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99, mappings.nof_site_residue_seq_id_between_50_and_90_site_id_between_70_and_99, mappings.Site_LDDT_seq_id_between_50_and_90_site_id_between_70_and_99,
                                mappings.Backbone_RMSD_seq_id_lower_50_site_id_between_70_and_99, mappings.All_atom_RMSD_seq_id_lower_50_site_id_between_70_and_99, mappings.nof_site_residue_seq_id_lower_50_site_id_between_70_and_99, mappings.Site_LDDT_seq_id_lower_50_site_id_between_70_and_99,
                                mappings.Backbone_RMSD_seq_id_greater_90_site_id_lower_70, mappings.All_atom_RMSD_seq_id_greater_90_site_id_lower_70, mappings.nof_site_residue_seq_id_greater_90_site_id_lower_70, mappings.Site_LDDT_seq_id_greater_90_site_id_lower_70,
                                mappings.Backbone_RMSD_seq_id_between_50_and_90_site_id_lower_70, mappings.All_atom_RMSD_seq_id_between_50_and_90_site_id_lower_70, mappings.nof_site_residue_seq_id_between_50_and_90_site_id_lower_70, mappings.Site_LDDT_seq_id_between_50_and_90_site_id_lower_70,
                                mappings.Backbone_RMSD_seq_id_lower_50_site_id_lower_70, mappings.All_atom_RMSD_seq_id_lower_50_site_id_lower_70, mappings.nof_site_residue_seq_id_lower_50_site_id_lower_70, mappings.Site_LDDT_seq_id_lower_50_site_id_lower_70
                            ))
                        ))
            #print(f'{u_ac} {pos} {[mappings.weighted_location, mappings.weighted_mainchain_location, mappings.weighted_sidechain_location, mappings.weighted_surface_value, mappings.weighted_mainchain_surface_value, mappings.weighted_sidechain_surface_value]}')

    return values


def insertComplexes(proteins, config):
    complex_list = proteins.get_complex_list()

    if len(complex_list) == 0:
        return

    smiles_path = config.smiles_path
    inchi_path = config.inchi_path
    pdb_path = config.pdb_path

    stored_ligands = {}
    results = select(config, ['Ligand_Id', 'Name'], 'Ligand')
    for row in results:
        ligand_name = row[1]
        stored_ligands[ligand_name] = row[0]

    values = []

    ligand_db = parseLigandDB(smiles_path, inchi_path)

    lig_values = []
    update_ligands = []
    new_ligands = set()

    for pdb_id in complex_list:
        if proteins.is_complex_stored(pdb_id):
            continue

        resolution = proteins.get_resolution(pdb_id)
        chains_str = proteins.get_chains_str(pdb_id)
        homomers_str = proteins.get_homomers_str(pdb_id)
        lig_str = proteins.get_lig_str(pdb_id)
        metal_str = proteins.get_metal_str(pdb_id)
        ion_str = proteins.get_ion_str(pdb_id)
        cc_str = proteins.get_chain_chain_str(pdb_id)

        values.append((pdb_id, resolution, chains_str, homomers_str, lig_str, metal_str, ion_str, cc_str))

        interaction_partners = proteins.get_interaction_partners(pdb_id)

        for iap in interaction_partners:
            ia_type = iap[0]
            if ia_type != "Ligand":
                continue
            name = iap[1]
            if name in stored_ligands:
                continue
            if name == "UNK" or name == "UNX":
                continue
            if name in new_ligands:
                continue
            res = iap[2]
            chain = iap[3]
            if name in ligand_db:
                (smiles, inchi) = ligand_db[name]
            elif len(name) < 3:
                if len(name) == 1:
                    smiles = "[%s]" % name
                else:
                    smiles = "[%s%s]" % (name[0], name[1].lower())
                inchi = "InChI=1S/%s" % smiles
            else:
                (smiles, inchi) = getSI(pdb_id, name, res, chain, pdb_path, config)
                update_ligands.append((name, smiles, inchi))
            new_ligands.add(name)
            lig_values.append((name, smiles, inchi))

    if len(update_ligands) > 0:
        updateLigandDB(update_ligands, smiles_path, inchi_path)

    if len(lig_values) > 0:
        try:
            insert('Ligand', ['Name', 'Smiles', 'Inchi'], lig_values, config)
        except:
            # There is a (or more than one) ligand that results in an sql error 1241. We need to find it, to figure out why.
            config.errorlog.add_error('Error trying inserting ligands\n%s' % '\n'.join([str(x) for x in lig_values]))

    if len(values) > 0:
        insert('Complex', ['PDB', 'Resolution', 'Chains', 'Homooligomers', 'Ligand_Profile', 'Metal_Profile', 'Ion_Profile', 'Chain_Chain_Profile'], values, config)

    stored_ligands = {}
    results = select(config, ['Ligand_Id', 'Name'], 'Ligand')
    for row in results:
        ligand_name = row[1]
        stored_ligands[ligand_name] = row[0]

    values = []
    results = select(config, ['Complex_Id', 'PDB'], 'Complex')
    for row in results:
        if not proteins.contains_complex(row[1]):
            continue
        if proteins.is_complex_stored(row[1]):
            continue
        proteins.set_complex_db_id(row[1], row[0])

        interaction_partners = proteins.get_interaction_partners(row[1])

        for iap in interaction_partners:
            ia_type = iap[0]
            if ia_type != "Ligand":
                continue
            name = iap[1]
            if name == "UNK" or name == "UNX" or name not in stored_ligands:
                continue
            lig_id = stored_ligands[name]
            res = iap[2]
            chain = iap[3]
            values.append((lig_id, row[0], chain, res))

    if len(values) > 0:
        insert('RS_Ligand_Structure', ['Ligand', 'Complex', 'Chain', 'Residue'], values, config)

