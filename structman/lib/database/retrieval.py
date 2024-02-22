import time

from structman.base_utils.base_utils import unpack
from structman.lib.sdsc import residue as residue_package
from structman.lib.sdsc import interface as interface_package
from structman.lib.database.database_core_functions import binningSelect

def getStoredResidues(proteins, config, custom_ids = None, retrieve_only_db_ids = False, exclude_interacting_chains = False):
    t0 = time.time()

    if custom_ids is None:
        stored_ids = proteins.getStoredStructureIds(exclude_interacting_chains = exclude_interacting_chains)
        if config.verbosity >= 3:
            print(f'Call of getStoredResidues, retrieve_only_db_ids {retrieve_only_db_ids}, exclude_interacting_chains {exclude_interacting_chains}: {len(stored_ids)}')
    else:
        stored_ids = custom_ids
        if config.verbosity >= 3:
            print(f'Call of getStoredResidues, retrieve_only_db_ids {retrieve_only_db_ids}, exclude_interacting_chains {exclude_interacting_chains} with custom_ids: {len(custom_ids)}')

    if config.verbosity >= 2:
        t1 = time.time()
        print("Time for getstoredresidues 1: %s" % str(t1 - t0))

    max_number_of_structures = 30000

    if len(stored_ids) > 0:
        if len(stored_ids) >= max_number_of_structures:
            number_of_packages = (len(stored_ids)//max_number_of_structures) + 1
            package_size = len(stored_ids) // number_of_packages
            if len(stored_ids) % number_of_packages != 0:
                package_size += 1
            id_packages = []
            complete_id_list = list(stored_ids.keys())
            for package_number in range(number_of_packages):
                left = package_number * package_size
                right = (package_number+1) * package_size
                id_packages.append(complete_id_list[left:right])
        else:
            id_packages = [list(stored_ids.keys())]

        """
        rows = ['Structure', 'Residue_Id', 'Number', 'Amino_Acid', 'Sub_Lig_Dist', 'Sub_Chain_Distances',
                'Relative_Surface_Access', 'Relative_Surface_Access_Main_Chain', 'Relative_Surface_Access_Side_Chain',
                'Secondary_Structure_Assignment', 'Homomer_Distances',
                'Interaction_Profile', 'Centralities', 'Modres', 'B_Factor', 'PHI', 'PSI', 'Intra_SSBOND', 'SSBOND_Length',
                'Intra_Link', 'Link_Length', 'CIS_Conformation', 'CIS_Follower', 'Inter_Chain_Median_KD',
                'Inter_Chain_Dist_Weighted_KD', 'Inter_Chain_Median_RSA',
                'Inter_Chain_Dist_Weighted_RSA', 'Intra_Chain_Median_KD', 'Intra_Chain_Dist_Weighted_KD',
                'Intra_Chain_Median_RSA', 'Intra_Chain_Dist_Weighted_RSA',
                'Inter_Chain_Interactions_Median', 'Inter_Chain_Interactions_Dist_Weighted',
                'Intra_Chain_Interactions_Median', 'Intra_Chain_Interactions_Dist_Weighted',
                'Interacting_Chains', 'Interacting_Ligands'
                ]
        """
        for id_package in id_packages:
            if not retrieve_only_db_ids:
                rows = ['Structure', 'Residue_Id', 'Number', 'Residue_Data']
            else:
                rows = ['Structure', 'Residue_Id', 'Number']
            table = 'Residue'
            results = binningSelect(id_package, rows, table, config)

            if config.verbosity >= 2:
                t10 = time.time()
                print("Time for getstoredresidues 2.1: %s" % str(t10 - t1))

            stored_res_db_ids = {}

            for row in results:
                try:
                    res_id = row[2]
                    if not retrieve_only_db_ids:
                        (one_letter, lig_dist_str, chain_dist_str, rsa, relative_main_chain_acc, relative_side_chain_acc,
                               ssa, homo_str, profile_str,
                               centrality_score_str, b_factor, modres, phi, psi, intra_ssbond, ssbond_length, intra_link, link_length,
                               cis_conformation, cis_follower, inter_chain_median_kd, inter_chain_dist_weighted_kd,
                               inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                               intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                               inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                               intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                               interacting_chains_str, interacting_ligands_str) = unpack(row[3])
                except:
                    config.errorlog.add_warning('Defective database entry in Residue table: %s %s' % (str(row[0]), str(row[1])))
                    continue

                s_id = row[0]
                pdb_id, chain = stored_ids[s_id]
                
                # Those residue inits include decoding of interaction profile and centrality score strings and thus takes some resources. For that a para function
                if not retrieve_only_db_ids:
                    residue = residue_package.Residue(res_id, aa=one_letter, lig_dist_str=lig_dist_str, chain_dist_str=chain_dist_str, RSA=rsa,
                                       relative_main_chain_acc=relative_main_chain_acc, relative_side_chain_acc=relative_side_chain_acc,
                                       SSA=ssa, homo_dist_str=homo_str, interaction_profile_str=profile_str, centrality_score_str=centrality_score_str,
                                       modres=modres, b_factor=b_factor, database_id=row[1], stored=True, phi=phi, psi=psi,
                                       intra_ssbond=intra_ssbond, ssbond_length=ssbond_length, intra_link=intra_link, link_length=link_length,
                                       cis_conformation=cis_conformation, cis_follower=cis_follower, inter_chain_median_kd=inter_chain_median_kd,
                                       inter_chain_dist_weighted_kd=inter_chain_dist_weighted_kd, inter_chain_median_rsa=inter_chain_median_rsa,
                                       inter_chain_dist_weighted_rsa=inter_chain_dist_weighted_rsa, intra_chain_median_kd=intra_chain_median_kd,
                                       intra_chain_dist_weighted_kd=intra_chain_dist_weighted_kd, intra_chain_median_rsa=intra_chain_median_rsa,
                                       intra_chain_dist_weighted_rsa=intra_chain_dist_weighted_rsa,
                                       inter_chain_interactions_median=inter_chain_interactions_median, inter_chain_interactions_dist_weighted=inter_chain_interactions_dist_weighted,
                                       intra_chain_interactions_median=intra_chain_interactions_median, intra_chain_interactions_dist_weighted=intra_chain_interactions_dist_weighted,
                                       interacting_chains_str=interacting_chains_str, interacting_ligands_str=interacting_ligands_str)
                    proteins.add_residue(pdb_id, chain, res_id, residue)
                    if config.verbosity >= 6:
                        print(f'Calling proteins.add_residue in getStoredResidues: {pdb_id} {chain} {res_id}')
                else:
                    if res_id not in proteins.structures[(pdb_id, chain)].residues:
                        residue = residue_package.Residue(res_id, database_id=row[1], stored=True)
                        proteins.add_residue(pdb_id, chain, res_id, residue)
                    else:
                        proteins.structures[(pdb_id, chain)].residues[res_id].database_id =  row[1]
                        proteins.structures[(pdb_id, chain)].residues[res_id].stored = True

                stored_res_db_ids[row[1]] = pdb_id, chain, res_id

            rows = ['Residue', 'Interacting_Residue']
            table = 'RS_Residue_Interface'
            results = binningSelect(stored_res_db_ids, rows, table, config)


            for row in results:
                residue_db_id = row[0]
                pdb_id, chain, res = stored_res_db_ids[residue_db_id]
                interacting_chain, i_res, score = row[1].split(',')

                if not (chain, interacting_chain) in proteins.complexes[pdb_id].interfaces:
                    proteins.complexes[pdb_id].interfaces[(chain, interacting_chain)] = interface_package.Interface(chain, interacting_chain, stored = True)
                if i_res != '':
                    proteins.complexes[pdb_id].interfaces[(chain, interacting_chain)].add_interaction(res, i_res, float(score))
                else:
                    proteins.complexes[pdb_id].interfaces[(chain, interacting_chain)].add_support(res)

    if config.verbosity >= 2:
        t2 = time.time()
        print("Time for getstoredresidues 2: %s" % str(t2 - t1))