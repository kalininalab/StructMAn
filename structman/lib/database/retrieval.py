import time
import ray
from structman.base_utils.base_utils import unpack, aggregate_times, print_times, add_to_times
from structman.lib.sdsc import residue as residue_package
from structman.lib.sdsc import interface as interface_package
from structman.lib.database.database_core_functions import binningSelect

@ray.remote(max_calls=1)
def remote_select_package_from_db(id_package, config, columns, pkg_number):
    table = 'Residue'
    results = binningSelect(id_package, columns, table, config)
    return (results, pkg_number)

def process_results_chunk(config, results, stored_ids, retrieve_only_db_ids, proteins):
    times = []
    t_init = time.time()

    stored_res_db_ids: dict[int, tuple[str, str, str | int]] = {}

    for row in results:
        try:
            try:
                res_id = int(row[2])
            except:
                res_id = row[2]
            
        except:
            config.errorlog.add_warning('Defective database entry in Residue table: %s %s' % (str(row[0]), str(row[1])))
            continue

        s_id = row[0]
        pdb_id, chain = stored_ids[s_id]
        
        # Those residue inits include decoding of interaction profile and centrality score strings and thus takes some resources. For that a para function
        if not retrieve_only_db_ids:
            
            residue = row[3]
            proteins.structures[pdb_id][chain].add_residue(res_id, (res_id, row[1], residue))
            if config.verbosity >= 6:
                print(f'Calling structures.add_residue in getStoredResidues: {pdb_id} {chain} {res_id}')
        else:
            if not proteins.structures[pdb_id][chain].residues.contains(res_id):
                residue = residue_package.Residue()
                residue.res_num = res_id
                residue.database_id=row[1]
                residue.stored=True
                proteins.structures[pdb_id][chain].add_residue(res_id, residue)
            else:
                proteins.structures[pdb_id][chain].residues.get_item(res_id).database_id =  row[1]
                proteins.structures[pdb_id][chain].residues.get_item(res_id).stored = True

        if not retrieve_only_db_ids:
            stored_res_db_ids[row[1]] = pdb_id, chain, res_id

    if config.verbosity >= 2:
        t_now = add_to_times(times, t_init)

    if not retrieve_only_db_ids:

        rows = ['Residue', 'Interacting_Residue']
        table = 'RS_Residue_Interface'
        results = binningSelect(stored_res_db_ids.keys(), rows, table, config)

        if config.verbosity >= 3:
            print(f'Get Interfaces: {len(stored_res_db_ids)=} {len(results)=}')

        total_subtimes = []

        for row in results:
            subtimes = []
            ti = time.time()

            residue_db_id = row[0]
            pdb_id, chain, res = stored_res_db_ids[residue_db_id]
            interacting_residue_strs = unpack(row[1])

            ti = add_to_times(subtimes, ti)

            for interacting_residue_str in interacting_residue_strs.split(';'):
                interacting_chain, i_res, score = interacting_residue_str.split(',')

                if not chain in proteins.complexes[pdb_id].interfaces:
                    proteins.complexes[pdb_id].interfaces[chain] = {}
                    proteins.complexes[pdb_id].interfaces[chain][interacting_chain] = interface_package.Interface(chain, interacting_chain, stored = True)
                elif interacting_chain not in proteins.complexes[pdb_id].interfaces[chain]:
                    proteins.complexes[pdb_id].interfaces[chain][interacting_chain] = interface_package.Interface(chain, interacting_chain, stored = True)

                if i_res != '':
                    proteins.complexes[pdb_id].interfaces[chain][interacting_chain].add_interaction(res, i_res, float(score))
                else:
                    proteins.complexes[pdb_id].interfaces[chain][interacting_chain].add_support(res)

            add_to_times(subtimes, ti)
            aggregate_times(total_subtimes, subtimes)

        if config.verbosity >= 2:
            times.append(total_subtimes)

    if config.verbosity >= 2:
        t_now = add_to_times(times, t_now)
        
    return times

def getStoredResidues(
        proteins, config,
        custom_ids = None,
        retrieve_only_db_ids = False,
        exclude_interacting_chains = False,
        para = False
        ):
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

    if retrieve_only_db_ids:
        max_number_of_residues = 10000000
    else:
        max_number_of_residues = 50000
        
    total_times = []
    if len(stored_ids) > 0:
        if len(stored_ids) >= max_number_of_residues:
            number_of_packages = (len(stored_ids)//max_number_of_residues) + 1
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

        remote_selects_process_ids = []
        conf_dump = ray.put(config)

        if not retrieve_only_db_ids:
            rows = ['Structure', 'Residue_Id', 'Number', 'Residue_Data']
        else:
            rows = ['Structure', 'Residue_Id', 'Number']

        for pkg_number, id_package in enumerate(id_packages):
            if config.verbosity >= 2:
                t10 = time.time()
            
            if len(id_packages) > 1 and para:
                remote_selects_process_ids.append(remote_select_package_from_db.remote(id_package, conf_dump, rows, pkg_number))
            else:
                table = 'Residue'
                results = binningSelect(id_package, rows, table, config)

            if config.verbosity >= 2:
                t11 = time.time()
                print(f"Time for getstoredresidues 2.1: {t11 - t10} {len(id_package)}")

        if config.verbosity >= 2:
            t11 = time.time()
    
        if len(id_packages) > 1 and para:
            count = 0
            returned_pkgs = set()
            while True:
                ready, not_ready = ray.wait(remote_selects_process_ids, timeout = 0.1)
                count += 1
                if len(ready) > 0:
                    count = 0
                    result_chunks = ray.get(ready)
                    for results, pkg_number in result_chunks:
                        returned_pkgs.add(pkg_number)
                        times = process_results_chunk(config, results, stored_ids, retrieve_only_db_ids, proteins)
                        total_times = aggregate_times(total_times, times)

                remote_selects_process_ids = not_ready
                if len(remote_selects_process_ids) == 0 and len(returned_pkgs) == len(id_packages):
                    break

                if count >= 50000:
                    if config.verbosity >= 2:
                        print(f'remote getStoredResidues froze: {len(id_packages)=} {len(returned_pkgs)=}')
                    for remote_process in remote_selects_process_ids:
                        ray.kill(remote_process)
                    if len(id_packages) - len(returned_pkgs) == 1:
                        results = binningSelect(id_package, rows, 'Residue', config)
                        times = process_results_chunk(config, results, stored_ids, retrieve_only_db_ids, proteins)
                        total_times = aggregate_times(total_times, times)
                        break
                    remote_selects_process_ids = []
                    for pkg_number, id_package in enumerate(id_packages):
                        if pkg_number not in returned_pkgs:
                            remote_selects_process_ids.append(remote_select_package_from_db.remote(id_package, conf_dump, rows, pkg_number))
                    
        else:
            if config.verbosity >= 2:
                t12 = time.time()
                print(f"Time for remote select: {t12 - t11}")

            times = process_results_chunk(config, results, stored_ids, retrieve_only_db_ids, proteins)
            total_times = aggregate_times(total_times, times)

    if config.verbosity >= 2:
        print_times(total_times, label = 'process residue data')
        t2 = time.time()
        print("Time for getstoredresidues 2: %s" % str(t2 - t1))