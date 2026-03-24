import time
import ray
import gzip
import zstd
from structman.base_utils.base_utils import unpack, aggregate_times, print_times, add_to_times, alphafold_model_id_to_db_key
from structman.base_utils.config_class import Config
from structman.lib.sdsc import residue as residue_package
from structman.lib.sdsc import interface as interface_package
from structman.lib.database.database_core_functions import binningSelect, select

@ray.remote(max_calls=1)
def remote_select_package_from_db(id_package, config, columns, pkg_number):
    table = 'Residue'
    results = binningSelect(id_package, columns, table, config)
    return (results, pkg_number)

def process_results_chunk(config: Config, results, stored_ids, retrieve_only_db_ids, proteins, locked, skip_interfaces, db_lock):
    times = []
    t_init = time.time()

    if config.verbosity >= 4:
        config.logger.info(f'process resutls chunk {len(stored_ids)=} {len(results)=} {(proteins is None)=} {locked=} {skip_interfaces=} {retrieve_only_db_ids=}')

    stored_res_db_ids: dict[int, tuple[str, str, str | int]] = {}

    if proteins is None:
        structures: dict[str, dict[str, residue_package.Residue_Map]] = {}
        complexes: dict[str, dict[str, dict[str, interface_package.Interface]]] = {}

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
        
        if not retrieve_only_db_ids:
            
            residue = row[3]
            if proteins is None:
                if pdb_id not in structures:
                    structures[pdb_id] = {}
                if chain not in structures[pdb_id]:
                    structures[pdb_id][chain] = residue_package.Residue_Map()
                structures[pdb_id][chain].add_item(res_id, [row[1], residue])
            else:
                proteins.structures[pdb_id][chain].add_residue(res_id, [res_id, row[1], residue])
            if config.verbosity >= 7:
                config.logger.info(f'Calling structures.add_residue in getStoredResidues: {pdb_id} {chain} {res_id}')
        else:
            if proteins is None:
                if pdb_id not in structures:
                    structures[pdb_id] = {}
                if chain not in structures[pdb_id]:
                    structures[pdb_id][chain] = residue_package.Residue_Map()

                if not structures[pdb_id][chain].contains(res_id):
                    structures[pdb_id][chain].add_item(res_id, [row[1], None])
                else:
                    structures[pdb_id][chain].get_item(res_id)[0] = row[1]
            else:
                if not proteins.structures[pdb_id][chain].residues.contains(res_id):
                    proteins.structures[pdb_id][chain].add_residue(res_id, [res_id, row[1], None])
                else:
                    proteins.structures[pdb_id][chain].residues.get_item(res_id)[1] = row[1]
                
        if not retrieve_only_db_ids:
            stored_res_db_ids[row[1]] = pdb_id, chain, res_id

    if config.verbosity >= 2:
        t_now = add_to_times(times, t_init)

    if not retrieve_only_db_ids and not skip_interfaces:

        rows = ['Residue', 'Interacting_Residue']
        table = 'RS_Residue_Interface'
        results = binningSelect(stored_res_db_ids.keys(), rows, table, config, locked=locked, db_lock=db_lock)

        if config.verbosity >= 3:
            config.logger.info(f'Get Interfaces: {len(stored_res_db_ids)=} {len(results)=}')

        total_subtimes = []

        n_interfaces = 0

        for row in results:
            subtimes = []
            ti = time.time()

            residue_db_id = row[0]
            pdb_id, chain, res = stored_res_db_ids[residue_db_id]
            interacting_residue_strs = unpack(row[1])

            ti = add_to_times(subtimes, ti)

            for interacting_residue_str in interacting_residue_strs.split(';'):
                interacting_chain, i_res, score = interacting_residue_str.split(',')

                if proteins is None:
                    if pdb_id not in complexes:
                        complexes[pdb_id] = {}
                    if chain not in complexes[pdb_id]:
                        n_interfaces += 1
                        complexes[pdb_id][chain] = {}
                        complexes[pdb_id][chain][interacting_chain] = interface_package.Interface(chain, interacting_chain, stored = True)
                    elif interacting_chain not in complexes[pdb_id][chain]:
                        complexes[pdb_id][chain][interacting_chain] = interface_package.Interface(chain, interacting_chain, stored = True)

                    if i_res != '':
                        complexes[pdb_id][chain][interacting_chain].add_interaction(res, i_res, float(score))
                    else:
                        complexes[pdb_id][chain][interacting_chain].add_support(res)

                else:
                    if chain not in proteins.complexes[pdb_id].interfaces:
                        n_interfaces += 1
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

        if config.verbosity >= 3:
            config.logger.info(f'Total number of received interfaces: {n_interfaces}')

        if config.verbosity >= 2:
            times.append(total_subtimes)

    if config.verbosity >= 2:
        t_now = add_to_times(times, t_now)
        
    if proteins is None:
        return times, structures, complexes
    return times

def getStoredResidues(
        proteins,
        config: Config,
        custom_ids = None,
        retrieve_only_db_ids = False,
        exclude_interacting_chains = False,
        para = False,
        locked = False,
        skip_interfaces = False,
        db_lock=None
        ):
    t0 = time.time()

    if custom_ids is None:
        stored_ids: dict[int, tuple[str, str]] = proteins.getStoredStructureIds(exclude_interacting_chains = exclude_interacting_chains)
        if config.verbosity >= 3:
            config.logger.info(f'Call of getStoredResidues, retrieve_only_db_ids {retrieve_only_db_ids}, exclude_interacting_chains {exclude_interacting_chains}: {len(stored_ids)}')
    else:
        stored_ids = custom_ids
        if config.verbosity >= 3:
            config.logger.info(f'Call of getStoredResidues, retrieve_only_db_ids {retrieve_only_db_ids}, exclude_interacting_chains {exclude_interacting_chains} with custom_ids: {len(custom_ids)}')

    if config.verbosity >= 2:
        t1 = time.time()
        config.logger.info(f"Time for getstoredresidues 1: {t1 - t0}")

    if retrieve_only_db_ids:
        max_number_of_residues = 10000000
    else:
        max_number_of_residues = 50000
        
    total_times = []
    total_structures: dict[str, dict[str, residue_package.Residue_Map]] = {}
    total_complexes: dict[str, dict[str, dict[str, interface_package.Interface]]] = {}
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

        remote_selects_process_ids = []
        if para:
            conf_dump = ray.put(config)

        if not retrieve_only_db_ids:
            rows = ['Structure', 'Residue_Id', 'Number', 'Residue_Data']
        else:
            rows = ['Structure', 'Residue_Id', 'Number']

        for pkg_number, id_package in enumerate(id_packages):

            t10 = time.time()
            times = []
            if len(id_packages) > 1 and para:
                remote_selects_process_ids.append(remote_select_package_from_db.remote(id_package, conf_dump, rows, pkg_number))
            else:
                table = 'Residue'
                results = binningSelect(id_package, rows, table, config, locked=locked, db_lock=db_lock)

                t11 = time.time()
                if config.verbosity >= 4:
                    config.logger.info(f"Time for getstoredresidues 2.1: {t11 - t10} {len(id_package)}")
                times.append(t11-t10)

                if proteins is not None:
                    proc_times = process_results_chunk(config, results, stored_ids, retrieve_only_db_ids, proteins, locked, skip_interfaces, db_lock)
                else:
                    proc_times, structures, complexes = process_results_chunk(config, results, stored_ids, retrieve_only_db_ids, None, locked, skip_interfaces, db_lock)
                    total_structures.update(structures)
                    total_complexes.update(complexes)
                times.extend(proc_times)
                
                t12 = time.time()
                if config.verbosity >= 4:
                    config.logger.info(f"Time for getstoredresidues 2.2: {t12 - t11} {len(id_package)}")
                times.append(t12-t11)

                total_times = aggregate_times(total_times, times)

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
                        times = process_results_chunk(config, results, stored_ids, retrieve_only_db_ids, proteins, skip_interfaces)
                        total_times = aggregate_times(total_times, times)

                remote_selects_process_ids = not_ready
                if len(remote_selects_process_ids) == 0 and len(returned_pkgs) == len(id_packages):
                    break

                if count >= 50000:
                    if config.verbosity >= 2:
                        config.logger.info(f'remote getStoredResidues froze: {len(id_packages)=} {len(returned_pkgs)=}')
                    for remote_process in remote_selects_process_ids:
                        ray.kill(remote_process)
                    if len(id_packages) - len(returned_pkgs) == 1:
                        results = binningSelect(id_package, rows, 'Residue', config, locked=locked, db_lock=db_lock)
                        times = process_results_chunk(config, results, stored_ids, retrieve_only_db_ids, proteins, skip_interfaces)
                        total_times = aggregate_times(total_times, times)
                        break
                    remote_selects_process_ids = []
                    for pkg_number, id_package in enumerate(id_packages):
                        if pkg_number not in returned_pkgs:
                            remote_selects_process_ids.append(remote_select_package_from_db.remote(id_package, conf_dump, rows, pkg_number))
                    
        else:
            if config.verbosity >= 2:
                t12 = time.time()
                config.logger.info(f"Time for remote select: {t12 - t11}")

    if config.verbosity >= 2:
        print_times(total_times, label = 'process residue data', logger=config.logger)
        t2 = time.time()
        config.logger.info(f"Time for getstoredresidues 2: {t2 - t1}")

    if proteins is None:
        return total_structures, total_complexes


def get_afdb_data(model_id: str, config: Config, columns_tuple: tuple[str, str]):
    afdb_id, main_db, chunk_id = alphafold_model_id_to_db_key(model_id)
    if main_db:
        subdb = 'main'
    else:
        subdb = 'collab'
    table_name = f'AFDB_{subdb}_{chunk_id}'

    if config.verbosity >= 4:
        config.logger.info(f'get_afdb_data: {model_id=} {columns_tuple=} {afdb_id=} {table_name=}')

    results = select(config, columns_tuple, table_name, equals_rows={'ID':afdb_id}, db_name=config.afdb)
    if len(results) == 0:
        return None

    pages = []
    for i in range(len(columns_tuple)//2):
        try:
            data, compr_type = results[0][(i*2):(i*2)+2]
        except ValueError as err: 
            raise ValueError(f'{err=} {columns_tuple=} {results=}')

        if compr_type == 0:
            page = gzip.decompress(data)
        elif compr_type == 1:
            page = zstd.decompress(data)
        else:
            page = None
        pages.append(page)
    if len(pages) == 1:
        return pages[0]
    return pages

def get_afdb_structure(model_id: str, config: Config):
    columns = ['Structure_data', 'Sd_compr_type']
    return get_afdb_data(model_id, config, columns)

def get_afdb_msa(model_id: str, config: Config):
    columns = ['MSA_data', 'Md_compr_type']
    return get_afdb_data(model_id, config, columns)

def get_afdb_rin_and_cent(model_id: str, config: Config):
    columns = ['RIN_data', 'Rd_compr_type', 'Centralities_data', 'Cd_compr_type']
    return get_afdb_data(model_id, config, columns)