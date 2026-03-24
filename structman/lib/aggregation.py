import ray
import time
import psutil
import gc
import math
import sys
from itertools import cycle

from ray.util.queue import Queue

from structman.scripts.createMMDB import mm_lookup
from structman.base_utils.base_utils import pack, unpack, aggregate_times, print_times, add_to_times
from structman.base_utils.config_class import Config
from structman.lib.sdsc import mappings as mappings_package
from structman.lib.sdsc import interface as interface_package
from structman.lib.sdsc import residue as residue_package
from structman.lib.sdsc import protein as protein_package
from structman.lib.sdsc.sdsc_utils import get_shortest_distances, triple_locate, classify
from structman.lib.database.insertion_lib import remote_insert_interfaces, remote_insert_classifications
from structman.lib.database.retrieval import getStoredResidues
from structman.lib.database.retrieval_2 import getBackmaps

def qualityScore(resolution: float, coverage: float, seq_id: float, resolution_wf=0.25, coverage_wf=0.5, seq_id_wf=1.0, pLDDT:float | None = None) -> float:

    if seq_id > 1.0:
        seq_id = seq_id / 100.0
    if coverage > 1.0:
        coverage = coverage / 100.0

    if pLDDT is None:
        pLDDT = 1.0

    # project the criteria to [0,1] via logistic regression
    resolution_value = (1 + math.exp((1.5 * resolution) - 4))**(-1)

    seq_value = (1 + math.exp(10 * (0.4 - seq_id)))**(-1)

    ws = sum((resolution_wf, coverage_wf, seq_id_wf))
    quality = pLDDT * (sum((resolution_value * resolution_wf, coverage * coverage_wf, seq_value * seq_id_wf)) / ws)
    return quality

@ray.remote(max_retries=100, max_calls = 1)
def para_classify_remote_wrapper(config: Config, complex_store, package, structure_db_dict, prot_db_id_map, packed_structures, backmap_store, com_queue, db_lock):
    t0 = time.time()
    #structs_from_db, interface_map = unpack(packed_structures)
    #del packed_structures
    complex_store = unpack(complex_store)
    structure_db_dict = unpack(structure_db_dict)
    package = unpack(package)
    t1 = time.time()

    return (para_classify(config,
                          complex_store,
                          package,
                          structure_db_dict,
                          prot_db_id_map,
                          stored_structures = packed_structures,
                          backmap_store = backmap_store,
                          com_queue = com_queue,
                          db_lock=db_lock, para=True, unpacking_time=(t1-t0)))


def para_classify(
        config: Config,
        complexes: dict[str, tuple[dict[str, str], float, int]],
        package: list[tuple[str, protein_package.Protein]],
        structure_db_dict: tuple[dict[int, tuple[str, str]], dict[int, tuple[str, str]]],
        prot_db_id_map: dict[int, list[int]],
        stored_structures = None,
        backmap_store = None,
        com_queue = None,
        db_lock= None, para=False, unpacking_time = 0.):
    ta = time.time()
    times = [unpacking_time]

    if config.verbosity >= 4:
        config.logger.info(f'Start of para_classify {len(package)=} {para=}')

    if backmap_store is None:
        backmap_store = {}

    structs_from_db: dict[str, dict[str, residue_package.Residue_Map[tuple[int, bytes]]]] = {}
    interface_map: dict[str, dict[str, dict[str, interface_package.Interface]]] = {}

    structs_from_db, interface_map = getStoredResidues(None, config,
        custom_ids = structure_db_dict[0],
        exclude_interacting_chains = True,
        locked = True,
        db_lock=db_lock)
    
    ta = add_to_times(times, ta)

    if stored_structures is not None:
        packed_structures, packed_interfaces = stored_structures
        for structure_id in packed_structures:
            if structure_id not in structs_from_db:
                structs_from_db[structure_id] = {}
            
            for chain in packed_structures[structure_id]:
                structs_from_db[structure_id][chain] = packed_structures[structure_id][chain]

        for structure_id in packed_interfaces:
            if structure_id not in interface_map:
                interface_map[structure_id] = packed_interfaces[structure_id]

    ta = add_to_times(times, ta)
        
    more_structs_from_db, _ = getStoredResidues(None, config,
            custom_ids = structure_db_dict[1],
            retrieve_only_db_ids=True,
            exclude_interacting_chains = True,
            skip_interfaces = True,
            locked = True,
            db_lock=db_lock)
    
    ta = add_to_times(times, ta)

    complex_to_db_ids: dict[str, dict[int, tuple[str, str]]] = {}
    structure_to_db_id: dict[str, dict[str, int]] = {}

    for db_id in structure_db_dict[0]:
        pdb_id, chain = structure_db_dict[0][db_id]
        if pdb_id not in complex_to_db_ids:
            complex_to_db_ids[pdb_id] = {}
        complex_to_db_ids[pdb_id][db_id] = (pdb_id, chain)
        if pdb_id not in structure_to_db_id:
            structure_to_db_id[pdb_id] = {}
        structure_to_db_id[pdb_id][chain] = db_id

    for db_id in structure_db_dict[1]:
        pdb_id, chain = structure_db_dict[1][db_id]
        if pdb_id not in complex_to_db_ids:
            complex_to_db_ids[pdb_id] = {}
        complex_to_db_ids[pdb_id][db_id] = (pdb_id, chain)
        if pdb_id not in structure_to_db_id:
            structure_to_db_id[pdb_id] = {}
        structure_to_db_id[pdb_id][chain] = db_id

    ta = add_to_times(times, ta)

    for pdb_id in more_structs_from_db:
        
        if pdb_id not in structs_from_db:
            structs_from_db[pdb_id] = more_structs_from_db[pdb_id]
        else:
            for chain in more_structs_from_db[pdb_id]:
                if chain not in structs_from_db[pdb_id]:
                    structs_from_db[pdb_id][chain] = more_structs_from_db[pdb_id][chain]

    del more_structs_from_db

    ta = add_to_times(times, ta)
    total_sub_times = []
    position_values = []

    for packnum, (protein_id, protein_obj) in enumerate(package):
        ts = time.time()
        sub_times = []
        annotation_list = protein_obj.get_annotation_list()

        position_value_dict = {}
        
        recommended_complexes: list[list[tuple[str, str, int | str]]] = [None]

        prot_db_id: int = protein_obj.database_id

        positions: list[int] = protein_obj.get_position_ids()

        structure_quality_measures: dict[str, dict[str, tuple[float, float]]] = {}

        every_second_time = True
        if len(annotation_list) > 1000:
            #flush_struct_store = True
            flush_struct_store = False

        else:
            flush_struct_store = False

        if config.verbosity >= 4:
            config.logger.info(f"In para_classify: {protein_id=} {len(annotation_list)=} {len(positions)=}")

        for pos in positions:
            
            if protein_obj.is_position_stored(pos):
                continue

            mappings_obj: mappings_package.Mappings = mappings_package.Mappings()
            aacbase = protein_obj.get_aac_base(pos)            

            aa1 = aacbase[0]

            for (pdb_id, chain) in annotation_list:

                try:
                    sub_info = protein_obj.structure_annotations[pdb_id][chain].get_sub_info(pos)

                except:
                    if config.verbosity >= 4:
                        config.logger.info(f'Skipped classification of {protein_id} due to pos {pos} was not in sub_infos')
                    continue

                if sub_info is None: #In this case the position is not mapped to any residue in the structure (gap in the alignment and getSubInfo called with ignore_gaps) 
                    continue

                res_nr = sub_info[0]

                if res_nr is None: #In this case the position is not mapped to any residue in the structure (gap in the alignment)
                    continue

                seq_id = protein_obj.structure_annotations[pdb_id][chain].get_sequence_id()

                if seq_id is None:
                    if config.verbosity >= 4:
                        config.logger.info(f'Skipped classification of {protein_id} due to {pdb_id} got no sequence identity')
                    continue

                cov = protein_obj.structure_annotations[pdb_id][chain].coverage                

                if pdb_id not in structure_quality_measures:
                    structure_quality_measures[pdb_id] = {}
                if chain not in structure_quality_measures[pdb_id]:
                    structure_quality_measures[pdb_id][chain] = (seq_id, cov)

                try:
                    (chains, resolution, _) = complexes[pdb_id]
                except:
                    (chains, resolution) = (None, None)

                """
                if pdb_id not in structs_from_db:
                    structs_from_db[pdb_id] = {}
                    interface_map[pdb_id] = {}
                if chain not in structs_from_db[pdb_id]:
                    structure_db_id = structure_to_db_id[pdb_id][chain]
                    custom_id = {structure_db_id: (pdb_id, chain)}
                    struct_from_db, struct_interface_map = getStoredResidues(None, config,
                        custom_ids = custom_id,
                        exclude_interacting_chains = True,
                        locked = True,
                        db_lock=db_lock)
                    structs_from_db[pdb_id][chain] = struct_from_db[pdb_id][chain]
                    try:
                        interface_map[pdb_id][chain] = struct_interface_map[pdb_id][chain]
                    except KeyError:
                        pass

                """
                if pdb_id not in structs_from_db:
                    continue
                if chain not in structs_from_db[pdb_id]:
                    continue

                if not flush_struct_store and every_second_time:
                    
                    if not isinstance(structs_from_db[pdb_id][chain], residue_package.Residue_Map):
                        structs_from_db[pdb_id][chain] = unpack(ray.get(structs_from_db[pdb_id][chain]))
                        every_second_time = False

                    if not structs_from_db[pdb_id][chain].contains(res_nr):
                        continue

                    _, packed_res_info = structs_from_db[pdb_id][chain].get_item(res_nr)
                    mapping = unpack(packed_res_info)
                    
                else:
                    
                    if not isinstance(structs_from_db[pdb_id][chain], residue_package.Residue_Map):
                        res_map = unpack(ray.get(structs_from_db[pdb_id][chain]))
                        every_second_time = True
                    else:
                        res_map = structs_from_db[pdb_id][chain]

                    if not res_map.contains(res_nr):
                        continue

                    _, packed_res_info = res_map.get_item(res_nr)
                    mapping = unpack(packed_res_info)

                if mapping is None:
                    continue

                (rsa, mc_rsa, sc_rsa, ssa, profile, centralities,
                 phi, psi, intra_ssbond, inter_ssbond, ssbond_length, intra_link, inter_link, link_length, cis_conformation, cis_follower,
                 inter_chain_median_kd, inter_chain_dist_weighted_kd,
                 inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                 intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                 inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                 intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                 b_factor, pLDDT, modres,
                 lig_dists, chain_distances, homomer_distances, res_aa) = mapping

                if res_aa is None:
                    continue

                identical_aa = res_aa == aa1

                gsd_return = get_shortest_distances(chains, lig_dists, chain_distances, homomer_distances)

                if gsd_return is None:
                    continue
                homo_dist, lig_dist, metal_dist, ion_dist, chain_dist, rna_dist, dna_dist, min_lig, min_metal, min_ion, iacs = gsd_return

                loc, mc_loc, sc_loc = triple_locate(rsa, mc_rsa, sc_rsa, config)

                rin_class, rin_simple_class = classify(profile, sc_loc)

                try:
                    qual: float = qualityScore(resolution, cov, seq_id, pLDDT = pLDDT)
                except TypeError:
                    qual = 0.01
                    resolution = 3.0

                structural_feature_dict = {
                    'b_factor' : b_factor, 'pLDDT': pLDDT, 'modres' : modres, 'ssa' : ssa, 'phi' : phi, 'psi' : psi, 'intra_ssbond' : intra_ssbond, 'inter_ssbond' : inter_ssbond, 'ssbond_length' : ssbond_length,
                    'intra_link' : intra_link, 'inter_link' : inter_link, 'link_length' : link_length, 'cis_conformation' : cis_conformation, 'cis_follower' : cis_follower,
                    'inter_chain_median_kd' : inter_chain_median_kd, 'inter_chain_dist_weighted_kd' : inter_chain_dist_weighted_kd,
                    'inter_chain_median_rsa' : inter_chain_median_rsa, 'inter_chain_dist_weighted_rsa' : inter_chain_dist_weighted_rsa,
                    'intra_chain_median_kd' : intra_chain_median_kd, 'intra_chain_dist_weighted_kd' : intra_chain_dist_weighted_kd,
                    'intra_chain_median_rsa' : intra_chain_median_rsa, 'intra_chain_dist_weighted_rsa' : intra_chain_dist_weighted_rsa,
                    'inter_chain_interactions_median' : inter_chain_interactions_median, 'inter_chain_interactions_dist_weighted' : inter_chain_interactions_dist_weighted,
                    'intra_chain_interactions_median' : intra_chain_interactions_median, 'intra_chain_interactions_dist_weighted' : intra_chain_interactions_dist_weighted,
                    'lig_dist' : lig_dist, 'metal_dist' : metal_dist, 'ion_dist' : ion_dist, 'chain_dist' : chain_dist, 'rna_dist' : rna_dist, 'dna_dist' : dna_dist, 'homo_dist' : homo_dist,
                    'surface_value' : rsa, 'mainchain_surface_value' : mc_rsa, 'sidechain_surface_value' : sc_rsa,
                    'rin_class' : rin_class, 'rin_simple_class' : rin_simple_class
                }

                rin_based_feature_dict = {
                    'profile' : profile,
                    'centralities' : centralities
                }

                mapping = (qual, seq_id, cov, structural_feature_dict, rin_based_feature_dict, identical_aa, resolution, res_aa)

                mappings_obj.add_mapping((pdb_id, chain, res_nr), mapping)
         	

            if config.verbosity >= 4:
                config.logger.info(f"In para_classify: {protein_id=} before weighting of {pos=}")

            err = mappings_obj.weight_all(config)
            
            if err is not None and config.verbosity >= 4:
                config.logger.info(f'{protein_id} {pos}: {err}')

            if config.verbosity >= 4:
                config.logger.info(f"In para_classify: {protein_id=} after weighting of {pos=}")

            recommendation_order: list[tuple[str, str, int | str]] = mappings_obj.recommendation_order
            mappings_obj.recommendation_order = None

            try:
                rec_res = mappings_obj.recommended_res
                max_res = mappings_obj.max_seq_res
            except AttributeError:
                rec_res = None
                max_res = None

            packed_recs = pack((rec_res, max_res))

            position_value_dict[pos] = [prot_db_id_map[prot_db_id][pos], packed_recs, pack(mappings_obj) ,None]
            del mappings_obj

            recommended_complexes.append(recommendation_order)

        ts = add_to_times(sub_times, ts)

        if config.verbosity >= 4:
            config.logger.info(f"In para_classify: {protein_id=} before calc_mm_features")

        mm_db_dict: dict[str, dict[str, dict[str | int, bytes]]]
        mm_db_dict, successful_mm_annotations = calc_mm_features(config, recommended_complexes)

        if config.verbosity >= 4:
            config.logger.info(f"In para_classify: {protein_id=} after calc_mm_features {len(mm_db_dict)=} {len(successful_mm_annotations)=}")

        ts = add_to_times(sub_times, ts)

        del recommended_complexes

        for pos in successful_mm_annotations:
            complex_id, chain_id, res = successful_mm_annotations[pos]
            if chain_id in mm_db_dict[complex_id]:
                if res in mm_db_dict[complex_id][chain_id]:
                    packed_mm_feats = mm_db_dict[complex_id][chain_id][res]
                    if config.verbosity >= 6:
                        config.logger.info(f'Setting mm_feat_vec {pos} {complex_id} {chain_id} {res}')
                    position_value_dict[pos][3] = packed_mm_feats

        ts = add_to_times(sub_times, ts)

        del mm_db_dict
        del successful_mm_annotations

        for pos in position_value_dict:
            pos_db_id, packed_recs, packed_mappings, packed_mm_feats = position_value_dict[pos]
            position_values.append((pos_db_id, packed_recs, pack((packed_mappings, packed_mm_feats))))

        ts = add_to_times(sub_times, ts)

        del position_value_dict

        ts = add_to_times(sub_times, ts)

        if packnum == (len(package) - 1):
            if com_queue is not None:
                com_queue.put(0)

        if config.compute_ppi:
            annotation_dict = {}
            for structure_id, chain in annotation_list:
                if structure_id not in annotation_dict:
                    annotation_dict[structure_id] = []
                annotation_dict[structure_id].append(chain)
            aggregated_interfaces: list[interface_package.Aggregated_interface] = interface_package.calculate_aggregated_interface_map(
                config,
                prot_db_id,
                annotation_dict,
                interface_map,
                structure_quality_measures,
                complex_to_db_ids,
                backmap_store,
                db_lock=db_lock)
            
            ts = add_to_times(sub_times, ts)

            if config.verbosity >= 5:
                config.logger.info(f'Len of aggregated_interfaces for {protein_id}: {len(aggregated_interfaces)}')

            if len(aggregated_interfaces) > 0:
                ins_times = remote_insert_interfaces(prot_db_id, aggregated_interfaces, structs_from_db, complexes, config, db_lock=db_lock)
            else:
                ins_times = [0.0]*10

            sub_times.extend(ins_times)

            del annotation_dict
            del aggregated_interfaces

            ts = add_to_times(sub_times, ts)

        del annotation_list
        del structure_quality_measures

        total_sub_times = aggregate_times(total_sub_times, sub_times)

    remote_insert_classifications(position_values, config)
    del position_values

    times.append(total_sub_times)

    del structs_from_db
    del interface_map 
    del complex_to_db_ids

    if config.verbosity >= 4:
        config.logger.info(f'End of para_classify {len(package)=} {para=}')

    return times


def get_res_info(structures, pdb_id, chain, res_nr):

    residue_obj: residue_package.Residue | list = structures[pdb_id][chain].residues.get_item(res_nr)
    if residue_obj is None:
        return None

    if isinstance(residue_obj, list):
        res_id, res_db_id, packed_res_data = residue_obj
        return packed_res_data
    if isinstance(residue_obj, tuple):
        res_id, res_db_id, packed_res_data = residue_obj
        return packed_res_data
    res_info = residue_obj.get_res_info()
    return res_info


def pack_packages(
        send_packages: int,
        classification_results: list,
        proteins,
        config,
        size_sorted,
        para,
        config_store,
        n_procs,
        prot_db_id_map: dict[str, tuple[int, list[int]]],
        cycle_locks,
        stored_structures,
        structure_store,
        backmap_store,
        com_queue
        ):
    
    if config.verbosity >= 5:
        config.logger.info(f'Call of pack_packages with: {len(size_sorted)=}, {send_packages=}, {n_procs=}')
        config.logger.info(f'Current CPU load: {psutil.cpu_percent()}')
    if config.verbosity >= 4:
        if len(size_sorted) == 0:
            config.logger.info('#############\ncalled pack packages with empty task\n##############')

    n_structs = 0
    structure_db_dict: tuple[dict[int, tuple[str, str]], dict[int, tuple[str, str]]] = ({}, {})
    complex_store: dict[str, tuple[dict[str, str], float, int]] = {}

    total_times = []
    package = []

    snap_size_sorted = [x for x in size_sorted]
    for u_ac, size in snap_size_sorted:
        ti = time.time()
        times = []

        protein_obj = proteins[u_ac]
        
        annotation_list = protein_obj.get_annotation_list()
        n_structs += len(annotation_list)

        ti = add_to_times(times, ti)

        for (structure_id, chain) in annotation_list:
            if structure_id not in complex_store:
                try:
                    complex_store[structure_id] = (proteins.complexes[structure_id].chains, proteins.complexes[structure_id].resolution, proteins.complexes[structure_id].database_id)
                except:
                    config.errorlog.add_warning(f'{structure_id} not in complexes. Proteins is None: {proteins is None}')
                    continue

            for c_chain in proteins.complexes[structure_id].chains:
                if proteins.complexes[structure_id].chains[c_chain] != 'Protein':
                    continue
                if c_chain not in proteins.structures[structure_id]:
                    continue
                
                structure_database_id = proteins.structures[structure_id][c_chain].database_id
                if structure_database_id in stored_structures:
                    continue
                #    structure_db_dict[0][structure_database_id] = (structure_id, c_chain)
                if c_chain == chain:
                    
                    if structure_database_id not in structure_db_dict[0]:
                        structure_db_dict[0][structure_database_id] = (structure_id, c_chain)
                        if structure_database_id in structure_db_dict[1]:
                            del structure_db_dict[1][structure_database_id]
                else:
                    if structure_database_id not in structure_db_dict[0] and structure_database_id not in structure_db_dict[1]:
                        structure_db_dict[1][structure_database_id] = (structure_id, c_chain)

        ti = add_to_times(times, ti)

        if n_structs >= 1000 and para:

            if config.verbosity >= 3:
                config.logger.info(f'Sending package: {len(structure_db_dict[0])=} {len(structure_db_dict[1])=}')

            package.append((u_ac, protein_obj))
            
            classification_results.append(para_classify_remote_wrapper.remote(config_store, pack(complex_store), pack(package), pack(structure_db_dict), prot_db_id_map, structure_store, backmap_store, com_queue, next(cycle_locks)))
            
            current_memory_load = psutil.virtual_memory()[2]

            send_packages += 1
            if config.verbosity >= 5:
                config.logger.info(f'Start remote classifcation {u_ac} Amount of different proteins in the package: {len(package)}')
                
                config.logger.info(f'Pending processes: {send_packages}')
                config.logger.info(f'RAM memory % used: {current_memory_load}')
                
            n_structs = 0
            package = []
            structure_db_dict = ({}, {})
            complex_store = {}

            if send_packages >= n_procs or current_memory_load >= 80.0:
                ti = add_to_times(times, ti)
        
                del size_sorted[0]

                ti = add_to_times(times, ti)
                total_times = aggregate_times(total_times, times)
                if config.verbosity >= 6:
                    print_times(total_times, label='pack_packages')
                return send_packages, classification_results, size_sorted
        else:
            package.append((u_ac, protein_obj))

        ti = add_to_times(times, ti)
        
        del size_sorted[0]

        ti = add_to_times(times, ti)
        total_times = aggregate_times(total_times, times)
        
    if config.verbosity >= 6:
        print_times(total_times, label='pack_packages')

    if len(package) > 0 and para:
        classification_results.append(para_classify_remote_wrapper.remote(config_store, pack(complex_store), pack(package), pack(structure_db_dict), prot_db_id_map, structure_store, backmap_store, com_queue, next(cycle_locks)))
        
        if config.verbosity >= 5:
            config.logger.info(f'Start final remote classifcation {u_ac} Amount of different proteins in the package: {len(package)=} {len(size_sorted)=}')
            
        send_packages += 1
        package = []
        structure_db_dict = ({},{})
        complex_store = {}

    if para:
        return send_packages, classification_results, size_sorted
    else:
        return complex_store, package, structure_db_dict



def group_proteins(config, size_sorted, proteins: protein_package.Proteins, max_number_of_groups = 10):
    groups = []
    group_order = []
    group_dict = {}
    count_dict = {}
    complex_count = {}
    for protein_id, size in size_sorted:
        struct_tuple_list = proteins[protein_id].get_annotation_list()
        for structure_id, chain in struct_tuple_list:
            if structure_id not in complex_count:
                complex_count[structure_id] = 1
            else:
                complex_count[structure_id] += 1
            if structure_id not in count_dict:
                count_dict[structure_id] = {}
            if chain not in count_dict[structure_id]:
                count_dict[structure_id][chain] = 1
            else:
                count_dict[structure_id][chain] += 1

        structures = set(struct_tuple_list)

        if len(groups) == 0:
            groups.append([[(protein_id, size)], structures])
            group_order.append([0, len(structures)])
            group_dict[protein_id] = 0
        else:
            broken = False
            for order, (group_number, _) in enumerate(group_order):
                for structure_tuple in structures:
                    if structure_tuple in groups[group_number][1]:
                        groups[group_number][0].append((protein_id, size))
                        groups[group_number][1] = groups[group_number][1].union(structures)
                        group_order[order][1] = len(groups[group_number][1])
                        group_dict[protein_id] = group_number
                        broken = True
                        break
                if broken:
                    break
            if not broken:
                if len(groups) == max_number_of_groups:
                    groups[group_number][0].append((protein_id, size))
                    groups[group_number][1] = groups[group_number][1].union(structures)
                    group_order[order][1] = len(groups[group_number][1])
                    group_dict[protein_id] = group_number
                else:
                    groups.append([[(protein_id, size)], structures])
                    group_order.append([(len(groups) - 1), len(structures)])
                    group_dict[protein_id] = len(groups) - 1
            group_order = sorted(group_order, reverse=True, key=lambda x: x[1])

    """
    count_distribution = {}
    for structure_id in count_dict:
        for chain in count_dict[structure_id]:
            count = count_dict[structure_id][chain]
            if count not in count_distribution:
                count_distribution[count] = 1
            else:
                count_distribution[count] += 1

    config.logger.info(count_distribution)
    """
    
    if config.verbosity >= 4:
        config.logger.info(f'Group order: {group_order}')

    return groups, group_dict, count_dict, complex_count

def group_clds(groups, proteins, config):
    clds = []
    for protein_ids, _ in groups:
        cld = generate_classification_dump(protein_ids, config, proteins = proteins)
        clds.append(cld)
    return clds

def group_structure_stores(groups, proteins):
    structure_stores = []
    #print_count = 0
    for _, structure_tuples in groups:
        structures = {}
        complex_store = {}
        for structure_id, chain in structure_tuples:
            if structure_id not in structures:
                structures[structure_id] = {}
            structures[structure_id][chain] = proteins.structures[structure_id][chain]
            if structure_id not in complex_store:
                complex_store[structure_id] = (proteins.complexes[structure_id].chains, proteins.complexes[structure_id].resolution, proteins.complexes[structure_id].database_id)

        structure_stores.append(ray.put(pack((structures, complex_store))))
    return structure_stores

def generate_classification_dump(size_sorted, config, proteins = None, protein_map = None, complexes = None, unpacked = False):
    t0 = time.time()
    complex_store = {}

    if proteins is not None:
        complexes = proteins.complexes


    total_number_of_backmaps = 0

    for siso in size_sorted:
        try:
            protein_id, _ = siso
        except:
            protein_id = siso

        if proteins is None:
            annotation_list = protein_map[protein_id].get_annotation_list()
        else:
            annotation_list = proteins.get_protein_annotation_list(protein_id)

        for (structure_id, chain) in annotation_list:

            if structure_id not in complex_store:
                if proteins is not None:
                    try:
                        complex_store[structure_id] = (complexes[structure_id].chains, complexes[structure_id].resolution, complexes[structure_id].database_id)
                    except:
                        config.errorlog.add_warning(f'{structure_id} not in complexes. Proteins is None: {proteins is None}')
                        continue
                else:
                    complex_store[structure_id] = (complexes[structure_id][0], complexes[structure_id][1], complexes[structure_id][2])


    t1 = time.time()
    if config.verbosity >= 5:
        size_in_gb_of_complex_store = sys.getsizeof(complex_store) / 1024 / 1024
        config.logger.info(f'Time for cld part 1: {t1-t0}')
        config.logger.info(f'Number of backmaps: {total_number_of_backmaps}')
        config.logger.info(f'Data sizes: {size_in_gb_of_complex_store}')

    if unpacked:
        return False, (config, complex_store)
    
    pcs = pack(complex_store)
    #packed_cld = pack((complex_store, structure_quality_measures, nested_structures))

    t2 = time.time()
    if config.verbosity >= 5:
        config.logger.info(f'Time for cld part 2: {t2-t1}')    
        size_in_gb_of_pcs = sys.getsizeof(pcs) / 1024 / 1024
        config.logger.info(f'packed data size: {size_in_gb_of_pcs}')

    store_reference =  ray.put((config, pcs))

    t3 = time.time()
    if config.verbosity >= 5:
        config.logger.info(f'Time for cld part 3: {t3-t2}')

    return True, store_reference

def process_classification_outputs(config: Config, outs, proteins, recommended_complexes):
    for protein_id, pos_outs in outs:
        if config.verbosity >= 5:
            config.logger.info(f'Processing aggregation results: {protein_id} {len(pos_outs)}')
        
        if protein_id not in recommended_complexes:
            recommended_complexes[protein_id] = {}
        for pos, packed_mappings_obj, recommendation_order, packed_res_recommends in pos_outs:
            if config.verbosity >= 7:
                config.logger.info(f'In process_classification_outputs: {protein_id} {pos}')
            proteins.protein_map[protein_id].positions[pos].mappings = packed_mappings_obj
            proteins.protein_map[protein_id].positions[pos].packed_res_recommends = packed_res_recommends
            recommended_complexes[protein_id][pos] = recommendation_order

@ray.remote(max_calls = 1)
def para_mm_lookup(package, store):
    config, n_sub_procs = store
    complex_feature_dict_dict = {}
    for complex_id in package:
        feature_dict, _ = mm_lookup(complex_id, config, n_sub_procs = n_sub_procs)
        complex_feature_dict_dict[complex_id] = {}
        for chain in feature_dict:
            if chain not in complex_feature_dict_dict[complex_id]:
                complex_feature_dict_dict[complex_id][chain] = {}
            for res_id in feature_dict[chain]:
                microminer_features = mappings_package.Microminer_features()
                microminer_features.set_values(feature_dict[chain][res_id])
                if res_id not in complex_feature_dict_dict[complex_id][chain]:
                    complex_feature_dict_dict[complex_id][chain][res_id] = pack(microminer_features)

    return complex_feature_dict_dict

def single_mm_lookup(package: set[str], config: Config)-> dict[str, dict[str, dict[str | int, bytes]]]:
    complex_feature_dict_dict: dict[str, dict[str, dict[str | int, bytes]]] = {}
    feature_dict: dict[str, dict[str | int, list[int | float | None]]]
    for complex_id in package:
        feature_dict, _ = mm_lookup(complex_id, config, n_sub_procs = 1)
        complex_feature_dict_dict[complex_id] = {}
        for chain in feature_dict:
            if chain not in complex_feature_dict_dict[complex_id]:
                complex_feature_dict_dict[complex_id][chain] = {}
            for res_id in feature_dict[chain]:
                microminer_features = mappings_package.Microminer_features()
                microminer_features.set_values(feature_dict[chain][res_id])
                if res_id not in complex_feature_dict_dict[complex_id][chain]:
                    complex_feature_dict_dict[complex_id][chain][res_id] = pack(microminer_features)
    return complex_feature_dict_dict

def classification(proteins, config, background_insert_residues_process, custom_structure_annotations=None):

    if config.verbosity >= 2:
        t0 = time.time()

    classification_results = []

    size_sorted = []

    if custom_structure_annotations is None:
        protein_ids = proteins.get_protein_ids()
    else:
        protein_ids = custom_structure_annotations.keys()

    total_mappings = 0

    prot_db_id_map: dict[int, list[int]] = {}

    for protein_id in protein_ids:
        prot_db_id_map[proteins[protein_id].database_id] = [None]
        for pos in proteins[protein_id].get_position_ids():
            prot_db_id_map[proteins[protein_id].database_id].append(proteins[protein_id].positions[pos].database_id)

        if proteins.protein_map[protein_id].stored:
            continue
        if len(proteins[protein_id].structure_annotations) == 0:
            continue

        size = proteins[protein_id].number_of_mappings
        size_sorted.append((protein_id, size))

        total_mappings += size

    if config.verbosity >= 3:
        config.logger.info(f'Total mappings: {total_mappings:_}')

    if config.verbosity >= 5:
        config.logger.info(f'{prot_db_id_map=}')
        if len(proteins.structures.keys()) <= 1000:
            config.logger.info(f'{proteins.structures.keys()=}')
        else:
            config.logger.info(f'Size of protein.structures: {len(proteins.structures)}')

    size_sorted = sorted(size_sorted, reverse=True, key=lambda x: x[1])

    para = (total_mappings > 20000) and (len(size_sorted) > 1)
    if config.proc_n == 1:
        para = False

    send_packages = 0

    if config.verbosity >= 3:
        config.logger.info(f'Begin of aggregation process, {para=}')

    
    if config.verbosity >= 2:
        t1 = time.time()
        config.logger.info(f'Aggregation part 1: {t1-t0} {len(size_sorted)=}')

    if para:
        config_store = ray.put(config)
        prot_db_id_map_store = ray.put(prot_db_id_map)

        number_of_locks = min([30 ,max([1, config.max_connections//2])])

        db_locks = [f'db_lock_{i}.lock' for i in range(number_of_locks)]
        cycle_locks = cycle(db_locks)

        if config.verbosity >= 2:
            t21 = time.time()
            config.logger.info(f'Aggregation part 2.1: {t21-t1}')

        groups, group_dict, count_dict, complex_count = group_proteins(config, size_sorted, proteins, max_number_of_groups = config.proc_n)

        if config.verbosity >= 2:
            t22 = time.time()
            config.logger.info(f'Aggregation part 2.2: {t22-t21} {len(groups)=}')

        count_thresh = max([min([config.proc_n, len(size_sorted)])//150, 1])
        complex_count_thresh = max([min([config.proc_n, len(size_sorted)])//150, 1])

        structure_store_db_dict: dict[int, tuple[str, str]] = {}
        stored_structures: set[int] = set()
        for structure_id in count_dict:
            for chain in count_dict[structure_id]:
                count = count_dict[structure_id][chain]
                if count > count_thresh:
                    struct_db_id = proteins.structures[structure_id][chain].database_id
                    structure_store_db_dict[struct_db_id] = (structure_id, chain)
                    stored_structures.add(struct_db_id)

        backmap_store_db_ids: dict[int, tuple[str, str]] = {}
        for structure_id in complex_count:
            count = complex_count[structure_id]
            if count > complex_count_thresh:
                for chain in proteins.complexes[structure_id].chains:
                    if proteins.complexes[structure_id].chains[chain] != 'Protein':
                        continue
                    if chain not in proteins.structures[structure_id]:
                        config.errorlog.add_warning(f'{chain=} not in proteins.structures[{structure_id=}], most likely Protein/Peptide mislabeling')
                        continue
                    struct_db_id = proteins.structures[structure_id][chain].database_id

                    backmap_store_db_ids[struct_db_id] = (structure_id, chain)

        if config.verbosity >= 2:
            t23 = time.time()
            config.logger.info(f'Aggregation part 2.3: {t23-t22} {len(stored_structures)=} {len(backmap_store_db_ids)=} {count_thresh=} {complex_count_thresh=}')

        structs_from_db, interface_map = getStoredResidues(None, config, custom_ids = structure_store_db_dict, exclude_interacting_chains = True)

        if config.verbosity >= 2:
            t24 = time.time()
            config.logger.info(f'Aggregation part 2.4: {t24-t23}')

        structure_backmaps: dict[str, dict[str, dict[int, bytes]]] = getBackmaps(
                backmap_store_db_ids,
                config,
                unpacked=False,
                locked = False)

        backmap_store = {}
        for structure_id in structure_backmaps:
            backmap_store[structure_id] = ray.put(structure_backmaps[structure_id])


        if config.verbosity >= 2:
            t25 = time.time()
            config.logger.info(f'Aggregation part 2.5: {t25-t24}')

        packed_structures = {}
        N = 0
        for structure_id in structs_from_db:
            packed_structures[structure_id] = {}
            for chain in structs_from_db[structure_id]:
                res_data = structs_from_db[structure_id][chain]
                packed_structures[structure_id][chain] = ray.put(pack(res_data))
                N += 1

        packed_interfaces = {}
        for structure_id in interface_map:
            packed_interfaces[structure_id]= ray.put(pack(interface_map[structure_id]))


        structure_store = (packed_structures, packed_interfaces)

        if config.verbosity >= 2:
            t26 = time.time()
            config.logger.info(f'Aggregation part 2.6: {t26-t25} Packed structures: {N=}')

        com_queue = Queue()

        start_with = 0
        for group_number, (size_sorted, _) in enumerate(groups[start_with:]):
            
            (send_packages, classification_results, size_sorted) = pack_packages(send_packages, classification_results,
                                                                                 proteins, config, size_sorted,
                                                                                 para, config_store, config.proc_n,
                                                                                 prot_db_id_map_store, cycle_locks,
                                                                                 stored_structures, structure_store,
                                                                                 backmap_store, com_queue)

            if send_packages >= config.proc_n:
                break
            start_with += 1

        if config.verbosity >= 2:
            t27 = time.time()
            config.logger.info(f'Aggregation part 2.7: {t27-t26} {send_packages=} {group_number=} {len(size_sorted)=}')
        

    else:

        (complex_store, package, structure_db_dict) = pack_packages(send_packages, classification_results, proteins,
                                            config, size_sorted, False, None, config.proc_n, prot_db_id_map, None, {}, None, None, None)
        
        agg_times = para_classify(config, complex_store, package, structure_db_dict, prot_db_id_map)
        if config.verbosity >= 3:
            print_times(agg_times, label='Serialized aggregation', logger=config.logger)


    if config.verbosity >= 2:
        t2 = time.time()
        config.logger.info(f'Time for aggregation part 2: {t2 - t1}, {para=}')

    if para:
        total_times = []

        total_agg_times = []

        loop_count = 0

        while True:
            loop_count += 1
            if loop_count >= 5000:
                if config.verbosity >= 3:
                    config.logger.info(f'5k waits {len(classification_results)=}')
                    loop_count = 0
            ready, not_ready = ray.wait(classification_results, timeout = 0.1)

            if len(ready) > 0:

                #send_packages -= len(ready)

                remote_returns = ray.get(ready)

                for returned_package in remote_returns:
                    agg_times = returned_package
                    total_agg_times = aggregate_times(total_agg_times, agg_times)
                    if config.verbosity >= 7:
                        print_times(agg_times, label='Para aggregation', logger=config.logger)

                classification_results = not_ready
            else:
                classification_results = not_ready


            if not com_queue.empty():
                ti = time.time()
                times = []
                
                com_queue.get()

                send_packages -= 1
                current_memory_load = psutil.virtual_memory()[2]
                if (len(size_sorted) > 0 or start_with < len(groups)) and current_memory_load <= 75.0:
                    sw = start_with
                    for group_number, (size_sorted, _) in enumerate(groups[start_with:]):
                        if len(size_sorted) > 0:
                            group_number += sw
                            (send_packages, classification_results, size_sorted) = pack_packages(send_packages, classification_results, proteins,
                                                                                                config, size_sorted, para, config_store,
                                                                                                config.proc_n, prot_db_id_map_store, cycle_locks,
                                                                                                stored_structures, structure_store, backmap_store, com_queue)
                            if config.verbosity >= 6:
                                config.logger.info(f'Packaged new aggregation process {group_number=} {send_packages=} {len(size_sorted)=}')
                            if send_packages >= config.proc_n:
                                break
                        start_with += 1
                ti = add_to_times(times, ti)
                total_times = aggregate_times(total_times, times)

            if len(classification_results) == 0:
                break

        if config.verbosity >= 3:
            print_times(total_times, label = 'Processing and packaging', logger=config.logger)
            print_times(total_agg_times, label='Para aggregation', logger=config.logger)


    else:
        gc.collect()


    if config.verbosity >= 2:
        t2 = time.time()
        config.logger.info(f'Time for classification part 2: {t2 - t1}')

    if background_insert_residues_process is not None:
        background_insert_residues_process.join()
        background_insert_residues_process.close()
        background_insert_residues_process = None

    if config.verbosity >= 2:
        t3 = time.time()
        config.logger.info(f'Time for classification part 3: {t3 - t2}')

def calc_mm_features(
        config,
        recommended_complexes: list[list[tuple[str, str, int | str]]]
        )-> tuple[dict[str, dict[str, dict[str | int, bytes]]], dict[int, tuple[str, str, str | int]]]:
    mm_db_dict: dict[str, dict[str, dict[str | int, bytes]]] = {}
    successful_mm_annotations: dict[int, tuple[str, str, str | int]] = {}
    if config.microminer_exe != '':
        mm_inputs: set[str] = set()
        processed_complexes = set()
        rec_complex_pos_dict = {}
        
        for pos, recommendation_order in enumerate(recommended_complexes):
            if recommendation_order is None:
                continue
            if len(recommendation_order) > 0:
                recommended_structure_triple = recommendation_order[0]
                complex_id = recommended_structure_triple[0]
                mm_inputs.add(complex_id)
                processed_complexes.add(complex_id)
                if complex_id not in rec_complex_pos_dict:
                    rec_complex_pos_dict[complex_id] = []
                rec_complex_pos_dict[complex_id].append((pos, 0))

        complex_feature_dict_dict: dict[str, dict[str, dict[str | int, bytes]]] = single_mm_lookup(mm_inputs, config)

        while len(mm_inputs) > 0:
            new_mm_inputs = set()
            n_empty = 0
            for complex_id in complex_feature_dict_dict:
                feature_dict = complex_feature_dict_dict[complex_id]

                if len(feature_dict) == 0:
                    n_empty += 1
                    for (pos, current_try) in rec_complex_pos_dict[complex_id]:
                        stay_in_loop = True
                        next_try = current_try + 1
                        while len(recommended_complexes[pos]) > next_try and stay_in_loop:
                            recommended_structure_triple = recommended_complexes[pos][next_try]
                            complex_id, chain_id, res = recommended_structure_triple

                            if complex_id in mm_db_dict:
                                if (chain_id, res) in mm_db_dict[complex_id]:
                                    successful_mm_annotations[pos] = complex_id, chain_id, res
                                    stay_in_loop = False
                                else:
                                    next_try += 1
                            elif complex_id not in processed_complexes:
                                new_mm_inputs.add(complex_id)
                                processed_complexes.add(complex_id)
                                if complex_id not in rec_complex_pos_dict:
                                    rec_complex_pos_dict[complex_id] = []
                                rec_complex_pos_dict[complex_id].append((pos, next_try))
                                stay_in_loop = False
                            else:
                                next_try += 1

                else:
                    mm_db_dict[complex_id] = feature_dict
                    for (pos, current_try) in rec_complex_pos_dict[complex_id]:
                        recommended_structure_triple = recommended_complexes[pos][current_try]
                        complex_id, chain_id, res = recommended_structure_triple
                        if chain_id in feature_dict:
                            if res in feature_dict[chain_id]:
                                successful_mm_annotations[pos] = complex_id, chain_id, res
                                success = True
                            else:
                                success = False
                        else:
                            success = False
                        if not success:
                            next_try = current_try + 1
                            stay_in_loop = True
                            while len(recommended_complexes[pos]) > next_try and stay_in_loop:
                                recommended_structure_triple = recommended_complexes[pos][next_try]
                                complex_id, chain_id, res = recommended_structure_triple
                                if complex_id in mm_db_dict:
                                    successful_mm_annotations[pos] = complex_id, chain_id, res
                                    stay_in_loop = False
                                elif complex_id not in processed_complexes:
                                    new_mm_inputs.add(complex_id)
                                    processed_complexes.add(complex_id)
                                    if complex_id not in rec_complex_pos_dict:
                                        rec_complex_pos_dict[complex_id] = []
                                    rec_complex_pos_dict[complex_id].append((pos, next_try))
                                    stay_in_loop = False
                                else:
                                    next_try += 1

            if config.verbosity >= 4:
                config.logger.info(f'MicroMiner-lookup succesful for {len(mm_inputs) - n_empty} Complexes.')
            
            mm_inputs = new_mm_inputs
            complex_feature_dict_dict = single_mm_lookup(mm_inputs, config)

    return mm_db_dict, successful_mm_annotations