import ray
import time
import psutil
import gc
import math
import sys
import traceback

from ray.util.queue import Queue

from structman.base_utils.ray_utils import ray_hack
from structman.scripts.createMMDB import mm_lookup
from structman.scripts.createMMDB import feature_names as microminer_feature_names
from structman.base_utils.base_utils import calculate_chunksizes, pack, unpack, is_alphafold_model
from structman.lib.sdsc import mappings as mappings_package
from structman.lib.sdsc import interface as interface_package
from structman.lib.rin import Interaction_profile, Centrality_scores, aggregate_IAmaps
from structman.lib.sdsc.sdsc_utils import get_shortest_distances, triple_locate, rin_classify
from structman.lib.database.insertion_lib import insert_interfaces

def qualityScore(resolution, coverage, seq_id, resolution_wf=0.25, coverage_wf=0.5, seq_id_wf=1.0):
    seq_id = float(seq_id)
    resolution = float(resolution)
    coverage = float(coverage)
    if seq_id > 1.0:
        seq_id = seq_id / 100.0
    if coverage > 1.0:
        coverage = coverage / 100.0

    # project the criteria to [0,1] via logistic regression
    resolution_value = (1 + math.exp((1.5 * resolution) - 4))**(-1)

    seq_value = (1 + math.exp(10 * (0.4 - seq_id)))**(-1)

    ws = sum((resolution_wf, coverage_wf, seq_id_wf))
    quality = sum((resolution_value * resolution_wf, coverage * coverage_wf, seq_value * seq_id_wf)) / ws
    return quality

@ray.remote(max_calls = 1)
def para_classify_remote_wrapper(classification_dump, package, package_size):
    ray_hack()

    return (para_classify(classification_dump, unpack(package), para=True), package_size)


def para_classify(classification_dump, package, para=False):
    t0 = time.time()
    if len(classification_dump) > 2:
        config, complexes, structure_quality_measures, mapped_positions, interface_map_store, backmap_store = classification_dump
    else:
        config, packed_dump = classification_dump
        complexes, structure_quality_measures, mapped_positions, interface_map_store, backmap_store = unpack(packed_dump)

    outs = []
    for protein_id, classification_inp, annotations in package:

        backmaps = {}
        interface_map = {}
        for (structure_id, chain) in annotations:
            backmaps[(structure_id, chain)] = backmap_store[(protein_id, structure_id, chain)]
            interface_map[(structure_id, chain)] = interface_map_store[(structure_id, chain)]

        pos_outs = []

        for pos, mappings, disorder_score, disorder_region in classification_inp:

            mappings_obj = mappings_package.Mappings()
            #print(f"In para_classify: {protein_id} {pos} {len(mappings)}")

            for mapp in mappings:

                if para:
                    id_triple, quality_measures, mapping = mapp

                    (pdb_id, chain, res_nr) = id_triple
                    (seq_id, cov, identical_aa) = quality_measures
                    try:
                        (chains, resolution) = complexes[pdb_id]
                    except:
                        (chains, resolution) = (None, None)
                        
                else:
                    id_triple, quality_measures, structure_infos = mapp

                    (pdb_id, chain, res_nr) = id_triple
                    (seq_id, cov, identical_aa) = quality_measures

                    structure_info, mapping = structure_infos
                    (chains, resolution) = structure_info

                (rsa, mc_rsa, sc_rsa, ssa, profile_or_str, centralities_or_str,
                 phi, psi, intra_ssbond, inter_ssbond, ssbond_length, intra_link, inter_link, link_length, cis_conformation, cis_follower,
                 inter_chain_median_kd, inter_chain_dist_weighted_kd,
                 inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                 intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                 inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                 intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                 b_factor, modres,
                 lig_dists, chain_distances, homomer_distances, res_aa) = mapping

                gsd_return = get_shortest_distances(chains, lig_dists, chain_distances, homomer_distances)

                if gsd_return is None:
                    continue
                homo_dist, lig_dist, metal_dist, ion_dist, chain_dist, rna_dist, dna_dist, min_lig, min_metal, min_ion, iacs = gsd_return

                if isinstance(profile_or_str, str):
                    profile = Interaction_profile(profile_str=profile_or_str)
                else:
                    profile = profile_or_str

                if isinstance(centralities_or_str, str):
                    centralities = Centrality_scores(code_str=centralities_or_str)
                else:
                    centralities = centralities_or_str

                loc, mc_loc, sc_loc = triple_locate(rsa, mc_rsa, sc_rsa, config)

                rin_class, rin_simple_class = rin_classify(profile, sc_loc)

                qual = qualityScore(resolution, cov, seq_id)

                structural_feature_dict = {
                    'b_factor' : b_factor, 'modres' : modres, 'ssa' : ssa, 'phi' : phi, 'psi' : psi, 'intra_ssbond' : intra_ssbond, 'inter_ssbond' : inter_ssbond, 'ssbond_length' : ssbond_length,
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
         	
            mappings_obj.weight_all(config, disorder_score, disorder_region)
            
            mapping_results = mappings_obj.get_raw_result()

            pos_outs.append((pos, mapping_results))

            mappings_obj.deconstruct()
            del mappings_obj

        if config.compute_ppi:
            aggregated_interfaces = interface_package.calculate_aggregated_interface_map(config, protein_id, interface_map, backmaps, structure_quality_measures, mapped_positions)
        else:
            aggregated_interfaces = []

        outs.append((protein_id, aggregated_interfaces, pos_outs))

    if para:
        outs = pack(outs)

    t1 = time.time()

    return outs, t1 - t0


def get_res_info_from_store(structure, res_nr):
    res_aa = structure.residues[res_nr].aa

    centralities_or_str = structure.residues[res_nr].get_centralities(get_whats_there=True)
    modres = structure.residues[res_nr].modres
    b_factor = structure.residues[res_nr].b_factor
    rsa, mc_rsa, sc_rsa = structure.get_residue_rsa_triple(res_nr)
    ssa = structure.residues[res_nr].SSA
    profile_or_str = structure.residues[res_nr].get_interaction_profile(get_whats_there=True)

    phi = structure.residues[res_nr].phi
    psi = structure.residues[res_nr].psi

    intra_ssbond, inter_ssbond, ssbond_length, intra_link, inter_link, link_length, cis_conformation, cis_follower = structure.get_residue_link_information(res_nr)

    (inter_chain_median_kd, inter_chain_dist_weighted_kd,
     inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
     intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa) = structure.get_residue_milieu(res_nr)

    (inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
     intra_chain_interactions_median, intra_chain_interactions_dist_weighted) = structure.residues[res_nr].get_interface_milieu()

    lig_dists = structure.get_residue_sld(res_nr)
    chain_distances = structure.get_residue_scd(res_nr)
    homomer_distances = structure.get_residue_homomer_dists(res_nr)

    res_info = (rsa, mc_rsa, sc_rsa, ssa, profile_or_str, centralities_or_str,
                phi, psi, intra_ssbond, inter_ssbond, ssbond_length, intra_link, inter_link, link_length, cis_conformation, cis_follower,
                inter_chain_median_kd, inter_chain_dist_weighted_kd,
                inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                b_factor, modres, lig_dists, chain_distances, homomer_distances, res_aa)
    return res_info


def get_res_info(structures, pdb_id, chain, res_nr):

    residue_obj = structures[(pdb_id, chain)].residues[res_nr]

    centralities_or_str = residue_obj.get_centrality_str()

    rsa, mc_rsa, sc_rsa = residue_obj.get_rsa(splitted = True)
    profile_or_str = residue_obj.get_interaction_profile_str()

    intra_ssbond, inter_ssbond, ssbond_length, intra_link, inter_link, link_length, cis_conformation, cis_follower = residue_obj.get_residue_link_information()

    (inter_chain_median_kd, inter_chain_dist_weighted_kd,
     inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
     intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa) = residue_obj.get_milieu()

    (inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
     intra_chain_interactions_median, intra_chain_interactions_dist_weighted) = residue_obj.get_interface_milieu()

    lig_dists = residue_obj.get_ligand_distances()
    chain_distances = residue_obj.get_chain_distances()
    homomer_distances = residue_obj.get_homomer_dists()

    res_info = (rsa, mc_rsa, sc_rsa, residue_obj.SSA, profile_or_str, centralities_or_str,
                residue_obj.phi, residue_obj.psi, intra_ssbond, inter_ssbond, ssbond_length, intra_link, inter_link, link_length, cis_conformation, cis_follower,
                inter_chain_median_kd, inter_chain_dist_weighted_kd,
                inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                residue_obj.b_factor, residue_obj.modres, lig_dists, chain_distances, homomer_distances, residue_obj.aa)
    return res_info

"""
def add_structure_to_store(structure_store, proteins, pdb_id, chain):
    chains = proteins.get_complex_chains(pdb_id)
    resolution = proteins.get_resolution(pdb_id)
    if resolution is None:
        if config.verbosity >= 4:
            print('Skipped classification of', u_ac, 'due to', pdb_id, 'got no resolution')
        return
    res_map = {}
    for res_nr in proteins.structures[(pdb_id, chain)].residues:
        res_map[res_nr] = get_res_info(proteins.structures, pdb_id, chain, res_nr)
      
    structure_store[(pdb_id, chain)] = ray.put((chains, resolution, res_map))
    del chains
    del resolution
    del res_map
"""

def pack_packages(package_size, package, send_packages, classification_results, protein_map, structures, config, size_sorted, processed_positions, para, classification_dump, max_package_size, n_procs, proteins = None):
    if config.verbosity >= 5:
        print('Call of pack_packages with:', package_size, len(package), send_packages)
        print('Current CPU load:', psutil.cpu_percent())
    if config.verbosity >= 5:
        if len(size_sorted) == 0:
            print('##############\n#############\ncalled pack packages with empty task\n##############\n#############')

    t_pack_0 = 0
    t_pack_01 = 0
    t_pack_02 = 0
    t_pack_03 = 0
    t_pack_1 = 0
    t_pack_11 = 0
    t_pack_12 = 0
    t_pack_13 = 0
    t_pack_14 = 0
    t_pack_15 = 0
    t_pack_16 = 0
    t_pack_21 = 0
    t_pack_22 = 0

    snap_size_sorted = [x for x in size_sorted]
    for u_ac, size in snap_size_sorted:
        classification_inp = []
        protein_obj = protein_map[u_ac]

        positions = protein_obj.get_position_ids()
        annotation_list = protein_obj.get_annotation_list()

        if not u_ac in processed_positions:
            processed_positions[u_ac] = set()

        for pos in positions:
            if pos in processed_positions[u_ac]:
                continue
            if protein_obj.is_position_stored(pos):
                continue
            t_pack_0 += time.time()

            mappings = []
            model_in_mappings = False
            aacbase = protein_obj.get_aac_base(pos)
            disorder_score = protein_obj.get_disorder_score(pos)
            disorder_region = protein_obj.get_disorder_region(pos)

            t_pack_01 += time.time()
            
            error_count = 0
            for (pdb_id, chain) in annotation_list:

                t_pack_11 += time.time()

                try:
                    sub_info = protein_obj.structure_annotations[(pdb_id, chain)].get_sub_info(pos)

                    t_pack_12 += time.time()
                except:
                    if config.verbosity >= 4:
                        print('Skipped classification of', u_ac, 'due to pos', pos, 'was not in sub_infos')

                    t_pack_12 += time.time()
                    continue

                if sub_info is None: #In this case the position is not mapped to any residue in the structure (gap in the alignment and getSubInfo called with ignore_gaps) 
                    continue

                res_nr = sub_info[0]

                if res_nr is None: #In this case the position is not mapped to any residue in the structure (gap in the alignment)
                    continue

                t_pack_13 += time.time()

                seq_id = protein_obj.structure_annotations[(pdb_id, chain)].get_sequence_id()

                t_pack_14 += time.time()

                if seq_id is None:
                    if config.verbosity >= 4:
                        print('Skipped classification of', u_ac, 'due to', pdb_id, 'got no sequence identity')
                    continue

                t_pack_15 += time.time()

                cov = protein_obj.structure_annotations[(pdb_id, chain)].get_coverage()

                t_pack_16 += time.time()

                try:
                    res_aa = structures[(pdb_id, chain)].get_residue_aa(res_nr)
                except:
                    error_count += 1
                    if error_count >= 2:
                        continue
                    [e, f, g] = sys.exc_info()
                    g = traceback.format_exc()
                    if config.verbosity >= 5:
                        print(f'Skipped annotation of {u_ac} {pos} to {pdb_id} {chain} {res_nr}:\n{e}\n{f}\n{g}')
                    continue

                if res_aa is None:
                    if config.verbosity >= 5:
                        print(f'Skipped annotation of {u_ac} {pos} to {pdb_id} {chain} {res_nr}: res_aa not defined')
                    continue

                identical_aa = res_aa == aacbase[0]

                t_pack_21 += time.time()

                if is_alphafold_model(pdb_id):
                    if len(mappings) > 0: #Mapped residues from models are only going into results aggregation, if it is the only mapped structure
                        continue
                    model_in_mappings = True
                elif model_in_mappings: #This happens only, if a model was the first mapped structure
                    mappings = []
                    model_in_mappings = False

                if not para:
                    try:
                        resolution = proteins.get_resolution(pdb_id)
                        if resolution is None:
                            if config.verbosity >= 4:
                                print('Skipped classification of', u_ac, 'due to', pdb_id, 'got no resolution')
                            continue
                        package_size += 1
                        chains = proteins.get_complex_chains(pdb_id)
                        mappings.append(((pdb_id, chain, res_nr), (seq_id, cov, identical_aa), ((chains, resolution), get_res_info(proteins.structures, pdb_id, chain, res_nr))))
                    except:
                        package_size += 1
                        mappings.append(((pdb_id, chain, res_nr), (seq_id, cov, identical_aa), get_res_info(structures, pdb_id, chain, res_nr)))

                else:
                    package_size += 1
                    mappings.append(((pdb_id, chain, res_nr), (seq_id, cov, identical_aa), get_res_info(structures, pdb_id, chain, res_nr)))

                t_pack_22 += time.time()


            if error_count >= 5 and config.verbosity >= 4:
                print(f'High number of error count: {error_count} in aggregation of {u_ac} {pos} in pack_packages')

            t_pack_02 += time.time()

            classification_inp.append((pos, mappings, disorder_score, disorder_region))

            t_pack_03 += time.time()

            processed_positions[u_ac].add(pos)

            t_pack_1 += time.time()

            if package_size >= max_package_size and para:

                package.append((u_ac, classification_inp, annotation_list))

                t_send_0 = time.time()
                classification_results.append(para_classify_remote_wrapper.remote(classification_dump, pack(package), package_size))
                t_send_1 = time.time()

                send_packages += 1
                if config.verbosity >= 5:
                    print('Start remote classifcation', u_ac, 'Package size:', package_size,'Amount of different proteins in the package:', len(package))
                    print('Amount of positions in last protein:', len(classification_inp))
                    print('Pending processes:', send_packages)
                    print('RAM memory % used:', psutil.virtual_memory()[2])
                    print('Time for sending the process:', (t_send_1 - t_send_0), 'Time for packaging:', (t_pack_1 - t_pack_0), (t_pack_01 - t_pack_0), (t_pack_02 - t_pack_01), (t_pack_03 - t_pack_02), (t_pack_1 - t_pack_03), (t_pack_22 - t_pack_21))
                    print('More times:', (t_pack_12 - t_pack_11), (t_pack_14 - t_pack_13), (t_pack_16 - t_pack_15), (t_pack_21 - t_pack_16))
                    t_pack_0 = 0
                    t_pack_01 = 0
                    t_pack_02 = 0
                    t_pack_03 = 0
                    t_pack_1 = 0
                    t_pack_11 = 0
                    t_pack_12 = 0
                    t_pack_13 = 0
                    t_pack_14 = 0
                    t_pack_15 = 0
                    t_pack_16 = 0
                    t_pack_21 = 0
                    t_pack_22 = 0
                package_size = 0
                package = []
                classification_inp = []

                if send_packages >= n_procs:
                    return package_size, package, send_packages, classification_results, size_sorted, processed_positions

        package.append((u_ac, classification_inp, annotation_list))
        del size_sorted[0]
        del processed_positions[u_ac]
    if len(package) > 0 and para:
        classification_results.append(para_classify_remote_wrapper.remote(classification_dump, pack(package), package_size))
        send_packages += 1
        package = []
        package_size = 0

    return package_size, package, send_packages, classification_results, size_sorted, processed_positions


@ray.remote(max_calls = 1)
def nested_classification_main_process(queue, config, size_sorted, max_package_size, n_procs, data_package):

    protein_map, structures, complexes, mapped_positions_map = unpack(data_package)

    if config.verbosity >= 5:
        print('Call of nested classification main process:', len(size_sorted), max_package_size, n_procs)

    if config.verbosity >= 6:
        print(f'\nMapped positions map:\n{mapped_positions_map}')

    package_size = 0
    package = []
    send_packages = 0
    processed_positions = {}
    classification_results = []

    if n_procs > 1:
        classification_dump = generate_classification_dump(size_sorted, config, protein_map = protein_map, complexes = complexes, structures = structures, mapped_positions_map = mapped_positions_map)

        package_size, package, send_packages, classification_results, size_sorted, processed_positions = pack_packages(package_size, package, send_packages, classification_results, protein_map, structures, config, size_sorted, processed_positions, True, classification_dump, max_package_size, n_procs)

        while True:
            ready, not_ready = ray.wait(classification_results)

            if len(ready) > 0:

                send_packages -= len(ready)

                queue.put(ready)

                classification_results = not_ready

                t_pack_0 = time.time()
                package_size, package, send_packages, classification_results, size_sorted, processed_positions = pack_packages(package_size, package, send_packages, classification_results, protein_map, structures, config, size_sorted, processed_positions, True, classification_dump, max_package_size, n_procs)
                t_pack_1 = time.time()

                if config.verbosity >= 5:
                    print('Time for packaging new processes:', (t_pack_1 - t_pack_0))
            else:
                classification_results = not_ready

            if len(classification_results) == 0 and len(package) == 0:
                queue.put('End')
                del classification_dump
                return
    else:
        classification_dump = generate_classification_dump(size_sorted, config, protein_map = protein_map, complexes = complexes, structures = structures, mapped_positions_map = mapped_positions_map, unpacked=True)

        package_size, package, send_packages, classification_results, size_sorted, processed_positions = pack_packages(package_size, package, send_packages, classification_results, protein_map, structures, config, size_sorted, processed_positions, False, classification_dump, max_package_size, 1)

        if len(package) > 0:
            classification_results = para_classify(classification_dump, package, para=True)
        del classification_dump
        return classification_results


def generate_classification_dump(size_sorted, config, proteins = None, protein_map = None, complexes = None, structures = None, mapped_positions_map = None, unpacked = False):
    complex_store = {}
    structure_quality_measures = {}
    interface_map_store = {}
    backmap_store = {}
    mapped_positions = {}


    if not proteins is None:
        complexes = proteins.complexes
        structures = proteins.structures

    for protein_id, size in size_sorted:
        if proteins is None:
            annotation_list = protein_map[protein_id].get_annotation_list()
        else:
            annotation_list = proteins.get_protein_annotation_list(protein_id)

        for (structure_id, chain) in annotation_list:

            if not structure_id in complex_store:
                try:
                    complex_store[structure_id] = (complexes[structure_id].chains, complexes[structure_id].resolution)
                except:
                    config.errorlog.add_warning(f'{structure_id} not in complexes. Proteins is None: {proteins is None}')
                    continue

            if proteins is None:
                structure_quality_measures[(protein_id, structure_id, chain)] = (protein_map[protein_id].structure_annotations[(structure_id, chain)].get_sequence_id(), protein_map[protein_id].structure_annotations[(structure_id, chain)].get_coverage())
            else:
                structure_quality_measures[(protein_id, structure_id, chain)] = (proteins[protein_id].structure_annotations[(structure_id, chain)].get_sequence_id(), proteins[protein_id].structure_annotations[(structure_id, chain)].get_coverage())
            if not (structure_id, chain) in interface_map_store:
                if proteins is None:
                    interface_map_store[(structure_id, chain)] = complexes[structure_id].interfaces
                else:
                    interface_map_store[(structure_id, chain)] = proteins.get_interfaces(structure_id)

            if proteins is None:
                backmap_store[(protein_id, structure_id, chain)] = protein_map[protein_id].get_backmap(structure_id, chain)
            else:
                backmap_store[(protein_id, structure_id, chain)] = proteins.get_backmap(protein_id, structure_id, chain)


    if proteins is None:
        for (structure_id, chain) in mapped_positions_map:
            if not structure_id in complex_store:
                continue
            if not (structure_id, chain) in mapped_positions:
                mapped_positions[(structure_id, chain)] = {}
                for res in mapped_positions_map[(structure_id, chain)]:
                    mapped_positions[(structure_id, chain)][res] = mapped_positions_map[(structure_id, chain)][res]

    else:
        for (structure_id, chain) in proteins.structures:
            if not structure_id in complex_store:
                continue
            if not (structure_id, chain) in mapped_positions:
                mapped_positions[(structure_id, chain)] = {}
                if config.verbosity >= 7:
                    print(f'{(structure_id, chain)} {structures[(structure_id, chain)].residues}')
                for res in structures[(structure_id, chain)].residues:
                    mapped_positions[(structure_id, chain)][res] = structures[(structure_id, chain)].get_mapped_positions(res, proteins)

        if config.verbosity >= 6:
            print(f'\nMapped positions in generate_classification_dump:\n{mapped_positions}\nState of structures:\n{proteins.structures}\nState of complex_store:\n{complex_store}')

    if unpacked:
        return (config, complex_store, structure_quality_measures, mapped_positions, interface_map_store, backmap_store)
    store_reference =  ray.put((config, pack((complex_store, structure_quality_measures, mapped_positions, interface_map_store, backmap_store))))
    del complex_store
    del structure_quality_measures
    del mapped_positions
    del interface_map_store
    del backmap_store

    return store_reference

def process_classification_outputs(outs, proteins, recommended_complexes):
    for protein_id, aggregated_interfaces, pos_outs in outs:
        proteins.set_aggregated_interface_map(protein_id, aggregated_interfaces)
        if protein_id not in recommended_complexes:
            recommended_complexes[protein_id] = {}
        for pos, mapping_results in pos_outs:
            #print(f'In process_classification_outputs: {protein_id} {pos} {mapping_results}')
            proteins.protein_map[protein_id].positions[pos].mappings = mappings_package.Mappings(raw_results=mapping_results)
            recommended_complexes[protein_id][pos] = proteins.protein_map[protein_id].positions[pos].mappings.recommendation_order

@ray.remote(max_calls = 1)
def para_mm_lookup(package, store):
    config, n_sub_procs = store
    complex_feature_dict_dict = {}
    for complex_id in package:
        feature_dict, _ = mm_lookup(complex_id, config, n_sub_procs = n_sub_procs)
        complex_feature_dict_dict[complex_id] = feature_dict
    return complex_feature_dict_dict


def classification(proteins, config, background_insert_residues_process, indel_analysis_follow_up=False, custom_structure_annotations=None):

    if config.verbosity >= 2:
        t0 = time.time()

    classification_results = []

    size_sorted = []

    if custom_structure_annotations is None:
        protein_ids = proteins.get_protein_ids()
    else:
        protein_ids = custom_structure_annotations.keys()

    #if config.verbosity >= 2:
    #    t_agg_cm_0 = 0
    #    t_agg_cm_1 = 0

    total_mappings = 0

    #structure_protein_backmap = {}
    #protein_interface_map = {}
    for protein_id in protein_ids:
        if proteins.protein_map[protein_id].stored:
            continue
        if custom_structure_annotations is None:
            annotation_list = proteins.get_protein_annotation_list(protein_id)
        else:
            annotation_list = custom_structure_annotations[protein_id]
        size = len(annotation_list) * len(proteins[protein_id].positions)
        size_sorted.append((protein_id, size))

        total_mappings += size

        #if config.verbosity >= 2:
        #    t_agg_cm_0 += time.time()

        #IAmaps = {}
        #backmaps = {}
        #interface_map = {}
        #for (structure_id, chain) in annotation_list:
        #    IAmaps[(structure_id, chain)] = proteins.get_IAmap(structure_id)
        #    backmaps[structure_id] = proteins.get_backmap(protein_id, structure_id, chain)
        #    interface_map[(structure_id, chain)] = proteins.get_interfaces(structure_id)
        #    if not (structure_id, chain) in structure_protein_backmap:
        #        structure_protein_backmap[(structure_id, chain)] = []
        #    structure_protein_backmap[(structure_id, chain)].append(protein_id)

        #proteins.set_aggregated_interface_map(protein_id, interface_package.calculate_aggregated_interface_map(config, proteins, protein_id, interface_map, backmaps, IAmaps))

        #proteins.set_aggregated_contact_matrix(protein_id, aggregate_IAmaps(IAmaps, backmaps))

        #protein_interface_map[protein_id] = interface_map

        #if config.verbosity >= 2:
        #    t_agg_cm_1 += time.time()

    #if background_insert_residues_process is not None:
    #    background_insert_residues_process.join()

    #insert_interfaces(proteins, config)

    #if config.verbosity >= 2:
    #    print(f'Time for calculating aggregated contact matrices: {t_agg_cm_1-t_agg_cm_0}')


    if config.verbosity >= 3:
        print('Total mappings:', total_mappings)

    if config.verbosity >= 5:
        if len(proteins.structures.keys()) <= 1000:
            print(proteins.structures.keys())
        else:
            print(f'Size of protein.structures: {len(proteins.structures)}')

    size_sorted = sorted(size_sorted, reverse=True, key=lambda x: x[1])

    packages = []
    package_counter = 0

    para_mappings_threshold = 20000
    para = (total_mappings > 20000)
    if config.proc_n == 1:
        para = False
    nested_procs = min([config.proc_n, max([len(size_sorted) // 10, 1])])
    #nested_procs = config.proc_n
    n_main_procs = config.proc_n // nested_procs
    if (config.proc_n % nested_procs) != 0:
        n_main_procs += 1
    nested_para = (total_mappings > para_mappings_threshold * nested_procs) and (len(size_sorted) >= n_main_procs)
    #nested_para = True

    max_package_size = min([50000, total_mappings // (config.proc_n * 4)])
    package_size = 0
    package = []

    send_packages = 0
    config_store = None
    recommended_complexes = {}

    if config.verbosity >= 3:
        print(f'Begin of aggregation process, para: {para}, nested_para: {nested_para}')

    if para:
        if not nested_para:

            classification_dump = generate_classification_dump(size_sorted, config, proteins = proteins)

        else:
            classification_dump = None
            t_big_put_0 = time.time()

            #proteins_in_store = ray.put((proteins, config))
            config_store = ray.put(config)

            t_big_put_1 = time.time()

            if config.verbosity >= 5:
                print('Time for putting proteins into store:', (t_big_put_1 - t_big_put_0))

            n_nested_proc_mappings = total_mappings // nested_procs

            nested_proc_size_sorted = []
            nested_proc_size_sorted_size = 0

            nested_protein_map = {}
            nested_structures = {}
            nested_complexes = {}
            nested_mapped_positions_map = {}

            queue = Queue(maxsize = (config.proc_n * 4))
            nested_process_ids = []

            truly_started = 0

            t_nest_0 = 0
            t_nest_1 = 0

            for protein_id, size in size_sorted:
                t_nest_0 += time.time()

                nested_proc_size_sorted.append((protein_id, size))
                nested_proc_size_sorted_size += size

                nested_protein_map[protein_id] = proteins.protein_map[protein_id]
                for (pdb, chain) in nested_protein_map[protein_id].get_annotation_list():
                    if (pdb, chain) in nested_structures:
                        continue
                    if (pdb, chain) not in proteins.structures:
                        if (pdb, chain) in config.removed_structures:
                            config.errorlog.add_warning(f'Structure {pdb}, {chain} not in structure dict and got removed in alignment, additional info: {protein_id}')
                        else:
                            config.errorlog.add_warning(f'Structure {pdb}, {chain} not in structure dict, additional info: {protein_id}')
                        continue
                    nested_structures[(pdb, chain)] = proteins.structures[(pdb, chain)]
                    if pdb in nested_complexes:
                        continue
                    try:
                        nested_complexes[pdb] = proteins.complexes[pdb]
                    except:
                        config.errorlog.add_warning(f'{pdb} not in proteins.complexes, additional_info: {chain} {protein_id}')

                for (structure_id, chain) in proteins.structures:
                    if not structure_id in nested_complexes:
                        continue
                    if not (structure_id, chain) in nested_mapped_positions_map:
                        nested_mapped_positions_map[(structure_id, chain)] = {}
                        for res in proteins.structures[(structure_id, chain)].residues:
                            nested_mapped_positions_map[(structure_id, chain)][res] = proteins.structures[(structure_id, chain)].get_mapped_positions(res, proteins)

                t_nest_1 += time.time()

                if nested_proc_size_sorted_size >= n_nested_proc_mappings:
                    #start nested main process

                    t_nest_2 = time.time()

                    nested_process_ids.append(nested_classification_main_process.remote(queue, config_store, nested_proc_size_sorted, max_package_size, n_main_procs, pack([nested_protein_map, nested_structures, nested_complexes, nested_mapped_positions_map])))

                    t_nest_3 = time.time()

                    if config.verbosity >= 5:
                        print('Time for preparing nested main process:', (t_nest_1 - t_nest_0), 'Time for starting the nested main process:', (t_nest_3 - t_nest_2))

                    t_nest_0 = 0
                    t_nest_1 = 0

                    truly_started += 1
                    nested_proc_size_sorted_size = 0
                    nested_proc_size_sorted = []
                    nested_protein_map = {}
                    nested_structures = {}
                    nested_complexes = {}
                    nested_mapped_positions_map = {}

            if len(nested_proc_size_sorted) > 0:

                max_package_size = min([50000, n_nested_proc_mappings // (n_main_procs * 4)])

                nested_process_ids.append(nested_classification_main_process.remote(queue, config_store, nested_proc_size_sorted, max_package_size, n_main_procs, pack([nested_protein_map, nested_structures, nested_complexes, nested_mapped_positions_map])))
                truly_started += 1

    else:
        classification_dump = generate_classification_dump(size_sorted, config, proteins = proteins, unpacked = True)

    if not nested_para:
        processed_positions = {}
        package_size, package, send_packages, classification_results, size_sorted, processed_positions = pack_packages(package_size, package, send_packages, classification_results, proteins.protein_map, proteins.structures, config, size_sorted, processed_positions, para, classification_dump, max_package_size, config.proc_n, proteins = proteins)

        if not para and len(package) > 0:
            classification_results.append(para_classify(classification_dump, package))
        elif not para and len(package) == 0:
            return

        if config.verbosity >= 2:
            t11 = time.time()
            print('Time for classification part 1.1:', t11 - t0)

        if config.verbosity >= 2:
            t1 = time.time()
            print('Time for classification part 1:', t1 - t0, para, max_package_size)

        if para:
            max_comp_time = 0
            total_comp_time = 0.
            max_comp_time_pos = None
            amount_of_positions = 0

            total_processed_mappings = 0

            while True:
                ready, not_ready = ray.wait(classification_results)

                if len(ready) > 0:

                    t_integrate_0 = 0.
                    t_integrate_1 = 0.
                    t_integrate_2 = 0.
                    send_packages -= len(ready)

                    t_get_0 = time.time()

                    remote_returns = ray.get(ready)

                    t_get_1 = time.time()

                    for returned_package, size_of_returned_package in remote_returns:
                        (outs, comp_time) = returned_package
                        total_processed_mappings += size_of_returned_package
                        t_integrate_0 += time.time()
                        outs = unpack(outs)
                        t_integrate_1 += time.time()
                        process_classification_outputs(outs, proteins, recommended_complexes)

                        total_comp_time += comp_time
                        amount_of_positions += 1

                        if comp_time > max_comp_time:
                            max_comp_time = comp_time
                            max_comp_time_pos = protein_id
                        t_integrate_2 += time.time()

                    classification_results = not_ready

                    t_pack_0 = time.time()
                    package_size, package, send_packages, classification_results, size_sorted, processed_positions = pack_packages(package_size, package, send_packages, classification_results, proteins.protein_map, proteins.structures, config, size_sorted, processed_positions, para, classification_dump, max_package_size, config.proc_n, proteins = proteins)
                    t_pack_1 = time.time()

                    if config.verbosity >= 5:
                        print('Time for getting:', (t_get_1 - t_get_0))
                        print('Time for unpacking:', (t_integrate_1 - t_integrate_0))
                        print('Time for result integration:', (t_integrate_2 - t_integrate_1))
                        print('Time for packaging new processes:', (t_pack_1 - t_pack_0))
                        print('Progress:', total_processed_mappings, 'of', total_mappings)
                else:
                    classification_results = not_ready

                if len(classification_results) == 0:
                    break
        else:
            gc.collect()
            for outs, comp_time in classification_results:
                process_classification_outputs(outs, proteins, recommended_complexes)

    else:
        if config.verbosity >= 2:
            t1 = time.time()

        if config.verbosity >= 5:
            print('Start collecting results from nested classification main processes:', truly_started)

        if n_main_procs > 1:

            finished = 0
            while finished < truly_started:

                try:
                    ready = queue.get()
                except:
                    continue
                if ready == 'End':
                    finished += 1
                    if config.verbosity >= 5:
                        print('Nested main process finished,', truly_started - finished, 'are left')
                    continue

                remote_returns = ray.get(ready)

                for returned_package, size_of_returned_package in remote_returns:
                    (outs, comp_time) = returned_package

                    outs = unpack(outs)
                    process_classification_outputs(outs, proteins, recommended_complexes)
        else:
            while True:
                ready, not_ready = ray.wait(nested_process_ids)

                if len(ready) > 0:
                    outs_list = ray.get(ready)
                    for outs, _ in outs_list:
                        outs = unpack(outs)
                        process_classification_outputs(outs, proteins, recommended_complexes)
                    nested_process_ids = not_ready
                else:
                    nested_process_ids = not_ready

                if len(nested_process_ids) == 0:
                    break
                

    del config_store
    del classification_dump

    if config.verbosity >= 6:
        print(recommended_complexes)

    mm_db_dict = {}
    successful_mm_annotations = {}
    if config.microminer_exe != '':
        mm_inputs = set()
        processed_complexes = set()
        rec_complex_pos_dict = {}
        for prot_id in recommended_complexes:
            for pos in recommended_complexes[prot_id]:
                if len(recommended_complexes[prot_id][pos]) > 0:
                    recommended_structure_triple = recommended_complexes[prot_id][pos][0][1]
                    complex_id = recommended_structure_triple[0]
                    mm_inputs.add(complex_id)
                    processed_complexes.add(complex_id)
                    if complex_id not in rec_complex_pos_dict:
                        rec_complex_pos_dict[complex_id] = []
                    rec_complex_pos_dict[complex_id].append((prot_id, pos, 0))


        while len(mm_inputs) > 0:
            if config.verbosity >= 2:
                print(f'Calling MicroMiner-lookup for {len(mm_inputs)} Complexes.')
            n_empty = 0
            new_mm_inputs = set()
            n_of_chunks = max([1, config.proc_n//4])
            small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(n_of_chunks, len(mm_inputs))
            real_n_of_chunks = n_of_small_chunks + n_of_big_chunks
            if real_n_of_chunks < config.proc_n//2:
                n_sub_procs = config.proc_n // real_n_of_chunks
            else:
                n_sub_procs = 1

            store = ray.put((config, n_sub_procs))
            para_mm_lookup_process_ids = []
            package = []
            for complex_id in mm_inputs:
                package.append(complex_id)
                if n_of_small_chunks > 0:
                    if len(para_mm_lookup_process_ids) < n_of_small_chunks:
                        if len(package) == small_chunksize:
                            para_mm_lookup_process_ids.append(para_mm_lookup.remote(package, store))
                            package = []
                            continue
                if n_of_big_chunks > 0:
                    if len(package) == big_chunksize:
                        para_mm_lookup_process_ids.append(para_mm_lookup.remote(package, store))
                        package = []

            if len(package) > 0:
                config.errorlog.add_error(f'Leftover mm_inputs: {package}')

            results = ray.get(para_mm_lookup_process_ids)
            for complex_feature_dict_dict in results:

                for complex_id in complex_feature_dict_dict:
                    feature_dict = complex_feature_dict_dict[complex_id]

                    if len(feature_dict) == 0:
                        n_empty += 1
                        for (prot_id, pos, current_try) in rec_complex_pos_dict[complex_id]:
                            stay_in_loop = True
                            next_try = current_try + 1
                            while len(recommended_complexes[prot_id][pos]) > next_try and stay_in_loop:
                                recommended_structure_triple = recommended_complexes[prot_id][pos][next_try][1]
                                complex_id, chain_id, res = recommended_structure_triple

                                if complex_id in mm_db_dict:
                                    if (chain_id, res) in mm_db_dict[complex_id]:
                                        successful_mm_annotations[(prot_id, pos)] = complex_id, chain_id, res
                                        stay_in_loop = False
                                    else:
                                        next_try += 1
                                elif not complex_id in processed_complexes:
                                    new_mm_inputs.add(complex_id)
                                    processed_complexes.add(complex_id)
                                    if complex_id not in rec_complex_pos_dict:
                                        rec_complex_pos_dict[complex_id] = []
                                    rec_complex_pos_dict[complex_id].append((prot_id, pos, next_try))
                                    stay_in_loop = False
                                else:
                                    next_try += 1

                    else:
                        mm_db_dict[complex_id] = feature_dict
                        for (prot_id, pos, current_try) in rec_complex_pos_dict[complex_id]:
                            recommended_structure_triple = recommended_complexes[prot_id][pos][current_try][1]
                            complex_id, chain_id, res = recommended_structure_triple
                            if (chain_id, res) in feature_dict:
                                successful_mm_annotations[(prot_id, pos)] = complex_id, chain_id, res
                            else:
                                next_try = current_try + 1
                                stay_in_loop = True
                                while len(recommended_complexes[prot_id][pos]) > next_try and stay_in_loop:
                                    recommended_structure_triple = recommended_complexes[prot_id][pos][next_try][1]
                                    complex_id, chain_id, res = recommended_structure_triple
                                    if complex_id in mm_db_dict:
                                        successful_mm_annotations[(prot_id, pos)] = complex_id, chain_id, res
                                        stay_in_loop = False
                                    elif not complex_id in processed_complexes:
                                        new_mm_inputs.add(complex_id)
                                        processed_complexes.add(complex_id)
                                        if complex_id not in rec_complex_pos_dict:
                                            rec_complex_pos_dict[complex_id] = []
                                        rec_complex_pos_dict[complex_id].append((prot_id, pos, next_try))
                                        stay_in_loop = False
                                    else:
                                        next_try += 1

            if config.verbosity >= 2:
                print(f'MicroMiner-lookup succesful for {len(mm_inputs) - n_empty} Complexes.')
            
            mm_inputs = new_mm_inputs

    prot_ids = proteins.get_protein_ids()

    for prot_id in prot_ids:
        positions = proteins[prot_id].get_position_ids()
        for pos in positions:
            if (prot_id, pos) in successful_mm_annotations:
                complex_id, chain_id, res = successful_mm_annotations[(prot_id, pos)]
                if (chain_id, res) in mm_db_dict[complex_id]:
                    mm_feat_vec = mm_db_dict[complex_id][(chain_id, res)]
                    proteins.protein_map[prot_id].positions[pos].mappings.microminer_features.set_values(mm_feat_vec)



    if background_insert_residues_process is not None:
        background_insert_residues_process.join()
        background_insert_residues_process.close()
        background_insert_residues_process = None

    if config.compute_ppi:
        insert_interfaces(proteins, config)

    if not indel_analysis_follow_up:
        if config.verbosity >= 3:
            print('Proteins object semi deconstruction')
        #for protein_id in protein_ids:
        #    proteins.remove_protein_annotations(protein_id)
        proteins.semi_deconstruct()


    if config.verbosity >= 2:
        t2 = time.time()
        print('Time for classification part 2:', t2 - t1)
        if para and not nested_para:
            print('Longest computation for:', max_comp_time_pos, 'with:', max_comp_time, 'In total', amount_of_positions, 'proteins', 'Accumulated time:', total_comp_time)