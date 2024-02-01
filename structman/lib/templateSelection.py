import multiprocessing

from structman.lib import pdbParser
from structman.lib.sdsc import structure as structure_package
from structman.lib.sdsc import complex as complex_package


# Called by serializedPipeline
def selectTemplates(structures, pdb_path):
    if len(structures) == 0:
        return {}
    filtered_structures = {}
    intermediate_results = {}
    for (pdb_id, chain) in structures:
        if pdb_id in intermediate_results:
            resolution, homomer_dict = intermediate_results[pdb_id]
        else:
            resolution, homomer_dict = pdbParser.getInfo(pdb_id, pdb_path)
            intermediate_results[pdb_id] = resolution, homomer_dict
        if resolution is None:
            continue
        filtered_structures[(pdb_id, chain)] = structures[(pdb_id, chain)]
        filtered_structures[(pdb_id, chain)]['Resolution'] = resolution
        filtered_structures[(pdb_id, chain)]['Oligo'] = homomer_dict[chain]

    return filtered_structures


def paraGetInfo(lock, input_queue, out_queue, pdb_path):
    with lock:
        input_queue.put(None)
    while True:
        with lock:
            inp = input_queue.get()
        if inp is None:
            return
        pdb_id = inp

        resolution, homomer_dict = pdbParser.getInfo(pdb_id, pdb_path)

        with lock:
            out_queue.put((pdb_id, resolution, homomer_dict))


# looks like an artifact
"""
# called by serializedPipeline
def filterRawStructureMap(raw_structure_map, pdb_ids, pdb_path, option_res_thresh, n_processes, proteins, manager, lock):

    input_queue = manager.Queue()
    out_queue = manager.Queue()

    info_map = {}

    for pdb_id in pdb_ids:
        input_queue.put(pdb_id)

    processes = {}
    for i in range(1, n_processes + 1):

        p = multiprocessing.Process(target=paraGetInfo, args=(lock, input_queue, out_queue, pdb_path))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()

    del input_queue

    out_queue.put(None)
    while True:
        out = out_queue.get()
        if out is None:
            break
        (pdb_id, resolution, homomer_dict) = out
        info_map[pdb_id] = (resolution, homomer_dict)

    del out_queue

    u_acs = proteins.get_protein_ids()
    structure_list = proteins.get_structure_list()
    complex_list = proteins.get_complex_list()

    for u_ac in raw_structure_map:
        for pdb_id, chain in raw_structure_map[u_ac]:
            if pdb_id not in info_map:
                continue
            resolution, homomer_dict = info_map[pdb_id]
            if resolution is None:
                continue
            if resolution > option_res_thresh:
                continue
            oligo = raw_structure_map[u_ac][(pdb_id, chain)]['Oligo']
            struct_anno = structure_package.StructureAnnotation(u_ac, pdb_id, chain)
            proteins.add_annotation(u_ac, pdb_id, chain, struct_anno)

            if not (pdb_id, chain) in structure_list:
                struct = structure_package.Structure(pdb_id, chain, oligo=oligo, mapped_proteins=[u_ac])
                proteins.add_structure(pdb_id, chain, struct)
            else:
                proteins.add_mapping_to_structure(pdb_id, chain, u_ac)

            if pdb_id not in complex_list:
                compl = complex_package.Complex(pdb_id, resolution, homomers=homomer_dict)
                proteins.add_complex(pdb_id, compl)
    return
"""
