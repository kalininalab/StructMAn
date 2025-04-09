from structman.lib import pdbParser


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

