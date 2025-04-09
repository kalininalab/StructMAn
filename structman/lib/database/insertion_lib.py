import time
import multiprocessing
from filelock import FileLock

from structman.lib.pdbParser import parseLigandDB, getSI, updateLigandDB
from structman.lib.database.retrieval import getStoredResidues
from structman.lib.sdsc import structure as structure_package
from structman.lib.sdsc import residue as residue_package
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
            print(f'In insert_interfaces: {protein_id} {len(proteins.protein_map[protein_id].aggregated_interface_map)}')
        prot_a_db_id = proteins.protein_map[protein_id].database_id

        pos_pos_map = {}

        for aggregated_interface in proteins.protein_map[protein_id].aggregated_interface_map:
            prot_db_ids.add(prot_a_db_id)
            structure_recommendation = f'{aggregated_interface.recommended_complex},{aggregated_interface.chain},{aggregated_interface.interacting_chain}'
            if config.verbosity >= 6:
                print(f'In insert_interfaces: {protein_id} {structure_recommendation} {len(aggregated_interface.pos_pos_interactions)}')
                aggregated_interface.print_interface(verbosity = config.verbosity)

            #structure_ids_for_residue_retrieval[proteins.structures[(aggregated_interface.recommended_complex, aggregated_interface.chain)].database_id] = (aggregated_interface.recommended_complex, aggregated_interface.chain)

            interface_hash_tuple = aggregated_interface.get_interface_hash()

            if aggregated_interface.database_id is None:
                mi, ma, me, le = interface_hash_tuple
                interface_values.append((prot_a_db_id, structure_recommendation, mi, ma, me, le))

            for pos_a in aggregated_interface.pos_pos_interactions:
                try:
                    pos_a_db_id = proteins.protein_map[protein_id].positions[pos_a].database_id
                except:
                    config.errorlog.add_warning(f'Position database id not found: {protein_id} {pos_a}')
                    continue

                if pos_a_db_id not in pos_pos_map:
                    pos_pos_map[pos_a_db_id] = {}

                for i_prot in aggregated_interface.pos_pos_interactions[pos_a]:
                    if config.verbosity >= 7:
                        print(f'Processing {i_prot=}({len(aggregated_interface.pos_pos_interactions[pos_a][i_prot])}) in Interface {protein_id} {interface_hash_tuple=}')

                    for i_pos in aggregated_interface.pos_pos_interactions[pos_a][i_prot]:
                        pos_pos_interaction = aggregated_interface.pos_pos_interactions[pos_a][i_prot][i_pos]

                        if pos_pos_interaction.protein_b is None or pos_pos_interaction.position_b is None:
                            config.errorlog.add_warning(f'Faulty pos_pos_interaction object: {protein_id} {pos_a} {i_prot} {i_pos} {pos_pos_interaction.protein_b=} {pos_pos_interaction.position_b=}')
                            continue
                        pos_b_db_id = proteins.protein_map[pos_pos_interaction.protein_b].positions[pos_pos_interaction.position_b].database_id
                        
                        (chain_a, res_a, chain_b, res_b) = pos_pos_interaction.recommended_interaction


                        #structure_ids_for_residue_retrieval[proteins.structures[(pos_pos_interaction.recommended_complex, chain_a)].database_id] = (pos_pos_interaction.recommended_complex, chain_a)
                        #structure_ids_for_residue_retrieval[proteins.structures[(pos_pos_interaction.recommended_complex, chain_b)].database_id] = (pos_pos_interaction.recommended_complex, chain_b)
                        if pos_b_db_id not in pos_pos_map[pos_a_db_id]:
                            pos_pos_pre_values.append((pos_a_db_id, pos_b_db_id, pos_pos_interaction.recommended_complex, chain_a, chain_b, res_a, res_b, pos_pos_interaction.interaction_score))
                            pos_pos_map[pos_a_db_id][pos_b_db_id] = (len(pos_pos_pre_values) - 1, pos_pos_interaction.interaction_score, (pos_pos_interaction.recommended_complex, chain_a, chain_b))

                            if not prot_a_db_id in prot_prot_map:
                                prot_prot_map[prot_a_db_id] = {}

                            prot_b_db_id = proteins.protein_map[pos_pos_interaction.protein_b].database_id

                            #if not prot_b_db_id in prot_prot_map:
                            #    pos_pos_pre_values.append((pos_b_db_id, pos_a_db_id, pos_pos_interaction.recommended_complex, chain_b, chain_a, res_b, res_a, pos_pos_interaction.interaction_score))

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

                        elif pos_pos_interaction.interaction_score > pos_pos_map[pos_a_db_id][pos_b_db_id][1]:
                            pos_pos_pre_values[pos_pos_map[pos_a_db_id][pos_b_db_id][0]] = (pos_a_db_id, pos_b_db_id, pos_pos_interaction.recommended_complex, chain_a, chain_b, res_a, res_b, pos_pos_interaction.interaction_score)
                            prot_prot_map[prot_a_db_id][prot_b_db_id][3] += (pos_pos_interaction.interaction_score - pos_pos_map[pos_a_db_id][pos_b_db_id][1])
                            
                            if not (pos_pos_interaction.recommended_complex, chain_a, chain_b) in prot_prot_map[prot_a_db_id][prot_b_db_id][2]:
                                prot_prot_map[prot_a_db_id][prot_b_db_id][2][(pos_pos_interaction.recommended_complex, chain_a, chain_b)] = 0
                            prot_prot_map[prot_a_db_id][prot_b_db_id][2][(pos_pos_interaction.recommended_complex, chain_a, chain_b)] += 1
                            
                            try:
                                prot_prot_map[prot_a_db_id][prot_b_db_id][2][pos_pos_map[pos_a_db_id][pos_b_db_id][2]] -= 1
                            except KeyError as e:
                                config.errorlog.add_warning(f'Unexpected KeyError: {prot_a_db_id} {prot_b_db_id}\n{e}')
                            pos_pos_map[pos_a_db_id][pos_b_db_id] = (pos_pos_map[pos_a_db_id][pos_b_db_id][0], pos_pos_interaction.interaction_score, (pos_pos_interaction.recommended_complex, chain_a, chain_b))


    for structure_id in proteins.complexes:
        for chain in proteins.complexes[structure_id].chains:
            if proteins.complexes[structure_id].chains[chain] == 'Protein' and structure_id in proteins.structures:
                if chain in proteins.structures[structure_id]:
                    if len(proteins.structures[structure_id][chain].residues) > 0: #If there are residues inside the structure datastructure ...
                        break_it = False
                        for res in proteins.structures[structure_id][chain].residues.get_keys():
                            residue = proteins.structures[structure_id][chain].residues.get_item(res)
                            if residue is None:
                                continue
                            if isinstance(residue, tuple):
                                break_it = True
                                break
                            elif residue.database_id != None: #look if the first residue has a database id ...
                                break_it = True
                                break
                        if break_it: #If yes, don't queue it for retrieval
                            continue
                    structure_ids_for_residue_retrieval[proteins.structures[structure_id][chain].database_id] = (structure_id, chain)

    getStoredResidues(proteins, config, custom_ids = structure_ids_for_residue_retrieval, retrieve_only_db_ids = True, exclude_interacting_chains=True)

    for (pos_a_db_id, pos_b_db_id, recommended_complex, chain_a, chain_b, res_a, res_b, interaction_score) in pos_pos_pre_values:
        if config.verbosity >= 7:
            print((pos_a_db_id, pos_b_db_id, recommended_complex, chain_a, chain_b, res_a, res_b, interaction_score))
        if not recommended_complex in proteins.structures:
            continue

        if not proteins.structures[recommended_complex][chain_a].residues.contains(res_a):
            continue

        if isinstance(proteins.structures[recommended_complex][chain_a].residues.get_item(res_a), tuple):
            res_a_db_id = proteins.structures[recommended_complex][chain_a].residues.get_item(res_a)[1]
        else:
            res_a_db_id = proteins.structures[recommended_complex][chain_a].residues.get_item(res_a).database_id

        if not chain_b in proteins.structures[recommended_complex]:
            continue

        if not proteins.structures[recommended_complex][chain_b].residues.contains(res_b):
            continue

        if isinstance(proteins.structures[recommended_complex][chain_b].residues.get_item(res_b), tuple):
            res_b_db_id = proteins.structures[recommended_complex][chain_b].residues.get_item(res_b)[1]
        else:
            res_b_db_id = proteins.structures[recommended_complex][chain_b].residues.get_item(res_b).database_id
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
        if len(prot_prot_map) < 1000:
            print(f'In insert_interfaces: contents of prot_prot_map before transformation to prot_prot_values:\n{prot_prot_map}')
        else:
            print(f'In insert_interfaces: size of prot_prot_map before transformation to prot_prot_values:{len(prot_prot_map)}')

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
            elif prot_b_db_id not in blank_interfaces:
                config.errorlog.add_warning(f'{prot_b_db_id} not in blank_interfaces')
                continue
            elif len(blank_interfaces[prot_b_db_id]) == 1:
                interface_b_db_id = list(blank_interfaces[prot_b_db_id].values())[0]
            else:
                #TODO retrace correct interface
                interface_b_db_id = list(blank_interfaces[prot_b_db_id].values())[0]

            prot_prot_values.append((interface_a_db_id, interface_b_db_id, proteins.complexes[recommended_complex].database_id, chain_a, chain_b, ppi_score))

            if prot_b_db_id not in prot_prot_map:
                prot_prot_values.append((interface_b_db_id, interface_a_db_id, proteins.complexes[recommended_complex].database_id, chain_b, chain_a, ppi_score))

    if config.verbosity >= 7:
        print(f'In insert_interfaces: contents of prot_prot_values:\n{prot_prot_values}')

    position_interface_values = []

    for protein_id in proteins.protein_map:
        prot_a_db_id = proteins.protein_map[protein_id].database_id
        if prot_a_db_id not in blank_interfaces:
            if config.verbosity >= 5:
                print(f'Protein {protein_id} ({prot_a_db_id}) omitted in position_interface insert')
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

    
def insert_interface_residues(proteins, config):
    residue_interface_values = []

    print_count = 0
 
    for structure_id in proteins.complexes:
        #print(structure_id, proteins.complexes[structure_id].interfaces)
        for chain in proteins.complexes[structure_id].interfaces:
            for i_chain in proteins.complexes[structure_id].interfaces[chain]:
                if proteins.complexes[structure_id].interfaces[chain][i_chain].stored:
                    continue
                try:
                    structure_db_id = proteins.structures[structure_id][chain].database_id
                except:
                    config.errorlog.add_warning(f'Structure missing: {structure_id} {chain} {i_chain}')
                    continue
                for res in proteins.complexes[structure_id].interfaces[chain][i_chain].residues:

                    if not proteins.structures[structure_id][chain].residues.contains(res):
                        if config.verbosity >= 4:
                            if print_count < 100:
                                print(f'Residue {res} filtered from {structure_id} {chain} {i_chain} in insert_interfaces: not in structures object')
                                print_count += 1
                        continue #Happens for ligands declared as part of the chain
                    res_db_id = proteins.structures[structure_id][chain].residues.get_item(res).database_id
                    if res_db_id is None: 
                        if config.verbosity >= 4:
                            print(f'Residue {res} filtered from {structure_id} {chain} {i_chain} in insert_interfaces: db_id is None')
                        continue
                    if res in proteins.complexes[structure_id].interfaces[chain][i_chain].interactions:
                        i_res_strings = []
                        for i_res in proteins.complexes[structure_id].interfaces[chain][i_chain].interactions[res]:
                            score = proteins.complexes[structure_id].interfaces[chain][i_chain].interactions[res][i_res]
                            i_res_strings.append(f'{i_chain},{i_res},{score}')
                        i_res_string = ';'.join(i_res_strings)
                    else:
                        i_res_string = f'{i_chain},,'
                    residue_interface_values.append((structure_db_id, res_db_id, pack(i_res_string)))

    insert('RS_Residue_Interface', ['Structure', 'Residue', 'Interacting_Residue'], residue_interface_values, config)


def background_insert_residues(values, config):

    columns = ['Structure', 'Number', 'Residue_Data']
    insert('Residue', columns, values, config)


def insertResidues(structural_analysis, proteins, config):

    if config.verbosity >= 2:
        t0 = time.time()

    values = []

    if len(structural_analysis) == 0:
        if config.verbosity >= 4:
            print('insertResidues was called with an empty structural_analysis dict')
        return

    t_0 = 0.
    t_1 = 0.
    t_2 = 0.


    for pdb_id in structural_analysis:
        for chain in structural_analysis[pdb_id]:
            
            if not proteins.contains_structure(pdb_id, chain):
                if config.verbosity >= 6:
                    config.errorlog.add_warning(f'{pdb_id} {chain} from structural_analysis not going into the database. contains_structure was False.')
                continue
            if proteins.contains_structure(pdb_id, chain) and proteins.is_structure_stored(pdb_id, chain) and not proteins.structures[pdb_id][chain].new_interacting_chain:  # all residues belonging to stored structures must not inserted twice
                if config.verbosity >= 6:
                    print(f'{pdb_id} {chain} from structural_analysis not going into the database. Is already stored.')
                continue

            analysis_map = structural_analysis[pdb_id][chain]
            
            s_id = proteins.get_structure_db_id(pdb_id, chain)
            if s_id is None:
                if config.verbosity >= 6:
                    config.errorlog.add_warning(f'{pdb_id} {chain} from structural_analysis not going into the database. s_id was None.')
                continue

            if config.verbosity >= 6:
                print(f'{pdb_id} {chain} from structural_analysis going into the database: {len(analysis_map)}')

            max_prints_per_prot = 10
            n_prints = 0

            for res_id in analysis_map.get_keys():
                t_0 += time.time()
                residue: residue_package.Residue = analysis_map.get_item(res_id)
                if residue is None:
                    t_2 += time.time()
                    continue

                packed_res_info = pack(residue.get_res_info())
                
                t_1 += time.time()

                values.append([s_id, res_id, packed_res_info])
                
                t_2 += time.time()

                if config.verbosity >= 6:
                    if n_prints < max_prints_per_prot:
                        print(f'Calling structures.add_residue in insertResidues: {pdb_id} {chain} {res_id}')
                        n_prints += 1
                proteins.structures[pdb_id][chain].add_residue(res_id, residue)
                t_13 += time.time()

    if config.verbosity >= 2:
        t1 = time.time()
        print(f'Time for insertResidues part 1: {t1 - t0} {config.low_mem_system}')
        print('Individual parts: 1 -', (t_1 - t_0), '2 -', (t_2 - t_1))

    if config.low_mem_system:
        process = None
        background_insert_residues(values, config)
    else:
        process = multiprocessing.Process(target=background_insert_residues, args=(values, config))
        process.start()

    return process

def remote_insertResidues(
        structural_analysis_list: list[tuple[str, dict[str, residue_package.Residue_Map]]],
        structure_ids: dict[str: dict[str, int]],
        config,
        locked: bool = True) -> None:

    if config.verbosity >= 4:
        t0 = time.time()

    values: list[tuple[int, int | str, bytes]] = []

    if len(structural_analysis_list) == 0:
        if config.verbosity >= 4:
            print('insertResidues was called with an empty structural_analysis dict')
        return

    t_0 = 0.
    t_1 = 0.
    t_2 = 0.

    for pdb_id, structural_analysis in structural_analysis_list:
        for chain in structural_analysis:
            analysis_map = structural_analysis[chain]
            try:
                s_id: int = structure_ids[pdb_id][chain]
            except KeyError as e:
                config.errorlog.add_warning(f'{pdb_id} {chain} not in structure ids ({e})')
                continue

            if s_id is None:
                config.errorlog.add_warning(f'{pdb_id} {chain} from structural_analysis not going into the database. s_id was None.')                
                continue

            if config.verbosity >= 6:
                print(f'{pdb_id} {chain} from structural_analysis going into the database: {len(analysis_map)}')

            for res_id in analysis_map.get_keys():
                t_0 += time.time()
                residue: residue_package.Residue = analysis_map.get_item(res_id)
                if residue is None:
                    t_2 += time.time()
                    continue

                packed_res_info = pack(residue.get_res_info())
                
                t_1 += time.time()

                values.append([s_id, res_id, packed_res_info])
                
                t_2 += time.time()

    if config.verbosity >= 4:
        t_1 = time.time()
        print(f'Time for insertResidues part 1: {t_2 - t0} {config.low_mem_system=}')
        print('Individual parts: 1 -', (t_1 - t_0), '2 -', (t_2 - t_1))

    if not locked:
        process = multiprocessing.Process(target=background_insert_residues, args=(values, config))
        process.start()
        return process
    else:
        with FileLock('insert_residues_lock.lock'):
            background_insert_residues(values, config)
    return None

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
        positions = proteins[u_ac].get_position_ids()

        for pos in positions:

            if proteins.is_position_stored(u_ac, pos):
                continue

            m = proteins.get_position_database_id(u_ac, pos)

            position = proteins.get_position(u_ac, pos)
            try:
                packed_mappings = position.mappings
            except AttributeError:
                packed_mappings = None

            try:
                packed_mm_features = position.microminer_features
            except AttributeError:
                packed_mm_features = None

            try:
                packed_res_recommends = position.packed_res_recommends
            except AttributeError:
                packed_res_recommends = None

            values.append((m, packed_res_recommends, pack((packed_mappings, packed_mm_features))))
            
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

        resolution = proteins.complexes[pdb_id].resolution
        chains_str = proteins.complexes[pdb_id].getChainStr()
        homomers_str = proteins.complexes[pdb_id].getHomomersStr()
        lig_str = proteins.complexes[pdb_id].getLigProfileStr()
        metal_str = proteins.complexes[pdb_id].getMetalProfileStr()
        ion_str = proteins.complexes[pdb_id].getIonProfileStr()
        cc_str = proteins.complexes[pdb_id].getChainChainProfileStr()

        values.append((pdb_id, resolution, chains_str, homomers_str, lig_str, metal_str, ion_str, cc_str))

        interaction_partners = proteins.get_interaction_partners(pdb_id)

        for iap in interaction_partners:
            ia_type = iap[0]
            if ia_type != "Ligand":
                continue
            name = iap[1]
            if isinstance(name, bytes):
                name = name.decode('ascii')
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
                    smiles = "[%s%s]" % (name[0], name[1:2].lower())
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


# called by indel_analysis
def insert_indel_results(proteins, config):
    table = 'Indel'
    columns = ['Indel_Id', 'Wildtype_Protein', 'Mutant_Protein', 'Indel_Notation', 'Analysis_Results']

    values = []
    for prot_id in proteins.indels:
        for indel_notation in proteins.indels[prot_id]:
            indel_obj = proteins.indels[prot_id][indel_notation]
            wt_prot_id = proteins.get_protein_database_id(proteins.indels[prot_id][indel_notation].wt_prot)
            mut_prot_id = proteins.get_protein_database_id(proteins.indels[prot_id][indel_notation].mut_prot)
            analysis_results = pack((indel_obj.size, indel_obj.delta_delta_classification, indel_obj.wt_aggregates, indel_obj.mut_aggregates,
                                        indel_obj.left_flank_wt_aggregates, indel_obj.left_flank_mut_aggregates,
                                        indel_obj.right_flank_wt_aggregates, indel_obj.right_flank_mut_aggregates))
            values.append((indel_obj.database_id, wt_prot_id, mut_prot_id, indel_notation, analysis_results))
    update(config, table, columns, values)
    return