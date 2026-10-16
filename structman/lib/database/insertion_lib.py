import time
import numpy as np
import multiprocessing

from structman.lib.pdbParser import parseLigandDB, getSI, updateLigandDB
from structman.lib.sdsc import interface as interface_package
from structman.lib.sdsc import residue as residue_package
from structman.base_utils.base_utils import pack, add_to_times
from structman.base_utils.config_class import Config
from structman.lib.database.database_core_functions import binningSelect, insert, select, update


def _group_rows(*columns: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Group equal rows of the given int columns.

    Returns the group of every row and the first (lowest) row of every group,
    groups numbered in key order, the first column being the primary key; what
    np.unique(axis=0, return_index=True, return_inverse=True) gives, but
    without its slow sort over void records.
    """
    order = np.lexsort(columns[::-1])
    starts = np.zeros(order.shape[0], dtype=bool)
    starts[0] = True
    for column in columns:
        sorted_column = column[order]
        starts[1:] |= sorted_column[1:] != sorted_column[:-1]
    inverse = np.empty(order.shape[0], dtype=np.int64)
    inverse[order] = np.cumsum(starts) - 1
    # lexsort is stable, so a run starts at its lowest row.
    return inverse, order[starts]


def build_prot_prot_rows(
        prot_a_db_id: int,
        interface_map: list[interface_package.Aggregated_interface]) -> tuple[list[tuple], int]:
    """The Protein_Protein_Interaction rows of one protein, plus the number of
    (interface, interacting protein) pairs for the timing log.

    Per interacting protein B, every recommended (complex, chain_a, chain_b) of
    the interactions with B is counted (score > 0) and its scores summed, over
    all interfaces; the most counted one is stored. This works directly on the
    interaction arrays of the interfaces instead of materialising every
    interaction as a Position_Position_Interaction. Only Pos_Pos_Data still
    needs the nested objects, and it takes them from the first interface that
    interacts with B.

    Results equal those of the former nested-dict walk exactly: ties go to the
    complex it would have met first, and scores are summed in its order. That
    walk visited, per interface and protein B, the interactions grouped by
    position_a in first-appearance order (the iteration order of the nested
    dict), and listed the proteins B in first-appearance order too.
    """
    chain_names: list[str] = []
    chain_lookup: dict[str, int] = {}

    iface_parts, key_parts, score_parts = [], [], []
    complex_parts, chain_a_parts, chain_b_parts = [], [], []

    for iface_number, aggregated_interface in enumerate(interface_map):
        n = aggregated_interface.interaction_count
        if n == 0:
            continue
        if not aggregated_interface.interaction_has_rec[:n].all():
            raise ValueError(f'Interaction without recommended complex in an interface of protein {prot_a_db_id}')

        # Chain indices point into the interface's own chain_ids list; re-intern
        # them so that the grouping below does not depend on all interfaces
        # sharing one list.
        chain_a = aggregated_interface.interaction_rec_chain_a[:n]
        chain_b = aggregated_interface.interaction_rec_chain_b[:n]
        local_ids = np.unique(np.concatenate((chain_a, chain_b)))
        global_ids = np.empty(local_ids.shape[0], dtype=np.int64)
        for i, local_id in enumerate(local_ids.tolist()):
            name = aggregated_interface.chain_ids[local_id]
            if name not in chain_lookup:
                chain_lookup[name] = len(chain_names)
                chain_names.append(name)
            global_ids[i] = chain_lookup[name]

        iface_parts.append(np.full(n, iface_number, dtype=np.int64))
        key_parts.append(aggregated_interface.interaction_key_array[:n])
        score_parts.append(aggregated_interface.interaction_scores[:n])
        complex_parts.append(aggregated_interface.interaction_rec_complex[:n])
        chain_a_parts.append(global_ids[np.searchsorted(local_ids, chain_a)])
        chain_b_parts.append(global_ids[np.searchsorted(local_ids, chain_b)])

    if len(iface_parts) == 0:
        return [], 0

    iface_col = np.concatenate(iface_parts)
    keys = np.concatenate(key_parts)
    prot_b = keys[:, 0]
    rows = np.arange(prot_b.shape[0])

    # Visit order of the old walk. Rows are in interface order already, so the
    # first row of a position_a inside its (interface, protein B) is enough to
    # order the groups, and the row itself orders inside a group.
    pos_a_group, first_of_pos_a = _group_rows(iface_col, prot_b, keys[:, 1])
    visit = np.lexsort((rows, first_of_pos_a[pos_a_group]))

    group_columns = (
        prot_b[visit],
        np.concatenate(complex_parts)[visit],
        np.concatenate(chain_a_parts)[visit],
        np.concatenate(chain_b_parts)[visit])
    scores = np.concatenate(score_parts)[visit]

    # On visit-ordered rows the first row of a group is its first visit, and
    # bincount accumulates in index order, so the sums match the old ones bit
    # for bit.
    group_inverse, first_visit = _group_rows(*group_columns)
    n_groups = first_visit.shape[0]
    positive = scores > 0.
    counts = np.bincount(group_inverse, weights=positive, minlength=n_groups)
    score_sums = np.bincount(group_inverse, weights=np.where(positive, scores, 0.), minlength=n_groups)
    group_keys = np.stack([column[first_visit] for column in group_columns], axis=1)

    # Per protein B the highest count, earliest visit first.
    order = np.lexsort((first_visit, -counts, group_keys[:, 0]))
    sorted_prot_b = group_keys[order, 0]
    run_starts = np.ones(order.shape[0], dtype=bool)
    run_starts[1:] = sorted_prot_b[1:] != sorted_prot_b[:-1]
    winners = order[run_starts]

    # Both winners and these are sorted by protein B.
    b_inverse, first_row_of_b = _group_rows(prot_b)
    total_interactions = int(_group_rows(iface_col, b_inverse)[1].shape[0])

    prot_prot_rows = []
    for b_index in np.argsort(first_row_of_b, kind='stable').tolist():
        group = winners[b_index]
        prot_b_db_id, recommended_complex, chain_a, chain_b = group_keys[group].tolist()
        source_interface = interface_map[int(iface_col[first_row_of_b[b_index]])]
        # The old walk summed onto an int 0, which stayed an int without any positive score.
        ppi_score = float(score_sums[group]) if counts[group] > 0 else 0
        prot_prot_rows.append((
            prot_a_db_id, prot_b_db_id, recommended_complex,
            chain_names[chain_a], chain_names[chain_b], ppi_score,
            pack(source_interface.nested_interactions_for(prot_b_db_id))))

    return prot_prot_rows, total_interactions


def remote_insert_interfaces(
        prot_a_db_id: int,
        interface_map: list[interface_package.Aggregated_interface],
        config: Config,
        db_lock=None) -> list[float]:
    interface_values = []
    prot_db_ids = set()

    if config.verbosity >= 5:
        config.logger.info(f'In remote_insert_interfaces: {prot_a_db_id} {len(interface_map)=}')

    times = []
    ti = time.time()

    prot_db_ids.add(prot_a_db_id)

    hashes_to_insert: dict[int, str] = {}
    interfaces_without_row: set[int] = set()

    for aggregated_interface in interface_map:
        
        structure_recommendation = f'{aggregated_interface.recommended_complex},{aggregated_interface.chain},{aggregated_interface.interacting_chain}'
        if config.verbosity >= 7:
            config.logger.info(f'In insert_interfaces: {prot_a_db_id} {structure_recommendation=} {list(aggregated_interface.interacting_protein_ids())=}')
            aggregated_interface.print_interface()

        #structure_ids_for_residue_retrieval[proteins.structures[(aggregated_interface.recommended_complex, aggregated_interface.chain)].database_id] = (aggregated_interface.recommended_complex, aggregated_interface.chain)

        interface_hash_tuple = aggregated_interface.get_interface_hash()

        if aggregated_interface.database_id is None:
            # The hash is the only handle on the row that is about to be
            # inserted blank. Two interfaces sharing one would take the same
            # Interface_Id below, merging their residues into one row and
            # orphaning the other; that used to pass unnoticed.
            if interface_hash_tuple in hashes_to_insert:
                interfaces_without_row.add(id(aggregated_interface))
                config.errorlog.add_warning(
                    f'Interface hash collision for protein {prot_a_db_id}: {interface_hash_tuple} given by '
                    f'{structure_recommendation} ({len(aggregated_interface.positions)} positions) and by '
                    f'{hashes_to_insert[interface_hash_tuple]}, which keeps the row. The interface residues '
                    f'of {structure_recommendation} are not stored.')
            else:
                hashes_to_insert[interface_hash_tuple] = structure_recommendation
                interface_values.append((prot_a_db_id, structure_recommendation, interface_hash_tuple))

    prot_prot_values, total_interactions = build_prot_prot_rows(prot_a_db_id, interface_map)

    ti = add_to_times(times, ti, f'build prot_prot_values {prot_a_db_id=} {len(interface_map)=} {len(prot_prot_values)=} {total_interactions=}')

    insert('Interface', ['Protein', 'Structure_Recommendation', 'Hash'], interface_values, config, locked=True)
    
    ti = add_to_times(times, ti, f'insert {len(interface_values)=} {len(prot_db_ids)=}')

    results = binningSelect(prot_db_ids, ['Protein', 'Interface_Id', 'Hash'], 'Interface', config, locked=True, db_lock=db_lock)

    ti = add_to_times(times, ti, f'select interfaces {len(results)}')

    blank_interfaces = {}
    for row in results:
        prot_db_id = row[0]
        interface_db_id = row[1]

        interface_hash_tuple = row[2]
        if prot_db_id not in blank_interfaces:
            blank_interfaces[prot_db_id] = {}

        blank_interfaces[prot_db_id][interface_hash_tuple] = interface_db_id

    ti = add_to_times(times, ti, f'build blank interfaces {len(blank_interfaces)=}')

    if config.verbosity >= 6:
        if len(prot_prot_values) < 1000:
            config.logger.info(f'In insert_interfaces: contents of prot_prot_values: {prot_prot_values=}')
        else:
            config.logger.info(f'In insert_interfaces: contents of prot_prot_values: {len(prot_prot_values)=}')

    if len(prot_prot_values) > 0:
        insert(f'Protein_Protein_Interaction_{prot_a_db_id % 20}', ['Protein_A', 'Protein_B', 'Complex', 'Chain_A', 'Chain_B', 'Interaction_Score', 'Pos_Pos_Data'], prot_prot_values, config, locked=True, buffer_factor=30.)

    ti = add_to_times(times, ti, f'insert Protein_Protein_Interactions {len(prot_prot_values)=}')

    position_interface_values = []

    if prot_a_db_id not in blank_interfaces:
        if config.verbosity >= 5:
            config.logger.info(f'Protein ({prot_a_db_id}) omitted in position_interface insert')
        return times
    for aggregated_interface in interface_map:
        if aggregated_interface.database_id is not None:
            continue

        # Lost a hash collision above, so it has no row of its own; its
        # residues would land on the interface that took the row.
        if id(aggregated_interface) in interfaces_without_row:
            continue

        interface_db_id = blank_interfaces[prot_a_db_id][aggregated_interface.get_interface_hash()]

        for position_db_id in aggregated_interface.positions:
            structure_id, chain, res = aggregated_interface.positions[position_db_id]
            position_interface_values.append((interface_db_id, position_db_id, f'{structure_id}:{chain} {res}'))

    ti = add_to_times(times, ti)

    insert('RS_Position_Interface', ['Interface', 'Position', 'Recommended_Residue'], position_interface_values, config, locked=True)

    ti = add_to_times(times, ti)

    return times

def remote_insert_interface_residues(
        interfaces: dict[str, dict[str, dict[str, interface_package.Interface]]],
        structure_ids: dict[str, dict[str, int]],
        config: Config,
        db_lock: None | str =None):
    structure_interface_values = []
 
    for structure_id in interfaces:
        if config.verbosity >= 5:
            config.logger.info(f'In insert interface residues: {structure_id}, {len(interfaces[structure_id])=}')
        for chain in interfaces[structure_id]:

            try:
                structure_db_id = structure_ids[structure_id][chain]
            except KeyError:
                config.errorlog.add_warning(f'Structure missing: {structure_id} {chain=}')
                continue

            structure_interface_values.append((structure_db_id, pack(interfaces[structure_id][chain])))

    insert('RS_Structure_Interface', ['Structure', 'Interacting_Residue_Data'], structure_interface_values, config, locked=True, db_lock=db_lock)


def background_insert_residues(value_lists, config: Config, locked=False, db_lock = None):

    for db_nr, values in enumerate(value_lists):
        columns = ['Structure_Id', 'Residue_Data', 'Interacting_Residue_Data', 'Simple_Class_Data']
        update(config, f'Structure_{db_nr}', columns, values, locked = locked, db_lock = db_lock)


def remote_insertResidues(
        structural_analysis_list: list[tuple[str, dict[str, list[np.ndarray[np.float32]]]]],
        structure_ids: dict[str: dict[str, int]],
        interfaces: dict[str, dict[str, list[tuple[str, np.ndarray[np.float32]]]]],
        config: Config,
        locked: bool = True,
        db_lock = None) -> None:

    if config.verbosity >= 4:
        t0 = time.time()

    if len(structural_analysis_list) == 0:
        if config.verbosity >= 4:
            config.logger.info('insertResidues was called with an empty structural_analysis dict')
        return

    value_lists = [[] for _ in range(20)]

    for structure_id, structural_analysis in structural_analysis_list:
        for chain in structural_analysis:
            analysis_list: list[np.ndarray[np.float32]] = structural_analysis[chain]
            try:
                internal_s_id: int = structure_ids[structure_id][chain]
            except KeyError as e:
                config.errorlog.add_warning(f'{structure_id} {chain} not in structure ids ({e})')
                continue

            if internal_s_id is None:
                config.errorlog.add_warning(f'{structure_id} {chain} from structural_analysis not going into the database. s_id was None.')                
                continue

            if config.verbosity >= 6:
                config.logger.info(f'{structure_id} {chain} from structural_analysis going into the database: {internal_s_id=} {len(analysis_list)=}')
            if config.verbosity >= 7:
                config.logger.info(f'{structure_id} {chain} from structural_analysis going into the database: {internal_s_id=} {analysis_list=}')
            
            try:
                # Explicit float32 to match Residue.encode_to_numpy(); an empty
                # analysis_list would otherwise fall back to float64.
                residue_data: np.ndarray[np.float32] = np.array(analysis_list, dtype=np.float32)
            except ValueError as err:
                for enc in analysis_list:
                    print(f'{len(enc)=}')
                raise ValueError(f'{analysis_list=} {err=}')

            # Derived before residue_data is replaced by its packed form. The
            # thresholds in effect now are baked into the location half of every
            # code, exactly as they are for the codes the aggregation path computes
            # on read; a database annotated under different thresholds would report
            # stale Core/Buried/Surface labels for non-interacting residues.
            simple_class_data: bytes = pack(residue_package.simple_class_codes_for_matrix(
                residue_data, config.surface_threshold, config.buried_threshold))

            residue_data: bytes = pack(residue_data)

            if_data = []
            if structure_id in interfaces:
                if chain in interfaces[structure_id]:
                    if_data = interfaces[structure_id][chain]

            if_data = pack(if_data)

            s_id = internal_s_id // 100
            db_nr = internal_s_id % 100

            value_lists[db_nr].append((s_id, residue_data, if_data, simple_class_data))

    if config.verbosity >= 4:
        t1 = time.time()
        config.logger.info(f'Time for insertResidues part 1: {t1 - t0} {locked=} {db_lock=} {len(value_lists[0])=}')

    background_insert_residues(value_lists, config, locked=True, db_lock=db_lock)
    


def remote_insert_classifications(value_lists, config: Config):
    if config.verbosity >= 2:
        t1 = time.time()

    for db_nr, values in enumerate(value_lists):
        if len(values) == 0:
            continue
        table = f'Position_{db_nr}'
        columns = ['Position_Id', 'Recommended_Structure_Data', 'Position_Data']

        update(config, table, columns, values, locked=True)

        if config.verbosity >= 4:
            t2 = time.time()
            config.logger.info(f'remote_insertClassifications part 1: {t2 - t1} Update Classifications of {len(values)} positions')



def remote_insert_complexes(config: Config, values, compl_interaction_partners, ligand_dump):
    stored_ligands, ligand_db = ligand_dump
    if len(values) > 0:
        insert('Complex', ['PDB', 'Resolution', 'Chains', 'Homooligomers', 'Chaintypes'], values, config, locked=True)

        lig_values = []
        update_ligands = []
        new_ligands = set()
        for value_list in values:
            structure_id = value_list[0]
            interaction_partners = compl_interaction_partners[structure_id]
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
                    (smiles, inchi) = getSI(structure_id, name, res, chain, config.pdb_path, config)
                    update_ligands.append((name, smiles, inchi))
                new_ligands.add(name)
                lig_values.append((name, smiles, inchi))

    if len(update_ligands) > 0:
        updateLigandDB(update_ligands, config.smiles_path, config.inchi_path)

    if len(lig_values) > 0:
        try:
            insert('Ligand', ['Name', 'Smiles', 'Inchi'], lig_values, config, locked=True)
        except BaseException as err:
            # There is a (or more than one) ligand that results in an sql error 1241. We need to find it, to figure out why.
            config.errorlog.add_error(f'Error trying inserting ligands {err=}\n {lig_values=}')

    stored_ligands = {}
    results = select(config, ['Ligand_Id', 'Name'], 'Ligand', locked=True)
    for row in results:
        ligand_name = row[1]
        stored_ligands[ligand_name] = row[0]

    values = []
    results = select(config, ['Complex_Id', 'PDB', 'Resolution'], 'Complex', locked=True)

    for row in results:
        structure_id = row[1]
        if structure_id not in compl_interaction_partners:
            continue
        interaction_partners = compl_interaction_partners[structure_id]

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
        insert('RS_Ligand_Structure', ['Ligand', 'Complex', 'Chain', 'Residue'], values, config, locked=True)

    return results

    
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