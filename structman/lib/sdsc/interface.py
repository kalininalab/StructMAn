import ray
from structman.lib.sdsc.sdsc_utils import doomsday_protocol, Slotted_obj
from structman.lib.sdsc.residue import Residue_Map
from structman.lib.database.retrieval_2 import getBackmaps
from structman.base_utils.base_utils import unpack

def calc_interface_annotation_score(seq_id, coverage, interactions):
    total_score = 0.
    for res_a in interactions:
        for res_b in interactions[res_a]:
            total_score += interactions[res_a][res_b]
    return total_score * seq_id * coverage

class Position_Position_Interaction(Slotted_obj):
    __slots__ = ['protein_a', 'position_a', 'protein_b',
                 'position_b', 'recommended_complex', 'recommended_interaction',
                 'interaction_score']
    
    slot_mask = [False, False, False,
                 False, True, True,
                 True]
    
    def __init__(self, protein_a = None, position_a = None, protein_b = None, position_b = None):
        self.protein_a = protein_a
        self.position_a = position_a
        self.protein_b = protein_b
        self.position_b = position_b
        self.interaction_score = 0.

    def deconstruct(self):
        del self.protein_a
        del self.position_a
        del self.protein_b
        del self.position_b
        doomsday_protocol(self)

    def set_recommended_complex(self, recommended_complex, recommended_interaction, residue_interaction_score):
        self.recommended_complex = recommended_complex
        self.recommended_interaction = recommended_interaction
        self.interaction_score = residue_interaction_score

class Aggregated_interface(Slotted_obj):
    __slots__ = [
        'protein',  'recommended_complex', 'chain',
        'interacting_chain', 'positions', 'interface_annotation_score',
        'prot_prot_interactions', 'database_id', 'interface_hash']
    def __init__(self, protein_id = None, recommended_complex = None, chain = None, interacting_chain = None, positions = None, interface_annotation_score = None, database_id = None):
        self.protein = protein_id
        self.recommended_complex = recommended_complex
        self.chain = chain
        self.interacting_chain = interacting_chain
        self.interface_annotation_score = interface_annotation_score
        if positions is None:
            self.positions = {}
        else:
            self.positions = positions
        self.prot_prot_interactions: dict[int, dict[int, dict[int, Position_Position_Interaction]]] = {}
        self.database_id = database_id
        self.interface_hash = None

    def deconstruct(self):
        del self.positions
        
        doomsday_protocol(self)

    def add_position(self, position, recommended_residue):
        self.positions[position] = recommended_residue

    def add_pos_pos_interaction(self, pos_pos_interaction: Position_Position_Interaction, score: float):
        if pos_pos_interaction.protein_b not in self.prot_prot_interactions:
            self.prot_prot_interactions[pos_pos_interaction.protein_b] = {}

        if pos_pos_interaction.position_a not in self.prot_prot_interactions[pos_pos_interaction.protein_b]:
            self.prot_prot_interactions[pos_pos_interaction.protein_b][pos_pos_interaction.position_a] = {}

        if pos_pos_interaction.position_b not in self.prot_prot_interactions[pos_pos_interaction.protein_b][pos_pos_interaction.position_a]:
            self.prot_prot_interactions[pos_pos_interaction.protein_b][pos_pos_interaction.position_a][pos_pos_interaction.position_b] = pos_pos_interaction
            return True
        elif score > self.prot_prot_interactions[pos_pos_interaction.protein_b][pos_pos_interaction.position_a][pos_pos_interaction.position_b].interaction_score:
            self.prot_prot_interactions[pos_pos_interaction.protein_b][pos_pos_interaction.position_a][pos_pos_interaction.position_b] = pos_pos_interaction
            return True
        return False

    def set_recommended_complex(self, recommended_complex, chain, interacting_chain, interface_annotation_score):
        self.recommended_complex = recommended_complex
        self.chain = chain
        self.interacting_chain = interacting_chain
        self.interface_annotation_score = interface_annotation_score

    def fuse(self, other_interface):
        if other_interface.interface_annotation_score > self.interface_annotation_score:
            return other_interface.fuse(self)
        for pos in other_interface.positions:
            self.positions[pos] = other_interface.positions[pos]
        for prot in other_interface.prot_prot_interactions:
            if prot not in self.prot_prot_interactions:
                self.prot_prot_interactions[prot] = {}
            for pos in other_interface.prot_prot_interactions[prot]:
                if pos not in self.prot_prot_interactions[prot]:
                    self.prot_prot_interactions[prot][pos] = {}
                for i_pos in other_interface.prot_prot_interactions[prot][pos]:
                    
                    self.prot_prot_interactions[prot][pos][i_pos] = other_interface.prot_prot_interactions[prot][pos][i_pos]
        return self

    def get_interface_hash(self) -> int:
        if self.interface_hash is not None:
            return self.interface_hash
        pos_list = list(self.positions.keys())
        if len(pos_list) == 0:
            return 0
        l = len(pos_list)
        mi = min(pos_list)
        ma = max(pos_list)
        self.interface_hash = 1_000_000_000*mi + 10_000*ma + l
        #print(f'{self.interface_hash:_}')
        return self.interface_hash

    def print_interface(self):
        max_prints = 100
        print('==============================================')
        print(f'Aggregated interface object: {self.protein}')
        print(f'Recommended complex: {self.recommended_complex} {self.chain} {self.interacting_chain}')
        if len(self.positions) < max_prints:
            print(f'Positions:\n{self.positions}')
        else:
            print(f'Positions, in total {len(self.positions)}')
        print(f'Interface score: {self.interface_annotation_score}')

            
        print('==============================================')

class Interface(Slotted_obj):
    __slots__ = ['chain', 'interacting_chain', 'residues', 'support_residues', 'interactions', 'stored']

    def __init__(self, chain = None, interacting_chain = None, stored = False):
        self.chain = chain
        self.interacting_chain = interacting_chain
        self.interactions = {}
        self.residues = set()
        self.support_residues = set()
        self.stored = stored

    def deconstruct(self):
        del self.residues
        del self.interactions
        del self.support_residues
        doomsday_protocol(self)

    def add_interaction(self, res_a, res_b, score):
        if res_a not in self.interactions:
            self.interactions[res_a] = {}
        self.interactions[res_a][res_b] = score
        self.residues.add(res_a)

    def add_support(self, res):
        self.support_residues.add(res)
        self.residues.add(res)

def calculate_aggregated_interface_map(
        config,
        protein_db_id: int,
        annotations: dict[str, list[str]],
        interface_map: dict[str, dict[str, dict[str, Interface]]],
        quality_measures: dict[str, dict[str, tuple[float, float]]],
        structure_to_db_ids: dict[str, dict[int, tuple[str, str]]],
        backmap_store: dict,
        db_lock = None) -> list[Aggregated_interface]:
    position_interface_map = {}
    aggregated_interfaces: list[Aggregated_interface] = []
    if config.verbosity >= 4:
        print(f'Call of calculate_aggregated_interface_map: {protein_db_id} {len(annotations)=} {len(interface_map)=}')

    n_interfaces_0 = 0
    n_interfaces_1 = 0
    n_interfaces_2 = 0
    count_skip_0 = 0
    count_skip_1 = 0
    count_skip_2 = 0
    for structure_id in annotations:
        #print(structure_to_db_ids[structure_id])
        if structure_id not in structure_to_db_ids:
            continue
        structure_backmaps: dict[str, dict[int, bytes | Residue_Map[int]]]
        if structure_id in backmap_store:
            structure_backmaps = ray.get(backmap_store[structure_id])
        else:
            structure_backmaps_ = getBackmaps(
                    structure_to_db_ids[structure_id],
                    config,
                    unpacked=False,
                    db_lock = db_lock)
            if structure_id not in structure_backmaps_:
                continue
            structure_backmaps = structure_backmaps_[structure_id]

        if structure_id not in interface_map:
            count_skip_0 += 1
            continue
        if not isinstance(interface_map[structure_id], dict):
            interface_map[structure_id] = unpack(ray.get(interface_map[structure_id]))
        interfaces = interface_map[structure_id]

        n_interfaces_0 += len(interfaces)
        for t_chain in annotations[structure_id]:
            
            if t_chain not in structure_backmaps:
                count_skip_1 += 1
                if config.verbosity >= 6:
                    print(f'skip in calculate_aggregated_interface_map: {protein_db_id}, {structure_id}, {t_chain}, {list(structure_backmaps.keys())=}')
                continue
            if protein_db_id not in structure_backmaps[t_chain]:
                count_skip_2 += 1
                if config.verbosity >= 6:
                    #This is quite normal for chains without interfaces
                    print(f'skip 2 in calculate_aggregated_interface_map: {protein_db_id}, {structure_id}, {t_chain}, {list(structure_backmaps[t_chain].keys())=}')
                continue
            if isinstance(structure_backmaps[t_chain][protein_db_id], bytes):
                backmap: Residue_Map[int] = unpack(structure_backmaps[t_chain][protein_db_id])
                structure_backmaps[t_chain][protein_db_id] = backmap
            else:
                backmap: Residue_Map[int] = structure_backmaps[t_chain][protein_db_id]

            if config.verbosity >= 6:
                print(f'\ncalculate_aggregated_interface_map: {protein_db_id}, {structure_id}, {t_chain}, {interfaces.keys()=}')

            if interfaces is None:
                config.errorlog.add_error(f'Interfaces is None: {protein_db_id} {structure_id} {t_chain}')

            for chain in interfaces:
                if chain != t_chain:
                    continue
                n_interfaces_1 += 1
                n_interfaces_2 += len(interfaces[chain])

                for i_chain in interfaces[chain]:
                    if config.verbosity >= 6:
                        print(f'At the beginning of interface aggregation process loop: {structure_id} {chain} {i_chain}')


                    other_interfaces = {}
                    back_mappable_positions = 0
                    for res in interfaces[chain][i_chain].residues:

                        position = backmap.get_item(res)
                        if position is None:
                            continue
                        back_mappable_positions += 1
                        if position in position_interface_map:
                            known_interfaces = position_interface_map[position]
                            for known_interface in known_interfaces:
                                if not known_interface in other_interfaces:
                                    other_interfaces[known_interface] = 1
                                else:
                                    other_interfaces[known_interface] += 1
                        else:
                            position_interface_map[position] = set()
                            #potential_new_interface = len(aggregated_interfaces)

                    if back_mappable_positions == 0:
                        continue

                    no_fusion_happened = True
                    fused_with = []
                    for known_interface_number in other_interfaces:
                        overlap = other_interfaces[known_interface_number]

                        if config.verbosity >= 7:
                            print(known_interface_number, chain, i_chain, overlap, back_mappable_positions, len(aggregated_interfaces[known_interface_number].positions), f'old recommended_complex: {aggregated_interfaces[known_interface_number].recommended_complex} {aggregated_interfaces[known_interface_number].chain} {aggregated_interfaces[known_interface_number].interacting_chain}')

                        #If the currently checked interface overlaps more than half with a known interface, fuse them
                        if overlap/back_mappable_positions > 0.5 or overlap/len(aggregated_interfaces[known_interface_number].positions) > 0.5:

                            if config.verbosity >= 7:
                                print(f'Fusion about to happen: {protein_db_id} {known_interface_number} {structure_id} {chain} {i_chain}')
                                print(f'Current state of the interface {known_interface_number}:')
                                aggregated_interfaces[known_interface_number].print_interface()

                            no_fusion_happened = False
                            seq_id, cov = quality_measures[structure_id][chain]
                            interface_annotation_score = calc_interface_annotation_score(seq_id, cov, interfaces[chain][i_chain].interactions)

                            if config.verbosity >= 7:
                                print(f'Interface annotation score: {interface_annotation_score}')

                            for res in interfaces[chain][i_chain].residues:
                                position = backmap.get_item(res)
                                if position is None:
                                    #print('position is none',protein_db_id, structure_id, chain, i_chain)
                                    continue
                                position_interface_map[position].add(known_interface_number)
                                if (not position in aggregated_interfaces[known_interface_number].positions) or (interface_annotation_score > aggregated_interfaces[known_interface_number].interface_annotation_score):
                                    aggregated_interfaces[known_interface_number].add_position(position, (structure_id, chain, res))
                                if res in interfaces[chain][i_chain].interactions:
                                    for i_res in interfaces[chain][i_chain].interactions[res]:
                                        residue_interaction_score = interfaces[chain][i_chain].interactions[res][i_res]

                                        if not i_chain in structure_backmaps:
                                            continue
                                        if config.verbosity >= 7:
                                            print(f'Fusing interaction: {chain} {res} <-> {i_chain} {i_res}: {len(structure_backmaps[i_chain])}')

                                        for i_protein in structure_backmaps[i_chain]:
                                            if isinstance(structure_backmaps[i_chain][i_protein], bytes):
                                                i_backmap = unpack(structure_backmaps[i_chain][i_protein])
                                                structure_backmaps[i_chain][i_protein] = i_backmap
                                            else:
                                                i_backmap = structure_backmaps[i_chain][i_protein]

                                            i_res_mapped_pos = i_backmap.get_item(i_res)

                                            if i_res_mapped_pos is None: #Unmapped residues
                                                continue

                                            pos_pos_i = Position_Position_Interaction(protein_db_id, position, i_protein, i_res_mapped_pos)
                                            new_pos_pos_i = aggregated_interfaces[known_interface_number].add_pos_pos_interaction(pos_pos_i, residue_interaction_score)

                                            if config.verbosity >= 7:
                                                print(new_pos_pos_i or (interface_annotation_score > aggregated_interfaces[known_interface_number].interface_annotation_score))

                                            if new_pos_pos_i or (interface_annotation_score > aggregated_interfaces[known_interface_number].interface_annotation_score):
                                                aggregated_interfaces[known_interface_number].prot_prot_interactions[i_protein][position][i_res_mapped_pos].set_recommended_complex(structure_id, (chain, res, i_chain, i_res), residue_interaction_score)
                                                if config.verbosity >= 7:
                                                    print(f'Overwrite of recommended complex of {position} -> {i_protein} {i_res_mapped_pos}: {new_pos_pos_i} {structure_id} {chain} {i_chain}')


                            fused_with.append(known_interface_number)
                            if interface_annotation_score > aggregated_interfaces[known_interface_number].interface_annotation_score:
                                aggregated_interfaces[known_interface_number].set_recommended_complex(structure_id, chain, i_chain, interface_annotation_score)
                                if config.verbosity >= 6:
                                    print(f'Overwrite of recommended_complex: {structure_id}, {chain}, {i_chain}')

                            if config.verbosity >= 6:
                                print(f'Fusion happend: {protein_db_id} {known_interface_number} {structure_id} {chain} {i_chain}, old recommended_complex: {aggregated_interfaces[known_interface_number].recommended_complex} {aggregated_interfaces[known_interface_number].chain} {aggregated_interfaces[known_interface_number].interacting_chain}')
                            if config.verbosity >= 7:
                                print(f'Current state of the fused interface {known_interface_number}:')
                                aggregated_interfaces[known_interface_number].print_interface()

                    if len(fused_with) > 1:
                        fused_with = sorted(fused_with)
                        if config.verbosity >= 6:
                            print(f'Superfusion of {protein_db_id} {fused_with}')
                        interface_number_1 = fused_with[0]
                        to_be_removed = set()
                        for interface_number_2 in fused_with[1:]:
                            aggregated_interfaces[interface_number_1] = aggregated_interfaces[interface_number_1].fuse(aggregated_interfaces[interface_number_2])
                            to_be_removed.add(interface_number_2)

                        new_aggregated_interfaces = []
                        new_position_interface_map = {}

                        diff_count = 0
                        interface_number_map = {}
                        for interface_number, agg_int in enumerate(aggregated_interfaces):
                            if not interface_number in to_be_removed:
                                new_aggregated_interfaces.append(agg_int)
                                interface_number_map[interface_number] = interface_number - diff_count
                            else:
                                diff_count += 1
                                interface_number_map[interface_number] = interface_number - diff_count

                        for position in position_interface_map:
                            new_position_interface_map[position] = set()
                            for known_interface in position_interface_map[position]:
                                new_position_interface_map[position].add(interface_number_map[known_interface])

                        aggregated_interfaces = new_aggregated_interfaces
                        position_interface_map = new_position_interface_map

                    if no_fusion_happened:
                        trigger = True

                        skip_0 = 0
                        skip_1 = 0
                        skip_2 = 0
                        skip_3 = 0
                        for res in interfaces[chain][i_chain].residues:

                            position = backmap.get_item(res)
                            if position is None:
                                #print('2 position is none',protein_db_id, structure_id, chain, i_chain)
                                continue

                            if trigger: #We need at least one position being backmapped
                                trigger = False
                                seq_id, cov = quality_measures[structure_id][chain]
                                interface_annotation_score = calc_interface_annotation_score(seq_id, cov, interfaces[chain][i_chain].interactions)
                                aggregated_interfaces.append(Aggregated_interface(protein_db_id, chain = chain, interacting_chain = i_chain, recommended_complex = structure_id, interface_annotation_score = interface_annotation_score))
                                if config.verbosity >= 6:
                                    print(f'Creation of new Aggregated_interface: {protein_db_id} {len(aggregated_interfaces) - 1} {chain} {i_chain} {structure_id}')

                            aggregated_interfaces[-1].add_position(position, (structure_id, chain, res))
                            position_interface_map[position].add(len(aggregated_interfaces) - 1)

                            if res not in interfaces[chain][i_chain].interactions:
                                skip_0 += 1
                                continue

                            for i_res in interfaces[chain][i_chain].interactions[res]:
                                residue_interaction_score = interfaces[chain][i_chain].interactions[res][i_res]
                                
                                if not i_chain in structure_backmaps:
                                    if config.verbosity >= 7:
                                        print(f'{i_chain} not in structure_back_maps[{structure_id}], {protein_db_id}, {chain}, {res}, {i_res}')
                                    skip_2 += 1
                                    continue
                                if config.verbosity >= 6:
                                    print(f'In residue aggregation of {res} - {i_res}')

                                for i_protein in structure_backmaps[i_chain]:
                                    if isinstance(structure_backmaps[i_chain][i_protein], bytes):
                                        i_backmap = unpack(structure_backmaps[i_chain][i_protein])
                                        structure_backmaps[i_chain][i_protein] = i_backmap
                                    else:
                                        i_backmap = structure_backmaps[i_chain][i_protein]
                                    i_protein_pos = i_backmap.get_item(i_res)
                                    if i_protein_pos is None: #Unmapped residues
                                        if config.verbosity >= 7:
                                            print(f'Unmapped residue: {structure_id} {i_chain} {i_protein} {i_res}')
                                        skip_3 += 1
                                        continue
                                    pos_pos_i = Position_Position_Interaction(protein_db_id, position, i_protein, i_protein_pos)
                                    new_pos_pos_i = aggregated_interfaces[-1].add_pos_pos_interaction(pos_pos_i, residue_interaction_score)
                                    aggregated_interfaces[-1].prot_prot_interactions[i_protein][position][i_protein_pos].set_recommended_complex(structure_id, (chain, res, i_chain, i_res), residue_interaction_score)
                        if config.verbosity >= 4 and len(aggregated_interfaces[-1].prot_prot_interactions) == 0:
                            print(f'Added Empty Agg_interfaces: {protein_db_id} {structure_id} {chain=} {i_chain=} {skip_0=} {skip_1=} {skip_2=} {skip_3=}')

    if config.verbosity >= 4:
        print(f'End of calculate_aggregated_interface_map: {protein_db_id} {n_interfaces_0=} {n_interfaces_1=} {n_interfaces_2=} {count_skip_0} {count_skip_1=} {count_skip_2=} {len(aggregated_interfaces)=}')

    return aggregated_interfaces




