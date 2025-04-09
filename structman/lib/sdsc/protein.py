import sys
from structman.lib.sdsc.consts import codons
from structman.lib.sdsc.indel import Indel
from structman.lib.sdsc.mutations import MultiMutation
from structman.lib.sdsc.sdsc_utils import process_recommend_structure_str, doomsday_protocol, Slotted_obj
from structman.lib.sdsc.position import Position


class Protein(Slotted_obj):
    __slots__ = [
        'primary_protein_id',       'u_ac',                         'u_id',
        'ref_id',                   'ref_nt_id',                    'other_ids',
        'pdb_id',                   'positions',                    'res_id_map',
        'sequence',                 'nucleotide_sequence',          'stored',
        'completely_stored',        'wildtype_protein',             'gene',
        'tags',                     'go_terms',                     'pathways',
        'database_id',              'structure_annotations',        'mapped_complexes',    
        'input_id',                 'multi_mutations',              'aggregated_contact_matrix',
        'aggregated_interface_map', 'mutant_type',                  'sav_positions',
        'insertion_positions',      'deletion_flanks',              'number_of_mappings'
        ]

    slot_mask = [
        True, True, True,
        True, True, True,
        True, True, True,
        False, True, True,
        True, True, True,
        True, True, True,
        True, True, True,
        True, True, True,
        True, True, True,
        True, True, True      
    ]

    def __init__(self, errorlog = None, primary_protein_id=None, u_ac=None, u_id=None, ref_id=None, ref_nt_id = None, wildtype_protein=None, gene = None,
                 pdb_id=None, positions=[], database_id=None, other_ids=[], input_id=None, sequence=None, aggregated_contact_matrix = {},
                 aggregated_interface_map = {}, mutant_type = None, sav_positions = None, insertion_positions = None, deletion_flanks = None,
                 tags = None
                ):
        self.primary_protein_id = primary_protein_id
        self.u_ac = u_ac  # UNIPROT accession number
        self.u_id = u_id  # UNIPROT ID
        self.ref_id = ref_id  # RefSeq ID
        self.ref_nt_id = ref_nt_id
        self.pdb_id = pdb_id
        self.input_id = input_id
        self.positions = [None]
        self.res_id_map = {}
        self.sequence = sequence
        self.nucleotide_sequence = None
        self.stored = False  # If the protein is stored in the database
        self.completely_stored = False  # Gets true if all corresponding input positions are stored
        self.database_id = database_id  # ID in the corresponding database table
        self.stored = (database_id is not None)
        self.go_terms = {}
        self.pathways = {}
        self.structure_annotations = {}
        self.mapped_complexes = set()
        self.multi_mutations = []
        self.wildtype_protein = wildtype_protein
        self.aggregated_contact_matrix = aggregated_contact_matrix
        self.aggregated_interface_map = aggregated_interface_map
        self.mutant_type = mutant_type # Can be None, 'SAV', 'indel', or 'multi'
        self.sav_positions = sav_positions
        self.insertion_positions = insertion_positions
        self.deletion_flanks = deletion_flanks
        self.number_of_mappings = 0

        self.gene = gene
        if tags is None:
            self.tags = set([])
        else:
            self.tags = tags
        if self.wildtype_protein is None and (self.mutant_type == 'indel' or self.mutant_type == 'multi'):
            self.wildtype_protein = self.primary_protein_id.split('_del')[0].split('_ins')[0]
        else:
            self.wildtype_protein = self.primary_protein_id

        if pdb_id is None:
            warns = self.add_positions(positions)
        else:
            warns = self.add_residues(positions)
        if warns is not None:
            for warn in warns:
                if errorlog is not None:
                    errorlog.add_warning(warn)
        self.other_ids = {}
        for id_id, other_id in other_ids:
            self.other_ids[id_id] = other_id

    def deconstruct(self):
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            pos_obj.deconstruct(completely = True)

        del self.positions
        del self.sequence
        for struct_id in self.structure_annotations:
            for chain in self.structure_annotations[struct_id]:
                self.structure_annotations[struct_id][chain].deconstruct()
        del self.structure_annotations
        for interface_obj in self.aggregated_interface_map:
            interface_obj.deconstruct()
        del self.aggregated_interface_map        
        doomsday_protocol(self)

    def print_state(self):
        print(f'----State of {self.primary_protein_id} ({len(self.positions)})----')
        print('Uniprot Id:', self.u_id)
        for pos_obj in self.positions:
            if pos_obj is None:
                print('None pos_obj')
                continue
            pos_obj.print_state()

    def print_all_snvs(self):
        print(f'---All SNVs of {self.primary_protein_id}:')
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            if len(pos_obj.mut_aas) > 0:
                pos_obj.print_state()
                

    def get_u_ac(self):
        if self.u_ac is not None:
            return self.u_ac
        else:
            return self.pdb_id

    def add_positions(self, positions):
        warns = []
        del_list = []
        for p, position in enumerate(positions):
            warn = False
            if position.pos >= len(self.positions):
                d = position.pos - len(self.positions)
                for _ in range(d):
                    self.positions.append(None)
                self.positions.append(position)
            elif self.positions[position.pos] is None:
                self.positions[position.pos] = position
            elif position.pos is not None:
                warn = self.positions[position.pos].fuse(position)
                del_list.append(p)
                if warn:
                    warns.append('Warning happened in add_positions fuse: %s %s' % (self.u_ac, self.pdb_id))

        del_list = sorted(del_list, reverse=True)

        for p in del_list:
            del positions[p]

        if len(warns) == 0:
            return None
        else:
            return warns

    def add_multi_mutation(self, multi_mutation, indel_map, mm_tags = None, mut_prot_id = None):
        corrected_multi_mutation = []
        for indel_or_snv in multi_mutation:
            if isinstance(indel_or_snv, tuple):
                pos_obj = indel_or_snv[0]
                aa2 = indel_or_snv[1]
                if isinstance(pos_obj, int):
                    correct_pos_obj = self.positions[pos_obj]
                else:
                    correct_pos_obj = self.positions[pos_obj.pos]
                corrected_snv = (correct_pos_obj, aa2)
                corrected_multi_mutation.append(corrected_snv)
            else:
                if indel_or_snv.get_notation() in indel_map:
                    indel = indel_map[indel_or_snv.get_notation()]
                else:
                    indel = indel_or_snv
                corrected_multi_mutation.append(indel)
        self.multi_mutations.append((corrected_multi_mutation, mut_prot_id, mm_tags))
        return

    def add_residues(self, positions):
        warns = []
        for position in positions:
            warn = False
            if position.pdb_res_nr not in self.res_id_map:
                self.res_id_map[position.pdb_res_nr] = position
            else:
                warn = self.res_id_map[position.pdb_res_nr].fuse(position)
                if warn:
                    warns.append('Warning happened in add_residues fuse: %s %s' % (self.u_ac, self.pdb_id))
        if len(warns) == 0:
            return None
        else:
            return warns

    def popNone(self):
        if self.pdb_id is None:
            if self.positions[0] is not None:
                tags = self.positions[0].pos_tags
                self.positions[0] = None
                return True, tags
            else:
                return True, None
        else:
            if None in self.res_id_map:
                tags = self.res_id_map[None].pos_tags
                del self.res_id_map[None]
                return True, tags
            else:
                return True, None

    def getRecommendedStructures(self):
        rec_structs = {}
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            recommended_structure = pos_obj.get_recommended_structure()
            if recommended_structure == '-':
                continue
            struct_tuple, res_info = recommended_structure.split()
            if struct_tuple not in rec_structs:
                rec_structs[struct_tuple] = [pos_obj.pos]
            else:
                rec_structs[struct_tuple].append(pos_obj.pos)
        return rec_structs

    def get_recommended_structure(self, scheme = 'functional_annotation', prefer_mutation_positions = False, custom_mutations = None, functionally_weighted = False):
        if scheme == 'functional_annotation':
            return self.get_major_recommend_structure(prefer_mutation_positions = prefer_mutation_positions, custom_mutations = custom_mutations, functionally_weighted = functionally_weighted)
        if scheme == 'most_similar':
            return self.get_most_similar_structure()
        return f'Unknown scheme: {scheme}'

    def get_most_similar_structure(self):
        max_score = 0
        rec_struct = None
        for structure_id in self.structure_annotations:
            for chain in self.structure_annotations[structure_id]:
                score = self.structure_annotations[structure_id][chain].sequence_identity * self.structure_annotations[structure_id][chain].coverage
                if score > max_score:
                    max_score = score
                    rec_struct = f'{structure_id}:{chain}'
        return rec_struct

    def get_major_recommend_structure(self, prefer_mutation_positions = False, custom_mutations = None, functionally_weighted = False):
        rec_structs = {}

        if custom_mutations is None:
            savs = self.sav_positions
            ins_pos = self.insertion_positions
            del_fl = self.deletion_flanks
        else:
            savs, ins_pos, del_fl = custom_mutations

        mut_pos = set([])
        if prefer_mutation_positions:
            prefer_mutation_positions = False
            if savs is not None:
                prefer_mutation_positions = True
                for pos in savs:
                    mut_pos.add(pos)
            if ins_pos is not None:
                prefer_mutation_positions = True
                for insert in ins_pos:
                    for pos in insert:
                        mut_pos.add(pos)
            if del_fl is not None:
                prefer_mutation_positions = True
                for lf, rf in del_fl:
                    if lf is not None:
                        mut_pos.add(lf)
                    if rf is not None:
                        mut_pos.add(rf)
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            if prefer_mutation_positions and pos_obj.pos not in mut_pos:
                continue
            recommended_structure = pos_obj.get_recommended_structure()
            if recommended_structure == '-':
                continue
            struct_tuple, res_info = recommended_structure.split()
            if struct_tuple not in rec_structs:
                rec_structs[struct_tuple] = 0
            weight = 1
            if functionally_weighted:
                weight = pos_obj.get_score()

            rec_structs[struct_tuple] += weight
            rec_structs[struct_tuple] += 1

        max_count = 0
        maj_rec_struct = None
        for struct_tuple in rec_structs:
            if rec_structs[struct_tuple] > max_count:
                max_count = rec_structs[struct_tuple]
                maj_rec_struct = struct_tuple
        return maj_rec_struct

    def get_mutation_position_inverts(self):
        if self.deletion_flanks is not None:
            wt_inserts = []
            for lf, rf in self.deletion_flanks:
                insert = []
                for i in range(lf, (rf+1)):
                    insert.append(i)
                wt_inserts.append(insert)
        else:
            wt_inserts = None

        if self.insertion_positions is not None:
            wt_del_flanks = []
            for insert in self.insertion_positions:
                lf = min(insert)
                rf = max(insert)
                wt_del_flanks.append((lf, rf))
        else:
            wt_del_flanks = None

        return wt_inserts, wt_del_flanks

    def invert_mutation_positions(self, mutant_protein):
        self.sav_positions = mutant_protein.sav_positions

        wt_inserts, wt_del_flanks = mutant_protein.get_mutation_position_inverts()

        self.insertion_positions = wt_inserts
        self.deletion_flanks = wt_del_flanks
        return

    def getAACList(self):
        aaclist = {}
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            aac_base = pos_obj.getAACBase()
            aaclist[aac_base] = pos_obj.get_database_id()
        return aaclist

    def get_aac_base(self, pos):
        return self.positions[pos].getAACBase()

    def status(self):
        stat = '\n'.join([self.u_ac, self.pdb_id, str(len(self.positions)), self.sequence, self.stored, self.database_id, str(len(self.structure_annotations))])
        return stat

    def is_sequence_set(self):
        return self.sequence is not None

    def set_sequence(self, value):
        self.sequence = value

    def get_sequence(self):
        return self.sequence

    def pop_sequence(self):
        seq = self.sequence
        self.sequence = None
        return seq

    def contains_position(self, pos):
        if pos < len(self.positions):
            if self.positions[pos] is None:
                return False
            else:
                return True
        else:
            return False

    def set_stored(self, value):
        self.stored = value

    def set_completely_stored(self):
        self.completely_stored = True
        return self.database_id

    def is_stored(self):
        return self.stored

    def add_other_ids(self, id_id, other_id):
        self.other_ids[id_id] = other_id

    def get_u_ac(self):
        return self.u_ac

    def get_ref_id(self):
        return self.ref_id

    def get_u_id(self):
        return self.u_id

    def get_go_terms(self):
        return self.go_terms

    def get_pathways(self):
        return self.pathways

    def is_position_stored(self, pos):
        if not self.contains_position(pos):
            return False
        else:
            return self.positions[pos].get_stored()

    def get_position(self, pos):
        if pos >= len(self.positions):
            return None
        return self.positions[pos]

    def get_position_ids(self):
        return list(range(1,len(self.positions)))

    def set_position_stored(self, pos, stored_value):
        self.positions[pos].set_stored(stored_value)

    def set_position_database_id(self, pos, value):
        self.positions[pos].set_database_id(value)

    def get_position_database_id(self, pos):
        return self.positions[pos].get_database_id()

    def get_pos_tags(self, pos):
        return self.positions[pos].get_pos_tags()

    def set_database_id(self, value):
        self.database_id = value

    def get_database_id(self):
        return self.database_id

    def get_res_id(self, pos):
        return self.positions[pos].get_res_id()

    def get_wt_aa(self, pos):
        return self.positions[pos].get_wt_aa()

    def get_mut_aas(self, pos):
        return self.positions[pos].get_mut_aas()

    def get_mut_tags(self, pos, new_aa):
        if pos >= len(self.positions):
            return None
        if self.positions[pos] is None:
            return None
        return self.positions[pos].get_mut_tags(new_aa)

    def retrieve_tag_value(self, tag_key, float_type = True, delimiter = '='):
        tag_value = None
        for tag in self.tags:
            split_result = tag.split(delimiter)
            if len(split_result) != 2:
                continue
            key, value = split_result
            if key != tag_key:
                continue
            elif float_type:
                try:
                    tag_value = float(value)
                except:
                    tag_value = None
        return tag_value

    def add_annotation(self, pdb_id, chain, anno_obj):
        self.mapped_complexes.add(pdb_id)
        if pdb_id not in self.structure_annotations:
            self.structure_annotations[pdb_id] = {}
        self.structure_annotations[pdb_id][chain] = anno_obj

    def remove_annotation(self, pdb_id, chain):
        if not pdb_id in self.structure_annotations:
            return
        elif not chain in self.structure_annotations[pdb_id]:
            return
        else:
            del self.structure_annotations[pdb_id][chain]
            if len(self.structure_annotations[pdb_id]) == 0:
                del self.structure_annotations[pdb_id]

    def remove_annotations(self):
        del self.structure_annotations

    def get_annotation_list(self) -> list[tuple[str, str]]:
        annotation_list: list[tuple[str, str]] = []
        for structure_id in self.structure_annotations:
            for chain in self.structure_annotations[structure_id]:
                annotation_list.append((structure_id, chain))
        return annotation_list

    def get_structure_annotations(self):
        return self.structure_annotations

    def is_annotation_stored(self, pdb_id, chain):
        return self.structure_annotations[pdb_id][chain].get_stored()

    def set_alignment(self, pdb_id, chain, value):
        self.structure_annotations[pdb_id][chain].set_alignment(value)

    def get_alignment(self, pdb_id, chain):
        if not pdb_id in self.structure_annotations:
            return 'Structure %s:%s not in annotation list of %s' % (pdb_id, chain, self.primary_protein_id)
        if not chain in self.structure_annotations[pdb_id]:
            return 'Structure %s:%s not in annotation list of %s' % (pdb_id, chain, self.primary_protein_id)

        return self.structure_annotations[pdb_id][chain].get_alignment()

    def pop_alignment(self, pdb_id, chain):
        return self.structure_annotations[pdb_id][chain].pop_alignment()

    def set_coverage(self, pdb_id, chain, value):
        self.structure_annotations[pdb_id][chain].set_coverage(value)
        self.number_of_mappings += int(len(self.positions) * value)

    def get_coverage(self, pdb_id, chain):
        return self.structure_annotations[pdb_id][chain].get_coverage()

    def set_sequence_id(self, pdb_id, chain, value):
        self.structure_annotations[pdb_id][chain].set_sequence_id(value)

    def get_sequence_id(self, pdb_id, chain):
        return self.structure_annotations[pdb_id][chain].get_sequence_id()

    def set_annotation_db_id(self, pdb_id, chain, value):
        self.structure_annotations[pdb_id][chain].set_database_id(value)

    def set_sub_infos(self, pdb_id, chain, value):
        self.structure_annotations[pdb_id][chain].set_sub_infos(value)

    def get_sub_infos(self, pdb_id, chain):
        return self.structure_annotations[pdb_id][chain].get_sub_infos()

    def set_backmap(self, structure_id, chain, backmap):
        self.structure_annotations[structure_id][chain].set_backmap(backmap)

    def get_backmap(self, structure_id, chain):
        if structure_id not in self.structure_annotations:
            return {}
        if chain in self.structure_annotations[structure_id]:
            return {}
        return self.structure_annotations[structure_id][chain].get_backmap()

    def set_aggregated_contact_matrix(self, aggregated_contact_matrix):
        self.aggregated_contact_matrix = aggregated_contact_matrix

    def set_aggregated_interface_map(self, aggregated_interface_map):
        self.aggregated_interface_map = aggregated_interface_map

    def classify(self, pos, config):
        self.positions[pos].classify(config)

    def get_ordered_classifications(self):
        classifications = []
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            c = pos_obj.get_classification()
            classifications.append(c)
        return classifications

    def mutate_snvs(self, proteins, config):
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            mut_seq = self.sequence
            for aa2 in pos_obj.mut_aas:
                protein_name = '%s_%s%s%s' % (self.primary_protein_id, pos_obj.wt_aa, str(pos_obj.pos), aa2)

                mut_seq = '%s%s%s' % (mut_seq[:(pos_obj.pos - 1)], aa2, mut_seq[pos_obj.pos:])
                mut_protein = Protein(config.errorlog, primary_protein_id=protein_name, sequence=mut_seq, wildtype_protein=self.primary_protein_id, mutant_type = 'SAV', sav_positions = [pos_obj.pos])
                proteins[protein_name] = mut_protein

                if config.verbosity >= 4:
                    print('Created new mutant:', protein_name)

                for (mpos, aa) in enumerate(mut_seq):
                    seq_pos = mpos + 1
                    position = Position(pos=seq_pos, wt_aa=aa, checked=True, recommended_structure = pos_obj.recommended_structure)
                    proteins[protein_name].positions[seq_pos] = position

    def create_multi_mutations(self, proteins, config):
        indel_position_order = {}
        multi_mutation_objects = []
        if config.verbosity >= 6:
            print('Starting multi mutation mutation of', self.primary_protein_id)
        for mm_nr, multi_mutation_tuple in enumerate(self.multi_mutations):
            multi_mutation, mut_prot_id, mm_tags = multi_mutation_tuple

            if config.verbosity >= 6:
                print(f'In create_multi_mutations, looping over new MM: {mut_prot_id}')
            name_parts = []

            mutant_positions = [[], [], []]
            for mut in multi_mutation:
                if not isinstance(mut, tuple): #In case of indels
                    name_parts.append(mut.get_notation())
                    insertion_pos, deletion_flanks = mut.get_mutant_positions()
                    if insertion_pos is not None:
                        mutant_positions[1].append(insertion_pos)
                    if deletion_flanks is not None:
                        mutant_positions[2].append(deletion_flanks)
                else:
                    position, aa2 = mut
                    mutant_positions[0].append(position.pos)
                    name_parts.append('%s%s%s' % (position.wt_aa, str(position.pos), aa2))
            
            if mut_prot_id is None:
                protein_name = f'{self.primary_protein_id}_{"_".join(name_parts)}'[:255]

                #print(f'Before shifts: {mutant_positions}')

                correction_map = []
                for insertion in mutant_positions[1]:
                    m_pos = max(insertion)
                    s = len(insertion)
                    correction_map.append((m_pos, s))
                for lf, rf in mutant_positions[2]:
                    if rf is not None:
                        if lf is None:
                            s = rf
                        else:
                            s = 1 + (rf - lf)
                        correction_map.append((rf, -s))

                for sav_pos, sav in enumerate(mutant_positions[0]):
                    total_shift = 0
                    for cor_pos, shift in correction_map:
                        if cor_pos < sav:
                            total_shift += shift
                    mutant_positions[0][sav_pos] = sav + total_shift

                for insertion_pos, insertion in enumerate(mutant_positions[1]):
                    for insertion_pos_pos, ins in enumerate(insertion):
                        total_shift = 0
                        for cor_pos, shift in correction_map:
                            if cor_pos < ins:
                                total_shift += shift
                        mutant_positions[1][insertion_pos][insertion_pos_pos] = ins + total_shift

                for deletion_pos, deletion_flanks in enumerate(mutant_positions[2]):
                    lf,rf = deletion_flanks
                    if lf is not None:
                        total_shift = 0
                        for cor_pos, shift in correction_map:
                            if cor_pos < lf:
                                total_shift += shift
                        shifted_lf = lf + total_shift
                    else:
                        shifted_lf = None
                    if rf is not None:
                        total_shift = 0
                        for cor_pos, shift in correction_map:
                            if cor_pos < rf:
                                total_shift += shift
                        shifted_rf = rf + total_shift
                    else:
                        shifted_rf = None

                    mutant_positions[2][deletion_pos] = (shifted_lf, shifted_rf)

                if len(mutant_positions[0]) == 0:
                    mutant_positions[0] = None

                if len(mutant_positions[1]) == 0:
                    mutant_positions[1] = None

                if len(mutant_positions[2]) == 0:
                    mutant_positions[2] = None

                #print(f'After shifts: {mutant_positions}')

                if mutant_positions[1] != None or mutant_positions[2] != None:
                    mut_protein = Protein(config.errorlog, primary_protein_id=protein_name, wildtype_protein=self.primary_protein_id, mutant_type = 'Multi', sav_positions = mutant_positions[0], insertion_positions = mutant_positions[1], deletion_flanks = mutant_positions[2])
                    proteins[protein_name] = mut_protein
                    if config.verbosity >= 4:
                        print('Created new mutant:', protein_name)
                    
            else:
                try:
                    if mut_prot_id not in proteins:
                        continue
                except:
                    if not proteins.contains(mut_prot_id):
                        continue
                protein_name = mut_prot_id

            multi_mutation_obj = MultiMutation(self.primary_protein_id, protein_name, multi_mutation, tags = mm_tags)
            multi_mutation_objects.append(multi_mutation_obj)

        return multi_mutation_objects
    
    def get_mutation_positions_by_mm(self):
        mutation_position_map = {}
        for multi_mutation_tuple in self.multi_mutations:
            multi_mutation, mut_prot_id, mm_tags = multi_mutation_tuple
            mutation_position_map[mut_prot_id] = [[], [], []]
            for mut in multi_mutation:
                if not isinstance(mut, tuple): #In case of indels
                    insertion_pos, deletion_flanks = mut.get_mutant_positions()
                    if insertion_pos is not None:
                        mutation_position_map[mut_prot_id][1].append(insertion_pos)
                    if deletion_flanks is not None:
                        mutation_position_map[mut_prot_id][2].append(deletion_flanks)
                else:
                    position, aa2 = mut
                    mutation_position_map[mut_prot_id][0].append(position.pos)

        return mutation_position_map


class Proteins(Slotted_obj):
    __slots__ = ['protein_map', 'stored_ids', 'stored_ids_mutant_excluded', 'completely_stored_ids',
                 'not_stored_ids', 'id_map', 'structures', 'complexes', 'indels', 'multi_mutations', 'multi_mutation_back_map']

    def __init__(self, proteins, indels, multi_mutation_objects):
        self.protein_map = proteins
        self.stored_ids = set()
        self.stored_ids_mutant_excluded = set()
        self.completely_stored_ids = set()
        self.not_stored_ids = set()
        self.id_map = {}
        self.structures = {}
        self.complexes = {}
        self.indels = {}
        self.multi_mutations = {}
        for multi_mutation_obj in multi_mutation_objects:
            if multi_mutation_obj.wt_prot not in self.multi_mutations:
                self.multi_mutations[multi_mutation_obj.wt_prot] = []
            self.multi_mutations[multi_mutation_obj.wt_prot].append(multi_mutation_obj)
        self.multi_mutation_back_map = {}
        for prot_id in indels:
            self.indels[prot_id] = {}
            for indel_notation in indels[prot_id]:
                self.indels[prot_id][indel_notation] = indels[prot_id][indel_notation]

    def __getitem__(self, key):
        return self.protein_map[key]

    def __setitem__(self, key, prot):
        self.protein_map[key] = prot
        if prot.stored:
            self.stored_ids.add(prot.get_database_id())
            self.id_map[prot.get_database_id()] = key

    def __delitem__(self, key):
        del self.protein_map[key]

    def semi_deconstruct(self):
        del self.stored_ids
        del self.stored_ids_mutant_excluded
        del self.not_stored_ids
        del self.id_map
        for struct in self.structures:
            for chain in self.structures[struct]:
                self.structures[struct][chain].deconstruct()
        del self.structures
        for comp in self.complexes:
            self.complexes[comp].deconstruct()
        del self.complexes

    def deconstruct(self):
        try:
            self.semi_deconstruct()
        except:
            pass
        for prot_id in self.protein_map:
            self.protein_map[prot_id].deconstruct()
        del self.protein_map
        for prot_id in self.indels:
            for indel_id in self.indels[prot_id]:
                self.indels[prot_id][indel_id].deconstruct()
        del self.indels
        for mm in self.multi_mutations:
            try:
                self.multi_mutations[mm].deconstruct()
            except:
                pass
        del self.multi_mutations
        doomsday_protocol(self)

    def remove_protein(self, prot_id):
        self.protein_map[prot_id].deconstruct()
        del self.protein_map[prot_id]
        try:
            del self.multi_mutations[prot_id]
        except:
            pass
        try:
            del self.indels[prot_id]
        except:
            pass

    def remove_protein_annotations(self, u_ac):
        self.protein_map[u_ac].remove_annotations()

    def print_protein_state(self):
        for key in self.protein_map:
            self.protein_map[key].print_state()

    def get_protein_map(self):
        return self.protein_map

    def get_protein_ids(self):
        return self.protein_map.keys()

    def get_position(self, u_ac, pos):
        return self.protein_map[u_ac].get_position(pos)

    def get_complexes(self):
        return self.complexes

    def get_complex_list(self):
        return self.complexes.keys()

    def add_complex(self, pdb_id, complex_obj):
        if pdb_id in self.complexes:
            complex_obj.set_atom_count(self.complexes[pdb_id].get_atom_count())
        self.complexes[pdb_id] = complex_obj

    def contains_complex(self, pdb_id):
        return pdb_id in self.complexes

    def is_complex_stored(self, pdb_id):
        return self.complexes[pdb_id].get_stored()

    def set_complex_db_id(self, pdb_id, value):
        self.complexes[pdb_id].set_database_id(value)

    def get_complex_chains(self, pdb_id, only_protein=False):
        return self.complexes[pdb_id].get_chains(only_protein=only_protein)

    def set_IAmap(self, structure_id, IAmap):
        self.complexes[structure_id].set_IAmap(IAmap)

    def get_IAmap(self, structure_id):
        return self.complexes[structure_id].get_IAmap()

    def set_interfaces(self, structure_id, interfaces):
        self.complexes[structure_id].set_interfaces(interfaces)

    def get_interfaces(self, structure_id):
        return self.complexes[structure_id].get_interfaces()

    def print_interfaces(self):
        for structure_id in self.complexes:
            print(structure_id)
            self.complexes[structure_id].print_interfaces()

    def set_backmap(self, protein_id, structure_id, chain, backmap):
        #self.protein_map[protein_id].set_backmap(structure_id, chain, backmap)
        self.structures[structure_id][chain].backmaps[protein_id] = backmap

    def get_backmap(self, protein_id, structure_id, chain):
        return self.protein_map[protein_id].get_backmap(structure_id, chain)

    def set_aggregated_contact_matrix(self, protein_id, aggregated_contact_matrix):
        self.protein_map[protein_id].set_aggregated_contact_matrix(aggregated_contact_matrix)

    def set_aggregated_interface_map(self, protein_id, aggregated_interface_map):
        self.protein_map[protein_id].set_aggregated_interface_map(aggregated_interface_map)


    def set_interaction_partners(self, pdb_id, value):
        self.complexes[pdb_id].set_interaction_partners(value)

    def get_interaction_partners(self, pdb_id):
        return self.complexes[pdb_id].get_interaction_partners()

    def set_chain_type_map(self, pdb_id, value, chainlist):
        self.complexes[pdb_id].set_chain_type_map(value, chainlist)

    def set_atom_count(self, pdb_id, atom_count):
        self.complexes[pdb_id].set_atom_count(atom_count)

    def get_atom_count(self, pdb_id):
        return self.complexes[pdb_id].get_atom_count()

    def get_complex_structures(self, pdb_id):
        chains = self.get_complex_chains(pdb_id)
        structures = {pdb_id : {}}
        for chain in chains:
            if not pdb_id in self.structures:
                continue
            if not chain in self.structures[pdb_id]:
                continue
            structures[pdb_id][chain] = self.structures[pdb_id][chain]
        return structures

    def get_structures(self):
        return self.structures

    def add_structure(self, pdb_id, chain, struct_obj):
        if pdb_id not in self.structures:
            self.structures[pdb_id] = {}
        self.structures[pdb_id][chain] = struct_obj

    def contains_structure(self, pdb_id, chain):
        if pdb_id in self.structures:
            if chain in self.structures[pdb_id]:
                return True

        return False

    def set_last_residue(self, pdb_id, chain, last_residue):
        self.structures[pdb_id][chain].set_last_residue(last_residue)
        return

    def set_first_residue(self, pdb_id, chain, first_residue):
        self.structures[pdb_id][chain].set_first_residue(first_residue)
        return

    def set_structure_db_id(self, pdb_id, chain, value):
        self.structures[pdb_id][chain].set_database_id(value)

    def get_structure_db_id(self, pdb_id, chain):
        return self.structures[pdb_id][chain].get_database_id()

    def set_structure_stored(self, pdb_id, chain, value):
        self.structures[pdb_id][chain].set_stored(value)

    def is_structure_stored(self, pdb_id, chain):
        return self.structures[pdb_id][chain].get_stored()

    def contains_residue(self, pdb_id, chain, res_nr):
        if not pdb_id in self.structures:
            return False
        if not chain in self.structures[pdb_id]:
            return False
        return self.structures[pdb_id][chain].contains_residue(res_nr)

    def get_protein_annotation_list(self, u_ac) -> list[tuple[str, str]]:
        return self.protein_map[u_ac].get_annotation_list()

    def get_protein_structure_annotations(self, u_ac):
        return self.protein_map[u_ac].get_structure_annotations()

    def add_mapping_to_structure(self, pdb_id, chain, u_ac):
        self.structures[pdb_id][chain].add_mapping(u_ac)

    def is_annotation_stored(self, pdb_id, chain, u_ac):
        return self.protein_map[u_ac].is_annotation_stored(pdb_id, chain)

    def set_alignment(self, u_ac, pdb_id, chain, alignment):
        self.protein_map[u_ac].set_alignment(pdb_id, chain, alignment)

    def set_alignment_by_db_id(self, prot_id, pdb_id, chain, value):
        self.getByDbId(prot_id).set_alignment(pdb_id, chain, value)

    def get_alignment(self, u_ac, pdb_id, chain):
        return self.protein_map[u_ac].get_alignment(pdb_id, chain)

    def pop_alignment(self, u_ac, pdb_id, chain):
        return self.protein_map[u_ac].pop_alignment(pdb_id, chain)

    def set_coverage_by_db_id(self, prot_id, pdb_id, chain, value):
        self.getByDbId(prot_id).set_coverage(pdb_id, chain, value)

    def set_coverage(self, u_ac, pdb_id, chain, value):
        self.protein_map[u_ac].set_coverage(pdb_id, chain, value)

    def get_coverage(self, u_ac, pdb_id, chain):
        return self.protein_map[u_ac].get_coverage(pdb_id, chain)

    def set_sequence_id_by_db_id(self, prot_id, pdb_id, chain, value):
        self.getByDbId(prot_id).set_sequence_id(pdb_id, chain, value)

    def set_sequence_id(self, u_ac, pdb_id, chain, value):
        self.protein_map[u_ac].set_sequence_id(pdb_id, chain, value)

    def get_sequence_id(self, u_ac, pdb_id, chain):
        return self.protein_map[u_ac].get_sequence_id(pdb_id, chain)

    def is_sequence_set(self, u_ac):
        return self.protein_map[u_ac].is_sequence_set()

    def set_annotation_db_id_by_db_id(self, prot_id, pdb_id, chain, value):
        self.getByDbId(prot_id).set_annotation_db_id(pdb_id, chain, value)

    def set_sub_infos(self, u_ac, pdb_id, chain, value):
        self.protein_map[u_ac].set_sub_infos(pdb_id, chain, value)

    def get_sub_infos(self, u_ac, pdb_id, chain):
        return self.protein_map[u_ac].get_sub_infos(pdb_id, chain)

    def set_not_stored_ids(self, value):
        self.not_stored_ids = value

    def get_not_stored_ids(self):
        return self.not_stored_ids

    def get_not_stored_acs(self):
        acs = []
        for db_id in self.not_stored_ids:
            acs.append(self.id_map[db_id])
        return acs

    def set_stored_ids(self, stored_ids, stored_ids_mutant_excluded):
        self.stored_ids = stored_ids
        self.stored_ids_mutant_excluded = stored_ids_mutant_excluded

    def get_stored_ids(self, exclude_indel_mutants=False, exclude_completely_stored=False):
        if not exclude_indel_mutants:
            if exclude_completely_stored:
                return self.stored_ids - self.completely_stored_ids
            return self.stored_ids
        else:
            if exclude_completely_stored:
                return self.stored_ids_mutant_excluded - self.completely_stored_ids
            return self.stored_ids_mutant_excluded

    def contains(self, prot_id):
        return prot_id in self.protein_map

    def isEmpty(self):
        return len(self.protein_map) == 0

    def is_protein_stored(self, u_ac):
        return self.protein_map[u_ac].stored

    def set_completely_stored(self, u_ac):
        prot_id = self.protein_map[u_ac].set_completely_stored()
        self.completely_stored_ids.add(prot_id)

    def is_completely_stored(self, u_ac):
        return self.protein_map[u_ac].completely_stored

    def isStored(self, prot_id):
        return prot_id in self.stored_ids

    def generate_id_map(self):
        for u_ac in self.protein_map:
            self.id_map[self.protein_map[u_ac].get_database_id()] = u_ac

    def getU_acByDbId(self, database_id):
        return self.id_map[database_id]

    def getByDbId(self, database_id):
        return self.protein_map[self.id_map[database_id]]

    def set_protein_stored(self, u_ac, value):
        self.protein_map[u_ac].set_stored(value)

    def set_protein_db_id(self, u_ac, value):
        self.protein_map[u_ac].set_database_id(value)

    def get_protein_db_id(self, u_ac):
        return self.protein_map[u_ac].get_database_id()

    def set_protein_sequence(self, u_ac, value):
        self.protein_map[u_ac].set_sequence(value)

    def get_sequence(self, u_ac):
        return self.protein_map[u_ac].get_sequence()

    def pop_sequence(self, u_ac):
        return self.protein_map[u_ac].pop_sequence()

    def get_u_ac(self, prot_id):
        return self.protein_map[prot_id].get_u_ac()

    def get_ref_id(self, prot_id):
        return self.protein_map[prot_id].get_ref_id()

    def get_u_id(self, prot_id):
        return self.protein_map[prot_id].get_u_id()

    def get_go_terms(self, prot_id):
        return self.protein_map[prot_id].get_go_terms()

    def get_pathways(self, prot_id):
        return self.protein_map[prot_id].get_pathways()

    def position_in_protein_by_db_id(self, database_id, pos):
        return self.getByDbId(database_id).contains_position(pos)

    def set_position_stored(self, database_id, pos, stored_value):
        self.getByDbId(database_id).set_position_stored(pos, stored_value)

    def set_position_database_id(self, database_id, pos, value):
        self.getByDbId(database_id).set_position_database_id(pos, value)

    def get_position_database_id(self, u_ac, pos):
        return self.protein_map[u_ac].get_position_database_id(pos)

    def get_pos_tags(self, u_ac, pos):
        return self.protein_map[u_ac].get_pos_tags(pos)

    def get_wt_aa(self, u_ac, pos):
        return self.protein_map[u_ac].get_wt_aa(pos)

    def get_mut_aas(self, u_ac, pos):
        return self.protein_map[u_ac].get_mut_aas(pos)

    def get_mut_tags(self, u_ac, pos, new_aa):
        if new_aa is None:
            return None
        return self.protein_map[u_ac].get_mut_tags(pos, new_aa)

    def get_classification(self, prot_id, pos):
        return self.protein_map[prot_id].positions[pos].get_classification()

    def get_structure_list(self):
        return set(self.structures.keys())

    def getStoredStructureIds(self, exclude_interacting_chains = False):
        stored_ids = {}
        for pdb_id in self.structures:
            for chain in self.structures[pdb_id]:
                if exclude_interacting_chains:
                    if self.structures[pdb_id][chain].interacting_structure:
                        continue
                if self.structures[pdb_id][chain].get_stored():
                    stored_ids[self.structures[pdb_id][chain].get_database_id()] = (pdb_id, chain)
        return stored_ids

    def get_protein_database_id(self, u_ac):
        try:
            return self.protein_map[u_ac].get_database_id()
        except:
            return None

    def is_position_stored(self, u_ac, pos):
        return self.protein_map[u_ac].is_position_stored(pos)

    def getAACList(self, u_ac):
        return self.protein_map[u_ac].getAACList()

    def get_res_id(self, u_ac, pos):
        return self.protein_map[u_ac].get_res_id(pos)

    def add_annotation(self, u_ac, pdb_id, chain, anno_obj):
        self.protein_map[u_ac].add_annotation(pdb_id, chain, anno_obj)

    def remove_annotation(self, u_ac, pdb_id, chain):
        try:
            self.protein_map[u_ac].remove_annotation(pdb_id, chain)
        except:
            return

    def remove_structures(self, structure_ids: list[tuple[str, list[str]]]):
        for (structure_id, chains) in structure_ids:
            for chain in chains:
                try:
                    del self.structures[structure_id][chain]
                    if len(self.structures[structure_id]) == 0:
                        del self.structures[structure_id]
                except:
                    continue

    def remove_complexes(self, pdb_ids):
        for pdb_id in pdb_ids:
            del self.complexes[pdb_id]

    def get_target_dict(self, pdb_id):
        target_dict = {}
        chains = self.get_complex_chains(pdb_id)
        for chain in chains:
            if chains[chain] != 'Protein':
                continue
            target_dict[chain] = {}
            if not pdb_id in self.structures:
                continue
            if not chain in self.structres[pdb_id]:
                continue
            mapped_proteins = self.structures[pdb_id][chain].get_mapped_proteins()
            for u_ac in mapped_proteins:
                try:
                    sub_infos = self.get_sub_infos(u_ac, pdb_id, chain)
                except:
                    continue
                for pos in sub_infos:
                    res_nr = sub_infos[pos][0]
                    if res_nr is None:
                        continue
                    if res_nr not in target_dict[chain]:
                        target_dict[chain][res_nr] = []
                    target_dict[chain][res_nr].append((u_ac, pos))
            if len(target_dict[chain]) == 0:
                del target_dict[chain]
        return target_dict

    def classify(self, u_ac, pos, config):
        self.protein_map[u_ac].classify(pos, config)
