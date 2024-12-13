from structman.lib import globalAlignment, pdbParser
from structman.lib.sdsc.sdsc_utils import process_alignment_data, doomsday_protocol


class Structure:
    __slots__ = ['pdb_id', 'chain', 'oligo', 'database_id', 'stored', 'mapped_proteins', 'residues', 'last_residue', 'first_residue', 'sequence', 'seq_len', 'new_interacting_chain', 'interacting_structure']

    def __init__(self, pdb_id, chain, oligo=set(), mapped_proteins = None, database_id=None, last_residue=None, first_residue=None, sequence=None, seq_len = None, new_interacting_chain = False):
        self.pdb_id = pdb_id
        self.chain = chain
        self.database_id = database_id
        if isinstance(oligo, str):
            self.oligo = set(oligo)
        else:
            self.oligo = oligo.copy()
        self.stored = (database_id is not None)
        if mapped_proteins is None:
            self.mapped_proteins = []
        else:
            self.mapped_proteins = mapped_proteins[:]
        self.residues = {}
        self.last_residue = last_residue
        self.first_residue = first_residue
        self.sequence = sequence
        self.seq_len = seq_len
        self.new_interacting_chain = new_interacting_chain
        self.interacting_structure = False

    def deconstruct(self):
        del self.mapped_proteins
        for res in self.residues:
            self.residues[res].deconstruct()
        del self.residues
        del self.sequence
        doomsday_protocol(self)

    def parse_page(self, page, config):
        seq_res_map, seq, last_residue, first_residue = globalAlignment.createTemplateFasta(page, self.pdb_id, self.chain, config, seqAndMap=True, could_be_empty=True)
        self.sequence = seq
        self.first_residue = first_residue
        self.last_residue = last_residue

    def status(self):
        stat = '\n'.join([self.pdb_id, self.chain, str(self.resolution), str(self.oligo), str(self.database_id), str(self.stored), str(len(self.residues))])
        return stat

    def add_mapping(self, u_ac):
        self.mapped_proteins.append(u_ac)

    def get_mapped_proteins(self):
        return self.mapped_proteins

    def get_mapped_positions(self, res_nr, proteins):

        #if len(self.mapped_proteins) > 0:
        #    print(f'Call of get_mapped_positions for {self.pdb_id} {self.chain} {res_nr}\nMapped proteins: {self.mapped_proteins}')

        mapped_positions = {}
        for protein_id in self.mapped_proteins:
            backmap = proteins.get_backmap(protein_id, self.pdb_id, self.chain)
            if res_nr in backmap:
                mapped_positions[protein_id] = backmap[res_nr]
        return mapped_positions

    def set_database_id(self, value):
        self.database_id = value

    def get_database_id(self):
        return self.database_id

    def set_oligo(self, value):
        self.oligo = value

    def set_last_residue(self, last_residue):
        self.last_residue = last_residue

    def set_first_residue(self, first_residue):
        self.first_residue = first_residue

    def get_oligo(self):
        return self.oligo

    def set_stored(self, value):
        self.stored = value

    def get_stored(self):
        return self.stored

    def add_residue(self, res_nr, residue_obj):
        self.residues[res_nr] = residue_obj

    def get_residue_list(self):
        return self.residues.keys()

    def get_last_residue(self):
        return self.last_residue

    def get_first_residue(self):
        return self.first_residue

    def set_residue_db_id(self, res_nr, value):
        self.residues[res_nr].set_database_id(value)

    def contains_residue(self, res_nr):
        return res_nr in self.residues

    def get_residue_db_id(self, res_nr):
        return self.residues[res_nr].get_database_id()

    def get_residue_aa(self, res_nr):
        try:
            return self.residues[res_nr].get_aa()
        except:
            #print(f'Length of self.residues in get_residue_aa({res_nr}): {len(self.residues)}')
            return None

    def get_residue_sld(self, res_nr):
        return self.residues[res_nr].get_ligand_distances()

    def get_residue_scd(self, res_nr):
        return self.residues[res_nr].get_chain_distances()

    def get_residue_homomer_dists(self, res_nr):
        return self.residues[res_nr].get_homomer_dists()

    def get_residue_centralities(self, res_nr, get_whats_there=False):
        return self.residues[res_nr].get_centralities(get_whats_there=get_whats_there)

    def get_residue_modres(self, res_nr):
        return self.residues[res_nr].get_modres()

    def get_residue_b_factor(self, res_nr):
        return self.residues[res_nr].get_b_factor()

    def get_residue_rsa(self, res_nr):
        return self.residues[res_nr].get_rsa()

    def get_residue_rsa_triple(self, res_nr):
        return self.residues[res_nr].get_rsa(splitted=True)

    def get_residue_ssa(self, res_nr):
        return self.residues[res_nr].get_ssa()

    def get_residue_phi(self, res_nr):
        return self.residues[res_nr].get_phi()

    def get_residue_psi(self, res_nr):
        return self.residues[res_nr].get_psi()

    def get_residue_link_information(self, res_nr):
        return self.residues[res_nr].get_residue_link_information()

    def get_residue_interaction_profile(self, res_nr, get_whats_there=False):
        return self.residues[res_nr].get_interaction_profile(get_whats_there=get_whats_there)

    def get_residue_interaction_profile_str(self, res_nr):
        return self.residues[res_nr].get_interaction_profile_str()

    def get_residue_milieu(self, res_nr):
        return self.residues[res_nr].get_milieu()

    def add_residue_classification(self, res_nr, Class, simpleClass):
        self.residues[res_nr].set_classification(Class, simpleClass)

    def getSequence(self, config, complex_obj=None, for_modeller=False):
        if self.sequence is not None and not for_modeller:
            return self.sequence
        else:
            if complex_obj is None:
                template_page, atom_count = pdbParser.standardParsePDB(self.pdb_id, config.pdb_path)
            else:
                template_page = complex_obj.getPage(config)
            seq_res_map, seq, last_residue, first_residue = globalAlignment.createTemplateFasta(template_page, self.pdb_id, self.chain, config, seqAndMap=True, for_modeller=for_modeller, could_be_empty=True, rare_residues=config.rare_residues)
            self.sequence = seq
        return self.sequence

    def get_seq_len(self):
        if self.seq_len is not None:
            return self.seq_len
        elif self.sequence is not None:
            self.seq_len = len(self.sequence)
            return self.seq_len
        else:
            return None


class StructureAnnotation:
    __slots__ = ['u_ac', 'pdb_id', 'chain', 'alignment', 'coverage', 'sequence_identity', 'sub_infos', 'stored', 'database_id', 'backmap']

    def __init__(self, u_ac, pdb_id, chain, alignment=None, stored=False, backmap = {}):
        self.u_ac = u_ac
        self.pdb_id = pdb_id
        self.chain = chain
        self.alignment = alignment
        self.coverage = None
        self.sequence_identity = None
        self.sub_infos = {}  # {pos:(res_nr,res_aa,structure_sequence_number)}
        self.stored = stored
        self.database_id = None
        self.backmap = backmap

    def deconstruct(self):
        del self.alignment
        del self.sub_infos
        del self.backmap
        doomsday_protocol(self)

    def set_alignment(self, value):
        self.alignment = value

    def get_alignment(self):
        if isinstance(self.alignment, tuple):
            return self.alignment
        return process_alignment_data(self.alignment)

    def pop_alignment(self):
        aln = self.alignment
        self.alignment = None
        if isinstance(aln, tuple):
            return aln
        return process_alignment_data(aln)

    def set_coverage(self, value):
        self.coverage = value

    def get_coverage(self):
        return self.coverage

    def set_sequence_id(self, value):
        self.sequence_identity = value

    def get_sequence_id(self):
        return self.sequence_identity

    def set_database_id(self, value):
        self.database_id = value

    def get_stored(self):
        return self.stored

    def set_sub_infos(self, value):
        self.sub_infos = value

    def get_sub_infos(self):
        return self.sub_infos

    def get_sub_info(self, pos):
        try:
            return self.sub_infos[pos]
        except:
            return None

    def set_backmap(self, backmap):
        self.backmap = backmap

    def get_backmap(self):
        return self.backmap

    def is_covered(self, pos):
        if pos not in self.sub_infos:
            return False
        if self.sub_infos[pos][0] is None:
            return False
        return True

    def adjacency_check(self, left_pos, right_pos):
        return ((self.sub_infos[right_pos][2] - self.sub_infos[left_pos][2]) == 1)

    def is_terminal(self, pos):
        if pos not in self.sub_infos:
            return False
        structure_sequence_number = self.sub_infos[pos][2]
        if structure_sequence_number == 1:
            return True
        if self.alignment is None:
            print('Alignment not found:', self.u_ac, self.pdb_id, self.chain)
        target_seq, template_seq = process_alignment_data(self.alignment)
        template_seq = template_seq.replace('-', '')
        if structure_sequence_number == len(template_seq):
            return True
        return False

    def model(self):
        return
