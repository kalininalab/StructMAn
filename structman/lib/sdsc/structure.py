from structman.lib import globalAlignment, pdbParser
from structman.lib.sdsc.sdsc_utils import process_alignment_data, doomsday_protocol, Slotted_obj
from structman.lib.sdsc.residue import Residue_Map
from structman.base_utils.base_utils import unpack

class Structure(Slotted_obj):
    __slots__ = [
        'pdb_id',       'chain',        'oligo',
        'database_id',  'stored',       'mapped_proteins',
        'residues',     'last_residue', 'first_residue',
        'sequence',     'seq_len',      'new_interacting_chain',
        'interacting_structure', 'backmaps'
        ]
    
    slot_mask = [
        True, True, True,
        True, True, True,
        True, True, True,
        False, True, True,
        True, True
    ]

    def __init__(self, pdb_id = None, chain = None, oligo = None, mapped_proteins = None, database_id=None, last_residue=None, first_residue=None, sequence=None, seq_len = None, new_interacting_chain = False):
        self.pdb_id = pdb_id
        self.chain = chain
        self.database_id = database_id

        if oligo is None:
            self.oligo: set[str] = set()
        elif isinstance(oligo, str):
            self.oligo = set()
            for ol in oligo:
                self.oligo.add(ol)
        else:
            self.oligo = oligo.copy()

        self.stored = (database_id is not None)
        if mapped_proteins is None:
            self.mapped_proteins = []
        else:
            self.mapped_proteins = mapped_proteins[:]
        self.residues = Residue_Map()
        self.last_residue = last_residue
        self.first_residue = first_residue
        self.sequence = sequence
        self.seq_len = seq_len
        self.new_interacting_chain = new_interacting_chain
        self.interacting_structure = False
        self.backmaps: dict[str, Residue_Map[int]] = {}

    def deconstruct(self):
        del self.mapped_proteins
        
        del self.residues
        del self.sequence
        del self.backmaps
        doomsday_protocol(self)

    def parse_page(self, page: bytes, config):
        seq_res_map, seq, last_residue, first_residue = globalAlignment.createTemplateFasta(page.split(b'\n'), self.pdb_id, self.chain.encode(), config, seqAndMap=True, could_be_empty=True)
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

    def set_database_id(self, value):
        self.database_id = value

    def get_database_id(self):
        return self.database_id

    def set_last_residue(self, last_residue):
        self.last_residue = last_residue

    def set_first_residue(self, first_residue):
        self.first_residue = first_residue

    def set_stored(self, value):
        self.stored = value

    def get_stored(self):
        return self.stored

    def add_residue(self, res_nr, residue_obj):
        self.residues.add_item(res_nr, residue_obj)


    def set_residue_db_id(self, res_nr, value):
        self.residues.get_item(res_nr).set_database_id(value)

    def contains_residue(self, res_nr):
        return self.residues.contains(res_nr)

    def getSequence(self, config, complex_obj=None, for_modeller=False):
        if self.sequence is not None and not for_modeller:
            return self.sequence
        else:
            if complex_obj is None:
                template_page, atom_count = pdbParser.standardParsePDB(self.pdb_id, config, return_bytes=True)
            else:
                template_page = complex_obj.getPage(config)
            seq_res_map, seq, last_residue, first_residue = globalAlignment.createTemplateFasta(template_page.split(b'\n'), self.pdb_id, self.chain.encode('ascii'), config, seqAndMap=True, for_modeller=for_modeller, could_be_empty=True, rare_residues=config.rare_residues)
            
            if config.verbosity >= 5:
                print(f'In getSequence of {self.pdb_id=} {self.chain=} {len(seq)=}')
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


class StructureAnnotation(Slotted_obj):
    __slots__ = [
        'u_ac',         'pdb_id',   'chain',
        'alignment',    'coverage', 'sequence_identity',
        'sub_infos',    'stored',   'database_id'
        ]

    slot_mask = [
        True, True, True,
        False, True, True,
        True, True, True
        ]

    def __init__(self, u_ac = None, pdb_id = None, chain = None, alignment=None, stored=False):
        self.u_ac = u_ac
        self.pdb_id = pdb_id
        self.chain = chain
        self.alignment: bytes | None = alignment
        self.coverage = None
        self.sequence_identity = None
        self.sub_infos = []  # {pos:(res_nr,res_aa,structure_sequence_number)}
        self.stored = stored
        self.database_id: int | None = None

    def deconstruct(self):
        del self.alignment
        del self.sub_infos
        doomsday_protocol(self)

    def get_alignment(self, unpacked = False) -> bytes | None | tuple[str, str]:   
        if isinstance(self.alignment, bytes) and unpacked:
            alignment = unpack(self.alignment)
            return process_alignment_data(alignment)
        return self.alignment

    def pop_alignment(self, unpacked = False) -> bytes | None | tuple[str, str]:
        aln = self.alignment
        self.alignment = None
        if isinstance(aln, bytes) and unpacked:
            aln = unpack(aln)
            return process_alignment_data(aln)
        return aln
        
    def set_coverage(self, value):
        self.coverage = value

    def set_sequence_id(self, value):
        self.sequence_identity = value

    def get_sequence_id(self):
        return self.sequence_identity

    def set_database_id(self, value):
        self.database_id = value

    def set_sub_infos(self, value):
        self.sub_infos = value

    def get_sub_infos(self):
        return self.sub_infos

    def get_sub_info(self, pos):
        try:
            return self.sub_infos[pos]
        except:
            return None

    def is_covered(self, pos):
        if self.sub_infos[pos] is None:
            return False
        if self.sub_infos[pos][0] is None:
            return False
        return True

    def adjacency_check(self, left_pos, right_pos):
        return ((self.sub_infos[right_pos][2] - self.sub_infos[left_pos][2]) == 1)

    def is_terminal(self, pos):
        if self.sub_infos[pos] is None:
            return False
        structure_sequence_number = self.sub_infos[pos][2]
        if structure_sequence_number == 1:
            return True
        if self.alignment is None:
            print('Alignment not found:', self.u_ac, self.pdb_id, self.chain)
        target_seq, template_seq = self.get_alignment()
        template_seq = template_seq.replace('-', '')
        if structure_sequence_number == len(template_seq):
            return True
        return False

    def model(self):
        return
