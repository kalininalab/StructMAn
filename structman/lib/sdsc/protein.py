import numpy as np

from structman.lib.sdsc.mutations import MultiMutation
import structman.lib.sdsc.structure as structure_package
from structman.lib.sdsc.complex import Complex
from structman.lib.sdsc.sdsc_utils import doomsday_protocol, is_mutant_ac
from structman.lib.sdsc.position import Position
from structman.base_utils.base_utils import Slotted_obj, sequence_hash


# The amino acid sequence is the identity of a protein: one sequence, one entry. The
# identifiers that mapped onto that sequence are all kept (see Protein.identifiers),
# but exactly one of them has to be the primary id, because it is the key of
# Proteins.protein_map and the string that a hundred or so places in the pipeline parse
# for isoform suffixes, RefSeq prefixes, PDB chains and mutant names.
#
# The order below picks it. Fixed rather than "whichever came first" so the surviving
# id cannot change between runs just because the chunking saw the proteins in a
# different order. PDB last: a chain id names a structure, not a protein.
ID_KIND_UNIPROT_AC = 0
ID_KIND_REFSEQ = 1
ID_KIND_ENSEMBL = 2
ID_KIND_OTHER = 3
ID_KIND_PDB = 4

# Kinds as they are stored, one joined string per kind. Keep in step with the columns
# of the Protein table.
ID_TYPE_UNIPROT_AC = 'uniprot_ac'
ID_TYPE_UNIPROT_ID = 'uniprot_id'
ID_TYPE_REFSEQ = 'refseq'
ID_TYPE_ENSEMBL = 'ensembl'
ID_TYPE_PDB = 'pdb'
ID_TYPE_OTHER = 'other'

ID_TYPES = (ID_TYPE_UNIPROT_AC, ID_TYPE_UNIPROT_ID, ID_TYPE_REFSEQ,
            ID_TYPE_ENSEMBL, ID_TYPE_PDB, ID_TYPE_OTHER)

# Separator for the joined strings that reach the database and the output columns.
ID_SEPARATOR = ','

# uniprot.integrate_protein fills Protein.other_ids with these keys. Until now that
# dict reached no table at all; mapping it onto the kinds above is what finally gets it
# stored.
_OTHER_ID_TYPES = {
    'ENSEMBL_ID': ID_TYPE_ENSEMBL,
    'HGNC_ID': ID_TYPE_OTHER,
}


def identifier_kind(identifier: str) -> str | None:
    """Which ID_TYPE a bare identifier string belongs to, or None when unrecognisable.

    Only shapes that are actually distinctive are claimed. None means "no idea", and
    callers must then not store the identifier at all rather than filing it under
    ID_TYPE_OTHER: a fasta entry id has no recognisable shape, and guessing put it in
    Other_Ids when it belongs in RS_Protein_Session.Input_Id and nowhere else.

    In particular a bare uniprot accession is NOT recognisable - nothing separates
    'P12345' from an arbitrary label - so it is only ever filed as one when a caller
    passes ID_TYPE_UNIPROT_AC explicitly, which every uniprot code path does.
    """
    if identifier is None or len(identifier) == 0:
        return None
    if identifier.count(':') == 1 and len(identifier) == 6:
        return ID_TYPE_PDB
    if identifier[:3] in ('NP_', 'NM_', 'XP_', 'XM_'):
        return ID_TYPE_REFSEQ
    if identifier[:3] == 'ENS':
        return ID_TYPE_ENSEMBL
    return None


def join_identifiers(identifiers) -> str | None:
    """Sorted, joined form of an identifier set, or None when it is empty.

    Sorted so the stored string is stable: the same protein must not produce a
    different column value just because the ids arrived in a different order.
    """
    if identifiers is None or len(identifiers) == 0:
        return None
    return ID_SEPARATOR.join(sorted(identifiers))


def split_identifiers(joined: str | None) -> set[str]:
    """Inverse of join_identifiers, for reading a stored column back."""
    if joined is None or len(joined) == 0:
        return set()
    return {part for part in (raw.strip() for raw in joined.split(ID_SEPARATOR)) if len(part) > 0}


def protein_id_kind(protein) -> int:
    """Which of the priority classes above a protein's identity falls into."""
    primary_protein_id = protein.primary_protein_id
    if protein.u_ac is not None:
        return ID_KIND_UNIPROT_AC
    if protein.pdb_id is not None:
        return ID_KIND_PDB
    if protein.ref_id is not None or primary_protein_id[:2] in ('NP', 'NM', 'XP', 'XM'):
        return ID_KIND_REFSEQ
    if primary_protein_id[:3] == 'ENS':
        return ID_KIND_ENSEMBL
    return ID_KIND_OTHER


def canonical_protein_id(protein_map, prot_ids) -> str:
    """The identifier that survives when these proteins are fused into one.

    Sorted on (kind, id) rather than taking the first match, so ties inside one kind
    resolve the same way on every run - "the uniprot accession, alphabetically first".
    """
    return sorted(prot_ids, key=lambda prot_id: (protein_id_kind(protein_map[prot_id]), prot_id))[0]


class Protein(Slotted_obj):
    __slots__ = [
        'primary_protein_id',       'u_ac',                         'u_id',
        'ref_id',                   'ref_nt_id',                    'other_ids',
        'pdb_id',                   'positions',
        'sequence',                 'nucleotide_sequence',          'stored',
        'completely_stored',        'wildtype_protein',             'gene',
        'tags',                     'go_terms',                     'pathways',
        'database_id',              '_structure_annotations',       'annotation_index',
        'input_id',                 'multi_mutations',
        'aggregated_interface_map', 'mutant_type',                  'sav_positions',
        'insertion_positions',      'deletion_flanks',              'number_of_mappings',
        'identifiers',              'input_ids',
        'structure_idx',            'structure_registry',           'annotations_built',
        'sequence_from_input'
        ]

    slot_mask = [
        True, True, True,
        True, True, True,
        True, True,
        False, False, True,
        True, True, False,
        False, False, False,
        True, False, False,
        True, False,
        True, True, True,
        True, True, True,
        True, True,
        False, False, False,
        True
    ]

    def __init__(self, errorlog = None, primary_protein_id=None, u_ac=None, u_id=None, ref_id=None, ref_nt_id = None, wildtype_protein=None, gene = None,
                 pdb_id=None, positions=[], database_id=None, other_ids=[], input_id=None, sequence=None,
                 aggregated_interface_map = {}, mutant_type = None, sav_positions = None, insertion_positions = None, deletion_flanks = None,
                 tags = None
                ):
        self.primary_protein_id: str | None = primary_protein_id
        self.u_ac: str | None = u_ac  # UNIPROT accession number
        self.u_id: str | None = u_id  # UNIPROT ID
        self.ref_id: str | None = ref_id  # RefSeq ID
        self.ref_nt_id: str | None = ref_nt_id
        self.pdb_id: str | None = pdb_id
        self.input_id: str | None = input_id
        self.positions: list[Position] = [None]
        self.sequence = sequence
        # True when the sequence came with the input (fasta), so the input id is only a
        # user label: protCheck then finds the protein by its sequence, never by that id.
        self.sequence_from_input: bool = False
        self.nucleotide_sequence = None
        self.stored = False  # If the protein is stored in the database
        self.completely_stored = False  # Gets true if all corresponding input positions are stored
        self.database_id = database_id  # ID in the corresponding database table
        self.stored = (database_id is not None)
        self.go_terms = {}
        self.pathways = {}
        self._structure_annotations: list[structure_package.StructureAnnotation] = []
        self.annotation_index = None
        # The mmseqs hits as indices into Proteins.annotated_structures, one per
        # (structure, chain), in hit order. See add_structure_idx.
        self.structure_idx: np.ndarray | None = None
        self.structure_registry: list[tuple[str, str]] | None = None
        self.annotations_built: int = 0
        self.multi_mutations = []
        self.wildtype_protein = wildtype_protein
        self.aggregated_interface_map = aggregated_interface_map
        self.mutant_type = mutant_type # Can be None, 'SAV', 'indel', or 'multi'
        self.sav_positions = sav_positions
        self.insertion_positions = insertion_positions
        self.deletion_flanks = deletion_flanks
        self.number_of_mappings = 0

        # Every identifier known to name this sequence, including the primary one, and
        # every input string that mapped onto it. Filled by register_own_identifiers at
        # the end of __init__, once other_ids exists.
        self.identifiers: dict[str, set[str]] = {}
        self.input_ids: set[str] = set()

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

        # Last, so every singular id slot it mirrors is already set.
        self.register_own_identifiers()

    def deconstruct(self):
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            pos_obj.deconstruct(completely = True)

        del self.positions
        del self.sequence
        # Neither the index array nor annotation objects hold reference cycles, so
        # refcounting frees them. Deconstructing each annotation used to be one Python
        # call per hit at the end of every chunk, and would now also build them first.
        self._structure_annotations = None
        self.structure_idx = None
        self.structure_registry = None
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
            if position.pos >= len(self.positions):
                d = position.pos - len(self.positions)
                for _ in range(d):
                    self.positions.append(None)
                self.positions.append(position)
            elif self.positions[position.pos] is None:
                self.positions[position.pos] = position
            elif position.pos is not None:
                err = self.positions[position.pos].fuse(position)
                del_list.append(p)
                if err is not None:
                    raise TypeError(f'Error happened in add_positions fuse: {self.primary_protein_id=} {self.pdb_id=}\n{err}')

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
        """File a PDB-based protein's positions by PDB residue number, not sequence
        position: the sequence position is not known until the structure itself is
        parsed (see serializedPipeline.map_pdb_sequences_onto_positions), so these are
        held in the meantime as a plain, short list past positions[0].

        positions[0] is reserved the same way a non-PDB protein's is: a
        Position(pdb_res_nr=None) is the "tag this chain, no specific residue" sentinel
        that popNone reads, and it goes through add_positions so both protein kinds
        share the same slot-0 convention and popNone needs no pdb_id branch.

        Real residue numbers are appended and, on a repeat of the same residue number -
        two mutations at the same position across separate input lines - fused into the
        existing one, exactly as a repeated sequence position is in add_positions. The
        list stays as short as the input actually names (typically empty for a plain
        structural run, one Protein.__init__ call's worth of explicit mutations
        otherwise), which is what makes a linear scan for that fine: the alternative,
        keeping a dict sized to the whole chain so map_pdb_sequences_onto_positions can
        write one entry per residue into it, was the actual cost - res_id_map used to
        take that for every residue of every chain, though only ever read back for the
        few the input named.
        """
        warns = []
        for position in positions:
            if position.pdb_res_nr is None:
                self.add_positions([position])
                continue
            existing = None
            for pending in self.positions[1:]:
                if pending is not None and pending.pdb_res_nr == position.pdb_res_nr:
                    existing = pending
                    break
            if existing is None:
                self.positions.append(position)
            else:
                err = existing.fuse(position)
                if err is not None:
                    raise TypeError(f'Error happened in add_positions fuse: {self.primary_protein_id=} {self.pdb_id=}\n{err}')
        if len(warns) == 0:
            return None
        else:
            return warns

    def popNone(self):
        if self.positions[0] is not None:
            tags = self.positions[0].pos_tags
            self.positions[0] = None
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
        for annotation in self.structure_annotations:
            score = annotation.sequence_identity * annotation.coverage
            if score > max_score:
                max_score = score
                rec_struct = f'{annotation.pdb_id}:{annotation.chain}'
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

    """
    def getAACList(self) -> dict[str, int]:
        aaclist = {}
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            aac_base = pos_obj.getAACBase()
            aaclist[aac_base] = pos_obj.get_database_id()
        return aaclist
    """
    def getAACList(self) -> list[str]:
        aaclist = []
        for pos_obj in self.positions:
            if pos_obj is None:
                continue
            aac_base = pos_obj.getAACBase()
            aaclist.append(aac_base)
        return aaclist

    def get_aac_base(self, pos):
        return self.positions[pos].getAACBase()

    def status(self):
        stat = '\n'.join([self.u_ac, self.pdb_id, str(len(self.positions)), self.sequence, self.stored, self.database_id, str(self.annotation_count())])
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
        self.add_identifier(other_id, id_type = _OTHER_ID_TYPES.get(id_id, ID_TYPE_OTHER))

    def add_identifier(self, identifier: str | None, id_type: str | None = None) -> None:
        """Record one more identifier for this sequence.

        Without an id_type the kind is taken from the identifier's shape, and an
        identifier whose shape says nothing is dropped rather than filed under
        ID_TYPE_OTHER. Everything worth keeping about such a string is already carried
        by Primary_Protein_Id or by RS_Protein_Session.Input_Id; adding it here as
        "other" only duplicated it into a column it does not belong in.
        """
        if identifier is None or len(identifier) == 0:
            return
        if id_type is None:
            id_type = identifier_kind(identifier)
            if id_type is None:
                return
        if id_type not in self.identifiers:
            self.identifiers[id_type] = set()
        self.identifiers[id_type].add(identifier)

    def add_input_id(self, input_id: str | None) -> None:
        if input_id is None or len(input_id) == 0:
            return
        self.input_ids.add(input_id)

    def register_own_identifiers(self) -> None:
        """Fold this protein's own singular id slots into the identifier set.

        The singular slots stay the primary representative of their kind - that is what
        keeps every id-parsing site in the pipeline working - so this only mirrors them
        into the set, it never changes them.
        """
        self.add_identifier(self.u_ac, ID_TYPE_UNIPROT_AC)
        self.add_identifier(self.u_id, ID_TYPE_UNIPROT_ID)
        self.add_identifier(self.ref_id, ID_TYPE_REFSEQ)
        self.add_identifier(self.pdb_id, ID_TYPE_PDB)
        # The primary id gets a column of its own (Primary_Protein_Id), so it is only
        # worth recording here when its shape names a kind it should be grouped under -
        # an Ensembl or RefSeq id that no singular slot carries. add_identifier drops it
        # when the shape says nothing, which is what keeps fasta entry ids out of
        # Other_Ids: those live in RS_Protein_Session.Input_Id and nowhere else.
        if self.primary_protein_id not in self.all_identifiers():
            self.add_identifier(self.primary_protein_id)
        # Only genuinely other identifiers reach Other_Ids - HGNC and anything future.
        for id_id, other_id in self.other_ids.items():
            self.add_identifier(other_id, _OTHER_ID_TYPES.get(id_id, ID_TYPE_OTHER))
        self.drop_redundant_other_identifiers()
        self.add_input_id(self.input_id)

    def drop_redundant_other_identifiers(self) -> None:
        """Remove anything from the catch-all kind that a real kind already claims.

        Rows written before identifier_kind stopped guessing carry accessions and fasta
        entry ids in Other_Ids. Reading such a row back merges them in under their stored
        kind, so this runs on every write-back and cleans the column as it goes - for the
        duplicates at least; an entry no other kind claims cannot be told apart from a
        genuine other identifier and is left alone.
        """
        other = self.identifiers.get(ID_TYPE_OTHER)
        if not other:
            return
        claimed_elsewhere = set()
        for id_type, identifiers in self.identifiers.items():
            if id_type != ID_TYPE_OTHER:
                claimed_elsewhere |= identifiers
        other -= claimed_elsewhere
        if len(other) == 0:
            del self.identifiers[ID_TYPE_OTHER]

    def get_identifiers(self, id_type: str) -> set[str]:
        return self.identifiers.get(id_type, set())

    def get_joined_identifiers(self, id_type: str) -> str | None:
        return join_identifiers(self.identifiers.get(id_type))

    def get_joined_input_ids(self) -> str | None:
        return join_identifiers(self.input_ids)

    def all_identifiers(self) -> set[str]:
        """Every identifier of every kind, flattened."""
        flat: set[str] = set()
        for identifiers in self.identifiers.values():
            flat |= identifiers
        return flat

    def get_sequence_hash(self):
        """CRC32 of this protein's sequence, or None when it has none yet.

        Recomputed rather than cached in a slot: `sequence` is assigned directly from
        a dozen places in the pipeline (sequenceScan, protCheck, mutate, ...), so a
        cached value would go stale silently. crc32 over one sequence is microseconds,
        and the dedup pass asks once per protein per chunk.
        """
        return sequence_hash(self.sequence)

    def is_duplicate_candidate(self):
        """Whether this protein takes part in sequence deduplication at all.

        The sequence is the identity of a protein, so there is no carve-out by input
        kind: PDB chains and proteins of different genes fuse like anything else.

        Two things still cannot take part. A protein with no sequence has no identity to
        compare - including the 0 a failed fetch leaves, which reaches protCheck now that
        the missing sequences are fetched before it. And mutants are out because their
        sequence only exists after multi_mutation.mutate(), which runs one step after
        protCheck.
        """
        if not isinstance(self.sequence, str) or len(self.sequence) == 0:
            return False
        return not is_mutant_ac(self.primary_protein_id)

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
            return self.positions[pos].stored

    def get_position(self, pos):
        if pos >= len(self.positions):
            return None
        return self.positions[pos]

    def get_position_ids(self) -> list[int]:
        return list(range(1,len(self.positions)))

    def set_position_stored(self, pos, stored_value):
        self.positions[pos].set_stored(stored_value)

    def set_position_database_id(self, pos, value):
        self.positions[pos].set_database_id(value)

    def get_position_database_id(self, pos):
        return self.positions[pos].get_database_id()

    def get_pos_tags(self, pos):
        return self.positions[pos].get_pos_tags()

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

    # structure_annotations is a plain list, not a map keyed by (pdb_id, chain). The
    # key was entirely redundant - StructureAnnotation stores its own u_ac, pdb_id and
    # chain - and paying for it dominated the construction loop in
    # serializedPipeline.autoTemplateSelection: measured on 3M annotations with gc
    # disabled, the dict plus its tuple keys cost 1.40s of a 2.36s loop and 126 of the
    # 272 bytes per annotation. On a whole-PDB chunk, with hundreds of millions of
    # annotations, that is most of the "create object datastructs" step and ~46 GiB.
    #
    # The keyed write also deduplicated repeat (pdb_id, chain) rows. Neither producer
    # can emit one: the mmseqs hit store has a single row per (pdb_id, chain) per
    # protein (MMseqs2._resolve_prot_collisions enforces it) and the database retrieval
    # path is keyed by structure db id, which is unique per protein.
    #
    # Nearly every reader already iterates one protein's annotations and reads the pair
    # back off the object, so they lost a dict lookup rather than gaining a scan. What
    # genuinely needs random access - output/massmodel.py, reaching a single annotation
    # named by a recommended structure - goes through get_annotation, which memoises
    # annotation_index. Nothing on the construction path calls it, so at whole-PDB
    # scale that map is never built.
    #
    # Template selection no longer builds these objects at all. A hit is only ever read
    # back as "which (pdb_id, chain)" - alignment, coverage and sequence identity are set
    # by nothing but database.getAlignments, which has no caller - so the hits are kept
    # as structure_idx, an int32 index into Proteins.annotated_structures. At whole-PDB
    # scale that is 4 bytes per hit instead of a ~120 byte object, and the pipeline's
    # per-hit loops (paraAlignment, the classification bookkeeping) became numpy
    # gathers. The objects are built from the index on first access to
    # structure_annotations, for the readers outside the pipeline (output/, massmodel).
    # Pipeline code must read structure_idx instead: touching structure_annotations
    # there builds one object per hit again.
    @property
    def structure_annotations(self) -> list:
        if self.structure_idx is not None and self.annotations_built < len(self.structure_idx):
            registry = self.structure_registry
            primary_protein_id = self.primary_protein_id
            for index in self.structure_idx[self.annotations_built:].tolist():
                pdb_id, chain = registry[index]
                self._structure_annotations.append(
                    structure_package.StructureAnnotation(primary_protein_id, pdb_id, chain))
            self.annotations_built = len(self.structure_idx)
            self.annotation_index = None
        return self._structure_annotations

    @structure_annotations.setter
    def structure_annotations(self, value: list) -> None:
        self._structure_annotations = value
        self.structure_idx = None
        self.structure_registry = None
        self.annotations_built = 0
        self.annotation_index = None

    @structure_annotations.deleter
    def structure_annotations(self) -> None:
        del self._structure_annotations
        self.structure_idx = None
        self.structure_registry = None
        self.annotations_built = 0

    def add_structure_idx(self, structure_idx: np.ndarray, registry: list[tuple[str, str]]) -> None:
        """Attach hits as indices into registry (Proteins.annotated_structures).

        Appended after any hits the protein already has, so several search databases
        keep the order the object path gave them.
        """
        if self.structure_idx is None or len(self.structure_idx) == 0:
            self.structure_idx = structure_idx
        else:
            self.structure_idx = np.concatenate([self.structure_idx, structure_idx])
        self.structure_registry = registry

    def annotation_count(self) -> int:
        """Number of annotations, without building objects for the index-backed ones."""
        n_indexed = 0 if self.structure_idx is None else len(self.structure_idx)
        return n_indexed - self.annotations_built + len(self._structure_annotations)

    def explicit_annotations(self) -> list:
        """The annotation objects that exist already, without building index-backed ones.

        In the pipeline these are only the stored annotations add_annotation received;
        once something has read structure_annotations it also holds the built ones,
        which are unstored.
        """
        return self._structure_annotations

    def add_annotation(self, anno_obj):
        self.structure_annotations.append(anno_obj)
        self.annotation_index = None

    def get_annotation(self, pdb_id, chain):
        """The annotation for one (pdb_id, chain), or None.

        Only for callers that are not already iterating this protein's annotations.
        The first call builds the index; add_annotation and remove_annotation drop it.
        """
        if self.annotation_index is None:
            self.annotation_index = {
                (annotation.pdb_id, annotation.chain): annotation
                for annotation in self.structure_annotations
            }
        return self.annotation_index.get((pdb_id, chain))

    def remove_annotation(self, pdb_id, chain):
        # A scan, but nothing in the pipeline calls this - the indel path that used to
        # delete annotations went through get_structure_annotations and the live map.
        for index, annotation in enumerate(self.structure_annotations):
            if annotation.pdb_id == pdb_id and annotation.chain == chain:
                del self.structure_annotations[index]
                self.annotation_index = None
                return

    def remove_annotations(self):
        del self.structure_annotations

    def get_structure_annotations(self) -> list:
        """The annotations as a list, NOT a map keyed by (pdb_id, chain).

        indel_analysis.py and sdsc/indel.py still index the result with a tuple and
        delete from it by key. Both are dead today - para_indel_analysis is only
        reached from a commented-out call - and were deliberately left untouched when
        this field stopped being a dict, so reviving either means porting them to the
        list first. get_annotation is the supported random-access path.
        """
        return self.structure_annotations

    def is_annotation_stored(self, pdb_id, chain):
        return self.get_annotation(pdb_id, chain).stored

    def set_coverage(self, pdb_id, chain, value):
        self.get_annotation(pdb_id, chain).set_coverage(value)
        self.number_of_mappings += int(len(self.positions) * value)

    def set_sequence_id(self, pdb_id, chain, value):
        self.get_annotation(pdb_id, chain).set_sequence_id(value)

    def set_annotation_db_id(self, pdb_id, chain, value):
        self.get_annotation(pdb_id, chain).set_database_id(value)

    def set_sub_infos(self, pdb_id, chain, value):
        self.get_annotation(pdb_id, chain).set_sub_infos(value)

    def get_sub_infos(self, pdb_id, chain):
        return self.get_annotation(pdb_id, chain).get_sub_infos()

    def get_backmap(self, structure_id, chain):
        # The map version tested `if chain in ...: return {}`, so it returned empty
        # exactly when the annotation existed. Correcting that exposed the real
        # problem underneath: backmaps live on Structure (Structure.backmaps), and
        # StructureAnnotation has no get_backmap at all, so the inverted guard was
        # the only reason this never raised. Left delegating rather than quietly
        # returning {}, so that a caller gets a loud failure naming the real issue
        # instead of a silently wrong empty map. Nothing calls it today - reach
        # backmaps through Proteins.structures[pdb_id][chain].backmaps instead.
        annotation = self.get_annotation(structure_id, chain)
        if annotation is None:
            return {}
        return annotation.get_backmap()

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

                if mutant_positions[1] is not None or mutant_positions[2] is not None:
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
                 'not_stored_ids', 'id_map', 'structures', 'complexes', 'max_internal_struct_db_id',
                 'indels', 'multi_mutations', 'multi_mutation_back_map', 'seq_store', 'database_ids',
                 'annotated_structures', 'annotated_structure_index']

    def __init__(self, proteins, indels, multi_mutation_objects):
        self.protein_map = proteins
        self.stored_ids = set()
        self.stored_ids_mutant_excluded = set()
        self.completely_stored_ids = set()
        self.not_stored_ids = set()
        self.id_map = {}
        self.structures: dict[str, dict[str, structure_package.Structure]] = {}
        self.complexes: dict[str, Complex] = {}
        # Every (pdb_id, chain) some protein's structure_idx points at, and its index.
        # Pairs rather than Structure objects, so a structure that structureCheck or the
        # complex parsing replaces is still resolved through self.structures.
        self.annotated_structures: list[tuple[str, str]] = []
        self.annotated_structure_index: dict[tuple[str, str], int] = {}
        self.indels = {}
        self.multi_mutations = {}
        self.seq_store = None
        self.database_ids = None
        for multi_mutation_obj in multi_mutation_objects:
            if multi_mutation_obj.wt_prot not in self.multi_mutations:
                self.multi_mutations[multi_mutation_obj.wt_prot] = []
            self.multi_mutations[multi_mutation_obj.wt_prot].append(multi_mutation_obj)
        self.multi_mutation_back_map = {}
        for prot_id in indels:
            self.indels[prot_id] = {}
            for indel_notation in indels[prot_id]:
                self.indels[prot_id][indel_notation] = indels[prot_id][indel_notation]

    def __getitem__(self, key: str) -> Protein:
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
        del self.annotated_structures
        del self.annotated_structure_index
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

    def _release_transferred_positions(self, canonical: Protein, duplicate: Protein) -> None:
        """Take every Position the survivor now holds off the duplicate's books.

        add_positions and add_residues store Position objects by reference, they do not
        copy them, so after a fusion the same objects hang in both proteins.

        remove_protein deconstructs the duplicate right after, and Protein.deconstruct
        runs pos_obj.deconstruct(completely = True) over its positions - a full
        doomsday_protocol that drops every slot, wt_aa included. Without this release
        the survivor is left holding gutted Positions: the next fusion into it dies on
        'Position' object has no attribute 'wt_aa', and anything that reads them earlier
        gets nonsense.

        What the duplicate keeps is exactly what the survivor did not take - the
        positions that were merged into a Position of its own, and any residue that was
        dropped - and those are still freed eagerly by remove_protein.
        """
        owned = {id(position) for position in canonical.positions if position is not None}

        if len(owned) == 0:
            return

        duplicate.positions = [
            None if (position is not None and id(position) in owned) else position
            for position in duplicate.positions
        ]

    @staticmethod
    def _position_conflict(canonical: Protein, duplicate: Protein) -> str | None:
        """The first position of duplicate that cannot be fused into canonical, or None."""
        for position in duplicate.positions:
            if position is None or position.pos is None or position.pos >= len(canonical.positions):
                continue
            existing = canonical.positions[position.pos]
            if existing is None:
                continue
            conflict = existing.conflict_with(position)
            if conflict is not None:
                return conflict
        return None

    def _merge_duplicate_into(self, canonical: Protein, duplicate: Protein, config) -> bool:
        """Move everything a duplicate asked for onto the survivor. Returns False, having
        changed nothing, when their positions disagree.

        Content only: references elsewhere in the chunk are left to
        _rekey_references, the survivor's identifier set is re-registered by the caller,
        and the duplicate is not removed. That split is what lets a whole pass of fusions
        share one reference sweep.
        """
        duplicate_id = duplicate.primary_protein_id
        canonical_id = canonical.primary_protein_id

        # All or nothing. Identical sequences cannot disagree on a wildtype, so a conflict
        # means the two are not the same protein - in practice a match by a shared
        # identifier onto a different sequence. Merging up to the conflict and dropping the
        # rest left the duplicate's multi mutations naming SNVs the survivor never got,
        # and insertMultiMutations failed on them hours later.
        conflict = self._position_conflict(canonical, duplicate)
        if conflict is not None:
            config.errorlog.add_warning(
                f'Not fusing {duplicate_id} into {canonical_id}: {conflict}. Both stay separate entries.'
            )
            return False

        # add_positions mutates and shortens the list it is handed, so it must never
        # see the duplicate's own list - it is still attached to the object that is
        # about to be deconstructed.
        transferred = [position for position in duplicate.positions if position is not None]
        if len(transferred) > 0:
            canonical.add_positions(transferred)

        canonical.tags = canonical.tags | duplicate.tags
        canonical.sequence_from_input = canonical.sequence_from_input or duplicate.sequence_from_input

        for id_id, other_id in duplicate.other_ids.items():
            canonical.add_other_ids(id_id, other_id)

        # Fill only what the survivor is missing; its own identity always wins, because
        # it is the primary id every id-parsing site in the pipeline reads.
        for id_slot in ('u_ac', 'u_id', 'ref_id', 'ref_nt_id', 'gene'):
            if getattr(canonical, id_slot) is None:
                setattr(canonical, id_slot, getattr(duplicate, id_slot))

        # The point of the whole exercise: the duplicate's identifiers are kept, not
        # thrown away with it. Its own singular slots were folded into its identifier
        # set at construction, so this one union carries all of them.
        for id_type, identifiers in duplicate.identifiers.items():
            for identifier in identifiers:
                canonical.add_identifier(identifier, id_type)
        canonical.input_ids |= duplicate.input_ids

        # Must come before the duplicate is removed: removal deconstructs it, and the
        # Positions just handed to the survivor are still reachable from it.
        self._release_transferred_positions(canonical, duplicate)
        return True

    def apply_identity_changes(
            self,
            renames: dict[str, str],
            fusions: dict[str, str],
            config,
            genes = None) -> dict[str, str]:
        """Apply a whole pass of renames and fusions with one reference sweep.

        renames: old primary id -> new primary id, for a protein taking over an id the
        chunk does not hold yet (a stored protein's). fusions: duplicate id -> survivor
        id. A rename onto an id the chunk already holds is applied as a fusion instead.

        Fusing one protein at a time scanned every protein, multi mutation, indel and
        gene for references to the dropped id, once per fusion. A whole-PDB chunk of
        300,000 proteins with 145,314 identical-sequence fusions spent over an hour in
        protCheck on that alone. Here the content merges run first, then one sweep
        rewrites every reference through a dict lookup, then the duplicates are removed,
        and each survivor's identifier set is registered once rather than once per
        protein fused into it.

        Returns the applied id -> final id mapping.
        """
        fusions = dict(fusions)
        mapping: dict[str, str] = {}
        touched: set[str] = set()

        for old_id, new_id in renames.items():
            if old_id == new_id or old_id not in self.protein_map:
                continue
            if new_id in self.protein_map:
                fusions[old_id] = new_id
                continue
            protein = self.protein_map.pop(old_id)
            # The old id is still a name for this sequence. It keeps whatever identifier
            # column its kind earned it, and reaches the output either way through
            # input_ids; the new id takes over Primary_Protein_Id.
            protein.add_identifier(old_id)
            protein.primary_protein_id = new_id
            self.protein_map[new_id] = protein
            mapping[old_id] = new_id
            touched.add(new_id)

        def resolve(identifier: str) -> str:
            # Follows rename and fusion chains to the id that survives. A cycle cannot
            # come out of the passes that call this; the guard just keeps a malformed
            # request from looping.
            seen = set()
            while identifier not in seen:
                seen.add(identifier)
                if identifier in mapping:
                    identifier = mapping[identifier]
                elif identifier in fusions and fusions[identifier] != identifier:
                    identifier = fusions[identifier]
                else:
                    break
            return identifier

        removed: list[str] = []
        for duplicate_id, canonical_id in list(fusions.items()):
            final_id = resolve(canonical_id)
            if final_id == duplicate_id:
                continue
            if duplicate_id not in self.protein_map or final_id not in self.protein_map:
                continue
            if not self._merge_duplicate_into(self.protein_map[final_id], self.protein_map[duplicate_id], config):
                # Refused: the duplicate stays an entry of its own, so nothing may resolve
                # through it to the survivor any more.
                del fusions[duplicate_id]
                continue
            mapping[duplicate_id] = final_id
            removed.append(duplicate_id)
            touched.add(final_id)

        for identifier in mapping:
            mapping[identifier] = resolve(mapping[identifier])

        self._rekey_references(mapping, genes = genes)

        for duplicate_id in removed:
            self.remove_protein(duplicate_id)

        for identifier in touched:
            if identifier in self.protein_map:
                # Idempotent; files anything the slot backfill moved onto the survivor.
                self.protein_map[identifier].register_own_identifiers()

        return mapping

    def fuse_protein(self, canonical_id: str, duplicate_id: str, config, genes = None) -> None:
        """Merge duplicate_id into canonical_id and remove it.

        Called when two proteins turned out to carry the same amino acid sequence. The
        duplicate's identifier disappears from the run entirely - no Protein row, no
        RS_Protein_Session row, absent from the output - but everything it *asked for*
        (positions, tags, alternative ids) is merged into the survivor first, so a
        position requested only under the duplicate id is still answered.

        One fusion costs one reference sweep over the chunk. For a whole pass use
        apply_identity_changes, which shares a single sweep between all of them.
        """
        if canonical_id == duplicate_id:
            return
        if duplicate_id not in self.protein_map or canonical_id not in self.protein_map:
            return
        self.apply_identity_changes({}, {duplicate_id: canonical_id}, config, genes = genes)

    def _rekey_references(self, mapping: dict[str, str], genes = None) -> None:
        """Point everything that named an id in mapping at the id it maps to, in one sweep."""
        if len(mapping) == 0:
            return

        # Mutant and indel objects of a dropped id become mutants of its survivor.
        # Both dicts are keyed by the wildtype's primary id.
        for store in (self.multi_mutations, self.indels):
            for old_id in [key for key in store if key in mapping]:
                moved = store.pop(old_id)
                new_id = mapping[old_id]
                if new_id not in store:
                    store[new_id] = moved
                elif isinstance(moved, dict):
                    store[new_id].update(moved)
                else:
                    store[new_id] += moved

        for multi_mutation_list in self.multi_mutations.values():
            for multi_mutation in multi_mutation_list:
                if multi_mutation.wt_prot in mapping:
                    multi_mutation.wt_prot = mapping[multi_mutation.wt_prot]
                if multi_mutation.mut_prot in mapping:
                    multi_mutation.mut_prot = mapping[multi_mutation.mut_prot]

        for indel_map in self.indels.values():
            for indel in indel_map.values():
                if indel.wt_prot in mapping:
                    indel.wt_prot = mapping[indel.wt_prot]
                if indel.mut_prot in mapping:
                    indel.mut_prot = mapping[indel.mut_prot]

        # Mutants carry their wildtype as a plain string.
        for protein in self.protein_map.values():
            if protein.wildtype_protein in mapping:
                protein.wildtype_protein = mapping[protein.wildtype_protein]

        # Gene.proteins is input_id -> primary_protein_id and would otherwise dangle.
        if genes is not None:
            for gene in genes.values():
                for input_id, primary_protein_id in gene.proteins.items():
                    if primary_protein_id in mapping:
                        gene.proteins[input_id] = mapping[primary_protein_id]

    def _rekey_protein_references(self, canonical_id: str, duplicate_id: str, genes = None) -> None:
        """Point everything that named duplicate_id at canonical_id instead."""
        self._rekey_references({duplicate_id: canonical_id}, genes = genes)

    def rename_protein(self, old_id: str, new_id: str, genes = None) -> None:
        """Give a protein a different primary id, keeping every reference to it valid.

        Used when a chunk protein turns out to duplicate one already in the database:
        it takes over the stored protein's primary id so its annotations are reported
        under the identifier the database already knows. For a whole pass use
        apply_identity_changes.
        """
        if old_id == new_id or old_id not in self.protein_map:
            return

        protein = self.protein_map.pop(old_id)
        # The old id is still a name for this sequence. It keeps whatever identifier
        # column its kind earned it, and reaches the output either way through
        # input_ids; the new id takes over Primary_Protein_Id.
        protein.add_identifier(old_id)
        protein.primary_protein_id = new_id
        protein.register_own_identifiers()
        self.protein_map[new_id] = protein

        self._rekey_references({old_id: new_id}, genes = genes)

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
        return self.complexes[pdb_id].stored

    def set_complex_db_id(self, pdb_id, value):
        self.complexes[pdb_id].set_database_id(value)

    def get_complex_chains(self, pdb_id, only_protein=False):
        return self.complexes[pdb_id].get_chains(only_protein=only_protein)

    def set_IAmap(self, structure_id, IAmap):
        self.complexes[structure_id].set_IAmap(IAmap)

    def get_IAmap(self, structure_id):
        return self.complexes[structure_id].get_IAmap()

    def get_interfaces(self, structure_id):
        return self.complexes[structure_id].get_interfaces()

    def print_interfaces(self):
        for structure_id in self.complexes:
            print(structure_id)
            self.complexes[structure_id].print_interfaces()

    def get_backmap(self, protein_id, structure_id, chain):
        return self.protein_map[protein_id].get_backmap(structure_id, chain)

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

    def get_complex_structures(self, pdb_id: str) -> dict[str, dict[str, structure_package.Structure]]:
        chains = self.get_complex_chains(pdb_id)
        structures = {pdb_id : {}}
        for chain in chains:
            if pdb_id not in self.structures:
                continue
            if chain not in self.structures[pdb_id]:
                continue
            structures[pdb_id][chain] = self.structures[pdb_id][chain]
        return structures

    def protein_chain_count(self, pdb_id: str) -> int:
        """The protein chains of a complex that have a stored Structure.

        These are the chains compl_struct_db_id_matrix holds for the complex during
        aggregation, i.e. the chains an interface can be found between, so a complex
        with fewer than two has no protein-protein interface.
        """
        compl_obj = self.complexes.get(pdb_id)
        if compl_obj is None:
            return 0
        struct_chains = self.structures.get(pdb_id, {})
        chain_count = 0
        for chain, chain_type in compl_obj.chains.items():
            if chain_type != 'Protein':
                continue
            struct = struct_chains.get(chain)
            if struct is None or struct.database_id is None:
                continue
            chain_count += 1
        return chain_count

    def get_structures(self):
        return self.structures

    def annotated_structure_objects(self) -> list[structure_package.Structure]:
        """The Structure behind every annotated_structures entry, in index order.

        Resolved now, not cached: structureCheck and the complex parsing set database
        ids on these objects, so per-structure columns must be read after them.
        """
        structures = self.structures
        return [structures[pdb_id][chain] for pdb_id, chain in self.annotated_structures]

    def add_structure(self, pdb_id, chain, struct_obj):
        if pdb_id not in self.structures:
            self.structures[pdb_id] = {}
        self.structures[pdb_id][chain] = struct_obj
        #self.complexes[pdb_id].chains[chain] = 'Protein'

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
        return self.structures[pdb_id][chain].stored

    def contains_residue(self, pdb_id, chain, res_nr):
        if pdb_id not in self.structures:
            return False
        if chain not in self.structures[pdb_id]:
            return False
        return self.structures[pdb_id][chain].contains_residue(res_nr)

    def get_protein_structure_annotations(self, u_ac):
        return self.protein_map[u_ac].get_structure_annotations()

    def add_mapping_to_structure(self, pdb_id, chain, u_ac):
        self.structures[pdb_id][chain].add_mapping(u_ac)

    def is_annotation_stored(self, pdb_id, chain, u_ac):
        return self.protein_map[u_ac].is_annotation_stored(pdb_id, chain)

    def set_coverage_by_db_id(self, prot_id, pdb_id, chain, value):
        self.getByDbId(prot_id).set_coverage(pdb_id, chain, value)

    def set_coverage(self, u_ac, pdb_id, chain, value):
        self.protein_map[u_ac].set_coverage(pdb_id, chain, value)

    def set_sequence_id_by_db_id(self, prot_id, pdb_id, chain, value):
        self.getByDbId(prot_id).set_sequence_id(pdb_id, chain, value)

    def set_sequence_id(self, u_ac, pdb_id, chain, value):
        self.protein_map[u_ac].set_sequence_id(pdb_id, chain, value)

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

    def getStoredStructureIds(self, exclude_interacting_chains: bool = False) -> dict[int, tuple[str, str]]:
        stored_ids: dict[int, tuple[str, str]] = {}
        pdb_id: str
        chain: str
        for pdb_id in self.structures:
            for chain in self.structures[pdb_id]:
                if exclude_interacting_chains:
                    if self.structures[pdb_id][chain].interacting_structure:
                        continue
                db_id = self.structures[pdb_id][chain].database_id
                if db_id is not None:
                    stored_ids[db_id] = (pdb_id, chain)
        return stored_ids

    def get_protein_database_id(self, u_ac):
        try:
            return self.protein_map[u_ac].get_database_id()
        except:
            return None

    def is_position_stored(self, u_ac, pos):
        return self.protein_map[u_ac].is_position_stored(pos)

    def get_res_id(self, u_ac, pos):
        return self.protein_map[u_ac].get_res_id(pos)

    def add_annotation(self, u_ac, anno_obj):
        self.protein_map[u_ac].add_annotation(anno_obj)

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
