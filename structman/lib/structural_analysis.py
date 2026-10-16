"""mmCIF structures for the structural analysis: header data, the parsed structure and DSSP.

One residue list per chain. The residues of a kept protein chain are the observed residues of
`_pdbx_poly_seq_scheme`, read with mmcifParser.select_protein_chains - the code that also built the
search DB, the PDB input sequences and the alignment templates (`Structure_N.Parse_Data`). Their order
is the contig order everywhere: row r of every per-residue column (DSSP, RIN profile, Residue_Data)
and contig r+1 of the alignments and interfaces are the same residue. Residue ids follow the
convention of the PDB format pipeline: the author number as int, or with its insertion code as str
('52A').

Chains are the author chain ids of the file, of any length (`A`, `AA`, `A-2` for a symmetry copy).
Like the search DB and the RINs, the analysis is about the kept instances - one per chain of the
deposited entry (mmcifParser.select_protein_chains, rin_engine.structure.kept_label_asyms). Symmetry
copies are only context: they keep their own chain id, point to the kept chain they copy
(`original_chain`) and are cut down to the residues near the kept instances.

The parsed structure is `Parsed_structure`: chain-number indexed columns as the PDB format parser
(templateFiltering.parsePDB) produced them, with chain numbers that are positions in the structure's
own chain list rather than in a fixed alphabet.
"""

from __future__ import annotations

import os
import subprocess
import uuid
from dataclasses import dataclass, field

import numpy as np

from structman.base_utils.base_utils import Residue_Map, Slotted_obj
from structman.base_utils.chain_consts import CHAIN_ID_ORDER
from structman.lib import mmcifParser
from structman.lib import structure_source

# ---------------------------------------------------------------------------------------------
# Chain types, as the PDB format parser named them

CHAIN_TYPE_NAMES = ('Protein', 'RNA', 'DNA', 'Peptide', 'MISC')
CHAIN_PROTEIN, CHAIN_RNA, CHAIN_DNA, CHAIN_PEPTIDE, CHAIN_MISC = range(5)
CHAIN_TYPE_CODES = {name: code for code, name in enumerate(CHAIN_TYPE_NAMES)}

# A polypeptide chain is a Peptide below this many observed residues - about the 300 heavy atoms
# below which the PDB format parser called a chain a peptide - or when at most two thirds of its
# residues are standard amino acids.
PEPTIDE_MAX_OBSERVED_RESIDUES = 37

STANDARD_AMINO_ACIDS = frozenset((
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
))
WATER_NAMES = frozenset(('HOH', 'DOD', 'WAT', 'H2O'))

STRUCTURE_HEADER_CATEGORIES = (
    'entity', 'entity_poly', 'pdbx_poly_seq_scheme', 'pdbx_nonpoly_scheme', 'pdbx_branch_scheme',
    'pdbx_chain_remapping', 'chem_comp',
)
_POLY_SCHEME_COLUMNS = ('asym_id', 'entity_id', 'seq_id', 'mon_id', 'pdb_mon_id', 'pdb_strand_id', 'pdb_seq_num', 'pdb_ins_code')
_NONPOLY_SCHEME_COLUMNS = ('asym_id', 'entity_id', 'mon_id', 'pdb_seq_num', 'pdb_ins_code', 'pdb_strand_id')
# Branched entities (glycans): pdb_asym_id is the chain _atom_site gives them.
_BRANCH_SCHEME_COLUMNS = ('asym_id', 'pdb_asym_id')
_REMAP_COLUMNS = ('label_asym_id', 'auth_asym_id', 'orig_label_asym_id', 'orig_auth_asym_id')

# Derived connectivity (the SSBOND / LINK / CISPEP records of the PDB format files).
SSBOND_MAX_DISTANCE = 2.5       # SG-SG, as the RIN engine's disulfide
# Longest covalent bond between two non-metal atoms of different residues (C-S ~1.8 A); closer
# contacts than hydrogen bonds but above this are left out.
COVALENT_LINK_MAX_DISTANCE = 1.95
METAL_LINK_MAX_DISTANCE = 2.8   # metal coordination, as the RIN engine's
# Atoms closer than this are a clash, not a bond (the clash rule of the distance features).
LINK_MIN_DISTANCE = 1.2
CIS_MAX_OMEGA = 30.0


def residue_id(number: int, ins_code: str) -> int | str:
    """The residue id of the PDB format pipeline: 52, or '52A' with an insertion code."""
    if ins_code in ('', '?', '.'):
        return int(number)
    return f'{int(number)}{ins_code}'


def residue_id_from_text(text: str) -> int | str:
    """'52' -> 52, '52A' -> '52A' (mmcifParser's residue ids)."""
    try:
        return int(text)
    except ValueError:
        return text


def chain_type_of_entity(poly_type: str | None, n_observed: int, n_irregular: int) -> str:
    if poly_type is None:
        return 'MISC'
    if poly_type == 'polypeptide(L)':
        if n_observed < PEPTIDE_MAX_OBSERVED_RESIDUES or (n_observed - n_irregular) <= 2 * n_irregular:
            return 'Peptide'
        return 'Protein'
    if poly_type == 'polypeptide(D)':
        return 'Peptide'
    if poly_type == 'polyribonucleotide':
        return 'RNA'
    if 'deoxyribonucleotide' in poly_type:
        return 'DNA'
    return 'MISC'


# ---------------------------------------------------------------------------------------------
# Header

@dataclass
class Instance:
    """A polymer instance of the file (one label_asym_id)."""
    label: str
    chain: str                      # author chain id, full length
    entity_id: str
    poly_type: str | None
    chain_type: str
    kept: bool
    # The author chain of the deposited entry it was made from (pdbx_chain_remapping), for names.
    orig_chain: str
    # For a symmetry copy: the chain of the kept instance it copies.
    original_chain: str | None = None
    residue_ids: list = field(default_factory = list)    # observed, in seq_id order
    residue_names: list = field(default_factory = list)


@dataclass
class StructureHeader:
    instances: dict[str, Instance]              # by label_asym_id
    chain_instance: dict[str, Instance]         # by author chain
    label_to_chain: dict[str, str]              # every label_asym_id (also non-polymers) -> author chain
    entity_types: dict[str, str]                # entity id -> _entity.type
    kept_labels: set[str]                       # all kept instances, polymers and non-polymers
    # Kept protein chains: chain -> (one letter sequence, residue ids in contig order), exactly as
    # in the search DB.
    kept_protein_residues: dict[str, tuple[str, list]]
    chain_origins: dict[str, str]               # copy chain -> kept chain it copies
    # (label, residue name, residue id, chain) of the kept non-polymer residues.
    nonpolymers: list[tuple[str, str, int | str, str]]
    sugar_names: set[str]
    observed_residues: int

    def chain_types(self) -> dict[str, str]:
        return {instance.chain: instance.chain_type for instance in self.instances.values()}

    def homomer_groups(self) -> list[list[str]]:
        """Kept Protein chains grouped by entity."""
        groups: dict[str, list[str]] = {}
        for instance in self.instances.values():
            if instance.kept and instance.chain_type == 'Protein':
                groups.setdefault(instance.entity_id, []).append(instance.chain)
        return list(groups.values())


def _category_rows(extract: mmcifParser.CategoryExtract, name: str, columns: tuple[str, ...]) -> list[tuple[str, ...]]:
    category_lines = extract.categories.get(name)
    if category_lines is None:
        return []
    try:
        return list(mmcifParser.iter_simple_loop_rows(category_lines, columns))
    except mmcifParser._NotASimpleLoop:
        return list(mmcifParser._biotite_rows(name, category_lines, columns))
    except mmcifParser.MissingCategoryError:
        return []


def _kept_from_remap(remap_rows, labels) -> tuple[set[str], dict[str, str]]:
    """(kept labels, original label -> kept label): the rule of rin_engine.structure.kept_label_asyms.

    Without remapping every instance is kept."""
    if not remap_rows:
        return set(labels), {label: label for label in labels}
    chosen: dict[str, tuple[str, bool]] = {}
    for label, auth, orig_label, orig_auth in remap_rows:
        is_original = auth == orig_auth
        current = chosen.get(orig_label)
        if current is None or (is_original and not current[1]):
            chosen[orig_label] = (label, is_original)
    kept = {label for label, _ in chosen.values()}
    covered = {row[0] for row in remap_rows}
    kept.update(label for label in labels if label not in covered)
    return kept, {orig_label: label for orig_label, (label, _) in chosen.items()}


def read_structure_header(lines) -> StructureHeader:
    extract = mmcifParser.extract_cif_categories(lines, STRUCTURE_HEADER_CATEGORIES, stop_category = 'atom_site')

    entity_types = {entity_id: entity_type for entity_id, entity_type in _category_rows(extract, 'entity', ('id', 'type'))}
    poly_types = {entity_id: poly_type for entity_id, poly_type in _category_rows(extract, 'entity_poly', ('entity_id', 'type'))}
    remap_rows = _category_rows(extract, 'pdbx_chain_remapping', _REMAP_COLUMNS)
    remap = {row[0]: row for row in remap_rows}

    # Every polymer instance with its observed residues, first observed row per seq_id.
    label_to_chain: dict[str, str] = {}
    instance_entity: dict[str, str] = {}
    residues_of: dict[str, dict[int, tuple[str, int | str]]] = {}
    observed_residues = 0
    for asym_id, entity_id, seq_id, mon_id, pdb_mon_id, strand_id, seq_num, ins_code in _category_rows(
            extract, 'pdbx_poly_seq_scheme', _POLY_SCHEME_COLUMNS):
        if asym_id not in label_to_chain:
            label_to_chain[asym_id] = strand_id
            instance_entity[asym_id] = entity_id
            residues_of[asym_id] = {}
        if pdb_mon_id in mmcifParser.UNOBSERVED:
            continue
        number = int(seq_id)
        chain_residues = residues_of[asym_id]
        if number not in chain_residues:
            chain_residues[number] = (mon_id, residue_id(int(seq_num), ins_code))
            observed_residues += 1

    nonpolymer_rows = _category_rows(extract, 'pdbx_nonpoly_scheme', _NONPOLY_SCHEME_COLUMNS)
    nonpolymer_labels = []
    for asym_id, entity_id, mon_id, seq_num, ins_code, strand_id in nonpolymer_rows:
        if asym_id not in label_to_chain:
            label_to_chain[asym_id] = strand_id
            nonpolymer_labels.append(asym_id)

    for asym_id, pdb_asym_id in _category_rows(extract, 'pdbx_branch_scheme', _BRANCH_SCHEME_COLUMNS):
        label_to_chain.setdefault(asym_id, pdb_asym_id)

    kept_labels, kept_label_of_orig = _kept_from_remap(remap_rows, list(label_to_chain))

    # Kept protein chains exactly as the search DB has them: select_protein_chains decides which
    # instance of a deposited chain is kept and which residues it has.
    protein_chains, _ = mmcifParser.get_observed_protein_residues(extract)
    kept_protein_residues: dict[str, tuple[str, list]] = {}
    for chain, (seq, res_pos_map) in protein_chains.items():
        ordered = sorted(res_pos_map.items(), key = lambda item: item[1])
        kept_protein_residues[chain] = (seq, [residue_id_from_text(text) for text, _ in ordered])
    kept_protein_chains = set(kept_protein_residues)

    instances: dict[str, Instance] = {}
    for label, chain in label_to_chain.items():
        if label not in residues_of:
            continue
        entity_id = instance_entity[label]
        poly_type = poly_types.get(entity_id)
        ordered = [residues_of[label][number] for number in sorted(residues_of[label])]
        names = [name for name, _ in ordered]
        n_irregular = sum(1 for name in names if name not in STANDARD_AMINO_ACIDS)
        if poly_type == 'polypeptide(L)':
            # The protein instances of select_protein_chains are the kept ones.
            kept = chain in kept_protein_chains
        else:
            kept = label in kept_labels
        remapped = remap.get(label)
        instances[label] = Instance(
            label = label, chain = chain, entity_id = entity_id, poly_type = poly_type,
            chain_type = chain_type_of_entity(poly_type, len(names), n_irregular),
            kept = kept, orig_chain = remapped[3] if remapped is not None else chain,
            residue_ids = [rid for _, rid in ordered], residue_names = names,
        )

    # Copies point to the kept instance of their original chain.
    kept_label_by_orig_label: dict[str, str] = dict(kept_label_of_orig)
    for instance in instances.values():
        if instance.kept and instance.label in remap:
            kept_label_by_orig_label[remap[instance.label][2]] = instance.label
    chain_origins: dict[str, str] = {}
    for instance in instances.values():
        if instance.kept or instance.label not in remap:
            continue
        original_label = kept_label_by_orig_label.get(remap[instance.label][2])
        if original_label is not None and original_label in instances and instances[original_label].kept:
            instance.original_chain = instances[original_label].chain
            chain_origins[instance.chain] = instance.original_chain

    nonpolymers = []
    for asym_id, entity_id, mon_id, seq_num, ins_code, strand_id in nonpolymer_rows:
        if asym_id not in kept_labels or mon_id in WATER_NAMES or entity_types.get(entity_id) == 'water':
            continue
        nonpolymers.append((asym_id, mon_id, residue_id(int(seq_num), ins_code), strand_id))

    sugar_names = {comp_id for comp_id, comp_type in _category_rows(extract, 'chem_comp', ('id', 'type'))
                   if 'saccharide' in comp_type.lower()}

    # Polymer instances are kept or not as decided above (select_protein_chains for proteins).
    kept_labels = ({label for label in kept_labels if label not in instances}
                   | {instance.label for instance in instances.values() if instance.kept})

    return StructureHeader(
        instances = instances,
        chain_instance = {instance.chain: instance for instance in instances.values()},
        label_to_chain = label_to_chain,
        entity_types = entity_types,
        kept_labels = kept_labels,
        kept_protein_residues = kept_protein_residues,
        chain_origins = chain_origins,
        nonpolymers = nonpolymers,
        sugar_names = sugar_names,
        observed_residues = observed_residues,
    )


# ---------------------------------------------------------------------------------------------
# What the complex parsing needs (serializedPipeline.parse_complexes_body)

@dataclass
class ComplexInfo:
    structure_id: str
    resolution: float | None
    chain_types: dict[str, str]                 # every polymer instance, copies included
    chain_names: dict[str, str]                 # kept chains -> UniProt name, where known
    homomer_groups: list[list[str]]
    chain_origins: dict[str, str]
    # Kept Protein chains -> (residue ids in contig order, one letter sequence): Parse_Data.
    protein_chains: dict[str, tuple[list, str]]
    interaction_partners: list[list]
    atom_count: int


def complex_info(structure_id: str, config) -> ComplexInfo | None:
    from structman.lib.sdsc.sdsc_utils import boring

    source = structure_source.mmcif_source(structure_id, config)
    if source is None:
        return None
    with source.open_lines() as lines:
        header = read_structure_header(lines)

    if source.kind == structure_source.KIND_PDB_ASSEMBLY:
        pdb_id, _ = mmcifParser.split_assembly_id(structure_id)
        entry = structure_source.entry_header(pdb_id, config)
        resolution = entry.resolution
        entry_names = entry.chain_names
    else:
        resolution = structure_source.MODEL_RESOLUTION
        entry_names = {}

    chain_types = header.chain_types()
    chain_names = {}
    for instance in header.instances.values():
        if instance.kept and instance.orig_chain in entry_names:
            chain_names[instance.chain] = entry_names[instance.orig_chain]

    protein_chains = {}
    for chain, (seq, residue_ids) in header.kept_protein_residues.items():
        if chain_types.get(chain) == 'Protein':
            protein_chains[chain] = (residue_ids, seq)

    interaction_partners: list[list] = []
    for label, name, rid, chain in header.nonpolymers:
        if chain_types.get(chain) == 'Peptide' or boring(name):
            continue
        interaction_partners.append(['Ligand', name, rid, chain])
    for instance in header.instances.values():
        if instance.kept:
            interaction_partners.append([instance.chain_type, instance.chain])

    # Only used to cost the annotation of the complex: about 8 heavy atoms per residue.
    atom_count = 8 * header.observed_residues + 10 * len(header.nonpolymers)

    return ComplexInfo(
        structure_id = structure_id, resolution = resolution, chain_types = chain_types,
        chain_names = chain_names, homomer_groups = header.homomer_groups(),
        chain_origins = header.chain_origins, protein_chains = protein_chains,
        interaction_partners = interaction_partners, atom_count = atom_count,
    )


# ---------------------------------------------------------------------------------------------
# The parsed structure

class Parsed_structure(Slotted_obj):
    """Chain-number indexed columns of a structure.

    chain_list holds the structure's chains; a chain number is a position in it. Without a chain
    list it is the 62 slot alphabet the PDB format parser (templateFiltering.parsePDB) uses.
    """

    __slots__ = [
        'chain_numbers', 'residue_names', 'residue_atoms', 'hetero_names',
        'hetero_atoms', 'box', 'chain_type', 'b_factors', 'res_contig', 'modres',
        'ssbond', 'link', 'cis_conformation', 'cis_follower', 'ligands', 'metals',
        'ions',
        # mmCIF only
        'chain_list', 'chain_index', 'kept', 'original_chain', 'entity_id', 'label_chain_number',
        'residue_ids', 'residue_rows', 'sugars', 'protein_atoms', 'header',
    ]

    def __init__(self, chain_list: list[str] | None = None):
        if chain_list is None:
            chain_list = list(CHAIN_ID_ORDER)
        self.chain_list = list(chain_list)
        self.chain_index = {chain: number for number, chain in enumerate(self.chain_list)}
        slots = len(self.chain_list)
        self.chain_numbers = np.empty(0, dtype=np.int64)
        self.residue_names = [None] * slots
        self.residue_atoms = [None] * slots
        self.hetero_names = [None] * slots
        self.hetero_atoms = [None] * slots
        self.box = np.full((slots, 6), np.nan, dtype=np.float32)
        self.chain_type = np.full(slots, -1, dtype=np.int8)
        self.b_factors = [None] * slots
        self.res_contig = [None] * slots
        self.modres = [None] * slots
        self.ssbond = [None] * slots
        self.link = [None] * slots
        self.cis_conformation = [None] * slots
        self.cis_follower = [None] * slots
        self.ligands = set()
        self.metals = set()
        self.ions = set()
        self.kept = np.zeros(slots, dtype=bool)
        self.original_chain = [None] * slots
        self.entity_id = [None] * slots
        self.label_chain_number = {}
        self.residue_ids = [None] * slots
        self.residue_rows = [None] * slots
        self.sugars = set()
        self.protein_atoms = None
        self.header = None

    @property
    def n_chains(self) -> int:
        return len(self.chain_list)

    def chain_id(self, chain_number: int) -> str:
        return self.chain_list[chain_number]

    def chain_number(self, chain_id: str) -> int:
        return self.chain_index.get(chain_id, -1)

    def chain_id_number_of(self, chain_id: str) -> int:
        return self.chain_number(chain_id)

    def chain_ids(self) -> list[str]:
        return [self.chain_list[number] for number in self.chain_numbers.tolist()]

    def is_protein(self, chain_number: int) -> bool:
        return self.chain_type[chain_number] == CHAIN_PROTEIN

    def type_name(self, chain_number: int) -> str | None:
        code = self.chain_type[chain_number]
        if code == -1:
            return None
        return CHAIN_TYPE_NAMES[code]

    def chain_type_map(self) -> dict[str, str]:
        return {self.chain_list[number]: CHAIN_TYPE_NAMES[self.chain_type[number]]
                for number in range(self.n_chains)
                if self.chain_type[number] != -1}

    def protein_chain_numbers(self) -> list[int]:
        return [number for number in range(self.n_chains) if self.chain_type[number] == CHAIN_PROTEIN]

    def res_contig_map(self) -> dict:
        return {self.chain_list[number]: self.res_contig[number]
                for number in self.chain_numbers.tolist()
                if self.res_contig[number] is not None}

    def chain_keyed(self, values: set) -> set:
        """A {(chain number, residue)} set back as {(chain id, residue)}."""
        return {(self.chain_list[number], residue) for number, residue in values}

    def analysed_chain_numbers(self) -> list[int]:
        """The chains the structural analysis is about: kept protein chains with residues."""
        return [number for number in self.chain_numbers.tolist()
                if self.kept[number] and self.is_protein(number) and self.res_contig[number] is not None]


def _residue_matrix(rows: list[np.ndarray]) -> np.ndarray:
    """(n residues, max atoms, 4) float32 of (atom serial, x, y, z), NaN padded."""
    max_atoms = max((len(row) for row in rows), default = 0)
    matrix = np.full((len(rows), max(max_atoms, 1), 4), np.nan, dtype = np.float32)
    for index, row in enumerate(rows):
        matrix[index, :len(row)] = row
    return matrix


def _read_atoms(source: structure_source.MmcifSource, region_radius: float):
    """The atoms of the first model (first alternative location, no waters or hydrogens)."""
    import biotite.structure.io.pdbx as pdbx
    from structman.lib.rin_engine import structure as rin_structure

    block = None
    if source.path is not None and os.path.getsize(source.path) >= rin_structure.STREAMING_MIN_BYTES:
        block, _ = rin_structure._streamed_region_block(source.path, radius = region_radius)
    if block is None:
        with source.open_lines() as f:
            block = pdbx.CIFFile.read(f).block
    atoms = pdbx.get_structure(
        block, model = 1, altloc = 'first',
        extra_fields = ['atom_id', 'label_asym_id', 'label_entity_id', 'label_seq_id', 'b_factor'])
    return atoms


def _region_mask(coords: np.ndarray, residue_index: np.ndarray, kept_atoms: np.ndarray, radius: float) -> np.ndarray:
    """Whole residues with an atom within `radius` of a kept atom."""
    if kept_atoms.all() or not kept_atoms.any():
        return np.ones(len(coords), dtype = bool)
    from scipy.spatial import cKDTree
    tree = cKDTree(coords[kept_atoms])
    distance, _ = tree.query(coords, k = 1, distance_upper_bound = radius)
    near = np.isfinite(distance) | kept_atoms
    return np.isin(residue_index, np.unique(residue_index[near]))


def _residue_starts(*columns) -> np.ndarray:
    n = len(columns[0])
    if n == 0:
        return np.zeros(0, dtype = np.int64)
    change = np.zeros(n, dtype = bool)
    change[0] = True
    for column in columns:
        change[1:] |= column[1:] != column[:-1]
    return np.flatnonzero(change)


def parse_structure(source: structure_source.MmcifSource, config, region_radius: float | None = None) -> Parsed_structure:
    """The parsed structure: kept instances plus everything within `region_radius` of them."""
    from structman.lib.sdsc.consts.residues import METAL_ATOMS, ION_ATOMS

    if region_radius is None:
        region_radius = max(float(config.ligand_interest_sphere), float(config.intra_milieu_threshold))

    with source.open_lines() as lines:
        header = read_structure_header(lines)
    atoms = _read_atoms(source, region_radius)

    label = np.asarray(atoms.label_asym_id).astype(str)
    entity = np.asarray(atoms.label_entity_id).astype(str)
    seq_id = np.asarray(atoms.label_seq_id).astype(str)
    res_name = np.asarray(atoms.res_name).astype(str)
    element = np.asarray(atoms.element).astype(str)
    water_entities = [entity_id for entity_id, entity_type in header.entity_types.items() if entity_type == 'water']
    keep = ~(np.isin(res_name, list(WATER_NAMES)) | np.isin(entity, water_entities) | np.isin(element, ('H', 'D')))
    atoms = atoms[keep]
    label = label[keep]
    entity = entity[keep]
    seq_id = seq_id[keep]
    res_name = res_name[keep]
    element = element[keep]

    res_ids = np.asarray(atoms.res_id)
    ins_codes = np.asarray(atoms.ins_code).astype(str)
    coords = np.asarray(atoms.coord, dtype = np.float32)
    starts = _residue_starts(label, res_ids, ins_codes, res_name)
    residue_index = np.zeros(len(label), dtype = np.int64)
    if len(starts) > 1:
        residue_index[starts[1:]] = 1
    residue_index = np.cumsum(residue_index)

    # Kept instances plus their surroundings.
    kept_atoms = np.isin(label, list(header.kept_labels)) if header.kept_labels else np.ones(len(label), dtype = bool)
    region = _region_mask(coords, residue_index, kept_atoms, region_radius)
    if not region.all():
        atoms = atoms[region]
        label, entity, res_name, element = label[region], entity[region], res_name[region], element[region]
        seq_id = seq_id[region]
        res_ids, ins_codes, coords = res_ids[region], ins_codes[region], coords[region]
        starts = _residue_starts(label, res_ids, ins_codes, res_name)

    atom_serial = np.asarray(atoms.atom_id, dtype = np.float32)
    b_factor = np.asarray(atoms.b_factor, dtype = np.float32)
    atom_name = np.asarray(atoms.atom_name).astype(str)
    # Full author chain ids: biotite keeps at most four characters of them.
    chain_of_label = dict(header.label_to_chain)
    biotite_chain = np.asarray(atoms.chain_id).astype(str)

    ends = np.append(starts[1:], len(label))
    residues = []       # (chain, label, residue id, name, start, end, is_polymer)
    chain_order: list[str] = []
    for start, end in zip(starts.tolist(), ends.tolist()):
        res_label = label[start]
        chain = chain_of_label.get(res_label, biotite_chain[start])
        instance = header.instances.get(res_label)
        # By label_seq_id: in converted models a ligand shares the label_asym_id of its chain.
        is_polymer = instance is not None and instance.poly_type is not None and seq_id[start] not in ('.', '?')
        residues.append((chain, res_label, residue_id(res_ids[start], ins_codes[start]), res_name[start], start, end, is_polymer))
        if chain not in chain_order:
            chain_order.append(chain)

    parsed = Parsed_structure(chain_order)
    parsed.header = header
    parsed.chain_numbers = np.arange(len(chain_order), dtype = np.int64)
    parsed.label_chain_number = {res_label: parsed.chain_number(chain_of_label.get(res_label, res_label))
                                 for res_label in np.unique(label).tolist()}

    def quads(start, end):
        block = np.empty((end - start, 4), dtype = np.float32)
        block[:, 0] = atom_serial[start:end]
        block[:, 1:4] = coords[start:end]
        return block

    polymer_of_chain: dict[int, list] = {number: [] for number in range(len(chain_order))}
    hetero_of_chain: dict[int, list] = {number: [] for number in range(len(chain_order))}
    for residue in residues:
        number = parsed.chain_number(residue[0])
        (polymer_of_chain if residue[6] else hetero_of_chain)[number].append(residue)

    protein_atom_ranges = []    # (chain number, row, start, end) of protein residues, for DSSP
    for number, chain in enumerate(chain_order):
        instance = header.chain_instance.get(chain)
        polymer_residues = polymer_of_chain[number]
        if instance is not None:
            parsed.kept[number] = instance.kept
            parsed.original_chain[number] = instance.original_chain
            parsed.entity_id[number] = instance.entity_id
            parsed.chain_type[number] = CHAIN_TYPE_CODES[instance.chain_type]
        elif hetero_of_chain[number]:
            parsed.chain_type[number] = CHAIN_MISC

        if polymer_residues:
            by_id = {}
            for residue in polymer_residues:
                by_id.setdefault(residue[2], residue)
            if instance is not None and instance.kept and chain in header.kept_protein_residues:
                # The search DB's residues, in their order.
                ordered_ids = header.kept_protein_residues[chain][1]
            else:
                ordered_ids = list(by_id)
            name_of = dict(zip(instance.residue_ids, instance.residue_names)) if instance is not None else {}
            rows, names, b_rows = [], [], []
            for rid in ordered_ids:
                residue = by_id.get(rid)
                if residue is None:
                    rows.append(np.zeros((0, 4), dtype = np.float32))
                    names.append(name_of.get(rid, 'UNK'))
                    b_rows.append(np.zeros(0, dtype = np.float32))
                    continue
                _, _, _, name, start, end, _ = residue
                rows.append(quads(start, end))
                names.append(name)
                b_rows.append(b_factor[start:end])
                if parsed.chain_type[number] in (CHAIN_PROTEIN, CHAIN_PEPTIDE):
                    protein_atom_ranges.append((number, len(rows) - 1, start, end))

            matrix = _residue_matrix(rows)
            parsed.residue_atoms[number] = matrix
            parsed.residue_names[number] = np.array(names, dtype = 'S5')
            parsed.residue_ids[number] = list(ordered_ids)
            parsed.residue_rows[number] = {rid: row for row, rid in enumerate(ordered_ids)}
            contig = Residue_Map()
            for row, (rid, name) in enumerate(zip(ordered_ids, names)):
                contig.add_item(rid, (row + 1, name))
            parsed.res_contig[number] = contig
            if parsed.chain_type[number] == CHAIN_PROTEIN:
                b_matrix = np.full(matrix.shape[:2], np.nan, dtype = np.float32)
                for row, values in enumerate(b_rows):
                    b_matrix[row, :len(values)] = values
                parsed.b_factors[number] = b_matrix
            modres = {row: name for row, name in enumerate(names)
                      if parsed.chain_type[number] == CHAIN_PROTEIN and name not in STANDARD_AMINO_ACIDS and name != 'UNK'}
            parsed.modres[number] = modres
        else:
            parsed.residue_atoms[number] = np.full((0, 1, 4), np.nan, dtype = np.float32)
            parsed.residue_names[number] = np.zeros(0, dtype = 'S5')
            parsed.residue_ids[number] = []
            parsed.residue_rows[number] = {}

        hetero = hetero_of_chain[number]
        if hetero:
            parsed.hetero_atoms[number] = _residue_matrix([quads(start, end) for _, _, _, _, start, end, _ in hetero])
            parsed.hetero_names[number] = np.array([name for _, _, _, name, _, _, _ in hetero], dtype = 'S5')
            for _, het_label, rid, name, start, end, _ in hetero:
                key = (number, rid)
                if name in METAL_ATOMS:
                    parsed.metals.add(key)
                elif name in ION_ATOMS:
                    parsed.ions.add(key)
                else:
                    parsed.ligands.add(key)
                if name in header.sugar_names or header.entity_types.get(entity[start]) == 'branched':
                    parsed.sugars.add(key)
        else:
            parsed.hetero_atoms[number] = np.full((0, 1, 4), np.nan, dtype = np.float32)
            parsed.hetero_names[number] = np.zeros(0, dtype = 'S5')

        chain_atoms = [parsed.residue_atoms[number].reshape(-1, 4), parsed.hetero_atoms[number].reshape(-1, 4)]
        chain_coords = np.concatenate(chain_atoms)[:, 1:4]
        chain_coords = chain_coords[~np.isnan(chain_coords[:, 0])]
        if len(chain_coords):
            # (x, z, y), the axis order the PDB format parser used for its boxes.
            parsed.box[number] = (chain_coords[:, 0].min(), chain_coords[:, 0].max(),
                                  chain_coords[:, 2].min(), chain_coords[:, 2].max(),
                                  chain_coords[:, 1].min(), chain_coords[:, 1].max())

    _derive_connectivity(parsed, residues, atom_name, element, coords, res_name)
    parsed.protein_atoms = (protein_atom_ranges, atom_name, element, res_name, coords, b_factor, atom_serial)
    return parsed


# ---------------------------------------------------------------------------------------------
# SSBOND / LINK / CISPEP from the coordinates

def _derive_connectivity(parsed: Parsed_structure, residues, atom_name, element, coords, res_name) -> None:
    from scipy.spatial import cKDTree
    from structman.lib.rin_engine.constants import METAL_ELEMENTS

    n_atoms = len(coords)
    for number in range(parsed.n_chains):
        parsed.ssbond[number] = {}
        parsed.link[number] = {}
        parsed.cis_conformation[number] = {}
        parsed.cis_follower[number] = {}
    if n_atoms == 0:
        return

    # Per atom: its residue (chain number, row or -1 for hetero residues, residue id, name).
    atom_chain = np.full(n_atoms, -1, dtype = np.int64)
    atom_row = np.full(n_atoms, -1, dtype = np.int64)
    atom_residue = np.full(n_atoms, -1, dtype = np.int64)
    residue_info = []
    for index, (chain, _, rid, name, start, end, is_polymer) in enumerate(residues):
        number = parsed.chain_number(chain)
        row = parsed.residue_rows[number].get(rid, -1) if is_polymer else -1
        atom_chain[start:end] = number
        atom_row[start:end] = row
        atom_residue[start:end] = index
        residue_info.append((number, row, rid, name))

    tree = cKDTree(coords)
    pairs = tree.query_pairs(METAL_LINK_MAX_DISTANCE, output_type = 'ndarray')
    if len(pairs):
        a, b = pairs[:, 0], pairs[:, 1]
        different = atom_residue[a] != atom_residue[b]
        a, b = a[different], b[different]
        distance = np.linalg.norm(coords[a] - coords[b], axis = 1)
        is_sg = (atom_name == 'SG') & np.isin(res_name, ('CYS',))
        ssbond = is_sg[a] & is_sg[b] & (distance <= SSBOND_MAX_DISTANCE)
        metal = np.isin(element, list(METAL_ELEMENTS))
        # A link connects a polymer residue to a non-polymer one (ligand, glycan, metal), as the
        # LINK records did: between two polymer residues a close pair is a backbone bond or, in
        # low resolution structures, a clash between symmetry copies.
        hetero_pair = (atom_row[a] < 0) | (atom_row[b] < 0)
        link = ~ssbond & hetero_pair & (distance >= LINK_MIN_DISTANCE) & (
            ((metal[a] | metal[b]) & (distance <= METAL_LINK_MAX_DISTANCE))
            | (~metal[a] & ~metal[b] & (distance <= COVALENT_LINK_MAX_DISTANCE)))

        for i, j, d in zip(a[ssbond].tolist(), b[ssbond].tolist(), distance[ssbond].tolist()):
            for x, y in ((i, j), (j, i)):
                number, row, _, _ = residue_info[atom_residue[x]]
                if row < 0:
                    continue
                other_number, _, other_rid, _ = residue_info[atom_residue[y]]
                parsed.ssbond[number][row] = (parsed.chain_id(other_number), other_rid, round(d, 2))

        best_link: dict[tuple[int, int], tuple] = {}
        for i, j, d in zip(a[link].tolist(), b[link].tolist(), distance[link].tolist()):
            for x, y in ((i, j), (j, i)):
                number, row, _, name = residue_info[atom_residue[x]]
                if row < 0:
                    continue
                other_number, _, other_rid, other_name = residue_info[atom_residue[y]]
                key = (number, row)
                if key not in best_link or d < best_link[key][6]:
                    best_link[key] = (str(atom_name[x]), str(name), str(atom_name[y]), str(other_name), other_rid,
                                      parsed.chain_id(other_number), round(d, 2))
        for (number, row), value in best_link.items():
            parsed.link[number][row] = value

    # Cis peptides: omega (CA_i, C_i, N_i+1, CA_i+1) within 30 degrees of 0, peptide bond present.
    backbone_of_chain: dict[int, dict[int, dict[str, np.ndarray]]] = {}
    for index, (chain_number, row, rid, name) in enumerate(residue_info):
        if row < 0 or parsed.chain_type[chain_number] not in (CHAIN_PROTEIN, CHAIN_PEPTIDE):
            continue
        start, end = residues[index][4], residues[index][5]
        backbone_of_chain.setdefault(chain_number, {})[row] = {
            atom_name[k]: coords[k] for k in range(start, end) if atom_name[k] in ('N', 'CA', 'C')}
    for number, backbone in backbone_of_chain.items():
        names = parsed.residue_names[number]
        rids = parsed.residue_ids[number]
        for row in sorted(backbone):
            nxt = backbone.get(row + 1)
            cur = backbone[row]
            if nxt is None or not all(atom in cur for atom in ('CA', 'C')) or not all(atom in nxt for atom in ('N', 'CA')):
                continue
            if np.linalg.norm(cur['C'] - nxt['N']) > 2.0:
                continue
            omega = _dihedral(cur['CA'], cur['C'], nxt['N'], nxt['CA'])
            if abs(omega) < CIS_MAX_OMEGA:
                chain = parsed.chain_id(number)
                parsed.cis_conformation[number][row] = (names[row + 1].decode(), chain, rids[row + 1], round(omega, 1))
                parsed.cis_follower[number][row + 1] = (names[row].decode(), chain, rids[row], round(omega, 1))


def _dihedral(p0, p1, p2, p3) -> float:
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


# ---------------------------------------------------------------------------------------------
# DSSP

# The atom_site columns the bundled DSSP (xssp 3.0.7, structure.cpp ReadmmCIF) reads.
_DSSP_ATOM_SITE_COLUMNS = (
    'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id',
    'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
    'B_iso_or_equiv', 'pdbx_formal_charge', 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id',
    'pdbx_PDB_model_num',
)


def _chain_codes(n: int) -> list[str]:
    """n distinct chain codes of at most 4 characters: single characters first."""
    alphabet = CHAIN_ID_ORDER
    codes = list(alphabet[:n])
    length = 2
    while len(codes) < n:
        from itertools import product
        for combination in product(alphabet, repeat = length):
            codes.append(''.join(combination))
            if len(codes) == n:
                break
        length += 1
    return codes


def _cif_token(value: str) -> str:
    if value == '' or any(character in value for character in ' \t'):
        return f'"{value}"'
    if "'" in value or '"' in value:
        return f'"{value}"' if '"' not in value else f"'{value}'"
    return value


def write_dssp_input(parsed: Parsed_structure, path: str) -> dict[str, int]:
    """The protein residues of the parsed structure as a minimal mmCIF that DSSP reads.

    Every chain gets a code of at most four characters (DSSP prints longer chain ids truncated) as
    its label and author chain id; author residue numbers and insertion codes stay. Returns
    {code: chain number}.
    """
    ranges, atom_name, element, res_name, coords, b_factor, atom_serial = parsed.protein_atoms
    chain_numbers = sorted({number for number, _, _, _ in ranges})
    codes = dict(zip(chain_numbers, _chain_codes(len(chain_numbers))))

    # No _entry.id: with one, DSSP 3.0.7 (smssp) looks for the entry's other header categories and
    # crashes (SIGSEGV, no output) when they are missing - every structure lost its DSSP that way.
    lines = ['data_dssp\n', 'loop_\n']
    lines.extend(f'_atom_site.{column}\n' for column in _DSSP_ATOM_SITE_COLUMNS)
    serial = 0
    for number, row, start, end in ranges:
        code = codes[number]
        rid = parsed.residue_ids[number][row]
        if isinstance(rid, int):
            seq_num, ins_code = rid, '?'
        else:
            digits = rid.rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz')
            seq_num, ins_code = int(digits), rid[len(digits):]
        name = _cif_token(res_name[start])
        for k in range(start, end):
            serial += 1
            x, y, z = coords[k]
            atom = _cif_token(atom_name[k])
            lines.append(
                f'ATOM {serial} {element[k]} {atom} . {name} {code} 1 {row + 1} {ins_code} '
                f'{x:.3f} {y:.3f} {z:.3f} 1.00 {b_factor[k]:.2f} ? {seq_num} {name} {code} {atom} 1\n')
    lines.append('#\n')
    with open(path, 'w') as f:
        f.writelines(lines)
    return {code: number for number, code in codes.items()}


def parse_dssp_output(out: str) -> tuple[dict[str, dict[int | str, tuple]], list[str]]:
    """{chain code: {residue id (52, '52A'): (rsa, rel mc acc, rel sc acc, ss, phi, psi, tco, kappa, alpha)}}.

    The patched DSSP format (dssp_n.patch: ACC MCA SCA), as templateFiltering.calcDSSP read it; the
    chain is the classic chain column, or for longer codes the label chain column at the end.
    """
    from structman.lib.sdsc.consts.residues import ONE_TO_THREE, RESIDUE_MAX_ACC

    results: dict[str, dict[int | str, tuple]] = {}
    errors: list[str] = []
    lines = out.split('\n')
    start = None
    for index, line in enumerate(lines):
        if line.startswith('  #  RESIDUE'):
            start = index + 1
            break
    if start is None:
        return results, ['No residue section in the DSSP output']

    for line in lines[start:]:
        if len(line) < 46 or line[13] == '!':
            continue
        if line[9] == ' ':
            continue
        try:
            rid = residue_id(int(line[5:10]), line[10].strip())
        except ValueError:
            errors.append(f'DSSP residue number not readable: {line}')
            continue
        code = line[11]
        if code == '>':
            code = line.split()[-2]
        aa = line[13]
        ssa = line[16]
        try:
            acc = float(line[34:38])
            mc_acc = float(line[38:42])
            sc_acc = float(line[42:46])
        except ValueError:
            errors.append(f'DSSP accessibility not readable: {line}')
            continue
        try:
            if aa == 'X':
                max_acc = 220.0
            elif aa.islower():  # cysteines in disulfide bonds
                max_acc = RESIDUE_MAX_ACC['Sander'][ONE_TO_THREE['C']]
            else:
                max_acc = RESIDUE_MAX_ACC['Sander'][ONE_TO_THREE[aa]]
            rsa, rel_mc, rel_sc = acc / max_acc, mc_acc / max_acc, sc_acc / max_acc
        except (KeyError, ZeroDivisionError):
            rsa = rel_mc = rel_sc = None
        phi = psi = tco = kappa = alpha = None
        try:
            if len(line) > 123:
                tco = float(line[93:99])
                kappa = float(line[99:105])
                alpha = float(line[105:111])
                phi = float(line[111:117])
                psi = float(line[117:123])
        except ValueError:
            pass
        results.setdefault(code, {})[rid] = (rsa, rel_mc, rel_sc, ssa, phi, psi, tco, kappa, alpha)
    return results, errors


def dssp_binary(config) -> str | None:
    path = getattr(config, 'dssp_path', '') or ''
    if path:
        return path
    import shutil
    return shutil.which('smssp')


def run_dssp(parsed: Parsed_structure, config) -> tuple[dict[str, list], list[str]]:
    """DSSP rows per protein chain, in residue order ({chain: [row or None, ...]}), and errors."""
    binary = dssp_binary(config)
    if binary is None or parsed.protein_atoms is None or not parsed.protein_atoms[0]:
        return {}, ['No DSSP binary' if binary is None else 'No protein residues']

    tmp_folder = getattr(config, 'temp_folder', None) or getattr(config, 'tmp_folder', None) or '.'
    path = os.path.join(tmp_folder, f'dssp_{uuid.uuid4().hex}.cif')
    codes = write_dssp_input(parsed, path)
    try:
        process = subprocess.run([binary, path], capture_output = True, text = True)
        out, err = process.stdout, process.stderr
    except OSError as e:
        out, err = '', str(e)
    finally:
        if getattr(config, 'verbosity', 0) < 6:
            try:
                os.remove(path)
            except OSError:
                pass

    if not out.strip():
        return {}, [f'DSSP produced no output: {err.strip()}']

    by_code, errors = parse_dssp_output(out)
    dssp_rows: dict[str, list] = {}
    for code, values in by_code.items():
        number = codes.get(code)
        if number is None:
            continue
        rows = []
        for rid in parsed.residue_ids[number]:
            rows.append(values.get(rid))
        dssp_rows[parsed.chain_id(number)] = rows
    return dssp_rows, errors
