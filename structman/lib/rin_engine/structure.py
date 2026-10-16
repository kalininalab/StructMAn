"""Loading a biological assembly file and preparing its atoms for interaction detection.

Preparation, in order:

1. Model 1, first alternative location, bonds from the CCD and the peptide links (biotite).
2. Waters and any hydrogens or deuteriums of the file are removed. Hydrogens are all placed by
   hydride later, so structures with and without deposited hydrogens are treated alike.
3. Kept chain instances: one instance per chain of the deposited entry, the same choice as the
   search DB makes for protein chains (mmcifParser.select_protein_chains). Interactions are only
   computed where at least one partner belongs to a kept instance; between two symmetry copies
   they would only repeat what the kept instances already show.
4. The atoms are cut down to the kept instances plus every residue within `relevant_radius` of
   them. For most assemblies that is everything; for a 60-mer capsid it is a small part.
5. Formal charges: the CCD's, with hydride's pH estimate for the 20 canonical amino acids. The
   assembly files carry no charges of their own.
6. Hydrogens by hydride, on the remaining atoms.
"""

from __future__ import annotations

import functools
import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure as struc
import biotite.structure.info as info
import biotite.structure.io.pdbx as pdbx

from structman.lib import mmcifParser
from structman.lib.rin_engine import constants as C

_EXTRA_FIELDS = ['label_asym_id', 'label_entity_id', 'label_seq_id', 'b_factor']
_AMINO_ACID_NAMES = frozenset(info.amino_acid_names())
_NUCLEOTIDE_NAMES = frozenset(info.nucleotide_names())


@dataclass
class PreparedAssembly:
    # Atoms of the relevant region with hydrogens. Annotations besides biotite's defaults:
    # charge, label_asym_id, label_entity_id, label_seq_id, kept (bool), kind (int8).
    atoms: struc.AtomArray
    stats: dict = field(default_factory = dict)


def read_block(path: str) -> pdbx.CIFBlock:
    with mmcifParser.open_cif_text(path) as f:
        return pdbx.CIFFile.read(f).block


def _columns(block, category_name: str, columns: tuple[str, ...]) -> list[list[str]] | None:
    if category_name not in block:
        return None
    category = block[category_name]
    if any(column not in category for column in columns):
        return None
    return [mmcifParser.category_column(category, column) for column in columns]


def entity_kinds(block) -> dict[str, int]:
    """entity id -> one of C.KIND_PROTEIN, KIND_NUCLEOTIDE, KIND_SUGAR, KIND_OTHER, or -1 for
    non-polymers (their kind depends on the residue) and -2 for water."""
    kinds: dict[str, int] = {}
    entity = _columns(block, 'entity', ('id', 'type'))
    if entity is not None:
        for entity_id, entity_type in zip(*entity):
            if entity_type == 'water':
                kinds[entity_id] = -2
            elif entity_type == 'branched':
                kinds[entity_id] = C.KIND_SUGAR
            elif entity_type == 'polymer':
                kinds[entity_id] = C.KIND_OTHER
            else:
                kinds[entity_id] = -1
    entity_poly = _columns(block, 'entity_poly', ('entity_id', 'type'))
    if entity_poly is not None:
        for entity_id, poly_type in zip(*entity_poly):
            if poly_type.startswith('polypeptide'):
                kinds[entity_id] = C.KIND_PROTEIN
            elif 'nucleotide' in poly_type:
                kinds[entity_id] = C.KIND_NUCLEOTIDE
            else:
                kinds[entity_id] = C.KIND_OTHER
    return kinds


_REMAP_COLUMNS = ('label_asym_id', 'auth_asym_id', 'orig_label_asym_id', 'orig_auth_asym_id')


def kept_label_asyms(block, label_asyms: list[str]) -> set[str]:
    """One instance per chain of the deposited entry (label asym ids of the assembly).

    The instance that carries the original author chain name if there is one, otherwise the
    first one in file order. Without _pdbx_chain_remapping every instance is kept.
    """
    return _kept_from_remap(_columns(block, 'pdbx_chain_remapping', _REMAP_COLUMNS), label_asyms)


def _kept_from_remap(remap: list[list[str]] | None, label_asyms: list[str]) -> set[str]:
    if remap is None:
        return set(label_asyms)
    chosen: dict[str, tuple[str, bool]] = {}
    for label, auth, orig_label, orig_auth in zip(*remap):
        is_original = auth == orig_auth
        current = chosen.get(orig_label)
        if current is None or (is_original and not current[1]):
            chosen[orig_label] = (label, is_original)
    kept = {label for label, _ in chosen.values()}
    # Label asyms of the atoms that the remapping does not mention are kept.
    covered = {label for label in remap[0]}
    kept.update(label for label in label_asyms if label not in covered)
    return kept


@functools.lru_cache(maxsize = None)
def _ccd_charges(res_name: str) -> dict[str, int] | None:
    # Only atom names and formal charges are read, so a component whose CCD entry lacks
    # coordinates for some atoms is as good as any other; biotite raises ValueError for those
    # unless told to allow it.
    try:
        residue = info.residue(res_name, allow_missing_coord = True)
    except (KeyError, ValueError):
        return None
    return {atom_name: int(charge) for atom_name, charge in zip(residue.atom_name, residue.charge)}


def assign_charges(atoms: struc.AtomArray, ph: float = C.PH) -> np.ndarray:
    """Formal charges: CCD values, hydride's estimate for the canonical amino acids."""
    charges = np.zeros(atoms.array_length(), dtype = np.int8)
    for res_name in np.unique(atoms.res_name):
        table = _ccd_charges(res_name)
        if not table:
            continue
        positions = np.nonzero(atoms.res_name == res_name)[0]
        charges[positions] = [table.get(atom_name, 0) for atom_name in atoms.atom_name[positions]]

    canonical = np.isin(atoms.res_name, list(C.CANONICAL_AMINO_ACIDS))
    if canonical.any():
        import hydride
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            estimated = hydride.estimate_amino_acid_charges(atoms, ph = ph)
        charges[canonical] = estimated[canonical]
    return charges


def residue_index(atoms: struc.AtomArray) -> np.ndarray:
    """Residue number 0..n-1 of every atom."""
    starts = struc.get_residue_starts(atoms)
    index = np.zeros(atoms.array_length(), dtype = np.int64)
    index[starts[1:]] = 1
    return np.cumsum(index)


def residue_kinds(atoms: struc.AtomArray, entity_kind: dict[str, int], res_index: np.ndarray) -> np.ndarray:
    """C.KIND_* per atom, constant within a residue. Without _entity categories (a file written
    without them) polymer residues are recognised by their CCD names."""
    n_residues = int(res_index[-1]) + 1 if len(res_index) else 0
    first = np.searchsorted(res_index, np.arange(n_residues))
    kinds = np.empty(n_residues, dtype = np.int8)
    heavy_counts = np.bincount(res_index, minlength = n_residues)
    for residue, start in enumerate(first):
        kind = entity_kind.get(atoms.label_entity_id[start], -1)
        if kind < 0 and not entity_kind:
            name = atoms.res_name[start]
            if heavy_counts[residue] > 1 and name in _AMINO_ACID_NAMES:
                kind = C.KIND_PROTEIN
            elif heavy_counts[residue] > 1 and name in _NUCLEOTIDE_NAMES:
                kind = C.KIND_NUCLEOTIDE
        if kind < 0:
            element = atoms.element[start]
            if heavy_counts[residue] == 1:
                kind = C.KIND_METAL if element in C.METAL_ELEMENTS else C.KIND_ION
            else:
                kind = C.KIND_LIGAND
        kinds[residue] = kind
    return kinds[res_index]


# Files from this size on are streamed and cut to the relevant region before biotite sees them: a
# 394 MB capsid (1m4x, 5040 chains, ~17M atoms) does not fit a worker's memory as a whole.
STREAMING_MIN_BYTES = 50 * 2**20
_HEADER_CATEGORIES = ('entity', 'entity_poly', 'pdbx_chain_remapping')
_STREAM_BATCH_ROWS = 200000


def _split_row(line: str) -> list[str]:
    """Tokens of an _atom_site row. Quoted values (atom names with primes) never contain blanks."""
    if '"' not in line and "'" not in line:
        return line.split()
    tokens = []
    i = 0
    n = len(line)
    while i < n:
        character = line[i]
        if character.isspace():
            i += 1
            continue
        if character in '"\'':
            j = line.find(character, i + 1)
            while j != -1 and j + 1 < n and not line[j + 1].isspace():
                j = line.find(character, j + 1)
            if j == -1:
                j = n
            tokens.append(line[i + 1:j])
            i = j + 1
        else:
            j = i
            while j < n and not line[j].isspace():
                j += 1
            tokens.append(line[i:j])
            i = j
    return tokens


class _AtomSiteStream:
    """The header categories of a file, then its _atom_site rows one by one."""

    def __init__(self, path: str):
        self.file = mmcifParser.open_cif_text(path)
        self.extract = mmcifParser.extract_cif_categories(self.file, _HEADER_CATEGORIES, required = (), stop_category = 'atom_site')
        if self.extract.stop_line is None:
            self.file.close()
            raise ValueError('No _atom_site category')
        self.tag_lines = [self.extract.stop_line]
        self.first_row = None
        for line in self.file:
            stripped = line.strip()
            if stripped.startswith('_atom_site.'):
                self.tag_lines.append(line)
                continue
            self.first_row = line
            break
        self.tags = [tag.strip().split()[0].split('.', 1)[1] for tag in self.tag_lines]
        self.is_loop = all(len(tag.split()) == 1 for tag in self.tag_lines)
        self.column = {tag: pos for pos, tag in enumerate(self.tags)}

    def rows(self):
        line = self.first_row
        while line is not None:
            stripped = line.strip()
            if stripped:
                if stripped[0] in '#_' or stripped.startswith(('loop_', 'data_')):
                    return
                yield line
            line = next(self.file, None)

    def close(self):
        self.file.close()


def _streamed_region_block(path: str, radius: float | None = None):
    """A CIF block holding the header categories and only the _atom_site rows of the kept chain
    instances and of the residues within `radius` of them (default: `relevant_radius`); None if
    nothing can be cut."""
    if radius is None:
        radius = C.THRESHOLDS['relevant_radius']

    stream = _AtomSiteStream(path)
    try:
        if not stream.is_loop:
            return None, {}
        remap_lines = stream.extract.categories.get('pdbx_chain_remapping')
        if remap_lines is None:
            return None, {}
        remap_category = mmcifParser.parse_category('pdbx_chain_remapping', remap_lines)
        remap = [mmcifParser.category_column(remap_category, column) for column in _REMAP_COLUMNS]
        kept = _kept_from_remap(remap, list(dict.fromkeys(remap[0])))
        if kept >= set(remap[0]):
            return None, {}
        i_asym = stream.column['label_asym_id']
        i_x, i_y, i_z = stream.column['Cartn_x'], stream.column['Cartn_y'], stream.column['Cartn_z']
        i_model = stream.column.get('pdbx_PDB_model_num')

        # Pass 1: coordinates of the kept instances.
        kept_coords = []
        n_rows = 0
        model = None
        for line in stream.rows():
            tokens = _split_row(line)
            if i_model is not None:
                if model is None:
                    model = tokens[i_model]
                elif tokens[i_model] != model:
                    continue
            n_rows += 1
            if tokens[i_asym] in kept:
                kept_coords.append((float(tokens[i_x]), float(tokens[i_y]), float(tokens[i_z])))
    finally:
        stream.close()
    if not kept_coords:
        return None, {}
    tree = cKDTree(np.array(kept_coords))

    # Pass 2: whole residues with an atom near a kept instance (rows of a residue are contiguous).
    stream = _AtomSiteStream(path)
    keep_lines: list[str] = []
    try:
        i_seq = stream.column.get('auth_seq_id', stream.column.get('label_seq_id'))
        i_ins = stream.column.get('pdbx_PDB_ins_code')
        i_comp = stream.column.get('label_comp_id')
        lines, coords, residues, is_kept = [], [], [], []
        previous = None
        residue = -1

        def flush(cut: int):
            if cut == 0:
                return
            near = np.isfinite(tree.query(np.array(coords[:cut]), k = 1, distance_upper_bound = radius)[0])
            residue_ids = np.array(residues[:cut])
            selected = np.isin(residue_ids, np.unique(residue_ids[near | np.array(is_kept[:cut])]))
            keep_lines.extend(line for line, take in zip(lines[:cut], selected) if take)
            del lines[:cut], coords[:cut], residues[:cut], is_kept[:cut]

        for line in stream.rows():
            tokens = _split_row(line)
            if i_model is not None and model is not None and tokens[i_model] != model:
                continue
            key = (tokens[i_asym], tokens[i_seq] if i_seq is not None else '',
                   tokens[i_ins] if i_ins is not None else '', tokens[i_comp] if i_comp is not None else '')
            if key != previous:
                residue += 1
                previous = key
                if len(lines) >= _STREAM_BATCH_ROWS:
                    flush(len(lines))
            lines.append(line)
            coords.append((float(tokens[i_x]), float(tokens[i_y]), float(tokens[i_z])))
            residues.append(residue)
            is_kept.append(tokens[i_asym] in kept)
        flush(len(lines))
        header = ''.join(stream.extract.categories[name].body() + '#\n'
                         for name in _HEADER_CATEGORIES if name in stream.extract.categories)
        text = ''.join(['data_region\n', header, 'loop_\n', *stream.tag_lines, *keep_lines, '#\n'])
    finally:
        stream.close()

    import io
    block = pdbx.CIFFile.read(io.StringIO(text)).block
    return block, {'atoms_in_file': n_rows, 'streamed': True, 'streamed_rows_kept': len(keep_lines)}


def prepare_assembly(path: str | None = None, add_hydrogens: bool = True, streaming_min_bytes: int = STREAMING_MIN_BYTES,
                     block = None) -> PreparedAssembly:
    """From a file (`path`), or from an already parsed CIF block (models read from SQL)."""
    import os
    stats: dict = {}
    if block is None and os.path.getsize(path) >= streaming_min_bytes:
        block, streamed_stats = _streamed_region_block(path)
        stats.update(streamed_stats)
    if block is None:
        block = read_block(path)
    atoms = pdbx.get_structure(block, model = 1, altloc = 'first', include_bonds = True, extra_fields = _EXTRA_FIELDS)
    stats.setdefault('atoms_in_file', atoms.array_length())

    entity_kind = entity_kinds(block)
    water_entities = [entity_id for entity_id, kind in entity_kind.items() if kind == -2]
    water = np.isin(atoms.res_name, list(C.WATER_NAMES)) | np.isin(atoms.label_entity_id, water_entities)
    hydrogen = np.isin(atoms.element, ('H', 'D'))
    atoms = atoms[~water & ~hydrogen]

    label_asyms = list(dict.fromkeys(atoms.label_asym_id.tolist()))
    kept_asyms = kept_label_asyms(block, label_asyms)
    kept = np.isin(atoms.label_asym_id, list(kept_asyms))

    # Cut away everything that cannot interact with a kept instance, whole residues at a time.
    res_index = residue_index(atoms)
    if not kept.all() and kept.any():
        tree = cKDTree(atoms.coord[kept])
        distance, _ = tree.query(atoms.coord, k = 1, distance_upper_bound = C.THRESHOLDS['relevant_radius'])
        near = np.isfinite(distance)
        region = np.isin(res_index, np.unique(res_index[near]))
        stats['region_restricted'] = True
        atoms = atoms[region]
        kept = kept[region]
        res_index = residue_index(atoms)
    else:
        stats['region_restricted'] = False
    stats['heavy_atoms'] = atoms.array_length()

    atoms.set_annotation('kept', kept)
    atoms.set_annotation('kind', residue_kinds(atoms, entity_kind, res_index))
    atoms.set_annotation('charge', assign_charges(atoms))

    stats['hydrogens_added'] = 0
    if add_hydrogens and atoms.array_length():
        import hydride
        with warnings.catch_warnings(record = True) as caught:
            warnings.simplefilter('always')
            hydrogenated, original = hydride.add_hydrogen(atoms)
        stats['hydride_warnings'] = len(caught)
        stats['hydrogens_added'] = int((~original).sum())
        atoms = hydrogenated
        # hydride fills annotations it does not know with defaults for the new atoms (kept False,
        # empty label ids). It places them behind the heavy atoms of their residue, so the first
        # atom of every residue still carries the residue's values.
        starts = struc.get_residue_starts(atoms)
        res_index = residue_index(atoms)
        for name in ('kept', 'kind', 'label_asym_id', 'label_entity_id', 'label_seq_id'):
            values = atoms.get_annotation(name)
            atoms.set_annotation(name, values[starts][res_index])

    return PreparedAssembly(atoms = atoms, stats = stats)
