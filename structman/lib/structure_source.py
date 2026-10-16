"""Where the structures of the mmCIF pipeline, their RINs and their header data come from.

Structure ids are the part of a search DB target id before the ':' (see MMseqs2.search_databases):

    1A04-2                      PDB biological assembly (mmcifParser.make_assembly_id): the
                                assembly file of the local mirror, else the RCSB's
    AF-P12345-F1                AlphaFold DB monomer          SQL database afdb_mmcif
    afdb_homo_<id>              AlphaFold homodimer           SQL database afdb_homo_mmcif
    afdb_hetero_<id>            AlphaFold heterodimer         SQL database afdb_hetero_mmcif
    P40200_S0__P49961_S0__AF    complex model                 <complex_model_db_path>_mmcif files
    <file name>                 custom DB structure           PDB format file in config.custom_db_path

The models of the SQL databases and the complex models were converted to mmCIF and got their RINs
in `structman update model_dbs` (structman.scripts.createModelDbs); the RINs of the assemblies are
the RIN DB of `structman update rindb` (structman.scripts.createAssemblyRinDb). A structure whose
RIN is missing gets one computed on the spot with structman.lib.rin_engine.

Worker processes import this module; it keeps structman.base_utils.base_utils out of its imports.
"""

from __future__ import annotations

import gzip
import io
import os
import urllib.error
import urllib.request
from dataclasses import dataclass, field

from structman.lib import mmcifParser

KIND_PDB_ASSEMBLY = 'pdb_assembly'
KIND_AFDB = 'afdb'
KIND_AFDB_HOMO = 'afdb_homo'
KIND_AFDB_HETERO = 'afdb_hetero'
KIND_COMPLEX_MODEL = 'complex_model'
KIND_CUSTOM = 'custom'

SQL_MODEL_KINDS = frozenset((KIND_AFDB, KIND_AFDB_HOMO, KIND_AFDB_HETERO))
MODEL_KINDS = frozenset((KIND_AFDB, KIND_AFDB_HOMO, KIND_AFDB_HETERO, KIND_COMPLEX_MODEL, KIND_CUSTOM))

# Models carry no experimental resolution; they are treated like a 2 A structure, as before.
MODEL_RESOLUTION = 2.0

CUSTOM_STRUCTURE_SUFFIXES = ('.pdb', '.pdb.gz', '.ent', '.ent.gz', '.cif', '.cif.gz')

RCSB_ENTRY_URL = 'https://files.rcsb.org/download/{pdb_id}.cif.gz'
DOWNLOAD_TIMEOUT = 60


# ---------------------------------------------------------------------------------------------
# Kinds of structure ids

def custom_structure_file(structure_id: str, config) -> str | None:
    """The custom DB file of a structure id: the first file of config.custom_db_path whose name
    starts with the id, as the PDB format pipeline chose it."""
    folder = getattr(config, 'custom_db_path', None)
    if not folder or not os.path.isdir(folder):
        return None
    for name in sorted(os.listdir(folder)):
        if name.startswith(structure_id) and name.endswith(CUSTOM_STRUCTURE_SUFFIXES):
            return os.path.join(folder, name)
    return None


def kind_of(structure_id: str, config = None) -> str | None:
    """One of the KIND_* constants, None for an id of none of the known forms.

    With a config, a structure of the custom DB is recognised as such first, since custom ids can
    take any form.
    """
    if config is not None and custom_structure_file(structure_id, config) is not None:
        return KIND_CUSTOM
    if mmcifParser.is_assembly_id(structure_id):
        return KIND_PDB_ASSEMBLY
    if structure_id.startswith('afdb_homo_'):
        return KIND_AFDB_HOMO
    if structure_id.startswith('afdb_hetero_'):
        return KIND_AFDB_HETERO
    if structure_id.startswith('AF-'):
        return KIND_AFDB
    if '__' in structure_id:
        return KIND_COMPLEX_MODEL
    return None


def is_model_kind(kind: str | None) -> bool:
    return kind in MODEL_KINDS


# ---------------------------------------------------------------------------------------------
# mmCIF

@dataclass
class MmcifSource:
    """An mmCIF structure, as a local file or as text."""
    structure_id: str
    kind: str
    path: str | None = None
    text: str | None = None
    # The path is the assembly file of the local PDB mirror (a RIN computed from it may be stored).
    mirror_file: bool = False

    def open_lines(self):
        """The file's lines, as a context manager."""
        if self.path is not None:
            return mmcifParser.open_cif_text(self.path)
        return io.StringIO(self.text)


def complex_model_mmcif_folder(config) -> str:
    return f'{config.complex_model_db_path.rstrip("/")}_mmcif'


def complex_model_stem(structure_id: str, config) -> str:
    """<complex_model_db_path>_mmcif/<A[-2:]>/<A[-4:-2]>/<id>, the layout of createModelDbs."""
    first = structure_id.split('_')[0]
    return os.path.join(complex_model_mmcif_folder(config), first[-2:], first[-4:-2], structure_id)


def _read_text(path: str) -> str:
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding = 'utf-8', errors = 'replace') as f:
        return f.read()


def mmcif_source(structure_id: str, config) -> MmcifSource | None:
    """The mmCIF of a structure; None if it cannot be found."""
    kind = kind_of(structure_id, config)

    if kind == KIND_PDB_ASSEMBLY:
        from structman.lib import assembly_retrieval
        pdb_id, assembly_nr = mmcifParser.split_assembly_id(structure_id)
        pdb_path = getattr(config, 'pdb_path', None)
        if pdb_path:
            path = mmcifParser.assembly_file_path(pdb_path, pdb_id, assembly_nr)
            if os.path.isfile(path):
                return MmcifSource(structure_id, kind, path = path, mirror_file = True)
        text = assembly_retrieval.download_assembly(pdb_id, assembly_nr, config)
        if text is None:
            return None
        return MmcifSource(structure_id, kind, text = text)

    if kind in SQL_MODEL_KINDS:
        from structman.lib.database import afdb_retrieval
        text = afdb_retrieval.get_model_mmcif(structure_id, config)
        if text is None:
            return None
        return MmcifSource(structure_id, kind, text = text)

    if kind == KIND_COMPLEX_MODEL:
        path = f'{complex_model_stem(structure_id, config)}.cif.gz'
        if not os.path.isfile(path):
            return None
        return MmcifSource(structure_id, kind, path = path)

    if kind == KIND_CUSTOM:
        path = custom_structure_file(structure_id, config)
        if path.endswith(('.cif', '.cif.gz')):
            return MmcifSource(structure_id, kind, path = path)
        from structman.lib import model_conversion
        cif_text, _ = model_conversion.pdb_to_mmcif(_read_text(path), structure_id)
        return MmcifSource(structure_id, kind, text = cif_text)

    return None


# ---------------------------------------------------------------------------------------------
# RINs

def assembly_rin_db_folder(config) -> str | None:
    """RINdb_assemblies next to config.rin_db_path (createAssemblyRinDb.rin_db_folder)."""
    rin_db_path = getattr(config, 'rin_db_path', None)
    if not rin_db_path:
        return None
    from structman.scripts import createAssemblyRinDb
    return createAssemblyRinDb.rin_db_folder(rin_db_path)


def stored_rin_document(structure_id: str, kind: str, config) -> dict | None:
    from structman.lib.rin_engine import storage

    if kind == KIND_PDB_ASSEMBLY:
        folder = assembly_rin_db_folder(config)
        if folder is None:
            return None
        path = storage.rin_file_path(folder, structure_id)
        return storage.read_rin(path) if os.path.isfile(path) else None

    if kind in SQL_MODEL_KINDS:
        from structman.lib.database import afdb_retrieval
        return afdb_retrieval.get_model_rin(structure_id, config)

    if kind == KIND_COMPLEX_MODEL:
        path = f'{complex_model_stem(structure_id, config)}{storage.SUFFIX}'
        return storage.read_rin(path) if os.path.isfile(path) else None

    return None


def compute_rin_document(source: MmcifSource) -> dict:
    from structman.lib import rin_engine
    if source.path is not None:
        stat = os.stat(source.path)
        return rin_engine.compute_rin(
            path = source.path, structure_id = source.structure_id,
            source = {'name': os.path.basename(source.path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns})
    return rin_engine.compute_rin(cif_text = source.text, structure_id = source.structure_id, source = {'name': source.structure_id})


def rin_document(source: MmcifSource, config, compute_missing: bool = True) -> dict | None:
    """The RIN document of a structure: the stored one, else one computed now.

    A RIN computed from an assembly file of the local mirror is written to the RIN DB, so the next
    run finds it.
    """
    document = stored_rin_document(source.structure_id, source.kind, config)
    if document is not None or not compute_missing:
        return document

    if getattr(config, 'verbosity', 0) >= 3 and getattr(config, 'logger', None) is not None:
        config.logger.info(f'No stored RIN for {source.structure_id}, computing it')
    document = compute_rin_document(source)

    if source.kind == KIND_PDB_ASSEMBLY and source.mirror_file:
        folder = assembly_rin_db_folder(config)
        if folder is not None and os.path.isdir(folder) and os.access(folder, os.W_OK):
            from structman.lib.rin_engine import storage
            try:
                storage.write_rin(storage.rin_file_path(folder, source.structure_id), document)
            except OSError as err:
                config.errorlog.add_warning(f'Could not store the RIN of {source.structure_id}: {err}')
    return document


# ---------------------------------------------------------------------------------------------
# Header data of PDB entries

ENTRY_HEADER_CATEGORIES = (
    'exptl', 'refine', 'reflns', 'em_3d_reconstruction', 'pdbx_coordinate_model', 'struct_ref', 'struct_ref_seq',
)

# Resolution from the first of these that has a value.
_RESOLUTION_ITEMS = (('refine', 'ls_d_res_high'), ('em_3d_reconstruction', 'resolution'), ('reflns', 'd_resolution_high'))

# The values the PDB format pipeline gave (pdbParser.getInfo).
NMR_RESOLUTION = 2.5
UNUSABLE_RESOLUTION = 100.0


@dataclass
class EntryHeader:
    # None when the header could not be read at all.
    resolution: float | None
    methods: list[str] = field(default_factory = list)
    # Author chain of the deposited entry -> UniProt entry name (or accession).
    chain_names: dict[str, str] = field(default_factory = dict)


def _columns(extract: mmcifParser.CategoryExtract, name: str, columns: tuple[str, ...]) -> list[list[str]] | None:
    category_lines = extract.categories.get(name)
    if category_lines is None:
        return None
    category = mmcifParser.parse_category(name, category_lines)
    if any(column not in category for column in columns):
        return None
    return [mmcifParser.category_column(category, column) for column in columns]


def resolution_from_header(extract: mmcifParser.CategoryExtract) -> float:
    """The resolution, with the rules of the PDB format pipeline: NMR and solution scattering 2.5,
    powder diffraction and C-alpha only models 100, no resolution at all 100."""
    exptl = _columns(extract, 'exptl', ('method',))
    methods = [method.upper() for method in exptl[0]] if exptl is not None else []
    if any('NMR' in method or method.startswith('SOLUTION SCATTERING') for method in methods):
        return NMR_RESOLUTION
    if any('POWDER' in method for method in methods):
        return UNUSABLE_RESOLUTION

    coordinate_model = _columns(extract, 'pdbx_coordinate_model', ('type',))
    if coordinate_model is not None and any('CA ATOMS ONLY' in value.upper() for value in coordinate_model[0]):
        return UNUSABLE_RESOLUTION

    for category_name, column_name in _RESOLUTION_ITEMS:
        values = _columns(extract, category_name, (column_name,))
        if values is None:
            continue
        for value in values[0]:
            try:
                return float(value)
            except ValueError:
                continue
    return UNUSABLE_RESOLUTION


def chain_names_from_header(extract: mmcifParser.CategoryExtract) -> dict[str, str]:
    """Author chain -> UniProt entry name, from the UNP references (DBREF of the PDB format)."""
    refs = _columns(extract, 'struct_ref', ('id', 'db_name', 'db_code', 'pdbx_db_accession'))
    ref_seqs = _columns(extract, 'struct_ref_seq', ('ref_id', 'pdbx_strand_id'))
    if refs is None or ref_seqs is None:
        return {}
    names: dict[str, str] = {}
    for ref_id, db_name, db_code, accession in zip(*refs):
        if db_name != 'UNP':
            continue
        name = db_code if db_code not in mmcifParser.UNOBSERVED else accession
        if name not in mmcifParser.UNOBSERVED:
            names[ref_id] = name
    chain_names: dict[str, str] = {}
    for ref_id, strand_id in zip(*ref_seqs):
        if ref_id in names and strand_id not in mmcifParser.UNOBSERVED:
            chain_names.setdefault(strand_id, names[ref_id])
    return chain_names


def _download_entry_header(pdb_id: str, config) -> mmcifParser.CategoryExtract | None:
    """The header categories of an entry from the RCSB, reading the download only up to _atom_site."""
    url = RCSB_ENTRY_URL.format(pdb_id = pdb_id.upper())
    try:
        with urllib.request.urlopen(urllib.request.Request(url), timeout = DOWNLOAD_TIMEOUT) as response:
            with gzip.GzipFile(fileobj = response) as unzipped:
                lines = io.TextIOWrapper(unzipped, encoding = 'utf-8', errors = 'replace')
                return mmcifParser.extract_cif_categories(lines, ENTRY_HEADER_CATEGORIES, stop_category = 'atom_site')
    except (urllib.error.URLError, OSError, ValueError) as err:
        if getattr(config, 'verbosity', 0) >= 3 and getattr(config, 'logger', None) is not None:
            config.logger.info(f'Could not download the header of {pdb_id} ({url}): {err}')
        return None


def entry_header(pdb_id: str, config) -> EntryHeader:
    """Resolution, methods and UniProt names of a PDB entry, from the header of its entry file."""
    extract = None
    pdb_path = getattr(config, 'pdb_path', None)
    if pdb_path:
        path = mmcifParser.entry_file_path(pdb_path, pdb_id)
        if os.path.isfile(path):
            try:
                extract = mmcifParser.read_cif_categories(path, ENTRY_HEADER_CATEGORIES, required = ())
            except (OSError, ValueError) as err:
                config.errorlog.add_warning(f'Could not read the header of {path}: {err}')
    if extract is None:
        extract = _download_entry_header(pdb_id, config)
    if extract is None:
        return EntryHeader(resolution = None)

    exptl = _columns(extract, 'exptl', ('method',))
    return EntryHeader(
        resolution = resolution_from_header(extract),
        methods = exptl[0] if exptl is not None else [],
        chain_names = chain_names_from_header(extract),
    )
