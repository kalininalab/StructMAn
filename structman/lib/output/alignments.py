"""Per-protein alignment matrices, gated by --fetch_alignments.

One TSV per input protein, MSA style: the top row is the protein sequence, every
following row is one annotated structure, and every column is a sequence position.
A cell holds the residue the position was aligned to, as
AA:ResNr:SimpleClass:Recommended, and a query cell holds AA:SimpleClass.

SimpleClass is a non-negative integer code (SIMPLE_CLASS_OUTPUT_CODE) rather than
the string, which would otherwise be repeated in nearly every cell; the first line of
the file is the legend. Recommended is 1 for the one residue per position that the
classification table names as its Recommended Structure, 0 otherwise. Structure rows
are ordered by how many recommended residues they hold, then by coverage.

Everything here is read back out of the database. The alignment itself is
Alignment_N.Alignment (sub_infos: protein position -> residue contig index), the
amino acid and the PDB residue number both come out of Structure_N.Parse_Data, and
the classification out of Structure_N.Simple_Class_Data. Residue_Data, the matrix
those classifications used to be recomputed from, is never selected. The query row's
classification and the recommended residues are the per-position values of the
classification table, handed over by classificationOutput (see position_annotation).
"""
import os

from structman.base_utils.config_class import Config
from structman.lib.database import database
from structman.lib.database.database_core_functions import select
from structman.lib.database.retrieval import get_alignment_data, get_structure_alignment_context
from structman.lib.rin import SIMPLE_CLASS_CODE_OF_STRING, SIMPLE_CLASS_STRING_OF_CODE
from structman.lib.sdsc import mappings as mappings_package
from structman.base_utils import base_utils

NULL_SYMBOL = '-'
CELL_SEPARATOR = ':'


def _safe_filename(protein_id: str) -> str:
    """A protein id that is safe to use as a file name.

    Primary protein ids are mostly uniprot accessions, but an input id can carry a
    path separator or worse, and this writes one file per protein.
    """
    safe = ''.join(char if (char.isalnum() or char in '.-_') else '_' for char in str(protein_id))
    return safe if len(safe) > 0 else 'unnamed_protein'


# The internal simple class codes (rin.SIMPLE_CLASS_CODE_OF_STRING) run from -4 to 7,
# with the locations negative. The file gets them renumbered from 0 in the same
# order, so Surface, Buried, Core come first and the interaction types follow.
# Internal code -1 (a None class, no location) has no string and gets no output code.
SIMPLE_CLASS_OUTPUT_CODE: dict[int, int] = {
    code: output_code for output_code, code in enumerate(
        sorted(code for code, text in SIMPLE_CLASS_STRING_OF_CODE.items() if text is not None))
}


def _class_code_str(code) -> str:
    """An internal simple class code as written into a cell.

    Code -1 stands for a None class (no location could be assigned); it and any
    code without a string render as the null symbol, as the strings used to.
    """
    if code is None:
        return NULL_SYMBOL
    output_code = SIMPLE_CLASS_OUTPUT_CODE.get(int(code))
    if output_code is None:
        return NULL_SYMBOL
    return str(output_code)


def simple_class_legend() -> str:
    """The legend line: every code a cell can hold, with its simple class string."""
    entries = [
        f'{output_code}={SIMPLE_CLASS_STRING_OF_CODE[code]}'
        for code, output_code in sorted(SIMPLE_CLASS_OUTPUT_CODE.items(), key=lambda item: item[1])
    ]
    return f'# SimpleClass codes: {", ".join(entries)}'


def position_annotation(simple_class: str | None, recommended_structure: str | None) -> tuple[int | None, tuple[str, str, int] | None]:
    """What the matrix needs of one classification table row.

    simple_class is the 'Simple Class' value (integrated_features.simple_class), and
    recommended_structure the recommended structure as getPositions returns it,
    'PDB:Chain ResContig:AA' or '-'. ResContig is the 1-based residue contig index,
    the same space the stored alignments point into, so it can be compared with
    sub_infos directly.

    Returns (internal simple class code or None, (pdb_id, chain, res_contig) or
    None); the code is renumbered for output by _class_code_str.
    """
    class_code = SIMPLE_CLASS_CODE_OF_STRING.get(simple_class) if simple_class is not None else None

    recommended = None
    if recommended_structure is not None and recommended_structure not in ('-', ''):
        try:
            struct_id_chain_id, res_contig_aa = recommended_structure.split()
            pdb_id, chain = struct_id_chain_id.split(':')
            res_contig = int(res_contig_aa.split(':')[0])
            if pdb_id != '':
                recommended = (pdb_id, chain, res_contig)
        except ValueError:
            recommended = None

    return class_code, recommended


def _residue_number(seq_res_map, contig_index: int) -> str:
    """The PDB residue number of a residue contig index, as a string.

    seq_res_map entries are ints, or bytes when the residue number carries an
    insertion code - see createTemplateFasta in globalAlignment.
    """
    res_nr = seq_res_map[contig_index]
    if isinstance(res_nr, bytes):
        return res_nr.decode('ascii', errors='replace')
    return str(res_nr)


def build_structure_row(
        sub_infos, seq_res_map, template_seq, simple_class_codes, seq_length: int,
        recommended_contigs: dict[int, int] | None = None) -> list[str]:
    """One structure's cells, one per sequence position (1-based).

    Mirrors the skip conditions of the aggregation mapping kernel: a position past
    the end of the alignment, a position on a template gap, and a contig index that
    runs off the end of the structure's own arrays all render as the null symbol.

    recommended_contigs maps a position to the residue contig index of this
    structure that the position recommends, for the positions that recommend it.
    """
    if recommended_contigs is None:
        recommended_contigs = {}

    cells: list[str] = []

    for pos in range(1, seq_length + 1):
        if pos >= len(sub_infos):
            cells.append(NULL_SYMBOL)
            continue

        res_contig = int(sub_infos[pos])
        if res_contig == -1:
            cells.append(NULL_SYMBOL)
            continue

        contig_index = res_contig - 1
        if contig_index < 0 or contig_index >= len(seq_res_map) or contig_index >= len(template_seq):
            cells.append(NULL_SYMBOL)
            continue

        res_aa = template_seq[contig_index]
        res_nr = _residue_number(seq_res_map, contig_index)

        # None for a structure annotated before Simple_Class_Data existed; there is
        # no fallback to Residue_Data, so the classification part stays empty while
        # the assignment itself is still reported.
        if simple_class_codes is None or contig_index >= len(simple_class_codes):
            simple_class = NULL_SYMBOL
        else:
            simple_class = _class_code_str(simple_class_codes[contig_index])

        recommended = '1' if recommended_contigs.get(pos) == res_contig else '0'

        cells.append(CELL_SEPARATOR.join((res_aa, res_nr, simple_class, recommended)))

    return cells


def count_recommended_cells(sub_infos, seq_res_map, template_seq, seq_length: int,
                            recommended_contigs: dict[int, int] | None) -> int:
    """How many cells of build_structure_row's row carry Recommended = 1.

    Counted from the recommendations rather than from the rendered cells, so rows
    can be ranked before any of them is built. Same conditions as a flagged cell:
    the position is in the sequence and the alignment, and its contig is the
    recommended one and lies within the structure's arrays.
    """
    if not recommended_contigs:
        return 0
    n_residues = min(len(seq_res_map), len(template_seq))
    count = 0
    for pos, rec_res_contig in recommended_contigs.items():
        if pos < 1 or pos > seq_length or pos >= len(sub_infos):
            continue
        if int(sub_infos[pos]) != rec_res_contig:
            continue
        if 1 <= rec_res_contig <= n_residues:
            count += 1
    return count


def write_alignment_matrix(
        outfile: str,
        primary_protein_id: str,
        input_id,
        sequence: str,
        alignments_for_protein: dict,
        structure_context: dict,
        position_annotations: dict[int, tuple] | None = None) -> bool:
    """Write one protein's matrix. False when it had nothing to report.

    position_annotations: {position number: position_annotation(...)} for the
    positions of this protein that the classification table reports.
    """
    seq_length = len(sequence)
    if position_annotations is None:
        position_annotations = {}

    # Grouped once per protein, so every row only looks up its own recommendations.
    recommended_by_structure: dict[tuple[str, str], dict[int, int]] = {}
    for pos, (_, recommended) in position_annotations.items():
        if recommended is None:
            continue
        rec_pdb_id, rec_chain, rec_res_contig = recommended
        recommended_by_structure.setdefault((rec_pdb_id, rec_chain), {})[pos] = rec_res_contig

    # Best templates first, so the top of the file is the most informative part.
    ranked_structures = []
    for structure_db_id, (sub_infos, coverage, seq_id) in alignments_for_protein.items():
        if structure_db_id not in structure_context:
            continue
        pdb_id, chain, seq_res_map, template_seq, _ = structure_context[structure_db_id]
        recommended_count = count_recommended_cells(
            sub_infos, seq_res_map, template_seq, seq_length, recommended_by_structure.get((pdb_id, chain)))
        ranked_structures.append((structure_db_id, sub_infos, coverage, seq_id, recommended_count))

    if len(ranked_structures) == 0:
        return False

    # The structures the classification table recommends most often first. Then
    # coverage, with sequence identity as the last tie breaker: within one protein
    # the identities are nearly all 1.000, so ranking on them leaves the rows in an
    # arbitrary order. How much of the sequence a template actually covers is what
    # separates them.
    ranked_structures.sort(
        key=lambda entry: (
            entry[4],
            entry[2] if entry[2] is not None else -1.0,
            entry[3] if entry[3] is not None else -1.0,
        ),
        reverse=True)

    positions = [str(pos) for pos in range(1, seq_length + 1)]

    query_cells = []
    for pos, aa in enumerate(sequence, start=1):
        class_code = position_annotations.get(pos, (None, None))[0]
        query_cells.append(f'{aa}{CELL_SEPARATOR}{_class_code_str(class_code)}')

    with open(outfile, 'w') as f:
        f.write(simple_class_legend() + '\n')
        f.write(f'# Protein: {primary_protein_id}\tInput ID: {input_id}\tLength: {seq_length}\n')
        f.write(
            f'# Cell format: AA{CELL_SEPARATOR}ResNr{CELL_SEPARATOR}SimpleClass{CELL_SEPARATOR}Recommended '
            f'(query: AA{CELL_SEPARATOR}SimpleClass; SimpleClass codes see first line, Recommended 1/0)\n'
        )
        f.write('\t'.join(['Structure', 'Seq_Id', 'Coverage'] + positions) + '\n')
        f.write('\t'.join(['<query>', NULL_SYMBOL, NULL_SYMBOL] + query_cells) + '\n')

        for structure_db_id, sub_infos, coverage, seq_id, _ in ranked_structures:
            pdb_id, chain, seq_res_map, template_seq, simple_class_codes = structure_context[structure_db_id]

            cells = build_structure_row(
                sub_infos, seq_res_map, template_seq, simple_class_codes, seq_length,
                recommended_contigs=recommended_by_structure.get((pdb_id, chain)))

            seq_id_str = NULL_SYMBOL if seq_id is None else f'{seq_id:.3f}'
            coverage_str = NULL_SYMBOL if coverage is None else f'{coverage:.3f}'

            f.write('\t'.join([f'{pdb_id}:{chain}', seq_id_str, coverage_str] + cells) + '\n')

    return True


def collect_position_annotations(config: Config, session_id) -> dict[int, dict[int, tuple]]:
    """{prot_db_id: {position number: position_annotation(...)}} for a session.

    Only the fallback for when classificationOutput did not run: it unpacks the
    same Position_Data the classification table is written from, with the same
    skip on an unpack failure, so both paths agree with the table.
    """
    results = select(config, ['Position'], 'RS_Position_Session', equals_rows={'Session': session_id})
    position_list, _ = database.getPositions([row[0] for row in results], config)

    position_annotations: dict[int, dict[int, tuple]] = {}
    for row in position_list:
        position_number = row[1]
        prot_db_id = row[2]
        recommended_structure = row[5]
        packed_pos_data = row[13]

        if packed_pos_data is None:
            simple_class = None
        else:
            try:
                packed_mappings, _ = base_utils.unpack(packed_pos_data)
                mappings = base_utils.unpack(packed_mappings)
            except Exception:
                continue
            if mappings is None:
                mappings = mappings_package.Mappings(init_none=True)
            simple_class = mappings.integrated_features.simple_class

        position_annotations.setdefault(prot_db_id, {})[position_number] = position_annotation(
            simple_class, recommended_structure)

    return position_annotations


def create_alignment_matrices(
        config: Config,
        session_id,
        outfolder: str,
        position_annotations: dict[int, dict[int, tuple]] | None = None) -> None:
    """position_annotations: what classificationOutput collected for the session,
    or None to collect it here."""
    results = select(config, ['Protein'], 'RS_Protein_Session', equals_rows={'Session': session_id})
    prot_db_ids = [row[0] for row in results]

    if len(prot_db_ids) == 0:
        if config.verbosity >= 2:
            config.logger.info('create_alignment_matrices: session has no proteins')
        return

    # A set, not the list: getProteinDict tests `prot_id in prot_id_list` per row.
    protein_dict, _ = database.getProteinDict(set(prot_db_ids), session_id, config, includeSequence=True)

    prot_structure_alignment_map = get_alignment_data(config, prot_db_ids)

    structure_db_ids = set()
    for alignments_for_protein in prot_structure_alignment_map.values():
        structure_db_ids.update(alignments_for_protein.keys())

    structure_context = get_structure_alignment_context(config, structure_db_ids)

    if position_annotations is None:
        position_annotations = collect_position_annotations(config, session_id)

    alignment_folder = os.path.join(outfolder, 'alignments')
    os.makedirs(alignment_folder, exist_ok=True)

    written = 0
    skipped = 0

    for prot_db_id in prot_db_ids:
        if prot_db_id not in prot_structure_alignment_map:
            skipped += 1
            continue
        if prot_db_id not in protein_dict:
            skipped += 1
            continue

        protein_entry = protein_dict[prot_db_id]
        primary_protein_id = protein_entry[0]
        input_id = protein_entry[6]
        sequence = protein_entry[9]

        if sequence is None or len(sequence) == 0:
            config.errorlog.add_warning(
                f'create_alignment_matrices: no sequence stored for {primary_protein_id}, skipped'
            )
            skipped += 1
            continue

        outfile = os.path.join(alignment_folder, f'{_safe_filename(primary_protein_id)}.tsv')

        if write_alignment_matrix(
                outfile, primary_protein_id, input_id, sequence,
                prot_structure_alignment_map[prot_db_id], structure_context,
                position_annotations=position_annotations.get(prot_db_id)):
            written += 1
        else:
            skipped += 1

    if config.verbosity >= 1:
        config.logger.info(
            f'create_alignment_matrices: wrote {written} alignment matrices to {alignment_folder} '
            f'({skipped} proteins skipped)'
        )
