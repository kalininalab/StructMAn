import gc
import multiprocessing
import re
import sys
import time
import traceback
from filelock import FileLock

import numpy as np
import ray

from structman.base_utils.base_utils import (
    calc_checksum, pack, unpack, calculate_chunksizes, 
    is_alphafold_model, struct_id_to_bin_nr, add_to_times, print_times,
    reset_logger_for_remotes, sequence_hash
    )
from structman.base_utils.config_class import Config
from structman.base_utils import ray_utils
from structman.base_utils.ray_utils import ray_init, parallel_db_access
from structman.lib import serializedPipeline
from structman.lib import uniprot
from structman.lib.lib_utils import fetch_msa
from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc import complex as complex_package
from structman.lib.sdsc import structure as structure_package
from structman.lib.sdsc import protein as protein_package
from structman.lib.sdsc import position as position_package
from structman.lib.sdsc import mappings as mappings_package
from structman.lib.sdsc.residue import string_decodes
from structman.lib.sdsc import snv as snv_package
from structman.lib.sdsc import indel as indel_package
from structman.lib.sdsc import interface as interface_package
from structman.base_utils import residue_consts
from structman.lib.database.retrieval import getStoredResidues, retrieve_structure_info_from_struct_ids
from structman.lib.database.database_core_functions import binningSelect, select, insert, update, set_null

# called by output.classification
def getProteinDict(prot_id_list, session_id, config, includeSequence=False):
    if len(prot_id_list) == 0:
        return {}, set()

    table = 'RS_Protein_Session'
    columns = ['Protein', 'Session', 'Input_Id', 'Tags']

    if session_id is not None:
        results = select(config, columns, table, equals_rows={'Session': session_id})

        prot_input_dict = {}

        for row in results:
            if row[1] != session_id:
                continue
            prot_id = row[0]
            if prot_id in prot_id_list:
                prot_input_dict[prot_id] = (row[2], row[3])
    else:
        prot_input_dict = set(prot_id_list)

    if not includeSequence:
        rows = ['Protein_Id', 'Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id', 'Error_Code', 'Error', 'Gene']
    else:
        rows = ['Protein_Id', 'Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id', 'Error_Code', 'Error', 'Gene', 'Sequence']
    table = 'Protein'
    results = binningSelect(prot_id_list, rows, table, config)

    gene_db_ids = set()

    protein_dict = {}
    for row in results:
        prot_id = row[0]
        gene_db_id = row[7]
        if gene_db_id is not None:
            gene_db_ids.add(gene_db_id)
        if prot_id in prot_input_dict:
            if session_id is not None:
                input_id, tags = prot_input_dict[prot_id]
            else:
                input_id = None
                tags = None
            if not includeSequence:
                protein_dict[prot_id] = (row[1], row[2], row[3], row[4], row[5], row[6], input_id, tags, gene_db_id, None)
            else:
                protein_dict[prot_id] = (row[1], row[2], row[3], row[4], row[5], row[6], input_id, tags, gene_db_id, row[8])

    return protein_dict, gene_db_ids


# called from serializedPipeline
def geneCheck(genes, config):
    results = select(config, ['Gene_Id', 'Primary_Gene_Id'], 'Gene')
    for row in results:
        gene_db_id = row[0]
        gene_id = row[1]
        if gene_id not in genes:
            continue
        genes[gene_id].database_id = gene_db_id

    if not config.read_only_mode:
        values = []
        for gene_id in genes:
            if genes[gene_id].database_id is not None:
                continue
            gene_name = genes[gene_id].gene_name
            values.append([gene_id, gene_name])

        if len(values) > 0:
            insert('Gene', ['Primary_Gene_Id' ,'Gene_Name'], values, config)

            results = select(config, ['Gene_Id', 'Primary_Gene_Id'], 'Gene')
            for row in results:
                gene_db_id = row[0]
                gene_id = row[1]
                if gene_id not in genes:
                    continue
                genes[gene_id].database_id = gene_db_id
    return

# called from serializedPipeline
# The identity columns of the Protein table, in the order protCheck reads them. The
# joined identifier columns are what lets a stored protein be recognised by any of its
# names, not only the one that happened to win the hierarchy when it was first stored.
PROTEIN_IDENTITY_COLUMNS = [
    'Protein_Id', 'Primary_Protein_Id', 'Sequence', 'Sequence_Hash',
    'Uniprot_Ac', 'Uniprot_Id', 'RefSeq_Ids', 'Ensembl_Ids', 'PDB_Ids', 'Other_Ids',
]

# Offsets of the joined identifier columns within a PROTEIN_IDENTITY_COLUMNS row.
_IDENTIFIER_COLUMN_OFFSETS = {
    protein_package.ID_TYPE_UNIPROT_AC: 4,
    protein_package.ID_TYPE_UNIPROT_ID: 5,
    protein_package.ID_TYPE_REFSEQ: 6,
    protein_package.ID_TYPE_ENSEMBL: 7,
    protein_package.ID_TYPE_PDB: 8,
    protein_package.ID_TYPE_OTHER: 9,
}


def stored_identifiers_of_row(row) -> set[str]:
    """Every identifier a stored Protein row carries, including its primary id."""
    identifiers: set[str] = set()
    if row[1] is not None:
        identifiers.add(row[1])
    for offset in _IDENTIFIER_COLUMN_OFFSETS.values():
        if offset < len(row):
            identifiers |= protein_package.split_identifiers(row[offset])
    return identifiers


def sequence_to_store(protein) -> str | None:
    """What goes into Protein.Sequence: the sequence, or None when there is none - which
    includes the 0 a failed fetch leaves, and a mutant sequence mutate() could not build."""
    sequence = protein.sequence
    if isinstance(sequence, str) and len(sequence) > 0:
        return sequence
    return None


def matchable_identifiers_of_row(row) -> set[str]:
    """The identifiers match_stored_by_identifier may recognise a stored row by: all of them
    for a row without a sequence, none for a row with one.

    A row with a sequence is matched by that sequence or not at all. A chunk protein only
    reaches the identifier match when the sequence match could not place it, so its sequence
    differs from the row's, and a shared name then joins two different proteins. It did:
    a crashed run wrote the identifiers its transcripts mapped to into stored PDB chain rows,
    and the next run merged the real UniProt isoforms of those names - other sequences - into
    the chains, which ended in a KeyError in insertMultiMutations.
    """
    stored_sequence = row[2]
    if stored_sequence is not None and len(stored_sequence) > 0:
        return set()
    return stored_identifiers_of_row(row)


def stored_identifiers_by_type(row) -> dict[str, set[str]]:
    """Same as stored_identifiers_of_row, but keeping the kind of each identifier.

    Needed when the row is written back: an identifier has to go into the column of its
    own kind, and the flattened set does not say which that is.
    """
    by_type: dict[str, set[str]] = {}
    for id_type, offset in _IDENTIFIER_COLUMN_OFFSETS.items():
        if offset >= len(row):
            continue
        identifiers = protein_package.split_identifiers(row[offset])
        if len(identifiers) > 0:
            by_type[id_type] = identifiers
    return by_type


IDENTIFIER_UPDATE_COLUMNS = [
    'Protein_Id', 'Uniprot_Ac', 'Uniprot_Id', 'RefSeq_Ids',
    'Ensembl_Ids', 'PDB_Ids', 'Other_Ids',
]

_IDENTIFIER_VALUE_ORDER = (
    protein_package.ID_TYPE_UNIPROT_AC, protein_package.ID_TYPE_UNIPROT_ID,
    protein_package.ID_TYPE_REFSEQ, protein_package.ID_TYPE_ENSEMBL,
    protein_package.ID_TYPE_PDB, protein_package.ID_TYPE_OTHER,
)


def protein_identifier_values(protein) -> tuple:
    """The joined identifier columns of a protein, in PROTEIN_IDENTITY_COLUMNS order."""
    return tuple(protein.get_joined_identifiers(id_type) for id_type in _IDENTIFIER_VALUE_ORDER)


def stored_identifier_values(stored_identity: dict[str, set[str]]) -> tuple:
    """The joined identifier columns of a stored row, for comparison with a protein's."""
    return tuple(
        protein_package.join_identifiers(stored_identity.get(id_type))
        for id_type in _IDENTIFIER_VALUE_ORDER
    )


def update_stored_identifiers(
        proteins: protein_package.Proteins,
        stored_identity_by_db_id: dict[int, dict[str, set[str]]],
        config: Config) -> int:
    """Write back identifiers the chunk knows but the stored Protein row does not.

    A protein already in the database is never re-inserted, and addProtInfos skips it
    too, so without this an identifier learned from a later input would stay in memory
    for one chunk and then be lost. Folding a duplicate onto a stored entry is exactly
    the case that produces new identifiers for an existing row.

    The stored identifiers are merged into the in-memory protein first. That direction
    matters: writing the column from a protein that had only seen this run's identifiers
    would delete the ones already recorded.
    """
    values = []

    for prot_id in proteins.get_protein_ids():
        protein = proteins[prot_id]
        if not protein.stored or protein.database_id is None:
            continue
        stored_identity = stored_identity_by_db_id.get(protein.database_id)
        if stored_identity is None:
            continue

        stored_flat: set[str] = set()
        for id_type, identifiers in stored_identity.items():
            stored_flat |= identifiers
            for identifier in identifiers:
                protein.add_identifier(identifier, id_type)

        # Rows written while identifier_kind still guessed carry accessions in Other_Ids
        # as well as in their own column. Pruning here makes the write-back clean them.
        protein.drop_redundant_other_identifiers()

        new_identifiers = protein.all_identifiers() - stored_flat
        # A row can also need writing when nothing is new but the pruning above moved an
        # identifier out of Other_Ids.
        redundant_dropped = protein_identifier_values(protein) != stored_identifier_values(stored_identity)
        if len(new_identifiers) == 0 and not redundant_dropped:
            continue

        if config.verbosity >= 2:
            config.logger.info(
                f'Updating stored protein {prot_id} (database id {protein.database_id}) '
                f'with newly seen identifiers: {sorted(new_identifiers)}'
            )

        values.append((protein.database_id,) + protein_identifier_values(protein))

    if len(values) > 0:
        update(config, 'Protein', IDENTIFIER_UPDATE_COLUMNS, values)
        if config.verbosity >= 1:
            config.logger.info(f'Updated the identifiers of {len(values)} stored proteins')

    return len(values)


def match_stored_by_identifier(
        proteins: protein_package.Proteins,
        stored_by_identifier: dict[str, tuple[int, str]],
        genes,
        config: Config) -> int:
    """Recognise still-unstored chunk proteins by any identifier a stored row carries.

    Only reached for proteins the sequence match could not place, so it never overrides
    the sequence - it fills the gap where a stored row has no sequence to match against.
    """
    if len(stored_by_identifier) == 0:
        return 0

    # Decide every match first and apply them in one batch; see
    # Proteins.apply_identity_changes for why fusing one at a time was quadratic.
    present = set(proteins.get_protein_ids())
    renames: dict[str, str] = {}
    fusions: dict[str, str] = {}
    stored_db_ids: dict[str, int] = {}

    for prot_id in list(proteins.get_protein_ids()):
        protein = proteins[prot_id]
        if protein.stored:
            continue
        if is_mutant_protein(protein):
            continue

        hit = None
        for identifier in sorted(protein.all_identifiers()):
            if identifier in stored_by_identifier:
                hit = stored_by_identifier[identifier]
                break
        if hit is None:
            continue

        stored_db_id, stored_prot_id = hit
        if stored_prot_id == prot_id:
            continue

        if config.verbosity >= 4:
            config.logger.info(
                f'Recognised {prot_id} as the stored protein {stored_prot_id} '
                f'(database id {stored_db_id}) by a shared identifier'
            )

        record_stored_match(prot_id, stored_prot_id, stored_db_id, present, renames, fusions, stored_db_ids)

    apply_stored_matches(proteins, renames, fusions, stored_db_ids, genes, config)

    return len(renames) + len(fusions)


def record_stored_match(
        prot_id: str,
        stored_prot_id: str,
        stored_db_id: int,
        present: set[str],
        renames: dict[str, str],
        fusions: dict[str, str],
        stored_db_ids: dict[str, int]) -> None:
    """Decide how a chunk protein joins the stored protein it matched.

    If the stored id is already in the chunk - as a protein of its own, or because an
    earlier match in this pass renamed onto it - the protein fuses into it. Otherwise it
    is renamed onto that id and becomes the stored protein. Exactly the choice the
    one-at-a-time loop made, with `present` standing in for the protein map it mutated.
    """
    present.discard(prot_id)
    if stored_prot_id in present:
        fusions[prot_id] = stored_prot_id
    else:
        renames[prot_id] = stored_prot_id
        present.add(stored_prot_id)
        stored_db_ids[stored_prot_id] = stored_db_id


def apply_stored_matches(proteins, renames, fusions, stored_db_ids, genes, config: Config) -> None:
    if len(renames) == 0 and len(fusions) == 0:
        return
    proteins.apply_identity_changes(renames, fusions, config, genes = genes)
    for stored_prot_id, stored_db_id in stored_db_ids.items():
        if stored_prot_id in proteins.get_protein_ids():
            proteins[stored_prot_id].database_id = stored_db_id
            proteins[stored_prot_id].stored = True


def is_mutant_protein(protein) -> bool:
    return sdsc_utils.is_mutant_ac(protein.primary_protein_id)


def group_duplicate_sequences(proteins: protein_package.Proteins, config: Config) -> list[list[str]]:
    """Classes of chunk proteins that carry byte-identical sequences.

    Bucketing by CRC32 first keeps this at one hash per protein instead of comparing
    every pair, but a 32 bit hash collides on different sequences often enough to
    matter, so a bucket is only a list of *candidates*: the sequences themselves decide.
    Returns one list per class of two or more, the rest are not duplicates.
    """
    buckets: dict[int, list[str]] = {}

    for prot_id in proteins.get_protein_ids():
        protein = proteins[prot_id]
        if not protein.is_duplicate_candidate():
            continue
        seq_hash = protein.get_sequence_hash()
        if seq_hash is None:
            continue
        if seq_hash not in buckets:
            buckets[seq_hash] = []
        buckets[seq_hash].append(prot_id)

    duplicate_classes: list[list[str]] = []

    for seq_hash, candidates in buckets.items():
        if len(candidates) < 2:
            continue

        # Split the bucket by actual sequence. This is what makes crc32 safe here.
        classes_by_sequence: dict[str, list[str]] = {}
        for prot_id in candidates:
            sequence = proteins[prot_id].sequence
            if sequence not in classes_by_sequence:
                classes_by_sequence[sequence] = []
            classes_by_sequence[sequence].append(prot_id)

        if config.verbosity >= 4 and len(classes_by_sequence) > 1:
            config.logger.info(
                f'Sequence hash collision on {seq_hash}: {len(candidates)} candidates '
                f'resolved into {len(classes_by_sequence)} distinct sequences'
            )

        for prot_ids in classes_by_sequence.values():
            if len(prot_ids) > 1:
                duplicate_classes.append(prot_ids)

    return duplicate_classes


def warn_on_ambiguous_fusion(proteins: protein_package.Proteins, prot_ids: list[str], config: Config) -> None:
    """Two distinct UniProt accessions on one sequence means something upstream went wrong.

    Fusing across id *kinds* - an accession and the Ensembl id of the same protein -
    is the expected case and stays silent.
    """
    accessions = {proteins[prot_id].u_ac for prot_id in prot_ids if proteins[prot_id].u_ac is not None}
    if len(accessions) > 1:
        config.errorlog.add_warning(
            f'Identical sequences for distinct Uniprot accessions: {sorted(accessions)}. '
            f'Fusing them into one entry, but this points at an upstream mapping problem.'
        )


def fuse_duplicate_sequences_in_chunk(proteins: protein_package.Proteins, genes, config: Config) -> int:
    """Pass 1: fuse proteins of this chunk that share a sequence. Returns how many were dropped."""
    # All fusions of the pass go through one apply_identity_changes call, which shares a
    # single reference sweep between them instead of one sweep per fusion.
    fusions: dict[str, str] = {}

    for prot_ids in group_duplicate_sequences(proteins, config):
        canonical_id = protein_package.canonical_protein_id(proteins.get_protein_map(), prot_ids)

        warn_on_ambiguous_fusion(proteins, prot_ids, config)

        for prot_id in prot_ids:
            if prot_id == canonical_id:
                continue
            if config.verbosity >= 4:
                config.logger.info(
                    f'Duplicate sequence: fusing {prot_id} into {canonical_id} (same chunk)'
                )
            fusions[prot_id] = canonical_id

    if len(fusions) > 0:
        proteins.apply_identity_changes({}, fusions, config, genes = genes)

    return len(fusions)


def fuse_duplicates_against_database(
        proteins: protein_package.Proteins,
        stored_by_hash: dict[int, list[tuple[int, str, str]]],
        genes,
        config: Config) -> int:
    """Pass 2: fold chunk proteins onto proteins already in the database.

    A chunk protein that duplicates a stored one takes over the stored protein's
    primary id, so its positions are annotated and reported under the identifier the
    database already knows and no second Protein row is created for the same sequence.
    """
    if len(stored_by_hash) == 0:
        return 0

    present = set(proteins.get_protein_ids())
    renames: dict[str, str] = {}
    fusions: dict[str, str] = {}
    stored_db_ids: dict[str, int] = {}
    matched_in_place = 0

    for prot_id in list(proteins.get_protein_ids()):
        protein = proteins[prot_id]
        if protein.stored:
            if config.verbosity >= 6:
                config.logger.info(
                    f' In fuse_duplicates_against_database: {prot_id} {protein.stored=}'
                )
            continue
        if not protein.is_duplicate_candidate():
            if config.verbosity >= 6:
                config.logger.info(
                    f' In fuse_duplicates_against_database: {prot_id} is not dup candidate'
                )
            continue

        seq_hash = protein.get_sequence_hash()
        if seq_hash not in stored_by_hash:
            if config.verbosity >= 6:
                config.logger.info(
                    f' In fuse_duplicates_against_database: {prot_id} {seq_hash not in stored_by_hash=} {seq_hash=}'
                )
            continue

        if config.verbosity >= 6:
            config.logger.info(
                f' In fuse_duplicates_against_database: {prot_id} {seq_hash=} {len(stored_by_hash[seq_hash])=}'
            )

        for stored_db_id, stored_prot_id, stored_sequence in stored_by_hash[seq_hash]:
            if config.verbosity >= 6:
                config.logger.info(
                    f' In fuse_duplicates_against_database inner loop: {prot_id=} {stored_db_id=} {stored_prot_id=} '
                )
            if stored_sequence != protein.sequence:
                continue
            if stored_prot_id == prot_id:
                # Same id and same sequence: only reachable for proteins the primary-id
                # match skips (fasta entries), so take the stored row as it is.
                protein.stored = True
                protein.database_id = stored_db_id
                matched_in_place += 1
                break

            if config.verbosity >= 4:
                config.logger.info(
                    f'Duplicate sequence: fusing {prot_id} into stored protein '
                    f'{stored_prot_id} (database id {stored_db_id})'
                )
            warn_on_ambiguous_fusion_against_stored(protein, stored_prot_id, stored_db_id, config)

            record_stored_match(prot_id, stored_prot_id, stored_db_id, present, renames, fusions, stored_db_ids)
            break

    apply_stored_matches(proteins, renames, fusions, stored_db_ids, genes, config)

    return len(renames) + len(fusions) + matched_in_place


def warn_on_ambiguous_fusion_against_stored(protein, stored_prot_id: str, stored_db_id: int, config: Config) -> None:
    """The database-side counterpart of warn_on_ambiguous_fusion.

    Only the stored primary id is known here, not the stored row's Uniprot_Ac, so the
    warning fires when both sides look like accessions and disagree.
    """
    if protein.u_ac is None or protein.u_ac == stored_prot_id:
        return
    if sdsc_utils.is_mutant_ac(stored_prot_id) or stored_prot_id.count(':') > 0:
        return
    if stored_prot_id[:2] in ('NP', 'NM', 'XP', 'XM') or stored_prot_id[:3] == 'ENS':
        return
    config.errorlog.add_warning(
        f'Identical sequences for distinct Uniprot accessions: {protein.u_ac} and the '
        f'stored {stored_prot_id} (database id {stored_db_id}). Fusing them into one '
        f'entry, but this points at an upstream mapping problem.'
    )


def protCheck(proteins: protein_package.Proteins, genes, session_id, config: Config) -> bool:
    # Scan the database for stored Proteins
    if proteins.isEmpty():
        return False

    # Pass 1 of the sequence deduplication: collapse proteins of this chunk that carry
    # the same sequence, before any of them can reach the Protein table. Wildtype
    # sequences are already known here: sequenceScan runs before core(), and core()
    # fetches whatever it left without one right before protCheck.
    fused_in_chunk = fuse_duplicate_sequences_in_chunk(proteins, genes, config)
    if fused_in_chunk > 0 and config.verbosity >= 1:
        config.logger.info(f'Sequence deduplication fused {fused_in_chunk} proteins within the chunk')

    if proteins.isEmpty():
        return False

    results = select(config, PROTEIN_IDENTITY_COLUMNS, 'Protein')

    if config.verbosity >= 3:
        config.logger.info(f'Just after protCheck selection {len(proteins.protein_map)=} {len(results)=}')

    prot_id_list = set()
    prot_ids_mutants_excluded = set()

    max_database_id = 0
    matched_by_primary_id = 0

    # Two indexes over the stored rows for pass 2, both built from the select above so
    # the database side costs no extra query: by sequence hash, which is the identity,
    # and by identifier, which catches stored rows whose Sequence was never written.
    stored_by_hash: dict[int, list[tuple[int, str, str]]] = {}
    stored_by_identifier: dict[str, tuple[int, str]] = {}
    # What each stored row already records, so the write-back below can tell a genuinely
    # new identifier from one the row already has.
    stored_identity_by_db_id: dict[int, dict[str, set[str]]] = {}

    for row in results:
        database_id = row[0]
        prot_id = row[1]
        stored_sequence = row[2]

        if stored_sequence is not None and len(stored_sequence) > 0:
            stored_hash = row[3]
            if stored_hash is None:
                # Written before Sequence_Hash existed, or by a path that only set the
                # sequence. Hashing it here keeps pass 2 working on such rows.
                stored_hash = sequence_hash(stored_sequence)
            else:
                stored_hash = int(stored_hash)
            if stored_hash not in stored_by_hash:
                stored_by_hash[stored_hash] = []
            stored_by_hash[stored_hash].append((database_id, prot_id, stored_sequence))

        for identifier in matchable_identifiers_of_row(row):
            stored_by_identifier[identifier] = (database_id, prot_id)

        stored_identity_by_db_id[database_id] = stored_identifiers_by_type(row)

        max_database_id = max(max_database_id, database_id)

        if not proteins.contains(prot_id):
            continue
        # A fasta entry id is a label the user chose, not an identity: two runs can put
        # different sequences under the same header. Such a protein is matched by its
        # sequence in pass 2, which also covers the rerun of the same entry.
        if proteins[prot_id].sequence_from_input:
            continue
        matched_by_primary_id += 1
        proteins.set_protein_stored(prot_id, True)
        proteins[prot_id].database_id = database_id
        proteins.set_protein_sequence(prot_id, stored_sequence)

        prot_id_list.add(database_id)
        if not sdsc_utils.is_mutant_ac(prot_id):
            prot_ids_mutants_excluded.add(database_id)

    if config.verbosity >= 4:
        config.logger.info(
            f' Before fuse_duplicates_against_database: {len(stored_by_hash)=}'
        )

    # Pass 2: fold the still-unstored proteins of this chunk onto stored proteins with
    # the same sequence but a different identifier.
    fused_against_db = fuse_duplicates_against_database(proteins, stored_by_hash, genes, config)

    # Pass 3: a stored row whose Sequence was never written cannot be matched by
    # sequence at all - addProtInfos never refreshes a stored protein, so a row inserted
    # when its sequence fetch failed keeps Sequence NULL forever. Recognising it by any
    # of its identifiers stops it collecting a duplicate on every later run.
    matched_by_identifier = match_stored_by_identifier(proteins, stored_by_identifier, genes, config)

    if fused_against_db > 0 or matched_by_identifier > 0:
        if config.verbosity >= 1:
            config.logger.info(
                f'Sequence deduplication fused {fused_against_db} proteins onto stored entries, '
                f'{matched_by_identifier} more were recognised by identifier'
            )
        for prot_id in proteins.get_protein_ids():
            protein = proteins[prot_id]
            if protein.stored and protein.database_id is not None:
                prot_id_list.add(protein.database_id)
                # indelCheck reads this set to find the stored indels of each wildtype.
                if not is_mutant_protein(protein):
                    prot_ids_mutants_excluded.add(protein.database_id)

    # Write back identifiers the chunk knows but the stored rows do not. Deliberately
    # not restricted to the proteins the passes above just fused: a protein matched by
    # its primary id can bring a new alias too, and it would be lost the same way.
    if not config.read_only_mode:
        update_stored_identifiers(proteins, stored_identity_by_db_id, config)

    del stored_by_hash
    del stored_by_identifier
    del stored_identity_by_db_id

    proteins.set_stored_ids(prot_id_list, prot_ids_mutants_excluded)

    if config.verbosity >= 3:
        config.logger.info(f'In protCheck  {len(prot_id_list)=} {len(proteins.stored_ids)=}')

    prot_ids = proteins.get_protein_ids()

    # Insert the new proteins into the database
    new_prots = set()
    
    for prot_id in prot_ids:
        if not proteins.is_protein_stored(prot_id):
            new_prots.add(prot_id)

    # Whether the input was found in the database at all, in every mode. Only an id count
    # at verbosity 3 said so before, which made a run that matched nothing look like one
    # that never looked.
    if config.verbosity >= 1:
        n_proteins = len(prot_ids)
        config.logger.info(
            f'protCheck: {n_proteins:_} proteins, {n_proteins - len(new_prots):_} already stored '
            f'({matched_by_primary_id:_} by id, {fused_against_db:_} by sequence, '
            f'{matched_by_identifier:_} by identifier), {len(new_prots):_} new')

    if config.read_only_mode:
        for prot_id in new_prots:
            proteins.remove_protein(prot_id)
        if len(proteins.stored_ids) == 0:
            return False
    else:
        if len(new_prots) > 0:
            values = []
            for prot_id in new_prots:
                u_ac = proteins.get_u_ac(prot_id)
                ref_id = proteins.get_ref_id(prot_id)
                u_id = proteins.get_u_id(prot_id)
                mut_type = proteins[prot_id].mutant_type
                gene_id = proteins[prot_id].gene
                if gene_id is not None:
                    if gene_id in genes:
                        gene_db_id = genes[gene_id].database_id
                    else:
                        gene_db_id = None
                else:
                    gene_db_id = None
                if mut_type is None:
                    mutant_positions_str = None
                else:
                    if proteins[prot_id].sav_positions is None:
                        sav_string = ''
                    else:
                        sav_string = ','.join([str(x) for x in proteins[prot_id].sav_positions])

                    if proteins[prot_id].insertion_positions is None:
                        insertion_string = ''
                    else:
                        insertion_parts = []
                        for ins in proteins[prot_id].insertion_positions:
                            insertion_parts.append(','.join([str(x) for x in ins]))
                        insertion_string = '!'.join(insertion_parts)

                    if proteins[prot_id].deletion_flanks is None:
                        deletion_string = ''
                    else:
                        deletion_parts = []
                        for deletion_fl in proteins[prot_id].deletion_flanks:
                            lf, rf = deletion_fl
                            if lf is None:
                                lfs = '_'
                            else:
                                lfs = str(lf)
                            if rf is None:
                                rfs = '_'
                            else:
                                rfs = str(rf)
                            deletion_parts.append(f'{lfs},{rfs}')
                        deletion_string = '!'.join(deletion_parts)

                    mutant_positions_str = f'{sav_string};{insertion_string};{deletion_string}'
                # The joined identifier columns carry every name this sequence is known
                # by; the singular Uniprot_Ac / Uniprot_Id / RefSeq_Ids values above are
                # subsumed by them, so they are written from the sets instead.
                joined_ac, joined_u_id, joined_refseq, joined_ensembl, joined_pdb, joined_other = (
                    protein_identifier_values(proteins[prot_id]))
                # The sequence goes in with the row. Every sequence is known by now - fetched
                # before protCheck, mutants built by mutate() before it too - so no row exists
                # without one, not even for a chunk that dies right after this insert.
                sequence = sequence_to_store(proteins[prot_id])
                values.append((
                    gene_db_id, prot_id, joined_ac, joined_refseq, joined_u_id,
                    session_id, mut_type, mutant_positions_str,
                    joined_ensembl, joined_pdb, joined_other,
                    sequence, sequence_hash(sequence)))

            insert('Protein', [
                'Gene', 'Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id',
                'Original_Session', 'Mutant_Type', 'Mutant_Positions',
                'Ensembl_Ids', 'PDB_Ids', 'Other_Ids',
                'Sequence', 'Sequence_Hash'], values, config)

            # Retrieve the Protein-Ids from the new added proteins

            db, cursor = config.getDB()

            sql = f"SELECT Protein_Id,Primary_Protein_Id FROM Protein WHERE Protein_Id > {max_database_id}"
            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except BaseException as err:
                [e, f, g] = sys.exc_info()
                raise NameError(f"Error in protCheck: {err=} {sql=}\n{e}\n{f}\n{g}")
            db.close()
        else:
            results = ()

        new_ids = set()
        for row in results:
            database_id = row[0]
            prot_id = row[1]
            if not proteins.contains(prot_id):
                continue
            proteins[prot_id].database_id = database_id

            if database_id not in prot_id_list:
                new_ids.add(database_id)

        proteins.set_not_stored_ids(new_ids)

        # The chunk's first new id is what 'structman database repair' needs when the chunk
        # dies before its alignments are written: from there on, its proteins cannot be told
        # apart from ones without templates.
        if len(new_ids) > 0 and config.verbosity >= 1:
            config.logger.info(
                f'{len(new_ids)} new proteins stored, Protein_Id {min(new_ids)}..{max(new_ids)}. If this chunk '
                f'does not finish: structman database repair --session_id {session_id} '
                f'--from_protein_id {min(new_ids)} --apply')

    proteins.generate_id_map()

    if config.isolate_ppi:
        db_id_set: set[int] = set()

    # Insert the Protein-Session-Connections into the database
    values = []
    for prot_id in prot_ids:
        database_id = proteins.get_protein_database_id(prot_id)
        # Every input identifier that fused onto this sequence, not just the one that
        # happened to be first. This is what the output's Input Protein ID column shows.
        input_id = proteins[prot_id].get_joined_input_ids()
        if input_id is None:
            input_id = proteins[prot_id].input_id
        tags_str = ','.join(proteins[prot_id].tags)
        values.append((database_id, session_id, input_id, tags_str))
        if config.isolate_ppi:
            db_id_set.add(database_id)

    insert('RS_Protein_Session', ['Protein', 'Session', 'Input_Id', 'Tags'], values, config)

    if config.isolate_ppi:
        proteins.database_ids = db_id_set

    return True


def getLastAutoInc(config):
    db, cursor = config.getDB()

    if config.db_address == '-':
        sql = "SELECT last_insert_rowid()"
    else:
        sql = "SELECT LAST_INSERT_ID()"

    cursor.execute(sql)
    results = cursor.fetchall()
    for row in results:
        auto_id = row[0]
    db.commit()

    return auto_id


# called by babel
def getUniprotAcFromId(Protein_Id, db, cursor):
    sql = "SELECT Uniprot_Ac FROM Protein WHERE Protein_Id = '%s'" % (str(Protein_Id))

    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("NameError in getUniprotFromId")
        db.rollback()

    if results == ():
        return 0
    row = results[0]
    return row[0]


# called by babel
def getAAChange(mutation_id, db, cursor):
    sql = "SELECT Amino_Acid_Change FROM Position WHERE Position_Id = '%s'" % str(mutation_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getAAChange: %s" % sql)
        db.rollback()

    if results == ():
        return ""
    for row in results:
        aac = row[0]
        return aac
    return ""


def _joined_ids(ids) -> str:
    """How a set of SNV or indel ids is stored in Multi_Mutation.SNVs / .Indels. Sorted, so
    the same set always gives the same string - the rows are matched back by it."""
    return ','.join(sorted(str(x) for x in ids))


def _joined_id_set(text) -> set[str]:
    """The ids in a Multi_Mutation.SNVs / .Indels value, as strings, without the empty one
    an empty set is stored as."""
    if text is None:
        return set()
    return {part for part in text.split(',') if part != ''}


def insertMultiMutations(proteins, genes, session, config: Config):

    if config.verbosity >= 3:
        config.logger.info(f'insertMultiMutations: {len(proteins.multi_mutations)}')

    t0 = time.time()
    db_ids = []
    for wt_prot in proteins.multi_mutations:
        db_ids.append(proteins.get_protein_database_id(wt_prot))

    if len(db_ids) > 0:
        t10 = time.time()
        columns = ['Wildtype_Protein', 'Multi_Mutation_Id', 'SNVs', 'Indels']
        table = 'Multi_Mutation'

        results = binningSelect(db_ids, columns, table, config)

        t11 = time.time()
        if config.verbosity >= 3:
            config.logger.info(f'insertMultiMutations part 1.1: {t11-t10}')

        for row in results:
            wt_database_id = row[0]
            wt_prot_id = proteins.getU_acByDbId(wt_database_id)
            if wt_prot_id not in proteins.multi_mutations:
                continue
            # Compared as strings on both sides: the row holds the ids as text, and comparing
            # that with the integer ids never matched, so every rerun stored them again.
            snv_db_ids = _joined_id_set(row[2])
            indel_db_ids = _joined_id_set(row[3])

            for multi_mutation in proteins.multi_mutations[wt_prot_id]:
                if snv_db_ids != {str(x) for x in multi_mutation.get_snv_db_ids(proteins)}:
                    continue

                if indel_db_ids != {str(x) for x in multi_mutation.get_indel_db_ids()}:
                    continue
                multi_mutation.database_id = row[1]
                multi_mutation.stored = True

        t12 = time.time()
        if config.verbosity >= 3:
            config.logger.info(f'insertMultiMutations part 1.2: {t12-t11}')

    t1 = time.time()
    if config.verbosity >= 3:
        config.logger.info(f'insertMultiMutations part 1: {t1-t0}')

    obj_map = {}
    if config.read_only_mode:
        wt_prot_db_ids = db_ids
    else:
        wt_prot_db_ids = []
        values = []
        
        for wt_prot in proteins.multi_mutations:
            wt_db_id = proteins.get_protein_database_id(wt_prot)
            wt_prot_db_ids.append(wt_db_id)
            for multi_mutation in proteins.multi_mutations[wt_prot]:
                if multi_mutation.stored:
                    continue
                if multi_mutation.mut_prot is not None:
                    if multi_mutation.mut_prot in proteins.protein_map:
                        mut_db_id = proteins[multi_mutation.mut_prot].database_id
                    else:
                        mut_db_id = None
                else:
                    mut_db_id = None
                snv_db_ids = _joined_ids(multi_mutation.get_snv_db_ids(proteins))
                indel_db_ids = _joined_ids(multi_mutation.get_indel_db_ids())

                mm_obj_id = (wt_db_id, snv_db_ids, indel_db_ids)

                # Identical multi mutations of one wildtype share a row. Sequence
                # deduplication produces them: identical transcripts carrying the same
                # variants end up on one survivor. Each used to add its own row, and all but
                # the last were left without an id.
                if mm_obj_id in obj_map:
                    obj_map[mm_obj_id].append(multi_mutation)
                    continue
                obj_map[mm_obj_id] = [multi_mutation]
                values.append((wt_db_id, mut_db_id, snv_db_ids, indel_db_ids))

        if len(values) > 0:
            columns = ['Wildtype_Protein', 'Mutant_Protein', 'SNVs', 'Indels']
            table = 'Multi_Mutation'
            insert(table, columns, values, config)

    t2 = time.time()
    if config.verbosity >= 3:
        config.logger.info(f'insertMultiMutations part 2: {t2-t1}')

    columns = ['Wildtype_Protein', 'Multi_Mutation_Id', 'SNVs', 'Indels']
    table = 'Multi_Mutation'

    results = binningSelect(wt_prot_db_ids, columns, table, config)

    t3 = time.time()
    if config.verbosity >= 3:
        config.logger.info(f'insertMultiMutations part 3: {t3-t2}')

    for row in results:
        wt_db_id = row[0]
        wt_prot_id = proteins.getU_acByDbId(wt_db_id)

        mm_obj_id = (wt_db_id, row[2], row[3])
        if mm_obj_id not in obj_map:
            continue

        for multi_mutation in obj_map[mm_obj_id]:
            multi_mutation.database_id = row[1]

    t4 = time.time()
    if config.verbosity >= 3:
        config.logger.info(f'insertMultiMutations part 4: {t4-t3}')

    columns = ['Multi_Mutation', 'Session', 'Tags']
    table = 'RS_Multi_Mutation_Session'

    values = []
    isoform_values = []
    for wt_prot in proteins.multi_mutations:
        wt_db_id = proteins.get_protein_database_id(wt_prot)
        for multi_mutation in proteins.multi_mutations[wt_prot]:
            if multi_mutation.tags is not None:
                mm_str = ','.join(multi_mutation.tags)
            else:
                mm_str = None
            values.append((multi_mutation.database_id, session, mm_str))

            if proteins[wt_prot].gene is not None:
                mut_db_id = proteins[multi_mutation.mut_prot].database_id
                gene_db_id = genes[proteins[wt_prot].gene].database_id
                isoform_values.append((wt_db_id, mut_db_id, gene_db_id, session, multi_mutation.database_id))

    t5 = time.time()
    if config.verbosity >= 3:
        config.logger.info(f'insertMultiMutations part 5: {t5-t4}')

    if len(values) > 0:
        insert(table, columns, values, config)

    t6 = time.time()
    if config.verbosity >= 3:
        config.logger.info(f'insertMultiMutations part 6: {t6-t5}')

    if len(isoform_values) > 0:
        insert('RS_Isoform', ['Protein_A', 'Protein_B', 'Gene', 'Session', 'Multi_Mutation'], isoform_values, config)

    t7 = time.time()
    if config.verbosity >= 3:
        config.logger.info(f'insertMultiMutations part 1: {t7-t6}')

    return


# called by serializedPipeline
def indelCheck(proteins, session, config):
    if len(proteins.indels) == 0:
        return
    stored_wt_ids = proteins.get_stored_ids(exclude_indel_mutants=True)

    if len(stored_wt_ids) > 0:
        rows = ['Wildtype_Protein', 'Indel_Notation', 'Indel_Id']
        table = 'Indel'

        results = binningSelect(stored_wt_ids, rows, table, config)

        for row in results:
            wt_prot_id = row[0]
            indel_notation = row[1]
            database_id = row[2]
            u_ac = proteins.getU_acByDbId(wt_prot_id)
            # stored_wt_ids holds every stored wildtype of the chunk, so the select also
            # returns indels that were stored by earlier runs but are not part of this one.
            if u_ac not in proteins.indels:
                continue
            if indel_notation not in proteins.indels[u_ac]:
                continue
            proteins.indels[u_ac][indel_notation].set_database_id(database_id)
            proteins.indels[u_ac][indel_notation].set_stored(True)

    # Wildtypes of the indels inserted below, whose new Indel_Ids still have to be read back.
    new_wt_ids = set()
    if not config.read_only_mode:
        values = []
        for u_ac in proteins.indels:
            for indel_notation in proteins.indels[u_ac]:
                if proteins.indels[u_ac][indel_notation].stored:
                    continue
                wt_prot_id = proteins.get_protein_database_id(proteins.indels[u_ac][indel_notation].wt_prot)
                mut_prot_id = proteins.get_protein_database_id(proteins.indels[u_ac][indel_notation].mut_prot)
                if wt_prot_id is None:
                    config.errorlog.add_error('Wildtype protein id is not allowed to be None')
                    continue
                values.append((wt_prot_id, mut_prot_id, indel_notation))
                new_wt_ids.add(wt_prot_id)

        if len(values) > 0:
            columns = ['Wildtype_Protein', 'Mutant_Protein', 'Indel_Notation']
            insert('Indel', columns, values, config)

    if len(new_wt_ids) > 0:
        rows = ['Wildtype_Protein', 'Indel_Notation', 'Indel_Id']
        table = 'Indel'
        results = binningSelect(new_wt_ids, rows, table, config)
        for row in results:
            indel_notation = row[1]
            wt_prot_id = row[0]
            u_ac = proteins.getU_acByDbId(wt_prot_id)
            if u_ac not in proteins.indels:
                continue
            if indel_notation not in proteins.indels[u_ac]:
                continue
            if proteins.indels[u_ac][indel_notation].stored:
                continue
            database_id = row[2]
            proteins.indels[u_ac][indel_notation].set_database_id(database_id)

    # One session row per indel of this run, stored or new. Built from the indels rather
    # than from the selected rows, which would also link every duplicate row an indel has.
    session_values = []
    for u_ac in proteins.indels:
        for indel_obj in proteins.indels[u_ac].values():
            if indel_obj.database_id is None:
                continue
            session_values.append((session, indel_obj.database_id, ','.join(indel_obj.tags)))

    if len(session_values) > 0:
        columns = ['Session', 'Indel', 'Tags']
        insert('RS_Indel_Session', columns, session_values, config)

@ray.remote(max_calls=1)
def remote_pos_check(*args, **kwargs):
    # Thin wrapper so the body can also be called directly, without ray, by the
    # serial branch of its dispatcher (see ray_utils.parallel_db_access).
    return pos_check_shard(*args, **kwargs)


def pos_check_shard(prot_db_ids: list[int], db_nr: int, config: Config):
    """The position ids of one shard's proteins.

    Recommended_Structure_Data used to come along and was unpacked per position to set
    Position.recommended_structure. Nothing in the pipeline reads that: the output stage
    rebuilds its proteins through proteinsFromDb, and Position.get_recommended_structure
    falls back to the mappings. On a whole-PDB chunk it meant unpacking tens of millions
    of blobs for nothing.
    """
    columns = ['Protein', 'Position_Number', 'Position_Id']
    table = f'Position_{db_nr}'

    results = binningSelect(prot_db_ids, columns, table, config)

    return db_nr, results


@ray.remote(max_calls=1, max_retries=10)
def remote_insert_positions(values, db_nr, config: Config):
    return insert_positions(values, db_nr, config)


def insert_positions(values, db_nr, config: Config):
    """Insert one shard's new positions. Reading the ids back is a separate step now:
    stored and new proteins are selected together, once, after every shard is inserted."""
    if len(values) > 0:
        columns = ['Protein', 'Position_Number', 'Residue_Id', 'Wildtype_Residue']
        insert(f'Position_{db_nr}', columns, values, config)


@ray.remote(max_calls=1, max_retries=10)
def remote_insert_position_sessions(values, config: Config):
    return insert_position_sessions(values, config)


def insert_position_sessions(values, config: Config):
    if len(values) > 0:
        insert('RS_Position_Session', ['Session', 'Position', 'Tag'], values, config)


# How many connections write RS_Position_Session at once. Plain appends with no
# uniqueness to enforce, so the slices are independent; kept small because they all hit
# the same table.
MAX_PARALLEL_POSITION_SESSION_INSERTS = 8

# called by serializedPipeline
def positionCheck(proteins: protein_package.Proteins, database_session: int, config: Config):
    """Store this chunk's positions and link them, and their SNVs, to the session.

    Rests on one property of the pipeline: Protein.popNone reports uni_pos for every
    protein, so sequenceScan gives each one all of its residues as positions and they are
    inserted together. A stored protein therefore has all of its positions stored and a
    new one has none of them, which is why nothing here decides that position by position
    any more.

    The walks over the position objects are what this step costs - 31.6M of them on a
    whole-PDB chunk, at ~29 us each on a loaded machine against 0.23 us on an idle one,
    so it is the touching that is expensive, not the work done per position. There used
    to be four such walks plus a row walk; now new proteins are walked once to build the
    insert and every protein's rows are walked once to do everything else.
    """
    times = []
    ta = time.time()

    protein_map = proteins.protein_map
    prot_ids = proteins.get_protein_ids()
    stored_prot_ids = proteins.get_stored_ids()

    # Every position of a stored protein is stored, so there is nothing left to check
    # for it - which is exactly what the old per-position loop concluded.
    for prot_id in prot_ids:
        if protein_map[prot_id].stored or config.read_only_mode:
            proteins.set_completely_stored(prot_id)

    ta = add_to_times(times, ta, 'mark stored proteins')

    if parallel_db_access(config):
        # ray.put would otherwise fall into ray's auto_init_hook and run a bare
        # ray.init() - every core on the machine, default object store - instead of the
        # sizing ray_init derives from the config; on a busy machine that times out.
        # ray_init returns immediately when ray is already up.
        ray_init(config, stop_existing = False)
        conf_store = ray.put(config)
    else:
        conf_store = None

    if not config.read_only_mode:
        # Only proteins that are not stored yet contribute positions; the stored ones are
        # not walked at all.
        value_bins = [[] for _ in range(20)]
        N_of_values = 0

        # These loops build lists of tens of millions of tuples while the process already
        # holds a very large live heap. Any generation-2 collection they trigger has to
        # walk that entire heap and can reclaim nothing, since none of these tuples can be
        # part of a cycle.
        gc_was_enabled = gc.isenabled()
        gc.disable()
        try:
            for prot_id in prot_ids:
                prot_obj = protein_map[prot_id]
                if prot_obj.stored:
                    continue
                prot_database_id = prot_obj.database_id
                value_bin = value_bins[prot_database_id % 20]
                # Index 0 is a placeholder and never a real position.
                for pos_obj in prot_obj.positions[1:]:
                    if pos_obj is None:
                        continue
                    value_bin.append(
                        (prot_database_id, pos_obj.pos, pos_obj.pdb_res_nr, pos_obj.wt_aa))
                    N_of_values += 1
        finally:
            if gc_was_enabled:
                gc.enable()

        ta = add_to_times(times, ta, f'prepare pos value bins {N_of_values=}')

        if parallel_db_access(config) and N_of_values > 10_000:
            ray.get([remote_insert_positions.remote(values, db_nr, conf_store)
                     for db_nr, values in enumerate(value_bins)])
        else:
            for db_nr, values in enumerate(value_bins):
                insert_positions(values, db_nr, config)

        del value_bins

        ta = add_to_times(times, ta, f'position insert {N_of_values=} {parallel_db_access(config)=}')

    # One id select for stored and new proteins together. The stored ones used to be
    # selected once before the insert and a second time after it.
    fused_prot_ids = proteins.get_not_stored_ids() | stored_prot_ids
    bins = [[] for _ in range(20)]
    for prot_db_id in fused_prot_ids:
        bins[prot_db_id % 20].append(prot_db_id)
    shards = [(bin_nr, prot_db_id_bin) for bin_nr, prot_db_id_bin in enumerate(bins) if len(prot_db_id_bin) > 0]

    if parallel_db_access(config):
        db_results = ray.get([remote_pos_check.remote(prot_db_id_bin, bin_nr, conf_store)
                              for bin_nr, prot_db_id_bin in shards])
    else:
        db_results = [pos_check_shard(prot_db_id_bin, bin_nr, config) for bin_nr, prot_db_id_bin in shards]

    ta = add_to_times(times, ta, f'select position ids {len(fused_prot_ids)=}')

    # The single walk over the returned rows: it sets every position's database id, builds
    # the session rows, and collects the few positions that carry SNVs. Separate walks
    # over all positions used to do each of those.
    id_map = proteins.id_map
    session_values = []
    append_session_value = session_values.append
    snv_positions: list = []
    append_snv_position = snv_positions.append
    n_rows = 0

    gc_was_enabled = gc.isenabled()
    gc.disable()
    try:
        for db_nr, results in db_results:
            for row in results:
                prot_database_id = row[0]
                pos = row[1]

                positions = protein_map[id_map[prot_database_id]].positions
                if pos >= len(positions):
                    continue
                pos_obj = positions[pos]
                if pos_obj is None:
                    continue

                n_rows += 1
                internal_pos_db_id = 100 * row[2] + db_nr
                pos_obj.database_id = internal_pos_db_id

                if not pos_obj.session_less:
                    append_session_value((
                        database_session, internal_pos_db_id,
                        None if pos_obj.pos_tags is None else ','.join(pos_obj.pos_tags)))

                if pos_obj.mut_aas:
                    append_snv_position(pos_obj)
    finally:
        if gc_was_enabled:
            gc.enable()

    del db_results

    ta = add_to_times(times, ta, f'process pos db results {n_rows=} {len(snv_positions)=}')

    if config.verbosity >= 4:
        # Costs a walk over every position, so verbosity 4 only. A position with no row
        # belongs to a stored protein whose sequence grew since it was stored; it gets no
        # session row. Normally zero.
        expected_rows = sum(
            1 for prot_id in prot_ids
            for pos_obj in protein_map[prot_id].positions[1:] if pos_obj is not None)
        if expected_rows != n_rows:
            config.logger.info(
                f'positionCheck: {expected_rows - n_rows} positions of this chunk have no '
                f'database row and are left out of the session')

    if len(snv_positions) > 0:
        # Only the positions that carry SNVs are queried, instead of every position id of
        # the chunk. full_select is gone with it: a handful of scattered ids turned into a
        # range over the whole SNV table. Resolving stored ids happens in read_only mode
        # too, where it is the only way the session rows below get their SNV ids.
        snv_position_ids = [pos_obj.database_id for pos_obj in snv_positions]
        stored_snvs = {}
        for row in binningSelect(snv_position_ids, ['Position', 'New_AA', 'SNV_Id'], 'SNV', config):
            stored_snvs[(row[0], row[1])] = row[2]

        missing = []
        for pos_obj in snv_positions:
            position_database_id = pos_obj.database_id
            for new_aa, snv_obj in pos_obj.mut_aas.items():
                stored_snv_id = stored_snvs.get((position_database_id, new_aa))
                if stored_snv_id is None:
                    missing.append((position_database_id, new_aa))
                else:
                    snv_obj.database_id = stored_snv_id

        ta = add_to_times(times, ta, f'select stored SNVs {len(stored_snvs)=}')

        # Only SNVs that are not stored yet. SNV has no unique key on (Position, New_AA),
        # so inserting the stored ones again duplicated them on every re-run.
        if len(missing) > 0 and not config.read_only_mode:
            insert('SNV', ['Position', 'New_AA'], missing, config)

            for row in binningSelect(snv_position_ids, ['Position', 'New_AA', 'SNV_Id'], 'SNV', config):
                stored_snvs[(row[0], row[1])] = row[2]

            for pos_obj in snv_positions:
                for new_aa, snv_obj in pos_obj.mut_aas.items():
                    if snv_obj.database_id is None:
                        snv_obj.database_id = stored_snvs.get((pos_obj.database_id, new_aa))

        ta = add_to_times(times, ta, f'insert SNVs {len(missing)=}')

    if len(snv_positions) > 0:
        values = [
            (database_session, snv_obj.database_id, ','.join(snv_obj.tags))
            for pos_obj in snv_positions
            for snv_obj in pos_obj.mut_aas.values()
        ]
        if len(values) > 0:
            insert('RS_SNV_Session', ['Session', 'SNV', 'Tag'], values, config)

        ta = add_to_times(times, ta, f'insert SNV sessions {len(snv_positions)=}')

    # The session rows are plain appends with no uniqueness to enforce, so slices of them
    # can go in parallel; one connection wrote all 31.6M of them before.
    n_session_values = len(session_values)
    if n_session_values > 0:
        if parallel_db_access(config) and n_session_values > 10_000:
            n_slices = max(1, min(MAX_PARALLEL_POSITION_SESSION_INSERTS, config.proc_n))
            slice_size = -(-n_session_values // n_slices)
            ray.get([
                remote_insert_position_sessions.remote(session_values[left:left + slice_size], conf_store)
                for left in range(0, n_session_values, slice_size)
            ])
        else:
            insert_position_sessions(session_values, config)

    del session_values

    ta = add_to_times(times, ta, f'insert pos sessions {n_session_values=}')

    if config.verbosity >= 2:
        print_times(times, 'positionCheck', logger=config.logger)

    return None


def backgroundInsertMS(values, config):
    insert('RS_Position_Session', ['Session', 'Position', 'Tag'], values, config)


# called by babel
def getPDB(template_id, db, cursor):
    sql = "SELECT Name FROM Template WHERE Template_Id = '%s'" % str(template_id)

    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        db.rollback()

    return results[0][0]


# called by babel
def checkMutationSession(mutation_id, session_id, db, cursor):
    sql = "SELECT Position FROM RS_Position_Session WHERE Position = '%s' AND Session = '%s'" % (str(mutation_id), str(session_id))
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in checkMutationSession")
        db.rollback()
    if results == ():
        return False
    else:
        return True

@ray.remote(max_calls=1)
def remote_getPositions(*args, **kwargs):
    # Thin wrapper so the body can also be called directly, without ray, by the
    # serial branch of its dispatcher (see ray_utils.parallel_db_access).
    return get_positions_shard(*args, **kwargs)


def get_positions_shard(db_nr: int, config: Config, pos_db_ids: list[int]):
    reset_logger_for_remotes(config)
    table = f'Position_{db_nr}'
    columns = ['Position_Id', 'Position_Number', 'Protein', 'Residue_Id', 'Wildtype_Residue', 'Recommended_Structure_Data', 'Position_Data']
    results: list[list] = binningSelect(pos_db_ids, columns, table, config)

    unpacked_results = []

    struct_ids = set()

    for row in results:
        mut_id = row[0]
        internal_pos_db_id = 100*mut_id + db_nr

        wt_aa = row[4]

        if row[5] is not None:
            (recommended_structure_str, max_seq_structure_str) = unpack(row[5])
        else:
            recommended_structure_str = None
            max_seq_structure_str = None

        if config.verbosity >= 5:
            config.logger.info(f'In remote_getPositions before process_rec_str {recommended_structure_str=}')

        recommended_structure, seq_id, cov, resolution = sdsc_utils.process_recommend_structure_str(recommended_structure_str)

        if config.verbosity >= 5:
            config.logger.info(f'In remote_getPositions before {recommended_structure=}')

        recommended_structure_words = recommended_structure.split(':')
        struct_ids.add(recommended_structure_words[0])
        if len(recommended_structure_words) == 3:
            struct_aa = string_decodes[35][int(float(recommended_structure_words[2]))]
            recommended_structure = ':'.join(recommended_structure_words[:2]) + f':{struct_aa}'

        if config.verbosity >= 5:
            config.logger.info(f'In remote_getPositions after {recommended_structure=}')

        max_seq_structure, max_seq_seq_id, max_seq_cov, max_seq_resolution = sdsc_utils.process_recommend_structure_str(max_seq_structure_str)

        max_seq_structure_words = max_seq_structure.split(':')
        struct_ids.add(max_seq_structure_words[0])
        if len(max_seq_structure_words) == 3:
            struct_aa = string_decodes[35][int(float(max_seq_structure_words[2]))]
            max_seq_structure = ':'.join(max_seq_structure_words[:2]) + f':{struct_aa}'
        
        unpacked_results.append(
            (
                internal_pos_db_id,
                row[1],
                row[2],
                row[3],
                wt_aa,
                recommended_structure,
                seq_id,
                cov,
                resolution,
                max_seq_structure,
                max_seq_seq_id,
                max_seq_cov,
                max_seq_resolution,
                row[6]
            ))
        
    return unpacked_results, struct_ids

# called by output.classification
def getPositions(
        position_id_list: list[int],
        config: Config
        ) -> tuple[list[tuple], set[str]]:
    if len(position_id_list) == 0:
        return [], []

    pos_db_bins = [[] for _ in range(20)]
    for internal_pos_db in position_id_list:
        db_nr = internal_pos_db % 100
        pos_db_id = internal_pos_db // 100
        pos_db_bins[db_nr].append(pos_db_id)

    if parallel_db_access(config):
        # Without this, ray.put falls into ray's auto_init_hook and starts an
        # unconfigured cluster. ray_init returns immediately if ray is already up.
        ray_init(config, stop_existing = False)
        conf_dump = ray.put(config)
        process_ids = [remote_getPositions.remote(db_nr, conf_dump, pos_db_bins[db_nr]) for db_nr in range(20)]
        db_results = ray.get(process_ids)
    else:
        db_results = [get_positions_shard(db_nr, config, pos_db_bins[db_nr]) for db_nr in range(20)]


    total_results = []
    total_struct_ids = set()
    for results, struct_ids in db_results:
        total_results += results
        total_struct_ids = total_struct_ids.union(struct_ids)

    return total_results, total_struct_ids


# called by serializedPipeline
def updateSession(session_id, time, config):
    update(config, 'Session', ['Session_Id', 'End'], [(str(session_id), str(time))])


# called by serializedPipeline
def insertSession(time, ori_filename, config):
    
    checksum = calc_checksum(ori_filename)

    values = [(ori_filename, checksum, str(time))]
   
    if len(values) > 0:
        insert('Session', ['Input_File', 'Checksum', 'Start'], values, config)

    session_id = getSessionId(ori_filename, config)

    return session_id

def getSessionId(infile: str, config: Config):
    checksum = calc_checksum(infile)

    table = 'Session'
    columns = ['Session_Id', 'Checksum']

    if config.verbosity >= 3:
        config.logger.info(f'Looking for infile {infile} with checksum {checksum} in getSessionId')

    results = select(config, columns, table, equals_rows={'Input_File': infile})

    if len(results) == 0:
        return None

    for row in results:
        if checksum == row[1]:
            return row[0]

    config.logger.info(f'Error in getSessionId:\nChecksum of infile {infile} not found in database:\n{results}')

    return None


# called by babel
# called by output
def getSessionFile(session_id, config):
    table = 'Session'
    columns = ['Input_File']

    results = select(config, columns, table, equals_rows={'Session_Id': session_id})

    if results == ():
        return None
    try:
        return results[0][0]
    except:
        return None

# called by serializedPipeline
def addProtInfos(proteins, config):
    """GO terms and pathways of the new proteins.

    The sequence and its hash used to be written here as well, in an update of the rows
    protCheck had inserted without them. protCheck now inserts them with the row.
    """
    go_term_ids = set()
    reac_ids = set()
    u_acs = proteins.get_protein_ids()

    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            go_term_ids.add(go_id)
        for reac_id in pathways:
            reac_ids.add(reac_id)

    stored_go_terms = {}
    if len(go_term_ids) > 0:
        results = select(config, ['Id', 'GO_Term_Id'], 'GO_Term', in_rows={'ID': go_term_ids})

        for row in results:
            stored_go_terms[row[0]] = row[1]

    stored_pathways = {}
    if len(reac_ids) > 0:
        results = select(config, ['Reactome_Id', 'Pathway_Id'], 'Pathway', in_rows={'Reactome_Id': reac_ids})

        for row in results:
            stored_pathways[row[0]] = row[1]

    new_go_terms = {}
    go_values = []
    new_pathways = {}
    pathway_values = []
    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            if go_id not in stored_go_terms:
                new_go_terms[go_id] = go_terms[go_id]
                go_values.append((go_terms[go_id].replace("'", "''"), go_id))
                stored_go_terms[go_id] = None
        for reac_id in pathways:
            if reac_id not in stored_pathways:
                new_pathways[reac_id] = pathways[reac_id]
                pathway_values.append((pathways[reac_id].replace("'", "''"), reac_id))
                stored_pathways[reac_id] = None

    if len(go_values) > 0:
        insert('GO_Term', ['Name', 'Id'], go_values, config)

    if len(pathway_values) > 0:
        insert('Pathway', ['Name', 'Reactome_Id'], pathway_values, config)

    stored_go_terms = {}
    if len(go_term_ids) > 0:
        results = select(config, ['Id', 'GO_Term_Id', ], 'GO_Term', in_rows={'Id': go_term_ids})

        for row in results:
            stored_go_terms[row[0]] = row[1]

    stored_pathways = {}
    if len(reac_ids) > 0:
        results = select(config, ['Reactome_Id', 'Pathway_Id'], 'Pathway', in_rows={'Reactome_Id': reac_ids})

        for row in results:
            stored_pathways[row[0]] = row[1]

    go_values = []
    pathway_values = []
    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        prot_id = proteins.get_protein_database_id(u_ac)
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            go_values.append((prot_id, stored_go_terms[go_id]))
        for reac_id in pathways:
            pathway_values.append((prot_id, stored_pathways[reac_id]))

    if len(go_values) > 0:
        insert('RS_Protein_GO_Term', ['Protein', 'GO_Term'], go_values, config)

    if len(pathway_values) > 0:
        insert('RS_Protein_Pathway', ['Protein', 'Pathway'], pathway_values, config)


def getComplexMap(config, pdb_ids=None):
    if pdb_ids is not None:
        if len(pdb_ids) == 0:
            return {}
    table = 'Complex'
    rows = ['Complex_Id', 'PDB', 'Resolution', 'Chains', 'Homooligomers', 'Ligand_Profile', 'Metal_Profile', 'Ion_Profile', 'Chain_Chain_Profile']

    if pdb_ids is None:
        results = select(config, rows, table)
    elif len(pdb_ids) > 50:
        results = select(config, rows, table)
    else:
        results = select(config, rows, table, in_rows = {'PDB': pdb_ids})

    complex_map = {}
    for row in results:
        pdb_id = row[1]
        if pdb_ids is not None:
            if pdb_id not in pdb_ids:
                continue
        comp_id = row[0]
        resolution = row[2]
        chains_str = row[3]
        homooligomers = row[4]
        lig_profile_str = row[5]
        metal_profile_str = row[6]
        ion_profile_str = row[7]
        cc_profile_str = row[8]

        lig_profile = {}
        if lig_profile_str != '':
            for lig_profile_part in lig_profile_str.split(','):
                part1, part2 = lig_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                lig_profile[(chain, res)] = int(deg), float(score)

        metal_profile = {}

        if metal_profile_str != '':
            for metal_profile_part in metal_profile_str.split(','):
                part1, part2 = metal_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                metal_profile[(chain, res)] = int(deg), float(score)

        ion_profile = {}

        if ion_profile_str != '':
            for ion_profile_part in ion_profile_str.split(','):
                part1, part2 = ion_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                ion_profile[(chain, res)] = int(deg), float(score)

        cc_profile = {}

        if cc_profile_str != '':
            for cc_profile_part in cc_profile_str.split(','):
                part1, part2 = cc_profile_part.split(':')
                chain, chain_b = part1.split('_')
                deg, score = part2.split('_')

                cc_profile[(chain, chain_b)] = int(deg), float(score)

        complex_map[pdb_id] = (comp_id, resolution, chains_str, homooligomers, lig_profile, metal_profile, ion_profile, cc_profile)

    return complex_map


@ray.remote(max_calls=1)
def remote_check_structures(*args, **kwargs):
    # Thin wrapper so the body can also be called directly, without ray, by the
    # serial branch of its dispatcher (see ray_utils.parallel_db_access).
    return check_structures_shard(*args, **kwargs)


def check_structures_shard(structure_ids: list[str], db_nr: int, config: Config):
    table = f'Structure_{db_nr}'
    columns = ['Structure_Id', 'PDB', 'Chain']
    if len(structure_ids) > 5000:
        return select(config, columns, table, locked=True), db_nr
    elif len(structure_ids) > 0:
        return select(config, columns, table, in_rows={'PDB': structure_ids}, locked=True), db_nr
    else:
        return [], db_nr

def structureCheck(proteins: protein_package.Proteins, config: Config):
    if config.verbosity >= 4:
        config.logger.info(f'Call of structureCheck {len(proteins.structures)=}')

    struct_bins = [[] for _ in range(20)]
    for structure_id in proteins.structures:
        bin_nr = struct_id_to_bin_nr(structure_id)
        struct_bins[bin_nr].append(structure_id)

    shards = [(db_nr, structure_ids) for db_nr, structure_ids in enumerate(struct_bins) if len(structure_ids) > 0]

    if parallel_db_access(config):
        # Without this, ray.put falls into ray's auto_init_hook and starts an
        # unconfigured cluster. ray_init returns immediately if ray is already up.
        ray_init(config, stop_existing = False)
        conf_dump = ray.put(config)
        process_ids = [remote_check_structures.remote(structure_ids, db_nr, conf_dump) for db_nr, structure_ids in shards]
        db_results = ray.get(process_ids)
    else:
        db_results = [check_structures_shard(structure_ids, db_nr, config) for db_nr, structure_ids in shards]

    internal_db_ids = []

    for results, db_nr in db_results:
        for row in results:
            s_id = row[0]
            pdb_id = row[1]
            chain = row[2]
            if config.verbosity >= 6:
                config.logger.info(f'{pdb_id} {chain} processed in structureCheck')
            if pdb_id not in proteins.structures:
                continue

            internal_structure_db_id = 100*s_id + db_nr
            internal_db_ids.append(internal_structure_db_id)

            if chain not in proteins.structures[pdb_id]:
                struct = structure_package.Structure(pdb_id, chain, database_id=internal_structure_db_id)
                proteins.add_structure(pdb_id, chain, struct)
            else:
                proteins.set_structure_db_id(pdb_id, chain, internal_structure_db_id)
                proteins.set_structure_stored(pdb_id, chain, True)  # all structures, mapped or not go into this dictionary, this is important for not reinserting residues from interacting structures
            if config.verbosity >= 6:
                config.logger.info(f'{pdb_id} {chain} set as stored in structureCheck')
            #proteins.complexes[pdb_id].stored = True
    if len(internal_db_ids) > 0:
        proteins.max_internal_struct_db_id = max(internal_db_ids)
    else:
        proteins.max_internal_struct_db_id = 0

def draw_complexes(config: Config, proteins, stored_complexes=[], draw_all=False, get_stable_objs = False, update_infos = False):
    results = select(config, ['Complex_Id', 'PDB', 'Resolution', 'Chains', 'Homooligomers', 'Chaintypes'], 'Complex')

    for row in results:
        pdb_id = row[1]
        if pdb_id not in stored_complexes and not draw_all:
            continue
        chaintype_map = unpack(row[5])

        try:
            resolution = float(row[2])
        except TypeError as err:
            config.errorlog.add_warning(f'{err=}, invalid resolution value in DB for: {row=}')
            resolution = 2.0

        compl: complex_package.Complex
        if get_stable_objs:
            compl = complex_package._Complex(pdb_id, chaintype_map=chaintype_map, resolution=resolution, chains_str=row[3], stored=True, database_id=row[0], homomers_str=row[4])

        else:
            
            compl = complex_package.Complex(pdb_id, chaintype_map=chaintype_map, resolution=resolution, chains_str=row[3], stored=True, database_id=row[0], homomers_str=row[4])

        if update_infos:
            compl.getPage(config, self_update=True)

        proteins.add_complex(pdb_id, compl)

def remoteInsertStructures(config: Config, structs_to_insert: dict[str, dict[str, tuple[tuple, str]]], compl_id_dict):

    struct_bins = [{} for _ in range(20)]
    for structure_id in structs_to_insert:
        bin_nr = struct_id_to_bin_nr(structure_id)
        struct_bins[bin_nr][structure_id] = structs_to_insert[structure_id]


    db_results: list[list[tuple[int, str, str]]] = []
    for db_nr, structs_to_insert in enumerate(struct_bins):

        db_lock_file = f'{config.tmp_folder}/remoteInsertStructures_{db_nr}_lock.lock'
        with FileLock(db_lock_file):
            results = retrieve_structure_info_from_struct_ids(config, structs_to_insert.keys(), db_nr, cols = ['PDB'])

        for row in results:
            struct_id = row[0]
            if struct_id in structs_to_insert:
                del structs_to_insert[struct_id]

        values = []
        for struct_id in structs_to_insert:
            for chain in structs_to_insert[struct_id]:
                parse_data = pack(structs_to_insert[struct_id][chain][0])
                homomers = structs_to_insert[struct_id][chain][1]
                compl_db_id, resolution = compl_id_dict[struct_id]
                values.append((struct_id, chain, parse_data, compl_db_id, resolution, homomers))

        with FileLock(db_lock_file):
            if len(values) > 0:
                insert(f'Structure_{db_nr}', ['PDB', 'Chain', 'Parse_Data', 'Complex', 'Resolution', 'Homooligomer'], values, config)

            results = retrieve_structure_info_from_struct_ids(config, structs_to_insert.keys(), db_nr, cols = ['Structure_Id', 'PDB', 'Chain'])
        db_results.append(results)

    return db_results

def get_parse_data_from_struct_ids(struct_ids: list[str], config: Config):
    struct_bins = [[] for _ in range(20)]
    for structure_id in struct_ids:
        bin_nr = struct_id_to_bin_nr(structure_id)
        struct_bins[bin_nr].append(structure_id)

    db_results = []
    for db_nr, struct_id_list in enumerate(struct_bins):
        results = retrieve_structure_info_from_struct_ids(config, struct_id_list, db_nr, cols = ['Structure_Id', 'PDB', 'Chain', 'Parse_Data'])
        db_results += results

    return db_results


# called by serializedPipeline
def insertStructures(structurelist: dict[str, set[str]], proteins, config):

    if len(structurelist) == 0:
        return

    t0 = time.time()

    table = 'Structure'
    rows = ['Structure_Id', 'PDB', 'Chain']
    results = select(config, rows, table, in_rows={'PDB': structurelist.keys()})

    t1 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f'insertStructures part 1: {t1-t0} {len(structurelist)=}')

    stored_complexes = set()

    for row in results:
        pdb_id = row[1]
        chain = row[2]
        stored_complexes.add(pdb_id)
        if not proteins.contains_structure(pdb_id, chain):
            continue
        
        s_id = row[0]
        proteins.set_structure_db_id(pdb_id, chain, s_id)
        proteins.set_structure_stored(pdb_id, chain, True)  # all structures, mapped or not go into this dictionary, this is important for not reinserting residues from interacting structures
        
        if pdb_id in structurelist:
            if chain in structurelist[pdb_id]:
                structurelist[pdb_id].remove(chain)
                if len(structurelist[pdb_id]) == 0:
                    del structurelist[pdb_id]

    t2 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f'insertStructures part 2: {t2-t1} {len(stored_complexes)=}')

    draw_complexes(config, proteins, stored_complexes=stored_complexes)

    t3 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f'insertStructures part 3: {t3-t2}')

    values = []

    for pdb_id in structurelist:
        for chain in structurelist[pdb_id]:
            oligos = proteins.structures[pdb_id][chain].oligo
            oligos = ''.join(oligos)
            values.append((pdb_id, chain, oligos))

    t4 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f'insertStructures part 4: {t4-t3} {len(values)=}')

    if len(values) > 0:
        insert('Structure', ['PDB', 'Chain', 'Homooligomer'], values, config)

    t5 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f'insertStructures part 5: {t5-t4}')

    if len(structurelist) > 0:
        table = 'Structure'
        rows = ['Structure_Id', 'PDB', 'Chain']
        results = select(config, rows, table, in_rows={'PDB': structurelist.keys()})

        for row in results:
            s_id = row[0]
            pdb_id = row[1]
            chain = row[2]

            if pdb_id not in structurelist:
                continue
            if chain not in structurelist[pdb_id]:
                continue

            proteins.set_structure_db_id(pdb_id, chain, s_id)

    t6 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f'insertStructures part 6: {t6-t5}')


# called by serializedPipeline
def insertAlignments(values, config: Config, db_lock = None):

    if config.verbosity >= 5:
        t0 = time.time()

    seperated_values = [[] for _ in range(20)]

    for value_list in values:
        prot_db_id = value_list[0]
        seperated_values[prot_db_id % 20].append(value_list)

    for db_nr, values in enumerate(seperated_values):
        if len(values) > 0:
            insert(f'Alignment_{db_nr}', ['Protein', 'Structure', 'Sequence_Identity', 'Coverage', 'Alignment', 'Backmap'], values, config, locked=True, db_lock=db_lock)
    if config.verbosity >= 5:
        t1 = time.time()
        config.logger.info(f'Time for insertAlignments: {t1 - t0} {len(values)=}')


# Structure ids per UPDATE statement: well inside SQLite's 32k parameter cap.
_PRUNE_IN_CHUNK = 10_000


@ray.remote(max_calls=1)
def remote_prune_backmaps_shard(*args, **kwargs):
    # Thin wrapper so the body can also be called directly, without ray, by the
    # serial branch of its dispatcher (see ray_utils.parallel_db_access).
    return prune_backmaps_shard(*args, **kwargs)


def prune_backmaps_shard(
        db_nr: int,
        prot_db_ids: list[int],
        multi_chain_struct_ids: np.ndarray,
        config: Config) -> tuple[int, int, int]:
    """Clears the backmaps of one shard's proteins that select_interface_annotations leaves out.

    multi_chain_struct_ids is the sorted array of the internal structure db ids whose
    complex has more than one protein chain. Rows of the same (protein, structure)
    collapse onto the last one, as in retrieval.get_alignment_data, which is where
    para_classify reads the same selection's inputs from.

    Returns (proteins with alignments, proteins with backmaps cleared, backmaps cleared).
    """
    table = f'Alignment_{db_nr}'
    limits = interface_package.interface_annotation_limits(config)
    results = binningSelect(prot_db_ids, ['Protein', 'Structure', 'Sequence_Identity'], table, config, locked=True)

    seq_ids_of_protein: dict[int, dict[int, float | None]] = {}
    for prot_db_id, structure_db_id, seq_id in results:
        if prot_db_id not in seq_ids_of_protein:
            seq_ids_of_protein[prot_db_id] = {}
        seq_ids_of_protein[prot_db_id][structure_db_id] = seq_id
    del results

    n_pruned_proteins = 0
    n_cleared = 0
    for prot_db_id, seq_id_of_structure in seq_ids_of_protein.items():
        n = len(seq_id_of_structure)
        struct_db_ids = np.fromiter(seq_id_of_structure.keys(), dtype=np.int64, count=n)
        seq_ids = np.fromiter(
            (np.nan if seq_id is None else seq_id for seq_id in seq_id_of_structure.values()),
            dtype=np.float64, count=n)
        if len(multi_chain_struct_ids) > 0:
            lookup = np.minimum(np.searchsorted(multi_chain_struct_ids, struct_db_ids),
                                len(multi_chain_struct_ids) - 1)
            multi_chain = multi_chain_struct_ids[lookup] == struct_db_ids
        else:
            multi_chain = np.zeros(n, dtype=bool)

        selected = interface_package.select_interface_annotations(
            struct_db_ids, seq_ids, multi_chain, *limits)

        # Single chain ones were written without a backmap in the first place.
        to_clear = struct_db_ids[multi_chain & ~selected].tolist()
        if len(to_clear) == 0:
            continue
        n_pruned_proteins += 1
        n_cleared += len(to_clear)
        for left in range(0, len(to_clear), _PRUNE_IN_CHUNK):
            set_null(config, table, 'Backmap',
                     equals_rows={'Protein': prot_db_id},
                     in_rows={'Structure': to_clear[left:left + _PRUNE_IN_CHUNK]},
                     locked=True)

    return len(seq_ids_of_protein), n_pruned_proteins, n_cleared


def prune_interface_backmaps(
        config: Config,
        prot_db_ids: list[int],
        multi_chain_struct_ids: np.ndarray) -> list:
    """Clears the stored backmaps that select_interface_annotations leaves out.

    align() writes a backmap for every annotation in a complex with more than one
    protein chain. Which of a protein's annotations keep theirs depends on the
    sequence identity ranking of all of them, and one protein's annotations can be
    spread over several alignment tasks, so this runs once they have all been
    written. prot_db_ids are the proteins with enough multi chain annotations to be
    cut at all; multi_chain_struct_ids is sorted.
    """
    times = []
    t0 = time.time()

    shards: dict[int, list[int]] = {}
    for prot_db_id in prot_db_ids:
        db_nr = prot_db_id % 20
        if db_nr not in shards:
            shards[db_nr] = []
        shards[db_nr].append(prot_db_id)

    if parallel_db_access(config):
        # Without this, ray.put falls into ray's auto_init_hook and starts an
        # unconfigured cluster. ray_init returns immediately if ray is already up.
        ray_init(config, stop_existing = False)
        conf_dump = ray.put(config)
        # Put once rather than serialised into each of the shard tasks.
        ids_store = ray.put(multi_chain_struct_ids)
        shard_results = ray.get([
            remote_prune_backmaps_shard.remote(db_nr, shard_ids, ids_store, conf_dump)
            for db_nr, shard_ids in shards.items()])
    else:
        shard_results = [prune_backmaps_shard(db_nr, shard_ids, multi_chain_struct_ids, config)
                         for db_nr, shard_ids in shards.items()]

    n_proteins = sum(result[0] for result in shard_results)
    n_pruned_proteins = sum(result[1] for result in shard_results)
    n_cleared = sum(result[2] for result in shard_results)

    t0 = add_to_times(times, t0, f'prune backmaps {n_proteins=} {n_pruned_proteins=} {n_cleared=}')

    if config.verbosity >= 2:
        config.logger.info(
            f'prune_interface_backmaps: {len(prot_db_ids)} candidate proteins, '
            f'{n_proteins} with alignments, {n_pruned_proteins} thinned out, '
            f'{n_cleared} backmaps cleared')

    return times


def get_interfaces(protein_db_ids, config: Config, proteins):
    if config.verbosity >= 3:    
        config.logger.info(f'Call of get_interfaces: {len(protein_db_ids)=}')
        t0 = time.time()

    rows = ['Protein', 'Interface_Id', 'Structure_Recommendation']
    table = 'Interface'
    results = binningSelect(protein_db_ids, rows, table, config)

    if config.verbosity >= 3:
        t1 = time.time()
        config.logger.info(f'Part 1 of get_interfaces: {t1-t0}')

    max_position_interface_select = 50000
    interface_map = {}
    interface_packages = []
    interface_key_package = []
    for row in results:
        protein_db_id = row[0]
        interface_id = row[1]
        interface_map[interface_id] = (protein_db_id, row[2])
        interface_key_package.append(interface_id)
        if len(interface_key_package) >= max_position_interface_select:
            interface_packages.append(interface_key_package)
            interface_key_package = []

    if config.verbosity >= 3:
        t2 = time.time()
        config.logger.info(f'Part 2 of get_interfaces: {t2-t1}')
    if len(interface_key_package) > 0: 
        interface_packages.append(interface_key_package)

    results = []
    for interface_key_package in interface_packages:
        results += binningSelect(interface_key_package, ['Interface', ' Position', 'Recommended_Residue'], 'RS_Position_Interface', config, full_select=True)

    if config.verbosity >= 3:
        t3 = time.time()
        config.logger.info(f'Part 3 of get_interfaces: {t3-t2}')
    interfaces = {}
    max_size = 500000
    position_packages = []
    positions = []
    for row in results:
        interface_db_id = row[0]
        pos_db_id = row[1]
        if interface_db_id not in interfaces:
            interfaces[interface_db_id] = {}
        interfaces[interface_db_id][pos_db_id] = row[2]
        positions.append(pos_db_id)
        if len(positions) >= max_size:
            position_packages.append(positions)
            positions = []

    if len(positions) > 0:
        position_packages.append(positions)

    if config.verbosity >= 3:
        t4 = time.time()
        config.logger.info(f'Part 4 of get_interfaces: {t4-t3}')
    pos_id_map = {}
    for positions in position_packages:
        if config.verbosity >= 3:
            config.logger.info(f'get_interfaces Position select: {len(positions)}')
        results = binningSelect(positions, ['Position_Id', 'Position_Number'], 'Position', config)

        for row in results:
            pos_id_map[row[0]] = row[1]

    if config.verbosity >= 3:
        t5 = time.time()
        config.logger.info(f'Part 5 of get_interfaces: {t5-t4}')
    protein_interface_map = {}
    for interface_db_id in interfaces:
        protein_db_id, structure_recommendation = interface_map[interface_db_id]
        recommended_complex, chain, interacting_chain = structure_recommendation.split(',')
        protein_id = proteins.getByDbId(protein_db_id).primary_protein_id
        positions = {}
        for pos_db_id in interfaces[interface_db_id]:
            structure_tuple, res = interfaces[interface_db_id][pos_db_id].split()
            structure_id, chain = structure_tuple.split(':')
            positions[(pos_id_map[pos_db_id])] = structure_id, chain, res
        interface_obj = interface_package.Aggregated_interface(protein_id, recommended_complex = recommended_complex, chain = chain, interacting_chain = interacting_chain, positions = positions, database_id = interface_db_id)
        if protein_id not in protein_interface_map:
            protein_interface_map[protein_id] = []
        protein_interface_map[protein_id].append(interface_obj)

    if config.verbosity >= 3:
        t6 = time.time()
        config.logger.info(f'Part 6 of get_interfaces: {t6-t5}')
    for protein_id in protein_interface_map:
        if config.verbosity >= 6:
            config.logger.info(f'Retrieved interfaces from database for {protein_id}:')
            for interface in protein_interface_map[protein_id]:
                interface.print_interface()
        proteins.set_aggregated_interface_map(protein_id, protein_interface_map[protein_id])

    del protein_interface_map

    if config.verbosity >= 3:
        t7 = time.time()
        config.logger.info(f'Part 7 of get_interfaces: {t7-t6}')

def prepare_alignments(proteins: protein_package.Proteins, config: Config):

    if config.verbosity >= 3:
        config.logger.info(f'Call fo prepare_Alignments: {len(proteins.stored_ids)=}')

    structureCheck(proteins, config)

    stored_structure_ids = proteins.getStoredStructureIds()

    stored_complexes = set()
    for structure_db_id in stored_structure_ids:
        structure_id = stored_structure_ids[structure_db_id][0]
        stored_complexes.add(structure_id)

    draw_complexes(config, proteins, stored_complexes=stored_complexes)

    
# called by serializedPipeline
def getAlignments(proteins, config: Config, get_all_alignments=False):
    if config.verbosity >= 2:
        t0 = time.time()

    if config.verbosity >= 3:
        config.logger.info(f'Call fo getAlignments: {get_all_alignments=} {len(proteins.stored_ids)=}')

    prot_db_ids = proteins.get_stored_ids(exclude_completely_stored=(not get_all_alignments))

    for prot_id in proteins.indels:
        if proteins[prot_id].stored:
            prot_db_ids.add(proteins[prot_id].database_id)

    rows = ['Protein', 'Structure', 'Sequence_Identity', 'Coverage', 'Alignment']
    table = 'Alignment'
    results = binningSelect(prot_db_ids, rows, table, config)

    if config.verbosity >= 2:
        t1 = time.time()
        config.logger.info(f'Time for part 1 in getAlignments: {t1 - t0}, number of stored proteins {len(prot_db_ids)}, number of stored alignments: {len(results)}')

    prot_structure_alignment_map = {}
    structure_ids = set()

    for row in results:

        prot_db_id = row[0]

        structure_id = row[1]
        seq_id = row[2]
        coverage = row[3]
        packed_alignment = row[4]

        structure_ids.add(structure_id)

        #target_seq, template_seq = sdsc_utils.process_alignment_data(alignment)

        if prot_db_id not in prot_structure_alignment_map:
            prot_structure_alignment_map[prot_db_id] = {}

        prot_structure_alignment_map[prot_db_id][structure_id] = (packed_alignment, coverage, seq_id)

    if config.verbosity >= 2:
        t2 = time.time()
        config.logger.info(f"Time for part 2 in getAlignments: {t2 - t1}")

    structureCheck(proteins, config)

    if config.verbosity >= 2:
        t3 = time.time()
        config.logger.info(f"Time for part 3 in getAlignments: {t3 - t2}")

    stored_structure_ids = proteins.getStoredStructureIds()

    if config.verbosity >= 2:
        t4 = time.time()
        config.logger.info(f"Time for part 4 in getAlignments: {t4 - t3}")

    structure_target_chain_map = {}
    stored_complexes = set()
    for structure_db_id in stored_structure_ids:
        (structure_id, chain) = stored_structure_ids[structure_db_id]
        structure_target_chain_map[structure_id] = chain
        stored_complexes.add(structure_id)

    draw_complexes(config, proteins, stored_complexes=stored_complexes)

    results = retrieve_structure_info_from_struct_ids(config, structure_target_chain_map.keys())

    interacting_structures = {}
    for row in results:
        structure_id = row[0]
        structure_db_id = row[1]
        chain = row[2]
        if structure_id in structure_target_chain_map:
            if chain != structure_target_chain_map[structure_id]:
                interacting_structures[structure_db_id] = (structure_id, chain)

    if config.verbosity >= 2:
        t5 = time.time()
        config.logger.info(f"Time for part 5 in getAlignments: {t5 - t4}")

    interacting_protein_db_ids = []

    if config.compute_ppi:
        rows = ['Structure', 'Protein', 'Sequence_Identity', 'Coverage', 'Alignment']
        table = 'Alignment'
        results = binningSelect(interacting_structures.keys(), rows, table, config)

        for row in results:

            prot_db_id = row[1]

            structure_id = row[0]
            seq_id = row[2]
            coverage = row[3]
            #alignment = unpack(row[4])
            packed_alignment = row[4]

            structure_ids.add(structure_id)

            #target_seq, template_seq = sdsc_utils.process_alignment_data(alignment)

            if prot_db_id not in prot_structure_alignment_map:
                prot_structure_alignment_map[prot_db_id] = {}
                interacting_protein_db_ids.append(prot_db_id)

            prot_structure_alignment_map[prot_db_id][structure_id] = (packed_alignment, coverage, seq_id)

    if config.verbosity >= 2:
        t6 = time.time()
        config.logger.info(f"Time for part 6 in getAlignments: {t6 - t5} {len(interacting_structures)=}")

    if config.compute_ppi:
        #pos_db_map = retrieve_stored_proteins(interacting_protein_db_ids, config, proteins, without_positions=True)
        _ = retrieve_stored_proteins(interacting_protein_db_ids, config, proteins, with_struct_recs=False, return_pos_db_map=False)

    if config.verbosity >= 2:
        t7 = time.time()
        config.logger.info(f'Time for part 7 in getAlignments: {str(t7 - t6)}, {len(interacting_protein_db_ids)}')

    #if config.compute_ppi:
    #    get_interfaces(interacting_protein_db_ids, config, proteins)

    if config.verbosity >= 2:
        t8 = time.time()
        config.logger.info(f"Time for part 8 in getAlignments: {t8 - t7}")

    structure_map: dict[str, dict[str, tuple[int, str]]]
    structure_map, id_structure_map = getStructure_map(structure_ids, config)

    if config.verbosity >= 2:
        t9 = time.time()
        config.logger.info(f"Time for part 9 in getAlignments: {t9 - t8}")

    #pdb_ids = set()
    #for pdb_id in structure_map:
    #    pdb_ids.add(pdb_id)

    complex_map = getComplexMap(config, pdb_ids=proteins.complexes.keys())

    if config.verbosity >= 2:
        t10 = time.time()
        config.logger.info(f"Time for part 10 in getAlignments: {t10 - t9} {len(prot_structure_alignment_map)=}")

    for prot_db_id in prot_structure_alignment_map:
        protein_obj = proteins.getByDbId(prot_db_id)
        prot_id = protein_obj.primary_protein_id
        for structure_id in prot_structure_alignment_map[prot_db_id]:
            (pdb_id, chain) = id_structure_map[structure_id]

            (packed_alignment, coverage, seq_id) = prot_structure_alignment_map[prot_db_id][structure_id]

            struct_anno = structure_package.StructureAnnotation(prot_id, pdb_id, chain, alignment=packed_alignment, stored=True)
            protein_obj.add_annotation(struct_anno)

            if not proteins.contains_structure(pdb_id, chain):
                oligo: str = structure_map[pdb_id][chain][1]
                struct = structure_package.Structure(pdb_id, chain, oligo=oligo, mapped_proteins=[prot_id], database_id=structure_id)
                if structure_id in interacting_structures:
                    struct.interacting_structure = True
                proteins.add_structure(pdb_id, chain, struct)
            else:
                proteins.add_mapping_to_structure(pdb_id, chain, prot_id)

            if not proteins.contains_complex(pdb_id) and (pdb_id in complex_map):
                (comp_id, resolution, chains_str, homooligomers, lig_profile, metal_profile, ion_profile, cc_profile) = complex_map[pdb_id]
                compl = complex_package.Complex(pdb_id, resolution=resolution, chains_str=chains_str, lig_profile=lig_profile,
                                     metal_profile=metal_profile, ion_profile=ion_profile,
                                     chain_chain_profile=cc_profile, stored=True, database_id=comp_id, homomers_str=homooligomers)
                proteins.add_complex(pdb_id, compl)
            elif (pdb_id in complex_map) and not proteins.complexes[pdb_id].stored:
                proteins.complexes[pdb_id].stored = True
                proteins.complexes[pdb_id].daatabes_id = complex_map[pdb_id][0]

            proteins.set_structure_db_id(pdb_id, chain, structure_id)

            # Set on the object just built rather than looked back up by
            # (pdb_id, chain): the three setters did one lookup each, and reaching
            # annotations by key now costs building the protein's annotation index.
            struct_anno.set_coverage(coverage)
            protein_obj.number_of_mappings += int(len(protein_obj.positions) * coverage)
            struct_anno.set_sequence_id(seq_id)
            struct_anno.set_database_id(True)

            if is_alphafold_model(pdb_id) and config.provide_msas:
                if config.verbosity >= 3:
                    config.logger.info(f'fetching msa for {prot_id}  {pdb_id}')
                seq = proteins[prot_id].sequence
                fetch_msa(config, pdb_id, prot_id, seq, seq_id, config.proc_n)

    if config.verbosity >= 2:
        t11 = time.time()
        config.logger.info(f"Time for part 11 in getAlignments: {t11 - t10}")


def getStructure_map(structure_ids, config):
    structure_map: dict[str, dict[str, tuple[int, str]]] = {}
    id_structure_map = {}
    if len(structure_ids) > 0:
        results = binningSelect(structure_ids, ['Structure_Id', 'PDB', 'Chain', 'Homooligomer'], 'Structure', config)

        for row in results:
            s_id = row[0]

            pdb_id = row[1]
            chain = row[2]
            oligo: str = row[3]
            
            if pdb_id not in structure_map:
                structure_map[pdb_id] = {}
            structure_map[pdb_id][chain] = (s_id, oligo)
            id_structure_map[s_id] = (pdb_id, chain)
    return structure_map, id_structure_map


def checkLigand(name, db, cursor):
    sql = "SELECT Ligand_Id FROM Ligand WHERE Name = '%s'" % name
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in checkLigand: %s" % sql)
        
    if results == ():
        return 0
    elif len(results) == 0:
        return 0
    else:
        row = results[0]
        return row[0]


def majority_vote(secs):
    class_dict = {}
    for (sec, qual) in secs:
        if sec not in class_dict:
            class_dict[sec] = qual
        else:
            class_dict[sec] += qual

    max_sec = None
    max_qual = 0.0
    for sec in class_dict:
        qual = class_dict[sec]
        if qual > max_qual:
            max_qual = qual
            max_sec = sec
    return max_sec


def getChemicalDistance(aac):
    try:
        if aac.count(',') < 1:
            aa1 = aac[0]
            aa2 = aac[-1]

            chemical_distance = residue_consts.CHEM_DIST_MATRIX[aa1][aa2]
        else:
            aa1 = aac[0]
            aa2s = aac.split(",")
            aa2s[0] = aa2s[0][-1]
            chem_dists = []
            for aa2 in aa2s:
                chem_dists.append(residue_consts.CHEM_DIST_MATRIX[aa1][aa2])
            chemical_distance = float(sum(chem_dists)) / float(len(chem_dists))
    except:
        return None
    return chemical_distance


def getBlosumValue(aac):
    if aac.count(',') < 1:
        try:
            try:
                blosum_value = residue_consts.BLOSUM62[(aac[0], aac[-1])]
            except:
                blosum_value = residue_consts.BLOSUM62[(aac[-1], aac[0])]
        except:
            blosum_value = 0.0
    else:
        aa1 = aac[0]
        aa2s = aac.split(",")
        aa2s[0] = aa2s[0][-1]
        bvs = []
        for aa2 in aa2s:
            try:
                try:
                    blosum_value = residue_consts.BLOSUM62[(aa1, aa2)]
                except:
                    blosum_value = residue_consts.BLOSUM62[(aa2, aa1)]
            except:
                blosum_value = 0.0
            bvs.append(blosum_value)
        blosum_value = float(sum(bvs)) / float(len(bvs))
    return blosum_value


'''
#called by babel
def getSLD(mutation_id,template_id,db,cursor):
    sql = "SELECT Sub_Lig_Dist FROM RS_Mutation_Template WHERE Template = '%s' AND Mutation = '%s'" % (str(template_id),str(mutation_id))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getSLD")
        # Rollback in case there is any NameError
        db.rollback()
    return results[0][0]
'''

'''
#called by babel
#structure of chosen_ones: {ligand_three_letter_code:tanimoto_score}
def getLigandAnnotation(chosen_ones,session_id,distance_threshold,db,cursor):
    anno_dict = {}

    ligand_map = {}
    ligand_ids = set()
    sql = "SELECT Ligand_Id,Name FROM Ligand"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in getLigandAnnotation: %s\n%s" % (sql, f))
    for row in results:
        ligand_name = row[1]
        if ligand_name not in chosen_ones:
            continue
        ligand_map[ligand_name] = row[0]
        ligand_ids.add(row[0])

    sql = "SELECT Template,Ligand FROM RS_Ligand_Template"
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getLigandAnnotation: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()

    ligand_template_map = {}
    template_ids = set()
    for row in results:
        if not row[1] in ligand_ids:
            continue
        template_ids.add(row[0])
        if not row[1] in ligand_template_map:
            ligand_template_map[row[1]] = set([row[0]])
        else:
            ligand_template_map[row[1]].add(row[0])

    sql = "SELECT Mutation,Template FROM RS_Annotation_Session WHERE Session = '%s'" % str(session_id)
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getLigandAnnotation: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()

    template_mutation_map = {}
    for row in results:
        if not row[1] in template_ids:
            continue
        if not row[1] in template_mutation_map:
            template_mutation_map[row[1]] = {}

        template_mutation_map[row[1]][row[0]] = None

    min_t_id = min(template_ids)
    max_t_id = max(template_ids)

    t = 50000
    old_max = None
    if max_t_id - min_t_id > t:
        old_max = max_t_id
        max_t_id = min_t_id + t

    while True:
        sql = "SELECT Mutation,Template,Sub_Lig_Dist FROM RS_Mutation_Template WHERE Template BETWEEN %s AND %s" % (str(min_t_id), str(max_t_id))
        try:
            # Execute the SQL command
            cursor.execute(sql)
            results = cursor.fetchall()
            # Commit your changes in the database
            db.commit()
        except:
            raise NameError("Error in getLigandAnnotation: %s" % sql)
            # Rollback in case there is any NameError
            db.rollback()

        for row in results:
            if not row[1] in template_mutation_map:
                continue
            if not row[0] in template_mutation_map[row[1]]:
                continue
            template_mutation_map[row[1]][row[0]] = row[2]

        if old_max is None:
            break
        else:
            min_t_id = max_t_id + 1
            max_t_id = old_max
            if max_t_id - min_t_id > t:
                old_max = max_t_id
                max_t_id = min_t_id + t
            else:
                old_max = None

    for lig in chosen_ones:
        if lig not in ligand_map:
            continue
        lig_id = ligand_map[lig]
        filtered_annotation_ids = {}
        if lig_id not in ligand_template_map:
            continue
        for template_id in ligand_template_map[lig_id]:
            if template_id not in template_mutation_map:
                continue
            for mutation_id in template_mutation_map[template_id]:

                sub_lig_distances = template_mutation_map[template_id][mutation_id]
                if sub_lig_distances is not None:
                    sub_lig_distances = sub_lig_distances.split(",")
                    good_dists = []
                    for sub_lig_distance in sub_lig_distances:
                        if not sub_lig_distance == "-":
                            infos = sub_lig_distance.split(":")
                            dist = float(infos[1])
                            lig_infos = infos[0].split("_")
                            lig_name = lig_infos[0]
                            if lig_name == lig and dist <= float(distance_threshold):
                                good_dists.append(sub_lig_distance)
                    if len(good_dists) > 0:
                        good_dists = sorted(good_dists, key=lambda x: float(x.rsplit(":")[1]))
                        if template_id in filtered_annotation_ids:
                            filtered_annotation_ids[template_id].append((mutation_id, good_dists))
                        else:
                            filtered_annotation_ids[template_id] = [(mutation_id, good_dists)]
        if len(filtered_annotation_ids) > 0:
            for template_id in filtered_annotation_ids:
                filtered_annotation_ids[template_id] = sorted(filtered_annotation_ids[template_id], key=lambda x: float(x[1][0].split(":")[1]))
            anno_dict[lig] = (chosen_ones[lig], filtered_annotation_ids)
    return anno_dict
'''


def getProtIdsFromSession(session_id, config: Config, prot_db_ids = None, filter_mutant_proteins=False) -> dict[int, tuple[str, None | set[str]]] :
    table = 'RS_Protein_Session'
    cols = ['Protein', 'Input_Id', 'Tags']
    eq_cols = {'Session': session_id}
    results = select(config, cols, table, equals_rows=eq_cols)
    prot_ids: dict[int, tuple[str, None | set[str]]] = {}
    for row in results:
        if config.verbosity >= 5:
            config.logger.info(f'In getProtIdsFromSession: {row=}')
        if filter_mutant_proteins and row[1] is None:
            continue
        if prot_db_ids is not None:
            if row[0] not in prot_db_ids:
                continue
        tags = row[2]
        if tags is not None:
            tags = set(tags.split(','))
        prot_ids[row[0]] = (row[1], tags)
    return prot_ids

def retrieve_gene_id_map(gene_db_ids, config):
    table = 'Gene'
    columns = ['Gene_Id', 'Primary_Gene_Id', 'Gene_Name']
    results = binningSelect(gene_db_ids, columns, table, config)

    gene_id_map = {}
    for row in results:
        gene_id_map[row[0]] = (row[1], row[2])
    return gene_id_map

def retrieve_stored_proteins(
        prot_db_ids: dict[int, tuple[str, None | set[str]]],
        config: Config,
        proteins,
        with_mappings = False,
        with_struct_recs = True,
        without_positions = False,
        return_pos_db_map = True,
        input_id_is_majored = False):
    if config.verbosity >= 2:
        t0 = time.time()
    cols = ['Protein_Id', 'Primary_Protein_Id', 'Sequence', 'Mutant_Type', 'Mutant_Positions', 'Gene']
    if not isinstance(prot_db_ids, dict):
        protein_ids_unknown = True
    else:
        protein_ids_unknown = False

    if protein_ids_unknown:
        results = binningSelect(prot_db_ids, cols, 'Protein', config)
    else:
        results = binningSelect(prot_db_ids.keys(), cols, 'Protein', config)

    if config.verbosity >= 2:
        t1 = time.time()
        config.logger.info(f'Retrieve stored proteins, part 1: {t1 - t0} {len(prot_db_ids)=} {input_id_is_majored=}')

    id_prot_id_map = {}

    prot_id_list = set()
    prot_ids_mutants_excluded = set()
    prot_gene_db_id_map = {}
    gene_db_ids = set()

    for row in results:
        protein_db_id: int = row[0]
        if config.verbosity >= 5:
            config.logger.info(f'{row=}')
        if protein_ids_unknown:
            input_id = row[1]
            tags = None
        else:
            input_id, tags = prot_db_ids[protein_db_id]
        mut_type = row[3]
        if row[4] is None:
            sav_positions = None
            insertion_positions = None
            deletion_flanks = None
        else:
            sav_string, insertion_string, deletion_string = row[4].split(';')
            if sav_string == '':
                sav_positions = None
            else:
                sav_positions = [int(x) for x in sav_string.split(',')]

            insertion_positions = []
            for ins_str in insertion_string.split('!'):
                if ins_str != '':
                    insertion_positions.append([int(x) for x in ins_str.split(',')])
            if len(insertion_positions) == 0:
                insertion_positions = None

            deletion_flanks = []
            for del_str in deletion_string.split('!'):
                if del_str != '':
                    lfs, rfs = del_str.split(',')
                    if lfs == '_':
                        lf = None
                    else:
                        lf = int(lfs)
                    if rfs == '_':
                        rf = None
                    else:
                        rf = int(rfs)
                    deletion_flanks.append((lf, rf))
            if len(deletion_flanks) == 0:
                deletion_flanks = None

        if mut_type == 'SAV':
            wt_prot = row[1]
        else:
            wt_prot = None

        gene_db_id = row[5]
        
        if input_id_is_majored:
            primary_id = input_id
        else:
            primary_id = row[1]

        id_prot_id_map[protein_db_id] = primary_id

        if gene_db_id is not None:
            prot_gene_db_id_map[primary_id] = gene_db_id
            gene_db_ids.add(gene_db_id)

        prot_obj = protein_package.Protein(config.errorlog, primary_protein_id=primary_id, database_id=protein_db_id, input_id=input_id, sequence=row[2], wildtype_protein = wt_prot, mutant_type = mut_type, sav_positions = sav_positions, insertion_positions = insertion_positions, deletion_flanks = deletion_flanks, tags = tags)
        proteins[primary_id] = prot_obj
        prot_id_list.add(protein_db_id)
        if mut_type is not None:
            prot_ids_mutants_excluded.add(protein_db_id)

    if config.verbosity >= 2:
        t2 = time.time()
        config.logger.info(f'Retrieve stored proteins, part 2: {t2 - t1}')

    proteins.set_stored_ids(prot_id_list, prot_ids_mutants_excluded)

    if config.verbosity >= 2:
        t3 = time.time()
        config.logger.info(f'Retrieve stored proteins, part 3: {t3 - t2}')

    gene_id_map = retrieve_gene_id_map(gene_db_ids, config)

    if config.verbosity >= 2:
        t4 = time.time()
        config.logger.info(f'Retrieve stored proteins, part 4: {t4 - t3}')

    for prot_id in prot_gene_db_id_map:
        proteins[prot_id].gene = gene_id_map[prot_gene_db_id_map[prot_id]][0]

    if config.verbosity >= 2:
        t5 = time.time()
        config.logger.info(f'Retrieve stored proteins, part 5: {t5 - t4} {without_positions=} {with_mappings=} {with_struct_recs=}')

    pos_db_map = {}

    error_count = 0
    error_message = None

    if not without_positions:
        if with_mappings:
            cols = ['Protein', 'Position_Number', 'Position_Id', 'Wildtype_Residue', 'Recommended_Structure_Data', 'Position_Data']
        elif with_struct_recs:
            cols = ['Protein', 'Position_Number', 'Position_Id', 'Wildtype_Residue', 'Recommended_Structure_Data']
        else:
            cols = ['Protein', 'Position_Number', 'Position_Id', 'Wildtype_Residue']
        table = 'Position'

        if protein_ids_unknown:
            results = binningSelect(prot_db_ids, cols, table, config)
        else:
            results = binningSelect(prot_db_ids.keys(), cols, table, config)

        if config.verbosity >= 2:
            t6 = time.time()
            config.logger.info(f'Retrieve stored proteins, part 6: {t6 - t5} {len(results)}')

        for row in results:
            p_id = row[0]

            m_id = row[2]
            pos = row[1]

            wt_aa = row[3]

            if with_struct_recs:
                if row[4] is not None:
                    (recommended_structure_str, max_seq_structure_str) = unpack(row[4])
                else:
                    recommended_structure_str = None
                    max_seq_structure_str = None

                recommended_structure, seq_id, cov, resolution = sdsc_utils.process_recommend_structure_str(recommended_structure_str)
                max_seq_structure, max_seq_seq_id, max_seq_cov, max_seq_resolution = sdsc_utils.process_recommend_structure_str(max_seq_structure_str)



            if with_mappings:

                if row[5] is None:
                    mappings_obj = mappings_package.Mappings()

                else:
                    try:
                        packed_mappings, _ = unpack(row[5])
                        mappings_obj = unpack(packed_mappings)
                        
                        mappings_obj.recommended_res = recommended_structure_str
                        mappings_obj.max_seq_res = max_seq_structure_str
                        
                    except:
                        error_count += 1
                        if error_message is None:
                            [e, f, g] = sys.exc_info()
                            g = traceback.format_exc()
                            error_message = f'Error in retrieving position data, couldnt unpack: {m_id}\n{e}\n{f}\n{g}'
                        continue

            if with_struct_recs:
                pos_obj = position_package.Position(pos=pos, checked=True, recommended_structure=recommended_structure, database_id=m_id, wt_aa = wt_aa)
            else:
                pos_obj = position_package.Position(pos=pos, checked=True, database_id=m_id, wt_aa = wt_aa)

            if with_mappings:
                pos_obj.mappings = mappings_obj
            prot_id = id_prot_id_map[p_id]
            proteins[prot_id].add_positions([pos_obj])
            if return_pos_db_map:
                pos_db_map[m_id] = (prot_id, pos)
    else:
        t6 = time.time()

    if config.verbosity >= 2:
        t7 = time.time()
        config.logger.info(f'Retrieve stored proteins, part 7: {t7 - t6}')

    if error_count > 0:
        config.logger.info(f'There were {error_count} number of errors retrieving position data, example error:\n{error_message}')
    return pos_db_map

# called by output

def proteinsFromDb(session, config: Config, with_residues=False, filter_mutant_proteins=False, with_mappings = False,
                   with_snvs=False, with_indels = False, with_multi_mutations = False, mutate_snvs=False, get_stable_objs=False,
                   with_alignments=False, with_complexes=False, keep_ray_alive=False, current_chunk = None, prot_db_id_dict = None,
                   input_id_is_majored = False, update_complexes = False):

    if config.verbosity >= 2:
        config.logger.info(f'Calling proteinsFromDb, {session=} {with_residues=} {with_mappings=} {input_id_is_majored=} {filter_mutant_proteins=}')

    if with_alignments:
        with_complexes = True

    prot_db_id_dict: dict[int, tuple[str, None | set[str]]] = getProtIdsFromSession(session, config, filter_mutant_proteins=filter_mutant_proteins, prot_db_ids = prot_db_id_dict)

    if config.verbosity >= 5:
        config.logger.info(f'In proteinsFromDb, {prot_db_id_dict=}')

    if current_chunk is not None:
        n_of_chunks = (len(prot_db_id_dict) // config.chunksize) + 1
        small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(n_of_chunks, len(prot_db_id_dict))

        if current_chunk <= n_of_big_chunks:
            l = (current_chunk - 1) * big_chunksize
            r = current_chunk * big_chunksize
        elif current_chunk <= n_of_big_chunks + n_of_small_chunks:
            l = (((current_chunk - 1) - n_of_big_chunks) * small_chunksize) + (big_chunksize * n_of_big_chunks)
            r = ((current_chunk - n_of_big_chunks) * small_chunksize) + (big_chunksize * n_of_big_chunks)
        else:
            return None, None

        try:
            prot_db_ids = list(prot_db_id_dict.keys())
            chunked_prot_db_id_dict = {}
            for prot_db_id in prot_db_ids[l:r]:
                chunked_prot_db_id_dict[prot_db_id] = prot_db_id_dict[prot_db_id]
            complete_prot_db_id_dict = prot_db_id_dict
            prot_db_id_dict = chunked_prot_db_id_dict
        except:
            prot_db_ids = list(prot_db_id_dict)
            chunked_prot_db_id_dict = []
            for prot_db_id in prot_db_ids[l:r]:
                chunked_prot_db_id_dict.append(prot_db_id)
            complete_prot_db_id_dict = prot_db_id_dict
            prot_db_id_dict = chunked_prot_db_id_dict

    proteins = protein_package.Proteins({}, {}, {})  # create empty Proteins object

    pos_db_map = retrieve_stored_proteins(prot_db_id_dict, config, proteins, with_mappings = with_mappings, input_id_is_majored = input_id_is_majored)

    if config.verbosity >= 5:
        config.logger.info(f'In proteinsFromDb, point 0 - {list(proteins.protein_map.keys())=}')


    if with_complexes:
        draw_complexes(config, proteins, draw_all=True, get_stable_objs=get_stable_objs, update_infos=update_complexes)

    if with_multi_mutations:
        with_snvs = True
        with_indels = True

    if with_snvs:
        cols = ['Position', 'New_AA', 'SNV_Id']
        table = 'SNV'

        results = binningSelect(pos_db_map.keys(), cols, table, config, full_select=True)

        snv_db_id_map = {}
        for row in results:
            snv = snv_package.SNV(row[1])
            (prot_id, pos) = pos_db_map[row[0]]
            proteins[prot_id].positions[pos].mut_aas[row[1]] = snv
            snv_db_id_map[row[2]] = (pos, snv)

    if config.verbosity >= 5:
        config.logger.info(f'In proteinsFromDb, point A - {list(proteins.protein_map.keys())=}')

    if with_indels:

        table = 'RS_Indel_Session'
        rows = ['Indel', 'Tags']
        eq_rows = {'Session': session}

        results = select(config, rows, table, equals_rows=eq_rows)
        indel_tags_map = {}
        for row in results:
            indel_tags_map[row[0]] = row[1]

        table = 'Indel'
        columns = ['Wildtype_Protein', 'Indel_Id', 'Indel_Notation', 'Mutant_Protein']
        try:
            prot_db_ids = prot_db_id_dict.keys()
        except:
            prot_db_ids = list(prot_db_id_dict)
            
        results = binningSelect(prot_db_ids, columns, table, config)

        indel_pre_map = {}
        indel_db_id_map = {}
        for row in results:
            mut_prot_db_id = row[3]

            indel_db_id = row[1]
            wt_prot_db_id = row[0]
            if wt_prot_db_id not in prot_db_id_dict:
                continue

            wt_prot = proteins.getByDbId(wt_prot_db_id).primary_protein_id

            indel_notation = row[2]

            if indel_db_id not in indel_tags_map:
                continue

            tags = indel_tags_map[indel_db_id]
            if isinstance(tags, str):
                tags = set(tags.split(','))

            if indel_notation.count('delins') == 1:
                indel = indel_package.Substitution(raw_str=indel_notation, tags=tags)

            elif indel_notation.count('del') == 1:
                indel = indel_package.Deletion(raw_str=indel_notation, tags=tags)

            elif indel_notation.count('ins') == 1:
                indel = indel_package.Insertion(raw_str=indel_notation, tags=tags)

            indel.database_id = indel_db_id

            indel_db_id_map[indel_db_id] = indel

            #if mut_prot_db_id is None:
            #    #Indels without mut_prot are part of a multi-mutation
            #    continue

            if wt_prot not in indel_pre_map:
                indel_pre_map[wt_prot] = []
            indel_pre_map[wt_prot].append(indel)

        indel_map = {}
        for wt_prot in indel_pre_map:
            uniprot.indel_insert(config, proteins, indel_map, indel_pre_map[wt_prot], wt_prot, proteins_is_object = True)

    

    if with_multi_mutations:
        table = 'RS_Multi_Mutation_Session'
        rows = ['Multi_Mutation', 'Tags']
        eq_rows = {'Session': session}

        results = select(config, rows, table, equals_rows=eq_rows)
        mm_tags_map = {}
        for row in results:
            mm_tags_map[row[0]] = row[1]

        table = 'Multi_Mutation'
        columns = ['Multi_Mutation_Id', 'SNVs', 'Indels', 'Wildtype_Protein', 'Mutant_Protein']
        results = binningSelect(mm_tags_map.keys(), columns, table, config)

        mm_db_map = {}
        for row in results:
            mm_db_id = row[0]
            snvs = row[1]
            indels = row[2]
            wt_prot_db_id = row[3]

            if wt_prot_db_id not in prot_db_id_dict:
                continue

            mm_tags = mm_tags_map[mm_db_id]

            mut_prot_db_id = row[4]

            wt_prot = proteins.getByDbId(wt_prot_db_id).primary_protein_id
            try:
                mut_prot = proteins.getByDbId(mut_prot_db_id).primary_protein_id
            except:
                mut_prot = '-'
            mut_list = []

            if snvs is not None:
                for snv_db_id in snvs.split(','):
                    try:
                        snv_db_id = int(snv_db_id)
                    except:
                        continue
                    pos, snv = snv_db_id_map[snv_db_id]
                    aa2 = snv.new_aa
                    position_obj = proteins[wt_prot].positions[pos]
                    mut_list.append((position_obj, aa2))
            if indels is not None:
                for indel_db_id in indels.split(','):
                    try:
                        indel_db_id = int(indel_db_id)
                    except:
                        continue
                    mut_list.append(indel_db_id_map[indel_db_id])

            if wt_prot not in indel_map:
                indel_map[wt_prot] = {}

            if config.verbosity >= 4:
                config.logger.info(f'In proteinsFromDB: Call of add_multi_mutation for {wt_prot} {mut_prot}')

            proteins[wt_prot].add_multi_mutation(mut_list, indel_map[wt_prot], mm_tags = mm_tags, mut_prot_id = mut_prot)

            if wt_prot not in mm_db_map:
                mm_db_map[wt_prot] = {}
            mm_db_map[wt_prot][mut_prot] = mm_db_id

        multi_mutation_objects = []
        for primary_protein_id in list(proteins.protein_map.keys()).copy():  # iterate over snapshot, since the mutant proteins get added to proteins in the loop
            if proteins[primary_protein_id].multi_mutations is not None:
                # creates mutation protein object
                multi_mutation_objects += proteins[primary_protein_id].create_multi_mutations(proteins, config)

        for multi_mutation_obj in multi_mutation_objects:
            if multi_mutation_obj.wt_prot not in proteins.multi_mutations:
                proteins.multi_mutations[multi_mutation_obj.wt_prot] = []
            proteins.multi_mutations[multi_mutation_obj.wt_prot].append(multi_mutation_obj)
            mm_db_id = mm_db_map[multi_mutation_obj.wt_prot][multi_mutation_obj.mut_prot]
            proteins.multi_mutation_back_map[mm_db_id] = multi_mutation_obj

    if config.verbosity >= 5:
        config.logger.info(f'In proteinsFromDb, point B - {list(proteins.protein_map.keys())=}')

    ray_utils.ray_init(config, stop_existing = False)

    if mutate_snvs:
        prot_ids = list(proteins.get_protein_ids()).copy()
        for prot_id in prot_ids:
            proteins[prot_id].mutate_snvs(proteins, config)

        serializedPipeline.autoTemplateSelection(config, proteins)

    if with_alignments:
        No_Errors = serializedPipeline.paraAlignment(config, proteins, skip_inserts=True, get_all_alignments=True)

    if not keep_ray_alive:
        ray.shutdown()

    if with_residues:
        getStoredResidues(proteins, config)

    if current_chunk is not None:
        return proteins, complete_prot_db_id_dict

    return proteins


# method for comparing/developing the new classifiction
def diffSurfs(mutation_surface_dict, g_u_dict, mutation_dict, outfile, config):
    class_dict = {}
    for m in mutation_surface_dict:
        (u_ac, u_id) = g_u_dict[mutation_dict[m][1]]
        aac = mutation_dict[m][0]
        DSC = 0.0  # decision sum core
        DSS = 0.0  # decision sum surface
        n = 0.0
        min_surf = 2.0
        lines = ["Coverage\trASA"]
        for (surface, qual, cov) in mutation_surface_dict[m]:
            if surface < config.surface_threshold:
                DSC += qual * (cov**5)
            else:
                DSS += qual * (cov**5)
            lines.append("%s\t%s" % (str(cov), str(surface)))
            n += 1.0

            if surface < min_surf:
                min_surf = surface

        weighted_surface_value = DSS - 2 * DSC
        if weighted_surface_value > 0:
            weighted_sc = "Surface"
            if min_surf < config.surface_threshold:
                f = open("%s.%s_%s.tsv" % (outfile, u_ac, aac), 'w')
                f.write("\n".join(lines))
                f.close()
        else:
            weighted_sc = "Core"

        conf_sc = (1.0 - 1.0 / (n + 1.0)) * abs(DSS - 2 * DSC) / (DSS + 2 * DSC)


def excludeFarClasses(c, sc):
    if c == "Surface" or c == "Core" or c == 'Disorder' or c is None:
        return c

    interactions = re.sub(r' Interaction$', '', re.sub(r'^[^:]*: ', '', c)).split(' and ')
    non_far_interactions = [x for x in interactions if ' far' not in x]

    if (len(non_far_interactions)) == 0:
        return sc

    clas = " and ".join(non_far_interactions)
    if len(non_far_interactions) == 4:
        clas = "Quadruple Interaction: " + clas
    elif len(non_far_interactions) == 3:
        clas = "Triple Interaction: " + clas
    elif len(non_far_interactions) == 2:
        clas = "Double Interaction: " + clas
    elif len(non_far_interactions) == 1:
        clas = clas + " Interaction"

    return clas


# called by output, needs update, currently no gene_score


def goTermAnalysis(session_id, outfile, db, cursor):
    sql = "SELECT Protein FROM RS_Protein_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)
    datasetsize = len(results)
    go_dict = {}
    gene_dict = {}
    Protein_Id_list = set([])
    for row in results:
        gene_score = row[1]
        if gene_score is None:
            gene_score = 0.0
        gene_dict[row[0]] = gene_score
        Protein_Id_list.add(row[0])

    sql = "SELECT GO_Term,Protein FROM RS_Protein_GO_Term"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)
    go_term_id_list = set([])
    for row in results:
        if not row[1] in Protein_Id_list:
            continue
        gene_score = gene_dict[row[1]]
        if not row[0] in go_dict:
            go_term_id_list.add(row[0])
            go_dict[row[0]] = ["", "", gene_score, 1.0, gene_score]
        else:
            i = go_dict[row[0]][3] + 1
            avg_score = go_dict[row[0]][2]
            go_dict[row[0]][2] = ((i - 1) / i) * avg_score + (1 / i) * gene_score
            go_dict[row[0]][3] = i
            go_dict[row[0]][4] += gene_score

    sql = "SELECT GO_Term_Id,Name,Id FROM GO_Term"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)

    for row in results:
        if not row[0] in go_term_id_list:
            continue
        go_dict[row[0]][0] = row[1]
        go_dict[row[0]][1] = row[2]

    go_list = sorted(list(go_dict.values()), key=lambda x: x[4], reverse=True)
    lines = ["GO-Term\tGO-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
    for go in go_list:
        normalized_score = go[4] / datasetsize
        lines.append("%s\t%s\t%s\t%s\t%s\t%s" % (str(go[0]), str(go[1]), str(go[4]), str(go[2]), str(int(go[3])), str(normalized_score)))
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()


# called by output, needs update, currently no gene_score
def pathwayAnalysis(session_id, outfile, db, cursor):
    sql = "SELECT Gene,Gene_Score FROM RS_Protein_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)
    datasetsize = len(results)
    path_dict = {}
    gene_dict = {}
    Protein_Id_list = set([])
    for row in results:
        gene_score = row[1]
        if gene_score is None:
            gene_score = 0.0
        gene_dict[row[0]] = gene_score
        Protein_Id_list.add(row[0])

    sql = "SELECT Pathway,Gene FROM RS_Protein_Pathway"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)
    pathway_id_list = set([])
    for row in results:
        if not row[1] in Protein_Id_list:
            continue
        gene_score = gene_dict[row[1]]
        if not row[0] in path_dict:
            pathway_id_list.add(row[0])
            path_dict[row[0]] = ["", "", gene_score, 1.0, gene_score]
        else:
            i = path_dict[row[0]][3] + 1
            avg_score = path_dict[row[0]][2]
            path_dict[row[0]][2] = ((i - 1) / i) * avg_score + (1 / i) * gene_score
            path_dict[row[0]][3] = i
            path_dict[row[0]][4] += gene_score

    if len(pathway_id_list) == 0:
        lines = ["Pathway\tReactome-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
        page = "\n".join(lines)
        f = open(outfile, "wb")
        f.write(page)
        f.close()
        return

    sql = "SELECT Pathway_Id,Name,Reactome_Id FROM Pathway"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)

    for row in results:
        if not row[0] in pathway_id_list:
            continue
        path_dict[row[0]][0] = row[1]
        path_dict[row[0]][1] = row[2]

    path_list = sorted(list(path_dict.values()), key=lambda x: x[4], reverse=True)
    lines = ["Pathway\tReactome-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
    for path in path_list:
        normalized_score = path[4] / datasetsize
        lines.append("%s\t%s\t%s\t%s\t%s\t%s" % (str(path[0]), str(path[1]), str(path[4]), str(path[2]), str(int(path[3])), str(normalized_score)))
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()


# called by postAnnoAnno
def updateAnnoAnno(anno_anno_map, session_id, db, cursor):
    valuestrs = []
    for (m_id1, m_id2) in anno_anno_map:
        (min_d, atom, atom2, t_id1, t_id2, chain1, chain2) = anno_anno_map[(m_id1, m_id2)]
        valuestrs.append("('%s','%s','%s','%s','%s','%s','%s','%s','%s %s')" % (str(t_id1), str(t_id2), str(m_id1), str(m_id2), chain1, chain2, str(session_id), str(min_d), str(atom), str(atom2)))
    sql = "INSERT IGNORE INTO RS_Annotation_Annotation(Template_1,Template_2,Mutation_1,Mutation_2,Chain_1,Chain_2,Session,Distance,Atompair) VALUES %s" % ','.join(valuestrs)
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in updateAnnoAnno: %s" % (f))


# called by calcAnnoRate
def calculateAnnotationRate(db, cursor, session_id):
    sql = "SELECT Position FROM RS_Position_Session WHERE Session = %s" % session_id
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in calculateAnnotationRate: %s" % sql)
    mut_id_list = set()
    for x in results:
        mut_id_list.add(str(x[0]))

    if len(mut_id_list) > 0:
        sql = "SELECT Position_Id,Gene FROM Position WHERE Position_Id in (%s)" % ",".join(mut_id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        gene_map = {}
        Protein_Id_list = set()
        for row in results:
            Protein_Id = str(row[1])
            mut_id = str(row[0])
            gene_map[mut_id] = Protein_Id
            Protein_Id_list.add(Protein_Id)

        sql = "SELECT Protein_Id,Error_Code FROM Protein WHERE Protein_Id in (%s)" % ",".join(Protein_Id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        error_map = {}
        for row in results:
            error_map[str(row[0])] = str(row[1])

        sql = "SELECT Position,Template,Error_Code FROM RS_Position_Template WHERE Position in (%s)" % ",".join(mut_id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        anno_mut_id_list = set()
        all_error_muts = set()
        fine_templates = set()

        for row in results:
            mut_id = str(row[0])
            if row[2] is None:
                anno_mut_id_list.add(mut_id)
                all_error_muts.discard(mut_id)
                fine_templates.add(row[1])
            elif mut_id not in anno_mut_id_list:
                all_error_muts.add(mut_id)

        sql = "SELECT Gene,Sequence_Identity,Template_Id FROM Template WHERE Gene in (%s)" % ",".join(Protein_Id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        number_of_templates = len(results)
        greater90 = set()
        greater95 = set()
        greater98 = set()
        greater99 = set()
        for row in results:
            if row[2] in fine_templates:
                g_id = str(row[0])
                s_id = float(row[1])
                if s_id > 99.0:
                    greater99.add(g_id)
                if s_id > 98.0:
                    greater98.add(g_id)
                if s_id > 95.0:
                    greater95.add(g_id)
                if s_id > 90.0:
                    greater90.add(g_id)

        muts_without_template = set()
        muts_with_unknown_gene_error = set()
        muts_with_anno = set()
        muts_with_temp_98 = 0
        for m in mut_id_list:
            if m in anno_mut_id_list:
                muts_with_anno.add(m)
                if gene_map[m] in greater98:
                    muts_with_temp_98 += 1
            if error_map[gene_map[m]] == '3':
                muts_without_template.add(m)
            elif error_map[gene_map[m]] == '4':
                muts_with_unknown_gene_error.add(m)

        return (len(mut_id_list), len(muts_with_anno), len(muts_without_template), len(muts_with_unknown_gene_error), len(all_error_muts), float(len(muts_with_anno)) / float(len(mut_id_list)), number_of_templates, len(greater90), len(greater95), len(greater98), len(greater99), muts_with_temp_98)
    else:
        return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
