import errno
import gc
import os
import shutil
import sys
import time
import traceback
import numpy as np
import psutil

import ray

#from structman.lib import annovar
from structman.lib import globalAlignment
from structman.lib import indel_analysis
from structman.lib import MMseqs2
from structman.lib import pdbParser
from structman.lib import templateFiltering
from structman.lib import uniprot

from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc import position as position_package
from structman.lib.sdsc.consts import codons as codon_consts
from structman.lib.sdsc import structure as structure_package
from structman.lib.sdsc import complex as complex_package
from structman.lib.sdsc import protein as protein_package
from structman.lib.sdsc import gene as gene_package
from structman.lib.sdsc import indel as indel_package
from structman.lib.sdsc import interface as interface_package

from structman.lib.database import database
from structman.lib.database.insertion_lib import remote_insert_complexes
from structman.lib.database.retrieval_2 import get_struct_base_infos, retrieve_stored_ligands, prepare_struct_base_info_store
from structman.base_utils.base_utils import (
    calculate_chunksizes,
    pack, unpack,
    is_alphafold_model, alphafold_model_id_to_file_path,
    is_model, complex_id_to_file_path,
    aggregate_times, print_times, add_to_times,
    is_seq_valid, reset_logger_for_remotes,
    build_ragged_int_store, ragged_row, build_sequence_store, sequence_from_store
    )
from structman.base_utils.config_class import Config
from structman.base_utils.ray_utils import ray_init, ray_hack, deref, parallel_db_access, read_object_store_stats
from structman.lib.lib_utils import process_mutations_str, generate_multi_mutation, extract_dedicated_tags, fetch_msa


# Inputs with at most this many proteins do not start ray up front; see the comment
# at the deferral in main() and the one at the lazy start in core().
SMALL_INPUT_RAY_DEFERRAL_THRESHOLD: int = 10


# Taken from https://stackoverflow.com/questions/2023608/check-what-files-are-open-in-python
def list_fds():
    """List process currently open FDs and their target """
    if not sys.platform.startswith('linux'):
        raise NotImplementedError('Unsupported platform: %s' % sys.platform)

    ret = {}
    base = '/proc/self/fd'
    for num in os.listdir(base):
        path = None
        try:
            path = os.readlink(os.path.join(base, num))
        except OSError as err:
            # Last FD is always the "listdir" one (which may be closed)
            if err.errno != errno.ENOENT:
                raise
        ret[int(num)] = path

    return ret


def kill_overhead_sockets(save_ids = None):
    if not sys.platform.startswith('linux'):
        raise NotImplementedError('Unsupported platform: %s' % sys.platform)

    if save_ids is None:
        ret_save_ids = set()

    base = '/proc/self/fd'
    for num in os.listdir(base):
        path = None
        try:
            path = os.readlink(os.path.join(base, num))
        except OSError as err:
            # Last FD is always the "listdir" one (which may be closed)
            if err.errno != errno.ENOENT:
                raise
        if save_ids is None:
            ret_save_ids.add(num)
        else:
            if num in save_ids:
                continue
            if path is None:
                continue
            if path[:6] == 'socket':
                try:
                    os.close(int(num))
                except OSError:
                    pass

    if save_ids is None:
        return ret_save_ids
    return save_ids


def _rss_breakdown_mb() -> tuple[float, float, float]:
    """(VmRSS, RssAnon, RssShmem) in MiB, all 0.0 where /proc/self/status is not readable.

    VmRSS alone is misleading for this pipeline. It counts the /dev/shm pages the driver
    has touched, and ray keeps its object store in /dev/shm - so the pages ray.put wrote
    are charged to the driver's RSS *and* sized again as the object store, and a reading
    like "477 GiB in the driver, 646 GiB object store" double-counts an unknown part of
    the total. RssAnon is what the driver's own data structures actually cost and is the
    number to compare against the cgroup limit; RssShmem is the store view.
    """
    anon = 0.0
    shmem = 0.0
    total = 0.0
    try:
        with open('/proc/self/status') as fh:
            for line in fh:
                if line.startswith('VmRSS:'):
                    total = int(line.split()[1]) / 1024.0
                elif line.startswith('RssAnon:'):
                    anon = int(line.split()[1]) / 1024.0
                elif line.startswith('RssShmem:'):
                    shmem = int(line.split()[1]) / 1024.0
    except Exception:
        pass
    return total, anon, shmem


def _object_store_note() -> str:
    """' store <used>/<capacity> GiB (<n>% primary)', or '' when ray cannot be asked.

    The spill watchdog only speaks up above 90% primary, so a run that died with the
    store merely large left no record of how large. Every stage boundary carrying it
    means the next crash can be reasoned about instead of reconstructed.
    """
    stats = read_object_store_stats()
    if stats is None or stats.get('capacity', 0) <= 0:
        return ''
    GiB = 1024 ** 3
    primary_share = 100.0 * stats['primary'] / stats['capacity']
    return (f', store {stats["used"] / GiB:.1f}/{stats["capacity"] / GiB:.1f} GiB '
            f'({primary_share:.0f}% primary)')


def _rss_mb() -> float:
    """Return current process RSS in MiB.

    On Linux reads /proc/self/status (VmRSS) which reflects the live RSS.
    Falls back to resource.getrusage ru_maxrss, which is only a peak
    high-water mark and should NOT be used to measure freed memory.
    """
    total, _anon, _shmem = _rss_breakdown_mb()
    if total:
        return total
    try:
        import resource
        # ru_maxrss is in KiB on Linux, bytes on macOS
        scale = 1024.0 if sys.platform.startswith('linux') else 1024.0 * 1024.0
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / scale
    except Exception:
        return 0.0


def clean_chunk_memory(config=None):
    """Free memory after each chunk.

    1. Three gc.collect() passes (gen 0→2) to break reference cycles in the
       driver process.
    2. ray._private.internal_api.global_gc() so every Ray worker also runs GC,
       releasing worker-side buffers pinned by finished remote tasks.
    3. malloc_trim(0) on Linux: forces glibc to return freed heap pages to the
       OS.  Python frees objects but the libc allocator holds onto arena pages
       for fast reuse; without trim, RSS grows monotonically even when Python's
       own heap is logically empty.

    Note: this intentionally does NOT call ray._private.internal_api.free() on
    tracked ray.put() objects, and does NOT restart the Ray cluster between
    chunks. Both were tried and caused native SIGSEGV crashes: free() forcibly
    evicts plasma objects regardless of whether something elsewhere (e.g. a
    zero-copy numpy view obtained via ray.get()) still references that memory,
    and Ray's docs explicitly warn it causes "severe undefined behavior" for
    objects still in use. Worker-side RSS growth is instead bounded via
    max_calls=1 on the @ray.remote functions that accumulate significant memory
    (filter_structures_wrapper, parse_complexes), which recycles worker processes
    after each task. para_parse_pdb also uses max_calls=1 to avoid reusing workers that previously
    ran annotate_wrapper, which leaves them in a bad/transitioning state. The initial
    task burst is throttled (proc_n//4) to prevent NFS thundering-herd on startup.
    """
    verbosity = config.verbosity if config is not None else 0

    if verbosity >= 2:
        rss_before = _rss_mb()

    # --- Step 1: Python GC cycles ---
    for generation in range(3):
        gc.collect(generation)

    # --- Step 2: Ray worker GC ---
    if ray.is_initialized():
        try:
            import ray._private.internal_api as _ray_api
            _ray_api.global_gc()
        except Exception as exc:
            if verbosity >= 2:
                config.logger.info(f'clean_chunk_memory: global_gc failed: {exc}')

    # --- Step 3: return freed glibc pages to the OS ---
    if sys.platform.startswith('linux'):
        try:
            import ctypes
            ctypes.cdll.LoadLibrary('libc.so.6').malloc_trim(0)
        except Exception as exc:
            if verbosity >= 2:
                config.logger.info(f'clean_chunk_memory: malloc_trim failed: {exc}')

    if verbosity >= 2:
        rss_after = _rss_mb()
        config.logger.info(f'clean_chunk_memory: RSS {rss_before:.0f} -> {rss_after:.0f} MiB (freed {rss_before - rss_after:.0f} MiB)')

    if verbosity >= 3:
        import tracemalloc
        if not tracemalloc.is_tracing():
            tracemalloc.start()
        snapshot = tracemalloc.take_snapshot()
        top = snapshot.statistics('lineno')[:10]
        config.logger.info('tracemalloc top-10 allocations after cleanup:')
        for stat in top:
            config.logger.info(f'  {stat}')


def SQLDateTime():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())


# For memory profiling, thanks to Fred Cirera
def sizeof_fmt(num, suffix='B'):
    ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def remove_sanity_filtered(config, proteins, indels, primary_protein_id):
    del_list = []
    for pos_obj in proteins[primary_protein_id].positions:
        if pos_obj is None:
            continue
        if not pos_obj.checked:
            del_list.append((pos_obj.pos, pos_obj.wt_aa))
    for pos, wt_aa in del_list:
        remove_position(config, proteins, indels, pos, wt_aa, primary_protein_id, degrade_to_session_less = True)


def remove_position(config, proteins, indels, pos, wt_aa, primary_protein_id, degrade_to_session_less = False):
    if pos == 0:
        return
    if degrade_to_session_less:
        warn_text = 'Degraded position through sanity check:%s,%s' % (primary_protein_id, pos)
    else:
        warn_text = 'Filtered position through sanity check:%s,%s' % (primary_protein_id, pos)

    config.errorlog.add_warning(warn_text, logger= config.logger)

    if proteins[primary_protein_id].multi_mutations is not None:
        mm_dels = []
        for mm_nr, multi_mutation in enumerate(proteins[primary_protein_id].multi_mutations):
            for mut in multi_mutation:
                if isinstance(mut, tuple):
                    position, aa2 = mut
                    if position.pos == pos and wt_aa == position.wt_aa:
                        mm_dels.append(mm_nr)
        for mm_nr in sorted(mm_dels, reverse=True):
            del proteins[primary_protein_id].multi_mutations[mm_nr]
        if len(proteins[primary_protein_id].multi_mutations) == 0:
            proteins[primary_protein_id].multi_mutations = None

    proteins[primary_protein_id].positions[pos].session_less = True

    if len(proteins[primary_protein_id].positions) == 0:
        if primary_protein_id in indels:
            if len(indels[primary_protein_id]) == 0:
                del proteins[primary_protein_id]
        else:
            del proteins[primary_protein_id]

@ray.remote(max_calls=1, max_retries=10)
def para_parse_pdb(config: Config, package: list[str], pdb_path: str) -> list[tuple[str, str, dict]]:
    reset_logger_for_remotes(config)
    if config.verbosity >= 4:
        rss0 = psutil.Process().memory_info().rss / 1024**2
        config.logger.info(f'[para_parse_pdb] start: {len(package)} PDB chains, rss={rss0:.0f} MB, avail={psutil.virtual_memory().available / 1024**2:.0f} MB')
    try:
        pdb_sequence_list: list[tuple[str, str, dict]] = pdbParser.getSequences(package, pdb_path, verbosity=config.verbosity)
    except Exception as e:
        config.logger.error(f'[para_parse_pdb] EXCEPTION {type(e).__name__}: {e}')
        raise
    if config.verbosity >= 4:
        rss1 = psutil.Process().memory_info().rss / 1024**2
        config.logger.info(f'[para_parse_pdb] done: {len(pdb_sequence_list)} results, rss={rss1:.0f} MB (+{rss1-rss0:.0f} MB)')
    return pdb_sequence_list

def map_pdb_sequences_onto_positions(pdb_sequence_list, sequenceScanPDB, proteins, removed_proteins):
    """Attach parsed PDB chain sequences to their Protein objects as Positions.

    Lifted out of sequenceScan so the hot path can be exercised directly by
    dev_tests/test_pdb_sequence_mapping.py, which pins it against the previous
    implementation. Mutates proteins and appends to removed_proteins.
    """
    # This was the slowest step of sequenceScan by a wide margin - 2989s for
    # 149773 chains, against 54s for the parallel PDB parsing that feeds it -
    # and almost none of that was the loop's own arithmetic. It allocates one
    # Position per residue, tens of millions of them, while the process already
    # holds a very large live heap, so every gen-2 collection walks that whole
    # heap and reclaims nothing: measured on a 20M object heap the collector
    # alone costs 4.6x. Nothing built here takes part in reference cycles, so
    # the collector has nothing to find. Same reasoning and same fix as the hit
    # parser in MMseqs2.parseHits.
    gc_was_enabled: bool = gc.isenabled()
    gc.disable()
    try:
        for pdb_chain_tuple, seq, res_pos_map in pdb_sequence_list:
            if seq == '':
                removed_proteins.append(pdb_chain_tuple)
                continue

            # Resolved once per chain. The nested loop below used to re-read
            # these out of proteins up to ten times per residue.
            protein_obj = proteins[pdb_chain_tuple]
            protein_obj.sequence = seq

            # add_residues (called from Protein.__init__/integrate_protein, before the
            # structure was parsed) filed any input residue-named position - a
            # PDB-style mutation such as '1ABC:A_123A>V' - past positions[0] by PDB
            # residue number, since the sequence position was not known yet.
            # positions[0] itself is the "tag this chain, no specific residue" sentinel,
            # already consumed by popNone earlier in sequenceScan. Pull the pending
            # residues into a local dict - built fresh per chain, so its size is exactly
            # how many residues the input actually named, typically none - and start
            # positions over for the real, sequence-indexed fill below. This used to be
            # a dict living on the protein itself (res_id_map), written to once per
            # residue of every chain by the uni_pos branch below even though nothing
            # ever read those entries back - every position it named was already
            # reachable through positions too.
            pending_by_res_id = {
                pending.pdb_res_nr: pending
                for pending in protein_obj.positions[1:] if pending is not None
            }
            positions = [protein_obj.positions[0]]
            protein_obj.positions = positions

            tags, uni_pos = sequenceScanPDB[pdb_chain_tuple]
            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False

            # pdbParser.parsePDBSequence builds this as
            # res_pos_map[res_nr] = len(res_pos_map) + 1, guarded by a used_res
            # set, so it maps distinct residue ids onto exactly 1..len(seq) with
            # no gaps, in ascending order. Walking it directly therefore visits
            # the same positions in the same order as enumerate(seq) did, without
            # first building the inverted backmap and without the per-residue
            # membership test that this invariant made unreachable.
            new_positions = []
            for res_id, seq_pos in res_pos_map.items():
                existing = pending_by_res_id.get(res_id)
                if existing is not None:
                    existing.pos = seq_pos
                    # Only reachable for a protein that already carried
                    # positions. The batch below has not grown the list yet, so
                    # pad here exactly as add_positions would have.
                    if seq_pos >= len(positions):
                        positions.extend([None] * (seq_pos + 1 - len(positions)))
                    positions[seq_pos] = existing
                    existing.add_tags(tags)
                elif uni_pos:
                    # contains_position(seq_pos), inlined against the local alias.
                    if seq_pos >= len(positions) or positions[seq_pos] is None:
                        position = position_package.Position(pos=seq_pos, pdb_res_nr=res_id, wt_aa=seq[seq_pos - 1], tags=tags, session_less = session_less)
                        new_positions.append(position)
                    # No else branch: a position already at positions[seq_pos] needs
                    # nothing more done to it, and nothing reads a resolved residue
                    # back by its PDB number again.

            # One call per chain rather than one per residue. add_positions
            # allocates a warning list and a deletion list and runs a sorted() on
            # every call, which made it 93% of this loop's own runtime.
            if new_positions:
                protein_obj.add_positions(new_positions)
    finally:
        if gc_was_enabled:
            gc.enable()



def sequenceScan(
        config: Config,
        proteins: dict[str, protein_package.Protein],
        indels:  dict[str, indel_package.Indel],
        genes: dict[str, gene_package.Gene] | None = None,
        skip_mm_calc: bool = False):
    def promote_uac_to_primary_id(u_ac, np_ref):
        primary_protein_id = u_ac
        proteins[primary_protein_id] = proteins[np_ref]

        if primary_protein_id in indels and np_ref in indels:
            indels[primary_protein_id].update(indels[np_ref])
        elif np_ref in indels:
            indels[primary_protein_id] = indels[np_ref]

        if primary_protein_id in indels:
            for indel_notation in indels[primary_protein_id]:
                indels[primary_protein_id][indel_notation].wt_prot = primary_protein_id

        if np_ref in indels:
            del indels[np_ref]
        del proteins[np_ref]
        proteins[primary_protein_id].primary_protein_id = primary_protein_id
        proteins[primary_protein_id].u_ac = u_ac
        return u_ac

    pdb_path = config.pdb_path

    sequenceScanProteins = {}
    sequenceScanPDB = {}
    sequenceScanNM = {}
    sequenceScanNP = {}
    sequenceScanEnsembl = {}
    fasta_inputs = {}

    removed_proteins = []

    if config.verbosity >= 2:
        if len(proteins) <= 100:
            config.logger.info(f'Protein IDs going into sequence scan: {proteins}, genes set: {genes is not None}')
        else:
            config.logger.info(f'{len(proteins)} number of Protein IDs going into sequence scan, genes set: {genes is not None}')
        t0 = time.time()

    for prot_id in proteins:
        uni_pos, tags = proteins[prot_id].popNone()

        if proteins[prot_id].is_sequence_set():
            if uni_pos:
                fasta_inputs[prot_id] = tags
            continue
        if prot_id.count(':') > 0:
            sequenceScanPDB[prot_id] = tags, uni_pos  # PDB inputs are always processed by the sequence scan, but the positions are only added if uni_pos is true
        elif not sdsc_utils.is_mutant_ac(prot_id):
            if proteins[prot_id].input_id[:2] == 'NM':
                if proteins[prot_id].input_id in indels:  # Indel proteins need complete classification
                    uni_pos = True
                sequenceScanNM[proteins[prot_id].input_id] = tags, uni_pos, proteins[prot_id].u_ac
            elif proteins[prot_id].input_id[:2] == 'NP' or proteins[prot_id].input_id[-2] == '.':
                if proteins[prot_id].input_id in indels:  # Indel proteins need complete classification
                    uni_pos = True
                sequenceScanNP[proteins[prot_id].input_id] = tags, uni_pos, proteins[prot_id].u_ac
            elif proteins[prot_id].input_id[:3] == 'ENS' and ((proteins[prot_id].input_id == proteins[prot_id].primary_protein_id) or genes is not None):
                sequenceScanEnsembl[proteins[prot_id].input_id] = tags, uni_pos
            else:
                if prot_id in indels:  # Indel proteins need complete classification
                    uni_pos = True
                sequenceScanProteins[proteins[prot_id].u_ac] = tags, uni_pos  # New: process everything to filter input by sanity checks

    if len(sequenceScanNP) > 0:
        if config.verbosity >= 1:
            config.logger.info(f"Amount of proteins going into NP sequenceScan: {len(sequenceScanNP)}")
        gene_sequence_map = uniprot.get_refseq_sequences(','.join(sequenceScanNP.keys()), config, seq_type='protein')
        stems = set()
        ref_stem_map = {}
        for np_ref in sequenceScanNP:
            (tags, uni_pos, u_ac) = sequenceScanNP[np_ref]
            if u_ac is None:
                continue
            if u_ac[:3] == 'NP_':
                continue
            u_ac_stem = u_ac.split('_')[0].split('-')[0]
            stems.add(u_ac_stem)
            ref_stem_map[np_ref] = u_ac_stem
        isoform_specific_id_map = uniprot.u_ac_isoform_search(gene_sequence_map, stems, ref_stem_map, config)
        for np_ref in sequenceScanNP:
            (tags, uni_pos, u_ac) = sequenceScanNP[np_ref]
            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False
            seq = gene_sequence_map[np_ref]
            if np_ref in isoform_specific_id_map:
                primary_protein_id = promote_uac_to_primary_id(isoform_specific_id_map[np_ref], np_ref)
            else:
                primary_protein_id = np_ref

            proteins[primary_protein_id].sequence = seq

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if proteins[primary_protein_id].contains_position(seq_pos):
                    proteins[primary_protein_id].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[primary_protein_id].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if not proteins[primary_protein_id].contains_position(seq_pos):
                        position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True, session_less = session_less)
                        proteins[primary_protein_id].add_positions([position])

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, primary_protein_id)

    if len(sequenceScanNM) > 0:
        if config.verbosity >= 1:
            config.logger.info(f"Amount of proteins going into NM sequenceScan: {len(sequenceScanNM)}")
        gene_sequence_map = uniprot.get_refseq_sequences(','.join(sequenceScanNM.keys()), config)
        ref_stem_map = {}
        stems = set()
        for nm_ref in sequenceScanNM:
            (tags, uni_pos, u_ac) = sequenceScanNM[nm_ref]
            if u_ac is None:
                continue
            if u_ac[:3] == 'NM_':
                continue
            u_ac_stem = u_ac.split('_')[0].split('-')[0]
            stems.add(u_ac_stem)
            ref_stem_map[nm_ref] = u_ac_stem

        isoform_specific_id_map = uniprot.u_ac_isoform_search(gene_sequence_map, stems, ref_stem_map, config)
        for nm_ref in sequenceScanNM:
            (tags, uni_pos, u_ac) = sequenceScanNM[nm_ref]
            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False
            if nm_ref not in gene_sequence_map:
                try:
                    proteins.remove_protein(nm_ref)
                except:
                    try:
                        del proteins[nm_ref]
                    except:
                        pass
                if config.verbosity >= 3:
                    config.logger.info(f'Removed {nm_ref} from Input')
                removed_proteins.append(nm_ref)
                continue

            seq = gene_sequence_map[nm_ref]
            if nm_ref in isoform_specific_id_map:
                primary_protein_id = promote_uac_to_primary_id(isoform_specific_id_map[nm_ref], nm_ref)
            else:
                primary_protein_id = nm_ref

            proteins[primary_protein_id].sequence = seq

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if proteins[primary_protein_id].contains_position(seq_pos):
                    proteins[primary_protein_id].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[primary_protein_id].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if not proteins[primary_protein_id].contains_position(seq_pos):
                        position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True, session_less = session_less)
                        proteins[primary_protein_id].add_positions([position])

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, primary_protein_id)

    if len(sequenceScanEnsembl) > 0:
        if config.verbosity >= 1:
            t01 = time.time()
            config.logger.info(f"Amount of proteins going into Ensembl sequenceScan: {len(sequenceScanEnsembl)}")
        
        gene_sequence_map, gene_id_map = uniprot.get_ensembl_seqs(config, list(sequenceScanEnsembl.keys()))

        if config.verbosity >= 2:
            t02 = time.time()
            config.logger.info(f"Size of gene_sequence_map: {len(gene_sequence_map)}, size of gene_id_map: {len(gene_id_map)}")
            config.logger.info(f'Time for uniprot.get_ensembl_seqs: {t02- t01}')

        mean_seq_len = 0
        n = 0

        for transcript_id in sequenceScanEnsembl:
            tags, uni_pos = sequenceScanEnsembl[transcript_id]
            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False

            if (transcript_id not in gene_sequence_map):# and (transcript_id not in gene_id_map):
                #This happens for problematic Ensembl transcript ids, for example for entries without an open reading frame (due to the transcript being non-coding)
                try:
                    proteins.remove_protein(transcript_id)
                except:
                    try:
                        del proteins[transcript_id]
                    except:
                        pass
                if config.verbosity >= 3:
                    config.logger.info(f'Removed {transcript_id} from Input')
                removed_proteins.append(transcript_id)
                continue
            
            if gene_sequence_map[transcript_id] is None:
                try:
                    proteins.remove_protein(transcript_id)
                except:
                    try:
                        del proteins[transcript_id]
                    except:
                        pass
                if config.verbosity >= 3:
                    config.logger.info(f'Removed {transcript_id} from Input')
                removed_proteins.append(transcript_id)
                continue

            if genes is not None:
                gene_id, gene_name = gene_id_map[transcript_id]
                primary_protein_id = proteins[transcript_id].primary_protein_id
                if gene_id not in genes:
                    gene_obj = gene_package.Gene(gene_id, proteins = {transcript_id:primary_protein_id}, gene_name = gene_name)
                    genes[gene_id] = gene_obj
                else:
                    genes[gene_id].proteins[transcript_id] = primary_protein_id
                    genes[gene_id].gene_name = gene_name

                proteins[transcript_id].gene = gene_id

            primary_protein_id = transcript_id

            if transcript_id in gene_sequence_map:
                seq = gene_sequence_map[transcript_id]

                prot_obj = proteins[primary_protein_id]

                prot_obj.sequence = seq
                n += 1
                mean_seq_len += len(seq)

                for (pos, aa) in enumerate(seq):
                    seq_pos = pos + 1
                    seq_pos_there = prot_obj.contains_position(seq_pos)
                    if seq_pos_there:
                        prot_obj.positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                        prot_obj.positions[seq_pos].add_tags(tags)
                    elif uni_pos:
                        if not seq_pos_there:
                            position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True, session_less = session_less)
                            prot_obj.add_positions([position])

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, primary_protein_id)
            
        if genes is not None:
            try:
                del genes['to_fetch']
            except:
                pass

        if n > 0:
            mean_seq_len = mean_seq_len/n
        if config.verbosity >= 2:
            config.logger.info(f'Retrieved {n} number of sequences in sequenceScanEnsembl with mean length: {mean_seq_len}')
            t03 = time.time()
            config.logger.info(f'Time for sequenceScan ensembl Part 2: {t03- t02}')


    if len(sequenceScanProteins) > 0:
        if config.verbosity >= 1:
            config.logger.info(f"Amount of proteins going into Uniprot sequenceScan: {len(sequenceScanProteins)}")

        gene_sequence_map = uniprot.getSequencesPlain(sequenceScanProteins.keys(), config)

        for u_ac in gene_sequence_map:
            if gene_sequence_map[u_ac] == 1 or gene_sequence_map[u_ac] == 0 or gene_sequence_map[u_ac] is None:
                config.errorlog.add_warning("Error in sequenceScan with gene: %s" % u_ac)
                continue
            seq = gene_sequence_map[u_ac]
            proteins[u_ac].sequence = seq

            tags, uni_pos = sequenceScanProteins[u_ac]

            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False
            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if proteins[u_ac].contains_position(seq_pos):
                    proteins[u_ac].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[u_ac].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if not proteins[u_ac].contains_position(seq_pos):
                        position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True, session_less = session_less)
                        proteins[u_ac].add_positions([position])

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, u_ac)
            del seq

        del gene_sequence_map
    
    if len(sequenceScanPDB) > 0:
        if config.verbosity >= 2:
            t_pdb_0 = time.time()
            config.logger.info(f"Amount of PDB chains going into PDB sequenceScan: {len(sequenceScanPDB)}")
        # parallel_db_access is False for read_only and proc_n == 1, where ray is never
        # started - the serial branch below is the only one that may run there.
        if not parallel_db_access(config) or len(sequenceScanPDB) < config.proc_n:
            pdb_sequence_list = pdbParser.getSequences(sequenceScanPDB.keys(), pdb_path, verbosity=config.verbosity)
        else:
            # Package size is capped by config.safe_task_chunksize (a RAM-based value, unaffected by a
            # user-supplied --chunksize) and packages are kept to at most config.proc_n in flight at once,
            # so a large --chunksize (which only controls how much gets read per outer batch) can't blow up
            # the amount of parsed PDB data resident in memory at any one time.
            task_size: int = max(1, min(config.safe_task_chunksize, -(-len(sequenceScanPDB) // config.proc_n)))
            pdb_ids: list[str] = list(sequenceScanPDB.keys())
            pending_packages: list[list[str]] = [pdb_ids[i:i + task_size] for i in range(0, len(pdb_ids), task_size)]

            parse_pdb_results: list[list[tuple[str, str, dict]]] = []
            in_flight: list = []
            next_package_idx = 0

            # One task per package, and every one of them was re-serialising the same
            # config into its task spec. ray resolves the ref back to the Config on the
            # worker, so para_parse_pdb's signature is unchanged.
            config_store = ray.put(config)

            # Warm up with a smaller initial burst to avoid NFS thundering-herd during
            # simultaneous worker startup. After the first wave workers are warm.
            initial_burst = max(1, config.proc_n // 4)
            while next_package_idx < len(pending_packages) and len(in_flight) < initial_burst:
                in_flight.append(para_parse_pdb.remote(config_store, pending_packages[next_package_idx], pdb_path))
                next_package_idx += 1

            while len(in_flight) > 0:
                ready, in_flight = ray.wait(in_flight)
                parse_pdb_results.extend(ray.get(ready))
                while next_package_idx < len(pending_packages) and len(in_flight) < config.proc_n:
                    in_flight.append(para_parse_pdb.remote(config_store, pending_packages[next_package_idx], pdb_path))
                    next_package_idx += 1

            pdb_sequence_list = []
            for pdb_sequence_sublist in parse_pdb_results:
                pdb_sequence_list.extend(pdb_sequence_sublist)

        if config.verbosity >= 2:
            t_pdb_1 = time.time()
            config.logger.info(f'Time for PDB file parsing (getSequences/para_parse_pdb): {t_pdb_1 - t_pdb_0}')

        map_pdb_sequences_onto_positions(pdb_sequence_list, sequenceScanPDB, proteins, removed_proteins)

        for pdb_chain_tuple in sequenceScanPDB:
            if pdb_chain_tuple not in proteins:
                removed_proteins.append(pdb_chain_tuple)

        if config.verbosity >= 2:
            t_pdb_2 = time.time()
            config.logger.info(f'Time for mapping PDB sequences onto positions: {t_pdb_2 - t_pdb_1}')

    if len(fasta_inputs) > 0:
        if config.verbosity >= 2:
            t_fasta_0 = time.time()
            config.logger.info(f"Amount of proteins going into fasta sequenceScan: {len(fasta_inputs)}")

        for prot_id in fasta_inputs:
            seq = proteins[prot_id].sequence

            tags = fasta_inputs[prot_id]

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if not proteins[prot_id].contains_position(seq_pos):
                    position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True)
                    proteins[prot_id].add_positions([position])

        if config.verbosity >= 2:
            t_fasta_1 = time.time()
            config.logger.info(f'Time for fasta sequenceScan: {t_fasta_1 - t_fasta_0}')

    if config.verbosity >= 2:
        t1 = time.time()
        config.logger.info(f'sequenceScan Part 1: {t1-t0}')

    if config.verbosity >= 3:
        config.logger.info(f'Before indel mutation with: {len(indels)} number of indels')

    for primary_protein_id in indels:
        for indel_notation in indels[primary_protein_id]:
            indels[primary_protein_id][indel_notation].mutate_sequence(proteins)

    if config.verbosity >= 2:
        t2 = time.time()
        config.logger.info(f'sequenceScan Part 2: {t2-t1}')

    multi_mutation_objects = []

    if (not config.only_wt) and (not skip_mm_calc):
        if config.verbosity >= 3:
            config.logger.info('Before multi mutations mutation')
        for primary_protein_id in list(proteins.keys()).copy():  # iterate over snapshot, since the mutant proteins get added to proteins in the loop
            if proteins[primary_protein_id].multi_mutations is not None:
                # creates mutation protein object
                multi_mutation_objects += proteins[primary_protein_id].create_multi_mutations(proteins, config)

    if config.verbosity >= 2:
        t3 = time.time()
        config.logger.info(f'sequenceScan Part 3: {t3-t2}')

    return proteins, indels, multi_mutation_objects, removed_proteins, genes


def parseFasta(config: Config, nfname: str):
    f = open(nfname, 'r')
    lines = f.readlines()
    f.close()

    seq_map = {}

    prot_map = {}

    del_map = set()
    positions = []
    multi_mutations = []
    prot_gene_map = {}
    prot_tags_map = {}

    t1_sum = 0.
    t2_sum = 0.
    t3_sum = 0.

    for line in lines:
        t0 = time.time()
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            pos_set = None
            words = line[1:].split()
            entry_id = words[0]

            seq_map[entry_id] = ''
            if len(words) > 1:
                aacs = words[1]
                aac_tag_tuples = aacs.split(';')
                mm = []
                multi_mutation_tags = None
                for aac_tag_tuple in aac_tag_tuples:
                    if aac_tag_tuple.count('<') == 1:
                        aac, tags = aac_tag_tuple.split('<')
                    else:
                        aac = aac_tag_tuple
                        tags = set()

                    pr_mut_str_result = process_mutations_str(config, aac, tags, pos_set = pos_set)

                    if pr_mut_str_result is None:
                        config.errorlog.add_warning(f'Couldnt process mutation_str: {aac}, {tags}')
                        continue

                    positions, new_multi_mutations, gene_id, protein_specific_tags, multi_mutation_tags, pos_set = pr_mut_str_result
                    mm += new_multi_mutations
                    if gene_id is not None:
                        prot_gene_map[entry_id] = gene_id

                    prot_tags_map[entry_id] = protein_specific_tags

                multi_mutations.append((mm, multi_mutation_tags))

            t1 = time.time()
            t1_sum += (t1-t0)
        elif line[0] == '<':
            words = line[1:].split()

            if len(words) == 2:
                aac, tags = words
            else:
                aac = words[0]
                tags = set()
            pr_mut_str_result = process_mutations_str(config, aac, tags, pos_set=pos_set)

            if pr_mut_str_result is None:
                config.errorlog.add_warning(f'Couldnt process mutation_str: {aac}, {tags}')
                continue

            new_positions, new_multi_mutations, gene_id, protein_specific_tags, multi_mutation_tags, pos_set = pr_mut_str_result

            positions += new_positions
            multi_mutations.append((new_multi_mutations, multi_mutation_tags))
            if gene_id is not None:
                prot_gene_map[entry_id] = gene_id
            #prot_tags_map[entry_id] = protein_specific_tags

            t2 = time.time()
            t2_sum += (t2-t0)
        else:
            if len(positions) == 0:
                position = position_package.Position()
                positions = [position]
            add_to_prot_map(prot_map, entry_id, positions, multi_mutations, config)
            positions = []
            multi_mutations = []
            cleaned_line = line.replace('\n', '').replace('/','').replace('*','').upper()
            for char_pos, char in enumerate(cleaned_line):
                if ord(char) < 65 or ord(char) > 90:
                    config.errorlog.add_warning(f'There is a not supported character in sequence at position {char_pos} for entry: {entry_id} --> Entry omitted!')
                    del_map.add(entry_id)
                    break
            seq_map[entry_id] += cleaned_line

            t3 = time.time()
            t3_sum += (t3-t0)

    if config.verbosity >= 3:
        config.logger.info(f'Parse loop times: {t1_sum} {t2_sum} {t3_sum}')

    for entry_id in seq_map:
        valid, reason = is_seq_valid(seq_map[entry_id])
        if not valid:
            config.errorlog.add_warning(f'Removed {entry_id} from input due to "{reason}"')
            del_map.add(entry_id)

    for entry_id in del_map:
        del seq_map[entry_id]

    proteins = {}
    indel_map = {}
    genes = {}

    for prot_id in seq_map:
        seq = seq_map[prot_id]

        if prot_id in prot_gene_map:
            gene_id = prot_gene_map[prot_id]
        else:
            gene_id = None

        if prot_id in prot_tags_map:
            protein_specific_tags = prot_tags_map[prot_id]
        else:
            protein_specific_tags = None
        if config.verbosity >= 3:
            config.logger.info(f'In parseFasta: Protein tags for {prot_id=} ({(prot_id in prot_tags_map)=}): {protein_specific_tags=}')
            t_0 = time.time()

        uniprot.integrate_protein(config, proteins, genes, indel_map, prot_id, prot_id, prot_map, gene_id = gene_id, protein_specific_tags = protein_specific_tags)
        proteins[prot_id].sequence = seq
        proteins[prot_id].sequence_from_input = True

        if config.verbosity >= 3:
            t_1 = time.time()
            config.logger.info(f'Time for integrate_protein: {t_1-t_0}')

    genes = gene_isoform_check(proteins, genes, indel_map, config)

    outlist: list[tuple[dict[str, protein_package.Protein], dict[str, indel_package.Indel]]] = input_chunking(config, proteins, indel_map)

    return outlist, genes


def add_to_prot_map(
        prot_map: dict[str, tuple[list[position_package.Position], list[tuple], list[tuple[list[tuple | indel_package.Indel]]]]],
        sp_id: str,
        positions: list[position_package.Position],
        multi_mutations: list[tuple | indel_package.Indel],
        config: Config):
    if sp_id not in prot_map:
        prot_map[sp_id] = [], [], []
    for multi_mutation, _ in multi_mutations:
        for indel_or_snv in multi_mutation:
            if not isinstance(indel_or_snv, tuple):
                prot_map[sp_id][1].append(indel_or_snv)
            #else:
            #    prot_map[sp_id][0].append(indel_or_snv[0])

    for position in positions:
        prot_map[sp_id][0].append(position)
        
    for multi_mutation, mm_tags in multi_mutations:
        if len(multi_mutation) > 1 and not config.only_snvs:
            prot_map[sp_id][2].append((multi_mutation, mm_tags))
    return prot_map


def input_chunking(
        config: Config,
        proteins: dict[str, protein_package.Protein],
        indels:dict[str, indel_package.Indel]
        ):
    
    outlist: list[tuple[dict[str, protein_package.Protein], dict[str, indel_package.Indel]]] = []

    s = len(proteins)
    if config.verbosity >= 1:
        config.logger.info(f"Total proteins: {s}")
    if s > config.chunksize:
        n_of_batches = s // config.chunksize
        if s % config.chunksize != 0:
            n_of_batches += 1
        batchsize = s // n_of_batches
        if s % n_of_batches != 0:
            batchsize += 1

        outlist = []

        no_more_indels = len(indels) == 0

        for i in range(0, n_of_batches):
            new_map = {}
            indel_map = {}
            prot_count = 0
            for j in range(0, batchsize):
                if prot_count >= batchsize:
                    break
                if len(proteins) == 0:
                    continue
                if no_more_indels:
                    (key, value) = proteins.popitem()
                    new_map[key] = value
                    prot_count += 1
                else:
                    key, indellist = indels.popitem()
                    new_map[key] = proteins[key]
                    del proteins[key]
                    prot_count += 1
                    for indel_notation in indellist:
                        indel = indellist[indel_notation]
                        if indel.mut_prot is not None:
                            new_map[indel.mut_prot] = proteins[indel.mut_prot]
                            del proteins[indel.mut_prot]
                            prot_count += 1
                    indel_map[key] = indellist
                    no_more_indels = len(indels) == 0
            outlist.append((new_map, indel_map))
        new_map = {}
        while len(proteins) > 0:
            (key, value) = proteins.popitem()
            new_map[key] = value
            if key in indels:
                for indel_notation in indels[key]:
                    indel = indels[key][indel_notation]
                    if indel.mut_prot is not None:
                        new_map[indel.mut_prot] = proteins[indel.mut_prot]
                        del proteins[indel.mut_prot]
                indel_map[key] = indels[key]
        if len(new_map) > 0:
            outlist.append((new_map, indel_map))
    else:
        outlist.append((proteins, indels))
    return outlist


def buildQueue(config: Config, filename: str, already_split: bool =False):
    t0 = time.time()

    proteins = {}
    u_ids = set()
    u_acs = set()
    id_map = {}
    ac_map: dict[str, tuple[list[position_package.Position], list[tuple], list[tuple[list[tuple | indel_package.Indel]]]]] = {}
    np_map: dict[str, tuple[list[position_package.Position], list[tuple], list[tuple[list[tuple | indel_package.Indel]]]]] = {}
    nm_map: dict[str, tuple[list[position_package.Position], list[tuple], list[tuple[list[tuple | indel_package.Indel]]]]] = {}
    hgnc_map: dict[str, tuple[list[position_package.Position], list[tuple], list[tuple[list[tuple | indel_package.Indel]]]]] = {}
    ensembl_map: dict[str, tuple[list[position_package.Position], list[tuple], list[tuple[list[tuple | indel_package.Indel]]]]] = {}

    pdb_map: dict[str, tuple[list[position_package.Position], list[tuple], list[tuple[list[tuple | indel_package.Indel]]]]] = {}

    if isinstance(filename, str):
        f = open(filename, "r")
        lines = f.read().split('\n')
        f.close()
    else:
        # In case of single line input
        lines = ['\t'.join(filename)]

    if config.low_mem_system and not already_split:
        prot_lines = {}
        prots = []
        for line in lines:
            # skip blank lines and remove c++ style `//` comments
            line = line.split('//')[0].strip()
            if line == '':
                continue
            if len(line) < 3:
                if config.verbosity >= 1:
                    config.logger.info(f"Skipped input line:\n{line}\nReason: Too short.\n")
                continue
            line = line.replace(' ', '\t')
            words = line.split("\t")
            if len(words) < 1:
                if config.verbosity >= 1:
                    config.logger.info(f"Skipped input line:\n{line}\nReason: Too few words.\n")
                continue
            sp_id = words[0]  # .replace("'","\\'")
            if sp_id not in prot_lines:
                prot_lines[sp_id] = [line]
                prots.append(sp_id)
            else:
                prot_lines[sp_id].append(line)

        total_num_of_raw_ids = len(prot_lines)
        if total_num_of_raw_ids > (10 * config.chunksize):
            temp_infiles = []
            num_of_infiles = total_num_of_raw_ids // (10 * config.chunksize)
            if total_num_of_raw_ids % (10 * config.chunksize) != 0:
                num_of_infiles += 1

            num_of_prots_per_file = total_num_of_raw_ids // num_of_infiles
            if total_num_of_raw_ids % num_of_infiles != 0:
                num_of_prots_per_file += 1

            for i in range(num_of_infiles):
                temp_file_lines = []
                for j in range(i * num_of_prots_per_file, (i + 1) * num_of_prots_per_file):
                    if j < len(prots):
                        temp_file_lines += prot_lines[prots[j]]
                if len(temp_file_lines) == 0:
                    continue
                temp_file_path = '%s/infile_split_%s.smlf' % (config.temp_folder, str(i))
                f = open(temp_file_path, 'w')
                f.write('\n'.join(temp_file_lines))
                f.close()
                temp_infiles.append(temp_file_path)
            return [], temp_infiles

    prot_tags_map = {}
    prot_gene_map = {}
    pos_sets = {}

    pdbs_to_expand: int = 0

    for line in lines:
        # skip blank lines and remove c++ style `//` comments
        line = line.split('//')[0].strip()
        if line == '':
            continue
        if len(line) < 3:
            if config.verbosity >= 1:
                config.logger.info(f"Skipped input line:\n{line}\nReason: Too short.\n")
            continue
        line = line.replace(' ', '\t')
        words = line.split("\t")
        if len(words) < 1:
            if config.verbosity >= 1:
                config.logger.info(f"Skipped input line:\n{line}\nReason: Too few words.\n")
            continue
        sp_id = words[0]  # .replace("'","\\'")

        if len(sp_id) < 2:
            if config.verbosity >= 1:
                config.logger.info(f"Skipped input line:\n{line}\nReason: ID too short.\n")
            continue
        elif (len(sp_id) == 6 and sp_id.count(':') == 1):
            sp_id = f'{sp_id[:4].upper()}:{sp_id[-1]}'
        elif (len(sp_id) == 4 and sp_id[0].isnumeric()):
            sp_id = sp_id.upper()

        tags = set()

        if config.verbosity >= 6:
            config.logger.info(f'Processing input line: {line}')

        if len(words) > 2:
            tags = set(words[2].split(','))

        try:
            if len(words) == 1 or words[1] == '':
                position = position_package.Position(tags=tags)
                
                multi_mutations = []
                positions: list[position_package.Position] = [position]

                tags, gene_id = extract_dedicated_tags(tags)
                if gene_id is not None:
                    prot_gene_map[sp_id] = gene_id
                prot_tags_map[sp_id] = tags

            else:
                mutation_str = words[1].replace("\n", "")
                multi_mutations = []
                positions: list[position_package.Position] = []
                if sp_id in pos_sets:
                    pos_set = pos_sets[sp_id]
                else:
                    pos_set = None
                if (not sp_id.count(':') == 1) or sp_id[0:5] == 'HGNC:':  # this means sp_id is not a pdb-id
                    positions, multi_mutation, gene_id, protein_specific_tags, multi_mutation_tags, pos_set = process_mutations_str(config, mutation_str, tags, pos_set=pos_set)

                    multi_mutations.append((multi_mutation, multi_mutation_tags))
                    if gene_id is not None:
                        prot_gene_map[sp_id] = gene_id
                    prot_tags_map[sp_id] = protein_specific_tags

                else:  # this means sp_id is a pdb-id
                    positions, multi_mutations, gene_id, protein_specific_tags, multi_mutation_tags, pos_set = process_mutations_str(config, mutation_str, tags, pos_set=pos_set, pdb_style=True)
                    multi_mutations.append((multi_mutation, multi_mutation_tags))
                    if gene_id is not None:
                        pdb_chain_tuple = '%s:%s' % (sp_id[:4].upper(), sp_id[-1])  # enforce uppercase pdb-id
                        prot_gene_map[pdb_chain_tuple] = gene_id
                    prot_tags_map[sp_id] = protein_specific_tags
                pos_sets[sp_id] = pos_set


        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            config.errorlog.add_error("File Format Error: %s\n%s\n%s\n%s" % (line, str(e), str(f), str(g)))

        if sp_id[2] == "_" or sp_id[-2] == '.':  # RefSeq Ids
            if sp_id[:2] == 'NP':
                np_map = add_to_prot_map(np_map, sp_id, positions, multi_mutations, config)
            elif sp_id[:2] == 'NM':
                nm_map = add_to_prot_map(nm_map, sp_id, positions, multi_mutations, config)
            else:
                #config.errorlog.add_warning('Unsupported Refseq Id:',line)
                # Try to put them all into the NP map
                np_map = add_to_prot_map(np_map, sp_id, positions, multi_mutations, config)
        else:
            if sp_id.count('_') > 0:
                u_ids.add(sp_id)
                id_map = add_to_prot_map(id_map, sp_id, positions, multi_mutations, config)

            elif sp_id[:3] == 'ENS' and len(sp_id) > 14:
                ensembl_map = add_to_prot_map(ensembl_map, sp_id, positions, multi_mutations, config)

            elif sp_id[:5] == 'HGNC:':
                hgnc_map = add_to_prot_map(hgnc_map, sp_id, positions, multi_mutations, config)

            elif len(sp_id) == 6 and sp_id.count(':') == 1:
                pdb_chain_tuple = '%s:%s' % (sp_id[:4].upper(), sp_id[-1])  # enforce uppercase pdb-id
                pdb_map = add_to_prot_map(pdb_map, pdb_chain_tuple, positions, multi_mutations, config)

            elif len(sp_id) == 4 and sp_id[0].isnumeric():
                pdb_chain_tuple = f'{sp_id.upper()}:-'
                pdb_map = add_to_prot_map(pdb_map, pdb_chain_tuple, positions, multi_mutations, config)
                pdbs_to_expand += 1
            else:
                u_acs.add(sp_id)
                ac_map = add_to_prot_map(ac_map, sp_id, positions, multi_mutations, config)

    t1 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f"buildQueue Part 1: {str(t1 - t0)} {len(ac_map)=} {len(prot_tags_map)=}")

    proteins: dict[str, protein_package.Protein]
    indel_map: dict[str, indel_package.Indel]
    genes: dict[str, gene_package.Gene]
    proteins, indel_map, genes = uniprot.IdMapping(config, ac_map, id_map, np_map, pdb_map, hgnc_map, nm_map, ensembl_map, prot_gene_map, prot_tags_map, pdbs_to_expand)

    t2 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f"buildQueue Part 2: {str(t2 - t1)}")

    if config.verbosity >= 3:
        config.logger.info(f'Size of proteins after IdMapping: {len(proteins)}')

    genes = gene_isoform_check(proteins, genes, indel_map, config)

    t3 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f"buildQueue Part 3: {str(t3 - t2)}")

    outlist: list[tuple[dict[str, protein_package.Protein], dict[str, indel_package.Indel]]] = input_chunking(config, proteins, indel_map)

    t4 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f"buildQueue Part 4: {str(t4 - t3)}")

    return outlist, [None], genes


def gene_isoform_check(
        proteins: dict[str, protein_package.Protein],
        genes: dict[str, gene_package.Gene],
        indel_map,
        config: Config) -> dict[str, gene_package.Gene]:
    if config.verbosity >= 2:
        t0 = time.time()

    extra_sequence_scan_proteins = {}
    prots_to_check = []
    for gene_id in genes:
        if len(genes[gene_id].proteins) < 2:
            continue

        if gene_id == 'to_fetch':
            for protein_id in genes[gene_id].proteins:
                extra_sequence_scan_proteins[protein_id] = proteins[genes[gene_id].proteins[protein_id]]
                prots_to_check.append(genes[gene_id].proteins[protein_id])
            continue

        isoform_list = list(genes[gene_id].proteins.keys())

        for iso1 in isoform_list:
            for iso2 in isoform_list:
                if iso1 == iso2:
                    continue
                extra_sequence_scan_proteins[iso1] = proteins[genes[gene_id].proteins[iso1]]
                extra_sequence_scan_proteins[iso2] = proteins[genes[gene_id].proteins[iso2]]

    if config.verbosity >= 2:
        t1 = time.time()
        config.logger.info(f'gene_isoform_check - Part 1: {t1-t0}')

    if config.verbosity >= 3:
        config.logger.info(f'In gene_isoform_check: Number of extra_sequence_scan_proteins: {len(extra_sequence_scan_proteins)}, Number of Genes: {len(genes)}')

    if config.verbosity >= 4:
        if len(genes) < 10:
            for gene_id in genes:
                genes[gene_id].print_content()

    extra_sequence_scan_proteins, indel_map, _, removed_proteins, genes = sequenceScan(config, extra_sequence_scan_proteins, indel_map, genes = genes, skip_mm_calc = True)

    if config.verbosity >= 2:
        t2 = time.time()
        config.logger.info(f'gene_isoform_check - Part 2: {t2-t1}')

    if config.verbosity >= 4:
        if len(genes) < 10:
            for gene_id in genes:
                genes[gene_id].print_content()

    n = 0
    mean_seq_len = 0
    for prot_id in proteins:
        if prot_id in extra_sequence_scan_proteins:
            proteins[prot_id] = extra_sequence_scan_proteins[prot_id]

            if proteins[prot_id].sequence is not None:
                seq_len = len(proteins[prot_id].sequence)
                n += 1
                mean_seq_len += seq_len
            else:
                config.errorlog.add_warning(f'Sequence is None for {prot_id}')
    if n > 0:
        mean_seq_len = mean_seq_len/n
    if config.verbosity >= 2:
        config.logger.info(f'Extra sequenceScan replace {n} number of protein objects with mean sequence lenght: {mean_seq_len}')

    if config.verbosity >= 2:
        t3 = time.time()
        config.logger.info(f'gene_isoform_check - Part 3: {t3-t2}')

    for removed_prot_id in removed_proteins:
        del proteins[removed_prot_id]

    if config.verbosity >= 2:
        t4 = time.time()
        config.logger.info(f'gene_isoform_check - Part 4: {t4-t3}')


    # Serial mode calls mm_generation directly, so config travels as itself.
    mm_parallel = parallel_db_access(config)
    pieces = []
    total_cost: int = 0

    for gene_id in genes:
        if len(genes[gene_id].proteins) < 2:
            continue
        isoform_list = list(genes[gene_id].proteins.values())

        sequence_map = {}
        iso_pairs = []

        cost: int = 0

        # Prune isoforms with no protein object, once, before pairing. Also dedupe:
        # Gene.proteins maps input_id -> primary_protein_id, and since identical
        # sequences fuse onto one primary id, several inputs of a gene can now point at
        # the same protein. Without this an isoform would be paired with itself.
        present_isoforms = []
        seen_isoforms = set()
        for iso in isoform_list:
            if iso not in proteins:
                if iso in genes[gene_id].proteins:
                    del genes[gene_id].proteins[iso]
                continue
            if iso in seen_isoforms:
                continue
            seen_isoforms.add(iso)
            present_isoforms.append(iso)

        # Each unordered pair once, both directions appended. Iterating the full list
        # twice instead produced every ordered pair twice over, which doubled the
        # alignments in mm_generation and, because add_multi_mutation does not dedupe,
        # attached each isoform multi mutation to its protein twice.
        for index, iso1 in enumerate(present_isoforms):
            for iso2 in present_isoforms[index + 1:]:
                if config.verbosity >= 6:
                    config.logger.info(f'Generating Isoform-Multimutation for {iso1} and {iso2}')

                if iso1 not in sequence_map:
                    sequence_map[iso1] = proteins[iso1].sequence
                if iso2 not in sequence_map:
                    sequence_map[iso2] = proteins[iso2].sequence
                iso_pairs.append((iso1, iso2))
                iso_pairs.append((iso2, iso1))
                cost += len(sequence_map[iso1]) * len(sequence_map[iso2])

        if len(iso_pairs) == 0:
            continue

        total_cost += cost

        pieces.append((sequence_map, iso_pairs, cost))

    # At least 1: with many cores and a little isoform work total_cost // proc_n floors
    # to zero, and the dispatch loop below then divides by it.
    package_cost: int = max(1, total_cost // config.proc_n)

    # config_dump is only read by the dispatch loop, so it is built here rather than
    # before the loop above - and only when there is something to dispatch.
    #
    # This is what made Part 5 cost 30-50 seconds on inputs with no isoform work at all.
    # gene_isoform_check runs inside buildQueue, which main() calls *before* its
    # ray_init, so on a small input ray is not up yet; an unconditional ray.put fell
    # into ray's auto_init_hook and cold started an unconfigured cluster, defeating
    # SMALL_INPUT_RAY_DEFERRAL_THRESHOLD.
    if len(pieces) == 0:
        # Nothing to distribute. Staying serial also keeps Part 7 off ray.get, which
        # raises check_connected() with no ray running, even on an empty list.
        mm_parallel = False
        config_dump = config
    elif mm_parallel:
        # Configured init, so the cluster matches config rather than ray's defaults.
        ray_init(config)
        config_dump = ray.put(config)
    else:
        config_dump = config

    if config.verbosity >= 2:
        t5 = time.time()
        config.logger.info(
            f'gene_isoform_check - Part 5: {t5-t4} {total_cost=} {package_cost=} '
            f'{len(pieces)=} {mm_parallel=}'
        )

    current_package_cost = 0
    package = []

    mm_generation_remote_ids = []
    for (sequence_map, iso_pairs, cost) in pieces:
        if cost > package_cost*2:
            # Never more splits than there are pairs to split, otherwise most of them
            # are empty tasks.
            splits = min(max(1, cost // package_cost), len(iso_pairs))
            # Ceiling division, so the splits together cover every pair.
            split_size = -(-len(iso_pairs) // splits)
            for i in range(splits):
                # The end of this slice used to read i*(split_size+1), which made split
                # 0 empty and later splits overlap: pairs were silently dropped and
                # others aligned twice.
                iso_pair_split = iso_pairs[i*split_size : (i+1)*split_size]
                if len(iso_pair_split) == 0:
                    continue
                if config.verbosity >= 3:
                    config.logger.info(f'Starting new para_mm_generatio process {cost=}')
                mm_generation_remote_ids.append(para_mm_generation.remote(config_dump, pack([(sequence_map, iso_pair_split)]))
                                                          if mm_parallel else mm_generation(config_dump, pack([(sequence_map, iso_pair_split)])))
            continue

        package.append((sequence_map, iso_pairs))
        current_package_cost += cost
        if current_package_cost >= package_cost:
            if config.verbosity >= 3:
                config.logger.info(f'Starting new para_mm_generatio process {current_package_cost=}')
            mm_generation_remote_ids.append(para_mm_generation.remote(config_dump, pack(package))
                                            if mm_parallel else mm_generation(config_dump, pack(package)))
            package = []
            current_package_cost = 0

    if len(package) > 0:
        mm_generation_remote_ids.append(para_mm_generation.remote(config_dump, pack(package))
                                        if mm_parallel else mm_generation(config_dump, pack(package)))

    if config.verbosity >= 2:
        t6 = time.time()
        config.logger.info(f'gene_isoform_check - Part 6: {t6-t5} {len(mm_generation_remote_ids)=}')

    # Serial dispatch already produced the values, so there is nothing to fetch -
    # and ray.get would raise check_connected() with no ray running, even on [].
    mm_generation_remote_results = ray.get(mm_generation_remote_ids) if mm_parallel else mm_generation_remote_ids

    if config.verbosity >= 2:
        t7 = time.time()
        config.logger.info(f'gene_isoform_check - Part 7: {t7-t6}')

    isoform_mm_matrix = {}

    for mm_generation_results in mm_generation_remote_results:
        for iso1, iso2, mm_results in unpack(mm_generation_results):
            if isinstance(mm_results, str):
                config.errorlog.add_warning(f'Multi mutation generation failed for {iso1} and {iso2} with error:\n{mm_results}')
                continue

            multi_mutation = mm_results

            for indel_or_snv in multi_mutation:
                if isinstance(indel_or_snv, tuple):
                    try:
                        proteins[iso1].add_positions([indel_or_snv[0]])
                    except TypeError as err:
                        raise TypeError(f'TypeError: {err=} {iso1=} {iso2=} {indel_or_snv=}')
                else:
                    uniprot.indel_insert(config, proteins, indel_map, [indel_or_snv], iso1, create_mut_protein_object = False)

            isoform_mm_matrix[(iso1, iso2)] = multi_mutation

            if iso1 not in indel_map:
                indel_map[iso1] = {}
            indel_map_entry = indel_map[iso1]
            proteins[iso1].add_multi_mutation(multi_mutation, indel_map_entry, mut_prot_id = iso2)

    #Remove the sequences for saving memory (sequences will be refetched in the individual chunk loops)
    for prot_id in prots_to_check:
        proteins[prot_id].sequence = None

    if config.verbosity >= 2:
        t8 = time.time()
        config.logger.info(f'gene_isoform_check - Part 8: {t8-t7}')

    return genes

@ray.remote(max_calls=1)
def para_mm_generation(*args, **kwargs):
    # Thin wrapper so the body can also run inline, without ray, when
    # parallel_db_access is False (read_only or proc_n == 1).
    return mm_generation(*args, **kwargs)


def mm_generation(config, packages):
    reset_logger_for_remotes(config)
    mm_results = []
    aligner_class = globalAlignment.init_bp_aligner_class()
    for sequence_map, iso_pairs in unpack(packages):
        for iso1, iso2 in iso_pairs:
            seq_iso1 = sequence_map[iso1]
            seq_iso2 = sequence_map[iso2]
            if len(seq_iso1) == 0:
                mm_results.append((iso1, iso2, f'Couldnt call generate_multi_mutation, sequence of {iso1} is empty'))
                continue
            if len(seq_iso2) == 0:
                mm_results.append((iso1, iso2, f'Couldnt call generate_multi_mutation, sequence of {iso2} is empty'))
                continue
            mm_result = generate_multi_mutation(seq_iso1, seq_iso2, config, aligner_class = aligner_class)
            if isinstance(mm_result, str):
                mm_results.append((iso1, iso2, mm_result))
            else:
                mm_results.append((iso1, iso2, mm_result[1]))
    return pack(mm_results)

def nToAA(seq):
    i = 1
    triple = ''
    aa_seq = ''
    for char in list(seq):
        triple += char
        if i % 3 == 0:
            if triple in codon_consts.CODONS:
                aa = codon_consts.CODONS[triple] if triple not in codon_consts.STOP_CODONS else ''
            aa_seq += aa
            triple = ''
        i += 1
    return aa_seq


def fetch_missing_sequences(proteins: protein_package.Proteins, config: Config):
    """Fetch what sequenceScan left without a sequence, before protCheck deduplicates by it.

    This used to happen after protCheck, as part of getSequences, so a protein that reached
    protCheck without its sequence was never deduplicated and got a new Protein row even next
    to a stored one of the same sequence.
    """
    t0 = time.time()

    n_missing = uniprot.fetch_missing_sequences(proteins, config)

    if config.verbosity >= 2:
        config.logger.info(f'Time for fetching {n_missing} missing sequences: {time.time() - t0}')


def add_protein_infos(proteins: protein_package.Proteins, config: Config):
    """GO terms and pathways, then the sequences and infos written to the new Protein rows.

    After protCheck, unlike the sequence fetch: the GO terms are looked up by the final primary
    ids, and addProtInfos needs the database ids.
    """
    t0 = time.time()

    uniprot.add_go_terms_and_pathways(proteins, config)

    t1 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f"Time for add_protein_infos Part 1: {t1 - t0}")

    if not config.read_only_mode:
        database.addProtInfos(proteins, config)

    t2 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f"Time for add_protein_infos Part 2: {t2 - t1}")


# Largest (pdb code x chain code) key space map_hits_to_structures builds dense lookup
# tables for: a bool and an int32 entry per key, ~320 MB at the cap. The PDB store needs
# ~200k pdb ids x ~60 chain letters; anything larger falls back to np.unique.
DENSE_HIT_KEY_TABLE_CAP = 1 << 26


def map_hits_to_structures(
        hit_store: MMseqs2.HitStore,
        pdb_table: list[str],
        proteins: protein_package.Proteins) -> np.ndarray:
    """Index into proteins.annotated_structures of every hit row of the store, int32.

    Creates the Structure of every (pdb_id, chain) and the Complex of every pdb_id the
    store names, if they do not exist yet, exactly as the per-hit object loop did. The
    Python work is one step per distinct (pdb_id, chain) of the store, none per hit.
    """
    pdb_codes = hit_store.pdb_codes
    chains = hit_store.chains
    if len(pdb_codes) == 0:
        return np.empty(0, dtype=np.int32)

    # Chain bytes compacted to the letters that occur, so the key space is
    # n_pdb x n_chain_letters rather than n_pdb x 256.
    chain_values = np.flatnonzero(np.bincount(chains, minlength=256))
    chain_code = np.zeros(256, dtype=np.int64)
    chain_code[chain_values] = np.arange(len(chain_values))
    n_chain_codes = len(chain_values)
    keys = pdb_codes.astype(np.int64) * n_chain_codes + chain_code[chains]

    key_space = len(pdb_table) * n_chain_codes
    if key_space <= DENSE_HIT_KEY_TABLE_CAP:
        present = np.zeros(key_space, dtype=bool)
        present[keys] = True
        unique_keys = np.flatnonzero(present)
        del present
        slot_of_key = np.empty(key_space, dtype=np.int32)
        slot_of_key[unique_keys] = np.arange(len(unique_keys), dtype=np.int32)
        local_index = slot_of_key[keys]
        del slot_of_key
    else:
        unique_keys, local_index = np.unique(keys, return_inverse=True)
        local_index = local_index.reshape(-1)
    del keys

    # The hit's target_len is mmseqs' length of the target sequence, so every row of one
    # (pdb_id, chain) carries the same value and any of them will do.
    lengths_of_key = np.empty(len(unique_keys), dtype=np.int32)
    lengths_of_key[local_index] = hit_store.target_len

    registry = proteins.annotated_structures
    registry_index = proteins.annotated_structure_index
    all_structures = proteins.structures
    all_complexes_map = proteins.complexes

    key_pdb_codes = (unique_keys // n_chain_codes).tolist()
    key_chains = [chr(c) for c in chain_values[unique_keys % n_chain_codes].tolist()]
    key_lengths = lengths_of_key.tolist()
    local_to_global = np.empty(len(unique_keys), dtype=np.int32)

    # Structure and Complex carry GC tracked containers, and a whole-PDB store creates
    # ~500k of each. Nothing here forms reference cycles; see the gc.freeze in
    # autoTemplateSelection.
    gc_was_enabled: bool = gc.isenabled()
    gc.disable()
    try:
        last_pdb_code = -1
        pdb_id = None
        known_chains = None
        # Keys are sorted, so the chains of one pdb are adjacent.
        for j, pdb_code in enumerate(key_pdb_codes):
            if pdb_code != last_pdb_code:
                last_pdb_code = pdb_code
                pdb_id = pdb_table[pdb_code]
                known_chains = all_structures.get(pdb_id)
                if pdb_id not in all_complexes_map:
                    all_complexes_map[pdb_id] = complex_package.Complex(pdb_id, None)

            chain = key_chains[j]
            if known_chains is None or chain not in known_chains:
                struct = structure_package.Structure(pdb_id, chain, seq_len = key_lengths[j])
                if known_chains is None:
                    known_chains = {}
                    all_structures[pdb_id] = known_chains
                known_chains[chain] = struct

            key = (pdb_id, chain)
            global_index = registry_index.get(key)
            if global_index is None:
                global_index = len(registry)
                registry_index[key] = global_index
                registry.append(key)
            local_to_global[j] = global_index
    finally:
        if gc_was_enabled:
            gc.enable()

    return local_to_global[local_index]


def annotated_structure_columns(proteins: protein_package.Proteins) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(sequence length, internal database id, multi chain) of every annotated_structures entry.

    Columns in registry order, so a protein's structure_idx gathers its hits' values.
    -1 marks an unknown length and 0 a missing database id, the value
    Structure.get_database_id reports for it. multi chain is whether the structure's
    complex has more than one stored protein chain, which an interface needs; read
    after the complex parsing, which is what fills in the chain types of new complexes.
    """
    structures = proteins.annotated_structure_objects()
    lengths = np.fromiter(
        (-1 if length is None else length for length in (structure.get_seq_len() for structure in structures)),
        dtype=np.int64, count=len(structures))
    database_ids = np.fromiter(
        (structure.get_database_id() for structure in structures),
        dtype=np.int64, count=len(structures))
    # Once per complex rather than once per chain of it.
    chain_counts: dict[str, int] = {}
    multi_chain = np.zeros(len(structures), dtype=bool)
    for index, (pdb_id, _) in enumerate(proteins.annotated_structures):
        chain_count = chain_counts.get(pdb_id)
        if chain_count is None:
            chain_count = proteins.protein_chain_count(pdb_id)
            chain_counts[pdb_id] = chain_count
        multi_chain[index] = chain_count > 1
    return lengths, database_ids, multi_chain


def autoTemplateSelection(config: Config, proteins: protein_package.Proteins):

    if config.verbosity >= 1:
        config.logger.info('Sequence search with mmseqs2')

    #filtering_dump = None

    t_sel = time.time()
    selection_times = []

    search_results: list[tuple[MMseqs2.HitStore, bool]] | None
    search_results = MMseqs2.search(proteins, config)

    t_sel = add_to_times(selection_times, t_sel, 'MMseqs2 search')

    if search_results is not None:
        total_loop_times= []

        # Consumed destructively. search_results holds one hit store per search
        # database and only one of them is ever being converted, but keeping the
        # list alive kept all of them alive for the whole function.
        # pop(0) rather than pop() so the databases are still processed in order.
        while len(search_results) > 0:
            search_result = search_results.pop(0)
            if len(search_result) == 0:
                continue
            loop_times = []
            t_loop = time.time()

            hit_store: MMseqs2.HitStore
            hit_store, is_model_db = search_result
            del search_result

            # The structure filtering that used to run here (resolution threshold, the
            # oligomer map from pdbParser) read the hits' seq_id, coverage, oligomers and
            # the pdb id set. The store carries none of them any more - only what is
            # still read, see MMseqs2.HitStore - so reviving that step means restoring
            # those columns in parseHits first.

            # One registry index per hit, and the Structure and Complex objects of every
            # (pdb, chain) the store names. No per-hit object is built: see
            # map_hits_to_structures and Protein.structure_annotations.
            pdb_table = MMseqs2.decode_pdb_table(hit_store)
            store_structure_idx = map_hits_to_structures(hit_store, pdb_table, proteins)

            registry = proteins.annotated_structures
            row_offsets = hit_store.row_offsets.tolist()
            log_annotations = config.verbosity >= 7

            for prot_index, prot_id in enumerate(MMseqs2.decode_prot_ids(hit_store)):
                protein = proteins.protein_map[prot_id]
                low = row_offsets[prot_index]
                high = row_offsets[prot_index + 1]
                if high == low:
                    continue
                # A copy, so the protein does not pin the whole store's array.
                protein.add_structure_idx(store_structure_idx[low:high].copy(), registry)

                if log_annotations:
                    for index in store_structure_idx[low:high].tolist():
                        pdb_id, chain = registry[index]
                        config.logger.info(f'Adding structural annotation: {prot_id} -> {pdb_id}:{chain}')

            del pdb_table, store_structure_idx, row_offsets
            del hit_store

            # What is built here stays alive for the rest of the chunk - one Structure per
            # chain and one Complex per pdb, each with GC tracked containers - and gen-2
            # collections during alignment and annotation would each walk all of it and
            # reclaim nothing. freeze moves it into the permanent generation so those
            # collections skip it: measured at 1.14s per 3M objects for a single gen-2
            # pass.
            #
            # Safe because the graph is acyclic - Structure and Complex reference strings
            # and numbers, never each other - so refcounting still reclaims all of it.
            # main() unfreezes at the end of the chunk, before clean_chunk_memory, whose
            # collections must be able to see this again.
            gc.freeze()

            t_loop = add_to_times(loop_times, t_loop, 'create object datastructs')

            total_loop_times = aggregate_times(total_loop_times, loop_times)

        selection_times.append(total_loop_times)

    if config.verbosity >= 2:
        print_times(selection_times, 'Template Selection', logger=config.logger)

    #del filtering_dump

@ray.remote(max_calls=1)
def filter_structures_wrapper(chunk, filtering_dump):

    pdb_path = filtering_dump
    
    return filter_structures(chunk, pdb_path)

def filter_structures(chunk, pdb_path):
    try:
        outs = []
        for pdb_id in chunk:
            resolution, homomer_dict, chain_names = pdbParser.getInfo(pdb_id, pdb_path)
            outs.append((pdb_id, resolution, homomer_dict, chain_names))
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        outs = f'Filter structure error:\n{e}\n{f}\n{g}'

    return outs

def package_alignment_processes(
        prots_todo: list[str],
        cost_map: dict[str, int],
        optimal_chunk_cost: int,
        proteins: protein_package.Proteins,
        protein_packages: dict[str, tuple[tuple, np.ndarray]],
        alignment_results,
        started_processes,
        chunk,
        chunk_cost,
        config: Config,
        mapping_dump: ray.ObjectRef[Config],
        sequence_storage,
        pos_db_storage,
        total_processes: int,
        db_id_vectors: np.ndarray[int],
        base_info_stores: list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]],
        structure_columns: tuple[np.ndarray, np.ndarray, np.ndarray],
        split_proteins: set[str] | None = None
        ):
    """split_proteins collects the proteins whose hits get spread over several alignment
    tasks; their package entries say so (complete=False), so align() leaves the interface
    limit on their backmaps to database.prune_interface_backmaps."""
    done = []
    package_cost = None
    if split_proteins is None:
        split_proteins = set()

    # From annotated_structure_columns: a package is a protein's structure_idx, and the
    # alignment task gets the database ids and multi chain flags gathered from it.
    structure_len, structure_db_id, structure_multi_chain = structure_columns

    if config.verbosity >= 4:
        config.logger.info(f'Call of package_alignment_processes: Number of proteins: {len(prots_todo)}, {optimal_chunk_cost=:_}, {chunk_cost=:_} {started_processes=} {total_processes=}')

    for pos, prot_id in enumerate(prots_todo):
        package_cost = cost_map[prot_id]
        if package_cost > optimal_chunk_cost:  #packages that are greater than the optimal_chunk_cost get split
            seq_len = len(proteins.get_sequence(prot_id))
            (prot_specific_mapping_dump, structure_idx) = protein_packages[prot_id]

            lengths = structure_len[structure_idx]
            costs = seq_len * np.where(lengths >= 0, lengths, seq_len)
            cumulative_costs = np.cumsum(costs)
            # A hit goes into this package while the cost collected before it is still
            # below optimal_chunk_cost, the rule of the per-hit loop this replaced.
            # Costs are never negative, so those hits are a prefix.
            n_split = int(np.searchsorted(cumulative_costs - costs, optimal_chunk_cost, side='left'))
            split_cost = int(cumulative_costs[n_split - 1]) if n_split > 0 else 0
            split_a = structure_idx[:n_split]
            split_b = structure_idx[n_split:]

            if config.verbosity >= 5:
                config.logger.info(f'Splitting protein alignment: {prot_id=}, {package_cost=:_} {split_cost=:_} {len(split_a)=} {len(split_b)=}')
            split_proteins.add(prot_id)
            if parallel_db_access(config):
                alignment_results.append(align_remote_wrapper.remote(
                    mapping_dump, sequence_storage,
                    pos_db_storage, db_id_vectors, base_info_stores,
                    [(prot_specific_mapping_dump, structure_db_id[split_a], structure_multi_chain[split_a], False)]))
                started_processes += 1
                total_processes += 1
            else:
                alignment_results.append(align(config, chunk, db_id_vectors, base_info_stores))
            if config.verbosity >= 5:
                config.logger.info(f'A, Start alignment package with cost: {split_cost} {len(alignment_results)=}')

            if len(split_b) > 0:
                protein_packages[prot_id] = (prot_specific_mapping_dump, split_b)
                cost_map[prot_id] = package_cost - split_cost
            else:
                done.append(pos)

        elif chunk_cost + package_cost <= optimal_chunk_cost:
            (prot_specific_mapping_dump, structure_idx) = protein_packages[prot_id]
            chunk.append((prot_specific_mapping_dump, structure_db_id[structure_idx], structure_multi_chain[structure_idx],
                          prot_id not in split_proteins))
            done.append(pos)
            chunk_cost += package_cost
            if chunk_cost >= optimal_chunk_cost:
                if parallel_db_access(config):
                    alignment_results.append(align_remote_wrapper.remote(
                        mapping_dump, sequence_storage, 
                        pos_db_storage, db_id_vectors, base_info_stores, 
                        chunk))
                    started_processes += 1
                    total_processes += 1
                else:
                    alignment_results.append(align(config, chunk, db_id_vectors, base_info_stores))
                if config.verbosity >= 5:
                    config.logger.info(f'B, Start alignment package with cost: {chunk_cost} {len(alignment_results)=}')
                chunk = []
                chunk_cost = 0
        
        if started_processes >= config.proc_n:
            break

    done.reverse()
    for pos in done:
        del prots_todo[pos]

    if (len(chunk) > 0) and len(prots_todo) == 0:
        if parallel_db_access(config):
            alignment_results.append(align_remote_wrapper.remote(
                mapping_dump, sequence_storage,
                pos_db_storage, db_id_vectors, base_info_stores,
                chunk))
            started_processes += 1
            total_processes += 1
        else:
            alignment_results.append(align(config, chunk, db_id_vectors, base_info_stores))
        if config.verbosity >= 5:
            config.logger.info(f'C, Start alignment package with cost: {chunk_cost} {len(alignment_results)=}')
        chunk = []
        chunk_cost = 0

    if (len(chunk) > 0) and package_cost is not None:
        if (chunk_cost + package_cost > optimal_chunk_cost)  and (started_processes < config.proc_n):
            if parallel_db_access(config):
                alignment_results.append(align_remote_wrapper.remote(
                    mapping_dump, sequence_storage, 
                    pos_db_storage, db_id_vectors, base_info_stores,
                    chunk))
                started_processes += 1
                total_processes += 1
            else:
                alignment_results.append(align(config, chunk, db_id_vectors, base_info_stores))
            if config.verbosity >= 5:
                config.logger.info(f'D, Start alignment package with cost: {chunk_cost} {len(alignment_results)=}')
            chunk = []
            chunk_cost = 0

    return alignment_results, started_processes, protein_packages, cost_map, chunk, chunk_cost, prots_todo, total_processes


def process_align_outs(
        align_outs,
        warn_map,
        started_processes,
        selection_totals: list[int] | None = None):
    """selection_totals, when given, sums what align() reports about the interface limit:
    [proteins selected in the alignment tasks, backmaps withheld]."""

    total_integration_times = []

    for out_chunks in align_outs:
        loop_times = []
        t_loop = time.time()

        package_runtimes, selection_counts = out_chunks
        if selection_totals is not None:
            selection_totals[0] += selection_counts[0]
            selection_totals[1] += selection_counts[1]

        t_loop = add_to_times(loop_times, t_loop, label = 'Unpack out_chunks')

        loop_times.append(package_runtimes)

        started_processes -= 1
        
        t_loop = add_to_times(loop_times, t_loop, label = 'process align outs loop')
        total_integration_times = aggregate_times(total_integration_times, loop_times)
        
    return (warn_map, total_integration_times, started_processes)


def paraAlignment(
        config: Config,
        proteins: protein_package.Proteins,
        skip_inserts: bool=False,
        indel_analysis_follow_up: bool=False,
        get_all_alignments: bool=False):
    indel_analysis_follow_up = indel_analysis_follow_up or get_all_alignments

    t_al = time.time()
    alignment_times = []

    #database.getAlignments(proteins, config, get_all_alignments=get_all_alignments)
    database.prepare_alignments(proteins, config)

    t_al = add_to_times(alignment_times, t_al, 'prepare_alignments') #0
        
    prot_ids = list(proteins.get_protein_ids())

    if parallel_db_access(config):
        mapping_dump = ray.put(config)
    else:
        mapping_dump = None

    #mapping_results = []

    #task_list: list[tuple[str, dict, int, list]] = []

    total_times = []

    #structure_insertion_list: dict[str, set[str]] = {}

    all_complexes = set()

    # The hits of unstored proteins, as one mask over annotated_structures rather than
    # a pass per hit: a structure that is not stored yet needs its complex parsed.
    needed_structures = np.zeros(len(proteins.annotated_structures), dtype=bool)

    for prot_id in prot_ids:

        if not proteins.is_protein_stored(prot_id):
            structure_idx = proteins.protein_map[prot_id].structure_idx
            if structure_idx is not None:
                needed_structures[structure_idx] = True
            continue

        if prot_id not in proteins.indels:
            if proteins.is_completely_stored(prot_id) and not get_all_alignments:
                continue

        times = []
        ti = time.time()
        # Only stored annotations count below, and index-backed hits never are.
        annotation_list = proteins.protein_map[prot_id].explicit_annotations()
        ti = add_to_times(times, ti)

        #prot_db_id = proteins.get_protein_db_id(prot_id)

        total_subtimes = []
        #struct_tasks: list[tuple[str, str, str, str, str]] = []
        for annotation in annotation_list:
            sub_times = []
            ts = time.time()

            pdb_id, chain = annotation.pdb_id, annotation.chain

            if not annotation.stored:
                continue
            else:
                if chain not in proteins.structures[pdb_id]:
                    all_complexes.add(pdb_id)
                elif not proteins.structures[pdb_id][chain].stored:
                    all_complexes.add(pdb_id)
                """
                if pdb_id not in structure_insertion_list:
                    structure_insertion_list[pdb_id] = set([chain])
                else:
                    structure_insertion_list[pdb_id].add(chain)
                """
                    
            #structure_db_id = proteins.get_structure_db_id(pdb_id, chain)
            ts = add_to_times(sub_times, ts)
            #try:
            #    struct_anno_obj = proteins[prot_id].structure_annotations[(pdb_id, chain)]
            #except KeyError as e:
            #    config.errorlog.add_error(f'KeyError {e}, no structural annotation object found: {prot_id} {pdb_id} {chain}')
            #    continue

            #if indel_analysis_follow_up:
            #    packed_alignment = struct_anno_obj.get_alignment()
            #else:
            #    packed_alignment = struct_anno_obj.pop_alignment()

            ts = add_to_times(sub_times, ts)
            #struct_tasks.append((pdb_id, chain, structure_db_id, packed_alignment))
            
            ts = add_to_times(sub_times, ts)

            total_subtimes = aggregate_times(total_subtimes, sub_times)
        
        #task_list.append((prot_id, prot_db_id, struct_tasks))

        times.append(total_subtimes)
        ti = add_to_times(times, ti)

        total_times = aggregate_times(total_times, times)

    if needed_structures.any():
        # After prepare_alignments, so structureCheck has marked the stored ones.
        structure_stored = np.fromiter(
            (structure.stored for structure in proteins.annotated_structure_objects()),
            dtype=bool, count=len(needed_structures))
        for index in np.flatnonzero(needed_structures & ~structure_stored).tolist():
            all_complexes.add(proteins.annotated_structures[index][0])
    del needed_structures

    t_al = add_to_times(alignment_times, t_al, f'{len(all_complexes)=}') #1
    if config.verbosity >= 2:
        print_times(total_times, label='process getAlignment prot_ids', logger=config.logger)

    stored_ligands = retrieve_stored_ligands(config)

    t_al = add_to_times(alignment_times, t_al, f'retrieve_stored_ligands {len(stored_ligands)=}')

    ligand_db = pdbParser.parseLigandDB(config.smiles_path, config.inchi_path)

    t_al = add_to_times(alignment_times, t_al, 'parseLigandDB')

    # An ObjectRef only when there is a ray to put it in; serial dispatch passes the
    # tuple itself, which parse_complexes_body unpacks identically.
    complex_parallel = parallel_db_access(config)
    ligand_dump = ray.put((stored_ligands, ligand_db)) if complex_parallel else (stored_ligands, ligand_db)

    t_al = add_to_times(alignment_times, t_al, 'put ligandDB')

    small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(config.proc_n, len(all_complexes))

    chunk = []
    n = 0
    parse_complex_process = []
    for struct_id in all_complexes:
        chunk.append(struct_id)
        if n < n_of_big_chunks:
            if len(chunk) == big_chunksize:
                parse_complex_process.append(parse_complexes.remote(mapping_dump, chunk, ligand_dump)
                                             if complex_parallel else parse_complexes_body(config, chunk, ligand_dump))

                chunk = []
                n += 1
        elif n <= n_of_big_chunks + n_of_small_chunks:
            if len(chunk) == small_chunksize:
                parse_complex_process.append(parse_complexes.remote(mapping_dump, chunk, ligand_dump)
                                             if complex_parallel else parse_complexes_body(config, chunk, ligand_dump))
                chunk = []
                n += 1
        else:
            config.logger.info('Warning: too few chunks')

    if len(chunk) != 0:
        parse_complex_process.append(parse_complexes.remote(mapping_dump, chunk, ligand_dump)
                                             if complex_parallel else parse_complexes_body(config, chunk, ligand_dump))

    # Serial dispatch already produced the values; ray.get would raise
    # check_connected() with no ray running, even on an empty list.
    complex_parse_results = ray.get(parse_complex_process) if complex_parallel else parse_complex_process

    t_al = add_to_times(alignment_times, t_al, 'parse complexes')

    internal_structure_db_ids = []
    for outs in complex_parse_results:
        structure_db_id_lists: list[list[tuple[int, str, str]]]
        structure_db_id_lists, atom_count_list, compl_id_struct_id_tuples = outs
        for db_nr, structure_db_id_list in enumerate(structure_db_id_lists):
            for s_id, pdb_id, chain in structure_db_id_list:
                internal_structure_db_id = 100*s_id + db_nr
                internal_structure_db_ids.append(internal_structure_db_id)
                #config.logger.info(f'{pdb_id=} {chain=} {internal_structure_db_id=}')
                
                if not proteins.contains_structure(pdb_id, chain):
                    
                    struct = structure_package.Structure(pdb_id, chain, database_id=internal_structure_db_id)
                    proteins.add_structure(pdb_id, chain, struct)
                else:
                    proteins.set_structure_db_id(pdb_id, chain, internal_structure_db_id)
                proteins.set_structure_stored(pdb_id, chain, True)

        for struct_id, atom_count, chain_type_map in atom_count_list:
            proteins.set_atom_count(struct_id, atom_count)
            proteins.complexes[struct_id].chains = chain_type_map

        for compl_db_id, struct_id, resolution in compl_id_struct_id_tuples:
            try:
                proteins.complexes[struct_id].database_id = compl_db_id
                proteins.complexes[struct_id].resolution = resolution
            except KeyError:
                continue

    if len(internal_structure_db_ids) > 0:
        misdi = max(internal_structure_db_ids)
        if misdi > proteins.max_internal_struct_db_id:
            proteins.max_internal_struct_db_id = misdi

    t_al = add_to_times(alignment_times, t_al, 'process complex parse results')

    total_cost = 0
    cost_map = {}

    # Per protein: its mapping dump and the structure_idx of the hits still to align.
    protein_packages: dict[str, tuple[tuple, np.ndarray]] = {}

    total_loop_times = []

    all_seqs = []
    pos_db_matrix = []
    current_seq_index = 0
    seq_index_dict: dict[int, tuple[int, int]] = {}

    weight_cost_by_structure_length = config.proc_n > 1

    # Every hit is an index into annotated_structures, so what the alignment needs per
    # hit - the structure's length for the cost, its database id for the package - is a
    # per-structure column, read once here now that structureCheck and the complex
    # parsing have set every database id. This used to be a python pass per hit: on a
    # whole-PDB chunk ~100M of them took ~30 min, most of it spent on the tuple and the
    # two list slots it kept for each hit.
    structure_len, structure_db_id, structure_multi_chain = annotated_structure_columns(proteins)
    structure_columns = (structure_len, structure_db_id, structure_multi_chain)

    # Multi chain hits per protein. align() applies the interface limit itself to a protein
    # whose hits all go into one task. Only a protein spread over several tasks
    # (split_proteins) has no task that sees its whole sequence identity ranking, so its
    # backmaps are thinned out once all alignments are written, by
    # database.prune_interface_backmaps - when it has more multi chain hits than
    # interface_annotations_keep_all. Stored proteins are never searched, so all Alignment
    # rows of these proteins come from this run, and the hit count bounds how many there are.
    multi_chain_hits: dict[str, int] = {}
    split_proteins: set[str] = set()
    selection_totals = [0, 0]
    keep_all_annotations, _, _ = interface_package.interface_annotation_limits(config)

    # The internal structure db ids some package needs, addressed by the id itself.
    # structureCheck and the complex parsing above bound them by
    # max_internal_struct_db_id; the max() only guards that invariant.
    max_structure_db_id = int(structure_db_id.max()) if len(structure_db_id) > 0 else 0
    structure_db_id_mask = np.zeros(max(proteins.max_internal_struct_db_id, max_structure_db_id) + 1, dtype=bool)

    for prot_id in prot_ids:
        loop_times = []
        t_loop = time.time()

        # Index-backed hits are never stored, which is what the per-annotation stored
        # check here used to select. The stored annotation objects a protein may carry
        # are exactly the ones that were skipped.
        structure_idx = proteins.protein_map[prot_id].structure_idx

        seq = proteins.get_sequence(prot_id)
        try:
            seq_len = len(seq)
        except TypeError:
            continue

        t_loop = add_to_times(loop_times, t_loop, label = 'get annotations')

        if structure_idx is not None and len(structure_idx) > 0:
            if weight_cost_by_structure_length:
                lengths = structure_len[structure_idx]
                protein_cost = int((seq_len * np.where(lengths >= 0, lengths, seq_len)).sum())
            else:
                protein_cost = 0
            cost_map[prot_id] = cost_map.get(prot_id, 0) + protein_cost
            structure_db_id_mask[structure_db_id[structure_idx]] = True

        t_loop = add_to_times(loop_times, t_loop, label = 'build cost_map')

        if prot_id not in cost_map:
            total_loop_times = aggregate_times(total_loop_times, loop_times)
            t_loop = add_to_times(loop_times, t_loop, label = 'build pos_db_ids')
            t_loop = add_to_times(loop_times, t_loop, label = 'numpyfy pos_db_ids')
            t_loop = add_to_times(loop_times, t_loop, label = 'make mapping dump')
            total_loop_times = aggregate_times(total_loop_times, loop_times)
            continue

        pos_db_ids: list[int] = [-1]
        for pos in proteins[prot_id].get_position_ids():
            pos_db_ids.append(proteins[prot_id].positions[pos].database_id)

        #config.logger.info(f'{prot_id=} {proteins.get_protein_db_id(prot_id)=} {pos_db_ids=}')

        t_loop = add_to_times(loop_times, t_loop, label = 'build pos_db_ids')

        pos_db_ids = np.array(pos_db_ids, dtype=np.int32)

        t_loop = add_to_times(loop_times, t_loop, label = 'numpyfy pos_db_ids')

        if parallel_db_access(config):
            prot_specific_mapping_dump = (prot_id, current_seq_index, proteins.get_protein_db_id(prot_id))
            seq_index_dict[proteins.get_protein_db_id(prot_id)] = current_seq_index, len(seq)
            all_seqs.append(seq)
            pos_db_matrix.append(pos_db_ids)
            current_seq_index += 1

        else:
            prot_specific_mapping_dump = (prot_id, seq, proteins.get_protein_db_id(prot_id), pos_db_ids)

        protein_packages[prot_id] = (prot_specific_mapping_dump, structure_idx)

        if structure_idx is not None:
            multi_chain_hits[prot_id] = int(structure_multi_chain[structure_idx].sum())

        total_cost += cost_map[prot_id]

        t_loop = add_to_times(loop_times, t_loop, label = 'make mapping dump')
        total_loop_times = aggregate_times(total_loop_times, loop_times)

    #config.logger.info(f'{protein_packages=}')
    #config.logger.info(f'{pos_db_matrix=}')

    alignment_times.append(total_loop_times)

    t_al = add_to_times(alignment_times, t_al, label = f'remote preparations {len(prot_ids)=}')

    if parallel_db_access(config):
        if len(pos_db_matrix) > 0:
            # Flat concatenation plus offsets rather than a matrix padded to the
            # longest protein: a single titin used to widen every row to 34k
            # entries, so the padding cost gigabytes that held nothing.
            sequence_storage = ray.put(build_sequence_store(all_seqs))
            pos_db_storage = ray.put(build_ragged_int_store(pos_db_matrix, dtype=np.int32))
        else:
            sequence_storage = None
            pos_db_storage = None
    else:
        sequence_storage = None
        pos_db_storage = None

    proteins.seq_store = (sequence_storage, seq_index_dict)

    t_al = add_to_times(alignment_times, t_al, label = 'Prepare Sequence Storage')

    db_id_vectors, base_info_stores, base_info_times = prepare_struct_base_info_store(structure_db_id_mask, config)
    del structure_db_id_mask

    if parallel_db_access(config):
        # One row per structure, so ~5.8 MB at the 180k structures a 500 protein run
        # reaches, and every alignment task needs all of it. Passed by value it was
        # serialised into each of the ~128 task specs and held as its own object store
        # entry per task - about 740 MB of copies of one array. Put once, the tasks
        # share a single copy. package_alignment_processes only sends this to a remote
        # call when proc_n > 1 and to a local align() otherwise, and ray resolves a
        # top level ObjectRef argument back to the array on the worker side, so the
        # callee sees the array either way (mapping_dump is already used like this).
        db_id_vectors = ray.put(db_id_vectors)

    alignment_times.append(base_info_times)

    t_al = add_to_times(alignment_times, t_al, label = 'Prepare Structure Storage')

    min_optimal_cost = 5_000_000 #too small jobs generate overhead
    max_optimal_cost = 10_000_000_000 #too large jobs can lead to ray.get to freeze
    optimal_chunk_cost = max([min_optimal_cost, min([(total_cost// config.proc_n) + 1, max_optimal_cost])])

    if config.verbosity >= 3:
        config.logger.info(f'Total cost for alignment: {total_cost=:_} {optimal_chunk_cost=:_}')

    alignment_results = []

    chunk = []
    chunk_cost = 0

    prots_todo: list[str] = list(protein_packages.keys())

    started_processes = 0
    total_processes = 0

    t_al = add_to_times(alignment_times, t_al, label = f'extend prots_todo {len(prots_todo)=}') #6

    while started_processes < config.proc_n and len(prots_todo) > 0:
        (alignment_results, started_processes,
         protein_packages, cost_map, chunk, 
         chunk_cost, prots_todo, total_processes) = package_alignment_processes(prots_todo, cost_map, optimal_chunk_cost,
                                                                                    proteins, protein_packages, alignment_results,
                                                                                    started_processes, chunk, chunk_cost, config,
                                                                                    mapping_dump, sequence_storage, pos_db_storage,
                                                                                    total_processes, db_id_vectors, base_info_stores,
                                                                                    structure_columns,
                                                                                    split_proteins
                                                                                )

    t_al = add_to_times(alignment_times, t_al, label = 'Initial pack alignment processes') #7

    warn_map = set()

    if parallel_db_access(config):

        total_loop_times = []
        while True:
            loop_times = []
            t_loop = time.time()

            if started_processes >= config.proc_n or len(prots_todo) == 0:
                if config.verbosity >= 5:
                    config.logger.info(f'Ray wait without timeout {started_processes=}')

                ready, not_ready = ray.wait(alignment_results)
            else:
                if config.verbosity >= 5:
                    config.logger.info('Ray wait with timeout')
                ready, not_ready = ray.wait(alignment_results, timeout = 0.01)

            t_loop = add_to_times(loop_times, t_loop, label = 'alignment wait')

            if len(ready) > 0:
                align_outs = ray.get(ready)
                
                t_loop = add_to_times(loop_times, t_loop, label = 'ray.get')

                (warn_map, package_runtimes, started_processes) = process_align_outs(align_outs, warn_map, started_processes, selection_totals)

                loop_times.append(package_runtimes)
                t_loop = add_to_times(loop_times, t_loop, label = 'process_align_outs')
            else:
                t_loop = add_to_times(loop_times, t_loop, label = 'ray.get')
                loop_times.append([])
                t_loop = add_to_times(loop_times, t_loop, label = 'process_align_outs')

            new_alignment_results = []
            if len(prots_todo) > 0 or len(chunk) > 0:
                (new_alignment_results, started_processes,
                 protein_packages, cost_map, chunk, 
                 chunk_cost, prots_todo, total_processes) = package_alignment_processes(
                                                                                            prots_todo, cost_map, optimal_chunk_cost,
                                                                                            proteins, protein_packages, new_alignment_results,
                                                                                            started_processes, chunk, chunk_cost, config,
                                                                                            mapping_dump, sequence_storage, pos_db_storage,
                                                                                            total_processes, db_id_vectors, base_info_stores,
                                                                                            structure_columns,
                                                                                            split_proteins
                                                                                        )
                
            t_loop = add_to_times(loop_times, t_loop, label = 'package new processes in wait loop')

            total_loop_times = aggregate_times(total_loop_times, loop_times)

            alignment_results = new_alignment_results + not_ready

            if config.verbosity >= 5:
                config.logger.info(f'Proteins to do: {len(prots_todo)}, Chunk: {len(chunk)}, Pending processes: {len(alignment_results)} {started_processes=} {total_processes=}')

            if len(alignment_results) == 0 and len(prots_todo) == 0 and len(chunk) == 0:
                if started_processes != 0 or len(prots_todo) != 0:
                    config.logger.info(f'\nWarning: left alignment loop while stuff was still to do {started_processes=} {prots_todo=}\n')
                break

        alignment_times.append(total_loop_times)

    else: #serial mode omitting ray

        while True:
            align_outs = alignment_results
            (warn_map, package_runtimes, started_processes) = process_align_outs(align_outs, warn_map, started_processes, selection_totals)

            n = len(prots_todo)
            if n == 0:
                break
            m = len(protein_packages[list(prots_todo)[0]][1])

            if started_processes < config.proc_n and len(prots_todo) > 0:
                alignment_results, started_processes, protein_packages, cost_map, chunk, chunk_cost, prots_todo, total_processes = package_alignment_processes(prots_todo, cost_map, optimal_chunk_cost, proteins, protein_packages, alignment_results, started_processes, chunk, chunk_cost, config, mapping_dump, sequence_storage, pos_db_storage, total_processes, db_id_vectors, base_info_stores, structure_columns, split_proteins)
            else:
                break

            if config.verbosity >= 5:
                config.logger.info(f'Proteins to do: {len(prots_todo)}, Chunk: {len(chunk)}, Pending processes: {len(alignment_results)}')

            if n == len(prots_todo) and n ==1 and m == len(protein_packages[list(prots_todo)[0]][1]):
                config.errorlog.add_error('Alignment todo list not decreased in serial alignment')
                break

    t_al = add_to_times(alignment_times, t_al, label = f'Alignment wait loop {total_processes=}') #8

    backmap_pruning_candidates = [proteins.get_protein_db_id(prot_id) for prot_id in split_proteins
                                  if multi_chain_hits.get(prot_id, 0) > keep_all_annotations]
    if config.verbosity >= 3:
        config.logger.info(
            f'Interface limit in the alignment tasks: {selection_totals[0]} proteins thinned out, '
            f'{selection_totals[1]} backmaps not written; {len(split_proteins)} proteins split over '
            f'several tasks, {len(backmap_pruning_candidates)} of them left to prune')

    # Every alignment task has returned, so each candidate's full ranking is in the
    # database now. Fast mode writes no backmaps at all.
    if not config.fast_pdb_annotation and len(backmap_pruning_candidates) > 0:
        multi_chain_struct_ids = np.unique(structure_db_id[structure_multi_chain])
        alignment_times.append(database.prune_interface_backmaps(
            config, backmap_pruning_candidates, multi_chain_struct_ids))
        del multi_chain_struct_ids

    t_al = add_to_times(alignment_times, t_al, label = f'prune interface backmaps {len(backmap_pruning_candidates)=}')

    for prot_id in list(protein_packages.keys()):
        del protein_packages[prot_id]
    del protein_packages

    t_al = add_to_times(alignment_times, t_al, label = 'protein_packages deletion') #9

    if len(config.rare_residues) > 0 and config.verbosity >= 3:
        config.logger.info(f'######### Detected some rare residues: {config.rare_residues}')

    del mapping_dump

    t_al = add_to_times(alignment_times, t_al, label = 'mapping dump deletion') #11

    if config.verbosity >= 2:
        print_times(alignment_times, label='Alignment main process', logger=config.logger)

    return

#called (also) by indel_analysis.py
@ray.remote(max_calls=1)
def align_remote_wrapper(align_dump, sequence_storage, pos_db_storage, db_id_vectors, base_info_stores, package):

    config = align_dump
    reset_logger_for_remotes(config)
    return align(config, package, db_id_vectors, base_info_stores, remote = True, sequence_storage=sequence_storage, pos_db_storage=pos_db_storage)


def align(
        config: Config,
        package: list[tuple[tuple[str, int, int], np.ndarray[int], np.ndarray[bool], bool]],
        db_id_vector: np.ndarray[int],
        base_info_stores: list[tuple[ray.ObjectRef, ray.ObjectRef, ray.ObjectRef]],
        remote=False,
        sequence_storage=None,
        pos_db_storage=None
        ):
    
    #results = []
    align_times = []
    t_align = time.time()
    aligner_class = globalAlignment.init_bp_aligner_class()

    alignment_insertion_list = []

    #all_struct_db_ids: np.ndarray[int] = package[0][1]
    #for _, struct_db_ids in package[1:]:
    #    all_struct_db_ids = np.union1d(all_struct_db_ids, struct_db_ids)

    resolutions_stores: list[np.ndarray[float]] = []
    bytes_matrix_stores: list[np.ndarray[bytes]] = []
    id_list_stores: list[np.ndarray[bytes]] = []

    for base_info_store in base_info_stores:
        try:
            p_id_list, p_bytes_matrix, p_resolutions = base_info_store
            # deref, not ray.get: these are ObjectRefs when the base info store was
            # built under ray and the plain arrays when it was built serially.
            resolutions_stores.append(deref(p_resolutions))
            bytes_matrix_stores.append(deref(p_bytes_matrix))
            id_list_stores.append(deref(p_id_list))
        except TypeError:
            resolutions_stores.append([])
            bytes_matrix_stores.append([])
            id_list_stores.append([])
            

    t_align = add_to_times(align_times, t_align, label = f'db id union {config.option_seq_thresh=}')

    #struct_info_dict: dict[int, tuple[str, str, tuple, float]]
    #struct_info_dict, retrieval_times = get_struct_base_infos(all_struct_db_ids, config)

    #align_times.append(retrieval_times)

    struct_store_dict = {row[0]: row[1:] for row in db_id_vector}

    t_align = add_to_times(align_times, t_align, label = 'retrieve struct infos')

    # The interface limit (interface.select_interface_annotations) is applied here to the
    # proteins whose hits all came in this package, so their left out backmaps are never
    # written; a protein split over several tasks is left to database.prune_interface_backmaps.
    # The identities are ranked as the database will hand them back to para_classify, which
    # makes the same selection again: a MySQL FLOAT column holds a float32.
    interface_limits = interface_package.interface_annotation_limits(config)
    rank_as_float32 = config.db_address != '-'
    n_proteins_selected = 0
    n_backmaps_withheld = 0

    total_loop_times = []
    for (prot_specific_mapping_dump, structure_infos, multi_chain_flags, complete) in package:

        loop_times = []
        t_loop = time.time()

        if remote:
            (prot_id, seq_index, prot_db_id) = prot_specific_mapping_dump
            seq = sequence_from_store(sequence_storage, seq_index)
            pos_db_ids = ragged_row(pos_db_storage, seq_index)
        else:
            (prot_id, seq, prot_db_id, pos_db_ids) = prot_specific_mapping_dump

        if config.verbosity >= 6:
            config.logger.info(f'Alignment of {prot_id} to {len(structure_infos)} number of Structures {pos_db_ids=}\n  {seq=}')

        t_loop = add_to_times(loop_times, t_loop, label='get mapping dump')

        counts = [0 for _ in range(5)]
        # [prot_db_id, structure_db_id, seq_id, coverage, packed_alignment, db_backmap or None];
        # the backmaps are packed once it is known which are kept.
        protein_rows: list[list] = []

        total_inner_times = []
        for structure_db_id, multi_chain in zip(structure_infos, multi_chain_flags):
            if structure_db_id == 0:
                counts[0] += 1
                continue

            #structure_id, chain, (seq_res_map, template_seq), resolution = struct_info_dict[structure_db_id]

            db_nr, store_index, trailing_count = struct_store_dict[structure_db_id]
            resolution = resolutions_stores[db_nr][store_index]

            if resolution > config.option_res_thresh:
                counts[1] += 1
                continue

            inner_times = []
            t_inner = time.time()

            fused_id = id_list_stores[db_nr][store_index].decode()
            structure_id = fused_id[:-1]
            chain = fused_id[-1]
            packed_parse_data: bytes = bytes_matrix_stores[db_nr][store_index].rstrip(b'\x00')
            if trailing_count > 0:
               packed_parse_data = packed_parse_data + (b'\x00' * trailing_count)
            try:
                (seq_res_map, template_seq) = unpack(packed_parse_data)
            except BaseException as err:
                raise NameError(f'{err=}\n{structure_db_id=} {fused_id=} {packed_parse_data=}')

            t_inner = add_to_times(inner_times, t_inner, label = 'unpack data from base_info_store')

            align_out = globalAlignment.alignBioPython(config,
                                                        prot_id,
                                                        seq,
                                                        structure_id,
                                                        None,
                                                        chain,
                                                        seq_info = (seq_res_map, template_seq),
                                                        aligner_class = aligner_class)
            
            t_inner = add_to_times(inner_times, t_inner, label = 'alignBioPython')

            
            if isinstance(align_out, str):
                total_inner_times = aggregate_times(total_inner_times, inner_times)
                counts[2] += 1
                continue
            backmap: np.ndarray
            (coverage, seq_id, sub_infos, backmap, _, bp_align_times, _, _) = align_out
            inner_times.append(bp_align_times)
            
            if sub_infos is None or seq_id is None:
                total_inner_times = aggregate_times(total_inner_times, inner_times)
                counts[3] += 1
                continue

            if config.verbosity >= 7:
                config.logger.info(f'Alignment of {prot_id} to {structure_id} {chain} successful: {seq_id=} {coverage=} {sub_infos=}')


            if 100.0 * seq_id >= config.option_seq_thresh:
                # Position db ids, which stay below 2**31. This is what gets
                # concatenated into the backmap store during aggregation, where
                # it is one of the biggest arrays held at once.
                if config.fast_pdb_annotation:
                    # No backmap at all in fast mode: aggregation skips getBackmaps
                    # there, so the column is never read back, and a NULL costs
                    # nothing where an array of -1 cost one int32 per residue. The
                    # alignment itself is still written - the alignment matrix output
                    # needs it, which is why this mode no longer drops the tables.
                    db_backmap = None
                elif not multi_chain:
                    # The only protein chain of its complex: no other chain can
                    # interact with it, so its backmap would never be read. See
                    # interface.select_interface_annotations.
                    db_backmap = None
                else:
                    # pos_db_ids carries a -1 sentinel at index 0; its real ids
                    # start at 1, matching backmap's 1-based positions. Two kinds
                    # of entry have no position to point at and are folded onto
                    # that sentinel, so they come out as the -1 the aggregation
                    # side already tests for:
                    #   - backmap holds -1 where a residue aligns to a gap, and a
                    #     raw -1 index would wrap onto the last position instead
                    #     of reaching the sentinel.
                    #   - a protein only has Position objects up to its highest
                    #     requested position, so residues past that have no id.
                    #     The dense padded store this replaces returned its zero
                    #     padding for them.
                    position_indices = np.where(backmap < 0, 0, backmap)
                    position_indices[position_indices >= len(pos_db_ids)] = 0
                    db_backmap = pos_db_ids[position_indices].astype(np.int32, copy = False)

                packed_alignment = pack(sub_infos)
                protein_rows.append([prot_db_id, structure_db_id, seq_id, coverage, packed_alignment, db_backmap])

                #if is_alphafold_model(structure_id) and config.provide_msas:
                #    fetch_msa(config, structure_id, prot_id, seq, seq_id, 1)
            else:
                counts[4] += 1

            t_inner = add_to_times(inner_times, t_inner, label = 'packing alignment outs')

            total_inner_times = aggregate_times(total_inner_times, inner_times)

        n_with_backmap = sum(1 for row in protein_rows if row[5] is not None)
        if complete and n_with_backmap > interface_limits[0]:
            struct_db_ids = np.array([row[1] for row in protein_rows], dtype=np.int64)
            seq_ids = np.array([row[2] for row in protein_rows], dtype=np.float64)
            if rank_as_float32:
                seq_ids = seq_ids.astype(np.float32).astype(np.float64)
            has_backmap = np.array([row[5] is not None for row in protein_rows], dtype=bool)
            selected = interface_package.select_interface_annotations(
                struct_db_ids, seq_ids, has_backmap, *interface_limits)
            for row, keep in zip(protein_rows, selected.tolist()):
                if row[5] is not None and not keep:
                    row[5] = None
                    n_backmaps_withheld += 1
            n_proteins_selected += 1

        for row in protein_rows:
            if row[5] is not None:
                row[5] = pack(row[5])
            alignment_insertion_list.append(tuple(row))
            if len(alignment_insertion_list) == 1000:
                database.insertAlignments(alignment_insertion_list, config)
                alignment_insertion_list = []
        del protein_rows

        if config.verbosity >= 6:
            config.logger.info(f'Alignment of {prot_id} {len(alignment_insertion_list)=} {counts=}')

        loop_times.append(total_inner_times)
        t_loop = add_to_times(loop_times, t_loop, label='inner alignment core loop')
        total_loop_times = aggregate_times(total_loop_times, loop_times)

    align_times.append(total_loop_times)
    t_align = add_to_times(align_times, t_align, label = 'outer alignment core loop')

    database.insertAlignments(alignment_insertion_list, config)

    t_align = add_to_times(align_times, t_align, label = 'insertAlignments')

    return align_times, (n_proteins_selected, n_backmaps_withheld)

@ray.remote(max_calls=1)
def parse_complexes(*args, **kwargs):
    # Thin wrapper so the body can also run inline, without ray, when
    # parallel_db_access is False (read_only or proc_n == 1).
    return parse_complexes_body(*args, **kwargs)


def parse_complexes_body(config, complexes, ligand_dump):
    reset_logger_for_remotes(config)
    static_model_path=None
    structure_insertion_dict = {}
    results = []
    compl_values = []
    compl_iap = {}
    for structure_id in complexes:
        model_path = None

        # find and sets model path in given custom db 
        if config.custom_db_path:
            custom_path_list = os.listdir(config.custom_db_path)
            for model in custom_path_list:
                if model.startswith(structure_id):
                    model_path = config.custom_db_path + f"/{model}"

        #is_model = False
        elif config.model_db_active:
            ret_of_is_model = is_model(structure_id)
            if isinstance(ret_of_is_model, str):
                config.errorlog.add_error(ret_of_is_model, logger=config.logger)
            if ret_of_is_model:
                if is_alphafold_model(structure_id):
                    model_path = alphafold_model_id_to_file_path(structure_id, config)
                else:
                    model_path = complex_id_to_file_path(structure_id, config)
                #is_model = True
            else:
                model_path = static_model_path
        else:
            model_path = static_model_path

        parse_out = pdbParser.getStandardizedPdbFile(structure_id, config, model_path=model_path)

        if parse_out is None or isinstance(parse_out, str):
            continue

        (interaction_partners, chain_type_map, _, atom_count, _, _, parse_times, seq_res_maps, _, _, seqs) = parse_out
        if model_path == 'complex_homo':
            homomer_dict = {'A' : set(['A','B']), 'B': set(['A','B'])}
            resolution = 2.0
            chain_names = {'A': 'I am a model', 'B': 'Also'}
        elif model_path == 'complex_hetero':
            homomer_dict = {'A' : set(['A']), 'B': set(['B'])}
            resolution = 2.0
            chain_names = {'A': 'I am a model', 'B': 'Also'}
        elif ret_of_is_model:
            homomer_dict = {'A' : set(['A'])}
            resolution = 2.0
            chain_names = {'A': 'I am a model'}
        else:
            resolution, homomer_dict, chain_names = pdbParser.getInfo(structure_id, config.pdb_path)

        if config.verbosity >= 5:
            print(f'In parse_complexes_body: {model_path=} {ret_of_is_model=} {resolution=} {homomer_dict=}')
        
        if chain_names is None:
            chain_str = ''
        else:
            chain_str_parts = [f'{chain}:{chain_names[chain]}' for chain in chain_names]
            chain_str = ','.join(chain_str_parts)

        del_list = []
        for chain in homomer_dict:
            if chain not in chain_type_map:
                del_list.append(chain)
            else:
                del_chains = []
                for h_chain in chain_type_map[chain]:
                    if h_chain not in chain_type_map:
                        del_chains.append(h_chain)
                for h_chain in del_chains:
                    try:
                        homomer_dict[chain].remove(h_chain)
                    except KeyError:
                        pass
        for chain in del_list:
            del homomer_dict[chain]

        used_chains = set()
        parts = []
        for chain in homomer_dict:
            if chain in used_chains:
                continue
            part = chain
            used_chains.add(chain)
            for chain2 in homomer_dict[chain]:
                if chain2 == chain:
                    continue
                used_chains.add(chain2)
                part += chain2
            parts.append(part)
        homomer_str = ','.join(parts)

        compl_values.append((structure_id, resolution, chain_str, homomer_str, pack(chain_type_map)))
        
        compl_iap[structure_id] = interaction_partners

        results.append((structure_id, atom_count, chain_type_map))

        for chain in chain_type_map:
            chain_type: str = chain_type_map[chain]
            if chain_type == 'Protein':
                
                if structure_id not in structure_insertion_dict:
                    structure_insertion_dict[structure_id] = {}
                
                if chain in homomer_dict:
                    if isinstance(homomer_dict[chain], set):
                        homomers = ''.join(homomer_dict[chain])
                else:
                    homomers = chain
                if homomers[:2] == 'I ':
                    homomers = chain
                structure_insertion_dict[structure_id][chain] = ((seq_res_maps[chain], seqs[chain]), homomers)

    compl_id_struct_id_tuples = remote_insert_complexes(config, compl_values, compl_iap, ligand_dump)

    compl_id_res_dict = {struct_id: (compl_db_id, resolution) for compl_db_id, struct_id, resolution in compl_id_struct_id_tuples}

    #compl_id_store = ray.put(compl_id_struct_id_tuples)
    struct_db_id_list: list[list[tuple[int, str, str]]] = database.remoteInsertStructures(config, structure_insertion_dict, compl_id_res_dict)
    return struct_db_id_list, results, compl_id_struct_id_tuples

def debug_freeze(config: Config, stop_number: int=1):
    config.logger.info(f'\n\n{config.ray_init_counter=}\n\n')
    if config.ray_init_counter == stop_number:
        time.sleep(999999999)
    return

def set_small_input_ray_deferral(config: Config, proteins_chunks, temp_infiles) -> None:
    """Decide whether this input is small enough to postpone the ray cold start.

    A cold start is not cheap - it prestarts one python worker per cpu and measured at
    ~33s of a ~60s run (see the note in ray_utils.parallel_db_access) - and an input of
    a handful of proteins that are all stored already never needs a cluster at all. For
    such inputs everything up to and including protCheck runs serially (defer_ray_init
    makes parallel_db_access report False) and core() hands the decision to
    resolve_deferred_ray_init, which starts ray only for an input that brought new
    proteins.

    temp_infiles != [None] means buildQueue split the infile because it did not fit into
    memory, which is by definition not a small input - proteins_chunks is empty in that
    case and would otherwise count as zero proteins.
    """
    input_proteins = sum(len(protein_list) for protein_list, _ in proteins_chunks)
    config.defer_ray_init = (
        temp_infiles == [None]
        and input_proteins <= SMALL_INPUT_RAY_DEFERRAL_THRESHOLD
    )

    if config.defer_ray_init and config.verbosity >= 1:
        config.logger.info(
            f'Small input ({input_proteins} proteins): running serially until protCheck, '
            'ray gets started only if the input contains proteins that are not stored yet'
        )


def resolve_deferred_ray_init(config: Config, proteins: protein_package.Proteins) -> bool:
    """Start the ray cluster a small input deferred, if protCheck found new proteins.

    Called right after protCheck, which has just filled not_stored_ids with the proteins
    it had to insert. If none were, the deferral stays in place and the whole run - this
    chunk and the output stage after it - keeps going serially without ever starting a
    cluster. Returns whether ray was started here.
    """
    if not config.defer_ray_init:
        return False

    not_stored = len(proteins.get_not_stored_ids())
    if not_stored == 0:
        return False

    config.defer_ray_init = False

    if config.verbosity >= 1:
        config.logger.info(f'{not_stored} of the input proteins were not stored yet, starting ray')

    if not parallel_db_access(config):
        return False

    ray_init(config)
    return True


def print_locals(items, config):
    for name, size in sorted(((name, sys.getsizeof(value)) for name, value in items), key=lambda x: -x[1])[:10]:
        config.logger.info(f"{name:>30}: {sizeof_fmt(size):>8}")

def core(
        protein_list: dict[str, protein_package.Protein],
        genes,
        indels,
        multi_mutation_objects,
        config: Config,
        session
        ):
    
    if len(protein_list) == 0:
        return 0

    if parallel_db_access(config):
        ray_init(config)
        config.ray_init_counter += 1    
        #monitor_ray_store()

    background_process_MS = None

    # transform the protein map into a Proteins object
    proteins: protein_package.Proteins = protein_package.Proteins(protein_list, indels, multi_mutation_objects)

    _rss_core_base = _rss_mb()
    def _rss_log(label: str):
        if config.verbosity >= 2:
            total, anon, shmem = _rss_breakdown_mb()
            # Off Linux, or wherever /proc/self/status cannot be read, the breakdown is
            # all zeros; _rss_mb still has the getrusage fallback for the total.
            if not total:
                total = _rss_mb()
            delta = total - _rss_core_base
            config.logger.info(
                f'  [RSS] {label}: {total:.0f} MiB (anon {anon:.0f} MiB, shmem {shmem:.0f} MiB) '
                f'(delta from core start: {delta:+.0f} MiB){_object_store_note()}'
            )

    try:
        # Before protCheck: its sequence deduplication only places proteins that have one.
        fetch_missing_sequences(proteins, config)

        # Also before protCheck, so the mutants have their sequences when protCheck writes
        # their Protein rows. mutate needs only the wildtype sequence and the mutation list,
        # and protCheck leaves mutants out of both deduplication passes.
        if not config.only_wt:
            for wt_prot_id in proteins.multi_mutations:
                if config.verbosity >= 6:
                    config.logger.info(f'Mutating {len(proteins.multi_mutations[wt_prot_id])} multi mutations for: {wt_prot_id}')
                for multi_mutation in proteins.multi_mutations[wt_prot_id]:
                    multi_mutation.mutate(proteins, config)

        t1 = time.time()

        if config.verbosity >= 1:
            config.logger.info("Before protCheck")
        # check for already stored genes and make all the necessary database interactions
        # sets the fields database_id and stored in the protein objects as well as the stored and non stored ids in the proteins object

        #if not config.fast_pdb_annotation:
        protCheck_successful = database.protCheck(proteins, genes, session, config)

        if not protCheck_successful:
            return 0

        # The ray start a small input deferred (see set_small_input_ray_deferral).
        # protCheck is the first point where it is known whether the input brought any
        # new work at all.
        if resolve_deferred_ray_init(config, proteins):
            config.ray_init_counter += 1

        t2 = time.time()
        if config.verbosity >= 2:
            config.logger.info(f"Time for protCheck: {t2 - t1}")
        _rss_log('after protCheck')

        # check for already stored mutations or position twins and make all the necessary database interactions
        # sets the fields database_id and stored in the position objects

        background_process_MS = database.positionCheck(proteins, session, config)
        
        t3 = time.time()
        if config.verbosity >= 2:
            config.logger.info(f"Time for positionCheck: {t3 - t2}")
        _rss_log('after positionCheck')

        database.indelCheck(proteins, session, config)

        t4 = time.time()
        if config.verbosity >= 2:
            config.logger.info(f"Time for indelCheck: {t4 - t3}")
        _rss_log('after indelCheck')

        database.insertMultiMutations(proteins, genes, session, config)

        t5 = time.time()
        if config.verbosity >= 2:
            config.logger.info(f"Time for insertMultiMutations: {t5 - t4}")
        _rss_log('after insertMultiMutations')

        if config.verbosity >= 1:
            config.logger.info(f"Before add_protein_infos {config.afdb=}")

        add_protein_infos(proteins, config)

        t6 = time.time()
        if config.verbosity >= 2:
            config.logger.info(f"Time for add_protein_infos: {t6 - t5}")
        _rss_log('after add_protein_infos')

        indel_analysis_follow_up = len(proteins.indels) > 0 and not config.skip_indel_analysis

        # read_only keeps only proteins that are already stored - protCheck deletes the
        # rest - so MMseqs2.search returns None ("all proteins are stored already") and
        # no new structure annotation can appear. With nothing annotated there is
        # nothing to align, nothing to analyse structurally, and nothing for the
        # aggregation that paraAnnotate calls at its end. The output stage reads what it
        # needs straight out of the database, so all three stages are skipped rather
        # than run to produce empty results.
        if config.read_only_mode:
            if config.verbosity >= 1:
                config.logger.info(
                    'read_only mode: skipping template selection, alignment and '
                    'structural analysis - every protein is already stored'
                )
            amount_of_structures = 0
            t7 = t8 = t9 = time.time()
        else:
            if config.verbosity >= 1:
                config.logger.info("Before autoTemplateSelection")

            autoTemplateSelection(config, proteins)

            t7 = time.time()
            if config.verbosity >= 2:
                config.logger.info(f"Time for Template selection: {t7 - t6}")
            _rss_log('after autoTemplateSelection')

            if config.verbosity >= 1:
                config.logger.info("Before paraAlignment")
            paraAlignment(config, proteins, indel_analysis_follow_up=indel_analysis_follow_up, get_all_alignments=config.provide_msas)

            t8 = time.time()
            if config.verbosity >= 2:
                config.logger.info(f"Time for Alignment: {t8 - t7}")
            _rss_log('after paraAlignment')

            if config.verbosity >= 1:
                config.logger.info("Before paraAnnotate")

            amount_of_structures = templateFiltering.paraAnnotate(config, proteins, indel_analysis_follow_up=indel_analysis_follow_up)

            t9 = time.time()
            if config.verbosity >= 2:
                config.logger.info(f"Time for paraAnnotate: {t9 - t8}")
            _rss_log('after paraAnnotate')

        if not config.skip_indel_analysis:
            #indel_analysis.para_indel_analysis(proteins, config)
            ...

        t10 = time.time()
        if config.verbosity >= 2:
            config.logger.info(f"Time for indel analysis: {t10 - t9}")

        proteins.deconstruct()
        del proteins

        t11 = time.time()
        if config.verbosity >= 2:
            config.logger.info(f"Time for deconstruction: {t11 - t10}")

        # Probe how much RSS survives proteins.deconstruct() inside core() itself,
        # before returning to the outer loop.  The gap between this reading and
        # the next chunk's _rss_core_base is the true persistent accumulation.
        if config.verbosity >= 2:
            for _g in range(3):
                gc.collect(_g)
            if sys.platform.startswith('linux'):
                try:
                    import ctypes as _ct
                    _ct.cdll.LoadLibrary('libc.so.6').malloc_trim(0)
                except Exception:
                    pass
            _rss_after_trim = _rss_mb()
            config.logger.info(
                f'  [RSS] after deconstruct+gc+trim: {_rss_after_trim:.0f} MiB '
                f'(persistent delta from core start: {_rss_after_trim - _rss_core_base:+.0f} MiB)'
            )
            if _rss_after_trim - _rss_core_base > 1000:
                import tracemalloc as _tm
                if _tm.is_tracing():
                    _snap = _tm.take_snapshot()
                    config.logger.info('  [tracemalloc] top allocations persisting after deconstruct:')
                    for _s in _snap.statistics('lineno')[:15]:
                        config.logger.info(f'    {_s}')

        # join the background inserts
        if background_process_MS is not None:
            background_process_MS.join()
            background_process_MS.close()
            background_process_MS = None

        t12 = time.time()
        if config.verbosity >= 2:
            config.logger.info(f"Time for background inserts: {t12 - t11}")

        if config.verbosity >= 3:
            for name, size in sorted(((name, sys.getsizeof(value)) for name, value in locals().items()), key=lambda x: -x[1])[:10]:
                config.logger.info("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
            config.logger.info(f'{list_fds()}')

    # Error-Handling for a whole input line
    except:

        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        errortext = '\n'.join([str(e), str(f), str(g)]) + '\n\n'
        config.errorlog.add_error(errortext)
        if background_process_MS is not None:
            background_process_MS.join()
            background_process_MS.close()
            background_process_MS = None

        amount_of_structures = 0
        sys.exit(1)

    if config.ray_init_counter >= 10000:
        ray.shutdown()
        #time.sleep(10)
        config.ray_init_counter = 0
        pass
    return amount_of_structures		


def main(filename: str, config: Config):
    # main function takes another argument custom db, default is False
    mrna_fasta = config.mrna_fasta

    session = 0  # This can later be used, structman_main.py could give a specific session id and the pipeline can then expand that session

    if config.verbosity >= 2:
        config.logger.info('==== Start main pipeline ====')
    if config.verbosity >= 1:
        # Name the code this run uses: a checkout does not change what is installed, and a
        # stale install of another branch has already cost several whole-PDB runs.
        from structman._build_stamp import describe_running
        config.logger.info(describe_running())

    config.errorlog.start(filename, config.outfolder)

    t0 = time.time()
    # Labels the rows this run writes into tables shared with other runs, such as the
    # classification memory samples.
    import socket
    config.run_id = f'{socket.gethostname()}:{time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t0))}'

    # The ray cold start used to happen here, before the input was even parsed. It
    # moved below buildQueue, where the number of input proteins is known, so that a
    # small input can defer it (see SMALL_INPUT_RAY_DEFERRAL_THRESHOLD). Nothing
    # between here and there touches ray.

    # Start tracemalloc before buildQueue so every allocation — including the
    # initial protein objects and any that accumulate across chunks — is tracked.
    if config.verbosity >= 2:
        import tracemalloc
        tracemalloc.start()
        config.logger.info(f'tracemalloc started; baseline RSS: {_rss_mb():.0f} MiB')

    '''
    mem_tracked_processes = []
    for proc in psutil.process_iter(['pid', 'name', 'username']):
        if proc.info['username'] == 'agr18':
            p_id = proc.info['pid']
            mem_tracked_processes.append(psutil.Process())
    '''

    if mrna_fasta is not None:
        if not os.path.exists(mrna_fasta):
            raise NameError("mRNA path not found: %s" % mrna_fasta)

    if isinstance(filename, str):
        # annovar-pipeline in case of vcf-file
        if filename.rsplit(".", 1)[1] == "vcf":
            """
            anno_db = "%s_annovar" % db_name.rsplit("_", 1)[0]  # BUG: undefined variable
            if config.verbosity >= 1:
                print('Convert vcf file format using Annovar')
            if mrna_fasta is not None:
                '... and using mrna file: ', mrna_fasta
            nfname = annovar.annovar_pipeline(filename, config.tax_id, config.annovar_path, config.db_address, config.db_user_name, config.db_password, anno_db, mrna_fasta, ref_id=config.ref_genome_id)
            """
            pass
        else:
            nfname = filename

    # Single line input check
    else:
        nfname = 'Single line input'
        single_line_inputs = filename
        if config.verbosity >= 1:
            config.logger.info('=== Single line input mode ===')

    t01 = time.time()

    if config.verbosity >= 2:
        config.logger.info(f"Time for preparation before buildQueue: {t01 - t0}")

    try:
        os.stat("%s/tmp_structman_pipeline" % (config.outfolder))
    except:
        os.mkdir("%s/tmp_structman_pipeline" % (config.outfolder))
    os.chdir("%s/tmp_structman_pipeline" % (config.outfolder))
    cwd = "%s/tmp_structman_pipeline" % (config.outfolder)

    config.temp_folder = cwd
    # Separate handle for the same path, because modelling.model rebinds
    # config.temp_folder to its own per-model process folder and runs during the output
    # stage. structman_main deletes this one once ray is down; deleting whatever
    # temp_folder happens to point at by then would take out a model folder instead.
    config.pipeline_temp_folder = cwd

    chunksize = config.chunksize

    if config.verbosity >= 1:
        config.logger.info(f"Call buildQueue with {chunksize=} and file: {nfname}")

    if nfname != 'Single line input':
        if config.fasta_input:
            proteins_chunks, genes = parseFasta(config, nfname)
            temp_infiles = [None]
        else:
            proteins_chunks, temp_infiles, genes = buildQueue(config, nfname)
    else:
        proteins_chunks, temp_infiles, genes = buildQueue(config, single_line_inputs)
        nfname = '/%s.' % (' '.join(single_line_inputs))

    t02 = time.time()
    if config.verbosity >= 2:
        config.logger.info(f"Time for buildQueue: {t02 - t01}")
    if config.verbosity >= 1:
        config.logger.info(f"Number of chunks: {len(proteins_chunks)}")

    set_small_input_ray_deferral(config, proteins_chunks, temp_infiles)

    if parallel_db_access(config):
        # need structman package path for ray
        ray_init(config)

        if config.verbosity >= 2:
            config.logger.info('ray init successful')

    newsession = False
    if session == 0:
        starttime = SQLDateTime()
        session = database.insertSession(starttime, nfname, config)
        newsession = True
    # The id 'structman database repair --session_id' takes, which still works once the
    # session tables are truncated and the input file can no longer be looked up.
    if config.verbosity >= 1:
        config.logger.info(f'Session {session} for {nfname}')
    session_name = (nfname.rsplit("/", 1)[1]).rsplit(".", 1)[0]

    total_amount_of_analyzed_structures = 0
    config.ray_init_counter = 0
    save_ids = None

    for nr_temp_file, temp_infile in enumerate(temp_infiles):
        if temp_infile is not None:
            if config.verbosity >= 1:
                config.logger.info(f'Infile splitting due to low memory system, processing infile split nr.: {nr_temp_file + 1} out of {len(temp_infiles)}')
            proteins_chunks: list[tuple[dict[str, protein_package.Protein], dict[str, indel_package.Indel]]]
            proteins_chunks, _, genes = buildQueue(config, temp_infile, already_split=True)
            os.remove(temp_infile)

        #The gene objects need to exist as a whole outside of the chunked core loop. Isoforms that correspond to the same gene can be chunked into different bins.
        database.geneCheck(genes, config)

        chunk_nr = 1
        for protein_list, indels in proteins_chunks:

            tc0 = time.time()
            config.indels_given_by_input = (len(indels) > 0)
            # Snapshot open FDs *before* this chunk so kill_overhead_sockets only
            # closes sockets that appeared and leaked *within* this chunk, never
            # Ray's persistent RPC connections from earlier chunks.
            ray_up_at_snapshot = ray.is_initialized()
            if sys.platform.startswith('linux'):
                save_ids = kill_overhead_sockets(save_ids=None)

            if config.verbosity >= 1:
                config.logger.info(f"Chunk {chunk_nr}/{len(proteins_chunks)}")
            config.current_chunk = chunk_nr
            chunk_nr += 1

            tc1 = time.time()
            protein_list, indels, multi_mutation_objects, removed_proteins, _ = sequenceScan(config, protein_list, indels)
            for removed_prot_id in removed_proteins:
                try:
                    del protein_list[removed_prot_id]
                except KeyError:
                    pass
            tc2 = time.time()
            if config.verbosity >= 2:
                config.logger.info(f'Time for sequenceScan: {tc2-tc1}')

            if config.verbosity >= 2:
                rss_pre_core = _rss_mb()
                config.logger.info(f'RSS before core: {rss_pre_core:.0f} MiB')

            amount_of_structures = core(protein_list, genes, indels, multi_mutation_objects, config, session)

            if config.verbosity >= 2:
                rss_post_core = _rss_mb()
                config.logger.info(f'RSS after core (before cleanup): {rss_post_core:.0f} MiB (delta: {rss_post_core - rss_pre_core:+.0f} MiB)')

            total_amount_of_analyzed_structures += amount_of_structures

            for prot_id in list(protein_list.keys()):
                del protein_list[prot_id]

            for indel_id in list(indels.keys()):
                del indels[indel_id]
            del multi_mutation_objects
            del removed_proteins

            if config.verbosity >= 3:
                try:
                    for name, size in sorted(((name, sys.getsizeof(value)) for name, value in globals().items()), key=lambda x: -x[1])[:10]:
                        config.logger.info("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
                    config.logger.info(f'{list_fds()}')
                except:
                    pass

            # A small input defers its ray cold start into core(), where
            # resolve_deferred_ray_init runs it right after protCheck - which is *after*
            # the snapshot above. Ray's sockets are then absent from save_ids and would
            # be closed here as if they had leaked, leaving the cluster on dead file
            # descriptors: the RPC layer starts reporting EBADF and the freed fd numbers
            # get recycled under files the output stage has already opened.
            #
            # Skipping this one round is enough. The next chunk snapshots with ray
            # already up, so its sockets are protected and cleaning resumes.
            ray_started_in_chunk = (not ray_up_at_snapshot) and ray.is_initialized()

            if sys.platform.startswith('linux') and not ray_started_in_chunk:
                kill_overhead_sockets(save_ids=save_ids)
            elif ray_started_in_chunk and config.verbosity >= 2:
                config.logger.info(
                    'Ray was started inside this chunk, skipping the socket cleanup so '
                    'its connections are not closed as leaks'
                )

            if config.verbosity >= 3:
                config.logger.info(f'save_ids: {save_ids}')
                config.logger.info(f'{list_fds()}')


            tc3 = time.time()
            if config.verbosity >= 1:
                config.logger.info(f'\nChunk time: {tc3-tc0}\n')

            # Undo the freeze autoTemplateSelection set, so the collections in
            # clean_chunk_memory can actually reach this chunk's objects, and the next
            # chunk starts from an empty permanent generation.
            gc.unfreeze()

            clean_chunk_memory(config)

    os.chdir(config.outfolder)

    if total_amount_of_analyzed_structures == 0:
        if config.verbosity >= 1:
            config.logger.info("=== Nothing got processed, probably the input didn't contain any usable protein ID. ===")

    if config.verbosity >= 2:
        t03 = time.time()

    # The temp folder is deliberately NOT removed here. ray_init starts the raylet while
    # this process sits inside it, so every worker inherits it as its working directory,
    # and the cluster is handed to the output stage still running. Deleting it now leaves
    # those workers with a cwd that no longer exists: any worker that starts afterwards
    # dies on import, because numba pulls in coverage, whose import calls
    # os.path.realpath('.') and raises FileNotFoundError.
    #
    # structman_main.main removes it once the output stage is done and ray is shut down.

    if config.verbosity >= 2:
        t04 = time.time()
        config.logger.info(f'Time for folder cleanup: {t04 - t03}')

    # Drop every ObjectRef still reachable, so the pipeline's share of the object
    # store is released before the output stage starts filling it again. This also
    # keeps the reason the collect was put here in the first place: an ObjectRef
    # finalizer reached during ray's teardown calls back into a core worker that is
    # being destroyed, so the finalizers have to run while ray is still up.
    gc.collect()

    # Ray deliberately stays up here. The output stage brings it back with
    # ray_init(stop_existing=False) when it finds it down, and that cold start is
    # not cheap - it prestarts one python worker per cpu, which measured at ~33s of
    # a ~60s run (see the note in ray_utils.parallel_db_access). Handing the running
    # cluster over instead saves the whole second start. The teardown moved to
    # structman_main.main, which does it once per infile after the output stage.

    returncode = config.errorlog.stop()

    if newsession:
        endtime = SQLDateTime()
        database.updateSession(session, endtime, config)

    tend = time.time()
    if config.verbosity >= 1:
        config.logger.info(f'{total_amount_of_analyzed_structures} structures in total got analyzed')
        config.logger.info(f'Total runtime of StructMAn: {tend - t0}')

    '''
    total_memory_peak = 0
    for p in mem_tracked_processes:
        total_memory_peak += p.memory_info().rss

    total_memory_peak = total_memory_peak/1024./1024./1024.

    print('Accumulated memory peak',total_memory_peak,'Gb')
    '''

    return session, config, returncode
