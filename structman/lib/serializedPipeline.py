import errno
import gc
import os
import shutil
import sys
import time
import traceback
import statistics

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

from structman.lib.database import database
from structman.base_utils.base_utils import calculate_chunksizes, pack, unpack, is_alphafold_model, alphafold_model_id_to_file_path, aggregate_times, print_times, add_to_times
from structman.base_utils.ray_utils import ray_init, ray_hack
from structman.lib.lib_utils import process_mutations_str, generate_multi_mutation, extract_dedicated_tags


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
    if config.verbosity >= 2:
        print(warn_text)
    config.errorlog.add_warning(warn_text)

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

@ray.remote(max_calls=1)
def para_parse_pdb(package: list[str], pdb_path: str) -> list[tuple[str, str, dict]]:
    pdb_sequence_list: list[tuple[str, str, dict]] = pdbParser.getSequences(package, pdb_path)
    return pdb_sequence_list

def sequenceScan(
        config,
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
            print(f'Protein IDs going into sequence scan: {proteins}, genes set: {genes is not None}')
        else:
            print(f'{len(proteins)} number of Protein IDs going into sequence scan, genes set: {genes is not None}')
        t0 = time.time()

    for prot_id in proteins:
        uni_pos, tags = proteins[prot_id].popNone()

        if proteins[prot_id].is_sequence_set():
            #print(f'proteins[prot_id].is_sequence_set(): {proteins[prot_id].is_sequence_set()}')
            if uni_pos:
                fasta_inputs[prot_id] = tags
                #print(f'fasta_inputs[prot_id]: {fasta_inputs[prot_id]}')
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
            print("Amount of proteins going into NP sequenceScan: ", len(sequenceScanNP))
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
            print("Amount of proteins going into NM sequenceScan: ", len(sequenceScanNM))
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
                    print(f'Removed {nm_ref} from Input')
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
            print("Amount of proteins going into Ensembl sequenceScan: ", len(sequenceScanEnsembl))
        
        gene_sequence_map, gene_id_map = uniprot.get_ensembl_seqs(config, list(sequenceScanEnsembl.keys()))

        if config.verbosity >= 2:
            t02 = time.time()
            print(f"Size of gene_sequence_map: {len(gene_sequence_map)}, size of gene_id_map: {len(gene_id_map)}")
            print(f'Time for uniprot.get_ensembl_seqs: {t02- t01}')

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
                    print(f'Removed {transcript_id} from Input')
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
                    print(f'Removed {transcript_id} from Input')
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
            print(f'Retrieved {n} number of sequences in sequenceScanEnsembl with mean length: {mean_seq_len}')
            t03 = time.time()
            print(f'Time for sequenceScan ensembl Part 2: {t03- t02}')


    if len(sequenceScanProteins) > 0:
        if config.verbosity >= 1:
            print("Amount of proteins going into Uniprot sequenceScan: ", len(sequenceScanProteins))

        gene_sequence_map = uniprot.getSequencesPlain(sequenceScanProteins.keys(), config)

        #print_locals(locals().items())

        for u_ac in gene_sequence_map:
            if gene_sequence_map[u_ac][0] == 1 or gene_sequence_map[u_ac][0] == 0 or gene_sequence_map[u_ac][0] is None:
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

        #print_locals(locals().items())

        del gene_sequence_map
        #print_locals(locals().items())
    
    if len(sequenceScanPDB) > 0:
        if len(sequenceScanPDB) < config.proc_n:
            pdb_sequence_list = pdbParser.getSequences(sequenceScanPDB.keys(), pdb_path)
        else:
            small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(config.proc_n, len(sequenceScanPDB))
            package: list[str] = []
            parse_pdb_processes = []
            
            for pdb_id in sequenceScanPDB.keys():
                package.append(pdb_id)
                if len(parse_pdb_processes) < n_of_big_chunks:
                    if len(package) == big_chunksize:
                        parse_pdb_processes.append(para_parse_pdb.remote(package, pdb_path))
                        package = []
                else:
                    if len(package) == small_chunksize:
                        parse_pdb_processes.append(para_parse_pdb.remote(package, pdb_path))
                        package = []
            if len(package) > 0:
                parse_pdb_processes.append(para_parse_pdb.remote(package, pdb_path))
                package = []

            parse_pdb_results = ray.get(parse_pdb_processes)
            pdb_sequence_list = []
            for pdb_sequence_sublist in parse_pdb_results:
                pdb_sequence_list.extend(pdb_sequence_sublist)

        if config.verbosity >= 6:
            print(sequenceScanPDB)
            print(pdb_sequence_list)

        for pdb_chain_tuple, seq, res_pos_map in pdb_sequence_list:
            if seq == '':
                removed_proteins.append(pdb_chain_tuple)
                continue
            proteins[pdb_chain_tuple].sequence = seq

            tags, uni_pos = sequenceScanPDB[pdb_chain_tuple]
            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False

            residue_id_backmap = {}
            for res_id in res_pos_map:
                residue_id_backmap[res_pos_map[res_id]] = res_id

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if seq_pos not in residue_id_backmap:
                    continue
                res_id = residue_id_backmap[seq_pos]
                if res_id in proteins[pdb_chain_tuple].res_id_map:
                    proteins[pdb_chain_tuple].res_id_map[res_id].pos = seq_pos
                    proteins[pdb_chain_tuple].positions[seq_pos] = proteins[pdb_chain_tuple].res_id_map[res_id]
                    if tags is None:
                        continue
                    proteins[pdb_chain_tuple].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if not proteins[pdb_chain_tuple].contains_position(seq_pos):
                        position = position_package.Position(pos=seq_pos, pdb_res_nr=res_id, wt_aa=aa, tags=tags, session_less = session_less)
                        proteins[pdb_chain_tuple].add_positions([position])
                        proteins[pdb_chain_tuple].res_id_map[res_id] = position
                    elif not res_id in proteins[pdb_chain_tuple].res_id_map:
                        proteins[pdb_chain_tuple].res_id_map[res_id] = proteins[pdb_chain_tuple].positions[seq_pos]

        for pdb_chain_tuple in sequenceScanPDB:
            if not pdb_chain_tuple in proteins:
                removed_proteins.append(pdb_chain_tuple) 
            

    if len(fasta_inputs) > 0:
        for prot_id in fasta_inputs:
            seq = proteins[prot_id].sequence

            tags = fasta_inputs[prot_id]

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if not proteins[prot_id].contains_position(seq_pos):
                    position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True)
                    proteins[prot_id].add_positions([position])

    if config.verbosity >= 2:
        t1 = time.time()
        print(f'sequenceScan Part 1: {t1-t0}')

    if config.verbosity >= 3:
        print('Before indel mutation with:', len(indels), 'number of indels')

    for primary_protein_id in indels:
        for indel_notation in indels[primary_protein_id]:
            indels[primary_protein_id][indel_notation].mutate_sequence(proteins)

    if config.verbosity >= 2:
        t2 = time.time()
        print(f'sequenceScan Part 2: {t2-t1}')

    multi_mutation_objects = []

    if (not config.only_wt) and (not skip_mm_calc):
        if config.verbosity >= 3:
            print('Before multi mutations mutation')
        for primary_protein_id in list(proteins.keys()).copy():  # iterate over snapshot, since the mutant proteins get added to proteins in the loop
            if proteins[primary_protein_id].multi_mutations is not None:
                # creates mutation protein object
                multi_mutation_objects += proteins[primary_protein_id].create_multi_mutations(proteins, config)

    if config.verbosity >= 2:
        t3 = time.time()
        print(f'sequenceScan Part 3: {t3-t2}')

    return proteins, indels, multi_mutation_objects, removed_proteins, genes


def parseFasta(config, nfname):
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
        line = line[:-1]
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
            prot_tags_map[entry_id] = protein_specific_tags

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
        print(f'Parse loop times: {t1_sum} {t2_sum} {t3_sum}')

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
            print(f'In parseFasta: Protein tags for {prot_id}: {protein_specific_tags}')
            t_0 = time.time()

        uniprot.integrate_protein(config, proteins, genes, indel_map, prot_id, prot_id, prot_map, gene_id = gene_id, protein_specific_tags = protein_specific_tags)
        proteins[prot_id].sequence = seq

        if config.verbosity >= 3:
            t_1 = time.time()
            print(f'Time for integrate_protein: {t_1-t_0}')

    genes = gene_isoform_check(proteins, genes, indel_map, config)

    outlist: list[tuple[dict[str, protein_package.Protein], dict[str, indel_package.Indel]]] = input_chunking(config, proteins, indel_map)

    return outlist, genes


def add_to_prot_map(
        prot_map: dict[str, tuple[list[position_package.Position], list[tuple], list[tuple[list[tuple | indel_package.Indel]]]]],
        sp_id: str,
        positions: list[position_package.Position],
        multi_mutations: list[tuple | indel_package.Indel],
        config):
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
        config,
        proteins: dict[str, protein_package.Protein],
        indels:dict[str, indel_package.Indel]
        ):
    
    outlist: list[tuple[dict[str, protein_package.Protein], dict[str, indel_package.Indel]]] = []

    s = len(proteins)
    if config.verbosity >= 1:
        print("Total proteins: ", s)
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


def buildQueue(config, filename, already_split=False):
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
                    print("Skipped input line:\n%s\nToo short.\n" % line)
                continue
            line = line.replace(' ', '\t')
            words = line.split("\t")
            if len(words) < 1:
                if config.verbosity >= 1:
                    print("Skipped input line:\n%s\nToo few words.\n" % line)
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
                print("Skipped input line:\n%s\nToo short.\n" % line)
            continue
        line = line.replace(' ', '\t')
        words = line.split("\t")
        if len(words) < 1:
            if config.verbosity >= 1:
                print("Skipped input line:\n%s\nToo few words.\n" % line)
            continue
        sp_id = words[0]  # .replace("'","\\'")

        if len(sp_id) < 2:
            if config.verbosity >= 1:
                print("Skipped input line:\n%s\nID too short.\n" % line)
            continue
        elif (len(sp_id) == 6 and sp_id.count(':') == 1):
            sp_id = f'{sp_id[:4].upper()}:{sp_id[-1]}'
        elif (len(sp_id) == 4 and sp_id[0].isnumeric()):
            sp_id = sp_id.upper()

        tags = set()

        if config.verbosity >= 6:
            print(f'Processing input line: {line}')

        if len(words) > 2:
            tags = set(words[2].split(','))

        try:
            if len(words) == 1 or words[1] == '':
                position = position_package.Position(tags=tags)
                indel = None
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
        print("buildQueue Part 1: ", str(t1 - t0), len(ac_map), len(prot_tags_map))

    proteins: dict[str, protein_package.Protein]
    indel_map: dict[str, indel_package.Indel]
    genes: dict[str, gene_package.Gene]
    proteins, indel_map, genes = uniprot.IdMapping(config, ac_map, id_map, np_map, pdb_map, hgnc_map, nm_map, ensembl_map, prot_gene_map, prot_tags_map, pdbs_to_expand)

    t2 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 2: ", str(t2 - t1))

    if config.verbosity >= 3:
        print(f'Size of proteins after IdMapping: {len(proteins)}')
    if config.verbosity >= 6:
        print(f'Proteins State after IdMapping:')
        for prot_id in proteins:
            proteins[prot_id].print_state()

    genes = gene_isoform_check(proteins, genes, indel_map, config)

    t3 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 3: ", str(t3 - t2))

    outlist: list[tuple[dict[str, protein_package.Protein], dict[str, indel_package.Indel]]] = input_chunking(config, proteins, indel_map)

    t4 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 4: ", str(t4 - t3))

    return outlist, [None], genes


def gene_isoform_check(
        proteins: dict[str, protein_package.Protein],
        genes: dict[str, gene_package.Gene],
        indel_map,
        config) -> dict[str, gene_package.Gene]:
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
        print(f'gene_isoform_check - Part 1: {t1-t0}')

    if config.verbosity >= 3:
        print(f'In gene_isoform_check: Number of extra_sequence_scan_proteins: {len(extra_sequence_scan_proteins)}, Number of Genes: {len(genes)}')

    if config.verbosity >= 4:
        if len(genes) < 10:
            for gene_id in genes:
                genes[gene_id].print_content()

    extra_sequence_scan_proteins, indel_map, _, removed_proteins, genes = sequenceScan(config, extra_sequence_scan_proteins, indel_map, genes = genes, skip_mm_calc = True)

    if config.verbosity >= 2:
        t2 = time.time()
        print(f'gene_isoform_check - Part 2: {t2-t1}')

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
        print(f'Extra sequenceScan replace {n} number of protein objects with mean sequence lenght: {mean_seq_len}')

    if config.verbosity >= 2:
        t3 = time.time()
        print(f'gene_isoform_check - Part 3: {t3-t2}')

    for removed_prot_id in removed_proteins:
        del proteins[removed_prot_id]

    if config.verbosity >= 2:
        t4 = time.time()
        print(f'gene_isoform_check - Part 4: {t4-t3}')


    config_dump = ray.put(config)
    pieces = []
    total_cost: int = 0

    for gene_id in genes:
        if len(genes[gene_id].proteins) < 2:
            continue
        isoform_list = list(genes[gene_id].proteins.values())

        sequence_map = {}
        iso_pairs = []

        cost: int = 0

        for iso1 in isoform_list:
            if iso1 not in proteins:
                if iso1 in genes[gene_id].proteins:
                    del genes[gene_id].proteins[iso1]
                continue
            for iso2 in isoform_list:
                if iso2 not in proteins:
                    if iso2 in genes[gene_id].proteins:
                        del genes[gene_id].proteins[iso2]
                    continue
                if iso1 == iso2:
                    continue

                if config.verbosity >= 6:
                    print(f'Generating Isoform-Multimutation for {iso1} and {iso2}')

                if iso1 not in sequence_map:
                    sequence_map[iso1] = proteins[iso1].sequence
                if iso2 not in sequence_map:
                    sequence_map[iso2] = proteins[iso2].sequence
                iso_pairs.append((iso1, iso2))
                cost += len(sequence_map[iso1]) * len(sequence_map[iso2])

        if len(iso_pairs) == 0:
            continue

        total_cost += cost

        pieces.append((sequence_map, iso_pairs, cost))

    package_cost: int = total_cost // config.proc_n

    if config.verbosity >= 2:
        t5 = time.time()
        print(f'gene_isoform_check - Part 5: {t5-t4} {total_cost=} {package_cost=}')

    current_package_cost = 0
    package = []

    mm_generation_remote_ids = []
    for (sequence_map, iso_pairs, cost) in pieces:
        if cost > package_cost*2:
            splits = cost // package_cost
            split_size = (len(iso_pairs) // splits) + 1
            for i in range(splits):
                iso_pair_split = iso_pairs[i*split_size : i*(split_size+1)]
                if config.verbosity >= 3:
                    print(f'Starting new para_mm_generatio process {cost=}')
                mm_generation_remote_ids.append(para_mm_generation.remote(config_dump, pack([(sequence_map, iso_pair_split)])))
            continue

        package.append((sequence_map, iso_pairs))
        current_package_cost += cost
        if current_package_cost >= package_cost:
            if config.verbosity >= 3:
                print(f'Starting new para_mm_generatio process {current_package_cost=}')
            mm_generation_remote_ids.append(para_mm_generation.remote(config_dump, pack(package)))
            package = []
            current_package_cost = 0

    if len(package) > 0:
        mm_generation_remote_ids.append(para_mm_generation.remote(config_dump, pack(package)))

    if config.verbosity >= 2:
        t6 = time.time()
        print(f'gene_isoform_check - Part 6: {t6-t5} {len(mm_generation_remote_ids)=}')

    mm_generation_remote_results = ray.get(mm_generation_remote_ids)

    if config.verbosity >= 2:
        t7 = time.time()
        print(f'gene_isoform_check - Part 7: {t7-t6}')

    isoform_mm_matrix = {}

    for mm_generation_results in mm_generation_remote_results:
        for iso1, iso2, mm_results in unpack(mm_generation_results):
            if isinstance(mm_results, str):
                config.errorlog.add_warning(f'Multi mutation generation failed for {iso1} and {iso2} with error:\n{mm_results}')
                continue

            multi_mutation = mm_results

            for indel_or_snv in multi_mutation:
                if isinstance(indel_or_snv, tuple):
                    if config.verbosity >= 6:
                        print(f'Add Position:')
                        indel_or_snv[0].print_state()
                        print(f'... to {iso1}')
                    proteins[iso1].add_positions([indel_or_snv[0]])
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
        print(f'gene_isoform_check - Part 8: {t8-t7}')

    return genes

@ray.remote(max_calls = 1)
def para_mm_generation(config, packages):
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


def getSequences(proteins, config):

    t0 = time.time()

    uniprot.getSequences(proteins, config)

    t1 = time.time()
    if config.verbosity >= 2:
        print("Time for getSequences Part 1: %s" % str(t1 - t0))

    if not config.read_only_mode:
        database.addProtInfos(proteins, config)

    t2 = time.time()
    if config.verbosity >= 2:
        print("Time for getSequences Part 2: %s" % str(t2 - t1))


def autoTemplateSelection(config, proteins):

    if config.verbosity >= 1:
        print('Sequence search with mmseqs2')

    filtering_dump = None

    t0 = time.time()

    # MMseqs search function takes one more argument now: custom_db
    search_results: list[tuple[dict[str, dict[str, dict[str, tuple[float, float, set[str], int, int]]]], set[str], bool]]
    search_results = MMseqs2.search(proteins, config)

    t1 = time.time()

    if search_results is not None:
        for search_result in search_results:
            if len(search_result) == 0:
                continue
            raw_structure_map: dict[str, dict[str, dict[str, tuple[float, float, set[str], int, int]]]]
            pdb_ids = set[str]
            raw_structure_map, pdb_ids, is_model_db = search_result

            info_map: dict[str, tuple[float, dict[str, set[str]]]] = {}
            filtering_results = []

            # if it is custom db, it creates fake oligo map
            if config.custom_db_path:
                for custom_id in pdb_ids:
                    fake_oligo_map = {'A' : set(['A'])}
                    info_map[custom_id] = (2.0, fake_oligo_map)

            elif not is_model_db:
                
                # calls filtering dump with custom db path rather than pdb path
                if config.proc_n > 1:
                    if config.custom_db_path:
                        filtering_dump = ray.put(config.custom_db_path)
                    else:
                        filtering_dump = ray.put(config.pdb_path)

                small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(config.proc_n, len(pdb_ids))

                n = 0
                chunk = []
                for pdb_id in pdb_ids:

                    if config.structure_limiter is not None:
                        if pdb_id[:4] != config.structure_limiter:
                            #print(f'Structure {pdb_id} got denied by limiter ({config.structure_limiter})')
                            continue

                    chunk.append(pdb_id)
                    if n < n_of_big_chunks:
                        if len(chunk) == big_chunksize:
                            if config.proc_n > 1:
                                filtering_results.append(filter_structures_wrapper.remote(chunk, filtering_dump))
                            else:
                                filtering_results.append(filter_structures(chunk, config.pdb_path))
                            chunk = []
                            n += 1
                    elif n <= n_of_big_chunks + n_of_small_chunks:
                        if len(chunk) == small_chunksize:
                            if config.proc_n > 1:
                                filtering_results.append(filter_structures_wrapper.remote(chunk, filtering_dump))
                            else:
                                filtering_results.append(filter_structures(chunk, config.pdb_path))
                            chunk = []
                            n += 1
                    else:
                        print('Warning: too few chunks')

                if len(chunk) != 0:
                    if config.proc_n > 1:
                        filtering_results.append(filter_structures_wrapper.remote(chunk, filtering_dump))
                    else:
                        filtering_results.append(filter_structures(chunk, config.pdb_path))

                if config.proc_n > 1:
                    filtering_out = ray.get(filtering_results)
                else:
                    filtering_out = filtering_results

                for out_chunk in filtering_out:

                    if isinstance(out_chunk, str):
                        config.errorlog.add_error(out_chunk)
                    else:
                        for (pdb_id, resolution, homomer_dict) in out_chunk:
                            info_map[pdb_id] = (resolution, homomer_dict)

            else:
                for model_id in pdb_ids:
                    fake_oligo_map = {'A' : set(['A'])}
                    info_map[model_id] = (2.0, fake_oligo_map)

            for prot_id in raw_structure_map:
                for pdb_id in raw_structure_map[prot_id]:
                    if config.fast_pdb_annotation:
                        #This filters every hit that is not identical to the input (pdb maps to pdb)
                        if prot_id[:4] != pdb_id[0:4]:
                            continue

                    for chain in raw_structure_map[prot_id][pdb_id]:
                        
                        if pdb_id not in info_map:
                            continue
                        resolution, homomer_dict = info_map[pdb_id]
                        if resolution is None:
                            continue
                        if resolution > config.option_res_thresh:
                            continue
                        
                        try:
                            oligo = homomer_dict[chain]
                        except:
                            oligo = raw_structure_map[prot_id][pdb_id][chain][2]

                        struct_anno = structure_package.StructureAnnotation(prot_id, pdb_id, chain)
                        proteins.add_annotation(prot_id, pdb_id, chain, struct_anno)

                        if config.verbosity >= 6:
                            print(f'Adding structural annotation: {prot_id} -> {pdb_id}:{chain}')

                        if not pdb_id in proteins.structures:
                            struct = structure_package.Structure(pdb_id, chain, oligo=oligo, mapped_proteins=[prot_id], seq_len = raw_structure_map[prot_id][pdb_id][chain][4])
                            proteins.add_structure(pdb_id, chain, struct)
                            
                        elif chain not in proteins.structures[pdb_id]:
                            struct = structure_package.Structure(pdb_id, chain, oligo=oligo, mapped_proteins=[prot_id], seq_len = raw_structure_map[prot_id][pdb_id][chain][4])
                            proteins.add_structure(pdb_id, chain, struct) 
                        else:
                            proteins.add_mapping_to_structure(pdb_id, chain, prot_id)

                        if pdb_id not in proteins.complexes:
                            compl = complex_package.Complex(pdb_id, resolution, homomers=homomer_dict)
                            proteins.add_complex(pdb_id, compl)

    t2 = time.time()
    if config.verbosity >= 2:
        print("Template Selection Part 1: %s" % (str(t1 - t0)))
        print("Template Selection Part 2: %s" % (str(t2 - t1)))

    del filtering_dump

@ray.remote(max_calls = 1)
def filter_structures_wrapper(chunk, filtering_dump):

    pdb_path = filtering_dump
    return filter_structures(chunk, pdb_path)

def filter_structures(chunk, pdb_path):
    try:
        outs = []
        for pdb_id in chunk:
            resolution, homomer_dict = pdbParser.getInfo(pdb_id, pdb_path)
            outs.append((pdb_id, resolution, homomer_dict))
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        outs = 'Filter structure error:\n%s\n%s\n%s' % (e, str(f), g)

    return outs

def package_alignment_processes(prots_todo, cost_map, optimal_chunk_cost, proteins, protein_packages, alignment_results, started_processes, chunk, chunk_cost, config, mapping_dump):
    done = []
    package_cost = None

    if config.verbosity >= 5:
        print(f'Call of package_alignment_processes: Number of proteins: {len(prots_todo)}, optimal chunk cost: {optimal_chunk_cost}, chunk cost: {chunk_cost}')

    for pos, prot_id in enumerate(prots_todo):
        package_cost = cost_map[prot_id]
        if package_cost > optimal_chunk_cost:  #packages that are greater than the optimal_chunk_cost get split
            structure_infos_split_a = []
            structure_infos_split_b = []
            seq = proteins.get_sequence(prot_id)
            (prot_specific_mapping_dump, structure_infos) = protein_packages[prot_id]
            split_cost = 0
            for structure_info in structure_infos:
                pdb_id, chain, oligo = structure_info
                str_len = proteins.structures[pdb_id][chain].get_seq_len()
                if str_len is not None:
                    cost = len(seq)*str_len
                else:
                    cost = (len(seq))**2
                if split_cost < optimal_chunk_cost:
                    structure_infos_split_a.append(structure_info)
                    split_cost += cost
                else:
                    structure_infos_split_b.append(structure_info)

            if config.verbosity >= 6:
                print('Splitting protein alignment:', prot_id, 'split_cost:', split_cost)
            if config.proc_n > 1:
                alignment_results.append(align_remote_wrapper.remote(mapping_dump, [(prot_specific_mapping_dump, structure_infos_split_a)]))
                started_processes += 1
            else:
                alignment_results.append(align(config, chunk))
            if config.verbosity >= 6:
                print('A, Start alignment package with cost:', split_cost)

            
            if len(structure_infos_split_b) > 0:
                protein_packages[prot_id] = (prot_specific_mapping_dump, structure_infos_split_b)
                cost_map[prot_id] = package_cost - split_cost
            else:
                done.append(pos)
        elif chunk_cost + package_cost <= optimal_chunk_cost:
            chunk.append(protein_packages[prot_id])
            done.append(pos)
            chunk_cost += package_cost
            if chunk_cost >= optimal_chunk_cost:
                if config.proc_n > 1:
                    alignment_results.append(align_remote_wrapper.remote(mapping_dump, chunk))
                    started_processes += 1
                else:
                    alignment_results.append(align(config, chunk))
                if config.verbosity >= 6:
                    print('B, Start alignment package with cost:', chunk_cost)
                chunk = []
                chunk_cost = 0
        
        if started_processes >= config.proc_n:
            break

    done.reverse()
    for pos in done:
        del prots_todo[pos]

    if (len(chunk) > 0) and len(prots_todo) == 0:
        if config.proc_n > 1:
            alignment_results.append(align_remote_wrapper.remote(mapping_dump, chunk))
            started_processes += 1
        else:
            alignment_results.append(align(config, chunk))
        if config.verbosity >= 6:
            print('C, Start alignment package with cost:', chunk_cost)
        chunk = []
        chunk_cost = 0

    if (len(chunk) > 0) and package_cost is not None:
        if (chunk_cost + package_cost > optimal_chunk_cost)  and (started_processes < config.proc_n):
            if config.proc_n > 1:
                alignment_results.append(align_remote_wrapper.remote(mapping_dump, chunk))
                started_processes += 1
            else:
                alignment_results.append(align(config, chunk))
            if config.verbosity >= 6:
                print('D, Start alignment package with cost:', chunk_cost)
            chunk = []
            chunk_cost = 0

    return alignment_results, started_processes, protein_packages, cost_map, chunk, chunk_cost, prots_todo

def process_align_outs(proteins, config, align_outs, warn_map, alignment_insertion_list, structure_insertion_list, indel_analysis_follow_up,
                        safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes, skip_inserts):


    total_integration_times = []
    total_package_runtimes = []

    for out_chunks in align_outs:
        if config.proc_n > 1:
            out_chunks, package_runtimes = unpack(out_chunks)
        else:
            out_chunks, package_runtimes = out_chunks

        total_package_runtimes = aggregate_times(total_package_runtimes, package_runtimes)

        for out in out_chunks:
            times = []
            t_init = time.time()
            if len(out) == 3:
                config.errorlog.add_error('Illegal alignment output: %s' % (str(out)))
                continue

            if len(out) == 4:
                (prot_id, pdb_id, chain, warn_text) = out
                proteins.remove_annotation(prot_id, pdb_id, chain)
                if not skip_inserts:
                    sus_complexes.add(pdb_id)
                    if pdb_id not in sus_structures:
                        sus_structures[pdb_id] = set([chain])
                    else:
                        sus_structures[pdb_id].add(chain)
                if prot_id not in warn_map:
                    config.errorlog.add_warning(warn_text)
                elif config.verbosity >= 4:
                    print('Alignment failed:', prot_id, pdb_id, chain, warn_text)
                warn_map.add(prot_id)
                continue

            if len(out) == 5:
                (prot_id, pdb_id, chain, sub_infos, align_info) = out
                proteins.remove_annotation(prot_id, pdb_id, chain)
                if not skip_inserts:
                    sus_complexes.add(pdb_id)
                    if pdb_id not in sus_structures:
                        sus_structures[pdb_id] = set([chain])
                    else:
                        sus_structures[pdb_id].add(chain)
                if config.verbosity >= 6:
                    seq_id, alignment_text = align_info
                    print('Alignment got filtered:', prot_id, pdb_id, chain, len(sub_infos), seq_id, alignment_text)
                continue

            (prot_id, pdb_id, chain, alignment, seq_id, coverage, interaction_partners, chain_type_map,
             oligo, sub_infos, backmap, atom_count, last_residue, first_residue, chainlist, rare_residues) = out

            if config.verbosity >= 6:
                print(f'Receiving alignment results from {prot_id} {pdb_id} {chain}: {seq_id} {coverage}')

            t_now = add_to_times(times, t_init)

            proteins.set_coverage(prot_id, pdb_id, chain, coverage)
            proteins.set_sequence_id(prot_id, pdb_id, chain, seq_id)
            proteins.set_sub_infos(prot_id, pdb_id, chain, sub_infos)
            proteins.set_atom_count(pdb_id, atom_count)

            proteins.set_backmap(prot_id, pdb_id, chain, backmap)

            if indel_analysis_follow_up:
                proteins.set_alignment(prot_id, pdb_id, chain, alignment)

            t_now = add_to_times(times, t_now)

            proteins.set_interaction_partners(pdb_id, interaction_partners)
            proteins.set_chain_type_map(pdb_id, chain_type_map, chainlist)

            t_now = add_to_times(times, t_now)

            proteins.structures[pdb_id][chain].oligo = oligo

            proteins.set_last_residue(pdb_id, chain, last_residue)
            proteins.set_first_residue(pdb_id, chain, first_residue)

            t_now = add_to_times(times, t_now)

            if pdb_id not in structure_insertion_list:
                structure_insertion_list[pdb_id] = set([chain])
            else:
                structure_insertion_list[pdb_id].add(chain)

            if config.verbosity >= 5:
                print(f'In process_align_outs: {prot_id=} {pdb_id=} {chain_type_map=}')

            for chain_id in chain_type_map:
                chain_type: str = chain_type_map[chain_id]
                if chain_type == 'Protein':
                    structure_insertion_list[pdb_id].add(chain_id)
                    if pdb_id not in safe_structures:
                        safe_structures[pdb_id] = set([chain_id])
                    else:
                        safe_structures[pdb_id].add(chain_id)
                    if chain_id not in proteins.structures[pdb_id]:
                        homomers = proteins.complexes[pdb_id].get_homomers(chain_id)
                        struct = structure_package.Structure(pdb_id, chain_id, oligo=homomers, new_interacting_chain = True)
                        proteins.add_structure(pdb_id, chain_id, struct)

            t_now = add_to_times(times, t_now)

            safe_complexes.add(pdb_id)
            if pdb_id not in safe_structures:
                safe_structures[pdb_id] = set([chain])
            else:
                safe_structures[pdb_id].add(chain)

            t_now = add_to_times(times, t_now)

            prot_db_id = proteins.get_protein_db_id(prot_id)
            alignment_insertion_list.append((prot_id, prot_db_id, pdb_id, chain, alignment))

            t_now = add_to_times(times, t_now)

            config.rare_residues.update(rare_residues)

            t_now = add_to_times(times, t_now)

            total_integration_times = aggregate_times(total_integration_times, times)

        started_processes -= 1

    return (warn_map, alignment_insertion_list, structure_insertion_list, total_package_runtimes, total_integration_times,
            safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes)


def paraAlignment(config, proteins, skip_inserts=False, indel_analysis_follow_up=False, get_all_alignments=False):
    indel_analysis_follow_up = indel_analysis_follow_up or get_all_alignments

    t0 = time.time()

    database.getAlignments(proteins, config, get_all_alignments=get_all_alignments)

    t1 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 1: %s" % (str(t1 - t0)))
    prot_ids = list(proteins.get_protein_ids())

    if config.proc_n > 1:
        mapping_dump = ray.put(config)
    else:
        mapping_dump = None

    mapping_results = []

    sus_complexes = set()
    sus_structures: dict[str, set[str]] = {}

    safe_structures: dict[str, set[str]] = {}
    safe_complexes = set()

    task_list: list[tuple[str, dict, int, list]] = []

    total_times = []

    for prot_id in prot_ids:
        if not proteins.is_protein_stored(prot_id):
            continue
        if prot_id not in proteins.indels:
            if proteins.is_completely_stored(prot_id) and not get_all_alignments:
                continue

        times = []
        ti = time.time()
        annotation_list = proteins.get_protein_annotation_list(prot_id)
        ti = add_to_times(times, ti)

        aaclist = proteins.getAACList(prot_id)
        prot_db_id = proteins.get_protein_db_id(prot_id)

        total_subtimes = []
        struct_tasks: list[tuple[str, str, str, str, str]] = []
        for (pdb_id, chain) in annotation_list:
            if not proteins.is_annotation_stored(pdb_id, chain, prot_id):
                continue

            sub_times = []
            ts = time.time()

            structure_id = proteins.get_structure_db_id(pdb_id, chain)
            ts = add_to_times(sub_times, ts)
            if indel_analysis_follow_up:
                get_alignment_out = proteins.get_alignment(prot_id, pdb_id, chain)
                if isinstance(get_alignment_out, str):
                    config.errorlog.add_error('Get alignment error: %s. %s %s %s %s' % (get_alignment_out, prot_id, pdb_id, chain, str(get_all_alignments)))
                    continue
                else:
                    target_seq, template_seq = get_alignment_out
            else:
                target_seq, template_seq = proteins.pop_alignment(prot_id, pdb_id, chain)
            ts = add_to_times(sub_times, ts)
            struct_tasks.append((pdb_id, chain, structure_id, target_seq, template_seq))
            
            ts = add_to_times(sub_times, ts)
            if pdb_id not in safe_structures:
                safe_structures[pdb_id] = set([chain])
            else:
                safe_structures[pdb_id].add(chain)
            ts = add_to_times(sub_times, ts)
            safe_complexes.add(pdb_id)
            ts = add_to_times(sub_times, ts)

            aggregate_times(total_subtimes, sub_times)
        
        task_list.append((prot_id, aaclist, prot_db_id, struct_tasks))

        times.append(total_subtimes)
        ti = add_to_times(times, ti)

        aggregate_times(total_times, times)

    t11 = time.time()
    if config.verbosity >= 2:
        print_times(total_times, label='process getAlignment prot_ids')
        print("Alignment Part 1.1: %s" % (str(t11 - t1)))

    if len(task_list) <= config.proc_n:
        for task in task_list:
            if config.proc_n > 1:
                mapping_results.append(paraMap_wrapper.remote(mapping_dump, pack([task])))
            else:
                mapping_results.append(paraMap(config, [task]))
    else:
        small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(config.proc_n, len(task_list))
        chunk = []
        n = 0
        for task in task_list:
            chunk.append(task)
            if n < n_of_big_chunks:
                if len(chunk) == big_chunksize:
                    if config.proc_n > 1:
                         mapping_results.append(paraMap_wrapper.remote(mapping_dump, pack(chunk)))
                    else:
                         mapping_results.append(paraMap(config, chunk))
                    chunk = []
                    n += 1
            elif n <= n_of_big_chunks + n_of_small_chunks:
                if len(chunk) == small_chunksize:
                    if config.proc_n > 1:
                        mapping_results.append(paraMap_wrapper.remote(mapping_dump, pack(chunk)))
                    else:
                        mapping_results.append(paraMap(config, chunk))
                    chunk = []
                    n += 1
            else:
                print('Warning: too few chunks')

        if len(chunk) != 0:
            if config.proc_n > 1:
                mapping_results.append(paraMap_wrapper.remote(mapping_dump, pack(chunk)))
            else:
                mapping_results.append(paraMap(config, chunk))

    t12 = time.time()
    if config.verbosity >= 2:
        print(f"Alignment Part 1.2: {t12 - t11}")

    if len(task_list) > 100:
        gc.collect()

    t2 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 2: %s" % (str(t2 - t12)))

    if config.proc_n > 1:
        mappings_outs = ray.get(mapping_results)
    else:
        mappings_outs = mapping_results

    t21 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 2.1: %s" % (str(t21 - t2)))

    for results in mappings_outs:
        for (prot_id, structure_id, chain, sub_infos, atom_count, last_residue, first_residue, backmap) in results:
            proteins.set_sub_infos(prot_id, structure_id, chain, sub_infos)
            proteins.set_atom_count(structure_id, atom_count)
            proteins.set_last_residue(structure_id, chain, last_residue)
            proteins.set_first_residue(structure_id, chain, first_residue)
            proteins.set_backmap(prot_id, structure_id, chain, backmap)

    t3 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 3: %s" % (str(t3 - t21)))

    total_cost = 0
    max_optimal_cost = 100000000
    cost_map = {}

    protein_packages = {}
    protein_bins: list[tuple[list[str], set[str]]] = []

    for prot_id in prot_ids:
        
        annotation_list: list[tuple[str, str]] = proteins.get_protein_annotation_list(prot_id)

        seq = proteins.get_sequence(prot_id)
        aaclist = proteins.getAACList(prot_id)
        structure_infos = []
        structure_set: set[str] = set()
        for (pdb_id, chain) in annotation_list:
            if proteins.is_annotation_stored(pdb_id, chain, prot_id):
                continue
            if prot_id not in cost_map:
                cost_map[prot_id] = 0
            if config.proc_n > 1:
                str_len = proteins.structures[pdb_id][chain].get_seq_len()
                if str_len is not None:
                    cost_map[prot_id] += len(seq)*str_len
                else:
                    cost_map[prot_id] += (len(seq))**2
            oligo = proteins.structures[pdb_id][chain].oligo
            structure_infos.append((pdb_id, chain, oligo))
            structure_set.add(pdb_id)

        if prot_id not in cost_map:
            continue

        binned = False
        for prot_bin, bin_structures in protein_bins:
            if pdb_id in bin_structures:
                prot_bin.append(prot_id)
                bin_structures.update(structure_set)
                binned = True
                break
        if not binned:
            protein_bins.append(([prot_id], structure_set))

        if config.proc_n > 1:
            prot_specific_mapping_dump = ray.put((prot_id, seq, aaclist))
            del seq
            del aaclist
        else:
            prot_specific_mapping_dump = (prot_id, seq, aaclist)
        protein_packages[prot_id] = (prot_specific_mapping_dump, structure_infos)

        total_cost += cost_map[prot_id]

    optimal_chunk_cost = min([(total_cost// (4 * config.proc_n)) + 1, max_optimal_cost])

    if config.verbosity >= 3:
        print('Total cost for alignment:', total_cost, optimal_chunk_cost)

    alignment_results = []

    chunk = []
    chunk_cost = 0

    prots_todo: list[str] = []
    for prot_bin, _ in protein_bins:
        prots_todo.extend(prot_bin)

    started_processes = 0

    while started_processes < config.proc_n and len(prots_todo) > 0:
        alignment_results, started_processes, protein_packages, cost_map, chunk, chunk_cost, prots_todo = package_alignment_processes(prots_todo, cost_map, optimal_chunk_cost, proteins, protein_packages, alignment_results, started_processes, chunk, chunk_cost, config, mapping_dump)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Alignment Part 4: %s" % (str(t4 - t3)))


    alignment_insertion_list = []
    structure_insertion_list: dict[str, set[str]] = {}
    warn_map = set()

    total_package_runtimes = []
    total_integration_times = []

    if config.proc_n > 1:

        while True:

            if started_processes >= config.proc_n or len(prots_todo) == 0:
                if config.verbosity >= 5:
                    print('Ray wait without timeout')

                ready, not_ready = ray.wait(alignment_results)
            else:
                if config.verbosity >= 5:
                    print('Ray wait with timeout')
                ready, not_ready = ray.wait(alignment_results, timeout = 0.01)

            if len(ready) > 0:

                if config.verbosity >= 2:
                    t_get_0 = time.time()
                try:
                    align_outs = ray.get(ready)
                except:
                    config.errorlog.add_error('Critical error in ray.get')
                    return
                if config.verbosity >= 2:
                    t_get_1 = time.time()

                if config.verbosity >= 6:
                    print('Times for get:', (t_get_1 - t_get_0))

                (warn_map, alignment_insertion_list, structure_insertion_list, package_runtimes, integration_times,
                    safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes) = process_align_outs(proteins, 
                            config, align_outs, warn_map, alignment_insertion_list, structure_insertion_list,
                            indel_analysis_follow_up, safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes, skip_inserts)

                total_package_runtimes = aggregate_times(total_package_runtimes, package_runtimes)
                total_integration_times = aggregate_times(total_integration_times, integration_times)
            t_new_0 = time.time()

            new_alignment_results = []
            if len(prots_todo) > 0 or len(chunk) > 0:
                new_alignment_results, started_processes, protein_packages, cost_map, chunk, chunk_cost, prots_todo = package_alignment_processes(prots_todo, cost_map, optimal_chunk_cost, proteins, protein_packages, new_alignment_results, started_processes, chunk, chunk_cost, config, mapping_dump)
                

            alignment_results = new_alignment_results + not_ready

            t_new_1 = time.time()
            if config.verbosity >= 5:
                print('Times for starting new processes:', (t_new_1 - t_new_0))
                print('Proteins to do:', len(prots_todo), 'Chunk:', len(chunk), 'Pending processes:', len(alignment_results))

            if len(alignment_results) == 0 and len(prots_todo) == 0 and len(chunk) == 0:
                if started_processes != 0 or len(prots_todo) != 0:
                    print('\nWarning: left alignment loop while stuff was still to do', started_processes, prots_todo,'\n')
                break

    else: #serial mode omitting ray

        while True:
            align_outs = alignment_results
            (warn_map, alignment_insertion_list, structure_insertion_list, package_runtimes, integration_times,
                safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes) = process_align_outs(proteins,
                    config, align_outs, warn_map, alignment_insertion_list, structure_insertion_list,
                    indel_analysis_follow_up, safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes, skip_inserts)
            total_package_runtimes = aggregate_times(total_package_runtimes, package_runtimes)
            total_integration_times = aggregate_times(total_integration_times, integration_times)

            n = len(prots_todo)
            if n == 0:
                break
            m = len(protein_packages[list(prots_todo)[0]][1])

            if started_processes < config.proc_n and len(prots_todo) > 0:
                alignment_results, started_processes, protein_packages, cost_map, chunk, chunk_cost, prots_todo = package_alignment_processes(prots_todo, cost_map, optimal_chunk_cost, proteins, protein_packages, alignment_results, started_processes, chunk, chunk_cost, config, mapping_dump)
            else:
                break

            if config.verbosity >= 5:
                print('Proteins to do:', len(prots_todo), 'Chunk:', len(chunk), 'Pending processes:', len(alignment_results))

            if n == len(prots_todo) and n ==1 and m == len(protein_packages[list(prots_todo)[0]][1]):
                config.errorlog.add_error('Alignment todo list not decreased in serial alignment')
                break

    for prot_id in list(protein_packages.keys()):
        del protein_packages[prot_id]
    del protein_packages

    if config.verbosity >= 2:
        print_times(total_package_runtimes, label='Alignment')
        print_times(total_integration_times, label='Align results integration')

    if config.verbosity >= 2:
        t5 = time.time()
        print("Alignment Part 5: %s" % (str(t5 - t4)))

    if len(config.rare_residues) > 0 and config.verbosity >= 3:
        print('######### Detected some rare residues:', config.rare_residues)

    len_of_structure_insertion_list = len(structure_insertion_list)
    if not skip_inserts:    
        database.insertStructures(structure_insertion_list, proteins, config)

    if config.verbosity >= 2:
        t6 = time.time()
        print(f"Alignment Part 6: {t6 - t5} {len_of_structure_insertion_list=}")

    if not skip_inserts:
        database.insertAlignments(alignment_insertion_list, proteins, config)
        
    if config.verbosity >= 2:
        t7 = time.time()
        print("Alignment Part 7: %s" % (str(t7 - t6)))

    # Due the removal of annotations in the previous loop, we might to remove some structures and complexes
    structures_to_remove: list[tuple[str, list[str]]] = []
    for pdb_id in sus_structures:
        if pdb_id not in safe_structures:
            structures_to_remove.append((pdb_id, sus_structures[pdb_id]))
        else:
            rchains: list[str] = []
            for chain in sus_structures[pdb_id]:
                if chain not in safe_structures[pdb_id]:
                    rchains.append(chain)
            if len(rchains) > 0:
                structures_to_remove.append((pdb_id, [x for x in rchains]))

    proteins.remove_structures(structures_to_remove)
    complexes_to_remove = sus_complexes - safe_complexes
    if config.verbosity >= 5:
        print('Len of sus_complexes:', len(sus_complexes), 'Len of safe complexes:', len(safe_complexes), 'Len of complexes_to_remove:', len(complexes_to_remove))
        print('Remove complexes:', complexes_to_remove)
    elif config.verbosity >= 4:
        print('Len of sus_complexes:', len(sus_complexes), 'Len of safe complexes:', len(safe_complexes), 'Len of complexes_to_remove:', len(complexes_to_remove))
        if len(complexes_to_remove) < 50:
            print('Remove complexes:', complexes_to_remove)
    proteins.remove_complexes(complexes_to_remove)

    if config.verbosity >= 2:
        t8 = time.time()
        print("Alignment Part 8: %s" % (str(t8 - t7)))

    del mapping_dump

    if config.verbosity >= 2:
        t9 = time.time()
        print(f"Alignment Part 9: {t9-t8}")

@ray.remote(max_calls = 1)
def paraMap_wrapper(mapping_dump, tasks, static_model_path = None):
    ray_hack()

    config = mapping_dump
    tasks = unpack(tasks)
    return paraMap(config, tasks, static_model_path = static_model_path)

def paraMap(config, tasks, static_model_path = None):
    pdb_path = config.pdb_path

    results = []

    parse_dict: dict[str, tuple] = {}
    structure_set: set[str] = set()
    multi_demand: dict[str, int] = {}

    for _, _, _, struct_tasks in tasks:
        for pdb_id, _, _, _, _, in struct_tasks:
            if pdb_id not in structure_set:
                structure_set.add(pdb_id)
            elif pdb_id not in multi_demand:
                multi_demand[pdb_id] = 1
            else:
                multi_demand[pdb_id] += 1

    for u_ac, aaclist, prot_id, struct_tasks in tasks:
        for pdb_id, chain, structure_id, target_seq, template_seq in struct_tasks:
            try:
                if config.model_db_active:
                    if is_alphafold_model(pdb_id):
                        model_path = alphafold_model_id_to_file_path(pdb_id, config)
                        #is_model = True
                    else:
                        model_path = static_model_path
                else:
                    model_path = static_model_path

                if pdb_id not in parse_dict:
                    template_page, _, _, _, atom_count, _, _ = pdbParser.getStandardizedPdbFile(pdb_id, pdb_path, verbosity=config.verbosity, obsolete_check=True, model_path=model_path)
                    if pdb_id in multi_demand:
                        parse_dict[pdb_id] = template_page, atom_count
                else:
                    template_page, atom_count = parse_dict[pdb_id]
                    multi_demand[pdb_id] -= 1
                    if multi_demand[pdb_id] == 1:
                        del parse_dict[pdb_id]

                seq_res_map, last_residue, first_residue = globalAlignment.createTemplateFasta(template_page, pdb_id, chain.encode('ascii'), config, onlySeqResMap=True)

                sub_out = globalAlignment.getSubPos(config, u_ac, target_seq, template_seq, aaclist, seq_res_map, ignore_gaps = True)

                if isinstance(sub_out, str):
                    config.errorlog.add_error("%s %s %s\n%s\n%s" % (sub_out, pdb_id, chain, template_seq.replace('-', ''), seq_res_map))
                    continue

                sub_infos, aaclist, backmap = sub_out

                results.append((u_ac, pdb_id, chain, sub_infos, atom_count, last_residue, first_residue, backmap))
            except:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                config.errorlog.add_error(f"Error in paraMap:\n{e}\n{f}\n{g}\n{u_ac} ({prot_id}), {pdb_id}:{chain} ({structure_id})\n{target_seq}\n{template_seq}\n{aaclist}")

    return results


#called (also) by indel_analysis.py
@ray.remote(max_calls = 1)
def align_remote_wrapper(align_dump, package, model_path=None):

    config = align_dump

    return pack(align(config, package, static_model_path=model_path, remote = True))


def align(config, package, static_model_path=None, remote=False):
    try:
        results = []
        total_times = []
        aligner_class = globalAlignment.init_bp_aligner_class()

        parse_dict: dict[str, tuple] = {}
        structure_set: set[str] = set()
        multi_demand: dict[str, int] = {}

        for _, structure_infos in package:
            for pdb_id, _, _ in structure_infos:
                if pdb_id not in structure_set:
                    structure_set.add(pdb_id)
                elif pdb_id not in multi_demand:
                    multi_demand[pdb_id] = 1
                else:
                    multi_demand[pdb_id] += 1

        for prot_specific_mapping_dump, structure_infos in package:

            if remote:
                (u_ac, seq, aaclist) = ray.get(prot_specific_mapping_dump, timeout = 600)
            else:
                (u_ac, seq, aaclist) = prot_specific_mapping_dump

            if config.verbosity >= 5:
                print(f'Alignment of {u_ac} to {len(structure_infos)} number of Structures')

            for pdb_id, chain, oligo in structure_infos:
                times = []
                t_init = time.time()

                pdb_path = config.pdb_path
                option_seq_thresh = config.option_seq_thresh

                try:
                    # find and sets model path in given custom db 
                    if config.custom_db_path:
                        custom_path_list = os.listdir(config.custom_db_path)
                        for model in custom_path_list:
                            if model.startswith(pdb_id):
                                model_path = config.custom_db_path + "/%s" %model

                    #is_model = False
                    elif config.model_db_active:
                        if is_alphafold_model(pdb_id):
                            model_path = alphafold_model_id_to_file_path(pdb_id, config)
                            #is_model = True
                        else:
                            model_path = static_model_path
                    else:
                        model_path = static_model_path

                    if pdb_id not in parse_dict:
                        parse_out = pdbParser.getStandardizedPdbFile(pdb_id, pdb_path, oligo=oligo, verbosity=config.verbosity, model_path=model_path)
                        if pdb_id in multi_demand:
                            parse_dict[pdb_id] = parse_out
                    else:
                        parse_out = parse_dict[pdb_id]
                        multi_demand[pdb_id] -= 1
                        if multi_demand[pdb_id] == 1:
                            del parse_dict[pdb_id]

                    t_now = add_to_times(times, t_init)
                    
                    if parse_out is None:
                        results.append((u_ac, pdb_id, chain, 'pdbParser failed for %s' % str(pdb_id)))
                        continue

                    if isinstance(parse_out, str):
                        results.append((u_ac, pdb_id, chain, f'pdbParser failed for {pdb_id}, error message: {parse_out}'))
                        continue

                    (template_page, interaction_partners, chain_type_map, oligo, atom_count, chainlist, rare_residues) = parse_out

                    align_out = globalAlignment.alignBioPython(config, u_ac, seq, pdb_id, template_page, chain, aaclist, rare_residues=rare_residues, aligner_class = aligner_class)
                    
                    t_now = add_to_times(times, t_now)

                    if isinstance(align_out, str):
                        align_out = '%s %s' % (align_out, str(model_path))
                        results.append((u_ac, pdb_id, chain, align_out))
                        continue

                    (coverage, seq_id, sub_infos, backmap, alignment_pir, align_times, aaclist, last_residue, first_residue) = align_out
                    times.append(align_times)

                    if sub_infos is None or seq_id is None:
                        results.append((u_ac, pdb_id, chain, sub_infos, seq_id))
                        continue

                    if 100.0 * seq_id >= option_seq_thresh:
                        results.append((u_ac, pdb_id, chain, alignment_pir, seq_id, coverage, interaction_partners, chain_type_map, oligo, sub_infos, backmap, atom_count, last_residue, first_residue, chainlist, rare_residues))
                    else:
                        if config.verbosity >= 6:
                            align_out = alignment_pir
                        else:
                            align_out = None
                        results.append((u_ac, pdb_id, chain, sub_infos, (seq_id, align_out)))

                except:
                    [e, f, g] = sys.exc_info()
                    g = traceback.format_exc()
                    results.append((u_ac, pdb_id, chain, '%s,%s %s\n%s\n%s\n%s' % (pdb_id, chain, str(model_path), e, str(f), g)))

                total_times = aggregate_times(total_times, times)

    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        results = [('', '', '', f'\nUnknown critical error in alignment: \n{e}\n{str(f)}\n{g}\n')]
        total_times = []

    return (results, total_times)


def debug_freeze(config, stop_number = 1):
    print('\n\n', config.ray_init_counter, '\n\n')
    if config.ray_init_counter == stop_number:
        time.sleep(999999999)
    return

def print_locals(items):
    for name, size in sorted(((name, sys.getsizeof(value)) for name, value in items), key=lambda x: -x[1])[:10]:
        print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))

def core(
        protein_list: dict[str, protein_package.Protein],
        genes,
        indels,
        multi_mutation_objects,
        config,
        session
        ):
    
    if len(protein_list) == 0:
        return None, 0

    if config.proc_n > 1:
        ray_init(config)
        config.ray_init_counter += 1    
        #monitor_ray_store()

    background_process_MS = None

    # transform the protein map into a Proteins object
    proteins: protein_package.Proteins = protein_package.Proteins(protein_list, indels, multi_mutation_objects)


    if config.verbosity >= 6:
        print('Proteins state after object initialization:')
        proteins.print_protein_state()

    try:
        t1 = time.time()

        if config.verbosity >= 1:
            print("Before protCheck")
        # check for already stored genes and make all the necessary database interactions
        # sets the fields database_id and stored in the protein objects as well as the stored and non stored ids in the proteins object

        database.protCheck(proteins, genes, session, config)

        if not config.only_wt:
            for wt_prot_id in proteins.multi_mutations:
                if config.verbosity >= 6:
                    print('Mutating', len(proteins.multi_mutations[wt_prot_id]), 'multi mutations for:', wt_prot_id)
                for multi_mutation in proteins.multi_mutations[wt_prot_id]:
                    #if proteins.contains(multi_mutation.mut_prot):
                    #    if config.verbosity >= 3:
                    #        print(f'Mutating for {multi_mutation.mut_prot} not necessary, skip.')
                    #    continue
                    multi_mutation.mutate(proteins, config)

        if config.verbosity >= 6:
            print('Proteins state after protCheck:')
            proteins.print_protein_state()

        t2 = time.time()
        if config.verbosity >= 2:
            print("Time for protCheck: %s" % (str(t2 - t1)))
        # check for already stored mutations or position twins and make all the necessary database interactions
        # sets the fields database_id and stored in the position objects

        
        background_process_MS = database.positionCheck(proteins, session, config)

        t3 = time.time()
        if config.verbosity >= 2:
            print(f"Time for positionCheck: {t3 - t2}")

        database.indelCheck(proteins, session, config)

        t4 = time.time()
        if config.verbosity >= 2:
            print(f"Time for indelCheck: {t4 - t3}")

        database.insertMultiMutations(proteins, genes, session, config)

        if config.verbosity >= 6:
            print('Proteins state after positionCheck:')
            proteins.print_protein_state()

        t5 = time.time()
        if config.verbosity >= 2:
            print(f"Time for insertMultiMutations: {t5 - t4}")

        if config.verbosity >= 1:
            print("Before getSequences")

        getSequences(proteins, config)

        t6 = time.time()
        if config.verbosity >= 2:
            print(f"Time for getSequences: {t6 - t5}")


        if config.verbosity >= 1:
            print("Before autoTemplateSelection")
        autoTemplateSelection(config, proteins)

        t7 = time.time()
        if config.verbosity >= 2:
            print(f"Time for Template selection: {t7 - t6}")


        indel_analysis_follow_up = len(proteins.indels) > 0 and not config.skip_indel_analysis
        
        if config.verbosity >= 1:
            print("Before paraAlignment")
        paraAlignment(config, proteins, indel_analysis_follow_up=indel_analysis_follow_up)

        t8 = time.time()
        if config.verbosity >= 2:
            print(f"Time for Alignment: {t8 - t7}")

        if config.verbosity >= 1:
            print("Before paraAnnotate")

        amount_of_structures = templateFiltering.paraAnnotate(config, proteins, indel_analysis_follow_up=indel_analysis_follow_up)
        
        t9 = time.time()
        if config.verbosity >= 2:
            print(f"Time for paraAnnotate: {t9 - t8}")


        if not config.skip_indel_analysis:
            indel_analysis.para_indel_analysis(proteins, config)

        t10 = time.time()
        if config.verbosity >= 2:
            print(f"Time for indel analysis: {t10 - t9}")

                
        proteins.deconstruct()
        del proteins

        t11 = time.time()
        if config.verbosity >= 2:
            print(f"Time for deconstruction: {t11 - t10}")

        # join the background inserts
        if background_process_MS is not None:
            background_process_MS.join()
            background_process_MS.close()
            background_process_MS = None

        t12 = time.time()
        if config.verbosity >= 2:
            print(f"Time for background inserts: {t12 - t11}")


        if config.verbosity >= 3:
            for name, size in sorted(((name, sys.getsizeof(value)) for name, value in locals().items()), key=lambda x: -x[1])[:10]:
                print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
            print(list_fds())

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

    if config.ray_init_counter >= 10000:
        ray.shutdown()
        #time.sleep(10)
        config.ray_init_counter = 0
        pass
    return amount_of_structures		


def main(filename, config):
    # main function takes another argument custom db, default is False
    mrna_fasta = config.mrna_fasta

    session = 0  # This can later be used, structman_main.py could give a specific session id and the pipeline can then expand that session

    if config.verbosity >= 2:
        print('==== Start main pipeline ====')

    config.errorlog.start(filename, config.outfolder)

    t0 = time.time()

    if config.proc_n > 1:
        # need structman package path for ray
        ray_init(config)

        if config.verbosity >= 2:
            print('ray init successful')

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
            print('=== Single line input mode ===')

    t01 = time.time()

    if config.verbosity >= 2:
        print("Time for preparation before buildQueue: %s" % (str(t01 - t0)))

    try:
        os.stat("%s/tmp_structman_pipeline" % (config.outfolder))
    except:
        os.mkdir("%s/tmp_structman_pipeline" % (config.outfolder))
    os.chdir("%s/tmp_structman_pipeline" % (config.outfolder))
    cwd = "%s/tmp_structman_pipeline" % (config.outfolder)

    config.temp_folder = cwd

    chunksize = config.chunksize

    if config.verbosity >= 1:
        print("Call buildQueue with chunksize: %s and file: %s" % (str(chunksize), nfname))

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
        print("Time for buildQueue: %s" % (str(t02 - t01)))
    if config.verbosity >= 1:
        print("Number of chunks: ", len(proteins_chunks))

    newsession = False
    if session == 0:
        starttime = SQLDateTime()
        session = database.insertSession(starttime, nfname, config)
        newsession = True
    session_name = (nfname.rsplit("/", 1)[1]).rsplit(".", 1)[0]

    total_amount_of_analyzed_structures = 0
    config.ray_init_counter = 0
    save_ids = None

    for nr_temp_file, temp_infile in enumerate(temp_infiles):
        if temp_infile is not None:
            if config.verbosity >= 1:
                print('Infile splitting due to low memory system, processing infile split nr.:', nr_temp_file + 1, 'out of', len(temp_infiles))
            proteins_chunks: list[tuple[dict[str, protein_package.Protein], dict[str, indel_package.Indel]]]
            proteins_chunks, _, genes = buildQueue(config, temp_infile, already_split=True)
            os.remove(temp_infile)

        #The gene objects need to exist as a whole outside of the chunked core loop. Isoforms that correspond to the same gene can be chunked into different bins.
        database.geneCheck(genes, config)

        chunk_nr = 1
        for protein_list, indels in proteins_chunks:

            tc0 = time.time()
            config.indels_given_by_input = (len(indels) > 0)

            if config.verbosity >= 1:
                print("Chunk %s/%s" % (str(chunk_nr), str(len(proteins_chunks))))
            chunk_nr += 1

            if config.verbosity >= 2:
                t0 = time.time()
            protein_list, indels, multi_mutation_objects, removed_proteins, _ = sequenceScan(config, protein_list, indels)
            for removed_prot_id in removed_proteins:
                del protein_list[removed_prot_id]
            if config.verbosity >= 2:
                t1 = time.time()
                print(f'Time for sequenceScan: {t1-t0}')

            amount_of_structures = core(protein_list, genes, indels, multi_mutation_objects, config, session)

            total_amount_of_analyzed_structures += amount_of_structures

            for prot_id in list(protein_list.keys()):
                del protein_list[prot_id]

            for indel_id in list(indels.keys()):
                del indels[indel_id]

            if config.verbosity >= 3:
                try:
                    for name, size in sorted(((name, sys.getsizeof(value)) for name, value in globals().items()), key=lambda x: -x[1])[:10]:
                        print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
                    print(list_fds())
                except:
                    pass

            save_ids = kill_overhead_sockets(save_ids = save_ids)

            if config.verbosity >= 3:
                print(f'save_ids: {save_ids}')
                print(list_fds())

            tc1 = time.time()
            if config.verbosity >= 1:
                print(f'\nChunk time: {tc1-tc0}\n')

    os.chdir(config.outfolder)

    if total_amount_of_analyzed_structures == 0:
        if config.verbosity >= 1:
            print("=== Nothing got processed, probably the input didn't contain any usable protein ID. ===")

    if config.verbosity >= 2:
        t03 = time.time()

    try:
        shutil.rmtree(cwd)
    except:
        pass

    if config.verbosity >= 2:
        t04 = time.time()
        print('Time for folder cleanup:', t04 - t03)

    ray.shutdown()

    returncode = config.errorlog.stop()

    if newsession:
        endtime = SQLDateTime()
        database.updateSession(session, endtime, config)

    tend = time.time()
    if config.verbosity >= 1:
        print(total_amount_of_analyzed_structures, 'structures in total got analyzed')
        print('Total runtime of StructMAn:', (tend - t0))

    '''
    total_memory_peak = 0
    for p in mem_tracked_processes:
        total_memory_peak += p.memory_info().rss

    total_memory_peak = total_memory_peak/1024./1024./1024.

    print('Accumulated memory peak',total_memory_peak,'Gb')
    '''

    return session, config, returncode
