import errno
import gc
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
import statistics

import ray

from structman.lib import annovar
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

from structman.lib.database import database
from structman.lib.output import output
from structman.base_utils.base_utils import calculate_chunksizes, pack, unpack, is_alphafold_model, alphafold_model_id_to_file_path
from structman.base_utils.ray_utils import ray_init, ray_hack
from structman.lib.lib_utils import process_mutations_str, generate_multi_mutation, extract_dedicated_tags

try:
    from memory_profiler import profile
except:
    pass


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
    for pos in proteins[primary_protein_id].positions:
        if not proteins[primary_protein_id].positions[pos].checked:
            del_list.append((pos, proteins[primary_protein_id].positions[pos].wt_aa))
    for pos, wt_aa in del_list:
        remove_position(config, proteins, indels, pos, wt_aa, primary_protein_id, degrade_to_session_less = True)


def remove_position(config, proteins, indels, pos, wt_aa, primary_protein_id, degrade_to_session_less = False):
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

# @profile


def sequenceScan(config, proteins, indels):
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
            print(f'Protein IDs going into sequence scan: {proteins}')

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
            elif proteins[prot_id].input_id[:3] == 'ENS':
                sequenceScanEnsembl[proteins[prot_id].input_id] = tags, uni_pos
            else:
                if prot_id in indels:  # Indel proteins need complete classification
                    uni_pos = True
                sequenceScanProteins[prot_id] = tags, uni_pos  # New: process everything to filter input by sanity checks

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
                if seq_pos in proteins[primary_protein_id].positions:
                    proteins[primary_protein_id].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[primary_protein_id].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if seq_pos not in proteins[primary_protein_id].positions:
                        position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True, session_less = session_less)
                        proteins[primary_protein_id].positions[seq_pos] = position

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
            seq = gene_sequence_map[nm_ref]
            if nm_ref in isoform_specific_id_map:
                primary_protein_id = promote_uac_to_primary_id(isoform_specific_id_map[nm_ref], nm_ref)
            else:
                primary_protein_id = nm_ref

            proteins[primary_protein_id].sequence = seq

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if seq_pos in proteins[primary_protein_id].positions:
                    proteins[primary_protein_id].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[primary_protein_id].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if seq_pos not in proteins[primary_protein_id].positions:
                        position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True, session_less = session_less)
                        proteins[primary_protein_id].positions[seq_pos] = position

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, primary_protein_id)

    if len(sequenceScanEnsembl) > 0:
        if config.verbosity >= 1:
            print("Amount of proteins going into Ensembl sequenceScan: ", len(sequenceScanEnsembl))
        
        gene_sequence_map = uniprot.get_ensembl_seqs(config, list(sequenceScanEnsembl.keys()))
        for transcript_id in sequenceScanEnsembl:
            tags, uni_pos = sequenceScanEnsembl[transcript_id]
            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False

            if transcript_id not in gene_sequence_map:
                #This happens for problematic Ensembl transcript ids, for example for entries without an open reading frame (due to the transcript being non-coding)
                try:
                    proteins.remove_protein(transcript_id)
                except:
                    del proteins[transcript_id]
                if config.verbosity >= 3:
                    print(f'Removed {transcript_id} from Input')
                removed_proteins.append(transcript_id)
                continue

            seq = gene_sequence_map[transcript_id]

            primary_protein_id = transcript_id

            proteins[primary_protein_id].sequence = seq

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if seq_pos in proteins[primary_protein_id].positions:
                    proteins[primary_protein_id].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[primary_protein_id].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if seq_pos not in proteins[primary_protein_id].positions:
                        position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True, session_less = session_less)
                        proteins[primary_protein_id].positions[seq_pos] = position

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, primary_protein_id)

    if len(sequenceScanProteins) > 0:
        if config.verbosity >= 1:
            print("Amount of proteins going into Uniprot sequenceScan: ", len(sequenceScanProteins))

        gene_sequence_map = uniprot.getSequencesPlain(sequenceScanProteins.keys(), config)

        #print_locals(locals().items())

        for u_ac in gene_sequence_map:
            if gene_sequence_map[u_ac][0] == 1 or gene_sequence_map[u_ac][0] == 0 or gene_sequence_map[u_ac][0] is None:
                config.errorlog.add_warning("Error in sequenceScan with gene: %s" % u_ac)
                continue
            seq, disorder_scores, disorder_regions_datastruct = gene_sequence_map[u_ac]
            proteins[u_ac].sequence = seq

            tags, uni_pos = sequenceScanProteins[u_ac]

            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False
            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if seq_pos in proteins[u_ac].positions:
                    proteins[u_ac].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[u_ac].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if seq_pos not in proteins[u_ac].positions:
                        position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True, session_less = session_less)
                        proteins[u_ac].positions[seq_pos] = position

            proteins[u_ac].set_disorder_scores(disorder_scores)
            proteins[u_ac].set_disorder_regions(disorder_regions_datastruct)

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, u_ac)
            del seq

        #print_locals(locals().items())

        del gene_sequence_map
        #print_locals(locals().items())
    
    if len(sequenceScanPDB) > 0:
        pdb_sequence_map, pdb_pos_map = pdbParser.getSequences(sequenceScanPDB.keys(), pdb_path)
        for u_ac in pdb_sequence_map:
            proteins[u_ac].sequence = pdb_sequence_map[u_ac][0]
            tags, uni_pos = sequenceScanPDB[u_ac]
            if tags is None:
                session_less = True
                tags = set()
            else:
                session_less = False
            residue_id_backmap = {}
            res_pos_map = pdb_pos_map[u_ac]

            for res_id in res_pos_map:
                residue_id_backmap[res_pos_map[res_id]] = res_id
            for (pos, aa) in enumerate(pdb_sequence_map[u_ac][0]):
                seq_pos = pos + 1
                if seq_pos not in residue_id_backmap:
                    continue
                res_id = residue_id_backmap[seq_pos]
                if res_id in proteins[u_ac].res_id_map:
                    proteins[u_ac].res_id_map[res_id].pos = seq_pos
                    proteins[u_ac].positions[seq_pos] = proteins[u_ac].res_id_map[res_id]
                    if tags is None:
                        continue
                    proteins[u_ac].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if seq_pos not in proteins[u_ac].positions:
                        position = position_package.Position(pos=seq_pos, pdb_res_nr=res_id, wt_aa=aa, tags=tags, session_less = session_less)
                        proteins[u_ac].positions[seq_pos] = position
                        proteins[u_ac].res_id_map[res_id] = position
                    elif not res_id in proteins[u_ac].res_id_map:
                        proteins[u_ac].res_id_map[res_id] = proteins[u_ac].positions[seq_pos]

    if len(fasta_inputs) > 0:
        for prot_id in fasta_inputs:
            seq = proteins[prot_id].sequence

            tags = fasta_inputs[prot_id]

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if seq_pos not in proteins[prot_id].positions:
                    position = position_package.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True)
                    proteins[prot_id].positions[seq_pos] = position

    if config.verbosity >= 3:
        print('Before indel mutation with:', len(indels), 'number of indels')

    for primary_protein_id in indels:
        for indel_notation in indels[primary_protein_id]:
            indels[primary_protein_id][indel_notation].mutate_sequence(proteins)

    multi_mutation_objects = []

    if not config.only_wt:
        if config.verbosity >= 3:
            print('Before multi mutations mutation')
        for primary_protein_id in list(proteins.keys()).copy():  # iterate over snapshot, since the mutant proteins get added to proteins in the loop
            if proteins[primary_protein_id].multi_mutations is not None:
                # creates mutation protein object
                multi_mutation_objects += proteins[primary_protein_id].create_multi_mutations(proteins, config)

    return proteins, indels, multi_mutation_objects, removed_proteins


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

    for line in lines:
        line = line[:-1]
        if len(line) == 0:
            continue
        if line[0] == '>':
            words = line[1:].split()
            entry_id = words[0]
            if entry_id.count('|') > 1:
                entry_id = entry_id.split('|')[1]
            seq_map[entry_id] = ''
            if len(words) > 1:
                aacs = words[1]
                aac_tag_tuples = aacs.split(';')
                mm = []
                for aac_tag_tuple in aac_tag_tuples:
                    if aac_tag_tuple.count('<') == 1:
                        aac, tags = aac_tag_tuple.split('<')
                    else:
                        aac = aac_tag_tuple
                        tags = set()

                    pr_mut_str_result = process_mutations_str(config, aac, tags)

                    if pr_mut_str_result is None:
                        config.errorlog.add_warning(f'Couldnt process mutation_str: {aac}, {tags}')
                        continue

                    positions, new_multi_mutations, gene_id, protein_specific_tags = pr_mut_str_result
                    mm += new_multi_mutations
                    if gene_id is not None:
                        prot_gene_map[entry_id] = gene_id

                    prot_tags_map[entry_id] = protein_specific_tags

                multi_mutations.append(mm)
        elif line[0] == '<':
            words = line[1:].split()

            if len(words) == 2:
                aac, tags = words
            else:
                aac = words[0]
                tags = set()
            pr_mut_str_result = process_mutations_str(config, aac, tags)

            if pr_mut_str_result is None:
                config.errorlog.add_warning(f'Couldnt process mutation_str: {aac}, {tags}')
                continue

            new_positions, new_multi_mutations, gene_id, protein_specific_tags = pr_mut_str_result

            positions += new_positions
            multi_mutations.append(new_multi_mutations)
            if gene_id is not None:
                prot_gene_map[entry_id] = gene_id
            prot_tags_map[entry_id] = protein_specific_tags

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

        uniprot.integrate_protein(config, proteins, genes, indel_map, prot_id, prot_id, prot_map, gene_id = gene_id, protein_specific_tags = protein_specific_tags)
        proteins[prot_id].sequence = seq

    gene_isoform_check(proteins, genes, indel_map, config)

    outlist = input_chunking(config, proteins, indel_map)

    return outlist, genes


def add_to_prot_map(prot_map, sp_id, positions, multi_mutations, config):
    if sp_id not in prot_map:
        prot_map[sp_id] = [], [], []
    for multi_mutation in multi_mutations:
        for indel_or_snv in multi_mutation:
            if not isinstance(indel_or_snv, tuple):
                prot_map[sp_id][1].append(indel_or_snv)
            else:
                prot_map[sp_id][0].append(indel_or_snv[0])

    for position in positions:
        prot_map[sp_id][0].append(position)
    for multi_mutation in multi_mutations:
        if len(multi_mutation) > 1 and not config.only_snvs:
            prot_map[sp_id][2].append(multi_mutation)
    return prot_map


def input_chunking(config, proteins, indels):
    outlist = []

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


# @profile
def buildQueue(config, filename, already_split=False):
    t0 = time.time()

    proteins = {}
    u_ids = set()
    u_acs = set()
    id_map = {}
    ac_map = {}
    np_map = {}
    nm_map = {}
    hgnc_map = {}
    ensembl_map = {}

    pdb_map = {}

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

        tags = set()

        if config.verbosity >= 5:
            print(f'Processing input line: {line}')

        if len(words) > 2:
            tags = set(words[2].split(','))

        try:
            if len(words) == 1 or words[1] == '':
                position = position_package.Position(tags=tags)
                indel = None
                multi_mutations = []
                positions = [position]

                tags, gene_id = extract_dedicated_tags(tags)
                if gene_id is not None:
                    prot_gene_map[sp_id] = gene_id
                prot_tags_map[sp_id] = tags

            else:
                mutation_str = words[1].replace("\n", "")
                multi_mutations = []
                positions = []
                if (not sp_id.count(':') == 1) or sp_id[0:5] == 'HGNC:':  # this means sp_id is not a pdb-id
                    positions, multi_mutation, gene_id, protein_specific_tags = process_mutations_str(config, mutation_str, tags)

                    multi_mutations.append(multi_mutation)
                    if gene_id is not None:
                        prot_gene_map[sp_id] = gene_id
                    prot_tags_map[sp_id] = protein_specific_tags

                else:  # this means sp_id is a pdb-id
                    positions, multi_mutations, gene_id, protein_specific_tags = process_mutations_str(config, mutation_str, tags, pdb_style=True)
                    multi_mutations.append(multi_mutation)
                    if gene_id is not None:
                        pdb_chain_tuple = '%s:%s' % (sp_id[:4].upper(), sp_id[-1])  # enforce uppercase pdb-id
                        prot_gene_map[pdb_chain_tuple] = gene_id
                    prot_tags_map[sp_id] = protein_specific_tags

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
            else:
                u_acs.add(sp_id)
                ac_map = add_to_prot_map(ac_map, sp_id, positions, multi_mutations, config)

    t1 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 1: ", str(t1 - t0), len(ac_map))

    proteins, indel_map, genes = uniprot.IdMapping(config, ac_map, id_map, np_map, pdb_map, hgnc_map, nm_map, ensembl_map, prot_gene_map, prot_tags_map)

    gene_isoform_check(proteins, genes, indel_map, config)

    t2 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 2: ", str(t2 - t1))

    outlist = input_chunking(config, proteins, indel_map)

    t3 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 3: ", str(t3 - t2))

    return outlist, [None], genes


def gene_isoform_check(proteins, genes, indel_map, config):
    extra_sequence_scan_proteins = {}
    for gene_id in genes:
        if len(genes[gene_id].proteins) < 2:
            continue
        isoform_list = list(genes[gene_id].proteins)

        for iso1 in isoform_list:
            for iso2 in isoform_list:
                if iso1 == iso2:
                    continue
                extra_sequence_scan_proteins[iso1] = proteins[iso1]
                extra_sequence_scan_proteins[iso2] = proteins[iso2]

    extra_sequence_scan_proteins, indel_map, _, removed_proteins = sequenceScan(config, extra_sequence_scan_proteins, indel_map)

    for prot_id in proteins:
        if prot_id in extra_sequence_scan_proteins:
            proteins[prot_id] = extra_sequence_scan_proteins[prot_id]

    for removed_prot_id in removed_proteins:
        del proteins[removed_prot_id]

    packages = []
    current_package_number = 0
    config_dump = ray.put(config)

    for gene_id in genes:
        if len(genes[gene_id].proteins) < 2:
            continue
        isoform_list = list(genes[gene_id].proteins)

        sequence_map = {}
        iso_pairs = []

        for iso1 in isoform_list:
            if iso1 not in proteins:
                if iso1 in genes[gene_id].proteins:
                    genes[gene_id].proteins.remove(iso1)
                continue
            for iso2 in isoform_list:
                if iso2 not in proteins:
                    if iso2 in genes[gene_id].proteins:
                        genes[gene_id].proteins.remove(iso2)
                    continue
                if iso1 == iso2:
                    continue

                if config.verbosity >= 4:
                    print(f'Generating Isoform-Multimutation for {iso1} and {iso2}')

                if iso1 not in sequence_map:
                    sequence_map[iso1] = proteins[iso1].sequence
                if iso2 not in sequence_map:
                    sequence_map[iso2] = proteins[iso2].sequence
                iso_pairs.append((iso1, iso2))

        if len(iso_pairs) == 0:
            continue

        if current_package_number >= len(packages):
            if len(packages) >= config.proc_n:
                current_package_number = 0
            else:
                packages.append([])

        packages[current_package_number].append((sequence_map, iso_pairs))
        current_package_number += 1

    mm_generation_remote_ids = []
    for package in packages:
        mm_generation_remote_ids.append(para_mm_generation.remote(config_dump, pack(package)))

    isoform_mm_matrix = {}

    for mm_generation_results in ray.get(mm_generation_remote_ids):
        for iso1, iso2, mm_results in unpack(mm_generation_results):
            if isinstance(mm_results, str):
                config.errorlog.add_warning(f'Multi mutation generation failed for {iso1} and {iso2} with error:\n{mm_results}')
                continue

            positions, multi_mutation = mm_results

            for indel_or_snv in multi_mutation:
                if isinstance(indel_or_snv, tuple):
                    if config.verbosity >= 4:
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

    return

@ray.remote(max_calls = 1)
def para_mm_generation(config, packages):
    mm_results = []
    aligner_class = globalAlignment.init_bp_aligner_class(open_gap_score = -2.0)
    for sequence_map, iso_pairs in unpack(packages):
        for iso1, iso2 in iso_pairs:
            mm_results.append((iso1, iso2, generate_multi_mutation(sequence_map[iso1], sequence_map[iso2], config, aligner_class = aligner_class)))
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


# @profile
def getSequences(proteins, config):

    t0 = time.time()

    uniprot.getSequences(proteins, config)

    t1 = time.time()
    if config.verbosity >= 2:
        print("Time for getSequences Part 1: %s" % str(t1 - t0))

    database.addProtInfos(proteins, config)

    t2 = time.time()
    if config.verbosity >= 2:
        print("Time for getSequences Part 2: %s" % str(t2 - t1))

    u_acs = proteins.get_protein_ids()

    store = None

    background_iu_process = None
    if config.iupred_path != '':
        iupred_path = config.iupred_path
        if iupred_path.count('mobidb-lite') > 0:
            mobi_lite = True
        else:
            mobi_lite = False

        t0 = time.time()

        if not mobi_lite:
            iupred_results = []
            sys.path.append(os.path.abspath(os.path.realpath(config.iupred_path)))
            if config.proc_n > 1:
                store = ray.put(config)
        else:
            mobi_list = []

        n_disorder = 0
        stored_disorder_ids = []

        if config.verbosity >= 3:
            print('Start disordered region calculations')

        for u_ac in u_acs:
            seq = proteins.get_sequence(u_ac)
            if seq == 0 or seq == 1 or seq is None or len(seq) <= 20:
                if config.verbosity >= 3:
                    print('Disordered regions calculations skipped for:', u_ac, 'since sequence is:', seq)
                continue
            disorder_scores = proteins.get_disorder_scores(u_ac)
            if disorder_scores is None:
                if not mobi_lite:
                    if config.verbosity >= 5:
                        print('Start disordered region calculation for:', u_ac)
                    if config.proc_n > 1:
                        iupred_results.append(para_iupred.remote(u_ac, seq, store))
                    else:
                        iupred_results.append(iupred(u_ac, seq, config))
                    n_disorder += 1
                else:
                    mobi_list.append('>%s\n%s\n' % (u_ac, seq))
            elif disorder_scores != 'Stored':
                proteins.set_disorder_tool(u_ac, 'MobiDB3.0')

        if not mobi_lite:
            if config.proc_n > 1:
                iupred_out = ray.get(iupred_results)
            else:
                iupred_out = iupred_results

            for iupred_parts in iupred_out:
                proteins.set_disorder_scores(iupred_parts[0], iupred_parts[2])
                proteins.set_disorder_regions(iupred_parts[0], iupred_parts[1])
                proteins.set_disorder_tool(iupred_parts[0], 'IUpred')

        else:
            mobi_tmp_file = 'mobi_tmp_file.fasta'
            f = open(mobi_tmp_file, 'w')
            f.write(''.join(mobi_list))
            f.close()
            mobi_bin_path = '%s/binx/' % iupred_path.rsplit('/', 1)[0]
            mobi_threads = min([7, config.proc_n])
            p = subprocess.Popen([iupred_path, mobi_tmp_file, '-t', str(mobi_threads), '-bin', mobi_bin_path, '-l'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            out, err = p.communicate()
            os.remove(mobi_tmp_file)
            if err != '':
                config.errorlog.add_warning('Warning: mobidb-lite threw an error: %s' % err)
            else:
                entries = []
                for line in out.split('}'):
                    if line.count('{') == 0:
                        continue
                    line = line + '}'
                    mobi_results = json.loads(line)
                    u_ac = mobi_results['acc']
                    raw_scores = mobi_results['p']
                    raw_regions = mobi_results['regions']
                    regions = []
                    for a, b in raw_regions:
                        regions.append([a, b, 'disorder'])
                    scores = {}
                    #seq = proteins.get_sequence(u_ac)
                    #disorder_scores = proteins.get_disorder_scores(u_ac)
                    #disorder_regions = proteins.get_disorder_regions(u_ac)
                    for pos, score in enumerate(raw_scores):
                        scores[pos + 1] = score
                    proteins.set_disorder_scores(u_ac, scores)
                    proteins.set_disorder_regions(u_ac, regions)
                    proteins.set_disorder_tool(u_ac, 'mobidb-lite')

        t1 = time.time()
        if config.verbosity >= 2:
            print(f"Time for addIupred Part 1: {t1 - t0}")

        background_iu_process = database.addIupred(proteins, config)

        t2 = time.time()
        if config.verbosity >= 2:
            print("Time for addIupred Part 2: %s" % str(t2 - t1))

    del store

    return background_iu_process


@ray.remote(max_calls=1)
def para_iupred(u_ac, seq, store):
    ray_hack()
    config = store
    return iupred(u_ac, seq, config)

def iupred(u_ac, seq, config):
    import iupred3_lib
    try:
        iupred2_result = iupred3_lib.iupred(seq, 'glob', 'medium')
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        config.errorlog.add_error("IUPred3 Error: %s\n%s\n%s\n%s\n%s" % (u_ac, seq, str(e), str(f), str(g)))
        return [u_ac, [], []]
    #anchor2_res = iupred3_lib.anchor2(seq)

    iupred_parts = [u_ac, [], {}]

    lines = iupred2_result[1].split('\n')
    for line in lines:
        words = line.split()
        if len(words) < 4:
            continue
        if words[0] == 'globular' and words[1] == 'domain':
            lower_bound, upper_bound = words[3].split('-')
            iupred_parts[1].append((int(lower_bound), int(upper_bound), 'globular'))

    for pos, residue in enumerate(seq):
        iupred_parts[2][pos+1] = iupred2_result[0][pos]

    return iupred_parts


# @profile
def autoTemplateSelection(config, proteins):

    if config.verbosity >= 1:
        print('Sequence search with mmseqs2')

    filtering_dump = None

    # Checks if custom_db path is given in config file, and sets custom_db boolean value
    if config.custom_db_path:
        custom_db = True
    else:
        custom_db = False
    t0 = time.time()

    # MMseqs search function takes one more argument now: custom_db
    search_results = MMseqs2.search(proteins, config, custom_db)

    t1 = time.time()

    if search_results is not None:
        for search_result in search_results:
            if len(search_result) == 0:
                continue
            raw_structure_map, pdb_ids, is_model_db = search_result

            info_map = {}
            filtering_results = []

            # if it is custom db, it creates fake oligo map
            if custom_db:
                for custom_id in pdb_ids:
                    fake_oligo_map = {'A' : ['A']}
                    info_map[custom_id] = (2.0, fake_oligo_map)

            elif not is_model_db:
                
                # calls filtering dump with custom db path rather than pdb path
                if config.proc_n > 1:
                    if custom_db:
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
                    fake_oligo_map = {'A' : ['A']}
                    info_map[model_id] = (2.0, fake_oligo_map)

            for u_ac in raw_structure_map:
                for pdb_id, chain in raw_structure_map[u_ac]:

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
                        oligo = raw_structure_map[u_ac][(pdb_id, chain)][2]

                    struct_anno = structure_package.StructureAnnotation(u_ac, pdb_id, chain)
                    proteins.add_annotation(u_ac, pdb_id, chain, struct_anno)

                    if config.verbosity >= 5:
                        print(f'Adding structural annotation: {u_ac} -> {pdb_id}:{chain}')

                    if not (pdb_id, chain) in proteins.structures:
                        struct = structure_package.Structure(pdb_id, chain, oligo=oligo, mapped_proteins=[u_ac], seq_len = raw_structure_map[u_ac][(pdb_id, chain)][4])
                        proteins.add_structure(pdb_id, chain, struct)
                    else:
                        proteins.add_mapping_to_structure(pdb_id, chain, u_ac)

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
    ray_hack()
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
    for prot_id in prots_todo:
        package_cost = cost_map[prot_id]
        if package_cost > optimal_chunk_cost:  #packages that are greater than the optimal_chunk_cost get split
            structure_infos_split_a = []
            structure_infos_split_b = []
            seq = proteins.get_sequence(prot_id)
            (prot_specific_mapping_dump, structure_infos) = protein_packages[prot_id]
            split_cost = 0
            for structure_info in structure_infos:
                pdb_id, chain, oligo = structure_info
                str_len = proteins.structures[(pdb_id, chain)].get_seq_len()
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
            else:
                alignment_results.append(align(config, chunk))
            if config.verbosity >= 6:
                print('A, Start alignment package with cost:', split_cost)

            started_processes += 1
            if len(structure_infos_split_b) > 0:
                protein_packages[prot_id] = (prot_specific_mapping_dump, structure_infos_split_b)
                cost_map[prot_id] = package_cost - split_cost
            else:
                done.append(prot_id)
        elif chunk_cost + package_cost <= optimal_chunk_cost:
            chunk.append(protein_packages[prot_id])
            done.append(prot_id)
            chunk_cost += package_cost
            if chunk_cost >= optimal_chunk_cost:
                if config.proc_n > 1:
                    alignment_results.append(align_remote_wrapper.remote(mapping_dump, chunk))
                else:
                    alignment_results.append(align(config, chunk))
                if config.verbosity >= 6:
                    print('B, Start alignment package with cost:', chunk_cost)
                chunk = []
                chunk_cost = 0
                started_processes += 1
        
        if started_processes >= config.proc_n:
            break

    if (len(chunk) > 0) and len(prots_todo) == 0:
        if config.proc_n > 1:
            alignment_results.append(align_remote_wrapper.remote(mapping_dump, chunk))
        else:
            alignment_results.append(align(config, chunk))
        if config.verbosity >= 6:
            print('C, Start alignment package with cost:', chunk_cost)
        chunk = []
        chunk_cost = 0
        started_processes += 1

    if (len(chunk) > 0) and package_cost is not None:
        if (chunk_cost + package_cost > optimal_chunk_cost)  and (started_processes < config.proc_n):
            if config.proc_n > 1:
                alignment_results.append(align_remote_wrapper.remote(mapping_dump, chunk))
            else:
                alignment_results.append(align(config, chunk))
            if config.verbosity >= 6:
                print('D, Start alignment package with cost:', chunk_cost)
            chunk = []
            chunk_cost = 0
            started_processes += 1
        
    for prot_id in done:
        prots_todo.remove(prot_id)

    return alignment_results, started_processes, protein_packages, cost_map, chunk, chunk_cost, prots_todo

def process_align_outs(proteins, config, align_outs, max_runtime, warn_map, alignment_insertion_list, structure_insertion_list,
                        package_runtimes, max_runtime_package, min_runtime_package, record_package, indel_analysis_follow_up,
                        safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes, skip_inserts):


    record_min_package = False
    if config.verbosity >= 4:
        record_min_package = True

    t_i_0 = 0.
    t_i_1 = 0.
    t_i_2 = 0.
    t_i_3 = 0.
    t_i_4 = 0.
    t_i_5 = 0.
    t_i_6 = 0.
    t_i_7 = 0.
    t_i_8 = 0.
    t_i_9 = 0.
    for out_chunks in align_outs:
        if config.proc_n > 1:
            out_chunks, package_runtime = unpack(out_chunks)
        else:
            out_chunks, package_runtime = out_chunks
        if config.verbosity >= 4:
            package_runtimes.append(package_runtime)
            if package_runtime > max_runtime:
                max_runtime = package_runtime
                record_package = True
                max_runtime_package = []

        for out in out_chunks:
            t_i_0 += time.time()
            if len(out) == 3:
                config.errorlog.add_error('Illegal alignment output: %s' % (str(out)))
                t_i_1 += time.time()
                continue

            if len(out) == 4:
                (prot_id, pdb_id, chain, warn_text) = out
                proteins.remove_annotation(prot_id, pdb_id, chain)
                if not skip_inserts:
                    sus_complexes.add(pdb_id)
                    sus_structures.add((pdb_id, chain))
                if prot_id not in warn_map:
                    config.errorlog.add_warning(warn_text)
                elif config.verbosity >= 4:
                    print('Alignment failed:', prot_id, pdb_id, chain, warn_text)
                warn_map.add(prot_id)
                t_i_1 += time.time()
                continue

            if len(out) == 5:
                (prot_id, pdb_id, chain, sub_infos, align_info) = out
                proteins.remove_annotation(prot_id, pdb_id, chain)
                if not skip_inserts:
                    sus_complexes.add(pdb_id)
                    sus_structures.add((pdb_id, chain))
                if config.verbosity >= 6:
                    seq_id, alignment_text = align_info
                    print('Alignment got filtered:', prot_id, pdb_id, chain, len(sub_infos), seq_id, alignment_text)
                t_i_1 += time.time()
                continue
            t_i_1 += time.time()
            t_i_2 += time.time()

            (prot_id, pdb_id, chain, alignment, seq_id, coverage, interaction_partners, chain_type_map,
             oligo, sub_infos, backmap, atom_count, last_residue, first_residue, chainlist, rare_residues) = out

            if config.verbosity >= 6:
                print(f'Receiving alignment results from {prot_id} {pdb_id} {chain}: {seq_id} {coverage}')

            if record_package:
                max_runtime_package.append((prot_id, pdb_id, chain, seq_id, coverage, len(alignment)))

            if record_min_package:
                min_runtime_package.append((prot_id, pdb_id, chain, seq_id, coverage, len(alignment)))

            proteins.set_coverage(prot_id, pdb_id, chain, coverage)
            proteins.set_sequence_id(prot_id, pdb_id, chain, seq_id)
            proteins.set_sub_infos(prot_id, pdb_id, chain, sub_infos)
            proteins.set_atom_count(pdb_id, atom_count)

            proteins.set_backmap(prot_id, pdb_id, chain, backmap)

            if indel_analysis_follow_up:
                proteins.set_alignment(prot_id, pdb_id, chain, alignment)

            t_i_3 += time.time()

            proteins.set_interaction_partners(pdb_id, interaction_partners)
            proteins.set_chain_type_map(pdb_id, chain_type_map, chainlist)

            t_i_4 += time.time()

            proteins.set_oligo(pdb_id, chain, oligo)
            proteins.set_last_residue(pdb_id, chain, last_residue)
            proteins.set_first_residue(pdb_id, chain, first_residue)

            t_i_5 += time.time()

            if (pdb_id, chain) not in safe_structures:
                structure_insertion_list.add((pdb_id, chain))

            t_i_6 += time.time()

            safe_complexes.add(pdb_id)
            safe_structures.add((pdb_id, chain))

            t_i_7 += time.time()

            prot_db_id = proteins.get_protein_db_id(prot_id)
            alignment_insertion_list.append((prot_id, prot_db_id, pdb_id, chain, alignment))

            t_i_8 += time.time()

            config.rare_residues.update(rare_residues)

            t_i_9 += time.time()

        started_processes -= 1

        record_package = False

    if config.verbosity >= 5:
        print('Times for data integration:', (t_i_1 - t_i_0), (t_i_3 - t_i_2), (t_i_4 - t_i_3), (t_i_5 - t_i_4), (t_i_6 - t_i_5), (t_i_7 - t_i_6), (t_i_8 - t_i_7), (t_i_9 - t_i_8))

    return (max_runtime, warn_map, alignment_insertion_list, structure_insertion_list, package_runtimes, max_runtime_package,
            min_runtime_package, record_package, safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes)


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
    sus_structures = set()

    safe_structures = set()
    safe_complexes = set()

    task_list = []

    for prot_id in prot_ids:
        if not proteins.is_protein_stored(prot_id):
            continue
        if prot_id not in proteins.indels:
            if proteins.is_completely_stored(prot_id) and not get_all_alignments:
                continue

        annotation_list = proteins.get_protein_annotation_list(prot_id)

        for (pdb_id, chain) in annotation_list:
            if not proteins.is_annotation_stored(pdb_id, chain, prot_id):
                continue

            structure_id = proteins.get_structure_db_id(pdb_id, chain)
            if indel_analysis_follow_up:
                get_alignment_out = proteins.get_alignment(prot_id, pdb_id, chain)
                if isinstance(get_alignment_out, str):
                    config.errorlog.add_error('Get alignment error: %s. %s %s %s %s' % (get_alignment_out, prot_id, pdb_id, chain, str(get_all_alignments)))
                    continue
                else:
                    target_seq, template_seq = get_alignment_out
            else:
                target_seq, template_seq = proteins.pop_alignment(prot_id, pdb_id, chain)
            aaclist = proteins.getAACList(prot_id)
            prot_db_id = proteins.get_protein_db_id(prot_id)

            task_list.append((prot_id, pdb_id, chain, structure_id, target_seq, template_seq, aaclist, prot_db_id))
            safe_structures.add((pdb_id, chain))
            safe_complexes.add(pdb_id)

    t11 = time.time()
    if config.verbosity >= 2:
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
                
    gc.collect()

    t2 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 2: %s" % (str(t2 - t11)))

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

    for prot_id in prot_ids:
        
        annotation_list = proteins.get_protein_annotation_list(prot_id)

        seq = proteins.get_sequence(prot_id)
        aaclist = proteins.getAACList(prot_id)
        structure_infos = []
        for (pdb_id, chain) in annotation_list:
            if proteins.is_annotation_stored(pdb_id, chain, prot_id):
                continue
            if prot_id not in cost_map:
                cost_map[prot_id] = 0
            str_len = proteins.structures[(pdb_id, chain)].get_seq_len()
            if str_len is not None:
                cost_map[prot_id] += len(seq)*str_len
            else:
                cost_map[prot_id] += (len(seq))**2
            oligo = proteins.get_oligo(pdb_id, chain)
            structure_infos.append((pdb_id, chain, oligo))

        if prot_id not in cost_map:
            continue

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

    prots_todo = set(protein_packages.keys())
    started_processes = 0

    while started_processes < config.proc_n and len(prots_todo) > 0:
        alignment_results, started_processes, protein_packages, cost_map, chunk, chunk_cost, prots_todo = package_alignment_processes(prots_todo, cost_map, optimal_chunk_cost, proteins, protein_packages, alignment_results, started_processes, chunk, chunk_cost, config, mapping_dump)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Alignment Part 4: %s" % (str(t4 - t3)))


    alignment_insertion_list = []
    structure_insertion_list = set()
    warn_map = set()

    database_structure_list = None

    package_runtimes = []
    max_runtime = 0
    min_runtime = sys.maxsize
    max_runtime_package = None
    min_runtime_package = []
    record_package = False

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

                (max_runtime, warn_map, alignment_insertion_list, structure_insertion_list, package_runtimes,
                    max_runtime_package, min_runtime_package, record_package,
                    safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes) = process_align_outs(proteins, 
                            config, align_outs, max_runtime, warn_map, alignment_insertion_list, structure_insertion_list, package_runtimes, max_runtime_package, min_runtime_package, record_package,
                            indel_analysis_follow_up, safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes, skip_inserts)

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
            (max_runtime, warn_map, alignment_insertion_list, structure_insertion_list, package_runtimes,
                max_runtime_package, min_runtime_package, record_package,
                safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes) = process_align_outs(proteins,
                    config, align_outs, max_runtime, warn_map, alignment_insertion_list, structure_insertion_list, package_runtimes, max_runtime_package, min_runtime_package, record_package,
                    indel_analysis_follow_up, safe_complexes, safe_structures, sus_complexes, sus_structures, started_processes, skip_inserts)
            
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

    if config.verbosity >= 5 and max_runtime_package is not None:
        if len(package_runtimes) > 0:
            print('Maximal package runtime:', max(package_runtimes))
            print('Minimal package runtime:', min(package_runtimes))
        if len(package_runtimes) > 0:
            print('Packages mean runtime:',  statistics.mean(package_runtimes))
            print('Packages runtime deviation:', statistics.stdev(package_runtimes))
        print('Max runtime package length:', len(max_runtime_package))
        if len(max_runtime_package) > 0:
            print('Max runtime package avg seq_id:', statistics.mean([x[3] for x in max_runtime_package]))
            print('Max runtime package avg coverage:', statistics.mean([x[4] for x in max_runtime_package]))
            print('Max runtime package total alignment length:', sum([x[5] for x in max_runtime_package]))
        print('All packages length:', len(min_runtime_package))
        if len(min_runtime_package) > 0:
            print('All packages avg seq_id:', statistics.mean([x[3] for x in min_runtime_package]))
            print('All packages avg coverage:', statistics.mean([x[4] for x in min_runtime_package]))
            print('All packages total alignment length:', sum([x[5] for x in min_runtime_package]))

    if config.verbosity >= 2:
        t5 = time.time()
        print("Alignment Part 5: %s" % (str(t5 - t4)))
        if len(config.rare_residues) > 0:
            print('######### Detected some rare residues:', config.rare_residues)

    if not skip_inserts:
        database_structure_list = database.insertStructures(structure_insertion_list,
                                                            proteins, config, results=database_structure_list,
                                                            return_results=config.low_mem_system)

    if config.verbosity >= 2:
        t6 = time.time()
        print("Alignment Part 6: %s" % (str(t6 - t5)))

    if not skip_inserts:
        if config.verbosity >= 2:
            t7 = time.time()
            print("Alignment Part 7: %s" % (str(t7 - t6)))

        database.insertAlignments(alignment_insertion_list, proteins, config)

        if config.verbosity >= 2:
            t8 = time.time()
            print("Alignment Part 8: %s" % (str(t8 - t7)))

    # Due the removal of annotations in the previous loop, we might to remove some structures and complexes
    structures_to_remove = sus_structures - safe_structures
    proteins.remove_structures(structures_to_remove)
    complexes_to_remove = sus_complexes - safe_complexes
    if config.verbosity >= 5:
        print('Len of sus_complexes:', len(sus_complexes), 'Len of safe complexes:', len(safe_complexes), 'Len of complexes_to_remove:', len(complexes_to_remove))
        if len(complexes_to_remove) < 50:
            print('Remove complexes:', complexes_to_remove)
    proteins.remove_complexes(complexes_to_remove)

    config.removed_structures = structures_to_remove

    del mapping_dump

@ray.remote(max_calls = 1)
def paraMap_wrapper(mapping_dump, tasks, static_model_path = None):
    ray_hack()

    config = mapping_dump
    tasks = unpack(tasks)
    return paraMap(config, tasks, static_model_path = static_model_path)

def paraMap(config, tasks, static_model_path = None):
    pdb_path = config.pdb_path

    results = []

    for u_ac, pdb_id, chain, structure_id, target_seq, template_seq, aaclist, prot_id in tasks:
        try:
            if config.model_db_active:
                if is_alphafold_model(pdb_id):
                    model_path = alphafold_model_id_to_file_path(pdb_id, config)
                    #is_model = True
                else:
                    model_path = static_model_path
            else:
                model_path = static_model_path

            #template_page, atom_count = pdbParser.standardParsePDB(pdb_id, pdb_path, obsolete_check=True, model_path=model_path)
            template_page, interaction_partners, chain_type_map, oligo, atom_count, chain_list, rare_residues = pdbParser.getStandardizedPdbFile(pdb_id, pdb_path, verbosity=config.verbosity, obsolete_check=True, model_path=model_path)

            seq_res_map, last_residue, first_residue = globalAlignment.createTemplateFasta(template_page, pdb_id, chain, config, onlySeqResMap=True)

            sub_out = globalAlignment.getSubPos(config, u_ac, target_seq, template_seq, aaclist, seq_res_map)

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
    ray_hack()

    config = align_dump

    return pack(align(config, package, static_model_path=model_path, remote = True))


def align(config, package, static_model_path=None, remote=False):
    try:

        t0 = time.time()

        results = []
        aligner_class = globalAlignment.init_bp_aligner_class()

        for prot_specific_mapping_dump, structure_infos in package:

            if remote:
                (u_ac, seq, aaclist) = ray.get(prot_specific_mapping_dump, timeout = 600)
            else:
                (u_ac, seq, aaclist) = prot_specific_mapping_dump

            if config.verbosity >= 5:
                print(f'Alignment of {u_ac} to {structure_infos}')

            for pdb_id, chain, oligo in structure_infos:

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

                    parse_out = pdbParser.getStandardizedPdbFile(pdb_id, pdb_path, oligo=oligo, verbosity=config.verbosity, model_path=model_path)

                    if parse_out is None:
                        results.append((u_ac, pdb_id, chain, 'pdbParser failed for %s' % str(pdb_id)))
                        continue

                    if isinstance(parse_out, str):
                        results.append((u_ac, pdb_id, chain, f'pdbParser failed for {pdb_id}, error message: {parse_out}'))
                        continue

                    (template_page, interaction_partners, chain_type_map, oligo, atom_count, chainlist, rare_residues) = parse_out

                    align_out = globalAlignment.alignBioPython(config, u_ac, seq, pdb_id, template_page, chain, aaclist, rare_residues=rare_residues, aligner_class = aligner_class)
                    
                    if isinstance(align_out, str):
                        align_out = '%s %s' % (align_out, str(model_path))
                        results.append((u_ac, pdb_id, chain, align_out))
                        continue

                    (coverage, seq_id, sub_infos, backmap, alignment_pir, times, aaclist, last_residue, first_residue) = align_out

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

        t1 = time.time()

        runtime = t1-t0
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        results = [('', '', '', f'\nUnknown critical error in alignment: \n{e}\n{str(f)}\n{g}\n')]
        runtime = 0

    return (results, runtime)


def debug_freeze(config, stop_number = 1):
    print('\n\n', config.ray_init_counter, '\n\n')
    if config.ray_init_counter == stop_number:
        time.sleep(999999999)
    return

def print_locals(items):
    for name, size in sorted(((name, sys.getsizeof(value)) for name, value in items), key=lambda x: -x[1])[:10]:
        print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))

def core(protein_list, genes, indels, multi_mutation_objects, config, session, outfolder, session_name, out_objects):
    if len(protein_list) == 0:
        return None, 0

    if config.proc_n > 1:
        ray_init(config)
        config.ray_init_counter += 1    
        #monitor_ray_store()

    background_process_MS = None

    # transform the protein map into a Proteins object
    proteins = protein_package.Proteins(protein_list, indels, multi_mutation_objects)


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
                if config.verbosity >= 3:
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

        database.indelCheck(proteins, session, config)

        database.insertMultiMutations(proteins, genes, session, config)

        if config.verbosity >= 6:
            print('Proteins state after positionCheck:')
            proteins.print_protein_state()

        t3 = time.time()
        if config.verbosity >= 2:
            print("Time for positionCheck: %s" % (str(t3 - t2)))

        if config.verbosity >= 1:
            print("Before getSequences")

        background_iu_process = getSequences(proteins, config)

        t4 = time.time()
        if config.verbosity >= 2:
            print("Time for getSequences: %s" % (str(t4 - t3)))

        if config.verbosity >= 1:
            print("Before autoTemplateSelection")
        autoTemplateSelection(config, proteins)

        t5 = time.time()
        if config.verbosity >= 2:
            print("Time for Template Selection: %s" % (str(t5 - t4)))

        indel_analysis_follow_up = len(proteins.indels) > 0 and not config.skip_indel_analysis
        
        if config.verbosity >= 1:
            print("Before paraAlignment")
        paraAlignment(config, proteins, indel_analysis_follow_up=indel_analysis_follow_up)

        t6 = time.time()
        if config.verbosity >= 2:
            print("Time for Alignment: %s" % (str(t6 - t5)))

        if background_iu_process is not None:
            background_iu_process.join()  # Disorder values have to finished before the classification happens
            background_iu_process.close()
            background_iu_process = None

        if config.verbosity >= 1:
            print("Before paraAnnotate")

        amount_of_structures = templateFiltering.paraAnnotate(config, proteins, indel_analysis_follow_up=indel_analysis_follow_up)
        
        t7 = time.time()
        if config.verbosity >= 2:
            print("Time for Annotation: %s" % (str(t7 - t6)))

        if not config.skip_indel_analysis:
            indel_analysis.para_indel_analysis(proteins, config)

        t1 = time.time()
        if config.verbosity >= 2:
            print("Time for Indelanalysis: %s" % (str(t1 - t7)))
                
        proteins.deconstruct()
        del proteins

        # join the background inserts
        if background_process_MS is not None:
            background_process_MS.join()
            background_process_MS.close()
            background_process_MS = None

        t2 = time.time()
        if config.verbosity >= 2:
            print('Resttime for background inserts: ', t2 - t1)

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

    if config.ray_init_counter >= 1000:
        ray.shutdown()
        #time.sleep(10)
        config.ray_init_counter = 0
        pass
    return out_objects, amount_of_structures		


# @profile
def main(filename, config, custom_db_mode=False):
    # main function takes another argument custom db, default is False
    mrna_fasta = config.mrna_fasta

    session = 0  # This can later be used, structman_main.py could give a specific session id and the pipeline can then expand that session

    if config.verbosity >= 2:
        print('==== Start main pipeline ====')

    config.errorlog.start(filename, config.outfolder)

    t0 = time.time()

    if config.iupred_path != '':
        sys.path.append(f'{os.path.abspath(os.path.realpath(config.iupred_path))}')
        #os.environ["PYTHONPATH"] = f'{os.path.abspath(os.path.realpath(config.iupred_path))}:{os.environ.get("PYTHONPATH", "")}'
        try:
            import iupred3_lib
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            config.errorlog.add_error(f'IUpred path was given, but import failed: {config.iupred_path}\n{e}\n{f}\n{g}')
            config.iupred_path = ''

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
            anno_db = "%s_annovar" % db_name.rsplit("_", 1)[0]  # BUG: undefined variable
            if config.verbosity >= 1:
                print('Convert vcf file format using Annovar')
            if mrna_fasta is not None:
                '... and using mrna file: ', mrna_fasta
            nfname = annovar.annovar_pipeline(filename, config.tax_id, config.annovar_path, config.db_address, config.db_user_name, config.db_password, anno_db, mrna_fasta, ref_id=config.ref_genome_id)
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

    out_objects = None

    total_amount_of_analyzed_structures = 0
    config.ray_init_counter = 0
    for nr_temp_file, temp_infile in enumerate(temp_infiles):
        if temp_infile is not None:
            if config.verbosity >= 1:
                print('Infile splitting due to low memory system, processing infile split nr.:', nr_temp_file + 1, 'out of', len(temp_infiles))
            proteins_chunks, _, genes = buildQueue(config, temp_infile, already_split=True)
            os.remove(temp_infile)

        #The gene objects need to exist as a whole outside of the chunked core loop. Isoforms that correspond to the same gene can be chunked into different bins.
        database.geneCheck(genes, config)

        chunk_nr = 1
        for protein_list, indels in proteins_chunks:

            config.indels_given_by_input = (len(indels) > 0)

            if config.verbosity >= 1:
                print("Chunk %s/%s" % (str(chunk_nr), str(len(proteins_chunks))))
            chunk_nr += 1

            protein_list, indels, multi_mutation_objects, _ = sequenceScan(config, protein_list, indels)

            out_objects, amount_of_structures = core(protein_list, genes, indels, multi_mutation_objects, config, session, config.outfolder, session_name, out_objects)

            total_amount_of_analyzed_structures += amount_of_structures

            if config.verbosity >= 3:
                for name, size in sorted(((name, sys.getsizeof(value)) for name, value in globals().items()), key=lambda x: -x[1])[:10]:
                    print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
                print(list_fds())

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

    config.errorlog.stop()

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

    return session, config
