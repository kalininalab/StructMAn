import os
import random
import shutil
import string
import subprocess
import time
from structman.scripts import createCustomDb


def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


def geneSeqMapToFasta(proteins, outfile, config, filtering_db=None):
    lines = []
    n = 0
    m = 0

    prot_ids = proteins.get_protein_ids()

    n_not_stored = 0

    for prot_id in prot_ids:

        # Previous if-statement:
        # if proteins[prot_id].stored and prot_id not in proteins.indels:
        # Didnt understood why proteins in the indels dict should be forced into the mmseqs search
        # Remove for now, retrieve it, if something comes up, then better comment it
        if proteins[prot_id].stored:
            continue
        n_not_stored += 1

        seq = proteins.get_sequence(prot_id)

        if seq == 0 or seq == 1 or seq == '' or seq == None:
            continue
        if filtering_db is not None:
            folder_key = prot_id.split('-')[0][-2:]
            filename = '%s/%s/%s_ref50_gpw.fasta.gz' % (filtering_db, folder_key, prot_id)
            if os.path.isfile(filename):
                n += 1
                continue

        lines.append('>%s' % prot_id)
        lines.append(seq)
        m += 1

    if n_not_stored == 0:
        return False

    if len(lines) > 0:
        if config.verbosity >= 2:
            print('Filtered ', n, ' Proteins before mmseqs2')
            print(m, ' sequences going into mmseqs2')
        f = open(outfile, 'w')
        f.write('\n'.join(lines))
        f.close()
        return True
    else:
        return 'Empty fasta file'


def parseHits(temp_outfile, option_seq_thresh, small_genes):
    f = open(temp_outfile, 'r')
    lines = f.read().split('\n')
    f.close()

    entries = {}

    if len(lines) == 0:
        return (entries)

    pdb_ids = set()

    for line in lines:
        if line == '':
            continue

        words = line.split()
        gene = words[0]
        hitlist = words[1].split(',')
        seq_id = 100.0 * float(words[2])
        aln_length = int(words[3])
        target_len = int(words[4])
        coverage = float(words[5])

        if aln_length < 50:
            if gene not in small_genes:
                continue
            elif seq_id < (option_seq_thresh * 2):
                continue

        if gene not in entries:
            entries[gene] = {}

        hits = {}

        for hit in hitlist:
            pdb, chain = hit.rsplit('-',1)
            if pdb not in hits:
                hits[pdb] = [chain, set([chain])]
            else:
                hits[pdb][1].add(chain)

        for hit in hits:
            pdb_id = hit
            pdb_ids.add(pdb_id)
            chain = hits[hit][0]
            oligos = hits[hit][1]
            if not len(chain) > 1:
                if not (pdb_id, chain) in entries[gene]:
                    entries[gene][(pdb_id, chain)] = [seq_id, coverage, oligos, aln_length, target_len]
                else:
                    if aln_length > entries[gene][(pdb_id, chain)][3]:
                        entries[gene][pdb_id] = [seq_id, coverage, oligos, aln_length, target_len]
                    entries[gene][(pdb_id, chain)][2].update(oligos)

    return entries, pdb_ids


def apply_mmseqs(mmseqs_tmp_folder, mmseqs2_path, temp_fasta, search_db, gigs_of_ram, errorlog, option_seq_thresh, small_proteins = None, verbosity = 1, number_of_returned_hits = '999999', is_model_db = False):

    if small_proteins is None:
        small_proteins = set()

    temp_outfile = '%s/tmp_outfile_%s.fasta' % (mmseqs_tmp_folder, randomString())

    if verbosity >= 2:
        print(mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder)

    out_format_str = 'query,target,fident,alnlen,tlen,qcov'

    debug_store = False

    if verbosity >= 4:
        if len(small_proteins) == 0:
            cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder, '--format-output', out_format_str,
                    '--max-seqs', number_of_returned_hits, '--min-aln-len', '50', '-s', '7.5', '--split-memory-limit', '%1.0fG' % (0.5 * gigs_of_ram)]
            p = subprocess.Popen(cmds)
        else:
            cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder, '--format-output', out_format_str,
                    '--max-seqs', number_of_returned_hits, '-s', '7.5', '--split-memory-limit', '%1.0fG' % (0.5 * gigs_of_ram)]
            p = subprocess.Popen(cmds)
        p.wait()
    else:
        FNULL = open(os.devnull, 'w')
        if len(small_proteins) == 0:
            cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder, '--format-output', out_format_str,
                    '--max-seqs', number_of_returned_hits, '--min-aln-len', '50', '-s', '7.5', '--split-memory-limit', '%1.0fG' % (0.5 * gigs_of_ram)]
            p = subprocess.Popen(cmds, stdout=FNULL)
        else:
            cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder, '--format-output', out_format_str,
                    '--max-seqs', number_of_returned_hits, '-s', '7.5', '--split-memory-limit', '%1.0fG' % (0.5 * gigs_of_ram)]
            p = subprocess.Popen(cmds, stdout=FNULL)
        p.wait()

        if not os.path.exists(temp_outfile):
            #Something went wrong, run again in not silenced mode
            errorlog.add_error(f'Error in MMSEQS2 call, see command line prints for possible bug log.\nCommand that called mmseqs:\n{cmds}\nDeletion of temp file: {temp_fasta} was suppressed')
            debug_store = True
            p = subprocess.Popen(cmds)


    if verbosity >= 3:
        print(f'MMseqs2 returned: {temp_outfile}')

    if verbosity >= 5:
        f = open(temp_outfile, 'r')
        print(f.read())
        f.close()

    hits, pdb_ids = parseHits(temp_outfile, option_seq_thresh, small_proteins)

    if verbosity >= 3:
        print(f'MMseqs2 results parsed: size of hit map: {len(hits)}')

    os.remove(temp_outfile)

    return hits, pdb_ids, debug_store


# called by serializePipeline
def search(proteins, config, custom_db=False):

    # accepts one more argument, custom db

    mmseqs_tmp_folder = '%s/mmseqs_tmp' % config.mmseqs_tmp_folder
    if not os.path.exists(mmseqs_tmp_folder):
        os.mkdir(mmseqs_tmp_folder)

    mmseqs2_path = config.mmseqs2_path

    # custom db mmseqs search processed here
    if custom_db:
        createCustomDb.main(config)

        mmseqs_custom_db_path = config.custom_db_path + "/custom_db"

        if not os.path.exists(mmseqs_custom_db_path):
            cmds = ["mkdir", "custom_db"]
            p = subprocess.Popen(cmds, cwd=config.custom_db_path)

        cmds = [mmseqs2_path, 'createdb', config.custom_db_fasta_path, "./structman_custom_db"]

        p = subprocess.Popen(cmds, cwd=mmseqs_custom_db_path)
        p.wait()

    # checks if custom db path is provided in config file
    if (config.custom_db_path is not None) and (config.custom_db_path != ''):
        search_dbs = ["%s/%s" % (mmseqs_custom_db_path, "structman_custom_db")]
        numbers_of_returned_hits = [None, '1']
        is_model_dbs = [True]

    elif config.model_db_active:
        search_dbs = [config.mmseqs2_db_path, config.mmseqs2_model_db_path]
        numbers_of_returned_hits = [None, '1']
        is_model_dbs = [False, True]
    else:
        search_dbs = [config.mmseqs2_db_path]
        numbers_of_returned_hits = [None]
        is_model_dbs = [False]
    option_seq_thresh = config.option_seq_thresh

    small_proteins = set()

    u_acs = proteins.get_protein_ids()
    for u_ac in u_acs:
        try:
            if len(proteins.get_sequence(u_ac)) < 100:
                small_proteins.add(u_ac)
        except:
            continue

    t0 = time.time()

    temp_fasta = '%s/tmp_%s.fasta' % (mmseqs_tmp_folder, randomString())
    to_fasta_out = geneSeqMapToFasta(proteins, temp_fasta, config)

    if isinstance(to_fasta_out, str):
        config.errorlog.add_warning('%s , mmseqs2 skipped, %s' % (to_fasta_out, str(list(u_acs)[:10])))
        return None
    if not to_fasta_out:
        # All proteins are stored, no need for mmseqs
        if config.verbosity >= 3:
            print('No need for sequence similarity search, all proteins are stored already.')
        return None

    default_number_of_returned_hits = '999999'

    search_results = []

    t1 = time.time()

    for pos, search_db in enumerate(search_dbs):

        if numbers_of_returned_hits[pos] is None:
            number_of_returned_hits = default_number_of_returned_hits
        else:
            number_of_returned_hits = numbers_of_returned_hits[pos]

        is_model_db = is_model_dbs[pos]

        hits, pdb_ids, debug_store = apply_mmseqs(mmseqs_tmp_folder, mmseqs2_path, temp_fasta, search_db, config.gigs_of_ram, config.errorlog, option_seq_thresh, small_proteins = small_proteins, verbosity = config.verbosity, number_of_returned_hits = number_of_returned_hits, is_model_db = is_model_db)

        for u_ac in small_proteins:
            if u_ac.count(':') == 1:
                if u_ac not in hits:
                    pdb_id, chain = u_ac.split(':')
                    hits[u_ac] = {(pdb_id, chain): [100.0, 1.0, [chain], len(proteins.get_sequence(u_ac)), len(proteins.get_sequence(u_ac))]}
                    pdb_ids.add(pdb_id)

        for fn in os.listdir(mmseqs_tmp_folder):
            subfolder_path = '%s/%s' % (mmseqs_tmp_folder, fn)
            if os.path.exists(subfolder_path):
                if os.path.getmtime(subfolder_path) > config.prog_start_time:
                    try:
                        shutil.rmtree(subfolder_path)
                    except:
                        if config.verbosity >= 4:
                            config.errorlog.add_warning('Tmp folder wipe failed for: %s' % subfolder_path)

        search_results.append((hits, pdb_ids, is_model_db))
    
    if not debug_store:
        os.remove(temp_fasta)

    t2 = time.time()

    if config.verbosity >= 2:
        print("MMseqs2 Part 1: %s" % (str(t1 - t0)))
        print("MMseqs2 Part 2: %s" % (str(t2 - t1)))

    return search_results
