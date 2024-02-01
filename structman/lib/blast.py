import os
import subprocess
import xml.etree.ElementTree as ET


# called by serializedPipeline
def blast(seq, name, blast_path, blast_db_path, nr=0, seq_thresh=0.35, cov_thresh=0.5, cwd=None):
    FNULL = open(os.devnull, 'w')
    if cwd is None:
        cwd = os.getcwd()
    seq = seq.replace(" ", "")
    target_length = len(seq)

    blast_in = "%s/blast_in.tmp.fasta" % cwd
    blast_out = "%s/blast_out.tmp" % cwd

    f = open(blast_in, "wb")
    f.write(">" + name + "\n" + seq)
    f.close()

    path = blast_path + "blastp"
    arg1 = "-db"
    arg2 = blast_db_path
    arg3 = "-evalue"
    arg4 = "1e-10"
    arg5 = "-outfmt"
    arg6 = "5"
    arg7 = "-query"
    arg8 = blast_in
    arg9 = "-out"
    arg10 = blast_out
    arg11 = "-num_alignments"
    arg12 = str(nr)
    if nr == 0:
        p = subprocess.Popen([path, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10])  # ,stdout=FNULL)
    else:
        p = subprocess.Popen([path, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12], stdout=FNULL)
    p.wait()
    f = open(blast_out, "r")
    lines = f.read()
    f.close()
    os.remove(blast_in)
    os.remove(blast_out)

    entries = {}

    if len(lines) == 0:
        return (entries)

    root = ET.fromstring(lines)
    blast_iter_hits = root[8][0][4]

    for child in blast_iter_hits:
        hit_defs = child[2].text.split("|")
        new_hit = False
        new_pdb = ''
        hits = {}
        for de in hit_defs:
            if de == "pdb":
                new_hit = True
            elif new_hit:
                new_pdb = de
                new_hit = False
            elif new_pdb != '':
                chain = de[0]
                if chain != ' ':
                    # In the blast-xml-output, if the number of chains exceeds the alphabet,
                    # the gi-identifier fails to support the correct chain id
                    if de[1] != " ":
                        words = de.split(" ")
                        found_true_chain_id = False
                        for w in words:
                            if w == 'Chain':
                                found_true_chain_id = True
                            elif found_true_chain_id:
                                chain = w[0]
                                break
                    if new_pdb not in hits:
                        hits[new_pdb] = [chain, set([chain])]
                    else:
                        hits[new_pdb][1].add(chain)
                new_pdb = ''
                chain = ''

        hit_hsp = child[5][0]

        aln_length = hit_hsp[13].text
        seq_id = 100 * float(hit_hsp[10].text) / float(aln_length)
        coverage = float(aln_length) / float(target_length)
        for hit in hits:
            pdb_id = hit
            chain = hits[hit][0]
            oligos = hits[hit][1]
            if not len(chain) > 1:
                if not (pdb_id, chain) in entries:
                    entries[(pdb_id, chain)] = {"Seq_Id": seq_id, "Coverage": coverage, "Oligo": oligos, "Length": aln_length}
                else:
                    if coverage > entries[(pdb_id, chain)]["Coverage"]:
                        entries[pdb_id] = {"Seq_Id": seq_id, "Coverage": coverage, "Oligo": oligos, "Length": aln_length}
                    entries[(pdb_id, chain)]["Oligo"].update(oligos)
    structures = {}

    for (pdb_id, chain) in entries:
        #print("\nstep1\n")
        if not (float(entries[(pdb_id, chain)]["Seq_Id"]) < seq_thresh):
            #print("\nstep2\n")
            if not ((entries[(pdb_id, chain)]["Length"] < 50) and (entries[(pdb_id, chain)]["Coverage"] < cov_thresh)):
                #print("\nstep3\n")
                structures[(pdb_id, chain)] = entries[(pdb_id, chain)]

    return structures
