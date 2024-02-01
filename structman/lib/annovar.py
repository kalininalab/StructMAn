import ftplib
import gzip
import os
import shutil
import subprocess
import sys

import pymysql as MySQLdb

from structman.lib.sdsc.consts import codons as codon_consts


# called by serializedPipeline
def annovar_pipeline(vcf_file, tax_id, annovar_path, host, user, pw, db_name, mrna, mail_address='', ref_id=None):
    if mrna is None:
        print(f"Start Annovar Pipeline with vcf_file: {vcf_file}, taxonomic ID: {str(tax_id)} and ref_id: {str(ref_id)}")
        print(annovar_path)
        (name, ucsc_id, ref_path) = searchInAnnovarDB(tax_id, host, user, pw, db_name, ref_id=ref_id)

        if name == "No Entry":
            (ucsc_id, name) = taxIdToUcscId(tax_id)
            name = name.lower()
            if ucsc_id == "not in ucsc":
                print("Did not found the species in the UCSC, search in NCBI-refSeq and construct a new reference ...")
                (ref_gene_path, ref_seq_path, name, ucsc_id) = buildRefFromNCBI(tax_id, annovar_path, mail_address)
            else:
                print("Found the species in the UCSC, constrcuting a new reference...")
                (ref_gene_path, ref_seq_path) = downloadRefFromUCSC(name, ucsc_id, annovar_path)
            ref_path = createRef(ref_gene_path, ref_seq_path, name, tax_id, ucsc_id, annovar_path, host, user, pw, db_name)
        else:
            print(f"Found the species in the local reference database: {ref_path}")
        smlf_file = vcfToSmlf(vcf_file, ref_path, ucsc_id, annovar_path)
        return smlf_file
    else:
        return ToSmlf(vcf_file, mrna)


def ToSmlf(vcf_file, fasta):
    f = open(fasta, 'r')
    lines = f.readlines()
    f.close()

    nucleic = True
    n_acids = set(['A', 'C', 'T', 'G'])

    gene_seq_map = {}
    for line in lines:
        line = line.replace('\n', '')
        if line[0] == '>':
            g_id = line.split()[0][1:]
            gene_name = line.split()[1]
            gene_seq_map[g_id] = ''
        else:
            gene_seq_map[g_id] += line
            """
            if nucleic:
                for char in line:
                    if not char in n_acids:
                        nucleic = False
            """

    stop_codons = set(codon_consts.STOP_CODONS)
    rev_stop_codons = set(codon_consts.REV_STOP_CODONS)

    inverse = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    inverse_strands = set()

    for g in gene_seq_map:
        seq = gene_seq_map[g]

        if seq[0:3] in rev_stop_codons and seq[-3:] not in stop_codons:
            seq = seq[::-1]
            inv_seq = ''

            for nuc in list(seq):
                inv_seq = f'{inv_seq}{inverse[nuc]}'
            gene_seq_map[g] = inv_seq
            inverse_strands.add(g)

    f = open(vcf_file, 'r')
    lines = f.readlines()
    f.close()

    smlf_lines = []

    for line in lines:
        if line[0] == '#':
            continue
        words = line.replace('\n', '').split()
        g_id = words[0]
        if g_id == '.':
            continue
        pos = int(words[1])
        pid = words[2]
        ref = words[3]
        alt = words[4]
        if ref == '-' or alt == '-':
            continue
        if len(ref) == 1 and len(alt) == 1:
            seq = gene_seq_map[g_id]

            inv_strand = g_id in inverse_strands

            if inv_strand:
                pos = (len(seq) - pos) + 1
                ref = inverse[ref]
                alt = inverse[alt]
            if ref != seq[pos - 1]:
                print(pos)
                print(seq)
                print(g_id)

            # if nucleic:
            if pos % 3 == 0:
                aa_pos = pos / 3
                triple = seq[pos - 3:pos]
                new_triple = f"{triple[0]}{triple[1]}{alt}"
            elif pos % 3 == 1:
                aa_pos = pos / 3 + 1
                triple = seq[pos - 1:pos + 2]
                new_triple = f"{alt}{triple[1]}{triple[2]}"
            elif pos % 3 == 2:
                aa_pos = pos / 3 + 1
                triple = seq[pos - 2:pos + 1]
                new_triple = f"{triple[0]}{alt}{triple[2]}"

            if new_triple in stop_codons:
                # nonsense mutation (Stop codon)
                continue

            if triple in stop_codons:
                continue

            old_aa = codon_consts.CODONS[triple]
            new_aa = codon_consts.CODONS[new_triple]

            if old_aa != new_aa:
                s_line = f"{g_id}\t{old_aa}{aa_pos}{new_aa}"
                smlf_lines.append(s_line)

    cwd = os.getcwd()
    if vcf_file.count("/") > 0:
        out_path = vcf_file.rsplit("/", 1)[0]
        name = vcf_file.rsplit("/", 1)[1].rsplit(".", 1)[0]
    else:
        out_path = cwd
        name = vcf_file.rsplit(".", 1)[0]
    filename = f"{out_path}/{name}.smlf"

    f = open(filename, 'w')
    f.write('\n'.join(smlf_lines))
    f.close()

    return filename


def buildRefFromNCBI(tax_id, annovar_path, mail_address):
    summary_path = f"{annovar_path}/Taxonomy/assembly_summary_refseq.txt"
    f = open(summary_path, "r")
    lines = f.readlines()
    f.close()
    for line in lines:
        words = line.split("\t")
        if words[5].replace("\n", "") == str(tax_id):
            ftp_path = words[19].replace("\n", "")
            reference_tag = words[4].replace("\n", "")
            version_status = words[11].replace("\n", "")
            sc_name = words[7].replace("\n", "")
            if version_status == "latest":
                break

    name = "_".join(sc_name.split())
    name = name.lower()
    ucsc_id = f"{name}_{str(tax_id)}"
    refdir = f"{annovar_path}/{name}db"
    if not os.path.isdir(refdir):
        os.mkdir(refdir)

    gfc_id = ftp_path.rsplit("/", 1)[1]

    # Connect to FTP
    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")

    ftp.login(user="anonymous", passwd=mail_address)
    ftp.set_pasv(True)
    ftp.cwd(ftp_path.replace("ftp://ftp.ncbi.nlm.nih.gov", ""))

    gff_gz_path = f"{refdir}/{ucsc_id}.gff.gz"
    f = open(gff_gz_path, 'wb')
    ftp.retrbinary(f'RETR {gfc_id}_genomic.gff.gz', f.write, 1024)
    f.close()

    if not os.path.isdir(f"{refdir}/{ucsc_id}_seq"):
        os.mkdir(f"{refdir}/{ucsc_id}_seq")
    fna_gz_path = f"{refdir}/{ucsc_id}_seq/{ucsc_id}.fna.gz"
    f = open(fna_gz_path, 'wb')
    ftp.retrbinary('RETR {gfc_id}_genomic.fna.gz', f.write)
    f.close()

    ftp.quit()

    # dezip files
    f = gzip.open(gff_gz_path, 'rb')
    page = f.read()
    f.close()

    gff_path = f"{refdir}/{ucsc_id}.gff"
    f = open(gff_path, "wb")
    f.write(page)
    f.close()

    os.remove(gff_gz_path)

    f = gzip.open(fna_gz_path, 'rb')
    page = f.read()
    f.close()

    fna_path = f"{refdir}/{ucsc_id}_seq/{ucsc_id}.fna"
    f = open(fna_path, "wb")
    f.write(page)
    f.close()

    os.remove(fna_gz_path)

    # use gff3toGenePred-Tool
    ref_gene_path = "%s/%s_refGene.txt" % (refdir, ucsc_id)

    subprocess.call(["%s/gff3togenepred/gff3ToGenePred" % annovar_path, gff_path, ref_gene_path])

    os.remove(gff_path)

    ref_seq_path = "%s/%s_seq" % (refdir, ucsc_id)

    return (ref_gene_path, ref_seq_path, name, ucsc_id)


def downloadRefFromUCSC(name, ucsc_id, annovar_path):
    if not os.path.isdir("%s/%sdb" % (annovar_path, name)):
        os.mkdir("%s/%sdb" % (annovar_path, name))
    name = name.lower()
    subprocess.call(["perl", "%s/annotate_variation.pl" % annovar_path, "--downdb", "refGene", "--buildver", ucsc_id, "%s/%sdb" % (annovar_path, name)])
    ref_gene_path = "%s/%sdb/%s_refGene.txt" % (annovar_path, name, ucsc_id)

    subprocess.call(["perl", "%s/annotate_variation.pl" % annovar_path, "--downdb", "seq", "--buildver", ucsc_id, "%s/%sdb/%s_seq" % (annovar_path, name, ucsc_id)])
    ref_seq_path = "%s/%sdb/%s_seq" % (annovar_path, name, ucsc_id)

    return (ref_gene_path, ref_seq_path)


def createRef(ref_gene_path, ref_seq_path, name, tax_id, ucsc_id, annovar_path, host, user, pw, db_name):
    seq_files = os.listdir(ref_seq_path)
    if "annovar_downdb.log" in seq_files:
        seq_files.remove("annovar_downdb.log")
    if len(seq_files) == 1:
        seqtype = "-seqfile"
        seq_path = "%s/%s" % (ref_seq_path, seq_files[0])
    elif len(seq_files) > 1:
        seqtype = "-seqdir"
        seq_path = ref_seq_path
    else:
        return "ERROR"

    folder_path = ref_gene_path.rsplit("/", 1)[0]
    out_path = "%s/%s_refGeneMrna.fa" % (folder_path, ucsc_id)
    subprocess.call(["perl", "%s/retrieve_seq_from_fasta.pl" % annovar_path, ref_gene_path, seqtype, seq_path, "-format", "refGene", "-out", out_path])

    shutil.rmtree(ref_seq_path)

    db = MySQLdb.connect(host, user, pw, db_name)
    cursor = db.cursor()

    sql = """INSERT INTO Organism(Name,Ucsc_Id,Mrna_Path,Taxon_Id)
            Values('%s','%s','%s','%s')""" % (name, ucsc_id, "/".join(out_path.rsplit("/", 2)[1:]), str(tax_id))
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in searchInAnnovarDB: %s,%f" % (sql, f))
    db.close()

    return out_path.rsplit("/", 1)[0]


def vcfToSmlf(vcf_file, ref_path, ucsc_id, annovar_path):
    FNULL = open(os.devnull, 'w')
    cwd = os.getcwd()
    if vcf_file.count("/") > 0:
        out_path = vcf_file.rsplit("/", 1)[0]
        name = vcf_file.rsplit("/", 1)[1].rsplit(".", 1)[0]
    else:
        out_path = cwd
        name = vcf_file.rsplit(".", 1)[0]

    f = open(vcf_file, "r")
    line = f.readline()
    f.close()
    info = line.split("=")
    if not len(info) == 2:
        raise NameError("first line of vcf-file must contain version-info")
    if not info[0] == "##fileformat":
        raise NameError("first line of vcf-file must contain version-info")

    vcf_version = info[1][4:]

    if int(vcf_version.split(".")[0]) > 3:
        subprocess.call(["perl", "%s/convert2annovar.pl" % annovar_path, "-format", "vcf4", "-allsample", "-withfreq", vcf_file, "-outfile", "%s/%s.avinput" % (out_path, name)])
    else:
        raise NameError("vcf-version is too old")

    # if os.path.isfile(vcf_file):
    #    os.remove(vcf_file)

    subprocess.call(["perl", "%s/annotate_variation.pl" % annovar_path, "-geneanno", "-buildver", ucsc_id, "%s/%s.avinput" %
                     (out_path, name), "%s/%s" % (annovar_path, ref_path), "-outfile", "%s/%s" % (cwd, name)])

    # if os.path.isfile("%s/%s.avinput" % (out_path, name)):
    #     os.remove("%s/%s.avinput" % (out_path, name))

    # if os.path.isfile("%s/%s.variant_function" % (cwd, name)):
    #     os.remove("%s/%s.variant_function" % (cwd, name))

    # if os.path.isfile("%s/%s.invalid_input" % (cwd, name)):
    #     os.remove("%s/%s.invalid_input" % (cwd, name))

    if os.path.isfile("%s/%s.exonic_variant_function" % (cwd, name)):
        f = open("%s/%s.exonic_variant_function" % (cwd, name), "r")
        lines = f.readlines()
        f.close()
    else:
        return ""

    outlines = []
    for line in lines:
        words = line.split("\t")
        if len(words) < 2:
            print(line)
        elif words[1] == "nonsynonymous SNV":
            muts = words[2].split(",")
            chromo = words[3]
            chr_pos = words[4]
            g_coord = '%s:%s' % (chromo, chr_pos)
            for mut in muts:
                infos = mut.split(":")
                if len(infos) > 3:
                    refseq = infos[1]
                    aac = infos[4].replace("p.", "")
                    outlines.append("%s\t%s\t\t%s" % (refseq, aac, g_coord))
    os.remove("%s/%s.exonic_variant_function" % (cwd, name))

    filename = "%s/%s.smlf" % (out_path, name)

    f = open(filename, "wb")
    f.write("\n".join(outlines))
    f.close()

    return filename


def searchInAnnovarDB(tax_id, host, user, pw, db_name, ref_id=None):
    db = MySQLdb.connect(host, user, pw, db_name)
    cursor = db.cursor()

    sql = "SELECT Name,Ucsc_Id,Mrna_Path FROM Organism WHERE Taxon_Id = '%s'" % str(tax_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in searchInAnnovarDB: %s,%f" % (sql, f))

    if results == ():
        return ("No Entry", "", "")

    if len(results) == 1:
        (name, ucsc_id, ref_path) = results[0]
    elif ref_id is None:
        raise NameError("Error in searchInAnnovarDB: %s" % str(ref_id))
    else:
        for row in results:
            if ref_id == row[1]:
                (name, ucsc_id, ref_path) = row

    return (name, ucsc_id, ref_path.split("/")[0])


def taxIdToUcscId(tax_id):
    host = "genome-mysql.cse.ucsc.edu"
    user = "genomep"
    pw = "password"
    db_name = "hgcentral"
    db = MySQLdb.connect(host, user, pw, db_name)
    cursor = db.cursor()

    sql = "SELECT name,active,organism FROM dbDb WHERE taxId = '%s'" % str(tax_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in taxIdToUcscId: %s,%f" % (sql, f))
    db.close()

    if results == ():
        return ("not in ucsc", "")

    ids = []
    for row in results:
        name = row[0]
        active = int(row[1])
        if active == 1:
            ids.append(name)
        orga = row[2]

    # Find the common substring in the id and the id-number
    cut_ids = ids
    found_new_common_char = True
    common = ""
    while found_new_common_char:
        c = cut_ids[0][0]
        for i in cut_ids:
            if not c == i[0]:
                found_new_common_char = False
        if found_new_common_char:
            common = "%s%s" % (common, c)
            cut_ids = [x[1:] for x in cut_ids]

    try:
        cut_ids = [int(x) for x in cut_ids]
    except:
        raise NameError("cut_ids did not contain only strings\nids: %s" % ",".join(ids))
    max_id_number = max(cut_ids)
    ucscId = "%s%s" % (common, str(max_id_number))
    return (ucscId, orga)
