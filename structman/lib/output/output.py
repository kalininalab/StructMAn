# Contains any classes/functions that are used externally
import os
import time

from structman.lib import babel
from structman.lib.database import database

try:
    from structman.lib import modelling
except:
    pass
from structman.lib.output import classification as classification_package



def makeDiffDict(f):
    f = open(f, "r")
    lines = f.readlines()[1:]
    f.close()
    diff_dict = {}
    for line in lines:
        words = line.split("\t")
        diff_dict[words[1]] = [words[0], words[5]]
    return diff_dict


def makeGeneDict(f):
    f = open(f, "r")
    lines = f.readlines()[1:]
    f.close()
    diff_dict = {}
    for line in lines:
        words = line.split("\t")
        diff_dict[words[0]] = words[2]
    return diff_dict


def godiffAna(fileA, fileB, papermode=False):
    go_dict_a = makeDiffDict(fileA)
    go_dict_b = makeDiffDict(fileB)
    result_list = []
    for Id in go_dict_a:
        if Id in go_dict_b:
            diff = float(go_dict_a[Id][1]) - float(go_dict_b[Id][1])
            result_list.append([Id, go_dict_a[Id][0], diff])
    result_list = sorted(result_list, key=lambda x: x[2], reverse=True)

    lines = ["GO-Term\tGO-ID\tScore-Difference"]
    for result in result_list:
        if not papermode:
            lines.append("%s\t%s\t%s" % (result[1], result[0], str(result[2])))
        elif result[1][:2] == ' P':
            lines.append("%s%s\t%s\t%s" % (result[1][3].upper(), result[1][4:], result[0], str(result[2])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    name_b = (fileB.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    outfile = "%s/GO_Diff_%s_%s.tsv" % (fileA.rsplit("/", 1)[0], name_a, name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()


def pathdiffAna(fileA, fileB, papermode=False):
    path_dict_a = makeDiffDict(fileA)
    path_dict_b = makeDiffDict(fileB)
    result_list = []
    for Id in path_dict_a:
        if Id in path_dict_b:
            diff = float(path_dict_a[Id][1]) - float(path_dict_b[Id][1])
            result_list.append([Id, path_dict_a[Id][0], diff])
    result_list = sorted(result_list, key=lambda x: x[2], reverse=True)

    lines = ["Pathway\tReactome-ID\tScore-Difference"]
    for result in result_list:
        lines.append("%s\t%s\t%s" % (result[1], result[0], str(result[2])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    name_b = (fileB.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    outfile = "%s/Path_Diff_%s_%s.tsv" % (fileA.rsplit("/", 1)[0], name_a, name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()


def genediffAna(fileA, fileB):
    gene_dict_a = makeGeneDict(fileA)
    gene_dict_b = makeGeneDict(fileB)
    result_list = []
    for Id in gene_dict_a:
        if Id in gene_dict_b:
            diff = float(gene_dict_a[Id]) - float(gene_dict_b[Id])
            result_list.append([Id, diff])
    result_list = sorted(result_list, key=lambda x: x[1], reverse=True)

    lines = ["Gene\tScore-Difference"]
    for result in result_list:
        lines.append("%s\t%s" % (result[0], str(result[1])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    name_b = (fileB.rsplit("/", 1)[1]).rsplit(".", 1)[0]
    outfile = "%s/Gene_Diff_%s_%s.tsv" % (fileA.rsplit("/", 1)[0], name_a, name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()


# called by structman
def main(sess_id, output_path, config):
    db_name = config.db_name
    db_address = config.db_address
    db_password = config.db_password
    db_user_name = config.db_user_name
    go = config.go
    godiff = config.godiff
    classification = config.classification
    path = config.path
    pathdiff = config.pathdiff
    do_modelling = config.do_modelling
    ligand_file = config.ligand_file
    multi_modelling = config.multi_modelling

    mod_per_gene = config.mod_per_gene
    infile = ''
    tanimoto_cutoff = config.tanimoto_cutoff
    distance_threshold = config.milieu_threshold

    if sess_id == 0 or sess_id is None:
        session_id = database.getSessionId(infile, config)
    else:
        session_id = sess_id

    db, cursor = config.getDB()

    if infile == '':
        infile = database.getSessionFile(session_id, config)

    db.close()

    session_name = (infile.rsplit("/", 1)[1]).rsplit(".", 1)[0]

    t0 = time.time()

    if classification and not config.skip_main_output_generation:
        t00 = time.time()
        classification_package.classificationOutput(config, output_path, session_name, session_id)
        t01 = time.time()
        if config.verbosity >= 2:
            print("Time for classificationOutput: ", t01 - t00)
        #for classfile in classfiles:
        #    out_utils.classDistributionFromFile(classfile, output_path, session_name, config)
        #    out_utils.classDistributionFromFile(classfile, output_path, session_name, config, rin_classes=True)
        t02 = time.time()
        if config.verbosity >= 2:
            print("Time for producing classification distributions: ", t02 - t01)


    t1 = time.time()
    if config.verbosity >= 2:
        print("Time for producing classification file: ", t1 - t0)

    #if config.indels_given_by_input:
        #if config.verbosity >= 2:
        #    print('Starting indel analysis output generation')
        #indel_package.create_indel_results_table(config, output_path, session_name, session_id)

    db, cursor = config.getDB()

    if go:
        database.goTermAnalysis(session_id, "%s/%s.goterm.tsv" % (output_path, session_name), db, cursor)
        if godiff:
            files = os.listdir(output_path)
            go_files = []
            for f in files:
                if '.goterm.tsv' in f:
                    go_files.append(f)
            print(go_files)
            print(files)
            if len(go_files) == 2:
                fileA = "%s/%s" % (output_path, go_files[0])
                fileB = "%s/%s" % (output_path, go_files[1])
                godiffAna(fileA, fileB)
    if path:
        database.pathwayAnalysis(session_id, "%s/%s.pathway.tsv" % (output_path, session_name), db, cursor)
        if pathdiff:
            files = os.listdir(output_path)
            path_files = []
            for f in files:
                if '.pathway.tsv' in f:
                    path_files.append(f)
            if len(path_files) == 2:
                fileA = "%s/%s" % (output_path, path_files[0])
                fileB = "%s/%s" % (output_path, path_files[1])
                pathdiffAna(fileA, fileB)
    if do_modelling:
        modelling.massModel(session_id, db, cursor, output_path, total_models=0, model_per_gene=int(mod_per_gene), multiple_mutations=multi_modelling)

    cursor.close()
    db.close()
    if ligand_file is None:
        try:
            ligand_file_names = os.listdir("%s/ligands" % infile.rsplit("/", 1)[0])
            ligand_files = []
            for ligand_file_name in ligand_file_names:
                ligand_files.append("%s/ligands/%s" % (infile.rsplit("/", 1)[0], ligand_file_name))
        except:
            ligand_files = []
    else:
        ligand_files = [ligand_file]
    for ligand_file in ligand_files:
        t0 = time.time()
        anno_dict = babel.ligandAnalyzer(ligand_file, session_id, db_name, db_address, db_user_name, db_password, cutoff=tanimoto_cutoff, distance_threshold=distance_threshold)
        t1 = time.time()
        babel.writeReport(anno_dict, "%s/Ligand_Report_%s_%s.tsv" % (output_path, ligand_file.rsplit('/', 1)[1].rsplit(".", 1)[0], session_name), db_name, db_address, db_user_name, db_password)
        t2 = time.time()
        if config.verbosity >= 2:
            print("Time for ligandAnalyzer: ", t1 - t0)
            print("Time for writeReport: ", t2 - t1)

    return 0
