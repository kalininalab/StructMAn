import os
import subprocess
import time

import pymysql as MySQLdb

from structman.lib.database import database


# Needs either a ligand_db or it creates one, but then it needs contact to the database
def findRelatives(config, infile, session_id, cutoff, ligand_db="", db=None, cursor=None):
    if ligand_db == "":
        if db is None:
            raise NameError("Illegal Input in findRelatives")
        session_file = database.getSessionFile(session_id, config)
        name = session_file.rsplit("/", 1)[1]
        name = name.rsplit(".", 1)[0]
        outname = "%s/%s.smi" % (infile.rsplit("/", 1)[0], name)
        database.createLigandDB(outname, session_id, db, cursor)
        ligand_db = outname

    p = subprocess.Popen(["babel", infile, ligand_db, "-ofpt"], stdout=subprocess.PIPE, universal_newlines=True)
    page = p.communicate()
    lines = page[0].split("\n")
    lines = lines[1:]

    chosen_ones = {}
    i = 0
    for line in lines:
        if line.count("=") > 0:
            lig_id = line.split()[0].replace(">", "")
            tanimoto = float(line.rsplit("=", 1)[1].replace(" ", ""))
            if tanimoto >= cutoff:
                chosen_ones[lig_id] = tanimoto
    return chosen_ones


#This is really oudated
# called by output
def ligandAnalyzer(config, infile, session_id, db_name, host, user, pw, ligand_db="", cutoff=0.05, distance_threshold=10.0):
    t0 = time.time()
    db = MySQLdb.connect(host, user, pw, db_name)
    cursor = db.cursor()
    if ligand_db == "":
        cwd = os.getcwd()
        session_file = database.getSessionFile(session_id, db, cursor)
        name = session_file.rsplit("/", 1)[1]
        name = name.rsplit(".", 1)[0]
        ligand_db = "%s/%s.smi" % (cwd, name)
        if os.path.isfile(ligand_db):
            chosen_ones = findRelatives(config, infile, session_id, float(cutoff), ligand_db=ligand_db)
        else:
            chosen_ones = findRelatives(config, infile, session_id, float(cutoff))
    else:
        chosen_ones = findRelatives(config, infile, session_id, cutoff, ligand_db=ligand_db)
    if len(chosen_ones) == 0:
        print("No relatives found")
        return {}

    t1 = time.time()
    anno_dict = database.getLigandAnnotation(chosen_ones, session_id, distance_threshold, db, cursor)

    db.close()

    t2 = time.time()

    print("Time for ligandAnalyzer Part1: ", t1 - t0)
    print("Time for ligandAnalyzer Part2: ", t2 - t1)

    return anno_dict


# called by output
def writeReport(anno_dict, outfile, db_name, host, user, pw):
    outlines = ["Ligand\tTanimoto score\tTemplate\tGene\tMutations\tDistance"]
    db = MySQLdb.connect(host, user, pw, db_name)
    cursor = db.cursor()
    for lig_name in sorted(anno_dict, key=lambda x: float(anno_dict[x][0]), reverse=True):
        temp_dict = anno_dict[lig_name][1]
        for template_id in sorted(temp_dict, key=lambda x: float(temp_dict[x][0][1][0].split(":")[1])):
            mutation_list = temp_dict[template_id]
            pdb = database.getPDB(template_id, db, cursor)
            gene_id = database.getGeneFromTemplate(template_id, db, cursor)
            sp_id = database.getUniprotAcFromId(gene_id, db, cursor)
            for (mutation_id, diststrs) in mutation_list:
                for diststr in diststrs:
                    dist = diststr.split(":")[1]
                    aac = database.getAAChange(mutation_id, db, cursor)
                    outlines.append("%s\t%s\t%s\t%s\t%s\t%s" % (lig_name, str(anno_dict[lig_name][0]), pdb, sp_id, aac, dist))
    db.close()
    if len(outlines) == 1:
        outlines.append("Sorry,\tdid\tnot\tfound\ta similar\tligand")
    print(len(outlines))
    f = open(outfile, "w")
    f.write("\n".join(outlines))
    f.close()
    return
