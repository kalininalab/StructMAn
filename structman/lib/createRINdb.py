#!/usr/bin/python3
import gzip
import os
import resource

import sys
import time
import traceback
import contextlib

from multiprocessing import Manager, Process

from structman import settings
from structman.lib import centrality
from structman.lib.rinerator import get_chains
from structman.lib.sdsc.consts import residues

#lowercase_order = {y:x for x,y in lowercase_map.iteritems()}

# forbidden symbols: : _ ( ) . -

# lowercase_map = {'a':'!','b':'@','c':'#'}#,'d':'$','e':'%','f':'^','g':'&','h':'*','i':'+','j':'=','k':'{','l':'}','m':'[','n':']','o':';','p':'<','q':'>','r':',','s':'?','t':'`','u':'~','v':'|','w':'/','x':'"','y':"'"}

chain_order = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz'

def parsePDB(page):

    try:
        page = page.decode('ascii')
    except:
        return None

    lines = page.split('\n')

    original_chains = set([])
    chain_ids = set()

    ligands = set()
    modres_map = {}

    new_lines = []

    not_crystal = False
    multi_model_mode = False
    multi_model_chain_dict = {}
    chain_order_pos_marker = 0
    asymetric_chain_number = None
    current_model_id = 0
    multi_model_chain_dict[current_model_id] = {}

    #changed_chains = []
    firstAltLoc = None
    for line in lines:
        record_name = line[0:6].replace(" ", "")

        if record_name == 'EXPDTA':
            words = line.split()
            if len(words) < 3:
                continue
            if words[2] == 'NMR':
                not_crystal = True

            if words[1] == 'SOLUTION' and (words[2] == 'SCATTERING' or words[2] == 'SCATTERING;' or words[2] == 'NMR;' or words[2] == 'NMR'):
                not_crystal = True

            if line.count('NMR') > 0:
                not_crystal = True
            new_lines.append(line)
            continue

        if len(line) > 26:
            if record_name == 'MODEL':
                current_model_id = line[10:14]
                if not current_model_id in multi_model_chain_dict:
                    multi_model_chain_dict[current_model_id] = {}
                if not multi_model_mode:
                    new_lines.append(line)
                elif len(multi_model_chain_dict) * asymetric_chain_number >= len(chain_order):
                    break
                continue

            elif record_name == 'ENDMDL':
                if not_crystal:
                    new_lines.append(line)
                    break
                elif not multi_model_mode:
                    multi_model_mode = True
                    asymetric_chain_number = len(chain_ids)
                continue

            atom_nr = line[6:11].replace(" ", "")
            atom_name = line[12:16].replace(" ", "")
            res_name = line[17:20].replace(" ", "")
            chain_id = line[21]

            # if chain_id in lowercase_map:
            #    chain_id = lowercase_map[chain_id]
            #    changed_chains.append(chain_id)

            res_nr = line[22:27].replace(" ", "")

            if record_name == 'MODRES':
                chain_id = line[16]
                res_nr = line[18:23].replace(" ", "")
                res_name = line[24:27].replace(" ", "")
                if not (chain_id, res_nr) in modres_map:
                    modres_map[(chain_id, res_nr)] = res_name
            if record_name == "ATOM" or record_name == "HETATM":

                if chain_id not in chain_ids:
                    chain_ids.add(chain_id)
                    multi_model_chain_dict[current_model_id][chain_id] = chain_id
                elif multi_model_mode:
                    if not chain_id in multi_model_chain_dict[current_model_id]:
                        new_chain_id = None
                        while new_chain_id is None:
                            if not chain_order[chain_order_pos_marker] in chain_ids:
                                new_chain_id = chain_order[chain_order_pos_marker]
                            chain_order_pos_marker += 1
                        multi_model_chain_dict[current_model_id][chain_id] = new_chain_id
                        chain_id = new_chain_id
                    else:
                        chain_id = multi_model_chain_dict[current_model_id][chain_id]

                altLoc = line[16]
                if firstAltLoc is None and altLoc != ' ':
                    firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                    continue
                if len(res_name) < 3:
                    ligands.add((res_name, chain_id, res_nr))
                else:
                    if record_name == "HETATM":
                        if not (chain_id, res_nr) in modres_map and res_name not in residues.THREE_TO_ONE:
                            ligands.add((res_name, chain_id, res_nr))
                        else:
                            original_chains.add(chain_id)
                    else:
                        original_chains.add(chain_id)

            line = '%s%s%s' % (line[:21], chain_id, line[22:])

        new_lines.append(line)

    page = '\n'.join(new_lines)
    # return original_chains,ligands,page,changed_chains
    return original_chains, ligands, page


def get_recently_modified_structures(time_of_last_update, status_path):
    if time_of_last_update is None:
        return None

    time_struct = time.strptime(time_of_last_update)

    recently_modified_structures = set()

    for date_sub_folder in os.listdir(status_path):
        date_sub_folder_path = '%s/%s' % (status_path, date_sub_folder)
        if not os.path.isdir(date_sub_folder_path):
            continue
        if date_sub_folder_path == 'latest':
            continue
        year = int(date_sub_folder[:4])
        month = int(date_sub_folder[4:6])
        day = int(date_sub_folder[6:8])
        if year < time_struct.tm_year:
            continue
        elif month < time_struct.tm_mon:
            continue
        elif day < time_struct.tm_mday:
            continue
        for status_file in date_sub_folder_path:
            if status_file[:8] != 'modified':
                continue
            f = open('%s/%s' % (date_sub_folder_path, status_file), 'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                line.replace('\n','')
                recently_modified_structures.add(line)
    return recently_modified_structures

'''
def changeBackChains(changed_chains,n_sif_file,n_intsc_file,n_nrint_file,n_res_file):
    for fn in [n_sif_file,n_intsc_file,n_nrint_file,n_res_file]:
        f = open(fn,'r')
        page = f.read()
        f.close()

        for special_chain_id in changed_chains:
            original_chain = lowercase_order[special_chain_id]
            page.replace(special_chain_id,original_chain)

        f = open(fn,'w')
        f.write(page)
        f.close()

    return
'''


def calcRIN(config, page, out_path, pdb_id, rinerator_path, remove_tmp_files, structure_path=None):

    if config.verbosity >= 4:
        print(f'Call of calcRIN: {pdb_id}, {out_path}, {structure_path}')

    parse_out = parsePDB(page)
    if parse_out is None:
        return
    original_chains, ligands, page = parse_out

    pdb_id = pdb_id.split('.')[0]
    
    if config.custom_db_path != '':
        custom_db_mode = True
        #custom_path_list = os.listdir(config.custom_db_path)
    else:
        custom_db_mode = False
                
    if len(original_chains) == 0:
        return
    if structure_path is None:
        tmp_pdb = "%s/%s.pdb" % (out_path, pdb_id)
        f = open(tmp_pdb, 'w')
        f.write(page)
        f.close()
    else:
        tmp_pdb = structure_path

    tmp_chain = "%s/chains_%s.txt" % (out_path, pdb_id)
    tmp_ligands = "%s/ligands_%s.txt" % (out_path, pdb_id)

    f = open(tmp_chain, 'w')
    f.write(",".join(original_chains))
    f.close()

    f = open(tmp_ligands, 'w')
    f.write("\n".join([','.join(x) for x in ligands]))
    f.close()
    rinerator_base_path = rinerator_path.rsplit('/', 1)[0]

    reduce_cmd = '%s/reduce' % rinerator_base_path
    probe_cmd = '%s/probe' % rinerator_base_path

    reduce_file = "%s/%s_h.ent" % (out_path, pdb_id)
    probe_file = "%s/%s_h.probe" % (out_path, pdb_id)

    try:
        if config.verbosity >= 4:
            get_chains.main(tmp_pdb, out_path, tmp_chain, tmp_ligands, True, reduce_cmd, probe_cmd, pdb_id, custom_db_mode)

        else:
            with contextlib.redirect_stdout(None), contextlib.redirect_stderr(None):
                get_chains.main(tmp_pdb, out_path, tmp_chain, tmp_ligands, True, reduce_cmd, probe_cmd, pdb_id, custom_db_mode)
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        print(f"\nRIN calc Error:\n{pdb_id} {original_chains}\n{structure_path}\n{tmp_pdb}\n{out_path}\n{e}\n{f}\n{g}")
        if config.verbosity >= 5:
            print(page)
        if os.path.exists(reduce_file):
            os.remove(reduce_file)
        if os.path.exists(probe_file):
            os.remove(probe_file)
        if remove_tmp_files:
            if structure_path is None:
                os.remove(tmp_pdb)
            os.remove(tmp_chain)
            os.remove(tmp_ligands)
        return

    '''
    cmd = ['python3',rinerator_path,'-d',tmp_pdb,out_path,tmp_chain,tmp_ligands]

    FNULL = open(os.devnull, 'w')
    if verbosity >= 2:
        p = subprocess.Popen(cmd,stdout=FNULL)
        outs,errs = p.communicate()
        print('Rinerator error output for ',pdb_id,'\n',errs)
    else:
        p = subprocess.Popen(cmd,stderr=FNULL,stdout=FNULL)
        p.communicate()
    '''

    if remove_tmp_files:
        if structure_path is None:
            os.remove(tmp_pdb)
        os.remove(tmp_chain)
        os.remove(tmp_ligands)
    else:
        print(tmp_pdb)
        print(tmp_chain)
        print(tmp_ligands)


    sif_file = "%s/%s_h.sif" % (out_path, pdb_id)
    n_sif_file = "%s/%s.sif" % (out_path, pdb_id)

    intsc_file = "%s/%s_h_intsc.ea" % (out_path, pdb_id)
    n_intsc_file = "%s/%s_intsc.ea" % (out_path, pdb_id)
    nrint_file = "%s/%s_h_nrint.ea" % (out_path, pdb_id)
    n_nrint_file = "%s/%s_nrint.ea" % (out_path, pdb_id)
    res_file = "%s/%s_h_res.txt" % (out_path, pdb_id)
    n_res_file = "%s/%s_res.txt" % (out_path, pdb_id)

    '''
    if p.returncode != 0:
        try:
            os.remove(reduce_file)
            os.remove(probe_file)
            os.remove(sif_file)
            os.remove(intsc_file)
            os.remove(nrint_file)
            os.remove(res_file)
        except:
            pass
        raise IOError("RINerator failed")
    '''

    if os.path.exists(reduce_file):
        os.remove(reduce_file)
    if os.path.exists(probe_file):
        os.remove(probe_file)

    os.rename(sif_file, n_sif_file)
    os.rename(intsc_file, n_intsc_file)
    os.rename(nrint_file, n_nrint_file)
    os.rename(res_file, n_res_file)

    # if changed_chains != []:
    #    changeBackChains(changed_chains,n_sif_file,n_intsc_file,n_nrint_file,n_res_file)

    if os.path.isfile("%s.gz" % n_sif_file):
        os.remove("%s.gz" % n_sif_file)
    if os.path.isfile("%s.gz" % n_intsc_file):
        os.remove("%s.gz" % n_intsc_file)
    if os.path.isfile("%s.gz" % n_nrint_file):
        os.remove("%s.gz" % n_nrint_file)
    if os.path.isfile("%s.gz" % n_res_file):
        os.remove("%s.gz" % n_res_file)

    centrality.main(n_sif_file, pdb_id=pdb_id)

    os.system("gzip %s" % n_sif_file)
    os.system("gzip %s" % n_intsc_file)
    os.system("gzip %s" % n_nrint_file)
    os.system("gzip %s" % n_res_file)


def createRinProc(config, in_queue, lock, i, remove_tmp_files, base_path, rinerator_path, errorlog, path_to_model_db):
    with lock:
        in_queue.put(None)
    while True:
        with lock:
            in_t = in_queue.get()
            if in_t is None:
                #print('Terminate Rinerator Process: ',i)
                return
            (pdbgz_path, pdb_id) = in_t
        try:

            if (pdb_id[:3] == 'AF-'): #Not a pdb structure, but an alphafold model
                uniprot_ac = pdb_id.split('-')[1]
                topfolder_id = uniprot_ac[-2:]
                subfolder_id = uniprot_ac[-4:]
                out_path = f'{path_to_model_db}/{topfolder_id}/{subfolder_id}'
            else:
                out_path = "%s/%s/%s" % (base_path, pdb_id[1:-1], pdb_id)
            if not os.path.exists(out_path):
                os.mkdir(out_path)

            n_sif_file = "%s/%s.sif" % (out_path, pdb_id)

            f = gzip.open(pdbgz_path, 'rb')
            page = f.read()
            f.close()
            #original_chains,ligands,page,changed_chains = parsePDB(page)

            calcRIN(config, page, out_path, pdb_id, rinerator_path, remove_tmp_files)

        except:
            [e, f, g] = sys.exc_info()
            #g = traceback.format_exc(g)
            print("Error: ", e, f, g)
            error = '\n'.join([pdb_id, pdbgz_path, str(e), str(f), str(g)])
            errortext = "###############################################################################\n%s\n###############################################################################\n" % (error)

            with lock:
                f = open(errorlog, 'a')
                f.write(errortext)
                f.close()


def test_single_file(config, pdb_id):
    calculateRINsFromPdbList(config, [pdb_id], remove_tmp_files=False)


def calculateRINsFromPdbList(config, pdbs, fromScratch=True, forceCentrality=True, remove_tmp_files=True, n_proc=32):

    pdbs = set([x.lower() for x in pdbs])

    lim = 100 * 1024 * 1024 * 1024

    resource.setrlimit(resource.RLIMIT_AS, (lim, lim))

    if not os.path.isfile(settings.REDUCE_HET_DICT):
        print("%s not found" % settings.REDUCE_HET_DICT)
        sys.exit(1)

    os.environ["REDUCE_HET_DICT"] = settings.REDUCE_HET_DICT

    num_of_proc = n_proc

    manager = Manager()
    lock = manager.Lock()

    in_queue = manager.Queue()

    bio_pdbs = set()

    total_structures = 0

    subfolders = os.listdir(bio_assembly_path)  # BUG: undefined variable
    for subfolder in subfolders:
        sub_path = "%s/%s" % (bio_assembly_path, subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (base_path, subfolder)):  # BUG: undefined variable
            os.mkdir("%s/%s" % (base_path, subfolder))

        for fn in files:
            if fn.count('.pdb1.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path, fn)
                if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                    continue
                pdb_id = fn.replace('.pdb1.gz', '')
                if pdb_id not in pdbs:
                    continue

                bio_pdbs.add(pdb_id)
                in_queue.put((pdbgz_path, pdb_id))
                total_structures += 1

    subfolders = os.listdir(AU_path)  # BUG: undefined variable
    for subfolder in subfolders:

        sub_path = "%s/%s" % (AU_path, subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (base_path, subfolder)):
            os.mkdir("%s/%s" % (base_path, subfolder))

        for fn in files:
            if fn.count('.ent.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path, fn)
                if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                    continue
                pdb_id = fn[3:7]
                if not '%s_au' % pdb_id in pdbs:
                    continue
                if pdb_id in bio_pdbs:
                    continue
                in_queue.put((pdbgz_path, pdb_id))
                total_structures += 1

    print('Amount of structures for RINerator: ', total_structures)

    processes = {}
    for i in range(1, num_of_proc + 1):
        p = Process(target=createRinProc, args=(config, in_queue, lock, fromScratch, i, forceCentrality, remove_tmp_files, base_path, rinerator_path, errorlog, config.path_to_model_db))  # BUG: undefined variable
        processes[i] = p
        print('Start RINerator Process: ', i)
        p.start()
    for i in processes:
        processes[i].join()



def main(fromScratch=False, pdb_p='', rin_db_path='', n_proc=32, rinerator_base_path='', process_model_db = False, config = None):
    bio_assembly_path = '%s/data/biounit/PDB/divided' % pdb_p
    AU_path = "%s/data/structures/divided/pdb" % pdb_p
    status_path = "%s/data/status" % pdb_p

    rinerator_path = "%s/get_chains.py" % rinerator_base_path
    het_dict_path = "%s/reduce_wwPDB_het_dict.txt" % rinerator_base_path

    base_path = rin_db_path

    errorlog = "createRINdb_errorlog.txt"

    lim = 100 * 1024 * 1024 * 1024

    #resource.setrlimit(resource.RLIMIT_AS, (lim, lim))

    if not os.path.isfile(het_dict_path):
        print("%s not found" % het_dict_path)
        sys.exit(1)

    meta_file_path = '%s/meta.txt' % base_path

    if os.path.isfile(meta_file_path):
        f = open(meta_file_path, 'r')
        time_of_last_update = f.read()
        f.close()
    else:
        time_of_last_update = None

    f = open(meta_file_path,'w')
    f.write(time.asctime(time.gmtime()))
    f.close()

    if not fromScratch:
        recently_modified_structures = get_recently_modified_structures(time_of_last_update, status_path)
    else:
        recently_modified_structures = None

    os.environ["REDUCE_HET_DICT"] = het_dict_path

    num_of_proc = n_proc

    manager = Manager()
    lock = manager.Lock()

    in_queue = manager.Queue()

    bio_pdbs = set()

    subfolders = os.listdir(bio_assembly_path)

    N = 0

    for subfolder in subfolders:

        sub_path = "%s/%s" % (bio_assembly_path, subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (base_path, subfolder)):
            os.mkdir("%s/%s" % (base_path, subfolder))

        for fn in files:
            if fn.count('.pdb1.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path, fn)
                if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                    continue

                pdb_id = fn.replace('.pdb1.gz', '')
                bio_pdbs.add(pdb_id)

                if recently_modified_structures is not None:
                    if not pdb_id in recently_modified_structures:
                        continue

                in_queue.put((pdbgz_path, pdb_id))
                N += 1

    subfolders = os.listdir(AU_path)
    for subfolder in subfolders:

        sub_path = "%s/%s" % (AU_path, subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (base_path, subfolder)):
            os.mkdir("%s/%s" % (base_path, subfolder))

        for fn in files:
            if fn.count('.ent.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path, fn)
                if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                    continue
                pdb_id = fn[3:7]

                if recently_modified_structures is not None:
                    if not pdb_id in recently_modified_structures:
                        continue

                if pdb_id in bio_pdbs:
                    continue

                in_queue.put((pdbgz_path, pdb_id))
                N += 1

    if process_model_db:
        for topfolder in os.listdir(config.path_to_model_db):
            if not os.path.isdir(f'{config.path_to_model_db}/{topfolder}'):
                continue
            for sub_folder in os.listdir(f'{config.path_to_model_db}/{topfolder}'):
                for fn in os.listdir(f'{config.path_to_model_db}/{topfolder}/{sub_folder}'):
                    if fn.count('.pdb.gz') == 1:
                        pdbgz_path = f'{config.path_to_model_db}/{topfolder}/{sub_folder}/{fn}'
                        if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                            continue

                        model_id = fn[:-7]
                        rin_file = f'{config.path_to_model_db}/{topfolder}/{sub_folder}/{model_id}_intsc.ea.gz'

                        if os.path.exists(rin_file):
                            continue

                        in_queue.put((pdbgz_path, model_id))
                        N += 1

    print('Creating RINs for', N, 'structures.')

    processes = {}
    for i in range(1, min([num_of_proc, N]) + 1):
        p = Process(target=createRinProc, args=(config, in_queue, lock, i, True, base_path, rinerator_path, errorlog, config.path_to_model_db))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()
    return recently_modified_structures


if __name__ == "__main__":
    #pdb_id = '1abr'
    # test_single_file(pdb_id)
    main(fromScratch=True)
