#!/usr/bin/python3
import gzip
import zstd
import os
import subprocess
import sys
import shutil
import traceback
import time

from multiprocessing import Manager, Process

from structman.base_utils.config_class import Config
from structman.base_utils.base_utils import fetch_af_msa, decode_zstd
from structman.lib.rin_utils import calcRIN, get_entry_path
from structman.scripts.updateAlphafoldModelDB import check_subfolder, call_af_model_rin_calc

#lowercase_order = {y:x for x,y in lowercase_map.iteritems()}

# forbidden symbols: : _ ( ) . -

# lowercase_map = {'a':'!','b':'@','c':'#'}#,'d':'$','e':'%','f':'^','g':'&','h':'*','i':'+','j':'=','k':'{','l':'}','m':'[','n':']','o':';','p':'<','q':'>','r':',','s':'?','t':'`','u':'~','v':'|','w':'/','x':'"','y':"'"}


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

def callRINmaker(
        infile_path: str,
        outfile_path: str,
        rinmaker_executable: str
        ) -> None:
    p: subprocess.Popen = subprocess.Popen([rinmaker_executable, '-i', infile_path, '-o', outfile_path, 'rin'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate()
    return


def createRinProc(config: Config, in_queue, lock, i, remove_tmp_files, path_to_rindb, rinerator_path, errorlog):
    with lock:
        in_queue.put(None)
    while True:
        with lock:
            in_t = in_queue.get()
            if in_t is None:
                #print('Terminate Rinerator Process: ',i)
                return
            
        if isinstance(in_t[0], int):
             (struct_compr, struct_file, afdb_id) = in_t
             call_af_model_rin_calc(config, struct_compr, struct_file, afdb_id)
             continue
        (pdbgz_path, structure_id, path_to_model_db) = in_t
        if structure_id is None:
            raise TypeError(f'Error in createRinProc: {structure_id=} {pdbgz_path=} {path_to_model_db=}')
        try:
            out_path, path_stem, structure_id = get_entry_path(structure_id, path_to_rindb, path_to_model_db=path_to_model_db)
            
            if not os.path.exists(out_path):
                os.mkdir(out_path)

            #n_sif_file = f"{path_stem}.sif"

            if pdbgz_path[-8:] == '.pdb.zst':
                f = open(pdbgz_path, 'rb')
                try:
                    page = decode_zstd(f.read())
                except zstd.Error as err:
                    print(f'Could not decode {pdbgz_path}, {err=}')
                    page = None
            else:
                f = gzip.open(pdbgz_path, 'rb')
                page = f.read()
            f.close()

            calcRIN(config, page, out_path, path_stem, structure_id, rinerator_path, remove_tmp_files)

        except (KeyError, FileNotFoundError) as err:
            raise KeyError(f'{structure_id=}\n{pdbgz_path=}\n{err=}')
        
        except BaseException as err:
            error = f'{structure_id=}\n{pdbgz_path=}\n{err=}'
            print(f'{error} {errorlog=}')
            errortext = f"###############################################################################\n{error}\n###############################################################################\n"

            with lock:
                f = open(errorlog, 'a')
                f.write(errortext)
                f.close()


def check_old_folder(structure_id, rin_db_path, path_to_model_db):
    out_path, _, _ = get_entry_path(structure_id, rin_db_path, path_to_model_db=path_to_model_db)
    topfolder, subfolder = out_path.rsplit('/',1)
    upper_out_path = f'{topfolder}/{subfolder.upper()}'
    if os.path.isdir(upper_out_path):
        shutil.rmtree(upper_out_path)                


def check_entry(structure_id, path_to_model_db, path_to_rindb):
    out_path, path_stem, _ = get_entry_path(structure_id, path_to_rindb, path_to_model_db=path_to_model_db)
    
    if not os.path.exists(out_path):
        return False

    n_sif_file = f"{path_stem}.sif.gz"
    return os.path.isfile(n_sif_file)


def process_loop_model_db(config, model_db_path, in_queue, N, fetch_msas = False, from_scratch=False):
    for topfolder in os.listdir(model_db_path):
        if len(topfolder) > 2:
            continue
        if not os.path.isdir(f'{model_db_path}/{topfolder}'):
            continue
        for sub_folder in os.listdir(f'{model_db_path}/{topfolder}'):
            N += check_subfolder(config, sub_folder, topfolder, fetch_msas= fetch_msas, inqueue=in_queue, from_scratch=from_scratch)

    return in_queue, N

def process_loop_complex_model_db(config, model_db_path, in_queue, N, from_scratch=False):
    for topfolder in os.listdir(model_db_path):
        if len(topfolder) > 2:
            continue
        if not os.path.isdir(f'{model_db_path}/{topfolder}'):
            continue
        for sub_folder in os.listdir(f'{model_db_path}/{topfolder}'):
            for fn in os.listdir(f'{model_db_path}/{topfolder}/{sub_folder}'):
                if fn.count('.pdb.gz') == 1 or fn[-8:] == '.pdb.zst':
                    pdbgz_path = f'{model_db_path}/{topfolder}/{sub_folder}/{fn}'
                    
                    if fn[-8:] == '.pdb.zst':
                        model_id = fn[:-8]
                    else:
                        model_id = fn[:-7]

                    if os.path.getsize(pdbgz_path) > 50 * 1024 * 1024:
                        continue

                    rin_file = f'{model_db_path}/{topfolder}/{sub_folder}/{model_id}_intsc.ea.gz'

                    if os.path.exists(rin_file) and not from_scratch:
                        continue

                    in_queue.put((pdbgz_path, model_id, model_db_path))
                    N += 1

    return in_queue, N

def main(
        fromScratch=False,
        pdb_p='',
        rin_db_path='',
        n_proc=32,
        rinerator_base_path='',
        update_alphafold_db = False,
        update_complex_model_db = False,
        config = None):
    

    bio_assembly_path = '%s/data/biounit/PDB/divided' % pdb_p
    AU_path = "%s/data/structures/divided/pdb" % pdb_p
    status_path = "%s/data/status" % pdb_p

    rinerator_path = "%s/get_chains.py" % rinerator_base_path
    het_dict_path = "%s/reduce_wwPDB_het_dict.txt" % rinerator_base_path

    base_path = rin_db_path

    errorlog = f"{rin_db_path}/createRINdb_errorlog.txt"

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

    print(f'Call of createRINdb.main: {fromScratch=} {rin_db_path=} {time_of_last_update=} {status_path=}')

    if not fromScratch:
        recently_modified_structures = get_recently_modified_structures(time_of_last_update, status_path)
    else:
        recently_modified_structures = None

    os.environ["REDUCE_HET_DICT"] = het_dict_path

    if config is None:
        num_of_proc = n_proc
    else:
        num_of_proc = config.proc_n

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

                structure_id = fn.replace('.pdb1.gz', '')
                bio_pdbs.add(structure_id)

                check_old_folder(structure_id, rin_db_path, config.path_to_model_db)                

                if recently_modified_structures is not None:
                    if structure_id not in recently_modified_structures and check_entry(structure_id, config.path_to_model_db, rin_db_path):
                        continue

                in_queue.put((pdbgz_path, structure_id, None))
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
                structure_id = f'{fn[3:7].lower()}_au'

                check_old_folder(structure_id, rin_db_path, config.path_to_model_db)

                if recently_modified_structures is not None:
                    if structure_id not in recently_modified_structures and check_entry(structure_id, config.path_to_model_db, rin_db_path):
                        continue

                in_queue.put((pdbgz_path, structure_id, None))
                N += 1

    if config is not None:
        process_model_db = config.model_db_active
        if (config.complex_model_db_path is not None) and (config.complex_model_db_path != ''):
            process_complex_model_db = True
        else:
            process_complex_model_db = False
    else:
        process_model_db = False
        process_complex_model_db = False

    if process_model_db and update_alphafold_db:
        in_queue, afdb_N = process_loop_model_db(config, config.path_to_model_db, in_queue, 0, from_scratch=fromScratch)
    else:
        afdb_N = 0
    if process_complex_model_db and update_complex_model_db:
        in_queue, complex_model_db_N = process_loop_complex_model_db(config, config.complex_model_db_path, in_queue, 0, from_scratch=fromScratch)
    else:
        complex_model_db_N = 0

    print(f'Creating RINs for {N=} {afdb_N=} {complex_model_db_N=} structures. {errorlog=} {num_of_proc=}')

    f = open(errorlog, 'w')
    f.write('')
    f.close()

    processes = {}
    for i in range(1, min([num_of_proc, N+afdb_N+complex_model_db_N]) + 1):
        p = Process(target=createRinProc, args=(config, in_queue, lock, i, True, base_path, rinerator_path, errorlog))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()
    return recently_modified_structures


if __name__ == "__main__":
    #pdb_id = '1abr'
    # test_single_file(pdb_id)
    main(fromScratch=True)
