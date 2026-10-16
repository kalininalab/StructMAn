import contextlib
import multiprocessing
import queue
import sys
import os
import traceback

from structman.base_utils.base_utils import is_alphafold_model
from structman.lib.rinerator import get_chains
from structman.lib import centrality
from structman.base_utils import residue_consts
from structman.base_utils.chain_consts import CHAIN_ID_ORDER as chain_order
from structman.base_utils.config_class import Config
try:
    from inflow_outflow_centrality import utilities as in_out_flow
except ImportError:
    pass

GET_CHAINS_TIMEOUT = 2 * 60 * 60  # 2 hours in seconds
GET_CHAINS_RESULT_TIMEOUT = 60  # seconds to wait for the result of an already exited worker


def _get_chains_worker(result_queue, tmp_pdb, out_path, tmp_chain, tmp_ligands, reduce_cmd, probe_cmd, structure_id, custom_db_mode, suppress_output):
    # BaseException, not Exception: get_chains.main calls sys.exit() on its failure paths
    try:
        if suppress_output:
            with contextlib.redirect_stdout(None), contextlib.redirect_stderr(None):
                get_chains.main(tmp_pdb, out_path, tmp_chain, tmp_ligands, True, reduce_cmd, probe_cmd, structure_id, custom_db_mode)
        else:
            get_chains.main(tmp_pdb, out_path, tmp_chain, tmp_ligands, True, reduce_cmd, probe_cmd, structure_id, custom_db_mode)
    except BaseException:
        # send the formatted traceback, an exception object itself might not be picklable
        result_queue.put(traceback.format_exc())
        return
    result_queue.put(None)


def _call_get_chains_with_timeout(tmp_pdb, out_path, tmp_chain, tmp_ligands, reduce_cmd, probe_cmd, structure_id, custom_db_mode, suppress_output):
    result_queue = multiprocessing.Queue()
    p = multiprocessing.Process(
        target=_get_chains_worker,
        args=(result_queue, tmp_pdb, out_path, tmp_chain, tmp_ligands, reduce_cmd, probe_cmd, structure_id, custom_db_mode, suppress_output)
    )
    p.start()
    p.join(timeout=GET_CHAINS_TIMEOUT)
    if p.is_alive():
        p.terminate()
        p.join()
        raise TimeoutError(f'get_chains.main timed out after 2 hours for structure {structure_id}')

    if p.exitcode != 0:
        # the worker died before it could report anything (killed by the OOM killer, segfault in reduce/probe, ...)
        return -1, f'get_chains.main process for {structure_id} died with exitcode {p.exitcode}'

    try:
        result = result_queue.get(timeout=GET_CHAINS_RESULT_TIMEOUT)
    except queue.Empty:
        return -1, f'get_chains.main process for {structure_id} exited without reporting a result'

    if result is not None:
        return -1, f'get_chains.main failed for {structure_id}:\n{result}'

    return 0, None

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

        if record_name == 'MODEL':
            current_model_id = line[10:14]
            if current_model_id not in multi_model_chain_dict:
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

        if len(line) > 26:
            

            #atom_nr = line[6:11].replace(" ", "")
            #atom_name = line[12:16].replace(" ", "")
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
                if (chain_id, res_nr) not in modres_map:
                    modres_map[(chain_id, res_nr)] = res_name
            if record_name == "ATOM" or record_name == "HETATM":

                if chain_id not in chain_ids:
                    chain_ids.add(chain_id)
                    multi_model_chain_dict[current_model_id][chain_id] = chain_id
                elif multi_model_mode:
                    if chain_id not in multi_model_chain_dict[current_model_id]:
                        new_chain_id = None
                        while new_chain_id is None:
                            if chain_order[chain_order_pos_marker] not in chain_ids:
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
                        if (chain_id, res_nr) not in modres_map and res_name not in residue_consts.THREE_TO_ONE:
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

def calcRIN(
        config: Config,
        page: bytes,
        out_path: str,
        path_stem: str,
        structure_id: str,
        rinerator_path: str,
        remove_tmp_files: bool,
        structure_path: str | None = None):

    if config.verbosity >= 4:
        config.logger.info(f'Call of calcRIN: {structure_id=}, {out_path=}, {path_stem=} {structure_path=} {remove_tmp_files=}')

    #print(f'Call of calcRIN: {structure_id=}, {out_path=}, {path_stem=} {structure_path=} {remove_tmp_files=}')

    parse_out = parsePDB(page)
    if parse_out is None:
        return
    original_chains, ligands, page = parse_out

    if config.custom_db_path != '':
        custom_db_mode = True
        #custom_path_list = os.listdir(config.custom_db_path)
    else:
        custom_db_mode = False
                
    if len(original_chains) == 0:
        return
    if structure_path is None:
        tmp_pdb = f"{path_stem}.pdb"
        f = open(tmp_pdb, 'w')
        f.write(page)
        f.close()
    else:
        tmp_pdb = structure_path

    tmp_chain = f"{out_path}/chains_{structure_id}.txt"
    tmp_ligands = f"{out_path}/ligands_{structure_id}.txt"

    f = open(tmp_chain, 'w')
    f.write(",".join(original_chains))
    f.close()

    f = open(tmp_ligands, 'w')
    f.write("\n".join([','.join(x) for x in ligands]))
    f.close()
    rinerator_base_path = rinerator_path.rsplit('/', 1)[0]

    reduce_cmd = '%s/reduce' % rinerator_base_path
    probe_cmd = '%s/probe' % rinerator_base_path

    reduce_file = f"{path_stem}_h.ent"
    probe_file = f"{path_stem}_h.probe"

    try:
        ret_code, ret_message = _call_get_chains_with_timeout(tmp_pdb, out_path, tmp_chain, tmp_ligands, reduce_cmd, probe_cmd, structure_id, custom_db_mode, config.verbosity < 4)
    except UnicodeDecodeError:
        tmp_pdb = f"{path_stem}.pdb"
        f = open(tmp_pdb, 'w')
        f.write(page)
        f.close()
        ret_code, ret_message = _call_get_chains_with_timeout(tmp_pdb, out_path, tmp_chain, tmp_ligands, reduce_cmd, probe_cmd, structure_id, custom_db_mode, config.verbosity < 4)
    except BaseException as err:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        err_message = f"\nRIN calc Error: \n{structure_id=} {original_chains=}\n{structure_path=}\n{tmp_pdb=}\n{out_path=}\n{err=}\n{e}\n{f}\n{g}"
        try:
            config.logger.info(err_message)
        except AttributeError:
            print(err_message)
        if config.verbosity >= 5:
            config.logger.info(page)
        if os.path.exists(reduce_file):
            os.remove(reduce_file)
        if os.path.exists(probe_file):
            os.remove(probe_file)
        if remove_tmp_files:
            if structure_path is None:
                if os.path.exists(tmp_pdb):
                    os.remove(tmp_pdb)
            if os.path.exists(tmp_chain):
                os.remove(tmp_chain)
            if os.path.exists(tmp_ligands):
                os.remove(tmp_ligands)
        return

    if ret_code == -1:
        config.logger.info(f'RIN calculation failed: {structure_id=} {tmp_pdb=} {ret_message}')
        if os.path.exists(reduce_file):
            os.remove(reduce_file)
        if os.path.exists(probe_file):
            os.remove(probe_file)
        if remove_tmp_files:
            if structure_path is None:
                if os.path.exists(tmp_pdb):
                    os.remove(tmp_pdb)
            if os.path.exists(tmp_chain):
                os.remove(tmp_chain)
            if os.path.exists(tmp_ligands):
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
            if os.path.exists(tmp_pdb):
                os.remove(tmp_pdb)
        if os.path.exists(tmp_chain):
            os.remove(tmp_chain)
        if os.path.exists(tmp_ligands):
            os.remove(tmp_ligands)
    else:
        config.logger.info(f'{tmp_pdb=}')
        config.logger.info(f'{tmp_chain=}')
        config.logger.info(f'{tmp_ligands=}')


    sif_file = f"{path_stem}_h.sif"
    n_sif_file = f"{path_stem}.sif"

    intsc_file = f"{path_stem}_h_intsc.ea"
    n_intsc_file = f"{path_stem}_intsc.ea"
    nrint_file = f"{path_stem}_h_nrint.ea"
    n_nrint_file = f"{path_stem}_nrint.ea"
    res_file = f"{path_stem}_h_res.txt"
    n_res_file = f"{path_stem}_res.txt"

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

    rin_files = (sif_file, intsc_file, nrint_file, res_file)
    missing_files = [rin_file for rin_file in rin_files if not os.path.isfile(rin_file)]
    if len(missing_files) > 0:
        config.logger.info(f'RINerator produced no complete output: {structure_id=} {tmp_pdb=} {missing_files=}')
        for rin_file in rin_files:
            if os.path.isfile(rin_file):
                os.remove(rin_file)
        return

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

    centrality.main(n_sif_file, pdb_id=structure_id)

    call_in_out_flow(n_intsc_file, f"{path_stem}_in_out_flow.csv")

    os.system("gzip %s" % n_sif_file)
    os.system("gzip %s" % n_intsc_file)
    os.system("gzip %s" % n_nrint_file)
    os.system("gzip %s" % n_res_file)

def call_in_out_flow(n_intsc_file: str, io_out_path: str):
    interactions_df = in_out_flow.read_interactions_ea(n_intsc_file)

    io_df = in_out_flow.flow_metric(interactions_df, None, 'both')

    io_df.to_csv(io_out_path,index=False)

def get_entry_path(
        structure_id: str, 
        path_to_rindb: str, 
        path_to_model_db: str | None = None, 
        custom_db_path: str | None = None, 
        model_path: str | None = None, 
        temp_folder: str | None = None
        ) -> tuple[str, str]:

    if model_path is not None: #First check, if this refering to a RIN of a custom model (structural analysis of models)
        #folder_path = temp_folder
        if model_path[-4:] == '.pdb':
            path_stem = model_path[:-4] # remove .pdb from the end of the model path
        else:
            path_stem = model_path[:-7] # .pdb.gz
        folder_path = os.path.dirname(path_stem)
        structure_id = path_stem.split('/')[-1].split('.')[0]
    elif custom_db_path is not None: #This is for custom structure DBs
        folder_path = custom_db_path
        #custom_files_list = os.listdir(config.temp_folder)     
        path_stem = f"{folder_path}/{structure_id}"
       
    elif path_to_model_db is not None and is_alphafold_model(structure_id): #Not a pdb structure, but a model
        try:
            uniprot_ac = structure_id.split('-')[1]
        except IndexError:
            uniprot_ac = structure_id.split('_')[0]
        topfolder_id = uniprot_ac[-2:]
        subfolder_id = uniprot_ac[-4:]
        folder_path = f'{path_to_model_db}/{topfolder_id}/{subfolder_id}'
        path_stem = f"{folder_path}/{structure_id}"
    elif path_to_rindb == '':
        folder_path = temp_folder
        path_stem = f"{folder_path}/{structure_id}"
    else:
        pdb_id: str = structure_id.replace('_AU', '').replace('_au','').lower().split('.')[0]
        pdb_id_middle_key = pdb_id[1:-1]
        folder_path = f"{path_to_rindb}/{pdb_id_middle_key}/{structure_id.lower()}"
        path_stem = f"{folder_path}/{structure_id.lower()}"

    return folder_path, path_stem, structure_id