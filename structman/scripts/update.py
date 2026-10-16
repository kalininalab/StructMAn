#!/usr/bin/python3
import os
import subprocess
import sys

# if running as script, add local structman package to path
if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__))))))
import structman

from structman.scripts import createPdbBaDb, updateMappingDB, updateAlphafoldModelDB, createCustomDb, updateComplexModelDB
from structman.scripts.createMMDB import update_microminer_db
from structman.base_utils.base_utils import create_mmseqs_index


def main(
            config, skipUpdatePDB=False, skip_rindb=False, rin_fromScratch=False, update_mapping_db = False,
            mapping_db_from_scratch = False, update_mapping_db_keep_raw_files = False, update_alphafold_db = False,
            update_source = None, update_microminer = False, microminer_from_scratch = False, check_search_db = False, update_complex_model_db = False
        ):
    skipStructureDBs = False

    verbose = True

    pdb_path = config.pdb_path
    pdb_update_script = config.pdb_sync_script
    print(f'================================\nStructMAn Update Script ({config.proc_n=}):\n Update PDB: {not skipUpdatePDB}\n Update RIN DB: {not skip_rindb}\n Update AF DB: {update_alphafold_db}\n Update Mapping DB: {update_mapping_db}\n Check Search DB: {check_search_db}\n Update Complex Model DB: {update_complex_model_db}\n================================')

    rinerator_base_path = config.rinerator_base_path
    rin_db_path = config.rin_db_path

    mmseqs2_db_path = config.mmseqs2_db_path
    search_db_base_path = mmseqs2_db_path.rsplit('/', 1)[0]
    mmseqs2_tmp = config.mmseqs_tmp_folder
    recently_modified_structures = set()
    
    if not skipStructureDBs:
        if not skipUpdatePDB:
            if not os.path.exists(pdb_path):
                print('ERROR: Did not found local pdb')
                return
            # Set the BASE_DIR variable in the sync script
            f = open(pdb_update_script, 'r')
            lines = f.readlines()
            f.close()
            newlines = []
            for line in lines:
                if line.count('BASE_DIR=') == 1:
                    line = 'BASE_DIR="%s"\n' % pdb_path
                newlines.append(line)
            f = open(pdb_update_script, 'w')
            f.write(''.join(newlines))
            f.close()

            # update local pdb
            if verbose:
                p = subprocess.Popen([pdb_update_script])
            else:
                FNULL = open(os.devnull, 'w')
                p = subprocess.Popen([pdb_update_script], stderr=FNULL, stdout=FNULL)
            p.wait()
            sync_returncode = p.returncode

            p = subprocess.Popen(['chmod', '777', '-R', '.'], cwd=pdb_path)
            p.wait()

            # The sync script exits non-zero only after it has already retried the
            # failing directories, so reaching here means the mirror could not be
            # completed. Everything below is built by walking that mirror, so
            # carrying on would spend hours producing a search database and a RIN
            # database that are silently short by exactly the entries that failed
            # to download - and the next run would have to redo all of it anyway.
            # Stop instead, with a non-zero exit code so a scripted update fails
            # loudly rather than looking like it worked.
            if sync_returncode != 0:
                print(
                    f'\n!!! ERROR !!!\nThe PDB sync script exited with '
                    f'{sync_returncode}: the local PDB mirror at {pdb_path} is '
                    f'INCOMPLETE (see the directory list above).\nStopping the '
                    f'update: anything built from a partial mirror - the MMseqs2 '
                    f'search database, the RIN database - would be missing the '
                    f'same entries.\nRe-run the update to resume the download; '
                    f'rsync continues where it stopped.\n'
                )
                sys.exit(1)

            print('Update PDB done')

        if update_alphafold_db:
            updateAlphafoldModelDB.main(config)

        if update_complex_model_db:
            updateComplexModelDB.main(config)

        if not skip_rindb:
            # update rin db
            recently_modified_structures = structman.lib.createRINdb.main(
                fromScratch=rin_fromScratch,
                pdb_p=pdb_path,
                rin_db_path=rin_db_path,
                update_alphafold_db=update_alphafold_db,
                update_complex_model_db=update_complex_model_db,
                n_proc=config.proc_n,
                rinerator_base_path=rinerator_base_path,
                config = config
                )

            p = subprocess.Popen(['chmod', '777', '-R', '.'], cwd=rin_db_path)
            p.wait()

            print('Update RIN db done')

            if recently_modified_structures is not None:
                print('Recently modified structures: ', len(recently_modified_structures), recently_modified_structures)

    if not skipUpdatePDB:
        pdb_fasta_name = 'pdbba_mmseqs2'
        config.pdb_fasta_path = f'{search_db_base_path}/{pdb_fasta_name}'
        model_db_fasta_name = 'model_db_mmseqs2'
        config.model_db_fasta_path = f'{search_db_base_path}/{model_db_fasta_name}'

        # update pdbba for mmseqs2
        if (config.custom_db_path is not None) and (config.custom_db_path != ''):
            createCustomDb.main(config)
        else:
            createPdbBaDb.main(config, update_source = update_source)

        print("Update search database fasta for MMseqs2 done")

        create_mmseqs_index(pdb_fasta_name, search_db_base_path, mmseqs2_tmp)

        print("Search database for MMseqs2 created!")

        if config.model_db_fasta_created:

            create_mmseqs_index(model_db_fasta_name, search_db_base_path, mmseqs2_tmp)

            print("Model search database for MMseqs2 created!")

    if update_alphafold_db:
        model_db_fasta_name = 'model_db_mmseqs2'
        config.model_db_fasta_path = f'{search_db_base_path}/{model_db_fasta_name}'
        createPdbBaDb.update_model_db(config, update_source = update_source, init_ray = True)
        if config.model_db_fasta_created:

            create_mmseqs_index(model_db_fasta_name, search_db_base_path, mmseqs2_tmp)

            print("Model search database for MMseqs2 created!")
    
    if update_complex_model_db:
        complex_model_db_fasta_name = 'complex_model_db_mmseqs2'
        config.complex_model_db_fasta_path = f'{search_db_base_path}/{complex_model_db_fasta_name}'
        createPdbBaDb.create_complex_model_db_fasta(config, update_source = update_source)
        
        if config.complex_model_db_fasta_created:

            create_mmseqs_index(complex_model_db_fasta_name, search_db_base_path, mmseqs2_tmp)

            print("Complex Model search database for MMseqs2 created!")

    if check_search_db:
        if os.path.exists(mmseqs2_db_path):
            pass
        else:
            compressed_pdb_fasta_file = structman.settings.PDB_FASTA
            pdb_fasta_name = 'pdbba_mmseqs2'
            config.pdb_fasta_path = f'{search_db_base_path}/{pdb_fasta_name}'
            if not os.path.exists(config.pdb_fasta_path):
                with open(config.pdb_fasta_path, 'w') as outfile:
                    p = subprocess.run(['gunzip', '-c', compressed_pdb_fasta_file], stdout=outfile)
            create_mmseqs_index(pdb_fasta_name, search_db_base_path, mmseqs2_tmp)


    # update the mapping database
    if update_mapping_db:
        updateMappingDB.main(config, fromScratch = mapping_db_from_scratch, update_mapping_db_keep_raw_files = update_mapping_db_keep_raw_files)

    if update_microminer:

       update_microminer_db(config, fromScratch = microminer_from_scratch)

    # update the human proteome mmseqs db, TODO if we want a simple mutation-calling for fasta inputs.

    print('Hurray!')


if __name__ == "__main__":
    config_path = sys.argv[1]
    if not os.path.isfile(config_path):
        print('ERROR: Need path to config file as second argument.')
        sys.exit(1)
    config = structman.Config(config_path, external_call=True)
    main(config)
