import os
import subprocess
import sys
import traceback
import gzip

from structman.lib.database.database_core_functions import update
from structman.base_utils.base_utils import pack

def check_instance(config, fromScratch = False):
    fresh_instance = False
    if config.mapping_db_is_set and fromScratch:
        db, cursor = config.getDB(mapping_db = True)
        sql = 'DROP DATABASE %s' % config.mapping_db
        cursor.execute(sql)
        db.close()
        fresh_instance = True
    if not config.mapping_db_is_set or fresh_instance:
        db, cursor = config.getDB(server_connection=True)

        sql = 'CREATE DATABASE %s' % config.mapping_db

        try:
            cursor.execute(sql)
            db.close()
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            print('\n'.join([str(e), str(f), str(g)]))
            if db is not None:
                db.close()
            sys.exit()

        if config.mapping_db_source_path[-3:] == '.gz':
            f = gzip.open(config.mapping_db_source_path, 'rb')
            binary = True
        else:
            f = open(config.mapping_db_source_path, 'r')
            binary = False
        lines = f.readlines()
        f.close()

        new_lines = []
        for line in lines:
            if binary:
                line = line.decode('ascii')
            if line[:4] == 'USE ':
                new_lines.append('USE `%s`;\n' % config.mapping_db)
            else:
                new_lines.append(line)

        db_file = 'tmp_mapping_database_file.sql'
        f = open(db_file, 'w')
        f.write(''.join(new_lines))
        f.close()

        cmds = ' '.join(['mysql', '-u', config.db_user_name, '-h', config.db_address, '--password=%s' % config.db_password, '<', db_file])

        p = subprocess.Popen(cmds, shell=True)
        p.wait()

        os.remove(db_file)
        
    return fresh_instance

def retrieve_raw_data(config, raw_files_folder_path, fromScratch = False):
    mapping_file_path = retrieve_data_from_uniprot(config, raw_files_folder_path, 'idmapping', 'idmapping.dat.gz', fromScratch = fromScratch)
    sequence_file_paths = []
    for seq_file_name in ['uniprot_sprot_varsplic.fasta.gz','uniprot_sprot.fasta.gz','uniprot_trembl.fasta.gz']:
        sequence_file_paths.append(retrieve_data_from_uniprot(config, raw_files_folder_path, 'complete', seq_file_name, fromScratch = fromScratch))
    return mapping_file_path, sequence_file_paths

def retrieve_data_from_uniprot(config, raw_files_folder_path, uniprot_sub_folder, uniprot_file_name, fromScratch = False):
    file_path = '%s/%s' % (raw_files_folder_path, uniprot_file_name)
    if os.path.isfile(file_path) and fromScratch:
        os.remove(file_path)
    if not os.path.isfile(file_path):
        cmds = ' '.join(['wget', 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/%s/%s' % (uniprot_sub_folder, uniprot_file_name)])
        p = subprocess.Popen(cmds, shell=True, cwd = config.tmp_folder)
        p.wait()
    return file_path

def put_seqs_to_database(seq_map, config):
    values = []
    for u_ac, seq in seq_map:
        values.append((u_ac, pack(seq)))

    update(config, 'UNIPROT', ['Uniprot_Ac','Sequence'], values, mapping_db = True)

def insert_sequence_data(config, seq_file_paths, skip_inserts = False):
    max_seqs_at_a_time = int(30000 * config.gigs_of_ram)

    for seq_file_number, seq_file in enumerate(seq_file_paths):
        if seq_file_number == 0: #The first sequence file contains the minor isoforms
            isoform_file = True
            iso_map = {}
        else:
            isoform_file = False
        with gzip.open(seq_file, 'rb') as f:
            seq_map = []
            current_seq = None
            for line in f:
                line = line.decode('ascii')
                if len(line) == 0:
                    continue
                line = line[:-1]
                if line[0] == '>':
                    if current_seq is not None:
                        seq_map.append((u_ac, current_seq))
                    if len(seq_map) == max_seqs_at_a_time:
                        if not skip_inserts:
                            put_seqs_to_database(seq_map, config)
                        seq_map = []

                    u_ac = line.split('|')[1]

                    current_seq = ''
                    if isoform_file: #save all isoform numbers
                        stem, iso_number = u_ac.split('-')
                        if not stem in iso_map:
                            iso_map[stem] = set()
                        iso_map[stem].add(iso_number)
                else:
                    current_seq += (line)
            seq_map.append((u_ac, current_seq))
            if len(seq_map) > 0:
                if not skip_inserts:
                    put_seqs_to_database(seq_map, config)
                seq_map = []
            if isoform_file:
                #Find all occasions of unusual major isoforms (major isoform number is not 1)
                unusual_major_isoforms = {}
                for stem in iso_map:
                    iso_numbers = iso_map[stem]
                    if not '1' in  iso_numbers:
                        continue
                    major_iso_number = None
                    for i in range(len(iso_numbers)):
                        if str(i+1) in iso_numbers:
                            continue
                        else:
                            major_iso_number = i+1
                    if major_iso_number is None:
                        unusual_major_isoforms[stem] = i+1
                    else:
                        unusual_major_isoforms[stem] = major_iso_number
            print(f'\nDatabase update of sequences done: {seq_file}\n')

    return unusual_major_isoforms

def main(config, fromScratch = False, update_mapping_db_keep_raw_files = False):
    #Step 1: Check if mapping SQL DB instance is there, create if not. Recreate for fromScratch mode
    fresh_instance = check_instance(config, fromScratch = fromScratch)

    print('Mapping DB instance checked, fresh_instance:', fresh_instance)

    if update_mapping_db_keep_raw_files:
        if config.container_version:
            raw_files_folder_path = '/structman/resources/'
        else:
            raw_files_folder_path = config.tmp_folder
    else:
        raw_files_folder_path = config.tmp_folder

    #Step 2: Check for raw files and download if necessary
    mapping_file_path, seq_file_paths = retrieve_raw_data(config, raw_files_folder_path, fromScratch = fromScratch)

    if update_mapping_db_keep_raw_files:
        print(f'\nDownloading all raw data files done. The files are stored in: {raw_files_folder_path}\n')
    else:
        print(f'\nDownloading all raw data files done. The files are temporarily stored in: {raw_files_folder_path}\n')

    #Step 3: Update the database
    #unusual_major_isoforms = insert_sequence_data(config, seq_file_paths, skip_inserts=True)

    unusual_major_isoforms = insert_sequence_data(config, seq_file_paths)

    ac_id_values = []
    ac_ref_values = []
    ac_ref_nt_values = []
    ac_ensembl_ids = []
    ensembl_acs = []

    max_values_at_a_time = int(100000 * config.gigs_of_ram)

    with gzip.open(mapping_file_path, 'rb') as f:
        previous_u_ac = None
        current_ensembl_ids = []
        for line in f:
            words = line.decode('ascii').split()
            if len(words) == 0:
                continue
            u_ac = words[0]
            id_name = words[1]
            id_value = words[2]

            if u_ac.count('-') == 1:
                stem, iso_number = u_ac.split('-')
                if stem in unusual_major_isoforms:
                    if int(iso_number) == unusual_major_isoforms[stem]:
                        u_ac = stem
                elif iso_number == '1':
                    u_ac = stem

            if id_name == 'UniProtKB-ID':
                ac_id_values.append((u_ac, id_value))

                if len(ac_id_values) == max_values_at_a_time:
                    update(config, 'UNIPROT', ['Uniprot_Ac', 'Uniprot_Id'], ac_id_values, mapping_db = True)
                    ac_id_values = []

            elif id_name == 'RefSeq':
                ac_ref_values.append((u_ac, id_value))

                if len(ac_ref_values) == max_values_at_a_time:
                    update(config, 'UNIPROT', ['Uniprot_Ac', 'RefSeq'], ac_ref_values, mapping_db = True)
                    ac_ref_values = []

            elif id_name == 'RefSeq_NT':
                ac_ref_nt_values.append((u_ac, id_value))

                if len(ac_ref_nt_values) == max_values_at_a_time:
                    update(config, 'UNIPROT', ['Uniprot_Ac', 'RefSeq_NT'], ac_ref_nt_values, mapping_db = True)
                    ac_ref_nt_values = []


            elif id_name == 'Ensembl_TRS' or id_name == 'Ensembl_PRO':
                id_value = id_value.split('.')[0]
                current_ensembl_ids.append(id_value)
                ensembl_acs.append((id_value, u_ac))

                if len(ensembl_acs) == max_values_at_a_time:
                    update(config, 'EMBL', ['EMBL_ID', 'Uniprot_Ac'], ensembl_acs, mapping_db = True)
                    ensembl_acs = []


            if u_ac != previous_u_ac:
                if len(current_ensembl_ids) > 0:
                    ac_ensembl_ids.append((u_ac, pack(current_ensembl_ids)))

                    current_ensembl_ids = []

                if len(ac_ensembl_ids) >= (max_values_at_a_time/5):
                    update(config, 'UNIPROT', ['Uniprot_Ac', 'Ensembl_Ids'], ac_ensembl_ids, mapping_db = True)
                    ac_ensembl_ids = []

            previous_u_ac = u_ac

    if not update_mapping_db_keep_raw_files:
        #Step 4: Remove the raw files
        os.remove(mapping_file_path)
        for seq_file in seq_file_paths:
            os.remove(seq_file)

        print('\nRemoving raw data files done.\n')

    print('\nParsing of mapping data file done.\n')

    if len(ac_id_values) > 0:
        update(config, 'UNIPROT', ['Uniprot_Ac', 'Uniprot_Id'], ac_id_values, mapping_db = True)
    print('\nDatabase update of Uniprot IDs done.\n')

    if len(ac_ref_values) > 0:
        update(config, 'UNIPROT', ['Uniprot_Ac', 'RefSeq'], ac_ref_values, mapping_db = True)
    print('\nDatabase update of RefSeq IDs done.\n')

    if len(ac_ref_nt_values) > 0:
        update(config, 'UNIPROT', ['Uniprot_Ac', 'RefSeq_NT'], ac_ref_nt_values, mapping_db = True)
    print('\nDatabase update of RefSeq NTs done.\n')


    if len(ensembl_acs) > 0:
        update(config, 'EMBL', ['EMBL_ID', 'Uniprot_Ac'], ensembl_acs, mapping_db = True)
    print('\nDatabase update of Ensembl IDs done.\n')


    if len(ac_ensembl_ids) > 0:
        update(config, 'UNIPROT', ['Uniprot_Ac', 'Ensembl_Ids'], ac_ensembl_ids, mapping_db = True)