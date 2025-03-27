import os
import sys
import time
import configparser
import traceback
import multiprocessing
import subprocess
import sqlite3

import pymysql as MySQLdb
import ray

from structman.base_utils.base_utils import resolve_path, is_alphafold_db_valid, Errorlog
from structman.lib import repairDB
from structman.settings import (
    PDB_SYNC_SCRIPT,
    STRUCTMAN_DB_SQL,
    STRUCTMAN_DB_SQLITE,
    UNIPROT_DB_SQL,
    RINERATOR_DIR,
    RINERATOR_PATH,
    REDUCE_HET_DICT,
    RESOURCES_DIR
    )

class Config:
    def __init__(self, config_path, num_of_cores=1, output_path='', basic_util_mode=False, configure_mode=False, local_db = None, db_mode = None,
                 util_mode=False, output_util=False, external_call=True, verbosity=None, dbname = None,
                 print_all_errors=False, print_all_warns=False, restartlog=False, compute_ppi = True):
        self.prog_start_time = time.time()
        # read config file, auto add section header so old config files work
        self.config_parser_obj = configparser.ConfigParser()
        try:
            self.config_parser_obj.read(config_path)
        except configparser.MissingSectionHeaderError:
            with open(config_path, 'r') as f:
                self.config_parser_obj.read_string(f'[user]\n{f.read()}')
        cfg = self.config_parser_obj['user']

        # fix typo `adress` -> `address`, but want old configs to still work, so also check for unfixed name
        old_address_name = cfg.get('db_adress', fallback='-')
        self.db_address = cfg.get('db_address', fallback=old_address_name)  # fixed address name takes precedence
        self.db_user_name = cfg.get('db_user_name', fallback='')
        self.db_password = cfg.get('db_password', fallback='')
        self.db_name = cfg.get('db_name', fallback='')

        self.mapping_db = cfg.get('mapping_db', fallback=None)

        #gets custom_db_path
        self.custom_db_path = cfg.get('custom_db_path', fallback='')
        self.custom_db_fasta_path = self.custom_db_path + "_custom_db.fasta"
        self.mmseqs2_custom_db_path = f'{self.custom_db_path }/custom_db_search_db_mmseqs2'

        if local_db is None:
            self.local_db_path = cfg.get('local_db_path', fallback = None)
        else:
            self.local_db_path = local_db

        self.ignore_local_mapping_db = False
        self.ray_local_mode = False
        self.python_env_expanded_for_ray = False
        self.user_mail = cfg.get('mail', fallback='')

        self.outfolder = resolve_path(output_path)
        self.num_of_cores = num_of_cores
        self.indels_given_by_input = True
        self.only_wt = False
        self.rare_residues = set()

        self.model_indel_structures = cfg.getboolean('model_indel_structures', fallback=False)
        self.force_modelling = cfg.getboolean('force_modelling', fallback = False)

        self.verbosity = cfg.getint('verbosity', fallback=1)

        self.test_low_mem_system = cfg.getboolean('test_low_mem_system', fallback=False)

        self.resources = cfg.get('resources', fallback='manu')
        self.container_version = cfg.getboolean('container_version', fallback = False)

        self.proc_n = 48

        # active options
        self.n_of_chain_thresh = cfg.getint('n_of_chain_thresh', fallback=15)  # Structures with n_of_chain_thresh or more chains get a nested paralellization
        self.option_seq_thresh = cfg.getfloat('seq_thresh', fallback=35.0)
        if self.option_seq_thresh <= 1.0:
            self.option_seq_thresh *= 100.0
        self.option_res_thresh = cfg.getfloat('res_thresh', fallback=4.5)
        self.option_ral_thresh = cfg.getfloat('cov_thresh', fallback=0.5)
        if self.option_ral_thresh > 1.0:
            self.option_ral_thresh *= 0.01

        self.tax_id = cfg.get('tax_id', fallback=None)
        self.ref_genome_id = cfg.get('ref', fallback=None)
        self.mrna_fasta = cfg.get('mrna', fallback=None)

        self.surface_threshold = cfg.getfloat('surface_threshold', fallback=0.16)
        self.buried_threshold = cfg.getfloat('buried_threshold', fallback=0.05)
        self.short_distance_threshold = cfg.getfloat('lig_short_dist_thresh', fallback=5.0)
        self.long_distance_threshold = cfg.getfloat('lig_long_dist_thresh', fallback=8.0)
        self.ligand_interest_sphere = cfg.getfloat('ligand_interest_sphere', fallback=25.0)

        # Paths to Helper Programs
        self.mmseqs2_path = cfg.get('mmseqs2_path', fallback='')
        self.output_path = cfg.get('output_path', fallback='')
        self.pdb_path = cfg.get('pdb_path', fallback='')
        self.path_to_model_db = cfg.get('path_to_model_db', fallback='')
        self.mmseqs2_db_path = cfg.get('mmseqs2_db_path', fallback='')
        self.mmseqs2_model_db_path = cfg.get('mmseqs2_model_db_path', fallback='')

        self.model_db_active = os.path.exists(self.mmseqs2_model_db_path)
        if not is_alphafold_db_valid(self):
            self.model_db_active = False

        self.read_only_mode = cfg.getboolean('read_only_mode', fallback=False)
        self.read_hybrid_mode = cfg.getboolean('read_hybrid_mode', fallback=False)
        self.fast_pdb_annotation = cfg.getboolean('fast_pdb_annotation', fallback=False)

        self.annovar_path = cfg.get('annovar_path', fallback='')
        self.dssp_path = cfg.get('dssp_path', fallback = '')
        self.rin_db_path = cfg.get('rin_db_path', fallback='')

        self.errorlog_path = cfg.get('errorlog_path', fallback=None)
        self.warnlog_path = cfg.get('warnlog_path', fallback=None)

        self.go = cfg.getboolean('do_goterm', fallback=False)
        self.anno = cfg.getboolean('do_anno', fallback=False)
        self.classification = cfg.getboolean('do_classification', fallback=True)
        self.gene = cfg.getboolean('do_genesort', fallback=False)
        self.path = cfg.getboolean('do_pathway', fallback=False)
        self.godiff = cfg.getboolean('do_godiff', fallback=False)
        self.pathdiff = cfg.getboolean('do_pathdiff', fallback=False)
        self.do_modelling = cfg.getboolean('do_modelling', fallback=False)
        self.multi_modelling = cfg.getboolean('multi_modelling', fallback=False)
        self.ligand_file = cfg.get('ligand_file', fallback=None)
        self.mod_per_mut = cfg.getint('mod_per_mut', fallback=0)
        self.mod_per_gene = cfg.getint('mod_per_gene', fallback=0)
        self.tanimoto_cutoff = cfg.getfloat('tanimoto_cutoff', fallback=0.05)
        self.milieu_threshold = cfg.getfloat('milieu_threshold', fallback=10.0)
        self.intra_milieu_threshold = cfg.getfloat('intra_milieu_threshold', fallback=15.0)
        self.proteome = cfg.getboolean('proteome', fallback=False)
        self.intertable_conf = cfg.getboolean('intertable', fallback=False)
        self.overwrite_incorrect_wt_aa = cfg.getboolean('overwrite_incorrect_wt_aa', fallback=False)

        self.microminer_exe = cfg.get('microminer_exe', fallback = '')
        self.microminer_db = cfg.get('microminer_db', fallback = '')
        self.microminer_identity = cfg.getfloat('microminer_identity', fallback = 0.6)
        self.microminer_kmer_matching_rate = cfg.getfloat('microminer_kmer_matching_rate', fallback = 0.4)
        self.microminer_fragment_length = cfg.getint('microminer_fragment_length', fallback=8)
        self.microminer_k = cfg.getint('microminer_k', fallback=7)
        self.update_microminer_active = cfg.getboolean('update_microminer_active', fallback=False)

        self.condition_1_tag = None
        self.condition_2_tag = None

        self.pdb_sync_script = cfg.get('pdb_sync_script', fallback=PDB_SYNC_SCRIPT)

        self.database_source_path = cfg.get('database_source_path', fallback=STRUCTMAN_DB_SQL)
        self.sqlite_database_source_path = cfg.get('sqlite_database_source_path', fallback=STRUCTMAN_DB_SQLITE)        
        self.mapping_db_source_path = cfg.get('mapping_db_source_path', fallback=UNIPROT_DB_SQL)
        self.ensembl_server = "https://rest.ensembl.org"

        self.rinerator_base_path = RINERATOR_DIR
        self.rinerator_path = RINERATOR_PATH
        os.environ["REDUCE_HET_DICT"] = REDUCE_HET_DICT

        self.mmseqs_tmp_folder = cfg.get('mmseqs_tmp_folder')
        self.tmp_folder = self.mmseqs_tmp_folder

        self.fasta_input = cfg.getboolean('fasta_input', fallback=False)

        self.compute_ppi = compute_ppi

        if verbosity is not None:
            self.verbosity = verbosity

        # Checking whether the given paths in config exist or not if it is given and not exist, system gives error message and exits

        if self.mmseqs2_path:
            if self.mmseqs2_path == 'mmseqs':
                cmd = "mmseqs"
                ret_val = subprocess.getstatusoutput(cmd)[0]
                #ret_val = os.system(cmd)
                if ret_val != 0:
                    print('MMSEQS2 path is not correct, please correct the path')
                    sys.exit()
            else:
                isExist = os.path.exists(self.mmseqs2_path)
                if not isExist:
                    print('MMSEQS2 path does not exist, please check the path')
                    sys.exit()

        if self.pdb_path:
            isExist = os.path.exists(self.pdb_path)
            if not isExist:
                print('PDB path does not exist, please check the path')
                #sys.exit()

        
        if self.rin_db_path:
            isExist = os.path.exists(self.rin_db_path)
            if not isExist:
                print('RIN DB path does not exist, please check the path')
                #sys.exit()


        if self.mmseqs_tmp_folder:
            isExist = os.path.exists(self.mmseqs_tmp_folder)
            if not isExist:
                print('MMSEQS temp folder path does not exist, please check the path')
                sys.exit()
            if not os.access(self.mmseqs_tmp_folder, os.R_OK):
                print('Need writing rights in the MMSEQS temp folder, please check the path')
                sys.exit()

        self.smiles_path = f'{RESOURCES_DIR}/Components-smiles-stereo-oe.smi'
        self.inchi_path = f'{RESOURCES_DIR}/inchi_base.tsv'

        if self.resources == 'auto' and self.num_of_cores is None:
            self.proc_n = multiprocessing.cpu_count() - 1
            if self.proc_n <= 0:
                self.proc_n = 1

        if self.num_of_cores is not None:
            self.proc_n = self.num_of_cores

        if self.proc_n > multiprocessing.cpu_count():
            if self.verbosity >= 1:
                print('More processes annotated (', self.proc_n, ') than cores registered in system (', multiprocessing.cpu_count(), ').')
            self.proc_n = multiprocessing.cpu_count()

        
        self.gigs_of_ram = ray._private.utils.get_system_memory() / 1024 / 1024 / 1024
        self.low_mem_system = self.gigs_of_ram < 40 or self.test_low_mem_system  # Less than 20Gb is a low memory system
        if self.test_low_mem_system:
            self.gigs_of_ram = 8
        if self.low_mem_system:
            self.chunksize = int(min([500, max([((self.gigs_of_ram * 40) // self.proc_n) - 120, 60 // self.proc_n, 1])]))
        else:
            self.chunksize = int(min([1500, max([((self.gigs_of_ram * 100) // self.proc_n) - 60, 120 // self.proc_n, 1])]))

        if not util_mode:
            if not external_call and not os.path.exists(self.outfolder):
                os.makedirs(self.outfolder)
            if self.verbosity >= 1:
                print(f'Using {self.proc_n} core(s), verbosity level: {self.verbosity}')
            if (not external_call) and (not print_all_warns):  # no need for errorlogs, when the config is generated not from the main script
                self.errorlog_path = os.path.join(self.outfolder, 'errorlogs', 'errorlog.txt')
                self.warnlog_path = os.path.join(self.outfolder, 'errorlogs', 'warnlog.txt')
                errorlog_dir = os.path.dirname(self.errorlog_path)
                if not os.path.exists(errorlog_dir):
                    os.makedirs(errorlog_dir)
                if restartlog and os.path.isfile(self.errorlog_path):
                    try:
                        os.remove(self.errorlog_path)
                        os.remove(self.warnlog_path)
                    except:
                        pass

        self.errorlog = Errorlog(path=self.errorlog_path, warn_path=self.warnlog_path, print_all_errors=print_all_errors, print_all_warns=print_all_warns)

        sqlite_max_package_size = 800000000
        if not configure_mode:
            # Determine maximal package size from database
            try:
                db, cursor = self.getDB(server_connection=True)
            except:

                db = None
                self.main_db_is_set = False
                self.max_package_size = None

            if db is not None and not basic_util_mode:
                if self.db_address != '-':
                    try:
                        cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_allowed_packet'")
                        self.max_package_size = int(cursor.fetchone()[1]) * 100 // 99
                        db.close()
                        db, cursor = self.getDB()
                        self.main_db_is_set = True
                        repairDB.check_structman_version(self)
                    except:
                        if self.verbosity >= 2:
                            print(f'Database instance ({self.db_address}:{self.db_name}) not found in config initialization. Trying to create the instance ...')
                        if self.verbosity >= 3:
                            [e, f, g] = sys.exc_info()
                            g = traceback.format_exc()
                            print('\n'.join([str(e), str(f), str(g)]))
                        repairDB.load(self)
                        try:
                            db, cursor = self.getDB()
                            self.main_db_is_set = True
                        except:
                            
                            self.main_db_is_set = False
                            self.max_package_size = None
                    try:
                        db.close()
                    except:
                        pass
                else:
                    self.main_db_is_set = True
                    self.max_package_size = sqlite_max_package_size
                    repairDB.check_structman_version(self)
                    db.close()
            elif db is not None and basic_util_mode:
                cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_allowed_packet'")
                self.max_package_size = int(cursor.fetchone()[1]) * 100 // 99
                db.close()
            elif self.db_address == '-':
                self.main_db_is_set = True
                self.max_package_size = sqlite_max_package_size


        else:
            self.main_db_is_set = False
            self.max_package_size = None
            db = None

        # Check for mapping DB instanceB
        if not configure_mode:
            self.check_mapping_db()

        if self.local_db_path is not None:
            if not os.path.exists(self.local_db_path) and db_mode != 'create':
                print(f"Couldn't find provided local structman database: {self.local_db_path}")

        if db is not None:
            try:
                db.close()
            except:
                pass

    def check_mapping_db(self):
        # Check for mapping DB instance
        try:
            db, cursor = self.getDB(mapping_db = True)
            db.close()
            self.mapping_db_is_set = True
        except:
            self.mapping_db_is_set = False

    def getDB(self, server_connection: bool = False, mapping_db: bool = False, try_again: int = 0, silent: bool = False, db_name: str | None = None, autocommit: bool = False):
        if self.db_address == '-':
            if mapping_db or server_connection:
                return None, None
            if self.local_db_path is None:
                self.local_db_path = os.path.abspath('local_structman_database.db')
            if not os.path.exists(self.local_db_path):
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                print(f'Database server connection not configured and provided local database path not found. Creating a new local database here:\n{self.local_db_path}\n{e}\n{f}\n{g}')
                db = sqlite3.connect(self.local_db_path)
            else:
                db = sqlite3.connect(self.local_db_path)
            cursor = db.cursor()
            return db, cursor

        try:
            if server_connection:
                db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, autocommit=autocommit)
                cursor = db.cursor()
            elif mapping_db:
                if self.ignore_local_mapping_db:
                    return None, None

                if self.mapping_db is None:
                    return None, None
                db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, database=self.mapping_db, autocommit=autocommit)
                cursor = db.cursor()
            elif db_name is None:
                db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, database=self.db_name, autocommit=autocommit)
                cursor = db.cursor()
                if cursor is None:
                    self.errorlog.add_warning(f'getDB misfunctioned: {self.db_address}, {self.db_user_name}, {self.db_password}, {self.db_name}')
            else:
                db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, database=db_name, autocommit=autocommit)
                cursor = db.cursor()
                if cursor is None:
                    self.errorlog.add_warning(f'getDB misfunctioned: {self.db_address}, {self.db_user_name}, {self.db_password}, {db_name}')
            return db, cursor
        except:
            if silent: #In database create we test the connection and expect it to fail
                return None, None
            if try_again < 10:
                time.sleep(0.1)
                return self.getDB(server_connection=server_connection, mapping_db=mapping_db, autocommit=autocommit, try_again = (try_again + 1))
            else:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                self.errorlog.add_error(f'Couldnt get database connection: {self.db_address}, {self.db_user_name}, {self.db_password}, {self.db_name}\n{e}\n{f}\n{g}')
                return None, None
