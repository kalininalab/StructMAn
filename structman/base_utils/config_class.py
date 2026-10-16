import os
import random
import sys
import time
import configparser
import traceback
import multiprocessing
import subprocess
import sqlite3
import logging

import structman.base_utils.mysql_driver as MySQLdb
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
                 util_mode=False, output_util=False, external_call=True, verbosity=None, dbname = None, force_db_name = None,
                 print_all_errors=False, print_all_warns=False, restartlog=False, compute_ppi = True):
        self.prog_start_time = time.time()
        # read config file, auto add section header so old config files work
        self.config_parser_obj = configparser.ConfigParser()
        try:
            self.config_parser_obj.read(config_path)
        except configparser.MissingSectionHeaderError:
            with open(config_path, 'r') as f:
                self.config_parser_obj.read_string(f'[user]\n{f.read()}')
        try:
            cfg = self.config_parser_obj['user']
        except KeyError:
            print(f'Got KeyError in config_parser: {config_path=}')
            return None

        main_logger = logging.getLogger(__name__)

        self.logfile: str = f'{os.path.dirname(config_path)}/default_log.log'
        logging.basicConfig(filename=self.logfile, encoding='utf-8', level=logging.DEBUG)
        self.logger: logging.Logger = main_logger

        # fix typo `adress` -> `address`, but want old configs to still work, so also check for unfixed name
        old_address_name = cfg.get('db_adress', fallback='-')
        self.db_address = cfg.get('db_address', fallback=old_address_name)  # fixed address name takes precedence
        self.db_user_name = cfg.get('db_user_name', fallback='')
        self.db_password = cfg.get('db_password', fallback='')
        self.db_name = cfg.get('db_name', fallback='')

        self.mapping_db = cfg.get('mapping_db', fallback=None)

        #gets custom_db_path
        self.custom_db_path = cfg.get('custom_db_path', fallback=None)
        if self.custom_db_path is not None:
            self.custom_db_fasta_path = self.custom_db_path + "_custom_db.fasta"
            self.mmseqs2_custom_db_path = f'{self.custom_db_path }/custom_db_search_db_mmseqs2'

        if local_db is None:
            self.local_db_path = cfg.get('local_db_path', fallback = None)
        else:
            self.local_db_path = local_db

        self.ignore_local_mapping_db = False
        self.ray_local_mode = False
        self.python_env_expanded_for_ray = False
        # Set by the pipeline for inputs small enough that the ray cold start is not
        # worth paying before it is known whether there is any new work at all. While
        # it is True, parallel_db_access reports False and nothing starts a cluster;
        # core() clears it right after protCheck if the input brought new proteins.
        self.defer_ray_init = False
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
        # Pre-filter for the mmseqs search, as a fraction (mmseqs --min-seq-id).
        # Hits below option_seq_thresh identity are discarded after the global
        # alignment anyway, so anything mmseqs drops here saves parsing, an
        # annotation object and a full alignment. The default keeps a wide margin
        # below option_seq_thresh because mmseqs measures identity on its own
        # local alignment rather than on the global one used for the decision.
        # Set to 0 to restore the unfiltered behaviour.
        self.mmseqs_min_seq_id = cfg.getfloat('mmseqs_min_seq_id', fallback=0.25)
        if self.mmseqs_min_seq_id > 1.0:
            self.mmseqs_min_seq_id *= 0.01
        # Retrieval of the stored residue matrices is the single most memory hungry
        # step of the aggregation: every running select task holds the packed rows,
        # the unpacked per-structure arrays and the concatenated result at the same
        # time. These two knobs cap that peak independently of proc_n, which is far
        # too large a parallelism for this particular workload on big inputs.
        self.max_structures_per_select_package = cfg.getint('max_structures_per_select_package', fallback=20000)
        self.max_parallel_structure_selects = cfg.getint('max_parallel_structure_selects', fallback=20)

        # Residue data holds nothing that needs 64 bit floats: accessibilities,
        # angles, distances, small interaction counts and the integer codes for
        # SSA/aa/modres, all exact in float32. It is written as float32 since
        # Residue.encode_to_numpy() was narrowed, so this only still matters for
        # rows written before that change - retrieved residue matrices are the
        # bulk of the ray object store during aggregation, and narrowing those
        # legacy rows on read halves that footprint. Set to false to keep them
        # as they were stored.
        self.residue_data_float32 = cfg.getboolean('residue_data_float32', fallback=True)

        # Every remote task opens its own connections, so with a few hundred
        # workers the server hands out 'Too many connections' regularly. That is
        # transient, so it is worth waiting it out rather than failing the task.
        self.db_connection_trials = cfg.getint('db_connection_trials', fallback=12)
        self.db_connection_max_backoff = cfg.getfloat('db_connection_max_backoff', fallback=5.0)

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
        self.path_to_af_complex_db = cfg.get('path_to_af_complex_db', fallback=None)
        self.mmseqs2_db_path = cfg.get('mmseqs2_db_path', fallback='')
        self.mmseqs2_model_db_path = cfg.get('mmseqs2_model_db_path', fallback='')
        self.complex_model_db_path = cfg.get('complex_model_db_path', fallback='')
        self.mmseqs2_complex_model_db_path = cfg.get('mmseqs2_complex_model_db_path', fallback='')

        self.model_db_active = os.path.exists(self.mmseqs2_model_db_path)
        self.afdb = cfg.get('afdb', fallback=None)
        self.afdb_hetero = cfg.get('afdb_hetero', fallback=None)
        self.afdb_homo = cfg.get('afdb_homo', fallback=None)
        self.afdb_hetero_mmseqs_path = cfg.get('afdb_hetero_mmseqs_path', fallback=None)
        self.afdb_homo_mmseqs_path = cfg.get('afdb_homo_mmseqs_path', fallback=None)

        if not is_alphafold_db_valid(self) and self.afdb is None:
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
        # Per-protein alignment matrices. Also settable as --fetch_alignments; defined
        # here so the attribute exists for programmatically built Config objects too.
        self.fetch_alignments = cfg.getboolean('fetch_alignments', fallback=False)
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
        self.isolate_ppi = False

        # How many of a protein's structure annotations take part in the interface
        # calculation: get a backmap stored, and go into calculate_aggregated_interface_map.
        # Annotations of single-protein-chain complexes never do. Of the rest, up to
        # interface_annotations_keep_all are all taken; above that the keep_all best by
        # sequence identity plus every one above interface_annotations_min_seq_id; above
        # interface_annotations_max only the max best. See
        # interface.select_interface_annotations.
        self.interface_annotations_keep_all = cfg.getint('interface_annotations_keep_all', fallback=100)
        self.interface_annotations_max = cfg.getint('interface_annotations_max', fallback=1000)
        # A fraction, like the stored Sequence_Identity.
        self.interface_annotations_min_seq_id = cfg.getfloat('interface_annotations_min_seq_id', fallback=0.5)
        if self.interface_annotations_min_seq_id > 1.0:
            self.interface_annotations_min_seq_id *= 0.01

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

        # Memory ray is allowed to size itself against, in GiB. Leave at 0 to
        # autodetect (cgroup limit if there is one, else system memory). Set it
        # when the effective ceiling is not visible in the cgroup files - for
        # example a systemd MemoryHigh, which cgroup v1 does not expose.
        self.ray_memory_limit_gigs = cfg.getfloat('ray_memory_limit_gigs', fallback=0.0)

        # Share of that memory handed to ray's object store, before the
        # total_memory_quantile (0.9) is applied - so the default 0.5 means the
        # store gets 45% of the effective memory, unchanged from before this was
        # configurable. The rest has to cover the driver and every worker heap,
        # so a store sized close to the working set leaves the node with nothing
        # to breathe with: once the store fills and nothing in it is evictable,
        # puts block and the raylet starts missing its heartbeats.
        self.ray_object_store_fraction = cfg.getfloat('ray_object_store_fraction', fallback=0.5)

        # Ceiling on the atom count paraAnnotate packs into one annotation chunk.
        # Without it the chunk is total_cost // (8 * proc_n), which grows with the
        # input: a whole-PDB run asked for 2.08M atoms per chunk, and annotate holds
        # a chunk's whole structural_analysis_list in RAM before inserting it. 255 of
        # those at once, on top of the driver and the object store, left the cgroup
        # with no room to start the next worker process, and ray gave up with
        # "Failed to startup worker after retrying 5 times". Capped, per-worker peak
        # memory no longer depends on how big the input is. 0 disables the cap.
        self.annotation_chunk_atom_cap = cfg.getint('annotation_chunk_atom_cap', fallback=250000)

        # How often one annotation chunk may be resubmitted after ray failed to start
        # a worker for it. Only that failure is retried - it means the task never ran,
        # so nothing was half-inserted.
        self.annotation_retry_limit = cfg.getint('annotation_retry_limit', fallback=3)

        # How many ray worker processes may start at once. Most stages start a fresh process per
        # task (max_calls=1), from an environment on NFS, and Ray's default of one per CPU made
        # 255 start together; 512 of them missed their registration deadline in one chunk's
        # structural analysis. -1: proc_n // 4 (at least 8), the throttle para_parse_pdb already
        # used for its first burst. 0: Ray's default.
        self.ray_worker_startup_concurrency = cfg.getint('ray_worker_startup_concurrency', fallback=-1)

        # Seconds a new ray worker has to register before its start counts as failed (Ray's
        # default: 60), so a slow but healthy start under load is not thrown away. 0: Ray's default.
        self.ray_worker_register_timeout_seconds = cfg.getint('ray_worker_register_timeout_seconds', fallback=300)

        # Ray puts its session directory, its logs and - most importantly - its
        # spilled objects under /tmp by default. Object spilling on large inputs
        # writes hundreds of gigabytes, which /tmp usually cannot take; the
        # resulting ENOSPC surfaces as an unrelated-looking IO worker error.
        # Point this at a directory on a large filesystem.
        self.ray_temp_dir = cfg.get('ray_temp_dir', fallback=None)
        if self.ray_temp_dir:
            self.ray_temp_dir = resolve_path(self.ray_temp_dir)
            if not os.path.isdir(self.ray_temp_dir):
                try:
                    os.makedirs(self.ray_temp_dir, exist_ok=True)
                except OSError as err:
                    print(f'Could not create ray_temp_dir {self.ray_temp_dir}: {err}')
                    sys.exit()
            if not os.access(self.ray_temp_dir, os.W_OK):
                print(f'Need writing rights in ray_temp_dir {self.ray_temp_dir}')
                sys.exit()

        # Whether ray may spill objects to disk when the object store is full. Off by
        # default: spilling goes to the session directory under /tmp/ray, which a large
        # run fills until it dies with ENOSPC hours later, after the run has already
        # crawled at spill speed. With it off, ray_utils' watchdog stops the run the
        # moment the store overflows (spilled or fallback-allocated bytes appear), so
        # the problem surfaces immediately and with a message that names it.
        self.ray_allow_spilling = cfg.getboolean('ray_allow_spilling', fallback=False)

        # ray_init splits the effective memory limit into the object store's appetite,
        # what the driver holds, and the pool ray admits task memory= requests against.
        # The pool used to be a hardcoded half of the machine, which on a 1435 GiB node
        # capped it at 646 GiB while the store was using 12 of the 646 GiB reserved for
        # it. These two say what to subtract instead.
        #
        # The store allowance is a ceiling on what the store is expected to *hold*, not
        # its reservation: a whole-PDB chunk measured 72.6 GiB of residue blobs plus
        # 60 GiB of backmaps.
        self.ray_store_allowance_gigs = cfg.getfloat('ray_store_allowance_gigs', fallback=200.0)
        # What the driver keeps out of reach of the workers. Its Proteins container and
        # annotation graph reached ~93 GiB on a whole-PDB chunk.
        self.ray_driver_allowance_gigs = cfg.getfloat('ray_driver_allowance_gigs', fallback=120.0)
        # Both allowances are ceilings that scale down with the machine: taking 120 GiB
        # literally on a 100 GiB node subtracts more than exists.
        self.ray_driver_allowance_share = cfg.getfloat('ray_driver_allowance_share', fallback=0.15)
        # Never hand ray a pool smaller than this, or a task that reserves anything at
        # all could never be scheduled on a small or heavily reserved machine.
        self.ray_task_memory_floor_gigs = cfg.getfloat('ray_task_memory_floor_gigs', fallback=8.0)

        # Share of the object store that the retrieval of stored residue data may hold.
        # It is the first big consumer of the store, and a whole-PDB run filled it to
        # 99%, leaving nothing for getBackmaps and the classification packages that
        # follow - so the run died later, somewhere that said nothing about the cause.
        # Passing this stops the run in the stage that consumed the store. 0 disables
        # the check.
        self.ray_residue_store_share = cfg.getfloat('ray_residue_store_share', fallback=0.7)

        # Per classification task, the budget for decoded structures that proteins of
        # the same package can reuse (aggregation.DecodedStructureCache). Every task has
        # one, so the worst case is this times the number of cores.
        self.aggregation_decode_cache_mb = cfg.getint('aggregation_decode_cache_mb', fallback=256)

        # Memory each para_classify task reserves from ray's memory resource, per
        # million units of package cost (positions x annotations). Ray only starts a
        # task once its reservation fits, which is what bounds the stage: its memory is
        # set by the heaviest proteins, not by the chunk, so chunksize cannot.
        #
        # Calibrate from the uss= field of the [para_classify][worker] done lines, not
        # rss=. A worker maps the backmap and residue stores out of the object store and
        # its rss grows by tens of GB as it touches them, but those pages are shared -
        # 255 workers touching them cost the node one copy. The 2026-09-23 run showed
        # rss growth of ~50 GB per task correlating with elapsed time (+0.78) and not
        # with cost (-0.16), and an earlier 0.46 GiB per Mcost read off rss was what
        # left the reservation ~10x too small and OOMed the node at 255 tasks.
        # 0 disables the reservation.
        #
        # This wants re-measuring: classification_batch_size took the residue matrix off
        # the per-task peak, so the peak_uss the current 1.0 was fitted to is no longer
        # the peak the same package reaches. Re-fit it from peak_uss against cost on a
        # verbosity-4 run and lower it, or the reservation keeps holding back tasks that
        # would now fit.
        self.classification_gib_per_mcost = cfg.getfloat('classification_gib_per_mcost', fallback=1.0)

        # Backstop independent of that calibration: para_classify holds dispatch when
        # the cgroup is this full, and resumes as tasks return. Ray's own OOM killer
        # fires at 95% and takes the raylet with it, so this stays well below.
        self.classification_dispatch_memory_percent = cfg.getfloat(
            'classification_dispatch_memory_percent', fallback=78.0)

        # Annotations whose residue data one para_classify batch decodes. The whole of a
        # heavy whole-PDB protein reaches 8.6 GiB, which is what held classification to
        # ~70 of 255 cores; a batch of a thousand structures is about 150 MB. Lower it if
        # a task still peaks too high, raise it if the per-batch python overhead shows up
        # (the position loop runs once per batch).
        self.classification_batch_size = cfg.getint('classification_batch_size', fallback=1000)

        # Smallest reservation of a para_classify task, in GiB: the worker's fixed cost -
        # interpreter, imports, numba, the decode cache - which every task pays afresh, since
        # max_calls=1 starts a new worker for each. Fitted as the intercept of peak USS over
        # 9,022 tasks of a whole-PDB chunk (2026-09-24): 1909 MB. The earlier 1.0 was read off
        # a single fresh worker and left 61% of tasks above their estimate.
        self.classification_min_task_gib = cfg.getfloat('classification_min_task_gib', fallback=1.9)

        # What a para_classify task's memory is actually made of. cost is
        # positions x annotations, which drives how long a task runs but not how much it
        # decodes - measured peak_rss per 1M cost spanned 186 MB to 4748 MB, a 25x spread,
        # and 7JWB:D reserved 12.8 GB against a 49 GB peak while 9CI3:A reserved 37 GB to
        # use 7. build_local_stores decodes one residue matrix per annotated structure and
        # the interface entries of every chain of their complexes, so those two counts are
        # what the footprint follows.
        #
        # Both coefficients are fitted, with the intercept above, to the peak USS the workers
        # of that chunk reported against their counts: peak_uss = 1909 MB + 0.1915 MB x
        # annotations + 0.035 MB x interface chains, R^2 0.79 over 9,022 tasks. USS rather than
        # RSS, because RSS also counts the object-store pages a worker maps, and the node holds
        # one copy of those however many workers map them. No margin is added: the budget sums
        # tasks whose peaks do not coincide, and the cgroup backstop covers the per-task tail
        # (p95 of actual/fitted is 1.6).
        #
        # Per annotated structure, in MB: its decoded residue matrix.
        self.classification_mb_per_annotation = cfg.getfloat('classification_mb_per_annotation', fallback=0.19)

        # Per interface chain, in MB: one chain's combined_keys (int64) and scores
        # (float32) entries.
        self.classification_mb_per_interface_chain = cfg.getfloat(
            'classification_mb_per_interface_chain', fallback=0.035)

        # The coefficients above are the static prior. Every task reports its peak, and the driver
        # refits the model on (1, annotations, chains, cost) from the finished tasks - the batched
        # rewrite of para_classify added a cost-driven term the prior lacks - then scales the fit
        # by this quantile of peak / fitted (at least 1). The result only ever raises the prior.
        # Higher: fewer tasks at once, fewer OOM kills.
        self.classification_estimate_quantile = cfg.getfloat('classification_estimate_quantile', fallback=0.9)

        # What lowers it instead: no task is estimated at more than this margin over the highest
        # peak seen, once at least classification_ceiling_min_samples finished tasks that were
        # estimated above that peaked below it (ClassificationMemoryModel.ceiling). Since the
        # batched para_classify and the interface limit, peaks plateau - 3.8 GiB at most on a
        # whole-PDB chunk whose heaviest tasks were estimated at 15 GiB. 0 switches the cap off.
        self.classification_ceiling_margin = cfg.getfloat('classification_ceiling_margin', fallback=1.5)
        self.classification_ceiling_min_samples = cfg.getint('classification_ceiling_min_samples', fallback=20)

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

        # RAM-based per-task size, kept separate from self.chunksize so a user-supplied --chunksize
        # (which only controls how much gets read per outer batch) can't inflate the memory footprint
        # of any single concurrently-running Ray task.
        self.safe_task_chunksize = self.chunksize

        if not util_mode:
            if not external_call and not os.path.exists(self.outfolder):
                os.makedirs(self.outfolder)
            if self.verbosity >= 1:
                print(f'Using {self.proc_n} core(s), verbosity level: {self.verbosity}')
                # mysqlclient escapes bulk-insert parameters in C; pymysql does it
                # in pure Python. Surfaced here so it is visible which one a run used.
                print(f'MySQL driver: {MySQLdb.DRIVER_NAME}')
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
        self.max_connections = 40
        if not configure_mode:
            # Determine maximal package size from database
            try:
                db, cursor = self.getDB(server_connection=True)
            except:

                db = None
                self.main_db_is_set = False
                self.max_package_size = None
                
            if db is not None and not basic_util_mode:
                if dbname is not None:
                    if not force_db_name:
                        if dbname[:len(self.db_user_name)] != self.db_user_name:
                            dbname = f'{self.db_user_name}_{dbname}'
                    self.db_name = dbname

                if self.db_address != '-':
                    try:
                        cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_allowed_packet'")
                        self.max_package_size = int(cursor.fetchone()[1]) * 100 // 99
                        cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_connections'")
                        self.max_connections = int(cursor.fetchone()[1])
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
                cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_connections'")
                self.max_connections = int(cursor.fetchone()[1])
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

    @property
    def ppi_data_available(self) -> bool:
        """Whether the PPI tables can actually hold data for this run.

        compute_ppi alone is not enough to decide that. fast_pdb_annotation skips
        getBackmaps, so aggregation skips the whole interface block that fills
        Interface, RS_Position_Interface and Protein_Protein_Interaction_N (see
        aggregation.para_classify). The output stage used to test compute_ppi on its
        own and then select those tables anyway - against position id sets that in a
        whole-PDB run are hundreds of millions of ids wide, for results that could
        only ever be empty.

        compute_ppi cannot simply be cleared for fast runs instead: it is forced back
        on by repairDB.check_structman_version whenever the database reports the PPI
        feature, and that runs during __init__, before structman_main gets to apply
        fast_pdb_annotation. So the two flags are combined here, at use time.
        """
        return self.compute_ppi and not self.fast_pdb_annotation

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

        if mapping_db and not server_connection:
            if self.ignore_local_mapping_db:
                return None, None
            if self.mapping_db is None:
                return None, None

        # 'Too many connections' is the normal answer when a few hundred ray
        # workers reach for the database at the same moment, and it is purely
        # transient. The old handler retried 10 times 0.1s apart - a whole second -
        # and then handed back (None, None), which every caller then used as a
        # cursor and died with a misleading AttributeError. Back off properly
        # instead, with jitter so the retries of all workers do not line up again.
        attempt = try_again
        last_exception_text = None
        while attempt <= self.db_connection_trials:
            try:
                if server_connection:
                    db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, autocommit=autocommit)
                elif mapping_db:
                    db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, database=self.mapping_db, autocommit=autocommit)
                elif db_name is None:
                    db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, database=self.db_name, autocommit=autocommit)
                else:
                    db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, database=db_name, autocommit=autocommit)

                cursor = db.cursor()
                if cursor is not None:
                    return db, cursor

                target_db = db_name if db_name is not None else self.db_name
                self.errorlog.add_warning(f'getDB misfunctioned: {self.db_address}, {self.db_user_name}, {target_db}')
                try:
                    db.close()
                except Exception:
                    pass
            except Exception as err:
                if silent:  # In database create we test the connection and expect it to fail
                    return None, None
                last_exception_text = f'{type(err).__name__}: {err}'

            attempt += 1
            if attempt > self.db_connection_trials:
                break
            # 0.2s, 0.4s, 0.8s ... capped, each with up to 50% extra jitter.
            backoff = min(0.1 * (2 ** attempt), self.db_connection_max_backoff)
            time.sleep(backoff * (1.0 + random.random() * 0.5))

        if silent:
            return None, None

        if last_exception_text is None:
            last_exception_text = 'connection succeeded but no cursor could be opened'
        target_db = db_name if db_name is not None else self.db_name
        self.errorlog.add_error(
            f'Couldnt get database connection after {self.db_connection_trials} trials: '
            f'{self.db_address}, {self.db_user_name}, {target_db}\n{last_exception_text}\n{traceback.format_exc()}'
        )
        return None, None
