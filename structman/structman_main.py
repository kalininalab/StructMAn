#!/usr/bin/python3
# This is the wrapper script, including and automatizing all functions StructMAn has to offer

import getopt
import os
import resource
import sys
import time

import pymysql as MySQLdb
import ray

# if running as script, add local structman package to path
if __name__ == "__main__":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__)))))

from structman.lib import repairDB, serializedPipeline, searchLargeStructures
from structman.lib.output import output
from structman.lib.output import ppi as ppi_output
from structman.lib.output import massmodel
from structman.lib.output import annotation
from structman.lib.output import genes as genes_out_package
from structman.lib.database import database
from structman.scripts import update
from structman.base_utils.base_utils import resolve_path, set_modeller_key
from structman.base_utils.config_class import Config

# Check if autopep8 is installed for developers. Disable warning for end users by utilizing conda environment variable `STRUCTMAN_WARN`
try:
    if os.environ.get('STRUCTMAN_WARN', 'yes') != 'no':
        import autopep8
except ImportError:
    print("Autopep8 is not installed. For developers, please set up autopep8 git hooks according to the internal wiki")
    user_input = input("If you are not a developer, please input `x` to disable this message and re-activate your conda environment.\n>>> ")
    if user_input.lower().startswith('x'):
        env_path = os.path.join(os.environ['CONDA_PREFIX'], 'etc/conda/')
        os.makedirs(os.path.join(env_path, 'activate.d'), exist_ok=True)
        os.makedirs(os.path.join(env_path, 'deactivate.d'), exist_ok=True)
        with open(os.path.join(env_path, 'activate.d', 'env_vars.sh'), 'a') as file:
            file.write('export STRUCTMAN_WARN=no\n')
        with open(os.path.join(env_path, 'deactivate.d', 'env_vars.sh'), 'a') as file:
            file.write('unset STRUCTMAN_WARN\n')



def limit_memory(soft_limit):
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (soft_limit, hard))

def main(infiles, out_folder, config):
    # resource.setrlimit(resource.RLIMIT_AS,(6442450944,8589934592))
    # resource.setrlimit(resource.RLIMIT_AS,(442450944,589934592))
    for infile in infiles:
        if isinstance(infile, str):
            if config.verbosity >= 1:
                print('Processing file: ', infile)
            # check if the infile got already processed
            session_id = database.getSessionId(infile, config)


            # create the output folder
            trunk, infilename = os.path.split(infile)

            if infilename[-6:] == '.fasta':
                config.fasta_input = True

            if out_folder is None:
                out_folder = resolve_path(os.path.join(trunk, "Output"))
        else:
            config.errorlog.add_error('Invalid input file path: {infile}')

        if not os.path.exists(out_folder):
            os.makedirs(out_folder)


        # run the main pipeline
        
        if session_id is None:
            session_id, config, main_returncode = serializedPipeline.main(infile, config)
        else:
            main_returncode = 0


        # run the output scripts
        session_name = infilename.rsplit('.', 1)[0]
        months = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun', 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
        time_struct = time.gmtime()
        year = str(time_struct[0])
        month = months[time_struct[1]]
        day = str(time_struct[2])
        date = f"{day}{month}{year}"
        outpath = os.path.join(out_folder, session_name)

        if not os.path.exists(outpath):
            os.makedirs(outpath)
        
        out_returncode = output.main(session_id, outpath, config)

        returncode = main_returncode | out_returncode

        if config.read_only_mode:
            repairDB.remove_sessions(config)
        elif config.read_hybrid_mode:
            repairDB.clear(config)

    return returncode

# needed for entry point in setup.py
def structman_cli():
    main_file_path = resolve_path(__file__)

    disclaimer = ''.join([
                
        ' .d8888b.  888                              888    888b     d888        d8888           .d8888b.       .d8888b. \n',
        'd88P  Y88b 888                              888    8888b   d8888       d88888          d88P  Y88b     d88P  Y88b\n', 
        'Y88b.      888                              888    88888b.d88888      d88P888                 888     888    888\n', 
        ' "Y888b.   888888 888d888 888  888  .d8888b 888888 888Y88888P888     d88P 888 88888b.       .d88P     888    888\n',
        '    "Y88b. 888    888P"   888  888 d88P"    888    888 Y888P 888    d88P  888 888 "88b  .od888P"      888    888\n', 
        '      "888 888    888     888  888 888      888    888  Y8P  888   d88P   888 888  888 d88P"          888    888\n',
        'Y88b  d88P Y88b.  888     Y88b 888 Y88b.    Y88b.  888   "   888  d8888888888 888  888 888"       d8b Y88b  d88P\n', 
        ' "Y8888P"   "Y888 888      "Y88888  "Y8888P  "Y888 888       888 d88P     888 888  888 888888888  Y8P  "Y8888P" \n',
        '\n',                                                                                                     
        '                    ======================================\n',
        '                    |                                    |\n',                                                                                                                                    
        '                    |  Usage: structman <-i input_file>  |\n',
        '                    |                                    |\n',           
        '                    ======================================\n\n',
        'more functionalities can be used giving a second key word:\n',
        'structman database    gives you more info about the database utility functions\n',
        'structman config      gives you more info about the config utility functions\n',
        'structman update      gives you more info about the different update options\n',
        'structman dev         gives you more info about options for developers\n\n\n',

        '----------------------------------------------------------------------------------------------------------------\n',
        '| Optional parameter  | Default value               | Description                                              |\n',
        '|---------------------|-----------------------------|----------------------------------------------------------|\n',
        '| <-n threads>        | All available -1            | Number of cores to be used                               |\n',
        '|---------------------|-----------------------------|----------------------------------------------------------|\n',
        '| <-o output_folder>  | cwd/Output/infilename       | Path to the output directory                             |\n',
        '|---------------------|-----------------------------|----------------------------------------------------------|\n',
        '| <--verbosity> [0-8] | 1                           | verbosity of command line messages                       |\n',
        '|                     |                             | 0: silence; 1: minimal progress reports;                 |\n',
        '|                     |                             | 2: more detailed progress reports;                       |\n',
        '|                     |                             | 3: 3 and greater are for debugging, more time stamps;    |\n',
        '|                     |                             | 4: subprocesses printing, best redirect into a log file; |\n',
        '|                     |                             | 5-8: more and more debug prints                          |\n',  
        '---------------------------------------------------------------------------------------------------------------|\n'
    ])

    update_util_disclaimer = ''.join([
        'Usage: structman.py update [commands] <-c config_file> <-p path_to_local_pdb>\n',
        'The update functionalities manage local database, that significantly increases the performance of StructMAn.\n',
        'They need a lot of disk space though. Thus the usage of local databases is only recommended when one wants to process larger inputs with StructMAn.\n'
        'Only needs the path to the local instance of the PDB (-p), if it is not already specified in the config file.\n\n'
        '#### Commands: ####\n\n',

        'pdb                       uses rsync to update all files of the local PDB that could potentially used by StructMAn.\n',
        '                          If the given path to the local PDB is empty, a large portion of the PDB will be downloaded,\n',
        '                          be sure that there is enought disk space available.\n\n',

        'rindb                     calculates and stores the RIN for each PDB structure in the local PDB.\n',
        '                          When called for the first time, this can take multiple hours to complete.\n',
        '                          Also requires a large amount of disk space.\n\n',

        'rindb_from_scratch        same as rindb, but overwrites all RINs currently in the RINdb.\n\n',

        'mapping_db                downloads protein identifier mapping information and sequences from UNIPROT.\n',
        '                          Builds a SQL database and removes the downloaded files again.\n',
        '                          Requires a configured SQL database server connection.\n',
        '                          Requires around 100Gb of temporal disk space and 35Gb of space on the SQL server.\n\n',

        'mapping_db_from_scratch   same as mapping_db but removes an existing instance of the mapping DB.\n\n'

        'alphafold_db              downloads all available protein structures models created by alphafold stored at: https://alphafold.ebi.ac.uk/\n'
    ])

    config_util_disclaimer = ''.join([
        'Usage: structman.py config [command] [value] <-c config_file>\n',
        'The config commands enable the expansion of StructMAn by giving it access to additional functionalities, which are too large to be included by default or require an external license.\n\n',
        '#### Commands: ####\n\n',
        'set_local_pdb_path <path_to_local_pdb> :               enables the usage of a local instance of the PDB. If this is not available, StructMAn has to download all structural information from the web.\n\n',
        '<any config variable name> <any value> :               modifies any variable in the config file. Warning: does not check if given input makes sense.',
    ])

    database_util_disclaimer = ''.join([
        'Usage: structman.py database [command] <-c config_file>\n\n',
        '#### Commands: ####\n\n',
        'reset :\n  deletes all content of the database\n\n',
        'clear :\n  deletes all content of the database, but keeps the results of all stored structures\n\n',
        'export <-p path_to_a_folder> :\n  exports the database to a .sql.gz format file\n\n',
        'destroy :\n  completely removes the database (database name is taken from config file)\n\n',
        'create  :\n  creates an instance of the database (database name is taken from config file)'
    ])

    dev_util_disclaimer = ''.join([
        'Hi, are you a dev? If yes, can you please write the disclaimer for the dev options? Just look in what miserable state they are currently in:\n'

        '| <-c config_file>    | standard installation:  | Path to the configuration file                               |\n',
        '|                     | config.txt              |                                                              |\n',
        '|                     | container installation: | create custom config.txt inside mounted                      |\n',
        '|                     | not mounted             | folder, or use <structman config>                            |\n',
        '|---------------------|-------------------------|-------------------------------------------                   |\n',
        '|---------------------|-------------------------|-------------------------------------------                   |\n',
        '| <--printerrors>     | False                   | prints error messages in the console                         |\n',
        '|                     |                         | instead of logging them. Warning messages                    |\n',
        '|                     |                         | still go into the log                                        |\n',
        '|---------------------|-------------------------|-------------------------------------------                   |\n',
        '| <--printwarnings>   | False                   | prints error and warning messages in                         |\n',
        '|                     |                         | the console instead of logging them                          |\n',
        '|---------------------|-------------------------|-------------------------------------------                   |\n',
        '| <--restartlog>      | False                   | wipes the log before starting the session                    |\n',
        '|---------------------|-----------------------------|----------------------------------------------------------|\n',
        '| <--custom_db>       | False                       | custom database usage                                    |\n',
    ])

    argv = sys.argv[1:]

    if len(argv) == 0:
        print(disclaimer)
        sys.exit(1)

    if argv[0] == 'dev':
        print(dev_util_disclaimer)
        sys.exit(0)

    database_util = False
    if argv[0] == 'database':
        database_util = True
        argv = argv[1:]

        possible_key_words = set(['reset', 'out', 'create', 'destroy', 'clear', 'export', 'reduce', 'remove_sessions', 'fix', 'clone'])

        if len(argv) == 0 or argv[0] == '-h' or argv[0] == '--help':
            print(database_util_disclaimer)
            sys.exit(1)
        if argv[0] in possible_key_words:
            db_mode = argv[0]
            argv = argv[1:]
        else:
            print(database_util_disclaimer)
            sys.exit(1)
    else:
        db_mode = None

    update_util = False
    update_pdb = False
    update_rindb = False
    update_rindb_from_scratch = False
    update_mapping_db = False
    update_mapping_db_from_scratch = False
    update_mapping_db_keep_raw_files = False
    update_alphafold_db = False
    structure_limiter = None
    update_microminer = False
    microminer_from_scratch = False
    check_search_db = False

    if len(argv) > 0:
        if argv[0] == 'update':
            update_util = True
            argv = argv[1:]

            if len(argv) == 0 or argv[0] == '-h' or argv[0] == '--help':
                print(update_util_disclaimer)
                sys.exit(1)
            if 'pdb' in argv:
                update_pdb = True
            if 'rindb' in argv:
                update_rindb = True
            if 'rindb_from_scratch' in argv:
                update_rindb = True
                update_rindb_from_scratch = True
            if 'mapping_db' in argv:
                update_mapping_db = True
            if 'mapping_db_from_scratch' in argv:
                update_mapping_db = True
                update_mapping_db_from_scratch = True
            if 'mapping_db_keep_raw_files' in argv:
                update_mapping_db = True
                update_mapping_db_keep_raw_files = True
            if 'mapping_db_from_scratch_keep_raw' in argv:
                update_mapping_db = True
                update_mapping_db_from_scratch = True
                update_mapping_db_keep_raw_files = True
            if 'alphafold_db' in argv:
                update_alphafold_db = True
            if 'mmdb' in argv:
                update_microminer = True
            if 'mmdb_from_scratch' in argv:
                update_microminer = True
                microminer_from_scratch = True

            if 'check_search_db' in argv:
                check_search_db = True

            if not (update_pdb or update_rindb or update_mapping_db or update_alphafold_db or update_microminer or check_search_db):
                print(update_util_disclaimer)
                sys.exit(1)

    configure_mode = False
    configure_help_mode = False
    if len(argv) > 0:
        if argv[0] == 'config':
            configure_mode = True
            argv = argv[1:]

            if len(argv) == 0 or argv[0] == '-h' or argv[0] == '--help':
                configure_help_mode = True
            else:
                conf_update_pdb_path = None

                config_variable_name = None

                if argv[0] == 'set_local_pdb_path':
                    if len(argv) == 1:
                        print(config_util_disclaimer)
                        sys.exit(1)
                    conf_update_pdb_path = argv[1]
                    if not os.path.exists(conf_update_pdb_path):
                        print('Did not found given path')
                        sys.exit(1)

                else:
                    if len(argv) == 1:
                        print(config_util_disclaimer)
                        sys.exit(1)
                    config_variable_name = argv[0]
                    config_variable_value = argv[1]

                argv = argv[2:]

    output_util = False

    if len(argv) > 0:
        if argv[0] == 'out':
            argv = argv[1:]
            possible_key_words = set(['suggest', 'PPI', 'LSS', 'FAT', 'MM', 'RAS', 'GIS', 'GIN', 'getsubgraph', 'isoform_relations', 'gene_report'])
            out_util_mode = None
            output_util = True
            if argv[0] in possible_key_words:
                out_util_mode = argv[0]
                argv = argv[1:]
            else:
                print(disclaimer)
                sys.exit(1)

    if len(argv) > 0:
        if argv[0] == 'set_m_key':
            if len(argv) != 2:
                print('Please provide the modeller key with "structman set_m_key [modeller key]"')
                sys.exit(1)
            modeller_key = argv[1]
            set_modeller_key(modeller_key)
            sys.exit(0)

    basic_util_mode = database_util or configure_mode
    util_mode = database_util or configure_mode or update_util or output_util

    # Custom single line input preparsing
    insert_flag_pos = None
    remove_args = []
    single_line_inputs = []
    for pos, arg in enumerate(argv):
        if insert_flag_pos is not None:
            # If the arg after the -i flag gives us a path, then it is not a single line input
            if pos == insert_flag_pos + 1:
                if os.path.exists(arg):
                    break
                else:
                    single_line_inputs.append(arg)
                    remove_args.append(insert_flag_pos)
                    remove_args.append(pos)
            else:
                if arg[0] == '-':
                    break
                else:
                    single_line_inputs.append(arg)
                    remove_args.append(pos)

        if arg == '-i':
            insert_flag_pos = pos
    # ignore --overwrrite by removing it, there will be a bug if someone combines single line input with overwrite
    for pos, arg in enumerate(argv):
        if arg == '--overwrite':
            remove_args.append(pos)
        if arg == 'pdb':
            remove_args.append(pos)
        if arg == 'rindb':
            remove_args.append(pos)
        if arg == 'mmdb':
            remove_args.append(pos)
        if arg == 'mmdb_from_scratch':
            remove_args.append(pos)            

    remove_args.reverse()
    for pos in remove_args:
        del argv[pos]

    try:
        long_paras = [
            'help', 'verbosity=', 'printerrors', 'printwarnings', 'chunksize=',
            'dbname=', 'fdbname=', 'restartlog', 'only_snvs', 'skip_indel_analysis', 'only_wt', 'mem_limit=',
            'model_indel_structures', 'ignore_local_pdb', 'ignore_local_rindb', 'ignore_local_mapping_db',
            'skip_main_output_generation', 'ray_local_mode', 'compute_ppi', 'structure_limiter=',
            'force_modelling', 'target=', 'custom_db', 'update_source=', 'condition_1=', 'condition_2=',
            'local_db=', 'read_only', 'read_hybrid', 'fast_pdb_annotation', 'target_db='
        ]
        opts, args = getopt.getopt(argv, "c:i:n:o:hp:", long_paras)

    except getopt.GetoptError:
        print("Illegal Input\n\n", disclaimer)
        sys.exit(2)

    infile = ''
    config_path = ''
    num_of_cores = None
    outfolder = ''
    verbosity = None

    print_all_errors = False
    print_all_warns = False
    chunksize = None
    dbname = None
    force_db_name = False
    restartlog = False
    only_snvs = False
    skip_indel_analysis = False
    only_wt = False
    model_indel_structures = None
    ignore_local_pdb = False
    ignore_local_rindb = False
    ignore_local_mapping_db = False
    skip_main_output_generation = False
    ray_local_mode = False
    compute_ppi = True
    overwrite_force_modelling = None
    condition_1_tag = None
    condition_2_tag = None
    target = None
    update_source = None
    local_db = None
    read_only = None
    read_hybrid = None
    structure_limiter = None
    fast_pdb_annotation = None
    target_db = None

    '''
    #mmcif mode flag is added
    mmcif_mode = False
    '''
    minus_p_path = None  # different modes can use this way a path given with -p
    mem_limit = None

    for opt, arg in opts:
        if opt == '-c':
            config_path = arg
        if opt == '-i':
            infile = arg
        if opt == '-n':
            num_of_cores = int(arg)
        if opt == '-o':
            outfolder = arg
        if opt == '-p':
            minus_p_path = arg
            if not os.path.exists(minus_p_path):
                try:
                    os.makedirs(minus_p_path)
                except:
                    print('Did not found given path', minus_p_path)
                    sys.exit(1)

        if opt == '-h' or opt == '--help':
            print(disclaimer)
            sys.exit(0)

        if opt == '--verbosity':
            verbosity = int(arg)

        if opt == '--printerrors':
            print_all_errors = True

        if opt == '--printwarnings':
            print_all_warns = True

        if opt == '--chunksize':
            chunksize = int(arg)

        if opt == '--dbname':
            dbname = arg
        if opt == '--fdbname':
            force_db_name = True
            dbname = arg
        if opt == '--restartlog':
            restartlog = True
        '''
        #mmcif option added to call mmcif mode while calling structman
        if opt == '--mmcif':
            mmcif_mode = True
        '''
        if opt == '--only_snvs':
            only_snvs = True

        if opt == '--skip_indel_analysis':
            skip_indel_analysis = True

        if opt == '--only_wt':
            only_wt = True
        
        if opt == '--mem_limit':
            mem_limit = float(arg)

        if opt == '--model_indel_structures':
            model_indel_structures = True

        if opt == '--ignore_local_pdb':
            ignore_local_pdb = True

        if opt == '--ignore_local_rindb':
            ignore_local_rindb = True

        if opt == '--ignore_local_mapping_db':
            ignore_local_mapping_db = True

        if opt == '--skip_main_output_generation':
            skip_main_output_generation = True

        if opt == '--ray_local_mode':
            ray_local_mode = True

        if opt == '--compute_ppi':
            compute_ppi = True

        if opt == '--structure_limiter':
            structure_limiter = arg

        if opt == '--force_modelling':
            overwrite_force_modelling = True

        if opt == '--condition_1':
            condition_1_tag = arg

        if opt == '--condition_2':
            condition_2_tag = arg

        if opt == '--target':
            target = arg

        if opt == '--update_source':
            update_source = arg

        if opt == '--local_db':
            local_db = arg

        if opt == '--read_only':
            read_only = True

        if opt == '--read_hybrid':
            read_hybrid = True

        if opt == '--fast_pdb_annotation':
            fast_pdb_annotation = True

        if opt == '--target_db':
            target_db = arg

    if not output_util and not util_mode:
        if infile == '' and len(single_line_inputs) == 0:
            input_folder = '/structman/input_data/'
            if not os.path.isdir(input_folder):
                print('Did not find the input path\n\n', disclaimer)
                sys.exit(2)
            filelist = os.listdir(input_folder)

            infiles = []
            for infile in filelist:
                if infile.split('.')[-1] == 'smlf' or infile.split('.')[-1] == 'vcf':
                    infiles.append(os.path.join(input_folder, infile))
            if infiles == []:
                print('Did not find any file inside the input folder\n\n', disclaimer)
                sys.exit(2)

        elif infile == '' and len(single_line_inputs) != 0:
            infiles = [single_line_inputs]
        else:
            # make infile to global path
            infile = resolve_path(infile)
            infiles = [infile]

    if not util_mode:
        if outfolder == '':
            outfolder = '/structman/results/'
            if not os.path.isdir(outfolder):
                if infile != '':
                    outfolder = resolve_path(os.path.join(os.path.dirname(infile), 'Output'))
                    if not os.path.exists(outfolder):
                        os.makedirs(outfolder)
                elif len(single_line_inputs) != 0:
                    outfolder = os.getcwd()
                else:
                    outfolder = None

    if config_path == '':
        # search highest to lowest priority
        search_paths = [
            os.path.dirname(infile),  # same directory as input file
            os.getcwd(),  # directory where structman command was triggered
            os.path.dirname(os.path.dirname(main_file_path)),  # top level of structman repo
            os.path.dirname(main_file_path),  # inside structman package (same level as structman_main.py)
            os.path.expanduser('~'),  # user's home folder
            '/structman/resources/', # container source folder
        ]
        search_names = ['config.txt', 'structman_config.txt']  # possible config filenames
        # search all possible paths with all possible filenames
        search_files = [os.path.join(path, name) for path in search_paths for name in search_names]
        for cfg_file in search_files:
            if os.path.exists(cfg_file):
                config_path = resolve_path(cfg_file)
                break

        if not config_path:
            print("No config file found, please use -c [Path to config]")
            sys.exit(2)


    config = Config(config_path, num_of_cores=num_of_cores, dbname = dbname,
                    output_path=outfolder, util_mode=util_mode, configure_mode=configure_mode, local_db = local_db, db_mode = db_mode,
                    basic_util_mode=basic_util_mode, output_util=output_util, external_call=False, verbosity=verbosity,
                    print_all_errors=print_all_errors, print_all_warns=print_all_warns, restartlog=restartlog, compute_ppi = compute_ppi)

    config.only_snvs = only_snvs
    config.skip_indel_analysis = skip_indel_analysis
    config.only_wt = only_wt
    config.ray_local_mode = ray_local_mode
    config.structure_limiter = structure_limiter
    config.condition_1_tag = condition_1_tag
    config.condition_2_tag = condition_2_tag
    
    if read_only is not None:
        config.read_only_mode = read_only

    if read_hybrid is not None:
        config.read_hybrid_mode = read_hybrid

    if fast_pdb_annotation is not None:
        config.fast_pdb_annotation = fast_pdb_annotation

    if overwrite_force_modelling is not None:
        config.force_modelling = True

    if model_indel_structures is not None:
        config.model_indel_structures = model_indel_structures

    if ignore_local_pdb:
        config.pdb_path = ''

    if ignore_local_rindb:
        config.rin_db_path = ''

    config.ignore_local_mapping_db = ignore_local_mapping_db
    if ignore_local_mapping_db:
        config.mapping_db_is_set = False

    config.skip_main_output_generation = skip_main_output_generation

    if mem_limit is not None:
        lim = ray._private.utils.get_system_memory() * mem_limit
        limit_memory(lim)
        print(f'Limit the memory consumption to {str(100 * mem_limit)}% ({str(lim / 1024. / 1024. / 1024.)} GB of RAM)')

    if chunksize is not None:
        config.chunksize = chunksize

    if config.verbosity >= 1:
        print(f"Using following config file: {config_path}")

    if dbname is not None:
        if not force_db_name:
            if dbname[:len(config.db_user_name)] != config.db_user_name:
                dbname = f'{config.db_user_name}_{dbname}'
        config.db_name = dbname

    if len(single_line_inputs) == 0:
        infile = resolve_path(infile)

    returncode = 0
    if database_util:
        if db_mode == 'reset':
            repairDB.empty(config)
        elif db_mode == 'clear':
            repairDB.clear(config)
        elif db_mode == 'out':
            db = MySQLdb.connect(config.db_address, config.db_user_name, config.db_password, config.db_name)
            cursor = db.cursor()
            # print(infile)
            session_id = database.getSessionId(infile, db, cursor)
            trunk, infilename = os.path.split(infile)
            session_name = infilename.rsplit('.', 1)[0]
            outfolder = resolve_path(os.path.join(trunk, 'Output'))
            outfile = os.path.join(outfolder, f'{session_name}_ligandlist.tsv')

            database.getLigandList(db, cursor, session_id, outfile)
            db.close()
        elif db_mode == 'create':
            if infile is not None and os.path.isfile(infile):
                config.database_source_path = infile

            repairDB.load(config)
        elif db_mode == 'destroy':
            repairDB.destroy(config)
        elif db_mode == 'export':
            if minus_p_path is None:
                print('database export needs a target_path given with -p [path to the target folder]')
                sys.exit(1)

            repairDB.export(config, minus_p_path)
        elif db_mode == 'reduce':
            if infile is None:
                sys.exit(1)
            repairDB.reduceToStructures(config, infile)
        elif db_mode == 'remove_sessions':
            repairDB.remove_sessions(config)

        elif db_mode == 'fix':
            repairDB.fix_database(config)

        elif db_mode == 'clone':
            if target_db is None:
                print('Database clone needs "--target_db" to be specified, Abort ...')
                sys.exit(1)
            repairDB.clone_db(config, target_db)

    elif output_util:
        if config.verbosity >= 1:
            print(infile)
        session_id = database.getSessionId(infile, config)
        trunk, infilename = os.path.split(infile)
        session_name = infilename.rsplit('.', 1)[0]
        outfolder = os.path.join(trunk, 'Output')

        if out_util_mode == 'suggest':
            outfile = os.path.join(outfolder, f'{session_name}_template_suggestions.tsv')
            anno_table = os.path.join(outfolder, f'{session_name}_full_annotation_table.tsv')
            if not os.path.exists(anno_table):
                annotation.create_annotation_table(session_id, config, anno_table)
            output.suggest(session_id, config, anno_table, outfile, ranking='consensus')
        elif out_util_mode == 'PPI':
            outfile = os.path.join(outfolder, f'{session_name}_ppi_network.tsv')
            ppi_output.create_ppi_network(session_id, config, outfile)
        elif out_util_mode == 'LSS':
            searchLargeStructures.search(config, chunksize, infile=infile)
        elif out_util_mode == 'FAT':
            outfile = os.path.join(outfolder, f'{session_name}_full_annotation_table.tsv')
            annotation.create_annotation_table(session_id, config, outfile)
        elif out_util_mode == 'MM':
            outfolder = os.path.join(outfolder, session_name)
            massmodel.mass_model(session_id, config, outfolder)
        elif out_util_mode == 'RAS':
            outfile = os.path.join(outfolder, f'{session_name}_all_annotated_structures.smlf')
            output.retrieve_annotated_structures(session_id, config, outfile)
        elif out_util_mode == 'GIS':
            out_f = f'{outfolder}/{session_name}'
            genes_out_package.generate_gene_isoform_scores(session_id, config, out_f, session_name)
        elif out_util_mode == 'GIN':
            out_f = f'{outfolder}/{session_name}'
            genes_out_package.generate_gene_interaction_network(session_id, config, out_f, session_name)
        elif out_util_mode == 'getsubgraph':
            if target is None:
                print(f'getsubgraph needs a target gene ID given with --target [gene ID]')
                sys.exit(1)
            genes_out_package.retrieve_sub_graph(infile, config, trunk, session_name, target)
        elif out_util_mode == 'isoform_relations':
            out_f = f'{outfolder}/{session_name}'
            genes_out_package.generate_isoform_relations(session_id, config, out_f, session_name)
        elif out_util_mode == 'gene_report':
            if target is None:
                print(f'gene_report needs a target gene ID given with --target [gene ID]')
                sys.exit(1)
            if session_id is None:
                print(f'Couldnt find input file {infile} in database. Is the correct database (--dbname) provided?')
                sys.exit(1)
            out_f = f'{outfolder}/{session_name}'
            genes_out_package.create_gene_report(session_id, target, config, out_f, session_name)

    elif update_util:
        if minus_p_path is not None:
            config.pdb_path = minus_p_path
            config.config_parser_obj.set('user', 'pdb_path', minus_p_path)
            f = open(config_path, 'w')
            config.config_parser_obj.write(f)
            f.close()
        elif config.container_version:
            if update_pdb:
                if not os.path.exists('/structman/resources/pdb'):
                    os.mkdir('/structman/resources/pdb')
                config.pdb_path = '/structman/resources/pdb'
                config.config_parser_obj.set('user', 'pdb_path', '/structman/resources/pdb')

            if update_alphafold_db:
                if not os.path.exists('/structman/resources/alphafold_db'):
                    os.mkdir('/structman/resources/alphafold_db')
                config.path_to_model_db = '/structman/resources/alphafold_db'
                config.config_parser_obj.set('user', 'path_to_model_db', '/structman/resources/alphafold_db')

            if update_rindb:
                if not os.path.exists('/structman/resources/rindb'):
                    os.mkdir('/structman/resources/rindb')
                config.rin_db_path = '/structman/resources/rindb'
                config.config_parser_obj.set('user', 'rin_db_path', '/structman/resources/rindb')

            if update_mapping_db:
                config.mapping_db = 'struct_man_db_mapping'
                config.config_parser_obj.set('user', 'mapping_db', 'struct_man_db_mapping')
                config.check_mapping_db()

            f = open(config_path, 'w')
            config.config_parser_obj.write(f)
            f.close()
        update.main(config, skipUpdatePDB=not update_pdb,
                            skip_rindb=not update_rindb,
                            rin_fromScratch=update_rindb_from_scratch,
                            update_mapping_db = update_mapping_db,
                            mapping_db_from_scratch = update_mapping_db_from_scratch,
                            update_mapping_db_keep_raw_files = update_mapping_db_keep_raw_files,
                            update_alphafold_db = update_alphafold_db,
                            update_source = update_source,
                            update_microminer = update_microminer,
                            microminer_from_scratch = microminer_from_scratch,
                            check_search_db = check_search_db
                    )

    elif configure_mode:
        if configure_help_mode:
            print(config_util_disclaimer)
            print(f'Current path to config: {config_path}')
            if os.path.isfile(config_path):
                with open(config_path, 'r') as f:
                    page = f.read()
                    print('\n===== Content of the config: ======\n')
                    print(page)
            print('\n===== Testing database configuration: ======\n')
            print(f'{config.db_address=} {config.db_name=}')
            db, cursor = config.getDB()
            if db is not None:
                db.close()
                print('DB connection was successful!')
            sys.exit(0)
        if conf_update_pdb_path is not None:
            config.config_parser_obj.set('user', 'pdb_path', conf_update_pdb_path)

        if config_variable_name is not None:
            config.config_parser_obj.set('user', config_variable_name, config_variable_value)

        f = open(config_path, 'w')
        config.config_parser_obj.write(f)
        f.close()

        f = open(config_path, 'r')
        lines = f.readlines()
        f.close()

        newlines = []
        for line in lines:
            if line == '' or line == '\n':
                continue
            newlines.append(line)

        f = open(config_path, 'w')
        f.write(''.join(newlines))
        f.close()

    else:        
        returncode = main(infiles, outfolder, config)

    sys.exit(returncode)

if __name__ == "__main__":
    structman_cli()
