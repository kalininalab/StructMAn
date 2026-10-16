import gzip

from structman.base_utils.base_utils import alphafold_model_id_to_db_key, decode_zstd, create_mmseqs_index
from structman.lib.database.database_core_functions import select
from structman.base_utils.config_class import Config
from structman.scripts.createPdbBaDb import parseSequences, great_seq_map_to_fasta

def get_afdb_data(model_id: str, config: Config, columns_tuple: list[str]):
    if model_id.startswith('afdb_hetero_'):
        model_id = model_id[12:]
        db_name = config.afdb_hetero
    elif model_id.startswith('afdb_homo_'):
        model_id = model_id[10:]
        db_name = config.afdb_homo
    else:
        db_name = config.afdb

    afdb_id, main_db, chunk_id = alphafold_model_id_to_db_key(model_id)
    if main_db:
        subdb = 'main'
    else:
        subdb = 'collab'
    table_name = f'AFDB_{subdb}_{chunk_id}'

    if config.verbosity >= 6:
        config.logger.info(f'get_afdb_data: {model_id=} {columns_tuple=} {afdb_id=} {table_name=}')

    try:
        results = select(config, columns_tuple, table_name, equals_rows={'ID':afdb_id}, db_name=db_name, locked=True)
    except (NameError, BaseException) as err:
        config.errorlog.add_warning(f'Could not retrieve {afdb_id=} from AFDB, {err=} {model_id=} {columns_tuple=} {afdb_id=} {table_name=}')
        return None
    
    if len(results) == 0:
        config.errorlog.add_warning(f'Could not retrieve {afdb_id=} from AFDB, {len(results)==0=} {model_id=} {columns_tuple=} {afdb_id=} {table_name=}')
        return None

    pages = []
    for i in range(len(columns_tuple)//2):
        try:
            data, compr_type = results[0][(i*2):(i*2)+2]
        except ValueError as err: 
            raise ValueError(f'{err=} {columns_tuple=} {results=}')

        if compr_type == 0:
            page = gzip.decompress(data)
        elif compr_type == 1:
            page = decode_zstd(data)
        else:
            page = data
        pages.append(page)
    if len(pages) == 1:
        return pages[0]
    return pages

def get_afdb_structure(model_id: str, config: Config):
    columns = ['Structure_data', 'Sd_compr_type']
    return get_afdb_data(model_id, config, columns)

def get_afdb_msa(model_id: str, config: Config):
    columns = ['MSA_data', 'Md_compr_type']
    return get_afdb_data(model_id, config, columns)

def get_afdb_rin_and_cent(model_id: str, config: Config):
    columns = ['RIN_data', 'Rd_compr_type', 'Centralities_data', 'Cd_compr_type', 'In_out_flow_data', 'Iod_compr_type']
    return get_afdb_data(model_id, config, columns)

def get_single_sequence(config: Config, prot_db_id: int):
    columns = ['Sequence']
    results = select(config, columns, 'Protein', equals_rows={'Protein_Id':prot_db_id})
    return results[0][0]

def afdb_create_mmseqs_index(config: Config, db_name: str, mmseqs_fasta_path: str):
    db, cursor = config.getDB(server_connection=True, db_name=db_name)
    cursor.execute(f"SHOW TABLES FROM {db_name}")
    
    tables = []

    for row in cursor.fetchall():
        tables.append(row[0])
        
    db.close()

    great_seq_map = {}

    for table_name in tables:
        columns_tuple = ['ID', 'Structure_data', 'Sd_compr_type']
        results = select(config, columns_tuple, table_name, db_name=db_name)

        for row in results:
            afdb_id, data, compr_type = row
            if compr_type == 0:
                page = gzip.decompress(data)
            elif compr_type == 1:
                page = decode_zstd(data)
            else:
                continue

            seqs = parseSequences(page.decode())

            for chain in seqs:
                seq = seqs[chain]
                if len(seq) > 10:
                    if seq not in great_seq_map:
                        great_seq_map[seq] = [(afdb_id, chain)]
                    else:
                        great_seq_map[seq].append((afdb_id, chain))

    great_seq_map_to_fasta(great_seq_map, mmseqs_fasta_path)

    search_db_base_path = mmseqs_fasta_path.rsplit('/',1)[0]

    create_mmseqs_index(mmseqs_fasta_path, search_db_base_path, config.mmseqs_tmp_folder)