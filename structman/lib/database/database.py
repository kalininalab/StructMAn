import cProfile
import multiprocessing
import re
import sys
import time
import traceback

import ray

from structman.base_utils.base_utils import calc_checksum, pack, unpack, calculate_chunksizes
from structman.base_utils import ray_utils
from structman.lib import pdbParser
from structman.lib import rin
from structman.lib import serializedPipeline
from structman.lib import uniprot
from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc import complex as complex_package
from structman.lib.sdsc import structure as structure_package
from structman.lib.sdsc import residue as residue_package
from structman.lib.sdsc import protein as protein_package
from structman.lib.sdsc import position as position_package
from structman.lib.sdsc import mappings as mappings_package
from structman.lib.sdsc import snv as snv_package
from structman.lib.sdsc import indel as indel_package
from structman.lib.sdsc import interface as interface_package
from structman.lib.sdsc.consts import residues as residue_consts
from structman.lib.database.retrieval import getStoredResidues
from structman.lib.database.database_core_functions import binningSelect, select, insert, update


# called by output.classification
def getProteinDict(prot_id_list, session_id, config, includeSequence=False):
    if len(prot_id_list) == 0:
        return {}, set()

    table = 'RS_Protein_Session'
    columns = ['Protein', 'Session', 'Input_Id']

    if session_id is not None:
        results = select(config, columns, table, equals_rows={'Session': session_id})

        prot_input_dict = {}

        for row in results:
            if row[1] != session_id:
                continue
            prot_id = row[0]
            if prot_id in prot_id_list:
                prot_input_dict[prot_id] = row[2]
    else:
        prot_input_dict = set(prot_id_list)

    if not includeSequence:
        rows = ['Protein_Id', 'Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id', 'Error_Code', 'Error', 'Gene']
    else:
        rows = ['Protein_Id', 'Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id', 'Error_Code', 'Error', 'Gene', 'Sequence']
    table = 'Protein'
    results = binningSelect(prot_id_list, rows, table, config)

    gene_db_ids = set()

    protein_dict = {}
    for row in results:
        prot_id = row[0]
        gene_db_id = row[7]
        if gene_db_id is not None:
            gene_db_ids.add(gene_db_id)
        if prot_id in prot_input_dict:
            if session_id is not None:
                input_id = prot_input_dict[prot_id]
            else:
                input_id = None
            if not includeSequence:
                protein_dict[prot_id] = (row[1], row[2], row[3], row[4], row[5], row[6], input_id, gene_db_id, None)
            else:
                protein_dict[prot_id] = (row[1], row[2], row[3], row[4], row[5], row[6], input_id, gene_db_id, row[8])

    return protein_dict, gene_db_ids


# called from serializedPipeline
def geneCheck(genes, config):
    results = select(config, ['Gene_Id', 'Primary_Gene_Id'], 'Gene')
    for row in results:
        gene_db_id = row[0]
        gene_id = row[1]
        if not gene_id in genes:
            continue
        genes[gene_id].database_id = gene_db_id

    if not config.read_only_mode:
        values = []
        for gene_id in genes:
            if genes[gene_id].database_id is not None:
                continue
            gene_name = genes[gene_id].gene_name
            values.append([gene_id, gene_name])

        if len(values) > 0:
            insert('Gene', ['Primary_Gene_Id' ,'Gene_Name'], values, config)

            results = select(config, ['Gene_Id', 'Primary_Gene_Id'], 'Gene')
            for row in results:
                gene_db_id = row[0]
                gene_id = row[1]
                if not gene_id in genes:
                    continue
                genes[gene_id].database_id = gene_db_id
    return

# called from serializedPipeline
def protCheck(proteins, genes, session_id, config):
    # Scan the database for stored Proteins
    if proteins.isEmpty():
        return {}, {}

    results = select(config, ['Protein_Id', 'Primary_Protein_Id', 'Sequence'], 'Protein')

    if config.verbosity >= 3:
        print('Just after protCheck selection')

    prot_id_list = set([])
    prot_ids_mutants_excluded = set()

    max_database_id = 0

    for row in results:
        database_id = row[0]
        prot_id = row[1]
        if not proteins.contains(prot_id):
            continue
        proteins.set_protein_stored(prot_id, True)
        proteins.set_protein_db_id(prot_id, database_id)
        proteins.set_protein_sequence(prot_id, row[2])

        prot_id_list.add(database_id)
        if not sdsc_utils.is_mutant_ac:
            prot_ids_mutants_excluded.add(database_id)

        if database_id > max_database_id:
            max_database_id = database_id

    proteins.set_stored_ids(prot_id_list, prot_ids_mutants_excluded)

    prot_ids = proteins.get_protein_ids()

    # Insert the new proteins into the database
    new_prots = set()
    
    for prot_id in prot_ids:
        if not proteins.is_protein_stored(prot_id):
            new_prots.add(prot_id)

    if config.read_only_mode:
        for prot_id in new_prots:
            proteins.remove_protein(prot_id)
    else:
        if len(new_prots) > 0:
            values = []
            for prot_id in new_prots:
                u_ac = proteins.get_u_ac(prot_id)
                ref_id = proteins.get_ref_id(prot_id)
                u_id = proteins.get_u_id(prot_id)
                mut_type = proteins[prot_id].mutant_type
                gene_id = proteins[prot_id].gene
                if gene_id is not None:
                    if gene_id in genes:
                        gene_db_id = genes[gene_id].database_id
                    else:
                        gene_db_id = None
                else:
                    gene_db_id = None
                if mut_type is None:
                    mutant_positions_str = None
                else:
                    if proteins[prot_id].sav_positions is None:
                        sav_string = ''
                    else:
                        sav_string = ','.join([str(x) for x in proteins[prot_id].sav_positions])

                    if proteins[prot_id].insertion_positions is None:
                        insertion_string = ''
                    else:
                        insertion_parts = []
                        for ins in proteins[prot_id].insertion_positions:
                            insertion_parts.append(','.join([str(x) for x in ins]))
                        insertion_string = '!'.join(insertion_parts)

                    if proteins[prot_id].deletion_flanks is None:
                        deletion_string = ''
                    else:
                        deletion_parts = []
                        for deletion_fl in proteins[prot_id].deletion_flanks:
                            lf, rf = deletion_fl
                            if lf is None:
                                lfs = '_'
                            else:
                                lfs = str(lf)
                            if rf is None:
                                rfs = '_'
                            else:
                                rfs = str(rf)
                            deletion_parts.append(f'{lfs},{rfs}')
                        deletion_string = '!'.join(deletion_parts)

                    mutant_positions_str = f'{sav_string};{insertion_string};{deletion_string}'
                values.append((gene_db_id, prot_id, u_ac, ref_id, u_id, session_id, mut_type, mutant_positions_str))

            insert('Protein', ['Gene', 'Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id', 'Original_Session', 'Mutant_Type', 'Mutant_Positions'], values, config)

            # Retrieve the Protein-Ids from the new added proteins

            db, cursor = config.getDB()

            sql = "SELECT Protein_Id,Primary_Protein_Id FROM Protein WHERE Protein_Id > %s" % str(max_database_id)
            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Error in protCheck: %s,\n%s" % (sql, f))
            db.close()
        else:
            results = ()

        new_ids = set()
        for row in results:
            database_id = row[0]
            prot_id = row[1]
            if not proteins.contains(prot_id):
                continue
            proteins.set_protein_db_id(prot_id, database_id)

            if database_id not in prot_id_list:
                new_ids.add(database_id)

        proteins.set_not_stored_ids(new_ids)

    proteins.generate_id_map()

    # Insert the Protein-Session-Connections into the database
    values = []
    for prot_id in prot_ids:
        database_id = proteins.get_protein_database_id(prot_id)
        input_id = proteins[prot_id].input_id
        tags_str = ','.join(proteins[prot_id].tags)
        values.append((database_id, session_id, input_id, tags_str))

    insert('RS_Protein_Session', ['Protein', 'Session', 'Input_Id', 'Tags'], values, config)


def getLastAutoInc(config):
    db, cursor = config.getDB()

    if config.db_address == '-':
        sql = "SELECT last_insert_rowid()"
    else:
        sql = "SELECT LAST_INSERT_ID()"

    cursor.execute(sql)
    results = cursor.fetchall()
    for row in results:
        auto_id = row[0]
    db.commit()

    return auto_id


# called by babel
def getUniprotAcFromId(Protein_Id, db, cursor):
    sql = "SELECT Uniprot_Ac FROM Protein WHERE Protein_Id = '%s'" % (str(Protein_Id))

    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("NameError in getUniprotFromId")
        db.rollback()

    if results == ():
        return 0
    row = results[0]
    return row[0]


# called by babel
def getAAChange(mutation_id, db, cursor):
    sql = "SELECT Amino_Acid_Change FROM Position WHERE Position_Id = '%s'" % str(mutation_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getAAChange: %s" % sql)
        db.rollback()

    if results == ():
        return ""
    for row in results:
        aac = row[0]
        return aac
    return ""


def insertMultiMutations(proteins, genes, session, config):

    if config.verbosity >= 3:
        print(f'insertMultiMutations: {len(proteins.multi_mutations)}')

    db_ids = []
    for wt_prot in proteins.multi_mutations:
        db_ids.append(proteins.get_protein_database_id(wt_prot))

    if len(db_ids) > 0:
        columns = ['Wildtype_Protein', 'Multi_Mutation_Id', 'SNVs', 'Indels']
        table = 'Multi_Mutation'

        results = binningSelect(db_ids, columns, table, config)

        for row in results:
            wt_database_id = row[0]
            wt_prot_id = proteins.getU_acByDbId(wt_database_id)
            if wt_prot_id not in proteins.multi_mutations:
                continue
            snv_db_ids = row[2].split(',')
            indel_db_ids = row[3].split(',')

            for multi_mutation in proteins.multi_mutations[wt_prot]:
                if snv_db_ids != multi_mutation.get_snv_db_ids():
                    continue

                if indel_db_ids != multi_mutation.get_indel_db_ids():
                    continue
                multi_mutation.database_id = row[1]
                multi_mutation.stored = True

    obj_map = {}
    if config.read_only_mode:
        wt_prot_db_ids = db_ids
    else:
        wt_prot_db_ids = []
        values = []
        
        for wt_prot in proteins.multi_mutations:
            wt_db_id = proteins.get_protein_database_id(wt_prot)
            wt_prot_db_ids.append(wt_db_id)
            for multi_mutation in proteins.multi_mutations[wt_prot]:
                if multi_mutation.stored:
                    continue
                if multi_mutation.mut_prot is not None:
                    if multi_mutation.mut_prot in proteins.protein_map:
                        mut_db_id = proteins[multi_mutation.mut_prot].database_id
                    else:
                        mut_db_id = None
                else:
                    mut_db_id = None
                snv_db_ids = ','.join([str(x) for x in multi_mutation.get_snv_db_ids()])
                indel_db_ids = ','.join([str(x) for x in multi_mutation.get_indel_db_ids()])
                values.append((wt_db_id, mut_db_id, snv_db_ids, indel_db_ids))

                mm_obj_id = (wt_db_id, snv_db_ids, indel_db_ids)

                obj_map[mm_obj_id] = multi_mutation

        if len(values) > 0:
            columns = ['Wildtype_Protein', 'Mutant_Protein', 'SNVs', 'Indels']
            table = 'Multi_Mutation'
            insert(table, columns, values, config)

    columns = ['Wildtype_Protein', 'Multi_Mutation_Id', 'SNVs', 'Indels']
    table = 'Multi_Mutation'

    results = binningSelect(wt_prot_db_ids, columns, table, config)
    for row in results:
        wt_db_id = row[0]
        wt_prot_id = proteins.getU_acByDbId(wt_db_id)

        mm_obj_id = (wt_db_id, row[2], row[3])
        if not mm_obj_id in obj_map:
            continue

        obj_map[mm_obj_id].database_id = row[1]

    columns = ['Multi_Mutation', 'Session', 'Tags']
    table = 'RS_Multi_Mutation_Session'

    values = []
    isoform_values = []
    for wt_prot in proteins.multi_mutations:
        wt_db_id = proteins.get_protein_database_id(wt_prot)
        for multi_mutation in proteins.multi_mutations[wt_prot]:
            if multi_mutation.tags is not None:
                mm_str = ','.join(multi_mutation.tags)
            else:
                mm_str = None
            values.append((multi_mutation.database_id, session, mm_str))

            if proteins[wt_prot].gene is not None:
                mut_db_id = proteins[multi_mutation.mut_prot].database_id
                gene_db_id = genes[proteins[wt_prot].gene].database_id
                isoform_values.append((wt_db_id, mut_db_id, gene_db_id, session, multi_mutation.database_id))

    if len(values) > 0:
        insert(table, columns, values, config)

    if len(isoform_values) > 0:
        insert('RS_Isoform', ['Protein_A', 'Protein_B', 'Gene', 'Session', 'Multi_Mutation'], isoform_values, config)

    return


# called by serializedPipeline
def indelCheck(proteins, session, config):
    if len(proteins.indels) == 0:
        return
    stored_wt_ids = proteins.get_stored_ids(exclude_indel_mutants=True)

    if len(stored_wt_ids) > 0:
        rows = ['Wildtype_Protein', 'Indel_Notation', 'Indel_Id']
        table = 'Indel'

        results = binningSelect(stored_wt_ids, rows, table, config)

        for row in results:
            wt_prot_id = row[0]
            indel_notation = row[1]
            database_id = row[2]
            u_ac = proteins.getU_acByDbId(wt_prot_id)
            proteins.indels[u_ac][indel_notation].set_database_id(database_id)
            proteins.indels[u_ac][indel_notation].set_stored(True)

    if config.read_only_mode:
        all_wt_ids = stored_wt_ids
    else:
        values = []
        all_wt_ids = set()
        for u_ac in proteins.indels:
            for indel_notation in proteins.indels[u_ac]:
                if proteins.indels[u_ac][indel_notation].stored:
                    continue
                wt_prot_id = proteins.get_protein_database_id(proteins.indels[u_ac][indel_notation].wt_prot)
                mut_prot_id = proteins.get_protein_database_id(proteins.indels[u_ac][indel_notation].mut_prot)
                if wt_prot_id is None:
                    config.errorlog.add_error('Wildtype protein id is not allowed to be None')
                    continue
                values.append((wt_prot_id, mut_prot_id, indel_notation))
                all_wt_ids.add(wt_prot_id)

        if len(values) > 0:
            columns = ['Wildtype_Protein', 'Mutant_Protein', 'Indel_Notation']
            insert('Indel', columns, values, config)

    session_values = []

    if len(all_wt_ids) > 0:
        rows = ['Wildtype_Protein', 'Indel_Notation', 'Indel_Id']
        table = 'Indel'
        results = binningSelect(all_wt_ids, rows, table, config)
        for row in results:
            indel_notation = row[1]
            wt_prot_id = row[0]
            u_ac = proteins.getU_acByDbId(wt_prot_id)
            if u_ac not in proteins.indels:
                continue
            if indel_notation not in proteins.indels[u_ac]:
                continue
            database_id = row[2]
            proteins.indels[u_ac][indel_notation].set_database_id(database_id)
            session_values.append((session, database_id, ','.join(proteins.indels[u_ac][indel_notation].tags)))

    if len(session_values) > 0:
        columns = ['Session', 'Indel', 'Tags']
        insert('RS_Indel_Session', columns, session_values, config)


# called by serializedPipeline
def positionCheck(proteins, database_session, config):
    # search for stored positions
    stored_positions = []
    pos_map = {}
    if len(proteins.get_stored_ids()) > 0:
        columns = ['Protein', 'Position_Number', 'Position_Id', 'Recommended_Structure_Data']
        table = 'Position'

        results = binningSelect(proteins.get_stored_ids(), columns, table, config)

        for row in results:
            prot_id = row[0]

            pos = row[1]
            position_database_id = row[2]
            if not proteins.position_in_protein_by_db_id(prot_id, pos):
                continue

            proteins.set_position_stored(prot_id, pos, True)
            proteins.set_position_database_id(prot_id, pos, position_database_id)
            stored_positions.append(position_database_id)
            pos_map[position_database_id] = (prot_id, pos)
            if row[3] is not None:
                recommended_structure_tuple = unpack(row[3])
                recommended_structure, seq_id, cov, resolution = sdsc_utils.process_recommend_structure_str(recommended_structure_tuple[0])
                proteins.getByDbId(prot_id).positions[pos].recommended_structure = recommended_structure

    # search for stored SNVs
    if len(stored_positions) > 0:
        columns = ['Position', 'New_AA', 'SNV_Id']
        table = 'SNV'
        results = binningSelect(stored_positions, columns, table, config)

        for row in results:
            prot_database_id, pos = pos_map[row[0]]
            prot_id = proteins.getU_acByDbId(prot_database_id)
            if not row[1] in proteins[prot_id].positions[pos].mut_aas:
                continue
            proteins[prot_id].positions[pos].mut_aas[row[1]].database_id = row[2]
            proteins[prot_id].positions[pos].mut_aas[row[1]].stored = True

    prot_ids = proteins.get_protein_ids()

    if not config.read_only_mode:
        # insert new positions
        
        values = []
        for prot_id in prot_ids:
            prot_database_id = proteins.get_protein_database_id(prot_id)
            positions = proteins.get_position_ids(prot_id)
            all_stored = True
            for pos in positions:
                if proteins.is_position_stored(prot_id, pos):
                    continue
                all_stored = False
                aac_base = proteins.get_aac_base(prot_id, pos)

                res_id = proteins.get_res_id(prot_id, pos)

                values.append((prot_database_id, int(aac_base[1:]), res_id, aac_base[0]))

            if all_stored:
                proteins.set_completely_stored(prot_id)

        if len(values) > 0:
            columns = ['Protein', 'Position_Number', 'Residue_Id', 'Wildtype_Residue']
            insert('Position', columns, values, config)

        # retrieve the database ids of the new positions
        columns = ['Protein', 'Position_Number', 'Position_Id']
        table = 'Position'

        fused_prot_ids = proteins.get_not_stored_ids() | proteins.get_stored_ids()

        results = binningSelect(fused_prot_ids, columns, table, config)

        position_database_ids = []
        prot_back_map = {}
        for row in results:
            prot_database_id = row[0]

            pos = row[1]

            position_database_id = row[2]
            if not proteins.position_in_protein_by_db_id(prot_database_id, pos):
                continue
            proteins.set_position_database_id(prot_database_id, pos, position_database_id)
            position_database_ids.append(position_database_id)
            prot_back_map[position_database_id] = (prot_database_id, pos)

        # insert the new SNVs
        if config.verbosity >= 6:
            for prot_id in prot_ids:
                proteins[prot_id].print_all_snvs()

        values = []
        for prot_id in prot_ids:
            for pos_obj in proteins[prot_id].positions:
                if pos_obj is None:
                    continue
                pos_database_id = pos_obj.database_id
                #proteins[prot_id].positions[pos].print_state()
                for new_aa in pos_obj.mut_aas:
                    values.append((pos_database_id, new_aa))

        if config.verbosity >= 6:
            print(f'SNV values inserted into Database: {values}')
        
        if len(values) > 0:
            table = 'SNV'
            columns = ['Position', 'New_AA']
            insert(table, columns, values, config)

        # retrieve the database ids of the new SNVs
        columns = ['Position', 'New_AA', 'SNV_Id']
        table = 'SNV'
        results = binningSelect(position_database_ids, columns, table, config)

        for row in results:
            position_database_id = row[0]
            if position_database_id not in prot_back_map:
                continue
            (prot_database_id, pos) = prot_back_map[position_database_id]
            prot_id = proteins.getU_acByDbId(prot_database_id)
            new_aa = row[1]
            if new_aa not in proteins[prot_id].positions[pos].mut_aas:
                continue
            proteins[prot_id].positions[pos].mut_aas[new_aa].database_id = row[2]

    # insert the the SNV session connections
    values = []
    for prot_id in prot_ids:
        for pos_obj in proteins[prot_id].positions:
            if pos_obj is None:
                continue
            for new_aa in pos_obj.mut_aas:
                tags = ','.join(pos_obj.mut_aas[new_aa].tags)
                snv_database_id = pos_obj.mut_aas[new_aa].database_id
                values.append((database_session, snv_database_id, tags))

    if len(values) > 0:
        table = 'RS_SNV_Session'
        columns = ['Session', 'SNV', 'Tag']
        insert(table, columns, values, config)

    # insert the position session connections
    values = []
    for prot_id in prot_ids:
        positions = proteins.get_position_ids(prot_id)
        for pos in positions:
            if proteins[prot_id].positions[pos].session_less:
                continue
            pos_id = proteins.get_position_database_id(prot_id, pos)
            pos_tags = proteins.get_pos_tags(prot_id, pos)
            if pos_tags is None:
                pos_tags_str = None
            else:
                pos_tags_str = ','.join(pos_tags)
            values.append((database_session, pos_id, pos_tags_str))

    """
    if len(values) > 10000:
        process = multiprocessing.Process(target=backgroundInsertMS, args=(values, config))
        try:
            process.start()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in PositionCheck: %s" % (f))
        return process
    else:
    """
    insert('RS_Position_Session', ['Session', 'Position', 'Tag'], values, config)
    return None


def backgroundInsertMS(values, config):
    insert('RS_Position_Session', ['Session', 'Position', 'Tag'], values, config)


# called by serializedPipeline
def addIupred(proteins, config):
    values = []
    u_acs = proteins.get_protein_ids()
    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        scores = proteins.get_disorder_scores(u_ac)
        regions = proteins.get_disorder_regions(u_ac)
        method = proteins.get_disorder_tool(u_ac)

        positions = proteins.get_position_ids(u_ac)

        for pos in positions:
            pos_id = proteins.get_position_database_id(u_ac, pos)

            if method == 'MobiDB3.0' or method == 'mobidb-lite':
                pos_region_type = 'globular'
            else:
                pos_region_type = 'disorder'
            if regions is None:
                continue
            for [a, b, region_type] in regions:
                if int(pos) > int(a) and int(pos) < int(b):
                    pos_region_type = region_type

            if scores is None:
                config.errorlog.add_warning('IUpred scores are None for: %s' % u_ac)
                break

            if pos not in scores:
                continue
            iupred_score = scores[pos]
            values.append((pos_id, iupred_score, pos_region_type))

    if config.verbosity >= 2:
        print(f'Before background add disorder process: {len(values)}')

    process = multiprocessing.Process(target=backgroundIU, args=(pack(values), config))
    process.start()
    #insert_IU(values, config)
    #process = None
    #process = None
    #if not len(values) == 0:
    #    update(config, 'Position', ['Position_Id', 'IUPRED', 'IUPRED_Glob'], values)
    return process


def backgroundIU(values, config):
    values = unpack(values)
    insert_IU(values, config)

def insert_IU(values, config): 
    if not len(values) == 0:
        update(config, 'Position', ['Position_Id', 'IUPRED', 'IUPRED_Glob'], values)
    return

# called by babel
def getPDB(template_id, db, cursor):
    sql = "SELECT Name FROM Template WHERE Template_Id = '%s'" % str(template_id)

    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        print("NameError in getPDB")
        db.rollback()

    return results[0][0]


# called by babel
def checkMutationSession(mutation_id, session_id, db, cursor):
    sql = "SELECT Position FROM RS_Position_Session WHERE Position = '%s' AND Session = '%s'" % (str(mutation_id), str(session_id))
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in checkMutationSession")
        db.rollback()
    if results == ():
        return False
    else:
        return True


# called by output.classification
def getMutationDict(mutation_id_list, config):
    if len(mutation_id_list) == 0:
        return {}

    rows = ['Position_Id', 'Position_Number', 'Protein', 'IUPRED', 'IUPRED_Glob', 'Residue_Id', 'Wildtype_Residue']
    table = 'Position'
    results = binningSelect(mutation_id_list, rows, table, config)

    mutation_dict = {}

    for row in results:
        mut_id = row[0]
        if mut_id in mutation_id_list:
            mutation_dict[mut_id] = (row[1], row[2], row[3], row[4], row[5], row[6])
    return mutation_dict


# called by serializedPipeline
def updateSession(session_id, time, config):
    update(config, 'Session', ['Session_Id', 'End'], [(str(session_id), str(time))])


# called by serializedPipeline
def insertSession(time, ori_filename, config):
    
    checksum = calc_checksum(ori_filename)

    values = [(ori_filename, checksum, str(time))]
   
    if len(values) > 0:
        insert('Session', ['Input_File', 'Checksum', 'Start'], values, config)

    session_id = getSessionId(ori_filename, config)

    return session_id

def getSessionId(infile, config):
    checksum = calc_checksum(infile)

    table = 'Session'
    columns = ['Session_Id', 'Checksum']

    if config.verbosity >= 3:
        print(f'Looking for infile {infile} with checksum {checksum} in getSessionId')

    results = select(config, columns, table, equals_rows={'Input_File': infile})

    if len(results) == 0:
        return None

    for row in results:
        if checksum == row[1]:
            return row[0]

    print(f'Error in getSessionId:\nChecksum of infile {infile} not found in database:\n{results}')

    return None


# called by babel
# called by output
def getSessionFile(session_id, config):
    table = 'Session'
    columns = ['Input_File']

    results = select(config, columns, table, equals_rows={'Session_Id': session_id})

    if results == ():
        return None
    try:
        return results[0][0]
    except:
        return None

# called by serializedPipeline
def addProtInfos(proteins, config):
    ref_value_strs = []
    go_term_ids = set()
    reac_ids = set()
    seq_values = []
    u_acs = proteins.get_protein_ids()

    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            go_term_ids.add(go_id)
        for reac_id in pathways:
            reac_ids.add(reac_id)
        prot_id = proteins.get_protein_database_id(u_ac)
        seq = proteins.get_sequence(u_ac)
        seq_values.append((prot_id, seq))

    stored_go_terms = {}
    if len(go_term_ids) > 0:
        results = select(config, ['Id', 'GO_Term_Id'], 'GO_Term', in_rows={'ID': go_term_ids})

        for row in results:
            stored_go_terms[row[0]] = row[1]

    stored_pathways = {}
    if len(reac_ids) > 0:
        results = select(config, ['Reactome_Id', 'Pathway_Id'], 'Pathway', in_rows={'Reactome_Id': reac_ids})

        for row in results:
            stored_pathways[row[0]] = row[1]

    new_go_terms = {}
    go_values = []
    new_pathways = {}
    pathway_values = []
    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            if go_id not in stored_go_terms:
                new_go_terms[go_id] = go_terms[go_id]
                go_values.append((go_terms[go_id].replace("'", "''"), go_id))
                stored_go_terms[go_id] = None
        for reac_id in pathways:
            if reac_id not in stored_pathways:
                new_pathways[reac_id] = pathways[reac_id]
                pathway_values.append((pathways[reac_id].replace("'", "''"), reac_id))
                stored_pathways[reac_id] = None

    if len(go_values) > 0:
        insert('GO_Term', ['Name', 'Id'], go_values, config)

    if len(pathway_values) > 0:
        insert('Pathway', ['Name', 'Reactome_Id'], pathway_values, config)

    stored_go_terms = {}
    if len(go_term_ids) > 0:
        results = select(config, ['Id', 'GO_Term_Id', ], 'GO_Term', in_rows={'Id': go_term_ids})

        for row in results:
            stored_go_terms[row[0]] = row[1]

    stored_pathways = {}
    if len(reac_ids) > 0:
        results = select(config, ['Reactome_Id', 'Pathway_Id'], 'Pathway', in_rows={'Reactome_Id': reac_ids})

        for row in results:
            stored_pathways[row[0]] = row[1]

    go_values = []
    pathway_values = []
    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        prot_id = proteins.get_protein_database_id(u_ac)
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            go_values.append((prot_id, stored_go_terms[go_id]))
        for reac_id in pathways:
            pathway_values.append((prot_id, stored_pathways[reac_id]))

    if len(go_values) > 0:
        insert('RS_Protein_GO_Term', ['Protein', 'GO_Term'], go_values, config)

    if len(pathway_values) > 0:
        insert('RS_Protein_Pathway', ['Protein', 'Pathway'], pathway_values, config)

    if len(seq_values) > 0:
        update(config, 'Protein', ['Protein_Id', 'Sequence'], seq_values)


def getComplexMap(config, pdb_ids=None):
    if pdb_ids is not None:
        if len(pdb_ids) == 0:
            return {}
    table = 'Complex'
    rows = ['Complex_Id', 'PDB', 'Resolution', 'Chains', 'Homooligomers', 'Ligand_Profile', 'Metal_Profile', 'Ion_Profile', 'Chain_Chain_Profile']

    if pdb_ids is None:
        results = select(config, rows, table)
    elif len(pdb_ids) > 50:
        results = select(config, rows, table)
    else:
        results = select(config, rows, table, in_rows = {'PDB': pdb_ids})

    complex_map = {}
    for row in results:
        pdb_id = row[1]
        if pdb_ids is not None:
            if pdb_id not in pdb_ids:
                continue
        comp_id = row[0]
        resolution = row[2]
        chains_str = row[3]
        homooligomers = row[4]
        lig_profile_str = row[5]
        metal_profile_str = row[6]
        ion_profile_str = row[7]
        cc_profile_str = row[8]

        lig_profile = {}
        if lig_profile_str != '':
            for lig_profile_part in lig_profile_str.split(','):
                part1, part2 = lig_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                lig_profile[(chain, res)] = int(deg), float(score)

        metal_profile = {}

        if metal_profile_str != '':
            for metal_profile_part in metal_profile_str.split(','):
                part1, part2 = metal_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                metal_profile[(chain, res)] = int(deg), float(score)

        ion_profile = {}

        if ion_profile_str != '':
            for ion_profile_part in ion_profile_str.split(','):
                part1, part2 = ion_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                ion_profile[(chain, res)] = int(deg), float(score)

        cc_profile = {}

        if cc_profile_str != '':
            for cc_profile_part in cc_profile_str.split(','):
                part1, part2 = cc_profile_part.split(':')
                chain, chain_b = part1.split('_')
                deg, score = part2.split('_')

                cc_profile[(chain, chain_b)] = int(deg), float(score)

        complex_map[pdb_id] = (comp_id, resolution, chains_str, homooligomers, lig_profile, metal_profile, ion_profile, cc_profile)

    return complex_map


def structureCheck(proteins, config):
    table = 'Structure'
    rows = ['Structure_Id', 'PDB', 'Chain']
    if len(proteins.structures) > 50:
        results = select(config, rows, table)
    elif len(proteins.structures) > 0:
        structure_ids = []
        for structure_id in proteins.structures:
            structure_ids.append(structure_id)
        results = select(config, rows, table, in_rows={'PDB': structure_ids})
    else:
        return

    stored_complexes = set()

    for row in results:
        s_id = row[0]
        pdb_id = row[1]
        chain = row[2]
        if config.verbosity >= 6:
            print(f'{pdb_id} {chain} processed in structureCheck')
        if not proteins.contains_structure(pdb_id, chain):
            continue
        proteins.set_structure_db_id(pdb_id, chain, s_id)
        proteins.set_structure_stored(pdb_id, chain, True)  # all structures, mapped or not go into this dictionary, this is important for not reinserting residues from interacting structures
        if config.verbosity >= 6:
            print(f'{pdb_id} {chain} set as stored in structureCheck')

def draw_complexes(config, proteins, stored_complexes=[], draw_all=False):
    results = select(config, ['Complex_Id', 'PDB', 'Resolution', 'Chains', 'Ligand_Profile', 'Metal_Profile', 'Ion_Profile', 'Chain_Chain_Profile', 'Homooligomers'], 'Complex')

    for row in results:
        pdb_id = row[1]
        if pdb_id not in stored_complexes and not draw_all:
            continue

        compl = complex_package.Complex(pdb_id, resolution=float(row[2]), chains_str=row[3], lig_profile_str=row[4], metal_profile_str=row[5], ion_profile_str=row[6], chain_chain_profile_str=row[7], stored=True, database_id=row[0], homomers_str=row[8])

        proteins.add_complex(pdb_id, compl)


# called by serializedPipeline
def insertStructures(structurelist, proteins, config, results=None, return_results=False):

    if len(structurelist) == 0:
        return []

    table = 'Structure'
    rows = ['Structure_Id', 'PDB', 'Chain']
    if results is None:
        results = select(config, rows, table)
        already_called = False
    else:
        already_called = True

    stored_complexes = set()
    del_list = []

    for pos, row in enumerate(results):
        pdb_id = row[1]
        chain = row[2]
        if not proteins.contains_structure(pdb_id, chain):
            if return_results:
                del_list.append(pos)
            continue
        if not already_called:
            s_id = row[0]
            proteins.set_structure_db_id(pdb_id, chain, s_id)
            proteins.set_structure_stored(pdb_id, chain, True)  # all structures, mapped or not go into this dictionary, this is important for not reinserting residues from interacting structures
            stored_complexes.add(pdb_id)
        if (pdb_id, chain) in structurelist:
            structurelist.remove((pdb_id, chain))

    if return_results:
        results = list(results)
        for pos in reversed(del_list):
            del results[pos]
        ret_results = results
    else:
        ret_results = None

    if not already_called:
        draw_complexes(config, proteins, stored_complexes=stored_complexes)

    values = []

    for (pdb_id, chain) in structurelist:
        oligos = proteins.get_oligo(pdb_id, chain)
        oligos = ''.join(oligos)
        values.append((pdb_id, chain, oligos))

    if len(values) > 0:
        insert('Structure', ['PDB', 'Chain', 'Homooligomer'], values, config)

    if len(structurelist) > 0:
        table = 'Structure'
        rows = ['Structure_Id', 'PDB', 'Chain']
        results = select(config, rows, table)

        for row in results:
            s_id = row[0]
            pdb_id = row[1]
            chain = row[2]

            if not (pdb_id, chain) in structurelist:
                continue

            if return_results:
                ret_results.append(row)

            proteins.set_structure_db_id(pdb_id, chain, s_id)

    return ret_results


# called by serializedPipeline
def insertAlignments(alignment_list, proteins, config):
    values = []
    if config.verbosity >= 2:
        t0 = time.time()
    for (u_ac, prot_id, pdb_id, chain, alignment_pir) in alignment_list:
        s_id = proteins.get_structure_db_id(pdb_id, chain)
        seq_id = proteins.get_sequence_id(u_ac, pdb_id, chain)
        coverage = proteins.get_coverage(u_ac, pdb_id, chain)
        if config.verbosity >= 6:
            print(f'Adding alignment to the database: {u_ac} {prot_id} {pdb_id} {chain} {s_id}')
        values.append((prot_id, s_id, seq_id, coverage, pack(alignment_pir)))
    if config.verbosity >= 2:
        t1 = time.time()
        print('Time for insertAlignments, part 1: ', t1 - t0)
    if len(values) > 0:
        insert('Alignment', ['Protein', 'Structure', 'Sequence_Identity', 'Coverage', 'Alignment'], values, config)
    if config.verbosity >= 2:
        t2 = time.time()
        print('Time for insertAlignments, part 2: ', t2 - t1)


def get_interfaces(protein_db_ids, config, proteins):
    if config.verbosity >= 3:    
        print('Call of get_interfaces:', len(protein_db_ids))
        t0 = time.time()

    rows = ['Protein', 'Interface_Id', 'Structure_Recommendation']
    table = 'Interface'
    results = binningSelect(protein_db_ids, rows, table, config)

    #print('In get_interfaces', results)
    if config.verbosity >= 3:
        t1 = time.time()
        print(f'Part 1 of get_interfaces: {t1-t0}')

    max_position_interface_select = 50000
    interface_map = {}
    interface_packages = []
    interface_key_package = []
    for row in results:
        protein_db_id = row[0]
        interface_id = row[1]
        interface_map[interface_id] = (protein_db_id, row[2])
        interface_key_package.append(interface_id)
        if len(interface_key_package) >= max_position_interface_select:
            interface_packages.append(interface_key_package)
            interface_key_package = []

    if config.verbosity >= 3:
        t2 = time.time()
        print(f'Part 2 of get_interfaces: {t2-t1}')
    if len(interface_key_package) > 0: 
        interface_packages.append(interface_key_package)

    results = []
    for interface_key_package in interface_packages:
        results += binningSelect(interface_key_package, ['Interface', ' Position', 'Recommended_Residue'], 'RS_Position_Interface', config)

    if config.verbosity >= 3:
        t3 = time.time()
        print(f'Part 3 of get_interfaces: {t3-t2}')
    interfaces = {}
    max_size = 500000
    position_packages = []
    positions = []
    for row in results:
        interface_db_id = row[0]
        pos_db_id = row[1]
        if not interface_db_id in interfaces:
            interfaces[interface_db_id] = {}
        interfaces[interface_db_id][pos_db_id] = row[2]
        positions.append(pos_db_id)
        if len(positions) >= max_size:
            position_packages.append(positions)
            positions = []

    if len(positions) > 0:
        position_packages.append(positions)

    if config.verbosity >= 3:
        t4 = time.time()
        print(f'Part 4 of get_interfaces: {t4-t3}')
    pos_id_map = {}
    for positions in position_packages:
        if config.verbosity >= 3:
            print(f'get_interfaces Position select: {len(positions)}')
        results = binningSelect(positions, ['Position_Id', 'Position_Number'], 'Position', config)

        for row in results:
            pos_id_map[row[0]] = row[1]

    if config.verbosity >= 3:
        t5 = time.time()
        print(f'Part 5 of get_interfaces: {t5-t4}')
    protein_interface_map = {}
    for interface_db_id in interfaces:
        protein_db_id, structure_recommendation = interface_map[interface_db_id]
        recommended_complex, chain, interacting_chain = structure_recommendation.split(',')
        protein_id = proteins.getByDbId(protein_db_id).primary_protein_id
        positions = {}
        for pos_db_id in interfaces[interface_db_id]:
            structure_tuple, res = interfaces[interface_db_id][pos_db_id].split()
            structure_id, chain = structure_tuple.split(':')
            positions[(pos_id_map[pos_db_id])] = structure_id, chain, res
        interface_obj = interface_package.Aggregated_interface(protein_id, recommended_complex = recommended_complex, chain = chain, interacting_chain = interacting_chain, positions = positions, database_id = interface_db_id)
        if not protein_id in protein_interface_map:
            protein_interface_map[protein_id] = []
        protein_interface_map[protein_id].append(interface_obj)

    #print('In get_interfaces', protein_interface_map)

    if config.verbosity >= 3:
        t6 = time.time()
        print(f'Part 6 of get_interfaces: {t6-t5}')
    for protein_id in protein_interface_map:
        if config.verbosity >= 6:
            print(f'Retrieved interfaces from database for {protein_id}:')
            for interface in protein_interface_map[protein_id]:
                interface.print_interface()
        proteins.set_aggregated_interface_map(protein_id, protein_interface_map[protein_id])

    del protein_interface_map

    if config.verbosity >= 3:
        t7 = time.time()
        print(f'Part 7 of get_interfaces: {t7-t6}')

# called by serializedPipeline
def getAlignments(proteins, config, get_all_alignments=False):
    if config.verbosity >= 2:
        t0 = time.time()

    prot_db_ids = proteins.get_stored_ids(exclude_completely_stored=(not get_all_alignments))

    for prot_id in proteins.indels:
        if proteins[prot_id].stored:
            prot_db_ids.add(proteins[prot_id].database_id)

    rows = ['Protein', 'Structure', 'Sequence_Identity', 'Coverage', 'Alignment']
    table = 'Alignment'
    results = binningSelect(prot_db_ids, rows, table, config)

    if config.verbosity >= 2:
        t1 = time.time()
        print('Time for part 1 in getAlignments:', (t1 - t0), ', number of stored proteins', len(prot_db_ids), ', number of stored alignments:', len(results))

    prot_structure_alignment_map = {}
    structure_ids = set()

    for row in results:

        prot_db_id = row[0]

        structure_id = row[1]
        seq_id = row[2]
        coverage = row[3]
        alignment = unpack(row[4])

        structure_ids.add(structure_id)

        target_seq, template_seq = sdsc_utils.process_alignment_data(alignment)

        if prot_db_id not in prot_structure_alignment_map:
            prot_structure_alignment_map[prot_db_id] = {}

        prot_structure_alignment_map[prot_db_id][structure_id] = (target_seq, template_seq, coverage, seq_id)

    if config.verbosity >= 2:
        t2 = time.time()
        print("Time for part 2 in getAlignments: %s" % (str(t2 - t1)))

    structureCheck(proteins, config)

    if config.verbosity >= 2:
        t3 = time.time()
        print("Time for part 3 in getAlignments: %s" % (str(t3 - t2)))

    stored_structure_ids = proteins.getStoredStructureIds()

    if config.verbosity >= 2:
        t4 = time.time()
        print("Time for part 4 in getAlignments: %s" % (str(t4 - t3)))

    structure_target_chain_map = {}
    for structure_db_id in stored_structure_ids:
        (structure_id, chain) = stored_structure_ids[structure_db_id]
        structure_target_chain_map[structure_id] = chain

    rows = ['PDB', 'Structure_Id', 'Chain']
    table = 'Structure'
    results = select(config, rows, table, in_rows={'PDB':structure_target_chain_map.keys()})

    interacting_structures = {}
    for row in results:
        structure_id = row[0]
        structure_db_id = row[1]
        chain = row[2]
        if chain != structure_target_chain_map[structure_id]:
            interacting_structures[structure_db_id] = (structure_id, chain)

    if config.verbosity >= 2:
        t5 = time.time()
        print("Time for part 5 in getAlignments: %s" % (str(t5 - t4)))

    interacting_protein_db_ids = []

    if config.compute_ppi:
        rows = ['Structure', 'Protein', 'Sequence_Identity', 'Coverage', 'Alignment']
        table = 'Alignment'
        results = binningSelect(interacting_structures.keys(), rows, table, config)

        for row in results:

            prot_db_id = row[1]

            structure_id = row[0]
            seq_id = row[2]
            coverage = row[3]
            alignment = unpack(row[4])

            structure_ids.add(structure_id)

            target_seq, template_seq = sdsc_utils.process_alignment_data(alignment)

            if prot_db_id not in prot_structure_alignment_map:
                prot_structure_alignment_map[prot_db_id] = {}
                interacting_protein_db_ids.append(prot_db_id)

            prot_structure_alignment_map[prot_db_id][structure_id] = (target_seq, template_seq, coverage, seq_id)
            #print('Adding interacting protein', prot_db_id, structure_id, seq_id)

    if config.verbosity >= 2:
        t6 = time.time()
        print("Time for part 6 in getAlignments: %s" % (str(t6 - t5)))

    if config.compute_ppi:
        #pos_db_map = retrieve_stored_proteins(interacting_protein_db_ids, config, proteins, without_positions=True)
        _ = retrieve_stored_proteins(interacting_protein_db_ids, config, proteins, with_struct_recs=False, return_pos_db_map=False)

    if config.verbosity >= 2:
        t7 = time.time()
        print(f'Time for part 7 in getAlignments: {str(t7 - t6)}, {len(interacting_protein_db_ids)}')

    if config.compute_ppi:
        get_interfaces(interacting_protein_db_ids, config, proteins)

    if config.verbosity >= 2:
        t8 = time.time()
        print("Time for part 8 in getAlignments: %s" % (str(t8 - t7)))

    structure_map, id_structure_map = getStructure_map(structure_ids, config)

    if config.verbosity >= 2:
        t9 = time.time()
        print("Time for part 9 in getAlignments: %s" % (str(t9 - t8)))

    pdb_ids = set()
    for pdb_id in structure_map:
        pdb_ids.add(pdb_id)

    complex_map = getComplexMap(config, pdb_ids=pdb_ids)

    if config.verbosity >= 2:
        t10 = time.time()
        print("Time for part 10 in getAlignments: %s" % (str(t10 - t9)))

    for prot_db_id in prot_structure_alignment_map:
        prot_id = proteins.getByDbId(prot_db_id).primary_protein_id
        for structure_id in prot_structure_alignment_map[prot_db_id]:
            (pdb_id, chain) = id_structure_map[structure_id]

            (target_seq, template_seq, coverage, seq_id) = prot_structure_alignment_map[prot_db_id][structure_id]

            struct_anno = structure_package.StructureAnnotation(prot_id, pdb_id, chain, alignment=(target_seq, template_seq), stored=True)
            proteins.add_annotation(prot_id, pdb_id, chain, struct_anno)

            #print('Add annotation', prot_id, pdb_id, chain)

            if not proteins.contains_structure(pdb_id, chain):
                oligo = structure_map[pdb_id][ chain][1]
                struct = structure_package.Structure(pdb_id, chain, oligo=oligo, mapped_proteins=[prot_id], database_id=structure_id)
                if structure_id in interacting_structures:
                    struct.interacting_structure = True
                proteins.add_structure(pdb_id, chain, struct)
            else:
                proteins.add_mapping_to_structure(pdb_id, chain, prot_id)

            if not proteins.contains_complex(pdb_id) and (pdb_id in complex_map):
                (comp_id, resolution, chains_str, homooligomers, lig_profile, metal_profile, ion_profile, cc_profile) = complex_map[pdb_id]
                compl = complex_package.Complex(pdb_id, resolution=resolution, chains_str=chains_str, lig_profile=lig_profile,
                                     metal_profile=metal_profile, ion_profile=ion_profile,
                                     chain_chain_profile=cc_profile, stored=True, database_id=comp_id, homomers_str=homooligomers)
                proteins.add_complex(pdb_id, compl)

            proteins.set_structure_db_id(pdb_id, chain, structure_id)

            proteins.set_coverage_by_db_id(prot_db_id, pdb_id, chain, coverage)
            proteins.set_sequence_id_by_db_id(prot_db_id, pdb_id, chain, seq_id)
            proteins.set_annotation_db_id_by_db_id(prot_db_id, pdb_id, chain, True)

    if config.verbosity >= 2:
        t11 = time.time()
        print("Time for part 11 in getAlignments: %s" % (str(t11 - t10)))


def getStructure_map(structure_ids, config):
    structure_map = {}
    id_structure_map = {}
    if len(structure_ids) > 0:
        results = binningSelect(structure_ids, ['Structure_Id', 'PDB', 'Chain', 'Homooligomer'], 'Structure', config)

        for row in results:
            s_id = row[0]

            pdb_id = row[1]
            chain = row[2]
            oligo = row[3]
            if pdb_id not in structure_map:
                structure_map[pdb_id] = {}
            structure_map[pdb_id][chain] = (s_id, oligo)
            id_structure_map[s_id] = (pdb_id, chain)
    return structure_map, id_structure_map


# called by indel_analysis
def insert_indel_results(proteins, config):
    table = 'Indel'
    columns = ['Indel_Id', 'Wildtype_Protein', 'Mutant_Protein', 'Indel_Notation', 'Analysis_Results']

    values = []
    for prot_id in proteins.indels:
        for indel_notation in proteins.indels[prot_id]:
            indel_obj = proteins.indels[prot_id][indel_notation]
            wt_prot_id = proteins.get_protein_database_id(proteins.indels[prot_id][indel_notation].wt_prot)
            mut_prot_id = proteins.get_protein_database_id(proteins.indels[prot_id][indel_notation].mut_prot)
            analysis_results = pack((indel_obj.size, indel_obj.delta_delta_classification, indel_obj.wt_aggregates, indel_obj.mut_aggregates,
                                        indel_obj.left_flank_wt_aggregates, indel_obj.left_flank_mut_aggregates,
                                        indel_obj.right_flank_wt_aggregates, indel_obj.right_flank_mut_aggregates))
            values.append((indel_obj.database_id, wt_prot_id, mut_prot_id, indel_notation, analysis_results))
    update(config, table, columns, values)
    return



def para_residue_init(rows):
    t0 = time.time()
    outs = []
    for row in rows:
        # Those residue inits include decoding of interaction profile and centrality score strings and thus takes some resources. For that a para function
        residue = residue_package.Residue(row[2], aa=row[3], lig_dist_str=row[4], chain_dist_str=row[5], RSA=row[6],
                               SSA=row[7], homo_dist_str=row[8], interaction_profile_str=row[9], centrality_score_str=row[10],
                               modres=row[11], b_factor=row[12], database_id=row[1], stored=True, phi=row[13], psi=row[14],
                               intra_ssbond=row[15], ssbond_length=row[16], intra_link=row[17], link_length=row[18],
                               cis_conformation=row[19], cis_follower=row[20], inter_chain_median_kd=row[21],
                               inter_chain_dist_weighted_kd=row[22], inter_chain_median_rsa=row[23],
                               inter_chain_dist_weighted_rsa=row[24], intra_chain_median_kd=row[25],
                               intra_chain_dist_weighted_kd=row[26], intra_chain_median_rsa=row[27],
                               intra_chain_dist_weighted_rsa=row[28])
        outs.append((row[0], row[2], residue))
    t1 = time.time()
    return outs, t1 - t0


def checkLigand(name, db, cursor):
    sql = "SELECT Ligand_Id FROM Ligand WHERE Name = '%s'" % name
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in checkLigand: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()
    if results == ():
        return 0
    elif len(results) == 0:
        return 0
    else:
        row = results[0]
        return row[0]


# called by Babel
def getLigandTemplates(name, db, cursor):
    ligand_id = checkLigand(name, db, cursor)
    sql = "SELECT Template FROM RS_Ligand_Template WHERE Ligand = '%s'" % str(ligand_id)
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getLigandTemplate: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()

    template_ids = []
    for row in results:
        template_ids.append(row[0])
    return template_ids


'''
#called by babel
def createLigandDB(outfile,session_id,db,cursor):
    sql = "SELECT Template,Mutation FROM RS_Annotation_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in createLigandDB: %s,\n%s" % (sql, f))
        db.rollback()

    template_ids = set()
    mutation_ids = set()
    for row in results:
        template_ids.add(row[0])
        mutation_ids.add(row[1])

    filtered_template_ids = set()
    if len(template_ids) > 0:
        sql = "SELECT Mutation,Template FROM RS_Mutation_Template WHERE Error IS NULL"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql, f))

        for row in results:
            if not row[1] in template_ids:
                continue
            if row[0] in mutation_ids:
                filtered_template_ids.add(row[1])

    ligand_ids = set()
    if len(filtered_template_ids) > 0:
        sql = "SELECT Ligand,Template FROM RS_Ligand_Template"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql, f))

        for row in results:
            if not row[1] in filtered_template_ids:
                continue
            ligand_ids.add(row[0])

    lines = []
    if len(ligand_ids) > 0:
        sql = "SELECT Name,Smiles,Ligand_Id FROM Ligand"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql, f))

        for row in results:
            if not row[2] in ligand_ids:
                continue
            lines.append("%s\t%s" % (row[1], row[0]))

    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()
'''


def getClassWT(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc, surface_t_id, chain_t_id, lig_t_id, config):
    if rel_sur_acc is not None:
        if rel_sur_acc > 1.0:
            rel_sur_acc = 0.01 * rel_sur_acc
    if lig_sub_dist is None:
        lig_sub_dist = "NONE"
    if chain_sub_dist is None:
        chain_sub_dist = "NONE"
    t_id = None
    try:
        mut_class = ""
        if (lig_sub_dist == "NONE" and chain_sub_dist == "NONE"):
            if rel_sur_acc is None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > config.surface_threshold:
                mut_class = "Surface isolated chain"
                t_id = surface_t_id
            else:
                mut_class = "Core isolated chain"
                t_id = surface_t_id
        elif lig_sub_dist == "NONE":
            if float(chain_sub_dist) > 5.0:
                if rel_sur_acc is None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > config.surface_threshold:
                    mut_class = "Surface"
                    t_id = surface_t_id
                else:
                    mut_class = "Core"
                    t_id = surface_t_id
            else:
                mut_class = "Contact %s" % chain_type
                t_id = chain_t_id
        elif chain_sub_dist == "NONE":
            if float(lig_sub_dist) > 5.0:
                if rel_sur_acc is None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > config.surface_threshold:
                    mut_class = "Surface"
                    t_id = surface_t_id
                else:
                    mut_class = "Core"
                    t_id = surface_t_id
            else:
                mut_class = "Contact Ligand"
                t_id = lig_t_id

        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) > 5.0):
            if rel_sur_acc is None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > config.surface_threshold:
                mut_class = "Surface"
                t_id = surface_t_id
            else:
                mut_class = "Core"
                t_id = surface_t_id
        elif (float(lig_sub_dist) < 5.0 and float(chain_sub_dist) > 5.0):
            mut_class = "Contact Ligand"
            t_id = lig_t_id
        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) < 5.0):
            mut_class = "Contact %s" % chain_type
            t_id = chain_t_id
        elif (float(lig_sub_dist) < 5.0 and float(chain_sub_dist) < 5.0):
            if (float(lig_sub_dist) <= float(chain_sub_dist)):
                mut_class = "Contact Ligand"
                t_id = lig_t_id
            elif (float(chain_sub_dist) < float(lig_sub_dist)):
                mut_class = "Contact %s" % chain_type
                t_id = chain_t_id

    except:
        mut_class = "Error"
        print("Error in class definition")
        print(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc)
        print(sys.exc_info())
    if mut_class == "":
        print(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc)

    return mut_class, t_id


def getClass(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc, config):
    if rel_sur_acc is not None:
        if rel_sur_acc > 1.0:
            rel_sur_acc = 0.01 * rel_sur_acc
    if lig_sub_dist is None:
        lig_sub_dist = "NONE"
    if chain_sub_dist is None:
        chain_sub_dist = "NONE"
    try:
        mut_class = ""
        if (lig_sub_dist == "NONE" and chain_sub_dist == "NONE"):
            if rel_sur_acc is None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > config.surface_threshold:
                mut_class = "Surface"
            else:
                mut_class = "Core"
        elif lig_sub_dist == "NONE":
            if float(chain_sub_dist) > 5.0:
                if rel_sur_acc is None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > config.surface_threshold:
                    mut_class = "Surface"
                else:
                    mut_class = "Core"
            else:
                mut_class = "Contact %s" % chain_type
        elif chain_sub_dist == "NONE":
            if float(lig_sub_dist) > 5.0:
                if rel_sur_acc is None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > config.surface_threshold:
                    mut_class = "Surface"
                else:
                    mut_class = "Core"
            else:
                mut_class = "Contact Ligand"

        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) > 5.0):
            if rel_sur_acc is None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > config.surface_threshold:
                mut_class = "Surface"
            else:
                mut_class = "Core"
        elif (float(lig_sub_dist) < 5.0 and float(chain_sub_dist) > 5.0):
            mut_class = "Contact Ligand"
        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) < 5.0):
            mut_class = "Contact %s" % chain_type
        elif (float(lig_sub_dist) < 5.0 and float(chain_sub_dist) < 5.0):
            if (float(lig_sub_dist) <= float(chain_sub_dist)):
                mut_class = "Contact Ligand"
            elif (float(chain_sub_dist) < float(lig_sub_dist)):
                mut_class = "Contact %s" % chain_type

    except:
        mut_class = "Error"
        print("Error in class definition")
        print(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc)
        print(sys.exc_info())
    if mut_class == "":
        print(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc)

    return mut_class


def majority_vote(secs):
    class_dict = {}
    for (sec, qual) in secs:
        if sec not in class_dict:
            class_dict[sec] = qual
        else:
            class_dict[sec] += qual

    max_sec = None
    max_qual = 0.0
    for sec in class_dict:
        qual = class_dict[sec]
        if qual > max_qual:
            max_qual = qual
            max_sec = sec
    return max_sec


def getChemicalDistance(aac):
    try:
        if aac.count(',') < 1:
            aa1 = aac[0]
            aa2 = aac[-1]

            chemical_distance = residue_consts.CHEM_DIST_MATRIX[aa1][aa2]
        else:
            aa1 = aac[0]
            aa2s = aac.split(",")
            aa2s[0] = aa2s[0][-1]
            chem_dists = []
            for aa2 in aa2s:
                chem_dists.append(residue_consts.CHEM_DIST_MATRIX[aa1][aa2])
            chemical_distance = float(sum(chem_dists)) / float(len(chem_dists))
    except:
        return None
    return chemical_distance


def getBlosumValue(aac):
    if aac.count(',') < 1:
        try:
            try:
                blosum_value = residue_consts.BLOSUM62[(aac[0], aac[-1])]
            except:
                blosum_value = residue_consts.BLOSUM62[(aac[-1], aac[0])]
        except:
            blosum_value = 0.0
    else:
        aa1 = aac[0]
        aa2s = aac.split(",")
        aa2s[0] = aa2s[0][-1]
        bvs = []
        for aa2 in aa2s:
            try:
                try:
                    blosum_value = residue_consts.BLOSUM62[(aa1, aa2)]
                except:
                    blosum_value = residue_consts.BLOSUM62[(aa2, aa1)]
            except:
                blosum_value = 0.0
            bvs.append(blosum_value)
        blosum_value = float(sum(bvs)) / float(len(bvs))
    return blosum_value


'''
#called by babel
def getSLD(mutation_id,template_id,db,cursor):
    sql = "SELECT Sub_Lig_Dist FROM RS_Mutation_Template WHERE Template = '%s' AND Mutation = '%s'" % (str(template_id),str(mutation_id))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getSLD")
        # Rollback in case there is any NameError
        db.rollback()
    return results[0][0]
'''

'''
#called by babel
#structure of chosen_ones: {ligand_three_letter_code:tanimoto_score}
def getLigandAnnotation(chosen_ones,session_id,distance_threshold,db,cursor):
    anno_dict = {}

    ligand_map = {}
    ligand_ids = set()
    sql = "SELECT Ligand_Id,Name FROM Ligand"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in getLigandAnnotation: %s\n%s" % (sql, f))
    for row in results:
        ligand_name = row[1]
        if ligand_name not in chosen_ones:
            continue
        ligand_map[ligand_name] = row[0]
        ligand_ids.add(row[0])

    sql = "SELECT Template,Ligand FROM RS_Ligand_Template"
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getLigandAnnotation: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()

    ligand_template_map = {}
    template_ids = set()
    for row in results:
        if not row[1] in ligand_ids:
            continue
        template_ids.add(row[0])
        if not row[1] in ligand_template_map:
            ligand_template_map[row[1]] = set([row[0]])
        else:
            ligand_template_map[row[1]].add(row[0])

    sql = "SELECT Mutation,Template FROM RS_Annotation_Session WHERE Session = '%s'" % str(session_id)
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getLigandAnnotation: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()

    template_mutation_map = {}
    for row in results:
        if not row[1] in template_ids:
            continue
        if not row[1] in template_mutation_map:
            template_mutation_map[row[1]] = {}

        template_mutation_map[row[1]][row[0]] = None

    min_t_id = min(template_ids)
    max_t_id = max(template_ids)

    t = 50000
    old_max = None
    if max_t_id - min_t_id > t:
        old_max = max_t_id
        max_t_id = min_t_id + t

    while True:
        sql = "SELECT Mutation,Template,Sub_Lig_Dist FROM RS_Mutation_Template WHERE Template BETWEEN %s AND %s" % (str(min_t_id), str(max_t_id))
        try:
            # Execute the SQL command
            cursor.execute(sql)
            results = cursor.fetchall()
            # Commit your changes in the database
            db.commit()
        except:
            raise NameError("Error in getLigandAnnotation: %s" % sql)
            # Rollback in case there is any NameError
            db.rollback()

        for row in results:
            if not row[1] in template_mutation_map:
                continue
            if not row[0] in template_mutation_map[row[1]]:
                continue
            template_mutation_map[row[1]][row[0]] = row[2]

        if old_max is None:
            break
        else:
            min_t_id = max_t_id + 1
            max_t_id = old_max
            if max_t_id - min_t_id > t:
                old_max = max_t_id
                max_t_id = min_t_id + t
            else:
                old_max = None

    for lig in chosen_ones:
        if lig not in ligand_map:
            continue
        lig_id = ligand_map[lig]
        filtered_annotation_ids = {}
        if lig_id not in ligand_template_map:
            continue
        for template_id in ligand_template_map[lig_id]:
            if template_id not in template_mutation_map:
                continue
            for mutation_id in template_mutation_map[template_id]:

                sub_lig_distances = template_mutation_map[template_id][mutation_id]
                if sub_lig_distances is not None:
                    sub_lig_distances = sub_lig_distances.split(",")
                    good_dists = []
                    for sub_lig_distance in sub_lig_distances:
                        if not sub_lig_distance == "-":
                            infos = sub_lig_distance.split(":")
                            dist = float(infos[1])
                            lig_infos = infos[0].split("_")
                            lig_name = lig_infos[0]
                            if lig_name == lig and dist <= float(distance_threshold):
                                good_dists.append(sub_lig_distance)
                    if len(good_dists) > 0:
                        good_dists = sorted(good_dists, key=lambda x: float(x.rsplit(":")[1]))
                        if template_id in filtered_annotation_ids:
                            filtered_annotation_ids[template_id].append((mutation_id, good_dists))
                        else:
                            filtered_annotation_ids[template_id] = [(mutation_id, good_dists)]
        if len(filtered_annotation_ids) > 0:
            for template_id in filtered_annotation_ids:
                filtered_annotation_ids[template_id] = sorted(filtered_annotation_ids[template_id], key=lambda x: float(x[1][0].split(":")[1]))
            anno_dict[lig] = (chosen_ones[lig], filtered_annotation_ids)
    return anno_dict
'''


def getProtIdsFromSession(session_id, config, filter_mutant_proteins=False):
    table = 'RS_Protein_Session'
    cols = ['Protein', 'Input_Id', 'Tags']
    eq_cols = {'Session': session_id}
    results = select(config, cols, table, equals_rows=eq_cols)
    prot_ids = {}
    for row in results:
        if filter_mutant_proteins and row[1] is None:
            continue
        tags = row[2]
        if tags is not None:
            tags = set(tags.split(','))
        prot_ids[row[0]] = row[1], tags
    return prot_ids

def retrieve_gene_id_map(gene_db_ids, config):
    table = 'Gene'
    columns = ['Gene_Id', 'Primary_Gene_Id', 'Gene_Name']
    results = binningSelect(gene_db_ids, columns, table, config)

    gene_id_map = {}
    for row in results:
        gene_id_map[row[0]] = (row[1], row[2])
    return gene_id_map

def retrieve_stored_proteins(prot_db_ids, config, proteins, with_mappings = False, with_struct_recs = True, without_positions = False, return_pos_db_map = True):
    if config.verbosity >= 2:
        t0 = time.time()
    cols = ['Protein_Id', 'Primary_Protein_Id', 'Sequence', 'Mutant_Type', 'Mutant_Positions', 'Gene']
    if not isinstance(prot_db_ids, dict):
        protein_ids_unknown = True
    else:
        protein_ids_unknown = False

    if protein_ids_unknown:
        results = binningSelect(prot_db_ids, cols, 'Protein', config)
    else:
        results = binningSelect(prot_db_ids.keys(), cols, 'Protein', config)

    if config.verbosity >= 2:
        t1 = time.time()
        print(f'Retrieve stored proteins, part 1: {t1 - t0} {len(prot_db_ids)}')

    id_prot_id_map = {}

    prot_id_list = set()
    prot_ids_mutants_excluded = set()
    prot_gene_db_id_map = {}
    gene_db_ids = set()

    for row in results:

        id_prot_id_map[row[0]] = row[1]
        if protein_ids_unknown:
            input_id = row[1]
            tags = None
        else:
            input_id, tags = prot_db_ids[row[0]]
        mut_type = row[3]
        if row[4] is None:
            sav_positions = None
            insertion_positions = None
            deletion_flanks = None
        else:
            sav_string, insertion_string, deletion_string = row[4].split(';')
            if sav_string == '':
                sav_positions = None
            else:
                sav_positions = [int(x) for x in sav_string.split(',')]

            insertion_positions = []
            for ins_str in insertion_string.split('!'):
                if ins_str != '':
                    insertion_positions.append([int(x) for x in ins_str.split(',')])
            if len(insertion_positions) == 0:
                insertion_positions = None

            deletion_flanks = []
            for del_str in deletion_string.split('!'):
                if del_str != '':
                    lfs, rfs = del_str.split(',')
                    if lfs == '_':
                        lf = None
                    else:
                        lf = int(lfs)
                    if rfs == '_':
                        rf = None
                    else:
                        rf = int(rfs)
                    deletion_flanks.append((lf, rf))
            if len(deletion_flanks) == 0:
                deletion_flanks = None

        if mut_type == 'SAV':
            wt_prot = row[1]
        else:
            wt_prot = None

        gene_db_id = row[5]
        if gene_db_id is not None:
            prot_gene_db_id_map[row[1]] = gene_db_id
            gene_db_ids.add(gene_db_id)

        prot_obj = protein_package.Protein(config.errorlog, primary_protein_id=row[1], database_id=row[0], input_id=input_id, sequence=row[2], wildtype_protein = wt_prot, mutant_type = mut_type, sav_positions = sav_positions, insertion_positions = insertion_positions, deletion_flanks = deletion_flanks, tags = tags)
        proteins[row[1]] = prot_obj
        prot_id_list.add(row[0])
        if mut_type is not None:
            prot_ids_mutants_excluded.add(row[0])

    if config.verbosity >= 2:
        t2 = time.time()
        print(f'Retrieve stored proteins, part 2: {t2 - t1}')

    proteins.set_stored_ids(prot_id_list, prot_ids_mutants_excluded)

    if config.verbosity >= 2:
        t3 = time.time()
        print(f'Retrieve stored proteins, part 3: {t3 - t2}')

    gene_id_map = retrieve_gene_id_map(gene_db_ids, config)

    if config.verbosity >= 2:
        t4 = time.time()
        print(f'Retrieve stored proteins, part 4: {t4 - t3}')

    for prot_id in prot_gene_db_id_map:
        proteins[prot_id].gene = gene_id_map[prot_gene_db_id_map[prot_id]][0]

    if config.verbosity >= 2:
        t5 = time.time()
        print(f'Retrieve stored proteins, part 5: {t5 - t4}')

    pos_db_map = {}

    error_count = 0
    error_message = None

    if not without_positions:
        if with_mappings:
            cols = ['Protein', 'Position_Number', 'Position_Id', 'Wildtype_Residue', 'Recommended_Structure_Data', 'IUPRED', 'IUPRED_Glob', 'Position_Data']
        elif with_struct_recs:
            cols = ['Protein', 'Position_Number', 'Position_Id', 'Wildtype_Residue', 'Recommended_Structure_Data']
        else:
            cols = ['Protein', 'Position_Number', 'Position_Id', 'Wildtype_Residue', 'Recommended_Structure_Data']
        table = 'Position'

        if protein_ids_unknown:
            results = binningSelect(prot_db_ids, cols, table, config)
        else:
            results = binningSelect(prot_db_ids.keys(), cols, table, config)

        if config.verbosity >= 2:
            t6 = time.time()
            print(f'Retrieve stored proteins, part 6: {t6 - t5} {len(results)}')

        for row in results:
            p_id = row[0]

            m_id = row[2]
            pos = row[1]

            wt_aa = row[3]

            if with_struct_recs:
                if row[4] is not None:
                    (recommended_structure_str, max_seq_structure_str) = unpack(row[4])
                else:
                    recommended_structure_str = None
                    max_seq_structure_str = None

                recommended_structure, seq_id, cov, resolution = sdsc_utils.process_recommend_structure_str(recommended_structure_str)
                max_seq_structure, max_seq_seq_id, max_seq_cov, max_seq_resolution = sdsc_utils.process_recommend_structure_str(max_seq_structure_str)



            if with_mappings:

                if row[7] is None:
                    mapping_results = None

                else:
                    try:
                        (
                            interaction_str,
                            amount_of_structures,
                            raw_structural_features,
                            raw_rin_features,
                            raw_microminer_features,
                            raw_integrated_features
                        ) = unpack(row[7])


                        mapping_results = ([], recommended_structure, max_seq_structure, None, amount_of_structures, raw_structural_features, raw_microminer_features, raw_rin_features, raw_integrated_features)

                    except:
                        error_count += 1
                        if error_message is None:
                            [e, f, g] = sys.exc_info()
                            g = traceback.format_exc()
                            error_message = f'Error in retrieving position data, couldnt unpack: {m_id}\n{e}\n{f}\n{g}'
                        continue
                mappings_obj = mappings_package.Mappings(raw_results = mapping_results)

            if with_struct_recs:
                pos_obj = position_package.Position(pos=pos, checked=True, recommended_structure=recommended_structure, database_id=m_id, wt_aa = wt_aa)
            else:
                pos_obj = position_package.Position(pos=pos, checked=True, database_id=m_id, wt_aa = wt_aa)

            if with_mappings:
                pos_obj.mappings = mappings_obj
            prot_id = id_prot_id_map[p_id]
            proteins[prot_id].add_positions([pos_obj])
            if return_pos_db_map:
                pos_db_map[m_id] = (prot_id, pos)
    else:
        t6 = time.time()

    if config.verbosity >= 2:
        t7 = time.time()
        print(f'Retrieve stored proteins, part 7: {t7 - t6}')

    if error_count > 0:
        print(f'There were {error_count} number of errors retrieving position data, example error:\n{error_message}')
    return pos_db_map

# called by output

def proteinsFromDb(session, config, with_residues=False, filter_mutant_proteins=False, with_mappings = False,
                   with_snvs=False, with_indels = False, with_multi_mutations = False, mutate_snvs=False,
                   with_alignments=False, with_complexes=False, keep_ray_alive=False, current_chunk = None, prot_db_id_dict = None):

    if config.verbosity >= 2:
        print(f'Calling proteinsFromDb, session {session}, with_residues {with_residues}, with_mappings {with_mappings}')

    if with_alignments:
        with_complexes = True

    if prot_db_id_dict is None:
        prot_db_id_dict = getProtIdsFromSession(session, config, filter_mutant_proteins=filter_mutant_proteins)

    if current_chunk is not None:
        n_of_chunks = (len(prot_db_id_dict) // config.chunksize) + 1
        small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(n_of_chunks, len(prot_db_id_dict))
        #print(config.chunksize, len(prot_db_id_dict),small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks)

        if current_chunk <= n_of_big_chunks:
            l = (current_chunk - 1) * big_chunksize
            r = current_chunk * big_chunksize
        elif current_chunk <= n_of_big_chunks + n_of_small_chunks:
            l = (((current_chunk - 1) - n_of_big_chunks) * small_chunksize) + (big_chunksize * n_of_big_chunks)
            r = ((current_chunk - n_of_big_chunks) * small_chunksize) + (big_chunksize * n_of_big_chunks)
        else:
            return None, None

        try:
            prot_db_ids = list(prot_db_id_dict.keys())
            chunked_prot_db_id_dict = {}
            for prot_db_id in prot_db_ids[l:r]:
                chunked_prot_db_id_dict[prot_db_id] = prot_db_id_dict[prot_db_id]
            complete_prot_db_id_dict = prot_db_id_dict
            prot_db_id_dict = chunked_prot_db_id_dict
        except:
            prot_db_ids = list(prot_db_id_dict)
            chunked_prot_db_id_dict = []
            for prot_db_id in prot_db_ids[l:r]:
                chunked_prot_db_id_dict.append(prot_db_id)
            complete_prot_db_id_dict = prot_db_id_dict
            prot_db_id_dict = chunked_prot_db_id_dict

    proteins = protein_package.Proteins({}, {}, {})  # create empty Proteins object

    pos_db_map = retrieve_stored_proteins(prot_db_id_dict, config, proteins, with_mappings = with_mappings)

    if with_complexes:
        draw_complexes(config, proteins, draw_all=True)

    if with_multi_mutations:
        with_snvs = True
        with_indels = True

    if with_snvs:
        cols = ['Position', 'New_AA', 'SNV_Id']
        table = 'SNV'

        results = binningSelect(pos_db_map.keys(), cols, table, config)

        snv_db_id_map = {}
        for row in results:
            snv = snv_package.SNV(row[1])
            (prot_id, pos) = pos_db_map[row[0]]
            proteins[prot_id].positions[pos].mut_aas[row[1]] = snv
            snv_db_id_map[row[2]] = (pos, snv)

    if with_indels:

        table = 'RS_Indel_Session'
        rows = ['Indel', 'Tags']
        eq_rows = {'Session': session}

        results = select(config, rows, table, equals_rows=eq_rows)
        indel_tags_map = {}
        for row in results:
            indel_tags_map[row[0]] = row[1]

        table = 'Indel'
        columns = ['Wildtype_Protein', 'Indel_Id', 'Indel_Notation', 'Mutant_Protein']
        try:
            prot_db_ids = prot_db_id_dict.keys()
        except:
            prot_db_ids = list(prot_db_id_dict)
            
        results = binningSelect(prot_db_ids, columns, table, config)

        indel_pre_map = {}
        indel_db_id_map = {}
        for row in results:
            mut_prot_db_id = row[3]

            indel_db_id = row[1]
            wt_prot_db_id = row[0]
            if wt_prot_db_id not in prot_db_id_dict:
                continue

            wt_prot = proteins.getByDbId(wt_prot_db_id).primary_protein_id

            indel_notation = row[2]

            if indel_db_id not in indel_tags_map:
                continue

            tags = indel_tags_map[indel_db_id]
            if isinstance(tags, str):
                tags = set(tags.split(','))

            if indel_notation.count('delins') == 1:
                indel = indel_package.Substitution(raw_str=indel_notation, tags=tags)

            elif indel_notation.count('del') == 1:
                indel = indel_package.Deletion(raw_str=indel_notation, tags=tags)

            elif indel_notation.count('ins') == 1:
                indel = indel_package.Insertion(raw_str=indel_notation, tags=tags)

            indel.database_id = indel_db_id

            indel_db_id_map[indel_db_id] = indel

            #if mut_prot_db_id is None:
            #    #Indels without mut_prot are part of a multi-mutation
            #    continue

            if not wt_prot in indel_pre_map:
                indel_pre_map[wt_prot] = []
            indel_pre_map[wt_prot].append(indel)

        indel_map = {}
        for wt_prot in indel_pre_map:
            uniprot.indel_insert(config, proteins, indel_map, indel_pre_map[wt_prot], wt_prot, proteins_is_object = True)

    if with_multi_mutations:
        table = 'RS_Multi_Mutation_Session'
        rows = ['Multi_Mutation', 'Tags']
        eq_rows = {'Session': session}

        results = select(config, rows, table, equals_rows=eq_rows)
        mm_tags_map = {}
        for row in results:
            mm_tags_map[row[0]] = row[1]

        table = 'Multi_Mutation'
        columns = ['Multi_Mutation_Id', 'SNVs', 'Indels', 'Wildtype_Protein', 'Mutant_Protein']
        results = binningSelect(mm_tags_map.keys(), columns, table, config)

        mm_db_map = {}
        for row in results:
            mm_db_id = row[0]
            snvs = row[1]
            indels = row[2]
            wt_prot_db_id = row[3]

            if wt_prot_db_id not in prot_db_id_dict:
                continue

            mm_tags = mm_tags_map[mm_db_id]

            mut_prot_db_id = row[4]

            wt_prot = proteins.getByDbId(wt_prot_db_id).primary_protein_id
            try:
                mut_prot = proteins.getByDbId(mut_prot_db_id).primary_protein_id
            except:
                mut_prot = '-'
            mut_list = []

            if snvs is not None:
                for snv_db_id in snvs.split(','):
                    try:
                        snv_db_id = int(snv_db_id)
                    except:
                        continue
                    pos, snv = snv_db_id_map[snv_db_id]
                    aa2 = snv.new_aa
                    position_obj = proteins[wt_prot].positions[pos]
                    mut_list.append((position_obj, aa2))
            if indels is not None:
                for indel_db_id in indels.split(','):
                    try:
                        indel_db_id = int(indel_db_id)
                    except:
                        continue
                    mut_list.append(indel_db_id_map[indel_db_id])

            if not wt_prot in indel_map:
                indel_map[wt_prot] = {}

            if config.verbosity >= 4:
                print(f'In proteinsFromDB: Call of add_multi_mutation for {wt_prot} {mut_prot}')

            proteins[wt_prot].add_multi_mutation(mut_list, indel_map[wt_prot], mm_tags = mm_tags, mut_prot_id = mut_prot)

            if wt_prot not in mm_db_map:
                mm_db_map[wt_prot] = {}
            mm_db_map[wt_prot][mut_prot] = mm_db_id

        multi_mutation_objects = []
        for primary_protein_id in list(proteins.protein_map.keys()).copy():  # iterate over snapshot, since the mutant proteins get added to proteins in the loop
            if proteins[primary_protein_id].multi_mutations is not None:
                # creates mutation protein object
                multi_mutation_objects += proteins[primary_protein_id].create_multi_mutations(proteins, config)

        for multi_mutation_obj in multi_mutation_objects:
            if multi_mutation_obj.wt_prot not in proteins.multi_mutations:
                proteins.multi_mutations[multi_mutation_obj.wt_prot] = []
            proteins.multi_mutations[multi_mutation_obj.wt_prot].append(multi_mutation_obj)
            mm_db_id = mm_db_map[multi_mutation_obj.wt_prot][multi_mutation_obj.mut_prot]
            proteins.multi_mutation_back_map[mm_db_id] = multi_mutation_obj

    ray_utils.ray_init(config)

    if mutate_snvs:
        prot_ids = list(proteins.get_protein_ids()).copy()
        for prot_id in prot_ids:
            proteins[prot_id].mutate_snvs(proteins, config)

        serializedPipeline.autoTemplateSelection(config, proteins)

    if with_alignments:
        No_Errors = serializedPipeline.paraAlignment(config, proteins, skip_inserts=True, get_all_alignments=True)

    if not keep_ray_alive:
        ray.shutdown()

    if with_residues:
        getStoredResidues(proteins, config)

    if current_chunk is not None:
        return proteins, complete_prot_db_id_dict

    return proteins


# method for comparing/developing the new classifiction
def diffSurfs(mutation_surface_dict, g_u_dict, mutation_dict, outfile, config):
    class_dict = {}
    for m in mutation_surface_dict:
        (u_ac, u_id) = g_u_dict[mutation_dict[m][1]]
        aac = mutation_dict[m][0]
        DSC = 0.0  # decision sum core
        DSS = 0.0  # decision sum surface
        n = 0.0
        min_surf = 2.0
        lines = ["Coverage\trASA"]
        for (surface, qual, cov) in mutation_surface_dict[m]:
            if surface < config.surface_threshold:
                DSC += qual * (cov**5)
            else:
                DSS += qual * (cov**5)
            lines.append("%s\t%s" % (str(cov), str(surface)))
            n += 1.0

            if surface < min_surf:
                min_surf = surface

        weighted_surface_value = DSS - 2 * DSC
        if weighted_surface_value > 0:
            weighted_sc = "Surface"
            if min_surf < config.surface_threshold:
                f = open("%s.%s_%s.tsv" % (outfile, u_ac, aac), 'w')
                f.write("\n".join(lines))
                f.close()
        else:
            weighted_sc = "Core"

        conf_sc = (1.0 - 1.0 / (n + 1.0)) * abs(DSS - 2 * DSC) / (DSS + 2 * DSC)


def createInterDict(mutation_inter_dict, chain_type='sc'):

    inter_dict = {}

    for m_id in mutation_inter_dict:
        profiles = mutation_inter_dict[m_id]
        total_qual = 0.0
        ion_qual = 0.0
        metal_qual = 0.0
        lig_qual = 0.0
        chain_qual = 0.0
        average_profile = [0.0] * 14
        for profile_str, qual in profiles:
            if None in (profile_str, qual):
                continue
            profile = rin.Interaction_profile(profile_str=profile_str)
            total_qual += qual

            Ion_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'ion')
            Ion_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'ion')
            Metal_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'metal')
            Metal_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'metal')
            Ligand_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'ligand')
            Ligand_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'ligand')
            Chain_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'interchain')
            Chain_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'interchain')
            Short_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'neighbor')
            Short_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'neighbor')
            Medium_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'short')
            Medium_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'short')
            Long_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'long')
            Long_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'long')

            if Ligand_Interaction_Degree > 0:
                lig_qual += qual
            if Metal_Interaction_Degree > 0:
                metal_qual += qual
            if Ion_Interaction_Degree > 0:
                ion_qual += qual
            if Chain_Interaction_Degree > 0:
                chain_qual += qual

            average_profile[0] += Ion_Interaction_Degree * qual
            average_profile[1] += Ion_Interaction_Score * qual
            average_profile[2] += Metal_Interaction_Degree * qual
            average_profile[3] += Metal_Interaction_Score * qual
            average_profile[4] += Ligand_Interaction_Degree * qual
            average_profile[5] += Ligand_Interaction_Score * qual
            average_profile[6] += Chain_Interaction_Degree * qual
            average_profile[7] += Chain_Interaction_Score * qual
            average_profile[8] += Short_Interaction_Degree * qual
            average_profile[9] += Short_Interaction_Score * qual
            average_profile[10] += Medium_Interaction_Degree * qual
            average_profile[11] += Medium_Interaction_Score * qual
            average_profile[12] += Long_Interaction_Degree * qual
            average_profile[13] += Long_Interaction_Score * qual

        if total_qual > 0.0:
            average_profile[8] = average_profile[8] / total_qual
            average_profile[9] = average_profile[9] / total_qual
            average_profile[10] = average_profile[10] / total_qual
            average_profile[11] = average_profile[11] / total_qual
            average_profile[12] = average_profile[12] / total_qual
            average_profile[13] = average_profile[13] / total_qual

        if ion_qual > 0.0:
            average_profile[0] = average_profile[0] / ion_qual
            average_profile[1] = average_profile[1] / ion_qual
        if metal_qual > 0.0:
            average_profile[2] = average_profile[2] / metal_qual
            average_profile[3] = average_profile[3] / metal_qual
        if lig_qual > 0.0:
            average_profile[4] = average_profile[4] / lig_qual
            average_profile[5] = average_profile[5] / lig_qual
        if chain_qual > 0.0:
            average_profile[6] = average_profile[6] / chain_qual
            average_profile[7] = average_profile[7] / chain_qual
        inter_dict[m_id] = average_profile

    return inter_dict


def excludeFarClasses(c, sc):
    if c == "Surface" or c == "Core" or c == 'Disorder' or c is None:
        return c

    interactions = re.sub(r' Interaction$', '', re.sub(r'^[^:]*: ', '', c)).split(' and ')
    non_far_interactions = [x for x in interactions if ' far' not in x]

    if (len(non_far_interactions)) == 0:
        return sc

    clas = " and ".join(non_far_interactions)
    if len(non_far_interactions) == 4:
        clas = "Quadruple Interaction: " + clas
    elif len(non_far_interactions) == 3:
        clas = "Triple Interaction: " + clas
    elif len(non_far_interactions) == 2:
        clas = "Double Interaction: " + clas
    elif len(non_far_interactions) == 1:
        clas = clas + " Interaction"

    return clas


def writeInterFile(outfile, inter_dict, mutation_dict, protein_dict, new_aa_map, tag_map, class_dict, header=True):
    startline = "Uniprot\tAAC\tSpecie\tTag\tLigand_Interaction_Degree\tLigand_Interaction_Score\tChain_Interaction_Degree\tChain_Interaction_Score\tShort_Interaction_Degree\tShort_Interaction_Score\tMedium_Interaction_Degree\tMedium_Interaction_Score\tLong_Interaction_Degree\tLong_Interaction_Score\tClass\tComplex class"
    if header:
        lines = [startline]
    else:
        lines = []
    for m in inter_dict:
        aac, Protein_Id = mutation_dict[m][0:2]

        new_aa = new_aa_map[m]
        aac = "%s%s" % (aac.split(',')[0], new_aa)
        #(u_ac,u_id,species) = g_u_dict[mutation_dict[m][1]]
        (u_ac, gpan, u_id, error_code, error, species, input_id) = protein_dict[Protein_Id]

        (Class, conf, weighted_sc, conf_sc, best_res, max_seq_res, amount_of_structures,
         weighted_c, conf_c,
         weighted_d, conf_d,
         weighted_r, conf_r,
         weighted_l, conf_l,
         weighted_m, conf_m,
         weighted_i, conf_i,
         weighted_h, conf_h,
         max_seq_id,
         weighted_raw, weighted_cent, weighted_norm,
         weighted_lig_degree, weighted_lig_score,
         weighted_metal_degree, weighted_metal_score,
         weighted_ion_degree, weighted_ion_score,
         weighted_prot_degree, weighted_prot_score,
         weighted_rna_degree, weighted_rna_score,
         weighted_dna_degree, weighted_dna_score,
         weighted_modres, modres_prop, b_factor,
         intra_ssbond_prop, inter_ssbond_prop,
         intra_link_prop, inter_link_prop,
         cis_prop, cis_follower_prop,
         weighted_inter_chain_median_kd, weighted_inter_chain_dist_weighted_kd,
         weighted_inter_chain_median_rsa, weighted_inter_chain_dist_weighted_rsa,
         weighted_intra_chain_median_kd, weighted_intra_chain_dist_weighted_kd,
         weighted_intra_chain_median_rsa, weighted_intra_chain_dist_weighted_rsa) = class_dict[m]
        simple_class = simplifyClass(Class, weighted_sc)  # BUG: undefined variable

        if m in inter_dict:
            interstr = '\t'.join([str(x) for x in inter_dict[m]])
        else:
            interstr = '\t'.join(([''] * 14))

        line = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac, aac, species, tag_map[m], interstr, simple_class, Class)
        lines.append(line)

    f = open(outfile, 'a')
    f.write("\n".join(lines))
    f.close()


def writeClassFile(outfile, mutation_surface_dict, mutation_sec_dict, mutation_dict, protein_dict, class_dict, tag_map, header=True):
    startline = "Uniprot-Ac\tUniprot Id\tRefseq\tPDB-ID (Input)\tResidue-Id\tAmino Acid\tPosition\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tIndividual Interactions\tConfidence Value\tSecondary Structure\tRecommended Structure\tSequence-ID\tCoverage\tResolution\tMax Seq Id Structure\tMax Sequence-ID\tMax Seq Id Coverage\tMax Seq Id Resolution\tAmount of mapped structures"

    if header:
        lines = [startline]
    else:
        lines = []
    for m in class_dict:
        if m in mutation_sec_dict:
            mv_sec_ass = majority_vote(mutation_sec_dict[m])
        else:
            mv_sec_ass = None

        (Class, conf, weighted_sc, conf_sc, best_res, max_seq_res, amount_of_structures,
         weighted_c, conf_c,
         weighted_d, conf_d,
         weighted_r, conf_r,
         weighted_l, conf_l,
         weighted_m, conf_m,
         weighted_i, conf_i,
         weighted_h, conf_h,
         max_seq_id,
         weighted_raw, weighted_cent, weighted_norm,
         weighted_lig_degree, weighted_lig_score,
         weighted_metal_degree, weighted_metal_score,
         weighted_ion_degree, weighted_ion_score,
         weighted_prot_degree, weighted_prot_score,
         weighted_rna_degree, weighted_rna_score,
         weighted_dna_degree, weighted_dna_score,
         weighted_modres, modres_prop, b_factor,
         intra_ssbond_prop, inter_ssbond_prop,
         intra_link_prop, inter_link_prop,
         cis_prop, cis_follower_prop,
         weighted_inter_chain_median_kd, weighted_inter_chain_dist_weighted_kd,
         weighted_inter_chain_median_rsa, weighted_inter_chain_dist_weighted_rsa,
         weighted_intra_chain_median_kd, weighted_intra_chain_dist_weighted_kd,
         weighted_intra_chain_median_rsa, weighted_intra_chain_dist_weighted_rsa) = class_dict[m]
        simple_class = simplifyClass(Class, weighted_sc)  # BUG: undefined variable

        if best_res is not None:
            [r_id, qual, res_aa, res_nr, pdb_id, chain, resolution, cov, seq_id, rsa, min_lig, min_metal, min_ion, iacs] = best_res

            recommended_structure = '%s:%s %s:%s' % (pdb_id, chain, res_nr, res_aa)
        else:
            resolution = '-'
            cov = '-'
            seq_id = '-'
            recommended_structure = '-'

        if max_seq_res is not None:
            [max_seq_r_id, max_seq_qual, max_seq_res_aa, max_seq_res_nr, max_seq_pdb_id, max_seq_chain, max_seq_resolution, max_seq_cov, max_seq_seq_id, max_seq_rsa, max_min_lig, max_min_metal, max_min_ion, max_iacs] = max_seq_res

            max_seq_structure = '%s:%s %s:%s' % (max_seq_pdb_id, max_seq_chain, max_seq_res_nr, max_seq_res_aa)
        else:
            max_seq_resolution = '-'
            max_seq_cov = '-'
            max_seq_seq_id = '-'
            max_seq_structure = '-'

        aac, Protein_Id = mutation_dict[m][0:2]
        input_res_id = mutation_dict[m][4]
        if input_res_id is None:
            input_res_id = ''

        (u_ac, gpan, u_id, error_code, error, species, input_id) = protein_dict[Protein_Id]

        input_pdb_id = ''
        if len(u_ac) == 6 and u_ac[4] == ':':
            input_pdb_id = u_ac
            u_ac = ''
        interaction_str = '-'  # TODO

        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac, u_id, gpan, input_pdb_id, input_res_id, aac[0], aac[1:], species, tag_map[m], weighted_sc, Class, simple_class, interaction_str, str(conf), mv_sec_ass, recommended_structure, str(seq_id), str(cov), str(resolution), max_seq_structure, str(max_seq_seq_id), str(max_seq_cov), str(max_seq_resolution), str(amount_of_structures)))
    f = open(outfile, 'a')
    f.write("\n".join(lines))
    f.close()

# called by output, needs update, currently no gene_score


def goTermAnalysis(session_id, outfile, db, cursor):
    sql = "SELECT Protein FROM RS_Protein_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)
    datasetsize = len(results)
    go_dict = {}
    gene_dict = {}
    Protein_Id_list = set([])
    for row in results:
        gene_score = row[1]
        if gene_score is None:
            gene_score = 0.0
        gene_dict[row[0]] = gene_score
        Protein_Id_list.add(row[0])

    sql = "SELECT GO_Term,Protein FROM RS_Protein_GO_Term"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)
    go_term_id_list = set([])
    for row in results:
        if not row[1] in Protein_Id_list:
            continue
        gene_score = gene_dict[row[1]]
        if not row[0] in go_dict:
            go_term_id_list.add(row[0])
            go_dict[row[0]] = ["", "", gene_score, 1.0, gene_score]
        else:
            i = go_dict[row[0]][3] + 1
            avg_score = go_dict[row[0]][2]
            go_dict[row[0]][2] = ((i - 1) / i) * avg_score + (1 / i) * gene_score
            go_dict[row[0]][3] = i
            go_dict[row[0]][4] += gene_score

    sql = "SELECT GO_Term_Id,Name,Id FROM GO_Term"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)

    for row in results:
        if not row[0] in go_term_id_list:
            continue
        go_dict[row[0]][0] = row[1]
        go_dict[row[0]][1] = row[2]

    go_list = sorted(list(go_dict.values()), key=lambda x: x[4], reverse=True)
    lines = ["GO-Term\tGO-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
    for go in go_list:
        normalized_score = go[4] / datasetsize
        lines.append("%s\t%s\t%s\t%s\t%s\t%s" % (str(go[0]), str(go[1]), str(go[4]), str(go[2]), str(int(go[3])), str(normalized_score)))
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()


# called by output, needs update, currently no gene_score
def pathwayAnalysis(session_id, outfile, db, cursor):
    sql = "SELECT Gene,Gene_Score FROM RS_Protein_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)
    datasetsize = len(results)
    path_dict = {}
    gene_dict = {}
    Protein_Id_list = set([])
    for row in results:
        gene_score = row[1]
        if gene_score is None:
            gene_score = 0.0
        gene_dict[row[0]] = gene_score
        Protein_Id_list.add(row[0])

    sql = "SELECT Pathway,Gene FROM RS_Protein_Pathway"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)
    pathway_id_list = set([])
    for row in results:
        if not row[1] in Protein_Id_list:
            continue
        gene_score = gene_dict[row[1]]
        if not row[0] in path_dict:
            pathway_id_list.add(row[0])
            path_dict[row[0]] = ["", "", gene_score, 1.0, gene_score]
        else:
            i = path_dict[row[0]][3] + 1
            avg_score = path_dict[row[0]][2]
            path_dict[row[0]][2] = ((i - 1) / i) * avg_score + (1 / i) * gene_score
            path_dict[row[0]][3] = i
            path_dict[row[0]][4] += gene_score

    if len(pathway_id_list) == 0:
        lines = ["Pathway\tReactome-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
        page = "\n".join(lines)
        f = open(outfile, "wb")
        f.write(page)
        f.close()
        return

    sql = "SELECT Pathway_Id,Name,Reactome_Id FROM Pathway"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)

    for row in results:
        if not row[0] in pathway_id_list:
            continue
        path_dict[row[0]][0] = row[1]
        path_dict[row[0]][1] = row[2]

    path_list = sorted(list(path_dict.values()), key=lambda x: x[4], reverse=True)
    lines = ["Pathway\tReactome-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
    for path in path_list:
        normalized_score = path[4] / datasetsize
        lines.append("%s\t%s\t%s\t%s\t%s\t%s" % (str(path[0]), str(path[1]), str(path[4]), str(path[2]), str(int(path[3])), str(normalized_score)))
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()


# called by postAnnoAnno
def updateAnnoAnno(anno_anno_map, session_id, db, cursor):
    valuestrs = []
    for (m_id1, m_id2) in anno_anno_map:
        (min_d, atom, atom2, t_id1, t_id2, chain1, chain2) = anno_anno_map[(m_id1, m_id2)]
        valuestrs.append("('%s','%s','%s','%s','%s','%s','%s','%s','%s %s')" % (str(t_id1), str(t_id2), str(m_id1), str(m_id2), chain1, chain2, str(session_id), str(min_d), str(atom), str(atom2)))
    sql = "INSERT IGNORE INTO RS_Annotation_Annotation(Template_1,Template_2,Mutation_1,Mutation_2,Chain_1,Chain_2,Session,Distance,Atompair) VALUES %s" % ','.join(valuestrs)
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in updateAnnoAnno: %s" % (f))


# called by calcAnnoRate
def calculateAnnotationRate(db, cursor, session_id):
    sql = "SELECT Position FROM RS_Position_Session WHERE Session = %s" % session_id
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in calculateAnnotationRate: %s" % sql)
    mut_id_list = set()
    for x in results:
        mut_id_list.add(str(x[0]))

    if len(mut_id_list) > 0:
        sql = "SELECT Position_Id,Gene FROM Position WHERE Position_Id in (%s)" % ",".join(mut_id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        gene_map = {}
        Protein_Id_list = set()
        for row in results:
            Protein_Id = str(row[1])
            mut_id = str(row[0])
            gene_map[mut_id] = Protein_Id
            Protein_Id_list.add(Protein_Id)

        sql = "SELECT Protein_Id,Error_Code FROM Protein WHERE Protein_Id in (%s)" % ",".join(Protein_Id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        error_map = {}
        for row in results:
            error_map[str(row[0])] = str(row[1])

        sql = "SELECT Position,Template,Error_Code FROM RS_Position_Template WHERE Position in (%s)" % ",".join(mut_id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        anno_mut_id_list = set()
        all_error_muts = set()
        fine_templates = set()

        for row in results:
            mut_id = str(row[0])
            if row[2] is None:
                anno_mut_id_list.add(mut_id)
                all_error_muts.discard(mut_id)
                fine_templates.add(row[1])
            elif mut_id not in anno_mut_id_list:
                all_error_muts.add(mut_id)

        sql = "SELECT Gene,Sequence_Identity,Template_Id FROM Template WHERE Gene in (%s)" % ",".join(Protein_Id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        number_of_templates = len(results)
        greater90 = set()
        greater95 = set()
        greater98 = set()
        greater99 = set()
        for row in results:
            if row[2] in fine_templates:
                g_id = str(row[0])
                s_id = float(row[1])
                if s_id > 99.0:
                    greater99.add(g_id)
                if s_id > 98.0:
                    greater98.add(g_id)
                if s_id > 95.0:
                    greater95.add(g_id)
                if s_id > 90.0:
                    greater90.add(g_id)

        muts_without_template = set()
        muts_with_unknown_gene_error = set()
        muts_with_anno = set()
        muts_with_temp_98 = 0
        for m in mut_id_list:
            if m in anno_mut_id_list:
                muts_with_anno.add(m)
                if gene_map[m] in greater98:
                    muts_with_temp_98 += 1
            if error_map[gene_map[m]] == '3':
                muts_without_template.add(m)
            elif error_map[gene_map[m]] == '4':
                muts_with_unknown_gene_error.add(m)

        return (len(mut_id_list), len(muts_with_anno), len(muts_without_template), len(muts_with_unknown_gene_error), len(all_error_muts), float(len(muts_with_anno)) / float(len(mut_id_list)), number_of_templates, len(greater90), len(greater95), len(greater98), len(greater99), muts_with_temp_98)
    else:
        return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
