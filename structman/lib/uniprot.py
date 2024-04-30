import os
import socket
import sys
import traceback
import time
import urllib.error
import urllib.parse
import urllib.request
import requests
import json

from Bio import Entrez

from structman.lib.database.database_core_functions import select, insert, update
from structman.lib.sdsc.sdsc_utils import translate
from structman.lib.sdsc import sdsc_utils
from structman.lib.sdsc import protein as protein_package
from structman.lib.sdsc import gene as gene_package
from structman.lib.pdbParser import standardParsePDB
from structman.lib.templateFiltering import parsePDB
from structman.lib.sdsc.consts import codons
from structman.base_utils import base_utils
from structman.settings import HUMAN_INFO_MAP

POLLING_INTERVAL = 3

def is_connected():
    try:
        # connect to the host -- tells us if the host is actually
        # reachable
        socket.create_connection(("1.1.1.1", 53))
        return True
    except OSError:
        pass
    return False


def connection_sleep_cycle(verbosity):
    while not is_connected():
        if verbosity >= 1:
            print('No connection, sleeping a bit and then try again')
        time.sleep(30)


def getUniprotId(query, querytype, verbosity=0):
    url = 'https://rest.uniprot.org/idmapping/run'
    params = {
        'from': '%s' % (querytype),
        'to': 'ID',
        'format': 'tab',
        'query': '%s' % (query)
    }
    connection_sleep_cycle(verbosity)
    data = urllib.parse.urlencode(params).encode('utf-8')
    request = urllib.request.Request(url, data)
    contact = ""  # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    try:
        response = urllib.request.urlopen(request, timeout=60)
    except:
        return None
    page = response.read(200000).decode('utf-8')
    try:
        lines = page.split("\n")
        line = lines[1].split()
        uniprot_id = line[1]
    # If unusable result, try without version number
    except:
        query = query.split(".")
        if len(query) > 1:
            query = query[0]
            return getUniprotId(query, querytype)
        else:
            return None

    return uniprot_id


def tag_update(tag_map, u_ac, new_entry):
    if u_ac not in tag_map:
        tag_map[u_ac] = new_entry
        return tag_map
    else:
        old_entry = tag_map[u_ac]

    for aac in new_entry:
        if aac not in old_entry:
            tag_map[u_ac][aac] = new_entry[aac]
        else:
            tags = set([])
            for tag in old_entry[aac].split(','):
                tags.add(tag)
            for tag in new_entry[aac].split(','):
                tags.add(tag)
            tag_map[u_ac][aac] = ','.join(tags)
    return tag_map

# deactive, needs update
def updateMappingDatabase(u_acs, db, config):
    cursor = db.cursor()
    ac_id_values = []
    ac_ref_values = []
    seq_values = []
    for u_ac in u_acs:
        seq_out = getSequence(u_ac, config, return_id=True)
        if seq_out is None:
            continue

        seq, refseqs, go_terms, pathways, u_id = seq_out
        if u_id is None:  # This can happen for uniprot entries, which got deleted from uniprot
            print("Warning: Uniprot entry:", u_ac, " not found, most probably the entry got deleted from uniprot")
            continue
        ac_id_values.append("('%s','%s','%s','%s')" % (u_ac, u_ac[-2:], u_id, u_id[:2]))
        for refseq in refseqs.split(','):
            ac_ref_values.append("('%s','%s','%s','%s')" % (u_ac, u_ac.split('-')[0][-2:], refseq, refseq[:2]))
        seq_values.append("('%s','%s')" % (u_ac, seq))

    try:

        if len(ac_id_values) > 0:

            sql = "INSERT IGNORE INTO AC_ID(Uniprot_Ac,Ac_Hash,Uniprot_Id,Id_Hash) VALUES %s " % (','.join(ac_id_values))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql, f))

        if len(ac_ref_values) > 0:

            sql = "INSERT IGNORE INTO AC_Refseq(Uniprot_Ac,Ac_Hash,Refseq,Refseq_Hash) VALUES %s" % (','.join(ac_ref_values))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql, f))

        if len(seq_values) > 0:

            sql = "INSERT IGNORE INTO Sequences(Uniprot_Ac,Sequence) VALUES %s " % (','.join(seq_values))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql, f))
    except:
        # this happens for Devs without insert rights to the mapping DB, just ignore for the moment
        pass


# called by serializedPipeline
def u_ac_isoform_search(gene_sequence_map, stems, ref_stem_map, config):
    isoform_map = get_all_isoforms(stems, config)
    isoform_specific_id_map = {}
    for nm_ref in gene_sequence_map:
        if nm_ref not in ref_stem_map:
            continue
        seq = gene_sequence_map[nm_ref]
        u_ac_stem = ref_stem_map[nm_ref]
        if u_ac_stem not in isoform_map:
            continue
        for u_ac in isoform_map[u_ac_stem]:
            if seq == isoform_map[u_ac_stem][u_ac][0]:
                isoform_specific_id_map[nm_ref] = u_ac
                break
    return isoform_specific_id_map

# called by serializedPipeline


def integrate_protein(config, proteins, genes, indel_map, primary_protein_id, input_id, prot_map, u_ac=None, u_id=None, ref=None, pdb_id = None, other_ids={}, is_pdb_input = False, gene_id = None, protein_specific_tags = None):

    if protein_specific_tags is None:
        protein_specific_tags = set([])
    if gene_id is not None:
        if gene_id not in genes:
            gene_obj = gene_package.Gene(gene_id, proteins = [primary_protein_id])
            genes[gene_id] = gene_obj
        else:
            genes[gene_id].proteins.add(primary_protein_id)

    if primary_protein_id not in proteins:
        protein = protein_package.Protein(config.errorlog, primary_protein_id=primary_protein_id, u_ac=u_ac, u_id=u_id, ref_id=ref, positions=prot_map[input_id][0], input_id=input_id, pdb_id = pdb_id, gene = gene_id, tags = protein_specific_tags)
        proteins[primary_protein_id] = protein
    else:
        proteins[primary_protein_id].u_id = u_id
        if is_pdb_input:
            proteins[primary_protein_id].add_residues(prot_map[input_id][0])
        else:
            proteins[primary_protein_id].add_positions(prot_map[input_id][0])

        proteins[primary_protein_id].tags = proteins[primary_protein_id].tags | protein_specific_tags

    for other_id_type in other_ids:
        proteins[primary_protein_id].add_other_ids(other_id_type, other_ids[other_id_type])

    if len(prot_map[input_id][1]) > 0:
        indel_insert(config, proteins, indel_map, prot_map[input_id][1], primary_protein_id)

    for multi_mutation in prot_map[input_id][2]:
        if len(multi_mutation) > 1:
            if primary_protein_id in indel_map:
                indel_map_entry = indel_map[primary_protein_id]
            else:
                indel_map_entry = {}
            proteins[primary_protein_id].add_multi_mutation(multi_mutation, indel_map_entry)
    return


def indel_insert(config, proteins, indel_map, indels, ac, proteins_is_object = False, create_mut_protein_object = True):
    def proteins_member_check(proteins, member, proteins_is_object):
        if proteins_is_object:
            return proteins.contains(member)
        else:
            return member in proteins

    if ac not in indel_map:
        indel_map[ac] = {}
    for indel in indels:
        if indel.get_notation() in indel_map[ac]:
            continue
        indel.wt_prot = ac
        if create_mut_protein_object:
            indel_protein_name = indel.create_protein_name(ac)
            if proteins_member_check(proteins, indel_protein_name, proteins_is_object):
                broken = False
                while proteins_member_check(proteins, indel_protein_name, proteins_is_object):
                    indel_protein_name = indel.create_other_protein_name(ac, indel_protein_name)
                    if indel_protein_name is None:  # Duplicate entry
                        #config.errorlog.add_error('All indel protein names are reserved, duplicate indels in input?')
                        broken = True
                        break
                if broken:
                    continue
            insertion_positions, deletion_flanks = indel.get_mutant_positions()
            if insertion_positions is not None:
                insertion_positions = [insertion_positions]
            if deletion_flanks is not None:
                deletion_flanks = [deletion_flanks]
            indel_mut_protein = protein_package.Protein(config.errorlog, primary_protein_id=indel_protein_name, mutant_type = 'indel', insertion_positions = insertion_positions, deletion_flanks = deletion_flanks)
            proteins[indel_protein_name] = indel_mut_protein
            indel.set_proteins(ac, indel_protein_name)
        indel_map[ac][indel.get_notation()] = indel
    return

# called by serializedPipeline


def IdMapping(config, ac_map, id_map, np_map, pdb_map, hgnc_map, nm_map, ensembl_map, prot_gene_map, prot_tags_map):

    if config.verbosity >= 4:
        print(f'prot_gene_map at the start of IdMapping: {prot_gene_map}')

    proteins = {}
    indel_map = {}
    genes = {}

    for ac in ac_map:
        if ac in prot_gene_map:
            gene_id = prot_gene_map[ac]
        else:
            gene_id = None
        integrate_protein(config, proteins, genes, indel_map, ac, ac, ac_map, u_ac=ac, gene_id = gene_id, protein_specific_tags = prot_tags_map[ac])

    # Step one: map everything to uniprot-ac
    if len(id_map) > 0:
        if config.mapping_db_is_set:

            results = select(config, ['Uniprot_Ac', 'Uniprot_Id'], 'UNIPROT', in_rows={'Uniprot_Id':id_map.keys()}, from_mapping_db = True)

            stored_ids = set()
            for row in results:
                u_ac = row[0]
                u_id = row[1]
                stored_ids.add(u_id)
                if u_id in prot_gene_map:
                    gene_id = prot_gene_map[u_id]
                else:
                    gene_id = u_id
                integrate_protein(config, proteins, genes, indel_map, u_ac, u_id, id_map, u_ac=u_ac, u_id=u_id, gene_id = gene_id, protein_specific_tags = prot_tags_map[u_id])

            unstored_ids = []
            for u_id in id_map:
                if u_id not in stored_ids:
                    unstored_ids.append(u_id)
            update_acs = []
            if len(unstored_ids) > 0:  # whenever ids are left in the dict, go to uniprot and download all unmapped entries (this happens for newer uniprot entries, which are not yet in the local mapping database)

                # This part is identical to the part, when no local database is used
                id_ac_map = getUniprotIds(config, unstored_ids, 'UniProtKB_AC-ID', target_sub_type = 'primaryAccession')
                if not id_ac_map is None:
                    for u_id in id_ac_map:
                        u_ac = id_ac_map[u_id]
                        if u_id in prot_gene_map:
                            gene_id = prot_gene_map[u_id]
                        else:
                            gene_id = u_id

                        integrate_protein(config, proteins, genes, indel_map, u_ac, u_id, id_map, u_ac=u_ac, u_id=u_id, gene_id = gene_id, protein_specific_tags = prot_tags_map[u_id])

                        # This part is different
                        if not sdsc_utils.is_mutant_ac(u_ac):
                            update_acs.append(u_ac)
                    # updateMappingDatabase(update_acs,db,config)

        else:
            id_ac_map = getUniprotIds(config, list(id_map.keys()), 'UniProtKB_AC-ID', target_sub_type = 'primaryAccession')
            for u_id in id_ac_map:
                u_ac = id_ac_map[u_id]
                if u_id in prot_gene_map:
                    gene_id = prot_gene_map[u_id]
                else:
                    gene_id = u_id
                integrate_protein(config, proteins, genes, indel_map, u_ac, u_id, id_map, u_ac=u_ac, u_id=u_id, gene_id = gene_id, protein_specific_tags = prot_tags_map[u_id])

    if len(np_map) > 0:
        if config.mapping_db_is_set:
            results = select(config, ['Uniprot_Ac', 'Refseq'], 'UNIPROT', in_rows={'Refseq': np_map.keys()}, from_mapping_db=True)

            ref_u_ac_map = {}  # different u_acs may have the same refseq, try to choose the right one, prefering u_acs containing '-'
            gene_id_snap = set(proteins.keys())  # snapshot of u_acs not originating from refseq-mapping

            stored_refs = set()

            for row in results:
                u_ac = row[0]
                ref = row[1]
                stored_refs.add(ref)
                if ref in prot_gene_map:
                    gene_id = prot_gene_map[ref]
                else:
                    gene_id = None
                if ref not in ref_u_ac_map:
                    ref_u_ac_map[ref] = u_ac

                    integrate_protein(config, proteins, genes, indel_map, ref, ref, np_map, u_ac=u_ac, gene_id = gene_id, protein_specific_tags = prot_tags_map[ref])

                elif u_ac in gene_id_snap:
                    if ref_u_ac_map[ref].count('-') == 0 and u_ac.count('-') > 0:
                        ref_u_ac_map[ref] = u_ac

                    integrate_protein(config, proteins, genes, indel_map, ref, ref, np_map, u_ac=u_ac, gene_id = gene_id, protein_specific_tags = prot_tags_map[ref])

                elif ref_u_ac_map[ref].count('-') == 0:
                    if u_ac.count('-') > 0:
                        # if the current u_ac does not contain a '-' and the new found u_ac contains a '-': replace the corresponding ids
                        old_ac = ref_u_ac_map[ref]
                        ref_u_ac_map[ref] = u_ac
                        del proteins[old_ac]
                        integrate_protein(config, proteins, genes, indel_map, ref, ref, np_map, u_ac=u_ac, gene_id = gene_id, protein_specific_tags = prot_tags_map[ref])
            # similar to uniprot-id mapping, we have to go to uniprot to get search for unstored refseq entries and if we find them, we have to update the local mapping database
            unstored_refs = []
            for ref in np_map:
                if ref not in stored_refs:
                    unstored_refs.append(ref)
            update_acs = []
            if len(unstored_refs) > 0:
                np_ac_map = getUniprotIds(config, unstored_refs, 'RefSeq_Protein', target_sub_type = 'primaryAccession')
                for ref in np_ac_map:
                    u_ac = np_ac_map[ref]
                    if ref in prot_gene_map:
                        gene_id = prot_gene_map[ref]
                    else:
                        gene_id = None
                    integrate_protein(config, proteins, genes, indel_map, ref, ref, np_map, u_ac=u_ac, gene_id = gene_id, protein_specific_tags = prot_tags_map[ref])
                    if not sdsc_utils.is_mutant_ac(u_ac):
                        update_acs.append(u_ac)
                # updateMappingDatabase(update_acs,db,config)

        else:
            np_ac_map = getUniprotIds(config, list(np_map.keys()), 'RefSeq_Protein', target_sub_type = 'primaryAccession')
            for ref in np_ac_map:
                u_ac = np_ac_map[ref]
                if ref in prot_gene_map:
                    gene_id = prot_gene_map[ref]
                else:
                    gene_id = None
                integrate_protein(config, proteins, genes, indel_map, ref, ref, np_map, u_ac=u_ac, gene_id = gene_id, protein_specific_tags = prot_tags_map[ref])

    if len(nm_map) > 0:
        nm_ac_map = getUniprotIds(config, list(nm_map.keys()), 'RefSeq_Nucleotide', target_sub_type = 'primaryAccession')
        for ref in nm_map:
            if ref in prot_gene_map:
                gene_id = prot_gene_map[ref]
            else:
                gene_id = None
            if ref in nm_ac_map:
                u_ac = nm_ac_map[ref]
                integrate_protein(config, proteins, genes, indel_map, ref, ref, nm_map, u_ac=u_ac, gene_id = gene_id, protein_specific_tags = prot_tags_map[ref])
            else:
                integrate_protein(config, proteins, genes, indel_map, ref, ref, nm_map, gene_id = gene_id, protein_specific_tags = prot_tags_map[ref])

    if len(hgnc_map) > 0:  # No support for mapping DB yet
        hgnc_ac_map = getUniprotIds(config, list(hgnc_map.keys()), 'HGNC', target_sub_type = 'primaryAccession')
        for hgnc in hgnc_ac_map:
            if hgnc in prot_gene_map:
                gene_id = prot_gene_map[hgnc]
            else:
                gene_id = None
            u_ac = hgnc_ac_map[hgnc]
            integrate_protein(config, proteins, genes, indel_map, u_ac, hgnc, hgnc_map, u_ac=u_ac, other_ids={'HGNC_ID': hgnc}, gene_id = gene_id, protein_specific_tags = prot_tags_map[hgnc])

    if len(ensembl_map) > 0:
        #Uniprot is doesnt seem to work currently for embl id mapping, thus skip for now
        #ensembl_ac_map = getUniprotIds(config, list(ensembl_map.keys()), 'EMBL', target_sub_type = 'primaryAccession')
        ensembl_ac_map = {}
        for ensembl_id in ensembl_map:

            #if ensembl_id in prot_gene_map:
            #    gene_id = prot_gene_map[ensembl_id]
            #else:
            gene_id = "to_fetch"

            if ensembl_id in ensembl_ac_map:
                u_ac = ensembl_ac_map[ensembl_id]
                integrate_protein(config, proteins, genes, indel_map, ensembl_id, ensembl_id, ensembl_map, other_ids={'ENSEMBL_ID': ensembl_id}, gene_id = gene_id, u_ac=u_ac, protein_specific_tags = prot_tags_map[ensembl_id])
            else:
                integrate_protein(config, proteins, genes, indel_map, ensembl_id, ensembl_id, ensembl_map, other_ids={'ENSEMBL_ID': ensembl_id}, gene_id = gene_id, protein_specific_tags = prot_tags_map[ensembl_id])


    # Step two: get uniprot-id and refseqs from uniprot-ac

    ac_iso_map = {}
    id_search = []
    for primary_protein_id in proteins:
        if proteins[primary_protein_id].u_id is not None:
            continue
        if primary_protein_id.count('_') > 0:
            continue
        split = primary_protein_id.split('-')
        if len(split) == 2:
            base, iso = split
            if base not in ac_iso_map:
                ac_iso_map[base] = [iso]
            else:
                ac_iso_map[base].append(iso)

            id_search.append(base)
        else:
            id_search.append(primary_protein_id)

    if len(id_search) > 0:
        if config.mapping_db_is_set:
            stored_u_acs = set()

            results = select(config, ['Uniprot_Ac', 'Uniprot_Id'], 'UNIPROT', in_rows={'Uniprot_Ac': id_search}, from_mapping_db=True)

            for row in results:
                u_ac = row[0]
                u_id = row[1]

                stored_u_acs.add(u_ac)
                if u_ac in proteins:
                    proteins[u_ac].u_id = u_id

                if u_ac in ac_iso_map:
                    for iso in ac_iso_map[u_ac]:
                        proteins['%s-%s' % (u_ac, iso)].u_id = u_id
                        stored_u_acs.add('%s-%s' % (u_ac, iso))

            unstored_u_acs = []
            for u_ac in id_search:
                if u_ac not in stored_u_acs:
                    if not sdsc_utils.is_mutant_ac(u_ac):
                        unstored_u_acs.append(u_ac)

            if len(unstored_u_acs) > 0:
                # updateMappingDatabase(unstored_u_acs, db, config)

                results = select(config, ['Uniprot_Ac', 'Uniprot_Id'], 'UNIPROT', in_rows={'Uniprot_Ac': unstored_u_acs}, from_mapping_db=True)

                for row in results:
                    u_ac = row[0]
                    u_id = row[1]

                    if u_ac in proteins:
                        proteins[u_ac].u_id = u_id
                    if u_ac in ac_iso_map:
                        for iso in ac_iso_map[u_ac]:
                            proteins['%s-%s' % (u_ac, iso)].u_id = u_id

        else:
            ac_id_map = getUniprotIds(config, id_search, 'UniProtKB_AC-ID', target_sub_type = 'uniProtkbId')
            if ac_id_map is not None:
                for u_ac in ac_id_map:
                    u_id = ac_id_map[u_ac]
                    if u_ac in proteins:
                        proteins[u_ac].u_id = u_id
                    if u_ac in ac_iso_map:
                        for iso in ac_iso_map[u_ac]:
                            proteins['%s-%s' % (u_ac, iso)].u_id = u_id

    if config.mapping_db_is_set:

        results = select(config, ['Uniprot_Ac', 'Refseq'], 'UNIPROT', in_rows={'Uniprot_Ac': proteins.keys()}, from_mapping_db=True)

        for row in results:
            u_ac = row[0]
            ref = row[1]
            proteins[u_ac].ref_id = ref
    else:
        ac_np_map = getUniprotIds(config, list(proteins.keys()), 'UniProtKB_AC-ID', target_type="RefSeq_Protein")
        for u_ac in ac_np_map:
            ref = ac_np_map[u_ac]
            proteins[u_ac].ref_id = ref

    expanded_pdbs = {}
    for pdb_tuple in pdb_map:
        if pdb_tuple in prot_gene_map:
            gene_id = prot_gene_map[pdb_tuple]
        else:
            gene_id = None
        if pdb_tuple[-1] == '-':
            page, _ = standardParsePDB(pdb_tuple[:4], config.pdb_path)
            if page == '':
                continue

            (_, _, _, _, _, _, _, chain_type_map, _, _, _, _, _, _, _) = parsePDB(page)

            chains = []
            for chain_id in chain_type_map:
                if chain_type_map[chain_id] == 'Protein':
                    chains.append(chain_id)

            if pdb_tuple in indel_map:
                for chain in chains:
                    new_pdb_tuple = f'{pdb_tuple[:4]}:{chain}'
                    indel_map[new_pdb_tuple] = indel_map[pdb_tuple]
                del indel_map[pdb_tuple]
            for chain in chains:
                new_pdb_tuple = f'{pdb_tuple[:4]}:{chain}'   
                expanded_pdbs[new_pdb_tuple] = pdb_map[pdb_tuple]
        else:
            integrate_protein(config, proteins, genes, indel_map, pdb_tuple, pdb_tuple, pdb_map, is_pdb_input = True, pdb_id = pdb_tuple, gene_id = gene_id, protein_specific_tags = prot_tags_map[pdb_tuple])

    for new_pdb_tuple in expanded_pdbs:
        integrate_protein(config, proteins, genes, indel_map, new_pdb_tuple, new_pdb_tuple, expanded_pdbs, is_pdb_input = True, pdb_id = new_pdb_tuple, protein_specific_tags = prot_tags_map[new_pdb_tuple[:4]])

    return proteins, indel_map, genes


def getUniprotIds(config, query_ids, querytype, target_type = "UniProtKB", target_sub_type = 'uniProtkbId', timeout = 60):
    if query_ids is None:
        config.errorlog.add_error(f'Query_ids is None, {querytype} {target_type}')
        return {}

    if len(query_ids) == 0:
        return {}

    if target_type != "UniProtKB":
        target_sub_type = None

    if config.verbosity >= 4:
        print('Call of getUniprotIds:', query_ids, querytype, target_type)

    query = ' '.join(query_ids)
    url = 'https://rest.uniprot.org/idmapping/run'
    connection_sleep_cycle(config.verbosity)
    params = {
        'from': '%s' % (querytype),
        'to': '%s' % (target_type),
        'ids': '%s' % (query)
    }
    request = requests.post(url, params)

    try:
        request.raise_for_status()
        job_id = request.json()['jobId']
        while True:
            request = requests.get(f"https://rest.uniprot.org/idmapping/status/{job_id}")
            request.raise_for_status()
            job = request.json()
            if "jobStatus" in job:
                if job["jobStatus"] == "RUNNING":
                    print(f"Retrying in {POLLING_INTERVAL}s")
                    time.sleep(POLLING_INTERVAL)
                else:
                    raise Exception(job["jobStatus"])
            else:
                results = job
                break
    except:
        e, f, g = sys.exc_info()
        config.errorlog.add_warning(f"Uniprot did not answer for {url}:\n{params}\n{e}\n{f}\n{g}")
        return {}

    if config.verbosity >= 4:
        print('Uniprot answer:\n',results)

    results_list = results['results']

    uniprot_ids = {}
    try:
        error = True
        for uniprot_entry_dict in results_list:
            error = False
            quer = uniprot_entry_dict['from']
            if quer in uniprot_ids:
                continue
            if target_sub_type is not None:
                target = uniprot_entry_dict['to'][target_sub_type]
            else:
                target = uniprot_entry_dict['to']

            if quer not in query_ids:
                error = True
                break
            uniprot_ids[quer] = target

        if error:

            if len(query_ids) == 1:
                quer = query_ids.pop()
                if querytype == 'UniProtKB_AC-ID' and quer.count('-') == 1:
                    int_quer = set([quer.split('-')[0]])
                    intermediate_map = getUniprotIds(config, int_quer, querytype, target_type=target_type)
                    return {quer: intermediate_map[quer.split('-')[0]]}
                elif target_sub_type == 'primaryAccession' and querytype == 'UniProtKB_AC-ID':
                    config.errorlog.add_error(f'Couldnt map Uniprot ID to AC: {quer}\nUniprot answer:\n{results}')
                    return {}
                else:
                    return {quer: None}
            # try a divide and conqer solution to isolate the problematic query
            set_A = set()
            set_B = set()
            while len(query_ids) > 0:
                set_A.add(query_ids.pop())
                if len(query_ids) > 0:
                    set_B.add(query_ids.pop())
            map_A = getUniprotIds(config, set_A, querytype, target_type=target_type)
            map_B = getUniprotIds(config, set_B, querytype, target_type=target_type)
            return map_A.update(map_B)

    # If unusable result, try without version number
    except:
        # This seems not up-to-date
        '''
        if query.count(".") > 0:
            ids = query.split("\t")
            nids = []
            for i in ids:
                ni = i.split(".")[0]
                nids.append(ni)
            nquery = "\t".join(nids)
            if nquery != query:
                return getUniprotId(nquery,querytype,target_type=target_type)
            else:
                return {}
        '''
        config.errorlog.add_warning(f'Parse error of uniprot ID-mapping webservice\n{results}')
        return {}
    return uniprot_ids


def retrieve_transcript_metadata(server, transcript_ids):
    rest_url = f'{server}/lookup/id'
    transcript_ids = list(transcript_ids)
    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    transcript_gene_id_dict = {}
    chunksize = 800
    chunk = transcript_ids[:chunksize]
    i = 0
    try_number = 0
    gene_ids = set()
    while len(chunk) > 0:
        data = json.dumps({'ids':chunk})

        try:
            r = requests.post(rest_url, headers = headers, data = data)
        except:
            try_number += 1
            if try_number == 5:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                print(f'POST request failed: {e}\n{f}\n{g}')
                break
            time.sleep(try_number**2)
            continue

        if not r.ok:
            r.raise_for_status()
            print(f'Transcript metadata Retrieval failed')
            return None

        decoded = r.json()
        try_number = 0

        for transcript_id in decoded:
            try:
                gene_id = decoded[transcript_id]['Parent']
                transcript_gene_id_dict[transcript_id] = gene_id
                gene_ids.add(gene_id)
            except:
                #print(transcript_id)
                continue
                
        i += 1
        chunk = transcript_ids[(chunksize*i):(chunksize*(i+1))]

    gene_name_dict = retrieve_gene_metadata(server, gene_ids)
    for transcript_id in transcript_gene_id_dict:
        gene_id = transcript_gene_id_dict[transcript_id]
        try:
            gene_name = gene_name_dict[gene_id]
        except:
            gene_name = None
        transcript_gene_id_dict[transcript_id] = (gene_id, gene_name)

    return transcript_gene_id_dict


def retrieve_gene_metadata(server, gene_ids, expand = False):
    rest_url = f'{server}/lookup/id'
    gene_ids = list(gene_ids)
    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    gene_name_dict = {}
    if expand:
        gene_transcript_dict = {}
    chunksize = 950
    if expand:
        chunksize = 50
    chunk = gene_ids[:chunksize]
    i = 0
    try_number = 0
    while len(chunk) > 0:
        data = json.dumps({'ids':chunk})

        try:
            r = requests.post(rest_url, headers = headers, data = data, params = {'expand' : int(expand)})
        except:
            try_number += 1
            if try_number == 5:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                print(f'POST request failed: {e}\n{f}\n{g}')
                break
            time.sleep(try_number**2)
            continue

        if not r.ok:
            r.raise_for_status()
            print(f'Transcript metadata Retrieval failed')
            return None

        decoded = r.json()
        try_number = 0

        for gene_id in decoded:
            try:
                gene_name = decoded[gene_id]['display_name']
                gene_name_dict[gene_id] = gene_name
                if expand:
                    gene_transcript_dict[gene_id] = []
                    for transcript_entry in decoded[gene_id]['Transcript']:
                        transcript_id = transcript_entry['id']
                        gene_transcript_dict[gene_id].append(transcript_id)
            except:
                continue
                
        i += 1
        chunk = gene_ids[(chunksize*i):(chunksize*(i+1))]

    if expand:
        return gene_name_dict, gene_transcript_dict
    return gene_name_dict


def retrieve_transcript_sequences(server, transcript_ids, recursed = False):
    rest_url = f'{server}/sequence/id'

    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    transcript_seq_dict = {}
    chunksize = 45
    transcript_ids = list(transcript_ids)
    chunk = transcript_ids[:chunksize]
    i = 0
    try_number = 0
    empty_seqs = {}
    while len(chunk) > 0:

        data = json.dumps({'ids':chunk})

        try:
            r = requests.post(rest_url, headers = headers, data = data, params = {'type' : 'protein'})
        except:
            try_number += 1
            if try_number == 5:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                print(f'POST request failed: {e}\n{f}\n{g}')
                break
            time.sleep(try_number**2)
            continue

        if not r.ok:
            print(f'Transcript Sequence Retrieval failed: {data}')
            if recursed:
                return {}, {transcript_ids[0]:None}
            for transcript_id in chunk:
                seq_sub_dict, sub_empty_seqs = retrieve_transcript_sequences(server, [transcript_id], recursed=True)
                transcript_seq_dict.update(seq_sub_dict)
                empty_seqs.update(sub_empty_seqs)
            i += 1
            chunk = transcript_ids[(chunksize*i):(chunksize*(i+1))]
            continue
            #r.raise_for_status()
            #return None

        decoded = r.json()
        try_number = 0

        for entry in decoded:
            transcript_id = entry['query']
            nt_seq = entry['seq']
            transcript_seq_dict[transcript_id] = nt_seq

        i += 1
        chunk = transcript_ids[(chunksize*i):(chunksize*(i+1))]

    return transcript_seq_dict, empty_seqs

def embl_database_lookup(config, transcript_ids):
    missing_ids = set(transcript_ids)
    results = select(config, ['Transcript_ID', 'Gene_ID', 'Name', 'Sequence'], 'EMBL_Transcripts', in_rows={'Transcript_ID': missing_ids}, from_mapping_db=True)
    sequence_map = {}
    gene_id_map = {}
    broken_ids = set()
    mean_seq_len = 0
    for row in results:
        transcript_id = row[0]
        gene_id = row[1]
        gene_name = row[2]
        try:
            seq = base_utils.unpack(row[3])
        except:
            broken_ids.add(transcript_id)
            continue

        if seq is None:
            missing_ids.remove(transcript_id)
            continue

        sequence_map[transcript_id] = seq

        mean_seq_len += len(seq)
        gene_id_map[transcript_id] = (gene_id, gene_name)
        missing_ids.remove(transcript_id)
    if len(sequence_map) > 0:
        mean_seq_len = mean_seq_len / len(sequence_map)
    if config.verbosity >= 2:
        print(f'Embl database lookup: {len(sequence_map)} sequences, mean length: {mean_seq_len}')
    return sequence_map, missing_ids, gene_id_map, broken_ids

def embl_database_update(config, sequence_map, transcript_id_gene_id_map, broken_ids, new_no_seqs):
    values = []
    for transcript_id in new_no_seqs:
        if transcript_id in broken_ids:
            continue
        try:
            gene_id, gene_name = transcript_id_gene_id_map[transcript_id]
        except:
            gene_id = None
            gene_name = None

        values.append((transcript_id, gene_id, gene_name, base_utils.pack(None)))

    for transcript_id in sequence_map:
        if transcript_id in broken_ids:
            continue
        try:
            gene_id, gene_name = transcript_id_gene_id_map[transcript_id]
        except:
            gene_id = None
            gene_name = None

        values.append((transcript_id, gene_id, gene_name, base_utils.pack(sequence_map[transcript_id])))

    insert('EMBL_Transcripts', ['Transcript_ID', 'Gene_ID', 'Name', 'Sequence'], values, config, mapping_db = True)

    values = []
    for transcript_id in broken_ids:
        try:
            gene_id, gene_name = transcript_id_gene_id_map[transcript_id]
        except:
            gene_id = None
            gene_name = None
        values.append((transcript_id, gene_id, gene_name, base_utils.pack(sequence_map[transcript_id])))

    update(config, 'EMBL_Transcripts', ['Transcript_ID', 'Gene_ID', 'Name', 'Sequence'], values, mapping_db = True)

def get_ensembl_seqs(config, full_transcript_ids):
    aa_seqs = {}
    total_gene_id_map = {}
    chunksize = 5000
    iter_number = 0
    while True:
        l = iter_number*chunksize
        r = (iter_number+1)*chunksize
        transcript_ids = full_transcript_ids[l:r]

        if config.verbosity >= 2:
            print(f'Get ensembl transcript sequences, iteration: {iter_number}, chunksize: {len(transcript_ids)}')

        if config.mapping_db_is_set:
            if config.verbosity >= 2:
                print(f'Calling embl_database_lookup with {len(transcript_ids)} number of transcript ids')
            sequence_map, missing_ids, gene_id_map, broken_ids = embl_database_lookup(config, transcript_ids)
            if config.verbosity >= 2:
                print(f'embl_database_lookup returned with {len(sequence_map)} number of sequences and {len(missing_ids)} number of still missing transcripts')
        else:
            sequence_map = {}
            gene_id_map = {}
            missing_ids = transcript_ids

        if len(missing_ids) > 0:
            fasta_sequences, new_no_seqs = retrieve_transcript_sequences(config.ensembl_server, missing_ids)
            transcript_id_gene_id_map = retrieve_transcript_metadata(config.ensembl_server, missing_ids)
            if config.mapping_db_is_set:
                try:
                    embl_database_update(config, fasta_sequences, transcript_id_gene_id_map, broken_ids, new_no_seqs)
                except:
                    [e, f, g] = sys.exc_info()
                    g = traceback.format_exc()
                    if config.verbosity >= 3:
                        print(f'Updating EMBL transcipt database failed:\n{e}\n{f}\n{g}')
        else:
            fasta_sequences = {}
            transcript_id_gene_id_map = {}

        fasta_sequences.update(sequence_map)
        aa_seqs.update(fasta_sequences)

        total_gene_id_map.update(gene_id_map)
        total_gene_id_map.update(transcript_id_gene_id_map)

        if r >= len(full_transcript_ids):
            break
        iter_number += 1
    return aa_seqs, total_gene_id_map

# First version of refseq sequence retrieval. Uses biopython.
# called by serializePipeline


def get_refseq_sequences(refseqs, config, seq_type='nucleotide'):

    """
    Fetches sequence data from NCBI refseq database.

    Parameters
    ----------

    refseqs
        comma-separated string of refseq identifiers

    config
        instance of structman config object

    seq_type
        type of sequence to be retrieved, either 'nucleotide' or 'protein'

    Returns
    -------

    seq_map
        dictionary of {refseq_identifier:sequence}
    """

    Entrez.email = config.user_mail

    ret_type = 'fasta_cds_aa'
    if seq_type == 'protein':
        ret_type = 'fasta'

    net_handle = Entrez.efetch(
        db=seq_type, id=refseqs, rettype=ret_type, retmode="text"
    )
    page = net_handle.read()
    net_handle.close()
    right_split = '_prot'
    left_split = '|'
    if seq_type == 'protein':
        right_split = ' '
        left_split = None
    seq_map = sdsc_utils.parseFasta(page=page, left_split=left_split, right_split=right_split)

    return seq_map


def translateGNSMap(gene_nuc_sequence_map):
    gene_sequence_map = {}
    for ref in gene_nuc_sequence_map:
        nuc_seq = gene_nuc_sequence_map[ref]
        aa_seq = translate(nuc_seq)
        gene_sequence_map[ref] = aa_seq
    return gene_sequence_map


def get_all_isoforms(u_ac_stems, config):

    could_find_more = True
    isoform_counter = 1
    isoform_counter_steps = 20
    isoform_map = {}
    while could_find_more:
        u_acs = []
        last_isoforms = {}
        for u_ac_stem in u_ac_stems:
            isoform_map[u_ac_stem] = {}
            stem_isoform_counter = isoform_counter
            for i in range(isoform_counter_steps):
                if stem_isoform_counter == 1:
                    u_acs.append(u_ac_stem)
                else:
                    isoform = '%s-%s' % (u_ac_stem, str(stem_isoform_counter))
                    u_acs.append(isoform)
                    if stem_isoform_counter == ((isoform_counter + isoform_counter_steps) - 1):
                        last_isoforms[u_ac_stem] = isoform
                stem_isoform_counter += 1
        isoform_counter += isoform_counter_steps

        gene_sequence_map = getSequencesPlain(u_acs, config, save_errors=False, skip_missing_routine=True)
        for u_ac in gene_sequence_map:
            u_ac_stem = u_ac.split('-')[0]
            isoform_map[u_ac_stem][u_ac] = gene_sequence_map[u_ac]
        could_find_more = False
        for u_ac_stem in last_isoforms:
            last_isoform = last_isoforms[u_ac_stem]
            if last_isoform in isoform_map[u_ac_stem]:
                could_find_more = True
            else:
                u_ac_stems.remove(u_ac_stem)

    return isoform_map

# called by serializedPipeline


def getSequencesPlain(u_acs, config, max_seq_len=None, filtering_db=None, save_errors=True, skip_missing_routine=False):
    gene_sequence_map = {}
    filtered_set = set()

    missing_set = set()

    if filtering_db is not None:
        filtering_db_path, dbs = filtering_db

    if config.mapping_db_is_set:
        if len(u_acs) == 0:
            return {}
        if config.verbosity >= 2:
            t0 = time.time()

        results = select(config, ['Uniprot_Ac', 'Sequence'], 'UNIPROT', in_rows={'Uniprot_Ac': u_acs}, from_mapping_db = True)

        if config.verbosity >= 2:
            t1 = time.time()
            print("getSequencesPlain Part 1: ", str(t1 - t0))

        n = 0
        for row in results:
            u_ac = row[0]
            if u_ac not in u_acs:
                continue
            if row[1] is None:
                config.errorlog.add_warning('Sequence field is None in Mapping DB, for: %s' % u_ac) 
                #gene_sequence_map[u_ac] = None, None, None
                continue
            try:
                seq = base_utils.unpack(row[1])
            except:
                config.errorlog.add_warning('Sequence field is defect in Mapping DB, for: %s' % u_ac)
                #gene_sequence_map[u_ac] = None, None, None
                continue

            if max_seq_len is not None:
                if len(seq) > max_seq_len:
                    filtered_set.add(u_ac)
                    n += 1
                    continue
            
            disorder_scores = None
            disorder_regions_datastruct = None

            gene_sequence_map[u_ac] = seq, disorder_scores, disorder_regions_datastruct

        if n > 0 and config.verbosity >= 2:
            print('Filtered ', n, ' Sequences due to max length: ', max_seq_len)

        if config.verbosity >= 2:
            t2 = time.time()
            print("getSequencesPlain Part 2: ", str(t2 - t1))

    else:
        missing_set = set(u_acs)

    if config.verbosity >= 2:
        t2 = time.time()

    in_db = set()
    for u_ac in u_acs:
        if filtering_db is not None:
            inside_all = True
            for db_name in dbs:
                folder_key = u_ac.split('-')[0][-2:]
                filename = '%s/%s/%s_%s_gpw.fasta.gz' % (filtering_db_path, folder_key, u_ac, db_name)
                if not os.path.isfile(filename):
                    inside_all = False
            if inside_all:
                in_db.add(u_ac)
                continue

        if u_ac in filtered_set:
            continue

        if u_ac not in gene_sequence_map:
            missing_set.add(u_ac)

    if config.verbosity >= 2:
        t3 = time.time()
        print("getSequencesPlain Part 3: ", str(t3 - t2))

    if not skip_missing_routine:

        if config.verbosity >= 2:
            print('Size of missing set: ', len(missing_set))
        if config.verbosity >= 3:
            print(missing_set)

        for u_ac in missing_set:
            seq_out = getSequence(u_ac, config)
            if seq_out is None and save_errors:
                config.errorlog.add_warning('getSequence output is None for %s' % u_ac)
                gene_sequence_map[u_ac] = 0, None, None
                continue
            elif seq_out is None:
                continue
            seq, refseqs, go_terms, pathways = seq_out
            gene_sequence_map[u_ac] = seq, None, None

        if config.verbosity >= 2:
            t4 = time.time()
            print("getSequencesPlain Part 4: ", str(t4 - t3))

    if filtering_db is not None:
        return gene_sequence_map, in_db
    return gene_sequence_map


# called by serializePipeline
def getSequences(proteins, config):

    t0 = time.time()
    u_acs = set()
    iso_map = {}
    protein_map = proteins.get_protein_map()
    for prot_id in protein_map:
        if protein_map[prot_id].sequence is None:
            u_acs.add(prot_id)

        if prot_id.count('-') == 1:
            [base, iso] = prot_id.split('-')
        else:
            base = prot_id
            iso = 'c'
        if base not in iso_map:
            iso_map[base] = [iso]
        else:
            iso_map[base].append(iso)

    gene_sequence_map = getSequencesPlain(u_acs, config)

    for u_ac in u_acs:
        protein_map[u_ac].sequence = gene_sequence_map[u_ac][0]
        proteins.set_disorder_scores(u_ac, gene_sequence_map[u_ac][1])
        proteins.set_disorder_regions(u_ac, gene_sequence_map[u_ac][2])

    info_map_path = HUMAN_INFO_MAP

    f = open(info_map_path, 'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        words = line[:-1].split('\t')
        u_ac = words[0]
        if u_ac in iso_map:
            go_terms = {}
            pathways = {}
            go_str = words[1]
            if go_str != '':
                go_doubles = go_str[:-1].split(';')
                for go_double in go_doubles:
                    go_id = (':'.join(go_double.split(':')[:2])).strip()
                    go_term = ':'.join(go_double.split(':')[2:])
                    go_terms[go_id] = go_term
            path_str = words[2]
            if path_str != '':
                path_doubles = path_str[:-1].split(';')
                for path_double in path_doubles:
                    [reac_id, pathway] = path_double.split(':', 1)
                    pathways[reac_id] = pathway
            for iso in iso_map[u_ac]:
                if iso == 'c':
                    iso_u_ac = u_ac
                else:
                    iso_u_ac = '%s-%s' % (u_ac, iso)
                protein_map[iso_u_ac].go_terms = go_terms
                protein_map[iso_u_ac].pathways = pathways

def get_last_version(u_ac):
    url = 'https://www.uniprot.org/uniprot/%s?version=*' % u_ac
    #connection_sleep_cycle(config.verbosity)
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        page = response.read(9000000).decode('utf-8')
    except:
        e, f, g = sys.exc_info()
        config.errorlog.add_error('Error trying to reach: %s\n%s\n%s' % (url, str(e), str(f)))
        return None
    print(page)

def get_obsolete_sequence(u_ac, config, return_id=False, tries = 0):
    uniparc_id_map = getUniprotIds(config, [u_ac], 'UniProtKB_AC', target_type="UniParc", timeout = 60 * (tries + 1))
    if not u_ac in uniparc_id_map:
        if tries < 3:
             return get_obsolete_sequence(u_ac, config, return_id=return_id, tries = (tries + 1))
        config.errorlog.add_warning(f'Couldnt find uniparc id for: {u_ac}, {uniparc_id_map}')
        return None
    uniparc_id = uniparc_id_map[u_ac]
    if config.verbosity >= 3:
        print('Uniparc ID of potential obsolete uniprot entry (', u_ac, ') is:', uniparc_id)
    if uniparc_id is None:
        config.errorlog.add_warning(f'Uniparc id was None for: {u_ac}, {uniparc_id_map}')
        return None
    return getSequence(uniparc_id, config, return_id=return_id, obsolete_try = True)

def getSequence(uniprot_ac, config, tries=0, return_id=False, obsolete_try = False):
    if config.verbosity >= 3:
        print('uniprot.getSequence for ', uniprot_ac)

    if uniprot_ac is None:
        config.errorlog.add_error(f'Uniprot Ac is None')
        return None

    if sdsc_utils.is_mutant_ac(uniprot_ac):
        config.errorlog.add_error('Cannot call getSequence with a mutant protein: %s' % uniprot_ac)
        return None

    # new part just for the sequence
    if len(uniprot_ac) < 2:
        return None

    if not uniprot_ac[0:3] == 'UPI':
        url = 'https://www.uniprot.org/uniprot/%s.fasta' % uniprot_ac
    else:
        url = 'https://www.uniprot.org/uniparc/%s.fasta' % uniprot_ac
    connection_sleep_cycle(config.verbosity)
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request, timeout=(tries + 1) * 10)
        page = response.read(9000000).decode('utf-8')
    except:
        if tries < 3:
            return getSequence(uniprot_ac, config, tries=tries + 1, return_id=return_id)
        else:
            return get_obsolete_sequence(uniprot_ac, config, return_id=return_id)

    if config.verbosity >= 3 and obsolete_try:
        print('Obsolete entry try:', uniprot_ac, url)
        if config.verbosity >= 4:
            print('Online return:\n', page)

    lines = page.split("\n")

    wildtype_sequences = []

    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            continue
        wildtype_sequences.append(line)

    wildtype_sequence = ("".join(wildtype_sequences)).replace(" ", "").replace("\n", "")

    if wildtype_sequence == '':
        if obsolete_try: #Unless it could come to an endless loop
            return None
        if config.verbosity >= 3:
            print('Try to search for obsolete uniprot entry:', uniprot_ac)
        return get_obsolete_sequence(uniprot_ac, config, return_id=return_id)

    if uniprot_ac[0:3] == 'UPI':
        if return_id:
            return (wildtype_sequence, {}, {}, {}, None)
        return (wildtype_sequence, {}, {}, {})

    # old part, now just for refseqs,go and reactome

    url = 'https://www.uniprot.org/uniprot/%s.txt' % uniprot_ac
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request, timeout=(tries + 1) * 10)
        page = response.read(9000000).decode('utf-8')
    except:
        if tries < 3:
            return getSequence(uniprot_ac, config, tries=tries + 1, return_id=return_id)

        else:
            return None

    lines = page.split("\n")

    refseqs = {}
    go_terms = {}
    pathways = {}
    u_id = None
    for line in lines:
        if line == '':
            continue
        words = line.split()
        if len(words) == 0:
            print(uniprot_ac)
            return None
        if words[0] == 'ID':
            u_id = words[1]
        if words[0] == "DR":
            if len(words) > 1:
                if words[1] == "GO;":
                    words_ = line.split(";")
                    go_id = words_[1].strip(";").strip()
                    go_name = words_[2].strip(";")
                    go_terms[go_id] = go_name
                if words[1] == "Reactome;":
                    split = line.split(";")
                    reac_id = split[1].replace(" ", "")
                    pathway = split[2][1:-1]
                    pathways[reac_id] = pathway
                # TODO filter out isoforms
                if words[1] == "RefSeq;":
                    if words[-1].count('[') > 0:
                        u_ac_iso = words[-1][1:-1]
                        refs = words[2:-1]
                    else:
                        u_ac_iso = uniprot_ac.split('-')[0]
                        refs = words[2:]
                    refs = [x[:-1] for x in refs]
                    if u_ac_iso not in refseqs:
                        refseqs[u_ac_iso] = []
                    refseqs[u_ac_iso] += refs
    if uniprot_ac in refseqs:
        refseqs = ",".join(refseqs[uniprot_ac])
    else:
        refseqs = ''
    if return_id:
        return (wildtype_sequence, refseqs, go_terms, pathways, u_id)
    return (wildtype_sequence, refseqs, go_terms, pathways)
