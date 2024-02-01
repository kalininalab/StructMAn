import os
from structman.lib.database import database

def retrieve_annotated_structures(session, config, outfile):
    if os.path.exists(outfile):
        os.remove(outfile)

    prot_db_ids = database.getProtIdsFromSession(session, config, filter_mutant_proteins=True)

    rows = ['Protein', 'Structure']
    table = 'Alignment'
    results = database.binningSelect(prot_db_ids, rows, table, config)

    structure_db_ids = []

    for row in results:
        structure_db_ids.append(row[1])

    rows = ['Structure_Id', 'PDB']
    table = 'Structure'
    results = database.binningSelect(structure_db_ids, rows, table, config)

    pdb_ids = set()

    for row in results:
        pdb_ids.add(row[1])

    complex_map = database.getComplexMap(config, pdb_ids=pdb_ids)

    outlines = []

    for pdb_id in complex_map:
        chains_str = complex_map[pdb_id][2]
        pdb_id = pdb_id.replace('_AU','')
        for c_tuple in chains_str.split(','):
            chain,chain_type = c_tuple.split(':')
            if chain_type == 'Protein':
                outlines.append('%s:%s\n' % (pdb_id,chain))

    f = open(outfile, 'w')
    f.write(''.join(outlines))
    f.close()
