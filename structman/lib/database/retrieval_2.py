from structman.lib.database.database_core_functions import binningSelect
from structman.base_utils.base_utils import unpack
from structman.lib.sdsc import residue as residue_package

def getBackmaps(
        structure_db_ids: dict[int, tuple[str, str]],
        config,
        unpacked = True,
        locked = True, db_lock = None
        ) -> dict[str, dict[str, dict[str, residue_package.Residue_Map[int] | bytes]]]:
    columns = ['Structure', 'Protein', 'Backmap']
    table = 'Alignment'
    results = binningSelect(structure_db_ids, columns, table, config, locked=locked, db_lock=db_lock)

    backmaps: dict[str, dict[str, dict[str, residue_package.Residue_Map[int] | bytes]]] = {}

    for row in results:
        structure_id, chain = structure_db_ids[row[0]]
        if structure_id not in backmaps:
            backmaps[structure_id] = {}
        if chain not in backmaps[structure_id]:
            backmaps[structure_id][chain] = {}
        if unpacked:
            backmap = unpack(row[2])
        else:
            backmap = row[2]
        backmaps[structure_id][chain][row[1]] = backmap

    return backmaps