import os

import structman.lib.database.database as database

def create_ppi_network(session, config, outfile):
    proteins = database.proteinsFromDb(session, config, with_residues = True, with_alignments = True)

    # Step 1: build a residue -> position and a position -> residue map
    r_m_map = {}

    u_acs = proteins.get_protein_ids()

    if config.verbosity >= 2:
        print('Amount of proteins:', len(u_acs))

    n_struct_anno = 0
    n_positions = 0
    n_sub_infos = 0

    for u_ac in u_acs:

        positions = proteins.get_position_ids(u_ac)

        annotation_list = proteins.get_protein_annotation_list(u_ac)

        n_struct_anno += len(annotation_list)
        n_positions += len(positions)

        for (pdb_id, chain) in annotation_list:

            sub_infos = proteins.get_sub_infos(u_ac, pdb_id, chain)
            n_sub_infos += len(sub_infos)
            for pos in positions:
                if pos not in sub_infos:
                    continue
                sub_info = sub_infos[pos]
                res_nr = sub_info[0]

                if not proteins.contains_residue(pdb_id, chain, res_nr):
                    continue
                r_id = proteins.get_residue_db_id(pdb_id, chain, res_nr)
                res_aa = proteins.get_residue_aa(pdb_id, chain, res_nr)
                if r_id not in r_m_map:
                    r_m_map[r_id] = []
                r_m_map[r_id].append((u_ac, pos))

    if config.verbosity >= 2:
        print('Size of r_m_map:', len(r_m_map))
        print('Amount of structure annotations:', n_struct_anno)
        print('Amount of positions:', n_positions)
        print('Amount of sub_infos:', n_sub_infos)

    # Step 2: build a residue,chain A - residue,chain B datastructure
    struct_res_res_dict = {}
    for (pdb_id, chain) in proteins.structures:
        struct_res_res_dict[(pdb_id, chain)] = {}
        for res_nr in proteins.structures[(pdb_id, chain)].residues:
            scd = proteins.get_residue_scd(pdb_id, chain, res_nr)
            r_id = proteins.get_residue_db_id(pdb_id, chain, res_nr)

            if scd is None:
                continue

            scds = scd.split(",")

            for sd in scds:
                sdi = sd.split(":")
                if len(sdi) > 1:
                    ichain = sdi[0][0]
                    res_nr_2 = sdi[0].split('.')[1]
                    mc_d = float(sdi[1])
                    if mc_d < 5.0:
                        if not proteins.contains_structure(pdb_id, ichain):
                            continue
                        if not (r_id, res_nr) in struct_res_res_dict[(pdb_id, chain)]:
                            struct_res_res_dict[(pdb_id, chain)][(r_id, res_nr)] = set([])
                        struct_res_res_dict[(pdb_id, chain)][(r_id, res_nr)].add((ichain, res_nr_2))

    if config.verbosity >= 2:
        print('Size of struct_res_res_dict:', len(struct_res_res_dict))

    if os.path.exists(outfile):
        os.remove(outfile)

    f = open(outfile, 'a')
    f.write('u_ac_1\taacbase_1\tpdb_id\tchain_1\tres_nr_1\tu_ac_2\taacbase_2\tchain_2\tres_nr_2\n')

    # Step 3 produce the table

    for pdb_id, chain_1 in struct_res_res_dict:
        for r_1, res_nr_1 in struct_res_res_dict[(pdb_id, chain_1)]:
            if r_1 not in r_m_map:
                continue
            for (chain_2, res_nr_2) in struct_res_res_dict[(pdb_id, chain_1)][(r_1, res_nr_1)]:

                r_2 = proteins.get_residue_db_id(pdb_id, chain_2, res_nr_2)

                if r_2 not in r_m_map:
                    continue

                for (u_ac_1, aacbase_1) in r_m_map[r_1]:
                    for (u_ac_2, aacbase_2) in r_m_map[r_2]:

                        f.write('\t'.join([u_ac_1, aacbase_1, pdb_id, chain_1, res_nr_1, u_ac_2, aacbase_2, pdb_id, chain_2, res_nr_2]) + '\n')
    f.close()