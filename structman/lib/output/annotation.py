import os
from structman.lib.database import database
from structman.lib.output import out_generator

# called by structman_main
def create_annotation_table(session, config, outfile):
    if os.path.exists(outfile):
        os.remove(outfile)
    proteins = database.proteinsFromDb(session, config, with_residues=True, filter_mutant_proteins=True, with_alignments=True)

    fat_output = out_generator.OutputGenerator()

    headers = ['Protein-ID', 'Input-ID', 'Position', 'PDB-ID', 'Chain', 'Residue-ID', 'Residue AA type',
               'Sequence identity', 'Coverage', 'Resolution',
               'Residue RIN simple classification',
               'Chain type interaction partners', 'Small molecules type interaction partners']

    fat_output.add_headers(headers)
    f = open(outfile, 'a')
    f.write(fat_output.get_header())

    prot_ids = proteins.get_protein_ids()
    if config.verbosity >= 2:
        print('Number of proteins:', len(prot_ids))

    n_annos = 0
    n_positions = 0

    for prot_id in prot_ids:
        annotation_list = proteins.get_protein_annotation_list(prot_id)
        input_id = proteins[prot_id].input_id
        positions = proteins.get_position_ids(prot_id)
        n_annos += len(annotation_list)
        n_positions += len(positions)
        for (pdb_id, chain) in annotation_list:
            sub_infos = proteins.get_sub_infos(prot_id, pdb_id, chain)
            resolution = proteins.complexes[pdb_id].resolution
            sequence_id = proteins[prot_id].structure_annotations[(pdb_id, chain)].sequence_identity
            coverage = proteins[prot_id].structure_annotations[(pdb_id, chain)].coverage

            for pos in positions:
                if sub_infos[pos] is None:
                    continue
                sub_info = sub_infos[pos]
                res_nr = sub_info[0]
                res_aa = sub_info[1]

                if not proteins.contains_residue(pdb_id, chain, res_nr):
                    continue
                fat_output.add_value('Protein-ID', prot_id)
                fat_output.add_value('Input-ID', input_id)
                fat_output.add_value('Position', pos)
                fat_output.add_value('PDB-ID', pdb_id)
                fat_output.add_value('Chain', chain)
                fat_output.add_value('Residue-ID', res_nr)
                fat_output.add_value('Residue AA type', res_aa)
                fat_output.add_value('Protein-ID', prot_id)

                fat_output.add_value('Sequence identity', sequence_id)
                fat_output.add_value('Coverage', coverage)
                fat_output.add_value('Resolution', resolution)

                rin_class, rin_simple_class = proteins.structures[pdb_id][chain].residues[res_nr].get_classification(config)

                fat_output.add_value('Residue RIN simple classification', rin_simple_class)

                interacting_chains, interacting_ligands = proteins.structures[pdb_id][chain].residues[res_nr].get_interaction_partners()

                fat_output.add_value('Chain type interaction partners', ','.join(interacting_chains))
                fat_output.add_value('Small molecules type interaction partners', ','.join(interacting_ligands))

                f.write(fat_output.pop_line())
    f.close()

    if config.verbosity >= 2:
        print('Number of annotated structures:', n_annos)
        print('Number of positions:', n_positions)


# structure of anno_anno_map: {(Mutation_id_1,Mutation_id_2):(distance,atom,atom2,template_id_1,template_id_2,chain1,chain2)} // Mutation_id_1 < Mutation_id_2
# structure of m_aac_map: {Mutation_id:(gene_id,AAC_base)}
# structure of gn_map: {Gene_id:Uniprot_Ac}ä
def annoAnnoNetwork(anno_anno_map, m_aac_map, gn_map, outfile, distance_threshold=None):
    lines = []
    for (m_id1, m_id2) in anno_anno_map:
        (min_d, atom, atom2, t_id1, t_id2, chain1, chain2) = anno_anno_map[(m_id1, m_id2)]
        if distance_threshold is not None:
            if min_d > distance_threshold:
                continue
        (g_id1, aac1) = m_aac_map[m_id1]
        (g_id2, aac2) = m_aac_map[m_id2]
        u_ac1 = gn_map[g_id1]
        u_ac2 = gn_map[g_id2]

        lines.append("%s:%s %s %s:%s" % (u_ac1, aac1, str(min_d), u_ac2, aac2))

    f = open(outfile, 'w')
    f.write('\n'.join(lines))
    f.close()
