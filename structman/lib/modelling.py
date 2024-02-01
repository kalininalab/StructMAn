
import contextlib
import os
import sys
import psutil
import ray
import time
import traceback

try:
    import modeller
    from modeller.automodel import AutoModel
except:
    pass

from structman.lib import globalAlignment, pdbParser, serializedPipeline, templateFiltering
from structman.lib.sdsc import structure as structure_package
from structman.lib import model as model_class
from structman.lib.sdsc.consts import residues as residue_consts

modeller_chains_order = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz'


def create_chain_map(chains):
    mod_chains = []
    chain_map = {}
    for pos, chain in enumerate(chains):
        if pos >= len(modeller_chains_order):
            return
        mod_chain = modeller_chains_order[pos]
        mod_chains.append(mod_chain)
        chain_map[chain] = mod_chain
    return chain_map, mod_chains


class CompleteModel(AutoModel):
    def special_restraints(self, aln):
        for homochains in self.homomer_sets:
            if len(homochains) < 2:
                continue
            selections = []
            for chain in homochains:
                selections.append(modeller.selection(self.chains[chain]).only_atom_types('CA'))
            for pos1, selection1 in enumerate(selections):
                for pos2, selection2 in enumerate(selections):
                    if pos2 <= pos1:
                        continue
                    if len(selection1) != len(selection2):
                        continue
                    self.restraints.symmetry.append(modeller.symmetry(selection1, selection2, 1.0))


def add_incomplete_residue_gaps(protein_aligned_sequence, template_aligned_sequence, template_sequence):
    ta_pos = 0
    for pos, aa in enumerate(template_sequence):
        if ta_pos < len(template_aligned_sequence):
            while template_aligned_sequence[ta_pos] == '-':
                if (ta_pos + 1) >= len(template_aligned_sequence) :
                    break
                ta_pos += 1
        # this is for incompletely resolved residues
        if aa == '-':
            protein_aligned_sequence = '%s-%s' % (protein_aligned_sequence[:ta_pos], protein_aligned_sequence[ta_pos:])
            template_aligned_sequence = '%s-%s' % (template_aligned_sequence[:ta_pos], template_aligned_sequence[(ta_pos):])
            #ta_pos += 1
        # this is for non-standard residues
        elif aa == '.':
            template_aligned_sequence = '%s.%s' % (template_aligned_sequence[:ta_pos], template_aligned_sequence[(ta_pos + 1):])
        elif aa == '?':
            template_aligned_sequence = '%s?%s' % (template_aligned_sequence[:ta_pos], template_aligned_sequence[(ta_pos + 1):])
        ta_pos += 1
    return protein_aligned_sequence, template_aligned_sequence


def createFullAlignment(config, pdb_id, prot_id, first_residue, last_residue, chains, sequences, target_chain, target_chains, ligand_counts, protein_aligned_sequence, template_aligned_sequence, removed_ligands, mutation_positions = None):

    if config.verbosity >= 5:
        print(f'Call of createFullAlignment: PDB ID: {pdb_id}:{target_chain} ({target_chains}), Prot ID: {prot_id}\nLigand Counts: {ligand_counts}\nRemoved ligands: {removed_ligands}')
        print(f'Protein aligned sequence before truncate:\n{protein_aligned_sequence}')

    if mutation_positions is not None:
        protein_aligned_sequence, template_aligned_sequence, remapped_positions = globalAlignment.truncateSequences(protein_aligned_sequence, template_aligned_sequence, flank_buffer = 5, remap_positions = mutation_positions)
    else:
        protein_aligned_sequence, template_aligned_sequence = globalAlignment.truncateSequences(protein_aligned_sequence, template_aligned_sequence, flank_buffer = 5)
        remapped_positions = None

    pdb_id = pdb_id.lower()

    if config.verbosity >= 5:
        print(f'Protein aligned sequence after truncate:\n{protein_aligned_sequence}')
        print('Template aligned sequence after truncate:\n', template_aligned_sequence)
        print('Template sequence after truncate:\n', sequences[target_chain])

    protein_aligned_sequence, template_aligned_sequence = add_incomplete_residue_gaps(protein_aligned_sequence, template_aligned_sequence, sequences[target_chain])

    if config.verbosity >= 5:
        print('Template aligned sequence after add_incomplete_residue_gaps:\n', template_aligned_sequence)

    truncated_prot_seq_len = len(protein_aligned_sequence) - protein_aligned_sequence.count('-')

    first_chain = chains[0]
    last_chain = chains[-1]

    temp_seq_lines = []
    prot_seq_lines = []
    total_prot_len = 0
    for chain in chains:
        lc = ligand_counts[chain]
        if chain in removed_ligands:
            lc -= removed_ligands[chain]
        if chain in target_chains:
            if template_aligned_sequence[-1] == '?':
                template_aligned_sequence = template_aligned_sequence[:-1]
                temp_seq_lines.append(template_aligned_sequence + ('.' * lc) + ('.'))
                prot_seq_lines.append(protein_aligned_sequence[:-1] + ('.' * lc) + ('.'))
            else:
                temp_seq_lines.append(template_aligned_sequence + ('.' * lc))
                prot_seq_lines.append(protein_aligned_sequence + ('.' * lc))
            total_prot_len += truncated_prot_seq_len + lc
        else:
            if len(sequences[chain]) == 0:
                continue
            if sequences[chain][-1] == '?':
                temp_seq = sequences[chain][:-1] + ('.' * lc) + ('.')
                temp_seq_lines.append(temp_seq)
                prot_seq = sequences[chain][:-1] + ('.' * lc) + ('.')
            else:
                temp_seq = sequences[chain] + ('.' * lc)
                temp_seq_lines.append(temp_seq)
                prot_seq = temp_seq
            prot_seq_lines.append(prot_seq)
            total_prot_len += len(prot_seq)

    template_header = '>P1;%s\nstructureX:%s:%s:%s:%s:%s::::\n' % (pdb_id, pdb_id, first_residue, first_chain, last_residue, last_chain)

    #template_header = '>P1;%s\nstructureX:%s:@:@:@:@::::\n' % (pdb_id,pdb_id)

    template_seq = '%s*\n' % '/\n'.join(temp_seq_lines)

    if config.verbosity >= 5:
        print('Template seq as is goes into the alignment pir:\n', template_seq)

    protein_header = '>P1;%s\nsequence:%s:1:%s:%s:%s::::\n' % (prot_id, prot_id, first_chain, str(total_prot_len), last_chain)

    protein_seq = '%s*\n' % '/\n'.join(prot_seq_lines)

    modeller_alignment = '%s%s\n%s%s' % (template_header, template_seq, protein_header, protein_seq)

    return modeller_alignment, remapped_positions


@ray.remote(num_cpus = 1, max_calls = 1)
def model_remote_wrapper(config, compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, skip_analysis=False, force_modelling=False, mutation_positions = None, label_add = '', target_path = None):
    return model(config, compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, skip_analysis=skip_analysis, force_modelling=force_modelling, mutation_positions = mutation_positions, label_add = label_add, target_path = target_path)

def model(config, compl_obj, structures, alignment_tuple, seq_id, cov, pdb_id, tchain, prot_id, skip_analysis=False, force_modelling=False, mutation_positions = None, label_add = '', target_path = None):
    def model_core():
        env = modeller.Environ()
        env.io.atom_files_directory = [process_folder]
        env.io.hetatm = True

        model_obj = CompleteModel(env, alnfile=tmp_aln_file, knowns=pdb_id.lower(), sequence=truncated_prot_id)

        homomers = compl_obj.getHomomersSets()
        mod_homomers = []
        for hom_chains in homomers:
            mod_hom_chains = []
            for chain in hom_chains:
                if chain not in modeller_chain_id_map:
                    continue
                mod_hom_chains.append(modeller_chain_id_map[chain])
            if len(mod_hom_chains) == 0:
                continue
            mod_homomers.append(mod_hom_chains)
        model_obj.homomer_sets = mod_homomers
        try:
            model_obj.make()
        except:
            if config.verbosity <= 2:
                sys.stdout = sys.__stdout__
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            errortext = '\n'.join([str(e), str(f), str(g)]) + '\n\n'
            config.errorlog.add_error('Modelling failed: %s %s %s\n%s' % (model_id, str(homomers), str(mod_homomers), errortext))

    t0 = time.time()

    serializedPipeline.ray_hack()

    if config.verbosity >= 3:
        print('Begin modelling of:', prot_id, pdb_id, tchain, label_add)

    # Step zero: create process-specific temporary cwd folder
    truncated_prot_id = prot_id[:50]
    model_id = f'{truncated_prot_id}_{pdb_id}_{tchain}{label_add}'

    model_id = model_id.replace('.', '_').replace('/', '_')

    if not os.path.exists('%s/model_base' % (config.tmp_folder)):
        try:
            os.mkdir('%s/model_base' % (config.tmp_folder))
        except:
            pass

    process_top_folder = '%s/model_base/%s' % (config.tmp_folder, model_id[:6])
    if not os.path.exists(process_top_folder):
        try:
            os.mkdir(process_top_folder)
        except:
            pass

    id_parts = model_id.split('_')
    if len(id_parts[0]) == 2:
        if len(id_parts[2]) == 1:
            wt_prot_id = '_'.join(id_parts[:3])
        else:
            wt_prot_id = '_'.join(id_parts[:2])
    else:
        wt_prot_id = id_parts[0]

    process_prot_folder = '%s/%s' % (process_top_folder, wt_prot_id)
    if not os.path.exists(process_prot_folder):
        try:
            os.mkdir(process_prot_folder)
        except:
            pass

    process_folder = '%s/%s' % (process_prot_folder, model_id)

    if not os.path.exists(process_folder):
        try:
            os.mkdir(process_folder)
        except:
            pass

    os.chdir(process_folder)
    config.temp_folder = process_folder

    model_path = None

    for fn in os.listdir(process_folder):
        full_path = '%s/%s' % (process_folder, fn)
        if fn[-4:] == '.pdb':
            if len(fn) <= 9:
                continue
            model_path = full_path

    # Step one: create full alignment: get sequences of all chains
    #template_page = compl_obj.getPage(config)

    if compl_obj.interaction_partners is None:
        template_page = compl_obj.getPage(config, self_update=True)
    else:
        template_page = compl_obj.getPage(config)

    (coordinate_map, siss_coord_map, res_contig_map, ligands, metals, ions, box_map, chain_type_map, chainlist,
     b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map) = templateFiltering.parsePDB(template_page)

    chains = compl_obj.chainlist
    reduced_chains = []
    for chain in chains:
        if chain == tchain:
            reduced_chains.append(tchain)
            continue

        if chain not in box_map:
            return 'Chain not in box_map %s %s %s %s' % (prot_id, pdb_id, tchain, chain)

        neighbors, center_dist = templateFiltering.box_check(box_map[tchain], box_map[chain], distance_threshold=config.short_distance_threshold)
        if neighbors:
            reduced_chains.append(chain)

    chains = reduced_chains
    chain_map_out = create_chain_map(chains)

    if chain_map_out is None:
        return 'Too many chains for %s %s %s %s' % (prot_id, pdb_id, tchain, str(chains))

    modeller_chain_id_map, modeller_chain_list = chain_map_out

    if model_path is None or force_modelling:

        if model_path is not None and os.path.exists(model_path):
            os.remove(model_path)

        # adding structure objects for not directly mapped chains
        for chain in chains:
            if not (pdb_id, chain) in structures:
                struct = structure_package.Structure(pdb_id, chain)
                struct.parse_page(template_page, config)
                structures[(pdb_id, chain)] = struct

        first_residue = structures[(pdb_id, chains[0])].first_residue
        last_residue = structures[(pdb_id, chains[-1])].last_residue
        sequences = {}
        ligand_counts = {}
        for chain in chains:
            sequences[chain] = structures[(pdb_id, chain)].getSequence(config, complex_obj=compl_obj, for_modeller=True)
            ligand_counts[chain] = compl_obj.countLigands(chain)

        if config.verbosity >=5:
            print(f'Complex IAP {compl_obj.interaction_partners}')

        # Unfortunately, we have to unmark homomeric chains, when their sequence does not match
        for base_chain in compl_obj.homomers:
            if base_chain not in modeller_chain_id_map:
                continue
            del_positions = []
            for pos, homo_chain in enumerate(compl_obj.homomers[base_chain]):
                if homo_chain not in modeller_chain_id_map:
                    del_positions.append(pos)
                elif len(sequences[homo_chain]) != len(sequences[base_chain]):
                    del_positions.append(pos)
                elif sequences[homo_chain] != sequences[base_chain]:
                    del_positions.append(pos)

            for pos in reversed(del_positions):
                if config.verbosity >= 2:
                    print(f'Had to remove homomer information for {pdb_id}: {base_chain} - {compl_obj.homomers[base_chain][pos]} due to non-identic sequences')
                del compl_obj.homomers[base_chain][pos]

        target_chains = compl_obj.homomers[tchain]

        if alignment_tuple is None:
            config.errorlog.add_error('Alignment not found: %s %s:%s' % (prot_id, pdb_id, tchain))

        protein_aligned_sequence, template_aligned_sequence = alignment_tuple

        # Step two: prepare environment: path to structure, alignment file, look for homooligo relationship, set symmetry restraints

        tmp_template_file = '%s/%s.pdb' % (process_folder, pdb_id.lower())

        if config.verbosity >= 3:
            print('Calling og relocate_hetatm with:', chains, residue_consts.NOT_IN_MODELLER)

        template_page, removed_ligands = pdbParser.relocate_hetatm(template_page, filter_chains=set(chains), filter_het=residue_consts.NOT_IN_MODELLER)

        if config.verbosity >= 6:
            print('Resulting page:\n', template_page)

        f = open(tmp_template_file, 'w')
        f.write(template_page)
        f.close()

        modeller_alignment, remapped_positions = createFullAlignment(config, pdb_id, truncated_prot_id, first_residue, last_residue, chains, sequences, tchain, target_chains, ligand_counts, protein_aligned_sequence, template_aligned_sequence, removed_ligands, mutation_positions = mutation_positions)

        if remapped_positions is not None:
            highlight_mutant_residue = {}
            for chain in target_chains:
                highlight_mutant_residue[chain] = remapped_positions
        else:
            highlight_mutant_residue = None

        if config.verbosity >= 3:
            print(f'{model_id} - Target chains: {target_chains}')
            print(f'{model_id} - Remapped highlight mutant residues:\n{remapped_positions}')

        tmp_aln_file = '%s/%s_%s.ali' % (process_folder, pdb_id.lower(), truncated_prot_id)

        f = open(tmp_aln_file, 'w')
        f.write(modeller_alignment)
        f.close()

        # Step three: model
        if config.verbosity <= 5:
            with contextlib.redirect_stdout(None), contextlib.redirect_stderr(None):
                model_core()
        else:
            model_core()

        if config.verbosity <= 2:
            try:
                os.remove(tmp_template_file)

            except:
                print(f'Could not remove: {tmp_template_file}')
            try:
                os.remove(tmp_aln_file)
            except:
                print(f'Could not remove: {tmp_aln_file}')
        else:
            print(f'Path to alignment tmp file: {tmp_aln_file}')

        model_path = None

        for fn in os.listdir(process_folder):
            full_path = '%s/%s' % (process_folder, fn)
            if full_path == tmp_template_file:
                continue
            if fn[-4:] == '.pdb':
                model_path = full_path

        if model_path is None:
            return 'No Model got produced for %s %s %s' % (prot_id, pdb_id, tchain)

    else:
        if config.verbosity >= 3:
            print('Found model:', model_path)
        
        if mutation_positions is not None:
            if alignment_tuple is None:
                config.errorlog.add_error('Alignment not found: %s %s:%s' % (prot_id, pdb_id, tchain))

            protein_aligned_sequence, template_aligned_sequence = alignment_tuple
            protein_aligned_sequence, template_aligned_sequence, remapped_positions = globalAlignment.truncateSequences(protein_aligned_sequence, template_aligned_sequence, flank_buffer = 5, remap_positions = mutation_positions)

            highlight_mutant_residue = {}
            target_chains = compl_obj.homomers[tchain]
            for chain in target_chains:
                highlight_mutant_residue[chain] = remapped_positions
        else:
            highlight_mutant_residue = None

    structman_model_obj = model_class.Model(model_id=model_id, path=model_path, tmp_folder=process_folder, template_structure=(pdb_id, tchain),
                                            target_protein=prot_id, chain_id_map=modeller_chain_id_map, truncated_prot_id = truncated_prot_id,
                                            template_resolution=compl_obj.resolution, sequence_identity=seq_id, coverage=cov, template_contig_map = res_contig_map, label_add = label_add)

    structman_model_obj.analyze(config, highlight_mutant_residue =  highlight_mutant_residue, target_path = f'{target_path}/{model_id}.pdb')

    if config.verbosity <= 2:
        for fn in os.listdir(process_folder):
            if fn[0] == '.':
                continue
            full_path = '%s/%s' % (process_folder, fn)
            if not (fn[-4:] == '.pdb' or fn[-3:] == '.gz'):
                if os.path.exists(full_path):
                    os.remove(full_path)

    t1 = time.time()

    if config.verbosity >= 3:
        print('Finished modelling of:', prot_id, pdb_id, tchain, ',computation time:', (t1 - t0))
    return (structman_model_obj)
