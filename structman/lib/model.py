import shutil
import os
import sys
import traceback
from structman.lib import templateFiltering


class Model:
    __slots__ = ['path', 'template_structure', 'target_protein', 'model_id', 'truncated_prot_id', 'structural_analysis_dict', 'ligand_profiles',
                 'metal_profiles', 'ion_profiles', 'chain_chain_profiles', 'chain_type_map', 'chainlist', 'chain_id_map',
                 'template_resolution', 'sequence_identity', 'coverage', 'tmp_folder', 'template_contig_map', 'label_add']

    def __init__(self, path='', template_structure=None, target_protein=None, model_id=None, chain_id_map=None, truncated_prot_id=None,
                 template_resolution=None, sequence_identity=None, coverage=None, tmp_folder=None, template_contig_map=None, label_add = ''):
        self.path = path
        self.template_structure = template_structure  # tuple of template pdb id and target chain
        self.target_protein = target_protein
        self.model_id = model_id
        self.truncated_prot_id = truncated_prot_id
        self.chain_id_map = chain_id_map  # maps the template chain ids to the model chain ids
        self.template_resolution = template_resolution
        self.sequence_identity = sequence_identity
        self.coverage = coverage
        self.tmp_folder = tmp_folder
        self.template_contig_map = template_contig_map
        self.label_add = label_add

    def analyze(self, config, highlight_mutant_residue = None, target_path= None):
        config.n_of_chain_thresh = 1000  # Hinder structuralAnalysis to spawn processes, since this function is already called by a remote
        model_target_chain = self.chain_id_map[self.template_structure[1]]

        if config.verbosity >= 3:
            print('Start model self analysis:', self.model_id, self.path, highlight_mutant_residue)

        if self.path[-12:] != '_refined.pdb' or target_path is not None:
            try:
                self.refine_model(highlight_mutant_residue = highlight_mutant_residue, target_path = target_path)
            except:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                config.errorlog.add_error(f'refinde_model failed in Model.analyze with: Path: {self.path},\ntarget_path: {target_path},\nModel ID: {self.model_id},\ntarget chain: {model_target_chain},\nchain_id_map:{self.chain_id_map}\n{e}\n{f}\n{g}')
                try:
                    os.remove(self.path)
                except:
                    pass
                return

        try:
            (structural_analysis_dict, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, _, _) = templateFiltering.structuralAnalysis(self.model_id, config, model_path=self.path, target_dict=[model_target_chain], keep_rin_files=True)
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            config.errorlog.add_error(f'structuralAnalysis failed in Model.analyze with: Path: {self.path}, Model ID: {self.model_id}, target chain: {model_target_chain}\n{e}\n{f}\n{g}')
            return

        self.structural_analysis_dict = structural_analysis_dict
        self.ligand_profiles = ligand_profiles
        self.metal_profiles = metal_profiles
        self.ion_profiles = ion_profiles
        self.chain_chain_profiles = chain_chain_profiles
        self.chain_type_map = chain_type_map
        self.chainlist = chainlist

        for error_text in errorlist:
            config.errorlog.add_warning(error_text)

    def clear_tmp(self):
        if self.tmp_folder is None:
            return
        try:
            shutil.rmtree(self.tmp_folder)
        except:
            print('Could not remove tmp folder:', self.tmp_folder)
        return

    def convert_highlight_map(self, highlight_map):
        sav_highlights = set()
        insertion_highlights = set()
        deletion_highlights = set()

        sav_pos, insertion_pos, deletion_flanks = highlight_map

        if sav_pos is not None:
            for pos in sav_pos:
                sav_highlights.add(str(pos))

        if insertion_pos is not None:
            for ins in insertion_pos:
                for i_pos in ins:
                    insertion_highlights.add(str(i_pos))

        if deletion_flanks is not None:
            for lf, rf in deletion_flanks:
                deletion_highlights.add(str(lf - 1))
                deletion_highlights.add(str(lf))

        return sav_highlights, insertion_highlights, deletion_highlights

    def check_highlights(self, res_nr, highlight_map):
        sav_highlights, insertion_highlights, deletion_highlights = highlight_map
        if res_nr in sav_highlights:
            return '  0.00'
        if res_nr in insertion_highlights:
            return ' 50.00'
        if res_nr in deletion_highlights:
            return ' 75.00'
        return '100.00'

    # modeller increases the residue ID continuously even when a new chain begun, this can lead to residue IDs > 9999
    # modeller solves this by using letters, for example: A000 = 10000
    def refine_model(self, highlight_mutant_residue = None, target_path = None):

        #print('\n\n',self.path,'\n',target_path,'\n\n')

        #print(highlight_mutant_residue)
        #print(self.template_structure, self.chain_id_map)

        if highlight_mutant_residue is not None:
            converted_highlight_map = {}
            for chain in highlight_mutant_residue:
                converted_highlight_map[self.chain_id_map[chain]] = self.convert_highlight_map(highlight_mutant_residue[chain])
            #print(converted_highlight_map)

        f = open(self.path, 'r')
        lines = f.readlines()
        f.close()

        new_lines = []
        current_chain = None
        current_res = None
        residue_nr_offset = 0
        for line in lines:
            if len(line) >= 21:
                record_name = line[0:6].rstrip()
                if record_name == "ATOM" or record_name == 'HETATM':
                    chain_id = line[21]
                    res_nr = line[22:26].strip()  # without insertion code
                    # reset res_id counter for every new chain
                    if current_chain != chain_id:
                        current_new_res_nr = 0
                        current_chain = chain_id
                    if current_res != res_nr:
                        current_res = res_nr
                        current_new_res_nr += 1
                        digit_res_str = str(current_new_res_nr)
                        current_res_str = '%s%s' % (' ' * (4 - len(digit_res_str)), digit_res_str)

                    newline = f'{line[:22]}{current_res_str}{line[26:]}'

                    if highlight_mutant_residue is not None:
                        if chain_id in converted_highlight_map:
                            b_factor = self.check_highlights(digit_res_str, converted_highlight_map[chain_id])
                        else:
                            b_factor = '25.00'
                        newline = f'{newline[:60]}{b_factor}{newline[66:]}'

                    new_lines.append(newline)

                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)

        stem = self.path[:-4].replace('.', '_')
        if stem[-12:] != '_refined.pdb':
            refined_path = f'{stem}_refined.pdb'
        else:
            refined_path = f'{stem}.pdb'

        if target_path is not None:
            f = open(target_path, 'w')
            f.write(''.join(new_lines))
            f.close()

        f = open(refined_path, 'w')
        f.write(''.join(new_lines))
        f.close()

        if self.path != refined_path:
            os.remove(self.path)

        self.path = refined_path
        return
