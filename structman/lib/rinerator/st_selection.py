import os.path

import Bio.PDB
from Bio.Data import PDBData

from . import sequence

# Copyright 2010 Max-Planck-Institut Informatik
#
#    This file is part of RINerator
#
#    RINerator is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RINerator is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RINerator.  If not, see <http://www.gnu.org/licenses/>.

# Table of 3 letter codes for post translation modifications
aa_trans_dic = PDBData.protein_letters_3to1.copy()
aa_trans_dic.update({'AP1': '1', 'AP2': '2', 'AP3': '3', 'AP4': '4', 'P1I': '5', 'EP2': '6'})

mc_atoms_dic = {'N': '', 'CA': '', 'C': '', 'O': '', 'H': '', 'HA': '', 'OXT': '', 'HA2': '', 'HA3': ''}


def set_str_sel(sel_id, pdb_file_bname, pdb_path, selection_lst):
    """
    Set the structure selection
    selection_lst = [[comp_name, comp_type, res_ranges],...]
    """
    pdb_file = os.path.join(pdb_path, pdb_file_bname)
    stsel_obj = StSel(sel_id, pdb_file)
    stsel_obj.get_pdb()
    for component_desc in selection_lst:
        [comp_name, comp_type, res_ranges] = component_desc
        stsel_obj.set_component(comp_name, comp_type, res_ranges)
    stsel_obj.set_res_id_lst()

    stsel_obj.set_subid_dic()
    print("INFO: Read PDB %s for selection %s" % (pdb_file_bname, sel_id))
    return stsel_obj


class StSel:
    """
    Structure Selection
    """

    def __init__(self, sel_id, pdb_file):
        self.sel_id = sel_id
        self.pdb_file = pdb_file
        self.st_obj = None  # Bio.PDB structure object
        self.seq_obj = None  # Sequence object
        self.protein_names_lst = []  # list of names of selected proteins
        self.protein_sel_dic = {}  # dict of selected proteins: dic[comp_name] = mcomp_dic
        self.ligand_sel_lst = []  # list of mcomp_dic of selected ligands
        # mcomp_dic['comp_name'] = comp_name ,label
        # mcomp_dic['comp_type'] = comp_type , protein, ligand, metal...
        # mcomp_dic['res_ranges'] = res_ranges , res_ranges = [range1,range2,...], range = [model, chain, [st_res_id_st, end_res_id]]
        # mcomp_dic['comp_res_id_lst'] = [res_id, res_id,...]
        self.res_id_dic = {}  # dict of res_id dic[res_id] = [comp_name,comp_type]
        self.res_id_lst = []  # sorted lsit of res_id, according to order of components in protein_names_lst, and then in ligand_names_lst
        self.res_subid_dic = {}  # dict of subset of residue identifiers, to interpret whatif data, ignores model and hetero flag
        #                      dic[(chain_id,res_seq,icode,resname)] = res_id
        self.atom_label_dic = {}  # dict of atom_label dic[atom_label] = [comp_name, comp_type] atom_label = (atom_id,atom_altloc)
        self.atom_subid_dic = {}  # dict of subset of atom identifiers, to interpret whatif data, ignores model hetero flag alternative_loc
        #                      dic[(chain_id,res_seq,icode,resname,atom_name)] = atom_label

    def get_pdb(self):
        parser_obj = Bio.PDB.PDBParser()
        self.st_obj = parser_obj.get_structure(self.sel_id, self.pdb_file)

    def set_component(self, comp_name, comp_type, res_ranges):
        mcomp_dic = self.set_mcomp(comp_name, comp_type)
        mcomp_dic = self.set_res(mcomp_dic, res_ranges)
        if comp_type == 'protein':
            self.protein_names_lst.append(comp_name)
            self.protein_sel_dic[comp_name] = mcomp_dic
        else:
            self.ligand_sel_lst.append(mcomp_dic)

    def set_res_id_lst(self):
        """
        set a ordered list of all residues in selection
        first the protein components ordered in protein_names_lst
        then the ligand components ordered in ligand_names_lst
        """
        self.res_id_lst = []
        comp_names_lst = []
        for protein_name in self.protein_names_lst:
            mcomp_dic = self.protein_sel_dic[protein_name]
            comp_res_id_lst = self.get_comp_res_id_lst(mcomp_dic)
            for res_id in comp_res_id_lst:
                self.res_id_lst.append(res_id)
        for mcomp_dic in self.ligand_sel_lst:
            comp_res_id_lst = self.get_comp_res_id_lst(mcomp_dic)
            for res_id in comp_res_id_lst:
                self.res_id_lst.append(res_id)

    def set_subid_dic(self):
        """
        set a dict of subset of residue identifiers and a dict of subset of atom identifiers
        only needed to interpret data from whatif
        """
        self.res_subid_dic = {}
        self.atom_subid_dic = {}
        for res_id in self.res_id_lst:
            chain_id = self.get_chain_from_res_id(res_id)
            res_seq = self.get_res_seq_from_res_id(res_id)
            icode = self.get_icode_from_res_id(res_id)
            resname = self.get_resname_from_res_id(res_id)
            dic_key = self.get_res_id_key(chain_id, res_seq, icode, resname)
            # res_subid_dic
            if dic_key in self.res_subid_dic:
                print("WARNING: Cannot identify residues. Repeated residue identifier: chain %s res_seq %s icode %s resname %s" % (chain_id, res_seq, icode, resname))
            else:
                self.res_subid_dic[dic_key] = res_id
            # atom_subid_dic
            res_obj = self.res_obj_from_res_id(res_id)
            atom_lst = res_obj.get_list()
            for atom_obj in atom_lst:
                atom_name = atom_obj.get_name()
                altloc = atom_obj.get_altloc()
                dic_key = self.get_atom_subid_key(chain_id, res_seq, icode, resname, atom_name, altloc)
                if dic_key in self.atom_subid_dic:
                    print("WARNING: Cannot identify atoms. Repeated atom identifier: chain %s res_seq %s icode %s resname %s atom_name %s altloc %s" % (chain_id, res_seq, icode, resname, atom_name, altloc))
                else:
                    atom_label = self.get_atom_label(atom_obj)
                    self.atom_subid_dic[dic_key] = atom_label

    def set_mcomp(self, comp_name, comp_type):
        """
        Initialise molecular component dictionary
        Can be a protein or ligand...
        """
        mcomp_dic = {}
        mcomp_dic['comp_name'] = comp_name  # label
        mcomp_dic['comp_type'] = comp_type  # protein, ligand, metal...
        mcomp_dic['res_ranges'] = {}  # res_ranges , res_ranges = [range1,range2,...], range = [model, chain, [st_res_id_st, end_res_id]]
        mcomp_dic['comp_res_id_lst'] = {}  # [res_id, res_id,...]
        return mcomp_dic

    def set_mcomp_res_ranges(self, mcomp_dic, res_ranges):
        mcomp_dic['res_ranges'] = res_ranges
        return mcomp_dic

    def set_mcomp_comp_res_id_lst(self, mcomp_dic, comp_res_id_lst):
        mcomp_dic['comp_res_id_lst'] = comp_res_id_lst
        return mcomp_dic

    def get_comp_name(self, mcomp_dic):
        return mcomp_dic['comp_name']

    def get_comp_type(self, mcomp_dic):
        return mcomp_dic['comp_type']

    def get_comp_res_id_lst(self, mcomp_dic):
        return mcomp_dic['comp_res_id_lst']

    def set_res(self, mcomp_dic, res_ranges):
        comp_res_id_lst = []
        # use the model of first segment, models of other segments ignored
        model_id = self.get_model_id_from_res_ranges(res_ranges)
        for seg_id in res_ranges:
            chain_id = self.get_chain_id_from_seg_id(seg_id)
            # If no chain is specified assume all possible chains
            if chain_id == " ":
                chain_lst = self.st_obj[model_id].get_list()
            else:
                chain_lst = [self.st_obj[model_id][chain_id]]
            for chain_obj in chain_lst:
                st_flag = False
                res_lst = chain_obj.get_list()
                for res_obj in res_lst:
                    if self.check_aa(res_obj) == False and self.get_comp_type(mcomp_dic) == 'protein':
                        continue
                    res_id = res_obj.get_full_id()
                    if not st_flag:
                        if self.check_st_res(res_id, seg_id):
                            comp_res_id_lst.append(res_id)
                            self.store_id_dic(res_id, mcomp_dic)
                            if self.get_end_res_seq_from_seg_id(seg_id) is None:
                                break
                            st_flag = True
                    else:
                        comp_res_id_lst.append(res_id)
                        self.store_id_dic(res_id, mcomp_dic)
                        if self.check_end_res(res_id, seg_id):
                            break
        mcomp_dic = self.set_mcomp_res_ranges(mcomp_dic, res_ranges)
        mcomp_dic = self.set_mcomp_comp_res_id_lst(mcomp_dic, comp_res_id_lst)
        return mcomp_dic

    def check_aa(self, res_obj):
        hetf = res_obj.get_id()[0]
        if hetf == ' ' or hetf[:2] == 'H_':
            res_name = res_obj.get_resname()
            if res_name in aa_trans_dic:
                return True
        print("INFO: Non amino acid: %s" % str(res_obj.get_full_id()))
        return False

    def check_st_res(self, res_id, seg_id):
        hetf = self.get_hetf_from_res_id(res_id)
        res_seq = self.get_res_seq_from_res_id(res_id)
        icode = self.get_icode_from_res_id(res_id)
        st_hetf = self.get_st_hetf_from_seg_id(seg_id)
        st_res_seq = self.get_st_res_seq_from_seg_id(seg_id)
        st_icode = self.get_st_icode_from_seg_id(seg_id)
        if st_res_seq == '_':
            return True
        if res_seq == st_res_seq:
            if icode == st_icode:
                return True
        return False

    def check_end_res(self, res_id, seg_id):
        hetf = self.get_hetf_from_res_id(res_id)
        res_seq = self.get_res_seq_from_res_id(res_id)
        icode = self.get_icode_from_res_id(res_id)
        end_res_seq = self.get_end_res_seq_from_seg_id(seg_id)
        end_icode = self.get_end_icode_from_seg_id(seg_id)
        end_hetf = self.get_end_hetf_from_seg_id(seg_id)
        if end_res_seq == '_':
            return False
        if res_seq == end_res_seq:
            if icode == end_icode:
                return True
        return False

    def store_id_dic(self, res_id, mcomp_dic):
        comp_name = self.get_comp_name(mcomp_dic)
        comp_type = self.get_comp_type(mcomp_dic)
        self.store_res_id_dic(res_id, comp_name, comp_type)
        self.store_atom_label_dic(res_id, comp_name, comp_type)

    def store_res_id_dic(self, res_id, comp_name, comp_type):
        self.res_id_dic[res_id] = [comp_name, comp_type]

    def store_atom_label_dic(self, res_id, comp_name, comp_type):
        res_obj = self.res_obj_from_res_id(res_id)
        atom_lst = res_obj.get_list()
        for atom_obj in atom_lst:
            atom_label = self.get_atom_label(atom_obj)
            self.atom_label_dic[atom_label] = [comp_name, comp_type]

    def res_obj_from_res_id(self, res_id):
        res_obj = self.st_obj[res_id[1]][res_id[2]][res_id[3]]
        return res_obj

    def get_model_id_from_res_ranges(self, res_ranges):
        """
        get model id from ranges
        use the model of first segment, models of other segments ignored
        """
        model_id = res_ranges[0][0]
        return model_id

    def get_res_id_key(self, chain_id, res_seq, icode, resname):
        dic_key = (chain_id, res_seq, icode, resname)
        return dic_key

    def get_atom_subid_key(self, chain_id, res_seq, icode, resname, atom_name, altloc):
        subid = (chain_id, res_seq, icode, resname, atom_name, altloc)
        return subid

    def get_atom_label_from_subid(self, chain_id, res_seq, icode, resname, atom_name, altloc):
        atom_key = self.get_atom_subid_key(chain_id, res_seq, icode, resname, atom_name, altloc)
        if atom_key in self.atom_subid_dic:
            atom_label = self.atom_subid_dic[atom_key]
        else:
            atom_label = None
        return atom_label

    def get_chain_id_from_seg_id(self, seg_id):
        chain_id = seg_id[1]
        return chain_id

    def get_st_hetf_from_seg_id(self, seg_id):
        return seg_id[2][0]

    def get_st_res_seq_from_seg_id(self, seg_id):
        return seg_id[2][1]

    def get_st_icode_from_seg_id(self, seg_id):
        return seg_id[2][2]

    def get_end_hetf_from_seg_id(self, seg_id):
        return seg_id[3][0]

    def get_end_res_seq_from_seg_id(self, seg_id):
        return seg_id[3][1]

    def get_end_icode_from_seg_id(self, seg_id):
        return seg_id[3][2]

    def get_chain_from_res_id(self, res_id):
        return res_id[2][0]

    def get_hetf_from_res_id(self, res_id):
        return res_id[3][0]

    def get_res_seq_from_res_id(self, res_id):
        return res_id[3][1]

    def get_icode_from_res_id(self, res_id):
        return res_id[3][2]

    def get_resname_from_res_id(self, res_id):
        res_obj = self.st_obj[res_id[1]][res_id[2]][res_id[3]]
        resname = res_obj.get_resname()
        return resname

    def get_atom_label(self, atom_obj):
        atom_id = atom_obj.get_full_id()
        atom_altloc = atom_obj.get_altloc()
        return (atom_id, atom_altloc)

    def get_atom_id_from_label(self, atom_label):
        return atom_label[0]

    def get_atom_name_from_label(self, atom_label):
        return atom_label[0][4][0]

    def check_mc_sc_ligand_from_label(self, atom_label):
        """
        return:
        mc
        sc
        ligand
        """
        comp_type = self.get_comp_type_from_label(atom_label)
        if comp_type is None:
            return None
        if comp_type != 'protein':
            return 'ligand'
        atom_name = self.get_atom_name_from_label(atom_label)
        if atom_name in mc_atoms_dic:
            return 'mc'
        return 'sc'

    def get_comp_type_from_label(self, atom_label):
        if atom_label in self.atom_label_dic:
            return self.atom_label_dic[atom_label][1]
        else:
            return None

    def atom_label_to_res_id(self, atom_label):
        atom_id = self.get_atom_id_from_label(atom_label)
        res_id = self.get_res_id_from_atom_id(atom_id)
        return res_id

    def get_res_id_from_atom_id(self, atom_id):
        return atom_id[:-1]

    def get_seq(self):
        char_lst = []
        for res_id in self.res_id_lst:
            res_obj = self.res_obj_from_res_id(res_id)
            resname = res_obj.get_resname()
            res_code = aa_trans_dic[resname]
            char_lst.append(res_code)
        sep = ''
        res_seq_str = sep.join(char_lst)
        seq_id = self.sel_id
        seq_obj = sequence.Sequence()
        seq_obj.set_seq(seq_id, res_seq_str)
        self.seq_obj = seq_obj

    def get_res_label_from_seq_idx(self, seq_idx):
        res_id = self.res_id_lst[seq_idx]
        chain_id = self.get_chain_from_res_id(res_id)
        res_seq = self.get_res_seq_from_res_id(res_id)
        icode = self.get_icode_from_res_id(res_id)
        if chain_id == ' ':
            chain_id = '_'
        if icode == ' ':
            icode = '_'
        res_label = "%s:%s:%s" % (chain_id, str(res_seq), icode)
        return res_label
