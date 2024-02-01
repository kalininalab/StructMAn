import copy

from . import nc_interactions, st_selection

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


def get_nci_res(nci_obj, stsel_obj, combine_int):
    nci_res = NCIResidues(nci_obj)
    nci_res.get_res_int_dic(stsel_obj, combine_int)
    return nci_res


class NCIResidues:
    def __init__(self, nci_obj):
        self.nci_obj = nci_obj
        self.res_int_dic = {}  # res_int_dic[res_int_key]= res_int_attr_obj
        #                       res_int_key = (res_id1,res_id2,int_type,int_subtype)
        self.res_int_lst = []  # res_int_lst=[res_int_key1,res_int_key2,...]

    def get_res_int_dic(self, stsel_obj, combine_int):
        for atom_int_key in self.nci_obj.atom_int_lst:
            atom_int_attr = self.nci_obj.atom_int_dic[atom_int_key]
            atom_label1 = atom_int_attr.id1
            # atom_label1 = self.nci_obj.atom_int_key_to_atom_label1(atom_int_key)
            res_id1 = stsel_obj.atom_label_to_res_id(atom_label1)
            atom_label2 = atom_int_attr.id2
            # atom_label2 = self.nci_obj.atom_int_key_to_atom_label2(atom_int_key)
            res_id2 = stsel_obj.atom_label_to_res_id(atom_label2)
            int_type = atom_int_attr.int_type
            # int_type = self.nci_obj.atom_int_key_to_int_type(atom_int_key)
            int_subtype = atom_int_attr.int_subtype
            res_int_key = self.get_res_int_key(res_id1, res_id2, int_type, int_subtype)
            res_int_rkey = self.get_res_int_key(res_id2, res_id1, int_type, int_subtype)
            if res_int_key in self.res_int_dic:
                self.update_res_int(res_int_key, atom_int_attr.int_score)
            elif res_int_rkey in self.res_int_dic:
                self.update_res_int(res_int_rkey, atom_int_attr.int_score)
            else:
                res_int_attr_obj = nc_interactions.IntAttr(res_id1, res_id2, int_type, int_subtype)
                self.res_int_dic[res_int_key] = res_int_attr_obj
                self.update_res_int_score(res_int_key, atom_int_attr.int_score)
                self.res_int_lst.append(res_int_key)
        if combine_int:
            self.get_combine_int()

    def get_res_int_key(self, res_id1, res_id2, int_type, int_subtype):
        return (res_id1, res_id2, int_type, int_subtype)

    def get_res_id1_from_key(self, res_int_key):
        return res_int_key[0]

    def get_res_id2_from_key(self, res_int_key):
        return res_int_key[1]

    def get_int_type_from_key(self, res_int_key):
        return res_int_key[2]

    def get_int_subtype_from_key(self, res_int_key):
        return res_int_key[3]

    def update_res_int(self, res_int_key, int_score):
        self.incr_res_int(res_int_key)
        self.update_res_int_score(res_int_key, int_score)

    def update_comb_int(self, comb_int_key, res_int_key):
        nr_int = self.res_int_dic[res_int_key].nr_int
        int_score = self.res_int_dic[res_int_key].int_score
        self.res_int_dic[comb_int_key].nr_int += nr_int
        self.res_int_dic[comb_int_key].int_score += int_score

    def new_comb_int(self, comb_int_key, res_int_key):
        res_id1 = self.get_res_id1_from_key(comb_int_key)
        res_id2 = self.get_res_id2_from_key(comb_int_key)
        int_type = self.get_int_type_from_key(comb_int_key)
        int_subtype = self.get_int_subtype_from_key(comb_int_key)
        nr_int = self.res_int_dic[res_int_key].nr_int
        int_score = self.res_int_dic[res_int_key].int_score
        comb_int_attr_obj = nc_interactions.IntAttr(res_id1, res_id2, int_type, int_subtype)
        self.res_int_dic[comb_int_key] = comb_int_attr_obj
        self.res_int_dic[comb_int_key].nr_int = nr_int
        self.res_int_dic[comb_int_key].int_score = int_score
        self.res_int_lst.append(comb_int_key)

    def update_res_int_score(self, res_int_key, new_int_score):
        res_int_attr_obj = self.res_int_dic[res_int_key]
        res_int_attr_obj.int_score += new_int_score
        self.res_int_dic[res_int_key] = res_int_attr_obj

    def incr_res_int(self, res_int_key):
        res_int_attr_obj = self.res_int_dic[res_int_key]
        res_int_attr_obj.incr_int()

    def get_combine_int(self):
        res_int_lst = copy.copy(self.res_int_lst)
        for res_int_key in res_int_lst:
            res_int_attr_obj = self.res_int_dic[res_int_key]
            int_score = res_int_attr_obj.int_score
            (comb_int_key, comb_int_rkey) = self.get_comb_int_keys(res_int_key)
            if comb_int_key in self.res_int_dic:
                self.update_comb_int(comb_int_key, res_int_key)
            elif comb_int_rkey in self.res_int_dic:
                self.update_comb_int(comb_int_rkey, res_int_key)
            else:
                self.new_comb_int(comb_int_key, res_int_key)

    def get_comb_int_keys(self, res_int_key):
        res_id1 = self.get_res_id1_from_key(res_int_key)
        res_id2 = self.get_res_id2_from_key(res_int_key)
        int_type = 'combi'
        int_sybtype = ('all', 'all')
        combi_int_key = self.get_res_int_key(res_id1, res_id2, int_type, int_sybtype)
        combi_int_rkey = self.get_res_int_key(res_id2, res_id1, int_type, int_sybtype)
        return (combi_int_key, combi_int_rkey)
