import math
import os
import os.path
import subprocess
import time
import gzip
import shutil

from . import st_selection

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


def get_ncint_probe(probe_filename, probe_path, stsel_obj, directed):
    probe_file = os.path.join(probe_path, probe_filename)
    nci_obj = NCI()
    nci_obj.get_ncint_probe(probe_file, stsel_obj, directed)
    print("INFO: Read non covalent interactions in PROBE file %s" % probe_filename)
    return nci_obj


def get_reduce_probe_rsl(pdb_filename, pdb_path, pdb_h_path, probe_path, reduce_cmd, probe_cmd, is_custom_db=False):
    #if is_custom_db:
    pdb_file = pdb_path
    
    #else:
    #    pdb_file = os.path.join(pdb_path, pdb_filename)
    pdb_h_filename = pdb_filename.split('.')[0] + '_h' + '.ent'
    pdb_h_file = os.path.join(pdb_h_path, pdb_h_filename)
    probe_filename = pdb_h_filename.split('.')[0] + '.probe'
    probe_file = os.path.join(probe_path, probe_filename)

    # temp_folder_path = pdb_h_path.split("/")[:-1]
    # print(f'temp_folder_path = {temp_folder_path}')

    # if pdb_file.endswith(".gz"):
    #     os.system(f"tar -xvf {pdb_file} -C ")

    run_reduce(pdb_file, pdb_h_file, reduce_cmd)
    run_probe(pdb_h_file, probe_file, probe_cmd)
    return (pdb_h_filename, probe_filename)


def get_probe_rsl(pdb_filename, pdb_path, pdb_h_path, probe_path, probe_cmd):
    pdb_file = os.path.join(pdb_path, pdb_filename)
    probe_filename = pdb_filename.split('.')[0] + '.probe'
    probe_file = os.path.join(probe_path, probe_filename)

    run_probe(pdb_file, probe_file, probe_cmd)
    return (pdb_filename, probe_filename)


def run_reduce(pdb_file, pdb_h_file, reduce_cmd):
    if pdb_file.endswith('.gz'):
        print(f'gzip pdb_file = {pdb_file}')
        in_pipe = gzip.open(pdb_file, 'rb')
    else:
        print(f'nongzip pdb_file = {pdb_file}')
        in_pipe = open(pdb_file, 'rb')
    
    cmd = ['reduce']
    if os.path.exists(reduce_cmd):
        cmd[0] = reduce_cmd
    cmd += ['-BUILD', '-']

    print("INFO: %s" % ' '.join(cmd))

    f = open(pdb_h_file, 'w')
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=f, stderr=open(os.devnull, 'w'))
    #shutil.copyfileobj(in_pipe, p.stdin)
    p.communicate(in_pipe.read())
    f.close()
    in_pipe.close()

    if p.returncode not in [0, 1]:
        raise IOError("reduce return code: %s" % str(p.returncode))


def run_probe(pdb_h_file, probe_file, probe_cmd):
    # additional options used in MolProbity: -4H -gap
    cmd = ['probe']
    if os.path.exists(probe_cmd):
        cmd[0] = probe_cmd
    cmd += ['-unformated', '-MC', '-MINOCCupancy0.0', '-self', 'all', pdb_h_file]

    print("INFO: %s" % ' '.join(cmd))

    f = open(probe_file, 'w')
    p = subprocess.Popen(cmd, stdout=f, stderr=open(os.devnull, 'w'))
    p.communicate()
    f.close()

    if p.returncode != 0:
        raise IOError("probe return code: %s" % str(p.returncode))


class NCI:
    def __init__(self):
        self.atom_int_dic = {}  # dictionary of atomic interactions, dic[atom_int_key] = atom_int_attr_obj
        #                      atom_int_key=(atom_label1,atom_label2,int_type)
        self.atom_int_lst = []  # list with keys of atom_int_dic, atom_int_lst = [atom_int_key,atom_int_key...]
        self.probe_radius = 0.25  # 0.25A is default parameter in probe

    def get_ncint_probe(self, probe_file, stsel_obj, directed):
        fobj = open(probe_file, 'r')
        lines = fobj.readlines()
        fobj.close()

        for line in lines:
            if self.check_line_probe(line) != True:
                continue
            line_lst = line.split(':')
            int_type_probe = line_lst[2]
            int_type = self.translate_probe_ncint_type(int_type_probe)
            if int_type is None:
                continue
            atom_desc1 = line_lst[3]
            atom_desc2 = line_lst[4]
            atom_label1 = self.atom_desc_to_atom_label(atom_desc1, stsel_obj)
            atom_label2 = self.atom_desc_to_atom_label(atom_desc2, stsel_obj)

            if (atom_label1 is None) or (atom_label2 is None):
                continue
            # raw score
            raw_score = float(line_lst[11])
            # get keys for atom_int_dic
            atom_int_key = self.get_atom_int_key(atom_label1, atom_label2, int_type)
            atom_int_rkey = self.get_atom_int_key(atom_label2, atom_label1, int_type)
            if atom_int_key in self.atom_int_dic:
                self.update_int(atom_int_key, raw_score)
                continue
            elif atom_int_rkey in self.atom_int_dic:
                self.update_int(atom_int_rkey, raw_score)
                continue
            # new atom_int
            self.new_int(atom_int_key, stsel_obj, directed)
            self.update_int(atom_int_key, raw_score)

    def new_int(self, atom_int_key, stsel_obj, directed):
        (atom_label1, atom_label2, int_type) = atom_int_key
        attr_obj = IntAttr(atom_label1, atom_label2, int_type, None)
        attr_obj.set_mc_sc_ligand(stsel_obj, directed)
        self.atom_int_dic[atom_int_key] = attr_obj
        self.atom_int_lst.append(atom_int_key)
        # combint_attr_obj = IntAttr(atom_label1,atom_label2,'combint', None)
        # combint_attr_obj.set_combint_subtype()
        # combint_key = self.get_combint_key(atom_int_key)
        # self.atom_int_dic[combint_key]=combint_attr_obj
        # self.atom_int_lst.append(combint_key)

    def check_line_probe(self, line):
        if len(line) >= 1:
            if line[0] == ':':
                return True
        return False

    def atom_desc_to_atom_label(self, atom_desc, stsel_obj):
        chain_id = atom_desc[1:2]
        res_seq = int(atom_desc[2:6])
        icode = atom_desc[6:7]
        resname = atom_desc[7:10]
        atom_name = atom_desc[11:15].strip()
        altloc = atom_desc[15:16]
        atom_label = stsel_obj.get_atom_label_from_subid(chain_id, res_seq, icode, resname, atom_name, altloc)
        return atom_label

    def translate_probe_ncint_type(self, ncint_type_probe):
        if ncint_type_probe == 'cc':
            ncint_type = 'cnt'  # close contact
        elif ncint_type_probe == 'wc':
            ncint_type = 'cnt'  # wide contact
        elif ncint_type_probe == 'so':
            ncint_type = 'cnt'  # small overlap as contact
        elif ncint_type_probe == 'bo':
            ncint_type = 'ovl'  # big overlap
        elif ncint_type_probe == 'hb':
            ncint_type = 'hbond'  # hydrogen bond
        else:
            ncint_type = None
        return ncint_type

    def update_int(self, atom_int_key, raw_score):
        attr_obj = self.atom_int_dic[atom_int_key]
        attr_obj.int_score = attr_obj.int_score + raw_score
        self.atom_int_dic[atom_int_key] = attr_obj
        # combint_key = self.get_combint_key(atom_int_key)
        # combint_attr_obj = self.atom_int_dic[combint_key]
        # combint_attr_obj.int_score = combint_attr_obj.int_score + raw_score
        # self.atom_int_dic[combint_key] = combint_attr_obj

    def get_combint_key(self, atom_int_key):
        (atom_label1, atom_label2) = self.get_atom_labels_from_key(atom_int_key)
        combint_key = self.get_atom_int_key(atom_label1, atom_label2, 'combint')
        return combint_key

    def get_atom_int_key(self, atom_label1, atom_label2, int_type):
        return (atom_label1, atom_label2, int_type)

    def get_atom_labels_from_key(self, atom_int_key):
        return (atom_int_key[0], atom_int_key[1])

    def atom_int_key_to_atom_label1(self, atom_int_key):
        return atom_int_key[0]

    def atom_int_key_to_atom_label2(self, atom_int_key):
        return atom_int_key[1]

    def atom_int_key_to_int_type(self, atom_int_key):
        return atom_int_key[2]

    def get_atom_int_attr(self, atom_int_key):
        return self.atom_int_dic[atom_int_key]


class IntAttr:
    def __init__(self, id1, id2, int_type, int_subtype):
        self.id1 = id1
        self.id2 = id2
        self.int_type = int_type  # cnt ovl hbond combint
        self.int_subtype = int_subtype
        self.nr_int = 1
        self.int_score = 0.0

    def incr_int(self):
        self.nr_int += 1

    def set_type(self, int_type):
        self.int_type = int_type

    def set_mc_sc_ligand(self, stsel_obj, directed=False):
        """
        undirected:
        ('mc','mc')
        ('mc','sc')
        ('sc','sc')
        ('mc','ligand')
        ('sc','ligand')
        ('ligand','ligand')
        ('all','all')
        directed:
        ('mc','mc')
        ('mc','sc')
        ('sc','mc')
        ('sc','sc')
        ('mc','ligand')
        ('ligand', 'mc')
        ('sc', 'ligand')
        ('ligand', 'sc')
        ('ligand','ligand')
        ('all','all')
        """
        atom_label1 = self.id1
        atom_label2 = self.id2
        subtype1 = stsel_obj.check_mc_sc_ligand_from_label(atom_label1)
        subtype2 = stsel_obj.check_mc_sc_ligand_from_label(atom_label2)

        if subtype1 is None or subtype2 is None:
            self.int_subtype = None
            return

        if directed:
            self.int_subtype = (subtype1, subtype2)
            return

        if subtype1 == 'mc':
            int_subtype = (subtype1, subtype2)
        elif subtype1 == 'sc':
            if subtype2 == 'mc':
                int_subtype = (subtype2, subtype1)
            else:
                int_subtype = (subtype1, subtype2)
        elif subtype1 == 'ligand':
            int_subtype = (subtype2, subtype1)
        self.int_subtype = int_subtype
