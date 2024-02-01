import os.path

from . import alig_set, sequence, st_selection

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


def get_prtset(pdb_selection_lst, pdb_path):
    prtset_obj = PrtSet()
    for pdb_selection_tpl in pdb_selection_lst:
        sel_id = pdb_selection_tpl[0]
        pdb_file_name = pdb_selection_tpl[1]
        selection_lst = pdb_selection_tpl[2]
        stsel_obj = st_selection.set_str_sel(sel_id, pdb_file_name, pdb_path, selection_lst)
        prtset_obj.get_stsel(stsel_obj)
    return prtset_obj


def get_stset_seqs(fasta_lst, seq_path, prtset_obj):

    for fasta_filename in fasta_lst:
        fasta_file = os.path.join(seq_path, fasta_filename)
        prtset_obj.get_seq(fasta_file)
        print("INFO: Read sequence %s" % fasta_file)
    return prtset_obj


def read_fasta_alig_write_alig(fasta_alig_file, alig_file, prtset_obj):
    prtset_obj.get_fasta_alig(fasta_alig_file)
    alig_set_obj = prtset_obj.alig_set_obj
    alig_set_obj.write_alig(alig_file, prtset_obj.seq_dic, prtset_obj.stsel_dic)
    return prtset_obj


# read structures, sequences
# write fasta

# read msa
# write alig


class PrtSet:
    def __init__(self):
        self.stsel_lst = []  # list of st_sel.sel_id
        self.stsel_dic = {}  # dic[st_sel.sel_id] = stsel_obj
        self.seq_lst = []  # list of sequence identifiers
        self.seq_dic = {}  # seq_dic[seq_id] = seq_obj
        self.alig_set_obj = None  # AligSet object

    def get_stsel(self, stsel_obj):
        if stsel_obj.sel_id in self.stsel_dic:
            print("WARNING: Skipping structure %s, already in set" % stsel_obj.sel_id)
            return
        if stsel_obj.seq_obj is None:
            stsel_obj.get_seq()
        self.stsel_lst.append(stsel_obj.sel_id)
        self.stsel_dic[stsel_obj.sel_id] = stsel_obj

    def get_seq(self, fasta_seq_file):
        seq_obj = sequence.Sequence()
        seq_obj.read_fasta(fasta_seq_file)
        self.seq_lst.append(seq_obj.seq_id)
        self.seq_dic[seq_obj.seq_id] = seq_obj

    def write_set_fasta(self, fasta_set_file):
        fobj = file(fasta_set_file, 'w')
        for seq_id in self.seq_lst:
            seq_obj = self.seq_dic[seq_id]
            fasta_seq_str = seq_obj.get_fasta_seq_str()
            fobj.write(fasta_seq_str)
        for sel_id in self.stsel_lst:
            stsel_obj = self.stsel_dic[sel_id]
            if stsel_obj.seq_obj is None:
                stsel_obj.get_seq()
            fasta_seq_str = stsel_obj.seq_obj.get_fasta_seq_str()
            fobj.write(fasta_seq_str)
        print("INFO: Wrote sequences in %s" % fasta_set_file)
        fobj.close()

    def get_fasta_alig(self, fasta_alig_file):
        entry_label_lst = []
        for entry_label in self.seq_lst:
            entry_label_lst.append(entry_label)
        for entry_label in self.stsel_lst:
            entry_label_lst.append(entry_label)
        alig_set_obj = alig_set.AligSet(entry_label_lst)
        alig_set_obj.read_fasta_alig(fasta_alig_file, self.seq_dic, self.stsel_dic)
        print("INFO: Read aliment from %s" % fasta_alig_file)
        self.alig_set_obj = alig_set_obj
