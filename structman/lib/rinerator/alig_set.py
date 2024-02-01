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


def get_alig_from_fasta(entry_label_lst, fasta_alig_file, stsel_dic, seq_dic):
    alig_set = AligSet(entry_label_lst)
    alig_set.read_fasta_alig(fasta_alig_file, stsel_dic, seq_dic)
    print("INFO: Read aliment from %s" % fasta_alig_file)
    return alig_set


class AligSet:
    def __init__(self, entry_label_lst):
        self.entry_label_lst = entry_label_lst  # list of entry_labels, entries are either sequences or structures
        self.alig_idx_lst = []  # alig_idx_lst[alig_idx] = map_alig_idx_dic
        #                    map_alig_idx_dic[entry_label] = seq_idx
        #                    entry_label is seq_id for sequences or
        #                    entry_label is st_sel.sel_id for structures
        self.alig_entry_dic = {}  # alig_entry_dic[entry_label] = map_seq_idx_lst
        #                          map_seq_idx_lst[seq_idx] = alig_idx
        self.alig_len = None

    def clear_alig(self):
        self.alig_idx_lst = []
        self.alig_entry_dic = {}
        self.alig_len = None

    def read_fasta_alig(self, fasta_alig_file, seq_dic, stsel_dic):
        gap_symbol = '-'
        fobj = file(fasta_alig_file, 'r')
        lines = fobj.readlines()
        fobj.close()

        current_entry = None
        current_seq_idx = 0
        current_alig_idx = 0
        current_seq_obj = None
        for line_str in lines:
            line_str = line_str.strip()
            if line_str == '':
                continue
            if line_str[0] == '>':
                entry_label = line_str[1:].strip().split()[0]
                entry_label = entry_label.split('/')[0]
                if entry_label in self.alig_entry_dic:
                    print("ERROR: Repeated entry %s in FASTA alignment %s" % (entry_label, fasta_alig_file))
                    self.clear_alig()
                    return
                if self.entry_label_lst.count(entry_label) == 1:
                    self.alig_entry_dic[entry_label] = []
                else:
                    print("WARNING: Skipping %s in file %s" % (entry_label, fasta_alig_file))
                    current_entry = 'Skip'
                    current_seq_idx = 0
                    current_alig_idx = 0
                    continue
                if entry_label in stsel_dic:
                    stsel_obj = stsel_dic[entry_label]
                    current_seq_obj = stsel_obj.seq_obj
                elif entry_label in seq_dic:
                    current_seq_obj = seq_dic[entry_label]
                else:
                    print("ERROR: No entry %s from alig %s in set" % (entry_label, fasta_alig_file))
                    self.clear_alig()
                    return
                current_entry = entry_label
                current_seq_idx = 0
                current_alig_idx = 0
            elif current_entry is None:
                print("ERROR: Cannot read FASTA alignment %s" % fasta_alig_file)
                self.clear_alig()
                return
            elif current_entry == 'Skip':
                continue
            else:
                for char_counter in range(len(line_str)):
                    if current_alig_idx >= len(self.alig_idx_lst):
                        map_alig_idx_dic = {}
                        self.alig_idx_lst.append(map_alig_idx_dic)
                    map_alig_idx_dic = self.alig_idx_lst[current_alig_idx]
                    if entry_label in map_alig_idx_dic:
                        print("ERROR: Entry %s, sequence mismatch in aligment %s, seq_idx %d, alig_idx %s" % (entry_label, fasta_alig_file, current_seq_idx, current_alig_idx))
                    res_code = line_str[char_counter]
                    if res_code != gap_symbol:
                        if res_code != current_seq_obj.seq_str[current_seq_idx]:
                            print("ERROR: Entry %s, sequence mismatch in aligment %s, seq_idx %d, alig_idx %s" % (entry_label, fasta_alig_file, current_seq_idx, current_alig_idx))
                        self.alig_entry_dic[entry_label].append(current_alig_idx)
                        map_alig_idx_dic[entry_label] = current_seq_idx
                        current_seq_idx += 1
                        current_alig_idx += 1
                    else:
                        map_alig_idx_dic[entry_label] = None
                        current_alig_idx += 1
        self.alig_len = current_alig_idx

    def write_alig(self, alig_file, seq_dic, stsel_dic):
        fobj = file(alig_file, 'w')
        for alig_idx in range(self.alig_len):
            alig_pos = alig_idx + 1
            map_alig_idx_dic = self.alig_idx_lst[alig_idx]
            for entry_label in self.entry_label_lst:
                seq_idx = map_alig_idx_dic[entry_label]
                if entry_label in seq_dic:
                    entry_type = 'seq'
                    res_label = '_'
                    if seq_idx is None:
                        seq_idx = '-'
                        seq_pos = '-'
                        res_name = '-'
                    else:
                        seq_pos = seq_idx + 1
                        res_name = seq_dic[entry_label].seq_str[seq_idx]
                    out_str = self.get_alig_line(alig_pos, entry_label, entry_type, seq_pos, res_label, res_name)
                    fobj.write(out_str)
                elif entry_label in stsel_dic:
                    stsel_obj = stsel_dic[entry_label]
                    entry_type = 'str'
                    if seq_idx is None:
                        seq_idx = '-'
                        seq_pos = '-'
                        res_name = '-'
                        res_label = '-'
                    else:
                        seq_pos = seq_idx + 1
                        res_name = stsel_obj.seq_obj.seq_str[seq_idx]
                        res_label = stsel_obj.get_res_label_from_seq_idx(seq_idx)
                    out_str = self.get_alig_line(alig_pos, entry_label, entry_type, seq_pos, res_label, res_name)
                    fobj.write(out_str)
        fobj.close()
        print("INFO: Wrote alignment %s" % alig_file)

    def get_alig_line(self, alig_pos, entry_label, entry_type, seq_pos, res_label, res_name):
        out_str = "%-7d > %-20s > %3s > %7s > %12s > %1s\n" % (alig_pos, entry_label, entry_type, str(seq_pos), res_label, res_name)
        return out_str
