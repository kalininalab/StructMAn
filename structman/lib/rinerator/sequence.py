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


class Sequence:
    def __init__(self):
        self.seq_id = None
        self.seq_str = None
        self.seq_len = None
        self.seq_file = None

    def set_seq(self, seq_id, seq_str):
        self.seq_id = seq_id
        self.seq_str = seq_str
        self.seq_len = len(self.seq_str)

    def read_fasta(self, seq_file):
        fobj = file(seq_file)
        lines = fobj.readlines()
        fobj.close()

        len_lines = len(lines)

        seq_label = lines[0]
        if seq_label[0] != '>':
            print("ERROR: Sequence not in FASTA format %s" % seq_file)
            return
        seq_id = seq_label[1:].strip().split()[0]
        seq_lst = []
        for seq_line in lines[1:]:
            seq_line = seq_line.strip()
            seq_lst.append(seq_line)
        sep = ''
        seq_str = sep.join(seq_lst)
        self.seq_id = seq_id
        self.seq_str = seq_str
        self.seq_len = len(self.seq_str)
        self.seq_file = seq_file

    def get_fasta_seq_str(self, line_length=60):
        """
        FASTA format the sequence string with line return every line_length
        """
        out_seq_str = ""
        out_seq_str = ">" + self.seq_id + '\n'
        seq_str = self.seq_str
        seq_len = self.seq_len
        nr_full_lines = int(seq_len / line_length)
        remain_length = seq_len % line_length
        for line_c in range(nr_full_lines):
            start_idx = line_c * line_length
            end_idx = start_idx + line_length
            out_seq_str = out_seq_str + seq_str[start_idx:end_idx] + '\n'
        start_idx = nr_full_lines * line_length
        out_seq_str = out_seq_str + seq_str[start_idx:] + '\n'
        return out_seq_str

    def write_fasta_seq(self, fasta_file):
        fasta_seq_str = self.get_fasta_seq_str()
        fobj = file(fasta_file)
        fobj.write(fasta_seq_str)
        fobj.close()
        print("INFO: Wrote %s" % fasta_file)
