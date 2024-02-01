#!/usr/bin/python

# Copyright 2014 Max-Planck-Institut Informatik
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

import math
import os
import string
import sys

alfa = 1.0
beta = 0.5
gamma = 2.0

# VALDAR ALGORITHM-------------------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------------------------------------------------
# w_i=1/(number of position) * sum_x(1/*(number of symbol types in position * frequency of i_th amino acidat the position x))
def calculate_wi(allignm, p, matrx):
    wi = 0.0
    for j in range(1, p + 1):
        kx = 0
        nxi = 0
        for k in range(24):
            if matrx[j][k] != 0:
                kx = kx + 1
            if matrx[0][k] == allignm[j - 1]:
                nxi = matrx[j][k]
        if kx * nxi != 0:
            wi = wi + 1.0 / float(kx * nxi)
    if p == 0:
        return 0
    else:
        return wi / p


# -----------------------------------------------------------------------------------------------------------------------------------------------
# according artickle "Scoring residue conservation", Valdar (p.238)
def calculate_scores(matrx, seq, p, lines, gap_format):

    C_trident = [0 for i5 in range(p)]

    # BLOSSOM32 matrix for calculating r(x)
    XX = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
          [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4],
          [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4],
          [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4],
          [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4],
          [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4],
          [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
          [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4],
          [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4],
          [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4],
          [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4],
          [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4],
          [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4],
          [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4],
          [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4],
          [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4],
          [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4],
          [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4],
          [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4],
          [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4],
          [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4],
          [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
          [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4],
          [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]]

    wi = [0 for i1 in range(seq)]
    s = 0
    kkk = 0
    number_position = 0  # position without "-"

    allignm = [[0 for i3 in range(p)] for i4 in range(seq)]
    allignm1 = [0 for j in range(p)]

    # ................................ Calculating scores (wi).............................

    for i in range(0, len(lines)):
        line = lines[i]
        if line[0] == '>':
            k = 0
            if kkk == 1:
                wi[s] = calculate_wi(allignm1, p, matrx)
                for iii in range(p):
                    allignm[s][iii] = allignm1[iii]
                s = s + 1
                allignm1 = [0 for j in range(p)]
            kkk = 0
        else:
            kkk = 1
            for j in range(0, len(line)):
                if line[j] in string.digits or line[j] in string.ascii_letters or line[j] in string.punctuation:
                    if line[j] == gap_format:
                        allignm1[k] = '-'
                    elif line[j] == 'J':
                        allignm1[k] = 'L'
                    else:
                        allignm1[k] = line[j]
                    k = k + 1
        if i == len(lines) - 1:
            wi[s] = calculate_wi(allignm1, p, matrx)
            for iii in range(p):
                allignm[s][iii] = allignm1[iii]

    for x in range(p):  # calculate score for each position of the allignments

        # ................................ Calculating palfa.............................
        palfa = [0.00 for i in range(24)]
        for i in range(24):
            for sss in range(seq):
                if matrx[0][i] == allignm[sss][x]:
                    palfa[i] = palfa[i] + wi[sss]

        # ................................ Calculating tx.............................
        # t(x)=lambda_t*sum_alfa(p_alfa*log_2(p_alfa))

        K = 24
        tx = 0.0
        for i in range(24):
            if palfa[i] != 0:
                tx = tx + palfa[i] * math.log(palfa[i], 2)
        if math.log(min(seq, K), 2) == 0:
            lt = 0
        else:
            lt = -1.0 * 1 / math.log(min(seq, K), 2)
        tx = float(tx) * float(lt)

        # ................................ Calculating gx.............................
        # g(x)=fraction of the symbol "-"

        gx = 0.0
        if seq != 0:
            gx = float(matrx[x + 1][23]) / float(seq)

        # ................................ Calculating rx.............................
        x_x = [0.0 for ii in range(24)]
        kx = 0.0
        mmax = -100
        mmin = 100
        for i2 in range(24):
            if matrx[x + 1][i2] != 0:
                kx = kx + 1
                for i3 in range(24):
                    x_x[i3] = x_x[i3] + XX[i2][i3]
                    if mmax < XX[i2][i3]:
                        mmax = XX[i2][i3]
                    if mmin > XX[i2][i3]:
                        mmin = XX[i2][i3]
        for i3 in range(24):
            if kx != 0:
                x_x[i3] = float(x_x[i3]) / float(kx)

        rx = 0.0
        for i2 in range(24):
            aaa = 0.0
            if matrx[x + 1][i2] != 0:
                for i3 in range(24):
                    aaa = aaa + float(math.pow((x_x[i3] - XX[i2][i3]), 2))
            if aaa > 0:
                rx = rx + float(math.sqrt(aaa))
        lr = math.sqrt(20 * math.pow(mmax - mmin, 2))
        if float(kx) * float(lr) != 0:
            rx = float(rx) / float((kx * lr))

        # ................................ Final calculating C_trident.............................

        C_trident[x] = math.pow(1 - tx, alfa) * math.pow(1 - rx, beta) * math.pow(1 - gx, gamma)
        if allignm[0][x] != "-":
            number_position = number_position + 1

    # ..............calculating vector with scores for the first(main) sequence without spaces.....
    C_trident_without_space = [0 for j in range(number_position)]

    pp = 0
    for x in range(p):
        if allignm[0][x] != "-":
            C_trident_without_space[pp] = C_trident[x]
            pp = pp + 1

    return C_trident, C_trident_without_space, number_position


# -----------------------------------------------------------------------------------------------------------------------------------------------
def matrixx(gap_format, output_format, path_alignment, path_out, path_log, path_id):
    path1 = path_alignment
    seq = 0
    KK = 0
    file_log = open(path_log, 'w')
    file_output = open(path_out, 'w')

    if os.path.exists(path1) is True:
        file1 = open(path1, 'rb')
        lines = file1.readlines()
        for i in range(1, len(lines)):
            line = lines[i]
            if line[0] == ">":
                break
            else:
                for j in range(0, len(line)):
                    if line[j] in string.digits or line[j] in string.ascii_letters or line[j] in string.punctuation:
                        KK = KK + 1
                    else:
                        break
        file1.close()
    else:
        file_log.write("Alignment file is not found\n")

    id_c = 0
    if os.path.exists(path_id) is True:
        file_id = open(path_id, 'r')
    else:
        if output_format == "name+score":
            id_c = 1
            file_log.write("File with identifiers is not found\n")

    # ..................................Calculating frequency matrix..................................................

    output_m = [[0 for c in range(24)] for r in range(KK + 1)]
    output_m[0] = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-']
    if os.path.exists(path1) is True and id_c != 1:
        file1 = open(path1, 'rb')
        lines = file1.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if line[0] == ">":
                seq = seq + 1
                k = 0  # new line
            else:
                for j in range(0, len(line)):
                    if line[j] in string.digits or line[j] in string.ascii_letters or line[j] in string.punctuation:
                        k = k + 1
                        if k > KK:
                            file_log.write("Something wrong with normalization of alignments: Lengths of sequences are not identical\n")
                        if line[j] == str(gap_format):
                            symbol = '-'
                        elif line[j] == 'J':
                            symbol = 'L'
                        else:
                            symbol = line[j]
                        check = 0
                        for r in range(24):
                            if symbol == output_m[0][r]:
                                output_m[k][r] = output_m[k][r] + 1
                                check = 1
                        if check == 0:
                            file_log.write("Something wrong with an alfabet - unexpected symbol was found: " + str(symbol) + "\n")
        file1.close()
        file_log.write("Number of sequences: " + str(seq) + '\n')
        file_log.write("Number of positions: " + str(k) + '\n')

        # ...................................Calculating scores..............................................................
        try:
            (C_trident, C_trident_2, ppp) = calculate_scores(output_m, seq, k, lines, str(gap_format))
            file_log.write("Number of nonempty positions in the first sequence: " + str(ppp) + '\n')
            if output_format == "name+score":
                lll_id = file_id.readlines()
                if len(lll_id) != ppp:
                    file_log.write("Length of the alignment and number of identifiers are not identical. Please check.\n")
                else:
                    for i in range(ppp):
                        file_output.write(str(lll_id[i]).strip() + ' ' + str(C_trident_2[i]) + '\n')
                    file_log.write("Results for the scores were written to the file " + str(path_out) + '\n')
                    file_log.write("###########Calculation is succesfully finished###############" + '\n')
                    file_id.close()
            elif output_format == "resid+score":
                for i in range(ppp):
                    file_output.write(str(int(i) + 1) + ' ' + str(C_trident_2[i]) + '\n')
                file_log.write("Results for the scores were written to the file " + str(path_out) + '\n')
                file_log.write("###########Calculation is succesfully finished###############" + '\n')
            elif output_format == "score":
                for i in range(ppp):
                    file_output.write(str(C_trident_2[i]) + '\n')
                file_log.write("Results for the scores were written to the file " + str(path_out) + '\n')
                file_log.write("###########Calculation is succesfully finished###############" + '\n')
            else:
                file_log.write("Parameter output_format wasn't defined correctly\n")
        except:
            file_log.write("Unexpected error: " + str(sys.exc_info()[1]) + '\n')

    file_log.close()
    file_output.close()

    # -----------------------------------------------------------------------------------------------------------------------------------------------
    # END VALDAR ALGORITHM---------------------------------------------------------------------------------------------------------------------------

    #matrixx(".", "name+score", "2ypi_dTIM.afasta", "out.txt", "log.txt", "2YPI_h_res.txt")
    '''
    matrixx(gap_format, output_format, path_alignment, path_out, path_log, path_id)
    gap_format: any symbol as "gap_format" will be considered as a gap
    output_format:
    - "name+score"
    - "resid+score"
    - "score"

    Inputs: 2 files     path_alignment, path_id
    - with multiple sequence alignment - FASTA format
    - with nodes identifiers (according to the NONgaps positions in the first sequence in alignment) - TXT format: shouldn't contain empty lines
    !!! if you don't have file with identifiers just put empty string as path_id

    Output: 2 files     path_out, path_log
    - file with 2 columns: identifire conservation_score (space separated)
    - log file

    Parameters of calculation:
        alfa=1
        beta=0.5
        gamma=2
    '''


if len(sys.argv) == 7:
    matrixx(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
elif len(sys.argv) == 6:
    matrixx(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], "")
else:
    print('''
    Usage: get_conservation.py gap_format output_format path_alignment path_out path_log [path_id]

    gap_format: any symbol as "gap_format" will be considered as a gap
    output_format:
    - "name+score": identifier conservation_score
    - "resid+score": index conservation_score
    - "score": conservation_score

    Input: 2 files
    path_alignment: file with multiple sequence alignment in FASTA format
    path_id: [optional] file with nodes identifiers (according to the NONgaps positions in the first sequence in alignment) in TXT format (shouldn't contain empty lines)

    Output: 2 files
    path_out: file with 1 or 2 columns (according to output_format) separated by space
    path_log: log file

    Parameters of calculation:
    alfa=1
    beta=0.5
    gamma=2
    ''')
    sys.exit()
