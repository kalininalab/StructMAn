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

import os
import sys

from . import ext_data

labels = {
    "combi unweighted": "number of 'combi' edges (residue interactions)",
    "hbond unweighted": "number of 'hbond' edges (residue interactions)",
    "ovl unweighted": "number of 'overlap' edges (residue interactions)",
    "cnt unweighted": "number of 'contact' edges (residue interactions)",
    "combi intsc": "sum of 'combi' interaction scores",
    "hbond intsc": "sum of 'hbond' interaction scores",
    "ovl intsc": "sum of 'overlap' interaction scores",
    "cnt intsc": "sum of 'contact' interaction scores",
    "combi nrint": "number of 'combi' edges (atomic interactions)",
    "hbond nrint": "number of 'hbond' edges (atomic interactions)",
    "ovl nrint": "number of 'overlap' edges (atomic interactions)",
    "cnt nrint": "number of 'contact' edges (atomic interactions)",
    "JANJ780101": "AA average accessible surface area",
    "GRAR740102": "AA polarity",
    "JURD980101": "AA modified Kyte-Doolittle hydrophobicity scale",
    "ISOY800101": "AA normalized relative frequency of alpha-helix",
    "GRAR740103": "AA volume",
    "JOND920101": "AA relative frequency of occurrence",
    "FAUJ880112": "AA negative charge",
    "FAUJ880111": "AA positive charge",
    "KLEP840101": "AA net charge",
    "ZIMJ680104": "AA isoelectric point",
    "JOND920102": "AA relative mutability",
    "ConSurfDB": "conservation scores from the ConSurfDB"
}

# urls
consurfdb_url = "http://bental.tau.ac.il/new_ConSurfDB/DB/"
aaindex_url = "http://www.genome.jp/dbget-bin/www_bget?aaindex:"

# AAindex parameters
aaindex_sel = "KLEP840101,ZIMJ680104,JOND920102,JANJ780101,GRAR740102,JURD980101,ISOY800101,GRAR740103,JOND920101,FAUJ880112,FAUJ880111"

# files to look for
consurf_suf = ".grades"
sif_suf = ".sif"
easc_suf = "intsc.ea"
eanr_suf = "nrint.ea"

# check input
if len(sys.argv) != 4 and len(sys.argv) != 5:
    print('''
    USAGE: get_data.py path_id path_output pdb_id [path_input]

    Parameters:
    path_id: list of residue identifiers

    path_output: file to save retrieved data

    path_input: path with additional input data, any of the following files will be considered
    - consurf.grades: conservation scores from ConSurfDB
    - *.sif: RIN file
    - *_nrint.ea
    - *_intsc.ea

    Short description of retrieved data:
    ''')
    for param in sorted(labels.keys()):
        print("    " + param + "\t" + labels[param])
    sys.exit()

input_file = sys.argv[1]
output_file = sys.argv[2]
pdb_id = sys.argv[3]
add_input = False

if not os.path.exists(input_file):
    print("File with residue identifiers not found. Exiting...")
    sys.exit()

# list all files in input dir
sif_filename = ""
easc_filename = ""
eanr_filename = ""
consurf_filenames = []
if len(sys.argv) == 5:
    input_dir = sys.argv[4]
    if os.path.exists(input_dir):
        for filename in os.listdir(input_dir):
            if filename.endswith(consurf_suf):
                consurf_filenames.append(filename)
            elif filename.endswith(sif_suf) and sif_filename == "":
                sif_filename = filename
            elif filename.endswith(easc_suf) and easc_filename == "":
                easc_filename = filename
            elif filename.endswith(eanr_suf) and eanr_filename == "":
                eanr_filename = filename
    else:
        print("Input directory not found. Exiting...")
        sys.exit()

# parse node identifiers
input_res = open(input_file, 'r').read().rstrip().split("\n")
chains = []
for res in input_res:
    res_parts = res.split(":")
    chain = ""
    if len(res_parts) == 4:
        chain = res_parts[0]
    elif len(res_parts) == 5:
        chain = res_parts[1]
    if len(chain) > 0 and chain not in chains:
        chains.append(chain)

# dictionary for saving the data
all_params_dict = {}

# get AA indices
aaindices = ext_data.getAAIndices(aaindex_url, aaindex_sel.split(","))
ext_data.getResidueAAIndices(input_res, aaindices, all_params_dict)

# get conservation
if add_input:
    for filename in consurf_filenames:
        ext_data.getConservation(input_dir + consurf_file, chains, all_params_dict)
else:
    ext_data.getConservation(consurfdb_url + pdb_id.upper() + "/", chains, all_params_dict)

# get number of interactions (edges)
if sif_filename != "":
    ext_data.getInteractionProperties(input_dir + sif_filename, "unweighted", chains, all_params_dict)

# get interaction score
if easc_filename != "":
    ext_data.getInteractionProperties(input_dir + easc_filename, "intsc", chains, all_params_dict)

# get number of atomic interactions
if eanr_filename != "":
    ext_data.getInteractionProperties(input_dir + eanr_filename, "nrint", chains, all_params_dict)

# save data
ext_data.saveData(output_file, input_res, [all_params_dict], all_params_dict.keys())
