#!/usr/bin/python3

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

import argparse
import os
import sys

from . import nc_interactions, nci_graph, nci_residues, sel_utils, st_selection


def main(input_path, output_path, chains_file, ligand_file, directed, reduce_cmd, probe_cmd, pdb_name = None, custom_db_path = False):

    # check if single or multiple files are given as input
    single_use = None
    if input_path.endswith("/") and os.access(input_path, os.R_OK):
        single_use = False
    elif os.path.exists(input_path):
        single_use = True

    if single_use is None:
        print("Cannot read input file or directory. Exiting...")
        sys.exit()

    ### read in chains

    file_content = open(chains_file, 'r').read()
    chains = [c.strip() for c in file_content.rstrip().split(",")]
    if len(chains) == 0:
        print("Cannot read chains file. Exiting...")
        sys.exit()

    ### read in ligands
    ligands = []

    if ligand_file:
        ligand_file = ligand_file
        file_content = open(ligand_file, 'r').read().rstrip().split("\n")
        for line in file_content:
            ligands.append([l.strip() for l in line.split(",")])

    # create RINs
    if single_use:
        # set file names
        if not pdb_name == None:
            pdb_filename = pdb_name
            # input_pdb_parts = input_path.split("/")
            # input_path = "/".join(input_pdb_parts[:len(input_pdb_parts) - 1])
            
        else:
            input_pdb_parts = input_path.split("/")
            input_path = "/".join(input_pdb_parts[:len(input_pdb_parts) - 1])
            pdb_filename = input_pdb_parts[len(input_pdb_parts) - 1]
        print("Processing single pdb file " + pdb_filename + "...")

        # create slection list
        if not pdb_name == None:
            sel_id = pdb_name
        else:
            sel_id = pdb_filename[:-4]

        selection_lst = sel_utils.chains_to_sel_list(sel_id, chains)

        ligand_lst = sel_utils.ligands_to_sel_list(ligands)

        if len(ligand_lst) > 0:
            selection_lst += ligand_lst
        
        # run reduce and probe
        if custom_db_path:
            is_custom_db = True
            (pdb_h_filename, probe_filename) = nc_interactions.get_reduce_probe_rsl(pdb_filename, input_path, output_path, output_path, reduce_cmd, probe_cmd, is_custom_db)
        
        else:
            (pdb_h_filename, probe_filename) = nc_interactions.get_reduce_probe_rsl(pdb_filename, input_path, output_path, output_path, reduce_cmd, probe_cmd)
        sif_file = pdb_h_filename[:-4] + ".sif"

        # get network sif file
        (stsel_obj, nci_obj, nci_res, nci_graph_) = nci_graph.pdb_to_sif(sel_id, pdb_h_filename, output_path, selection_lst, probe_filename, output_path, sif_file, directed=directed)
        
    else:
        print("Processing multiple pdb files in directory " + input_path)
        for pdb_filename in os.listdir(input_path):
            if pdb_filename.endswith(".pdb") or pdb_filename.endswith(".ent"):
                # create selection list
                sel_id = pdb_filename[:-4]
                selection_lst = sel_utils.chains_to_sel_list(sel_id, chains)
                ligand_lst = sel_utils.ligands_to_sel_list(ligands)
                if len(ligand_lst) > 0:
                    selection_lst += ligand_lst

                # check if probe and reduce were already run on this file
                pdb_h_filename = pdb_filename[:-4] + "_h.ent"
                probe_filename = pdb_filename[:-4] + "_h.probe"
                if not os.path.exists(output_path + pdb_h_filename) or not os.path.exists(output_path + probe_filename):
                    # run reduce and probe
                    nc_interactions.get_reduce_probe_rsl(pdb_filename, input_path, output_path, output_path, reduce_cmd, probe_cmd)

                # get network sif file
                sif_file = pdb_h_filename[:-4] + ".sif"
                if not os.path.exists(sif_file):
                    nci_graph.pdb_to_sif(sel_id, pdb_h_filename, output_path, selection_lst, probe_filename, output_path, sif_file, directed=directed)


if __name__ == "__main__":
    # set program paths
    reduce_cmd = os.path.dirname(sys.argv[0]) + '/reduce'  # reduce command
    probe_cmd = os.path.dirname(sys.argv[0]) + '/probe'  # probe command

    parser = argparse.ArgumentParser(description='Create residue interaction network for chains', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_path', metavar='path_pdb', help='PDB file or directory containing PDB files of the same protein structure')
    parser.add_argument('output_path', metavar='path_output', help='directory to save generated files')
    parser.add_argument('chains_file', metavar='path_chains', help='''file with chain identifiers separated by commas
    - chain identifier might be any letter or an empty string
    - examples: "A,B,I" or "" or "A,"''')
    parser.add_argument('ligand_file', nargs='?', metavar='path_ligands', help='''file with ligand identifiers separated by new lines
    - each ligand is listed with a name, a chain identifier and a residue number separated by commas
    - examples: "NOA,I,1" or "NOA,,1"''')
    parser.add_argument('-d', '--directed', action="store_true", help='generate directed interaction types.\n\'mc,sc\' with this option indicates an interaction of the main chain of the first with the side chain of the second residue instead of the direction being ambiguous')
    args = parser.parse_args()

    print("Processing input...")

    # set paths for files to read/write
    input_path = args.input_path
    output_path = args.output_path

    directed = args.directed
    chains_file = args.chains_file
    ligand_file = args.ligand_file

    main(input_path, output_path, chains_file, ligand_file, directed, reduce_cmd, probe_cmd)
