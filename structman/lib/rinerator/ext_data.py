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
from urllib.request import Request, URLError, urlopen

from . import st_selection


def getAAIndices(url, indices, aaindices={}):
    for index in indices:
        index = index.strip()
        print("Retrieve data for index " + index + "..."),
        req = Request(url + index)
        try:
            response = urlopen(req)
        except URLError as e:
            print(e)
            break
        else:
            html = response.read()
            lines = html.split("\n")
            line_num = 0
            for i in range(len(lines)):
                if lines[i].startswith("I"):
                    line_num = i
            aas = lines[line_num].split()[1:]
            indices1 = lines[line_num + 1].split()
            indices2 = lines[line_num + 2].split()
            aaindex_dict = {}
            for a in range(len(aas)):
                (aa1, aa2) = aas[a].split('/')
                aaindex_dict[aa1] = float(indices1[a])
                aaindex_dict[aa2] = float(indices2[a])
        aaindices[index] = aaindex_dict
        print("done.")
    return aaindices


def getResidueAAIndices(nodes, aaindices, aaindices_res={}):
    for aaindex in aaindices:
        aaindex_dict = {}
        for node in nodes:
            node_name = node.split(":")
            if len(node_name) != 4:
                print("format wrong for residue " + node)
                continue
            if node_name[3] not in st_selection.aa_trans_dic:
                print("missing aa mapping for " + node_name[3])
                continue
            aa = st_selection.aa_trans_dic[node_name[3]]
            if aa not in aaindices[aaindex]:
                print("aa " + aa + " not found.")
                continue
            aaindex_dict[node] = aaindices[aaindex][aa]
        aaindices_res[aaindex] = aaindex_dict
    return aaindices_res


def getConservation(url, chains, params={}):
    print("Retrieve data from ConSurfDB..."),
    conserv = {}
    check_file = False
    for chain in chains:
        req = Request(url + chain + "/consurf.grades")
        try:
            response = urlopen(req)
            html = response.read()
            lines = html.split("\n")
        except URLError as e:
            check_file = True
        except Exception as e:
            check_file = True
        if check_file:
            if os.path.exists(url):
                f = open(url, "r")
                lines = f.readlines()
                f.close()
            else:
                print("Could not get ConSurfDB data")
                continue
        line_num = 0
        for line in lines:
            words = line.split("\t")
            if len(words) != 14:
                continue
            score = float(words[3].strip())
            res_parts = words[2].split(":")
            resid = ""
            if len(res_parts) == 2 and len(res_parts[1]) == 1:
                res_chain = res_parts[1].strip()
                res = res_parts[0].strip()
                res_type = '_'
                res_number = '_'
                res_icode = '_'
                if len(res) > 3:
                    res_type = res[0:3]
                    res_number = res[3:]
                    if res_number.isdigit() != True and len(res_number) > 0:
                        res_number = res[3:-1]
                        res_icode = res[-1:]
                resid = res_chain + ":" + res_number + ":" + res_icode + ":" + res_type
            if len(resid) == 0:
                continue
            conserv[resid] = score
    if len(conserv) > 0:
        params["ConSurfDB"] = conserv
    print("done.")
    return params


def getInteractionProperties(filename, edge_weight, chains, int_params={}):
    print("Retrieve '" + edge_weight + "' interaction properties..."),
    network_file = open(filename, "r")
    line_number = 0
    nodes = set()
    edge_params = set()
    for line in network_file:
        line_number += 1
        line = line.strip("\n")
        if line == "" or (edge_weight != "unweighted" and line_number == 1):
            continue
        line_parts_tmp = line.split("\t")
        line_parts = line_parts_tmp[0].split()
        if len(line_parts) != 3:
            continue
        node1 = line_parts[0]
        node2 = line_parts[2]
        if len(line_parts_tmp) == 1:
            edge_type = line_parts[1].split(":")[0]
            weight = 1
        elif len(line_parts_tmp) == 2:
            edge_type = line_parts[1][1:-1].split(":")[0]
            weight = float(line_parts_tmp[1])
        else:
            continue
        chain1 = node1.split(":")[0]
        chain2 = node2.split(":")[0]
        if chain1 not in chains or chain2 not in chains:
            continue
        if node1 not in nodes:
            nodes.add(node1)
        if node2 not in nodes:
            nodes.add(node2)
        edge_param = edge_type + " " + edge_weight
        if edge_param not in edge_params:
            edge_params.add(edge_param)
            int_params[edge_param] = {}
        if node1 not in int_params[edge_param]:
            int_params[edge_param][node1] = 0
        if node2 not in int_params[edge_param]:
            int_params[edge_param][node2] = 0
        int_params[edge_param][node1] += weight
        int_params[edge_param][node2] += weight
    network_file.close()
    for node in nodes:
        for param in edge_params:
            if node not in int_params[param]:
                int_params[param][node] = 0
    print("done.")
    return int_params


# save a dictionary of residue properties, indexed by the property name and
# containing dictionaries of residues and their respective values
def saveData(filename, nodes, all_params, param_names):
    print("Save all parameters... "),
    results = []
    results.append("Node\t")
    append_nodes = True
    for params in all_params:
        results[0] += "\t".join(param_names) + "\t"
        for i in range(0, len(nodes)):
            node = nodes[i]
            if append_nodes:
                results.append(node)
            for param in param_names:
                if node in params[param]:
                    results[i + 1] += "\t" + str(params[param][node])
                else:
                    results[i + 1] += "\t" + "null"
        append_nodes = False
    results_file = open(filename, "w")
    for line in results:
        results_file.write(line + "\n")
    results_file.close()
    print("done.")
