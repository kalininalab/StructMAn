import os.path

from . import nc_interactions, nci_residues, st_selection

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


def pdb_to_sif(sel_id, pdb_file_name, pdb_path, selection_lst, probe_filename, probe_path, sif_file, combine_int=True, directed=False):
    print("Running RINerator version 0.5")

    stsel_obj = st_selection.set_str_sel(sel_id, pdb_file_name, pdb_path, selection_lst)

    nci_obj = nc_interactions.get_ncint_probe(probe_filename, probe_path, stsel_obj, directed)

    nci_res = nci_residues.get_nci_res(nci_obj, stsel_obj, combine_int)

    # save sif file in the same path as pdb
    nci_graph = get_sif(os.path.join(probe_path, sif_file), nci_res, stsel_obj)

    return (stsel_obj, nci_obj, nci_res, nci_graph)


def get_sif(sif_file, nci_res, stsel_obj):
    nci_graph = NCIGraph(nci_res)
    nci_graph.get_edges(stsel_obj)
    file_names_tpl = nci_graph.write_sif(sif_file)
    print("INFO: Wrote SIF: %s and attribute files: %s %s" % (nci_graph.sif_file, nci_graph.ea_nr_int_file, nci_graph.ea_int_score_file))
    return nci_graph


class NCIGraph:
    def __init__(self, nci_res_obj):
        self.nci_res_obj = nci_res_obj
        self.edge_dic = {}  # edge_dic[edge_key] = edge_attr
        #                    edge_key = (node_label1,node_label2,edge_type)
        #                    node_label = chain:resseq:icode:resname    node_label=ligand_name
        #                    edge_type = int_type:nodesubtype1_nodesubtype2 nodesubtype= mc_mc or mc_sc or mc_ligand or sc_sc or ...
        #                    other edge types are possible
        self.edge_lst = []  # list of edge keys edge_lst=[edge_key1,edge_key2,...]
        self.not_connected_nodes_dic = {}  # dict of nodes that are not connected
        self.not_connected_nodes_lst = []  # list of nodes that are not connected
        self.node_lst = []  # list of nodes
        self.sif_file = None
        self.na_res_file = None
        self.ea_nr_int_file = None
        self.ea_int_score_file = None

    def get_edges(self, stsel_obj):
        """
        other related functions can be constructed where
        edge_type and edge_key can be changed to group interactions in different way
        or some interactions types can be filtered out
        """
        node_dic = {}
        for res_int_key in self.nci_res_obj.res_int_lst:
            res_int_attr = self.nci_res_obj.res_int_dic[res_int_key]
            res_id1 = res_int_attr.id1
            res_id2 = res_int_attr.id2
            if res_id1 == res_id2:
                # exclude edges with self
                continue
            int_type = res_int_attr.int_type
            int_subtype = res_int_attr.int_subtype
            nr_int = res_int_attr.nr_int
            int_score = res_int_attr.int_score
            node_label1 = self.get_node_label(res_id1, stsel_obj)
            node_label2 = self.get_node_label(res_id2, stsel_obj)
            edge_type = "%s:%s_%s" % (int_type, int_subtype[0], int_subtype[1])
            edge_key = (node_label1, node_label2, edge_type)
            edge_rkey = (node_label2, node_label1, edge_type)
            if edge_key in self.edge_dic:
                self.update_edge(edge_key, nr_int, int_score)
            elif edge_rkey in self.edge_dic:
                self.update_edge(edge_rkey, nr_int, int_score)
            else:
                self.new_edge(edge_key, edge_type, node_label1, node_label2, nr_int, int_score)
                self.edge_lst.append(edge_key)
                node_dic[node_label1] = True
                node_dic[node_label2] = True
        self.get_not_connected(node_dic, stsel_obj)

    def get_not_connected(self, node_dic, stsel_obj):
        for res_id in stsel_obj.res_id_lst:
            node_label = self.get_node_label(res_id, stsel_obj)
            self.node_lst.append(node_label)
            if node_label in node_dic:
                continue
            self.not_connected_nodes_dic[node_label] = True
            self.not_connected_nodes_lst.append(node_label)

    def write_sif(self, sif_file):
        self.sif_file = sif_file
        self.ea_nr_int_file = os.path.splitext(self.sif_file)[0] + "_nrint.ea"
        self.ea_int_score_file = os.path.splitext(self.sif_file)[0] + "_intsc.ea"
        self.na_res_file = os.path.splitext(self.sif_file)[0] + "_res.txt"
        sif_fobj = open(sif_file, 'w')
        nr_int_fobj = open(self.ea_nr_int_file, 'w')
        nr_int_fobj.write('EdgeID\tNrInteractions\n')
        int_score_fobj = open(self.ea_int_score_file, 'w')
        int_score_fobj.write('EdgeID\tInteractionScore\n')
        res_fobj = open(self.na_res_file, 'w')
        # res_fobj.write('NodeID\n')
        for edge_key in self.edge_lst:
            edge_attr = self.edge_dic[edge_key]
            node_label1 = edge_attr.node_label1
            node_label2 = edge_attr.node_label2
            edge_type = edge_attr.edge_type
            nr_int = edge_attr.nr_int
            int_score = edge_attr.int_score
            sif_fobj.write("%s %s %s\n" % (node_label1, edge_type, node_label2))
            nr_int_fobj.write("%s (%s) %s\t%d\n" % (node_label1, edge_type, node_label2, nr_int))
            int_score_fobj.write("%s (%s) %s\t%f\n" % (node_label1, edge_type, node_label2, int_score))
        for node_label in self.not_connected_nodes_lst:
            sif_fobj.write("%s\n" % node_label)
        for res_label in self.node_lst:
            res_fobj.write("%s\n" % res_label)
        int_score_fobj.close()
        nr_int_fobj.close()
        sif_fobj.close()

    def get_node_label(self, res_id, stsel_obj):
        if res_id in stsel_obj.res_id_dic:
            (comp_name, comp_type) = stsel_obj.res_id_dic[res_id]
        else:
            print("ERROR: No residue %s" % str(res_id))
            return None
        if comp_type != 'protein':
            # check if it works always?
            chain_id = stsel_obj.get_chain_from_res_id(res_id)
            if chain_id == ' ':
                chain_id = '_'
            res_seq = stsel_obj.get_res_seq_from_res_id(res_id)
            icode = stsel_obj.get_icode_from_res_id(res_id)
            if icode == ' ':
                icode = '_'
            node_label = "%s:%s:%s:%s" % (chain_id, str(res_seq), icode, comp_name)
        else:
            chain_id = stsel_obj.get_chain_from_res_id(res_id)
            if chain_id == ' ':
                chain_id = '_'
            res_seq = stsel_obj.get_res_seq_from_res_id(res_id)
            icode = stsel_obj.get_icode_from_res_id(res_id)
            if icode == ' ':
                icode = '_'
            resname = stsel_obj.get_resname_from_res_id(res_id)
            node_label = "%s:%s:%s:%s" % (chain_id, str(res_seq), icode, resname)
        return node_label

    def update_edge(self, edge_key, nr_int, int_score):
        edge_attr = self.edge_dic[edge_key]
        edge_attr.incr_nr_int(nr_int)
        edge_attr.incr_int_score(int_score)

    def new_edge(self, edge_key, edge_type, node_label1, node_label2, nr_int, int_score):
        edge_attr = EdgeAttr(node_label1, node_label2, edge_type, nr_int, int_score)
        self.edge_dic[edge_key] = edge_attr


class EdgeAttr:
    def __init__(self, node_label1, node_label2, edge_type, nr_int, int_score):
        self.node_label1 = node_label1
        self.node_label2 = node_label2
        self.edge_type = edge_type
        self.nr_int = nr_int
        self.int_score = int_score

    def incr_nr_int(self, nr_int):
        self.nr_int += nr_int

    def incr_int_score(self, int_score):
        self.int_score += int_score
