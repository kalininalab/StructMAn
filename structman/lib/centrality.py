#!/usr/bin/python3
import argparse
import gzip
import os
from collections import defaultdict

import igraph

from structman.lib.sdsc.consts import residues as residue_consts
from structman.lib.sdsc.consts import ligands as ligand_consts


class ChainGraph():
    def __init__(self):
        self.res_to_vertex = defaultdict(lambda: len(self.res_to_vertex))
        self.edges = []
        self.weights = []
        self.weights2 = []


def normalize_centrality(btw_cent, min_btw, max_btw, num_residues):
    # Betweenness centrality only makes sense if there is at least one vertex connected to two edges, otherwise everything will be 0
    if max_btw == 0:
        return (None, None, None)

    btw_cent_norm_len = btw_cent / ((num_residues - 1) * (num_residues - 2) / 2.0)

    if min_btw == max_btw:
        btw_cent_norm_max = 1.0
    else:
        btw_cent_norm_max = (btw_cent - min_btw) / (max_btw - min_btw)

    return (btw_cent, btw_cent_norm_len, btw_cent_norm_max)


def main(path, pdb_id=None):
    dir_name = os.path.dirname(path)
    if pdb_id is None:
        pdb_id = os.path.basename(path)[0:4]

    scores_path = dir_name + '/' + pdb_id + '_intsc.ea'
    if os.path.exists(scores_path):
        scores_file = open(scores_path, 'r')
    else:
        scores_path = scores_path + ".gz"
        scores_file = gzip.open(scores_path, 'rt')

    out_path = dir_name + '/' + pdb_id + '_btw_cent.txt.gz'
    out_file = gzip.open(out_path, 'wt')

    chains = defaultdict(ChainGraph)
    all_chains = ChainGraph()

    # Skip header
    next(scores_file)

    for line in scores_file:
        res1, int_type, res2, score = line.strip().split()

        if int_type != "(combi:all_all)":
            continue

        # Don't include clashing residues in our network
        score = float(score)
        score2 = score
        if score < 0:
            score = 0

        score2 = abs(score2)

        chain1 = res1.split(":")[0]
        chain2 = res2.split(":")[0]

        res_type1 = res1.split(":")[3]
        res_type2 = res2.split(":")[3]
        if res_type1 in ligand_consts.BORING_LIGANDS and res_type1 not in residue_consts.THREE_TO_ONE:
            continue
        if res_type2 in ligand_consts.BORING_LIGANDS and res_type2 not in residue_consts.THREE_TO_ONE:
            continue

        vertex_id1 = chains[chain1].res_to_vertex[res1]
        vertex_id2 = chains[chain2].res_to_vertex[res2]

        # Ignore contacts between chains unless we are calculating complex centralities
        if chain1 == chain2:
            chains[chain1].edges.append((vertex_id1, vertex_id2))
            chains[chain1].weights.append(float('inf') if score == 0 else 1.0 / score)
            chains[chain1].weights2.append(float('inf') if score2 == 0 else 1.0 / score2)

        complex_vertex_id1 = all_chains.res_to_vertex[res1]
        complex_vertex_id2 = all_chains.res_to_vertex[res2]

        all_chains.edges.append((complex_vertex_id1, complex_vertex_id2))
        all_chains.weights.append(float('inf') if score == 0 else 1.0 / score)
        all_chains.weights2.append(float('inf') if score2 == 0 else 1.0 / score2)

    scores_file.close()

    out_file.write('#Residue\t'
                   'AbsoluteCentrality\tLengthNormalizedCentrality\tMinMaxNormalizedCentrality\t'
                   'AbsoluteCentralityWithNegative\tLengthNormalizedCentralityWithNegative\tMinMaxNormalizedCentralityWithNegative\t'
                   'AbsoluteComplexCentrality\tLengthNormalizedComplexCentrality\tMinMaxNormalizedComplexCentrality\t'
                   'AbsoluteComplexCentralityWithNegative\tLengthNormalizedComplexCentralityWithNegative\tMinMaxNormalizedComplexCentralityWithNegative\n')

    num_complex_residues = len(all_chains.res_to_vertex)
    complex_graph = igraph.Graph(n=len(all_chains.res_to_vertex), edges=all_chains.edges, directed=False)
    if num_complex_residues < 2 or len(all_chains.edges) < 1:
        return

    try:
        complex_betweenness = complex_graph.betweenness(directed=False, weights=all_chains.weights)
        complex_betweenness_with_negative = complex_graph.betweenness(directed=False, weights=all_chains.weights2)
    except:
        print('Centrality error:', pdb_id, num_complex_residues, len(all_chains.edges))

    min_complex_btw_cent = min(complex_betweenness)
    max_complex_btw_cent = max(complex_betweenness)
    min_complex_btw_cent_with_negative = min(complex_betweenness_with_negative)
    max_complex_btw_cent_with_negative = max(complex_betweenness_with_negative)

    for chain, chain_graph in sorted(chains.items()):
        residues = [r for r, v in sorted(chain_graph.res_to_vertex.items(), key=lambda x: x[1])]
        num_residues = len(residues)
        if num_residues < 2 or len(chain_graph.edges) < 1:
            continue
        graph = igraph.Graph(n=num_residues, edges=chain_graph.edges, directed=False)

        try:
            betweenness = graph.betweenness(directed=False, weights=chain_graph.weights)
            betweenness_with_negative = graph.betweenness(directed=False, weights=chain_graph.weights2)
        except:
            print('Centrality error:', pdb_id, chain, num_residues, len(chain_graph.edges))

        min_btw_cent = min(betweenness)
        max_btw_cent = max(betweenness)
        min_btw_cent_with_negative = min(betweenness_with_negative)
        max_btw_cent_with_negative = max(betweenness_with_negative)

        res_names = []
        res_centralities = []

        for res, btw_cent, btw_cent_with_negative in sorted(zip(residues, betweenness, betweenness_with_negative), key=lambda x: x[1], reverse=True):
            centralities = []
            centralities.extend(normalize_centrality(btw_cent, min_btw_cent, max_btw_cent, num_residues))
            centralities.extend(normalize_centrality(btw_cent_with_negative, min_btw_cent_with_negative, max_btw_cent_with_negative, num_residues))

            complex_vertex_id = all_chains.res_to_vertex[res]
            centralities.extend(normalize_centrality(complex_betweenness[complex_vertex_id],
                                                     min_btw_cent,
                                                     max_btw_cent,
                                                     num_complex_residues))
            centralities.extend(normalize_centrality(complex_betweenness_with_negative[complex_vertex_id],
                                                     min_btw_cent_with_negative,
                                                     max_btw_cent_with_negative,
                                                     num_complex_residues))
            if any(centralities) is not None:
                res_names.append(res)
                res_centralities.append(centralities)

        for entries in zip(res_names, res_centralities):
            out_file.write(entries[0] + '\t' + '\t'.join(map(str, entries[1])) + '\n')

    out_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate betweenness centrality for a given RIN')
    parser.add_argument('in_file', metavar='input.sif', help='input .sif or .sif.gz file')
    args = parser.parse_args()

    main(args.in_file)
