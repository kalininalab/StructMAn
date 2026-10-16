"""Graph metrics of a residue network: betweenness centralities and in/out flow.

Both run on the contact graph - residue pairs weighted by their total contact area - with boring
ligands (crystallisation additives, as in the legacy RIN DB) left out. Betweenness uses the legacy
implementation in structman/lib/centrality.py; contact areas are never negative, so its
"WithNegative" variants would only repeat the plain ones and are not computed. Above
MAX_COMPLEX_BETWEENNESS_NODES residues the complex-wide betweenness is skipped (it is
O(nodes x edges)); per-chain betweenness is always computed.
"""

from __future__ import annotations

import numpy as np

from structman.lib.rin_engine.network import ResidueNetwork

# Complex-wide weighted betweenness took 113 s for the 19850 residues of 4v6x and grows with
# nodes x edges; this bound keeps the worst cases at a few minutes.
MAX_COMPLEX_BETWEENNESS_NODES = 25000

CENTRALITY_NAMES = (
    'AbsoluteCentrality', 'LengthNormalizedCentrality', 'MinMaxNormalizedCentrality',
    'AbsoluteComplexCentrality', 'LengthNormalizedComplexCentrality', 'MinMaxNormalizedComplexCentrality',
)
FLOW_NAMES = (
    'InflowValue', 'OutflowValue', 'Degree', 'WeightSum', 'WeightedInflowValue', 'WeightedOutflowValue',
    'InflowValueFirstPass', 'OutflowValueFirstPass', 'WeightedInflowValueFirstPass', 'WeightedOutflowValueFirstPass',
)


def _graph(network: ResidueNetwork):
    from structman.base_utils import ligand_consts, residue_consts
    a, b, area = network.contact_graph()
    boring = np.array([
        name in ligand_consts.BORING_LIGANDS and name not in residue_consts.THREE_TO_ONE for name in network.res_name
    ], dtype = bool)
    keep = ~boring[a] & ~boring[b] & (area > 0)
    return a[keep], b[keep], area[keep]


def betweenness(network: ResidueNetwork) -> tuple[np.ndarray, str]:
    """(n_nodes x 6 float32 with NaN for residues without values, mode)."""
    from structman.lib import centrality

    values = np.full((network.n_nodes, len(CENTRALITY_NAMES)), np.nan, dtype = np.float32)
    a, b, area = _graph(network)
    if len(a) == 0:
        return values, 'empty'
    n_graph_nodes = len(np.union1d(a, b))
    compute_complex = n_graph_nodes <= MAX_COMPLEX_BETWEENNESS_NODES
    chain = network.chain
    # Chain-level values only for the kept instances; the copies around them are cut to the
    # region near the kept ones, their chain graphs would be fragments.
    kept_chains = set(chain[network.kept].tolist())
    results = centrality.compute_betweenness(
        zip(a.tolist(), b.tolist(), area.tolist()), lambda node: chain[node],
        with_negative = False, compute_complex = compute_complex, report_chains = kept_chains,
    )
    for node, row in results:
        # row: chain (3), chain WithNegative (3), complex (3), complex WithNegative (3)
        values[node] = [np.nan if v is None else v for v in (row[0:3] + row[6:9])]
    return values, 'exact' if compute_complex else 'complex_skipped'


def flow(network: ResidueNetwork) -> tuple[np.ndarray, str]:
    """(n_nodes x 10 float32 with NaN, mode) from inflow_outflow_centrality, if installed."""
    values = np.full((network.n_nodes, len(FLOW_NAMES)), np.nan, dtype = np.float32)
    try:
        from inflow_outflow_centrality.utilities import flow_metric
    except ImportError:
        return values, 'unavailable'
    import pandas as pd

    a, b, area = _graph(network)
    if len(a) == 0:
        return values, 'empty'
    table = flow_metric(pd.DataFrame({'a': a, 'b': b, 'w': area.astype(float)}), None, 'both')
    nodes = table['Node'].to_numpy(dtype = np.int64)
    for column, name in enumerate(FLOW_NAMES):
        if name in table.columns:
            values[nodes, column] = table[name].to_numpy(dtype = np.float64)
    return values, 'computed'
