"""Views of stored RINs: an interactive HTML page (Cytoscape.js), a text summary, a networkx graph.

Used by `structman show` (structman/scripts/showRin.py). Everything here works on storage
documents (storage.read_rin) and touches neither the config nor a database.

The page is one self-contained file: the network is embedded as JSON, Cytoscape.js and the fcose
layout come from public CDNs (the page says so when they cannot be loaded). Nodes are residues
(colour by chain or by a metric, shape by kind, size by a metric), edges are residue pairs per
interaction family with the atom class pairs merged; the class breakdown is shown on click.
"""

from __future__ import annotations

import glob
import json
import math
import os
import re

import numpy as np

from structman.lib.rin_engine import constants as C
from structman.lib.rin_engine import storage

# Above this many residues the page starts with no chain selected.
LARGE_NETWORK_NODES = 5000
# Above this many residues of kept chains a page would be too large for a browser (8OTZ: 184k
# residues, 1.9M residue pairs, 145 MB); show asks for chains instead.
MAX_PAGE_NODES = 20000

_QUERY_RE = re.compile(r'^(?P<pdb>[0-9][A-Za-z0-9]{3}|[Pp][Dd][Bb]_[0-9A-Za-z]{8})(-(?P<nr>[1-9][0-9]*))?(:(?P<chain>.+))?$')


# ---------------------------------------------------------------------------------------------
# Finding RIN files

def parse_query(query: str) -> tuple[str | None, int | None, str | None]:
    """'1abc' -> ('1ABC', None, None), '1ABC-2:A' -> ('1ABC', 2, 'A'); (None, None, None) otherwise."""
    match = _QUERY_RE.match(query.strip())
    if match is None:
        return None, None, None
    return match.group('pdb').upper(), int(match.group('nr')) if match.group('nr') else None, match.group('chain')


def _assembly_number(path: str) -> int:
    name = os.path.basename(path)[:-len(storage.SUFFIX)]
    try:
        return int(name.rsplit('-', 1)[1])
    except (IndexError, ValueError):
        return 0


def find_rin_files(db_folder: str | None, query: str) -> tuple[list[tuple[str, str]], str | None]:
    """([(structure id, path)], chain) for a PDB id, an assembly id, either with ':chain', or a path."""
    if query.endswith(storage.SUFFIX) and os.path.isfile(query):
        return [(os.path.basename(query)[:-len(storage.SUFFIX)], os.path.abspath(query))], None
    pdb_id, assembly_nr, chain = parse_query(query)
    if pdb_id is None or not db_folder:
        return [], chain
    if assembly_nr is not None:
        structure_id = f'{pdb_id}-{assembly_nr}'
        path = storage.rin_file_path(db_folder, structure_id)
        return ([(structure_id, path)] if os.path.isfile(path) else []), chain
    folder = os.path.dirname(storage.rin_file_path(db_folder, f'{pdb_id}-1'))
    paths = sorted(glob.glob(os.path.join(folder, f'{glob.escape(pdb_id)}-*{storage.SUFFIX}')), key = _assembly_number)
    paths = [path for path in paths if _assembly_number(path) > 0]
    return [(os.path.basename(path)[:-len(storage.SUFFIX)], path) for path in paths], chain


# ---------------------------------------------------------------------------------------------
# The network as shown

def _clean(value) -> float | None:
    value = float(value)
    if not math.isfinite(value):
        return None
    return float(f'{value:.4g}')


def _residue_number(res_id, ins_code) -> str:
    ins = str(ins_code)
    return f'{int(res_id)}{ins if ins not in ("", "?", ".", "_") else ""}'


def _select_chains(kept_chains: list[str], wanted: list[str] | None) -> set[str] | None:
    """Exact chain names; a name that matches none also matches its kept symmetry copies ('A' -> 'A-2')."""
    if not wanted:
        return None
    selected = set()
    for name in wanted:
        if name in kept_chains:
            selected.add(name)
        else:
            selected.update(chain for chain in kept_chains if chain.startswith(f'{name}-'))
    return selected


def network_payload(document: dict, chains: list[str] | None = None, families: list[str] | None = None) -> dict:
    """What the page shows of one RIN.

    Core nodes are the residues of the kept chain instances (of the chosen chains); partner nodes
    are residues outside of them with an interaction to a core node. Edges are residue pairs per
    family, from the chosen families, with at least one core end.
    """
    nodes, edges = document['nodes'], document['edges']
    n_nodes = len(nodes['res_id'])
    family_names = list(document.get('families') or C.FAMILY_NAMES)
    class_names = list(document.get('classes') or C.CLASS_NAMES)
    kind_names = list(document.get('kinds') or C.KIND_NAMES)
    chain_of = np.asarray(nodes['chain']).astype(str)
    kept = np.asarray(nodes['kept']).astype(bool)

    kept_chains = list(dict.fromkeys(chain_of[kept].tolist()))
    selected = _select_chains(kept_chains, chains)
    if chains and not selected:
        raise ValueError(f'No kept chain named {", ".join(chains)}; kept chains: {", ".join(kept_chains) or "-"}')
    core = kept.copy() if selected is None else (kept & np.isin(chain_of, sorted(selected)))

    unknown = sorted(set(families or ()) - set(family_names))
    if unknown:
        raise ValueError(f'Unknown families {", ".join(unknown)}; known: {", ".join(family_names)}')
    family_codes = None if not families else np.array([family_names.index(name) for name in families])

    a = np.asarray(edges['a']).astype(np.int64)
    b = np.asarray(edges['b']).astype(np.int64)
    family = np.asarray(edges['family']).astype(np.int64)
    use = core[a] | core[b] if len(a) else np.zeros(0, dtype = bool)
    if family_codes is not None and len(a):
        use &= np.isin(family, family_codes)
    order = np.flatnonzero(use)

    shown = core.copy()
    shown[a[order]] = True
    shown[b[order]] = True
    node_index = np.full(n_nodes, -1, dtype = np.int64)
    shown_ids = np.flatnonzero(shown)
    node_index[shown_ids] = np.arange(len(shown_ids))

    metric_names = list(document.get('metrics', {}).get('names') or [])
    metric_values = np.asarray(document.get('metrics', {}).get('values', np.zeros((n_nodes, 0))), dtype = np.float64)
    if metric_values.shape[0] != n_nodes:
        metric_values = np.full((n_nodes, len(metric_names)), np.nan)
    extra = [('SAS area', nodes.get('sas_area')), ('B-factor / pLDDT', nodes.get('b_factor'))]
    for name, column in extra:
        if column is not None and len(column) == n_nodes:
            metric_names.append(name)
            metric_values = np.concatenate([metric_values, np.asarray(column, dtype = np.float64)[:, None]], axis = 1)

    node_rows = []
    for i in shown_ids:
        node_rows.append([
            str(chain_of[i]), _residue_number(nodes['res_id'][i], nodes['ins_code'][i]), str(nodes['res_name'][i]),
            int(nodes['kind'][i]), int(kept[i]), int(core[i]), [_clean(v) for v in metric_values[i]],
        ])

    # Residue pairs per family; the atom class pairs become the breakdown.
    edge_rows = []
    counts = {name: 0 for name in family_names}
    if len(order):
        ea, eb, ef = node_index[a[order]], node_index[b[order]], family[order]
        key = (ea * len(shown_ids) + eb) * len(family_names) + ef
        sort = np.argsort(key, kind = 'stable')
        key_sorted = key[sort]
        starts = np.flatnonzero(np.r_[True, key_sorted[1:] != key_sorted[:-1]])
        ends = np.r_[starts[1:], len(sort)]
        rows = order[sort]
        class_a = np.asarray(edges['class_a']).astype(np.int64)[rows]
        class_b = np.asarray(edges['class_b']).astype(np.int64)[rows]
        count = np.asarray(edges['count']).astype(np.int64)[rows]
        score = np.asarray(edges['score']).astype(np.float64)[rows]
        distance = np.asarray(edges['min_distance']).astype(np.float64)[rows]
        for start, end in zip(starts, ends):
            first = sort[start]
            fam = int(ef[first])
            breakdown = [[int(class_a[k]), int(class_b[k]), int(count[k]), _clean(score[k]), _clean(distance[k])] for k in range(start, end)]
            edge_rows.append([int(ea[first]), int(eb[first]), fam, int(count[start:end].sum()), _clean(score[start:end].sum()),
                              _clean(distance[start:end].min()), breakdown])
            counts[family_names[fam]] += 1

    stats = document.get('stats', {})
    return {
        'id': document.get('structure_id', ''),
        'engine': document.get('engine_version', ''),
        'mode': document.get('mode', {}),
        'source': {key: value for key, value in (document.get('source') or {}).items() if isinstance(value, (str, int, float))},
        'families': family_names,
        'classes': class_names,
        'kinds': kind_names,
        'metrics': metric_names,
        'chains': list(dict.fromkeys(chain_of[core].tolist())),
        'all_chains': kept_chains,
        'nodes': node_rows,
        'edges': edge_rows,
        'counts': counts,
        'filter': {'chains': sorted(selected) if selected else None, 'families': list(families) if families else None},
        'totals': {'nodes': int(stats.get('nodes', n_nodes)), 'kept_nodes': int(stats.get('kept_nodes', int(kept.sum()))),
                   'edges': int(stats.get('edges', len(a)))},
    }


def kept_chain_sizes(document: dict) -> dict[str, int]:
    """Residues per kept chain instance, largest first."""
    chain_of = np.asarray(document['nodes']['chain']).astype(str)
    kept = np.asarray(document['nodes']['kept']).astype(bool)
    names, counts = np.unique(chain_of[kept], return_counts = True)
    order = np.argsort(-counts, kind = 'stable')
    return {str(names[i]): int(counts[i]) for i in order}


def core_size(document: dict, chains: list[str] | None = None) -> int:
    sizes = kept_chain_sizes(document)
    selected = _select_chains(list(sizes), chains)
    return sum(n for chain, n in sizes.items() if selected is None or chain in selected)


def _chain_list(chains: list[str], limit: int = 16) -> str:
    if len(chains) <= limit:
        return ', '.join(chains) or '-'
    return f'{", ".join(chains[:limit])}, ... ({len(chains)} chains)'


# ---------------------------------------------------------------------------------------------
# Text summary

def summary_text(payload: dict, top: int = 10) -> str:
    lines = [f'{payload["id"]}  (engine {payload["engine"]})']
    core = sum(1 for node in payload['nodes'] if node[5])
    lines.append(f'  residues: {core} in kept chain instances ({_chain_list(payload["chains"])}), '
                 f'{len(payload["nodes"]) - core} partners outside')
    counts = ', '.join(f'{name} {n}' for name, n in payload['counts'].items() if n)
    lines.append(f'  residue pairs: {sum(payload["counts"].values())} ({counts or "none"})')
    mode = payload.get('mode') or {}
    flags = [f'{key}={value}' for key, value in mode.items()]
    if flags:
        lines.append(f'  mode: {", ".join(flags)}')
    if 'AbsoluteCentrality' in payload['metrics']:
        column = payload['metrics'].index('AbsoluteCentrality')
        ranked = sorted((node for node in payload['nodes'] if node[5] and node[6][column] is not None),
                        key = lambda node: -node[6][column])[:top]
        if ranked:
            lines.append('  highest betweenness: ' + ', '.join(f'{n[0]}:{n[1]} {n[2]} ({n[6][column]:g})' for n in ranked))
    return '\n'.join(lines)


# ---------------------------------------------------------------------------------------------
# networkx (Cytoscape desktop: cyjs, graphml)

def to_networkx(payload: dict):
    import networkx as nx

    graph = nx.MultiGraph(name = payload['id'])
    for i, (chain, number, res_name, kind, kept, core, metrics) in enumerate(payload['nodes']):
        attributes = {'name': f'{chain}:{number}:{res_name}', 'chain': chain, 'residue': number, 'res_name': res_name,
                      'kind': payload['kinds'][kind], 'kept': bool(kept), 'core': bool(core)}
        attributes.update({name: value for name, value in zip(payload['metrics'], metrics) if value is not None})
        graph.add_node(i, **attributes)
    for a, b, family, count, score, distance, breakdown in payload['edges']:
        attributes = {'family': payload['families'][family], 'count': count,
                      'classes': ','.join(f'{payload["classes"][ca]}_{payload["classes"][cb]}' for ca, cb, *_ in breakdown)}
        if score is not None:
            attributes['score'] = score
        if distance is not None:
            attributes['min_distance'] = distance
        graph.add_edge(a, b, key = attributes['family'], **attributes)
    return graph


def write_export(payload: dict, path: str, kind: str) -> None:
    import networkx as nx

    graph = to_networkx(payload)
    if kind == 'graphml':
        nx.write_graphml(graph, path)
    elif kind == 'cyjs':
        with open(path, 'w') as f:
            json.dump(nx.cytoscape_data(graph), f)
    else:
        raise ValueError(f'Unknown export format {kind}')


# ---------------------------------------------------------------------------------------------
# HTML

def render_html(payloads: list[dict], title: str) -> str:
    data = json.dumps({'title': title, 'large': LARGE_NETWORK_NODES, 'assemblies': payloads}, allow_nan = False, separators = (',', ':'))
    data = data.replace('</', '<\\/')
    return _TEMPLATE.replace('__TITLE__', _escape(title)).replace('__DATA__', data)


def _escape(text: str) -> str:
    return text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')


_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ RIN</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.30.2/cytoscape.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/layout-base@2.0.1/layout-base.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cose-base@2.2.0/cose-base.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cytoscape-fcose@2.2.0/cytoscape-fcose.js"></script>
<style>
:root {
  --bg: #f7f7f5; --panel: #ffffff; --text: #1d1d1f; --muted: #6b6b70; --border: #dedee2;
  --accent: #2b6cb0; --canvas: #ffffff; --chip: #eef0f3; --label: #333333;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #17181b; --panel: #202226; --text: #e8e8ea; --muted: #9a9ba1; --border: #34363c;
    --accent: #6aa5e8; --canvas: #1b1c20; --chip: #2a2c31; --label: #d8d8dc;
  }
}
:root[data-theme="dark"] {
  --bg: #17181b; --panel: #202226; --text: #e8e8ea; --muted: #9a9ba1; --border: #34363c;
  --accent: #6aa5e8; --canvas: #1b1c20; --chip: #2a2c31; --label: #d8d8dc;
}
* { box-sizing: border-box; }
html, body { margin: 0; height: 100%; background: var(--bg); color: var(--text);
  font: 13px/1.4 system-ui, -apple-system, "Segoe UI", sans-serif; }
header { display: flex; flex-wrap: wrap; align-items: center; gap: 8px 16px; padding: 10px 16px;
  border-bottom: 1px solid var(--border); background: var(--panel); }
header h1 { font-size: 16px; margin: 0; font-weight: 600; }
.tabs { display: flex; flex-wrap: wrap; gap: 4px; }
.tabs button, .btn { border: 1px solid var(--border); background: var(--chip); color: var(--text);
  border-radius: 6px; padding: 3px 9px; cursor: pointer; font: inherit; }
.tabs button.active { background: var(--accent); border-color: var(--accent); color: #fff; }
.meta { color: var(--muted); font-size: 12px; }
main { display: grid; grid-template-columns: 250px 1fr 300px; height: calc(100% - 52px); min-height: 0; }
aside { overflow-y: auto; padding: 12px 14px; background: var(--panel); min-height: 0; }
aside.left { border-right: 1px solid var(--border); }
aside.right { border-left: 1px solid var(--border); }
h2 { font-size: 12px; text-transform: uppercase; letter-spacing: .04em; color: var(--muted); margin: 14px 0 6px; }
h2:first-child { margin-top: 0; }
label.row { display: flex; align-items: center; gap: 6px; padding: 1px 0; cursor: pointer; }
label.row .count { margin-left: auto; color: var(--muted); font-variant-numeric: tabular-nums; }
.swatch { width: 12px; height: 12px; border-radius: 3px; flex: none; }
select, input[type=text] { width: 100%; font: inherit; padding: 4px 6px; border: 1px solid var(--border);
  border-radius: 6px; background: var(--bg); color: var(--text); }
input[type=range] { width: 100%; }
.buttons { display: flex; gap: 4px; margin: 4px 0; flex-wrap: wrap; }
#cy { width: 100%; height: 100%; background: var(--canvas); min-height: 300px; }
.center { position: relative; min-height: 0; }
#notice { position: absolute; top: 12px; left: 12px; right: 12px; padding: 8px 12px; border-radius: 6px;
  background: var(--chip); border: 1px solid var(--border); display: none; z-index: 2; }
table { border-collapse: collapse; width: 100%; font-size: 12px; table-layout: fixed; }
td { padding: 2px 4px; border-bottom: 1px solid var(--border); vertical-align: top; overflow-wrap: anywhere; }
td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; width: 38%; }
.partner { color: var(--accent); cursor: pointer; }
.empty { color: var(--muted); }
@media (max-width: 900px) {
  main { grid-template-columns: 1fr; grid-template-rows: auto 60vh auto; height: auto; }
  aside.left, aside.right { border: none; border-bottom: 1px solid var(--border); }
}
</style>
</head>
<body>
<header>
  <h1 id="title"></h1>
  <div class="tabs" id="tabs"></div>
  <div class="meta" id="meta"></div>
  <div style="margin-left:auto" class="buttons">
    <button class="btn" id="fit">Fit</button>
    <button class="btn" id="png">PNG</button>
    <button class="btn" id="theme">Theme</button>
  </div>
</header>
<main>
  <aside class="left">
    <h2>Find residue</h2>
    <input type="text" id="search" placeholder="A:123, 123 or LYS">
    <h2>Chains</h2>
    <div class="buttons"><button class="btn" id="chains-all">all</button><button class="btn" id="chains-none">none</button></div>
    <div id="chains"></div>
    <label class="row"><input type="checkbox" id="others" checked> partners outside the selection</label>
    <h2>Interactions</h2>
    <div id="families"></div>
    <div id="area-box"><div class="meta">minimum contact area: <span id="area-value">0</span> &#8491;&sup2;</div>
      <input type="range" id="area" min="0" max="40" step="0.5" value="5"></div>
    <h2>Node size</h2>
    <select id="size"></select>
    <h2>Node colour</h2>
    <select id="color"></select>
    <h2>Layout</h2>
    <select id="layout">
      <option value="fcose">force directed</option>
      <option value="sequence">sequence (chains as rows)</option>
      <option value="concentric">concentric by node size metric</option>
      <option value="circle">circle</option>
    </select>
    <div class="buttons"><button class="btn" id="relayout">Run layout again</button></div>
    <h2>Shapes</h2>
    <div class="meta">circle protein &middot; rectangle nucleotide &middot; diamond ligand &middot;
      triangle metal / ion &middot; hexagon sugar &middot; faded: outside the kept chain instances</div>
  </aside>
  <div class="center"><div id="notice"></div><div id="cy"></div></div>
  <aside class="right"><h2>Details</h2><div id="details" class="empty">Click a residue or an interaction.</div></aside>
</main>
<script type="application/json" id="rin-data">__DATA__</script>
<script>
(function () {
const DATA = JSON.parse(document.getElementById('rin-data').textContent);
const FAMILY_COLORS = {
  contact: '#a0a4ab', hbond: '#1f77b4', weak_hbond: '#82b4e0', polar: '#17becf', weak_polar: '#93d8e0',
  ionic: '#d62728', pipi_parallel: '#9467bd', pipi_t: '#c3a8de', pication: '#e377c2', pihydrogen: '#f2a8cf',
  disulfide: '#bcbd22', metal: '#8c564b', halogen: '#2ca02c', hydrophobic: '#ff7f0e'
};
const CHAIN_COLORS = ['#4e79a7', '#f28e2b', '#59a14f', '#e15759', '#76b7b2', '#edc948', '#b07aa1', '#ff9da7',
  '#9c755f', '#bab0ac', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
  '#7f7f7f', '#bcbd22', '#17becf'];
const KIND_SHAPES = {protein: 'ellipse', nucleotide: 'round-rectangle', ligand: 'diamond', metal: 'triangle',
  ion: 'triangle', sugar: 'hexagon', other: 'star'};
const KIND_COLORS = {protein: '#4e79a7', nucleotide: '#59a14f', ligand: '#f28e2b', metal: '#8c564b',
  ion: '#b07aa1', sugar: '#edc948', other: '#7f7f7f'};
const $ = (id) => document.getElementById(id);
const state = {tab: 0, chains: new Set(), families: new Set(), minArea: 5, others: true,
  size: null, color: 'chain', layout: 'fcose', positions: {}};
let cy = null;
let adjacency = null;

function notice(text) { const box = $('notice'); box.textContent = text || ''; box.style.display = text ? 'block' : 'none'; }
function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
function assembly() { return DATA.assemblies[state.tab]; }
function nodeLabel(n) { return n[0] + ':' + n[1] + ' ' + n[2]; }
function fmt(v) { return v === null || v === undefined ? '–' : (Math.abs(v) >= 1000 || Number.isInteger(v) ? String(v) : v.toPrecision(3)); }
function escapeHtml(t) { return String(t).replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c])); }

function heat(t) {
  const stops = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
  t = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(Math.floor(t), stops.length - 2), f = t - i;
  const c = stops[i].map((v, k) => Math.round(v + f * (stops[i + 1][k] - v)));
  return 'rgb(' + c.join(',') + ')';
}

function setupTabs() {
  $('title').textContent = DATA.title;
  const tabs = $('tabs');
  tabs.innerHTML = '';
  if (DATA.assemblies.length < 2) return;
  DATA.assemblies.forEach((a, i) => {
    const b = document.createElement('button');
    b.textContent = a.id;
    b.className = i === state.tab ? 'active' : '';
    b.onclick = () => { state.tab = i; setupTabs(); setupControls(); build(true); };
    tabs.appendChild(b);
  });
}

function setupControls() {
  const a = assembly();
  adjacency = null;
  const core = a.nodes.filter(n => n[5]).length;
  $('meta').textContent = core + ' residues in kept chains, ' + (a.nodes.length - core) + ' partners · ' +
    a.edges.length + ' residue pairs · engine ' + a.engine + (a.mode && a.mode.region_restricted ? ' · region restricted' : '');
  const large = core > DATA.large;
  state.chains = new Set(large ? [] : a.chains);
  const chainBox = $('chains');
  chainBox.innerHTML = '';
  a.chains.forEach((c, i) => {
    const n = a.nodes.filter(x => x[5] && x[0] === c).length;
    const row = document.createElement('label');
    row.className = 'row';
    row.innerHTML = '<input type="checkbox"' + (state.chains.has(c) ? ' checked' : '') + '><span class="swatch" style="background:' +
      CHAIN_COLORS[i % CHAIN_COLORS.length] + '"></span>' + escapeHtml(c) + '<span class="count">' + n + '</span>';
    row.querySelector('input').onchange = (e) => { e.target.checked ? state.chains.add(c) : state.chains.delete(c); build(false); };
    chainBox.appendChild(row);
  });
  const famBox = $('families');
  famBox.innerHTML = '';
  if (state.families.size === 0) a.families.forEach(f => { if (f !== 'contact') state.families.add(f); });
  a.families.forEach(f => {
    if (!a.counts[f]) return;
    const row = document.createElement('label');
    row.className = 'row';
    row.innerHTML = '<input type="checkbox"' + (state.families.has(f) ? ' checked' : '') + '><span class="swatch" style="background:' +
      (FAMILY_COLORS[f] || '#888') + '"></span>' + f.replace('_', ' ') + '<span class="count">' + a.counts[f] + '</span>';
    row.querySelector('input').onchange = (e) => { e.target.checked ? state.families.add(f) : state.families.delete(f); build(false); };
    famBox.appendChild(row);
  });
  $('area-box').style.display = a.counts.contact ? 'block' : 'none';
  const size = $('size');
  size.innerHTML = '<option value="">uniform</option>' + a.metrics.map(m => '<option>' + escapeHtml(m) + '</option>').join('');
  if (state.size === null) state.size = a.metrics.includes('AbsoluteCentrality') ? 'AbsoluteCentrality' : '';
  size.value = a.metrics.includes(state.size) ? state.size : '';
  const color = $('color');
  color.innerHTML = '<option value="chain">chain</option><option value="kind">kind</option>' +
    a.metrics.map(m => '<option value="m:' + escapeHtml(m) + '">' + escapeHtml(m) + '</option>').join('');
  color.value = state.color.startsWith('m:') && !a.metrics.includes(state.color.slice(2)) ? 'chain' : state.color;
  state.color = color.value;
  notice(large ? 'Large network (' + core + ' residues): pick chains on the left to show them.' : '');
}

function visibleElements() {
  const a = assembly();
  const chainIndex = {};
  a.chains.forEach((c, i) => chainIndex[c] = i);
  const coreVisible = a.nodes.map(n => n[5] && state.chains.has(n[0]));
  const nodeVisible = coreVisible.slice();
  const edges = [];
  a.edges.forEach((e, i) => {
    const fam = a.families[e[2]];
    if (!state.families.has(fam)) return;
    if (fam === 'contact' && (e[4] === null || e[4] < state.minArea)) return;
    const ca = coreVisible[e[0]], cb = coreVisible[e[1]];
    if (!(ca && cb) && !(state.others && (ca || cb))) return;
    nodeVisible[e[0]] = nodeVisible[e[1]] = true;
    edges.push(i);
  });
  return {nodes: nodeVisible, core: coreVisible, edges: edges, chainIndex: chainIndex};
}

function metricColumn(name) { return assembly().metrics.indexOf(name); }

function nodeStyleData(n, range, colorRange, isCore) {
  const a = assembly();
  let size = 16;
  if (state.size) {
    const v = n[6][metricColumn(state.size)];
    size = v === null ? 10 : 10 + 30 * (range[1] > range[0] ? (v - range[0]) / (range[1] - range[0]) : 0.5);
  }
  let color;
  const kind = a.kinds[n[3]];
  if (state.color === 'chain') {
    const i = a.chains.indexOf(n[0]);
    color = i < 0 ? '#9a9a9a' : CHAIN_COLORS[i % CHAIN_COLORS.length];
  } else if (state.color === 'kind') {
    color = KIND_COLORS[kind] || '#888';
  } else {
    const v = n[6][metricColumn(state.color.slice(2))];
    color = v === null ? '#9a9a9a' : heat(colorRange[1] > colorRange[0] ? (v - colorRange[0]) / (colorRange[1] - colorRange[0]) : 0.5);
  }
  return {size: (kind === 'ion' ? 0.7 : 1) * size, color: color, shape: KIND_SHAPES[kind] || 'ellipse', op: isCore ? 1 : 0.35};
}

function range(name, visible) {
  if (!name) return [0, 1];
  const col = metricColumn(name);
  let lo = Infinity, hi = -Infinity;
  assembly().nodes.forEach((n, i) => { if (visible[i] && n[6][col] !== null) { lo = Math.min(lo, n[6][col]); hi = Math.max(hi, n[6][col]); } });
  return lo === Infinity ? [0, 1] : [lo, hi];
}

function build(relayout) {
  if (!cy) return;
  const a = assembly();
  const vis = visibleElements();
  const sizeRange = range(state.size, vis.nodes);
  const colorRange = range(state.color.startsWith('m:') ? state.color.slice(2) : '', vis.nodes);
  const positions = state.positions[state.tab] || (state.positions[state.tab] = {});
  const elements = [];
  let missing = 0;
  a.nodes.forEach((n, i) => {
    if (!vis.nodes[i]) return;
    const d = nodeStyleData(n, sizeRange, colorRange, vis.core[i]);
    d.id = 'n' + i; d.label = n[1] + ' ' + n[2]; d.idx = i;
    const el = {group: 'nodes', data: d};
    if (positions[i]) el.position = Object.assign({}, positions[i]); else missing++;
    elements.push(el);
  });
  vis.edges.forEach(i => {
    const e = a.edges[i];
    const fam = a.families[e[2]];
    const w = fam === 'contact' ? 0.5 + Math.min(4, Math.log1p(e[4] || 0) / 1.2) : 1.2 + Math.min(4, Math.log2(e[3]));
    elements.push({group: 'edges', data: {id: 'e' + i, source: 'n' + e[0], target: 'n' + e[1], color: FAMILY_COLORS[fam] || '#888', w: w, idx: i}});
  });
  cy.batch(() => {
    cy.elements().remove();
    cy.add(elements);
    cy.style().selector('edge').style('curve-style', vis.edges.length > 4000 ? 'haystack' : 'bezier').update();
  });
  if (relayout || missing > elements.filter(e => e.group === 'nodes').length / 3) {
    runLayout();
  } else if (missing) {
    // New nodes next to a placed neighbour, the rest stays where it was.
    cy.nodes().forEach(node => {
      if (positions[node.data('idx')]) return;
      const placed = node.neighborhood('node').filter(m => positions[m.data('idx')]);
      const p = placed.length ? positions[placed[0].data('idx')] : {x: 0, y: 0};
      node.position({x: p.x + (Math.random() - 0.5) * 60, y: p.y + (Math.random() - 0.5) * 60});
      positions[node.data('idx')] = Object.assign({}, node.position());
    });
  }
}

function sequencePositions() {
  const a = assembly();
  const rows = {};
  cy.nodes().forEach(node => {
    const n = a.nodes[node.data('idx')];
    (rows[n[0]] = rows[n[0]] || []).push(node);
  });
  const order = a.all_chains.filter(c => rows[c]).concat(Object.keys(rows).filter(c => !a.all_chains.includes(c)).sort());
  const pos = {};
  order.forEach((c, r) => {
    rows[c].sort((x, y) => {
      const nx = a.nodes[x.data('idx')], ny = a.nodes[y.data('idx')];
      const px = a.kinds[nx[3]] === 'protein' || a.kinds[nx[3]] === 'nucleotide' ? 0 : 1;
      const py = a.kinds[ny[3]] === 'protein' || a.kinds[ny[3]] === 'nucleotide' ? 0 : 1;
      return px - py || parseInt(nx[1]) - parseInt(ny[1]) || x.data('idx') - y.data('idx');
    });
    rows[c].forEach((node, k) => { pos[node.id()] = {x: k * 18, y: r * 160 + (k % 2) * 14}; });
  });
  return pos;
}

function runLayout() {
  if (!cy || cy.nodes().length === 0) return;
  const n = cy.nodes().length;
  let options;
  if (state.layout === 'sequence') {
    const pos = sequencePositions();
    options = {name: 'preset', positions: node => pos[node.id()]};
  } else if (state.layout === 'concentric') {
    options = {name: 'concentric', concentric: node => node.data('size'), levelWidth: () => 4, minNodeSpacing: 4, animate: false};
  } else if (state.layout === 'circle') {
    options = {name: 'circle', animate: false};
  } else {
    // 'draft' (spectral only) collapses large, loosely connected networks into lines.
    options = {name: 'fcose', quality: 'default', randomize: true, animate: false, packComponents: true,
      nodeRepulsion: 5000, idealEdgeLength: 45, numIter: n > 3000 ? 1000 : 2500};
  }
  const positions = state.positions[state.tab] = {};
  let layout;
  try { layout = cy.layout(options); } catch (err) { layout = cy.layout({name: 'cose', animate: false}); }
  layout.one('layoutstop', () => {
    cy.nodes().forEach(node => { positions[node.data('idx')] = Object.assign({}, node.position()); });
    cy.fit(undefined, 30);
  });
  layout.run();
}

function neighbours(i) {
  if (!adjacency) {
    adjacency = {};
    assembly().edges.forEach((e, k) => { (adjacency[e[0]] = adjacency[e[0]] || []).push(k); (adjacency[e[1]] = adjacency[e[1]] || []).push(k); });
  }
  return adjacency[i] || [];
}

function showNode(i) {
  const a = assembly();
  const n = a.nodes[i];
  let html = '<table>' +
    '<tr><td>residue</td><td class="num">' + escapeHtml(nodeLabel(n)) + '</td></tr>' +
    '<tr><td>kind</td><td class="num">' + a.kinds[n[3]] + '</td></tr>' +
    '<tr><td>instance</td><td class="num">' + (n[4] ? 'kept' : 'symmetry copy / other') + '</td></tr>';
  a.metrics.forEach((m, k) => { html += '<tr><td>' + escapeHtml(m) + '</td><td class="num">' + fmt(n[6][k]) + '</td></tr>'; });
  html += '</table>';
  const byFamily = {};
  neighbours(i).forEach(k => { const e = a.edges[k]; (byFamily[a.families[e[2]]] = byFamily[a.families[e[2]]] || []).push(e); });
  a.families.forEach(f => {
    const list = byFamily[f];
    if (!list) return;
    html += '<h2><span class="swatch" style="display:inline-block;background:' + (FAMILY_COLORS[f] || '#888') + '"></span> ' +
      f.replace('_', ' ') + ' (' + list.length + ')</h2><table>';
    list.sort((x, y) => (y[4] || 0) - (x[4] || 0)).forEach(e => {
      const other = e[0] === i ? e[1] : e[0];
      const classes = e[6].map(b => (e[0] === i ? a.classes[b[0]] + '-' + a.classes[b[1]] : a.classes[b[1]] + '-' + a.classes[b[0]])).join(', ');
      html += '<tr><td><span class="partner" data-idx="' + other + '">' + escapeHtml(nodeLabel(a.nodes[other])) + '</span><br><span class="meta">' +
        classes + '</span></td><td class="num">' + (f === 'contact' ? fmt(e[4]) + ' &#8491;&sup2;' : '&times;' + e[3]) +
        '<br><span class="meta">' + fmt(e[5]) + ' &#8491;</span></td></tr>';
    });
    html += '</table>';
  });
  if (!neighbours(i).length) html += '<p class="empty">No interactions in the selected families.</p>';
  $('details').className = '';
  $('details').innerHTML = html;
  $('details').querySelectorAll('.partner').forEach(el => el.onclick = () => focusNode(parseInt(el.dataset.idx)));
}

function showEdge(k) {
  const a = assembly();
  const e = a.edges[k];
  const f = a.families[e[2]];
  let html = '<table><tr><td>family</td><td class="num">' + f.replace('_', ' ') + '</td></tr>' +
    '<tr><td>residues</td><td class="num"><span class="partner" data-idx="' + e[0] + '">' + escapeHtml(nodeLabel(a.nodes[e[0]])) +
    '</span><br><span class="partner" data-idx="' + e[1] + '">' + escapeHtml(nodeLabel(a.nodes[e[1]])) + '</span></td></tr>' +
    '<tr><td>' + (f === 'contact' ? 'atom contacts' : 'atom pairs') + '</td><td class="num">' + e[3] + '</td></tr>' +
    '<tr><td>' + (f === 'contact' ? 'contact area (&#8491;&sup2;)' : 'score') + '</td><td class="num">' + fmt(e[4]) + '</td></tr>' +
    '<tr><td>shortest distance (&#8491;)</td><td class="num">' + fmt(e[5]) + '</td></tr></table>' +
    '<h2>By atom classes</h2><table>';
  e[6].forEach(b => { html += '<tr><td>' + a.classes[b[0]] + ' &ndash; ' + a.classes[b[1]] + '</td><td class="num">&times;' + b[2] +
    '</td><td class="num">' + fmt(b[3]) + '</td><td class="num">' + fmt(b[4]) + ' &#8491;</td></tr>'; });
  html += '</table><p class="meta">mc main chain, sc side chain, nbb / nba nucleotide backbone / base, lig ligand</p>';
  $('details').className = '';
  $('details').innerHTML = html;
  $('details').querySelectorAll('.partner').forEach(el => el.onclick = () => focusNode(parseInt(el.dataset.idx)));
}

function focusNode(i) {
  const node = cy.getElementById('n' + i);
  showNode(i);
  if (node.empty()) { notice(nodeLabel(assembly().nodes[i]) + ' is not shown with the current selection.'); return; }
  notice('');
  cy.elements().unselect();
  node.select();
  cy.animate({center: {eles: node}, zoom: Math.max(cy.zoom(), 1.5)}, {duration: 300});
}

function search(text) {
  text = text.trim().toUpperCase();
  if (!text) return;
  const a = assembly();
  const hits = [];
  a.nodes.forEach((n, i) => {
    const full = (n[0] + ':' + n[1]).toUpperCase();
    if (full === text || n[1].toUpperCase() === text || n[2].toUpperCase() === text) hits.push(i);
  });
  if (!hits.length) { notice('No residue ' + text + ' in this network.'); return; }
  if (hits.length === 1) { focusNode(hits[0]); return; }
  cy.elements().unselect();
  const shown = cy.collection(hits.map(i => cy.getElementById('n' + i)).filter(x => !x.empty()));
  shown.select();
  notice(hits.length + ' matches' + (shown.length < hits.length ? ', ' + shown.length + ' shown' : '') + '.');
  if (shown.length) cy.animate({fit: {eles: shown, padding: 60}}, {duration: 300});
}

function start() {
  setupTabs();
  if (typeof cytoscape === 'undefined') {
    setupControls();
    notice('Cytoscape.js could not be loaded from the CDN (no internet access?). The page needs it to draw the network.');
    return;
  }
  try { if (typeof cytoscapeFcose !== 'undefined') cytoscape.use(cytoscapeFcose); } catch (err) { /* already registered */ }
  cy = cytoscape({
    container: $('cy'), wheelSensitivity: 0.3, minZoom: 0.02, maxZoom: 6,
    style: [
      {selector: 'node', style: {'width': 'data(size)', 'height': 'data(size)', 'background-color': 'data(color)', 'shape': 'data(shape)',
        'opacity': 'data(op)', 'label': 'data(label)', 'font-size': 7, 'color': css('--label'), 'text-valign': 'center',
        'text-halign': 'center', 'min-zoomed-font-size': 7, 'border-width': 0}},
      {selector: 'node:selected', style: {'border-width': 3, 'border-color': css('--accent')}},
      {selector: 'edge', style: {'line-color': 'data(color)', 'width': 'data(w)', 'opacity': 0.75, 'curve-style': 'bezier'}},
      {selector: 'edge:selected', style: {'opacity': 1, 'width': 5}},
    ],
  });
  cy.on('tap', 'node', ev => showNode(ev.target.data('idx')));
  cy.on('tap', 'edge', ev => showEdge(ev.target.data('idx')));
  cy.on('dragfree', 'node', ev => { (state.positions[state.tab] || {})[ev.target.data('idx')] = Object.assign({}, ev.target.position()); });
  setupControls();
  build(true);
}

$('chains-all').onclick = () => { state.chains = new Set(assembly().chains); setupChainsChecked(); build(false); };
$('chains-none').onclick = () => { state.chains = new Set(); setupChainsChecked(); build(false); };
function setupChainsChecked() {
  $('chains').querySelectorAll('input').forEach((box, i) => { box.checked = state.chains.has(assembly().chains[i]); });
}
$('others').onchange = (e) => { state.others = e.target.checked; build(false); };
$('area').oninput = (e) => { state.minArea = parseFloat(e.target.value); $('area-value').textContent = e.target.value; };
$('area').onchange = () => build(false);
$('area-value').textContent = $('area').value;
$('size').onchange = (e) => { state.size = e.target.value; build(false); };
$('color').onchange = (e) => { state.color = e.target.value; build(false); };
$('layout').onchange = (e) => { state.layout = e.target.value; runLayout(); };
$('relayout').onclick = () => runLayout();
$('fit').onclick = () => cy && cy.fit(undefined, 30);
$('search').onkeydown = (e) => { if (e.key === 'Enter') search(e.target.value); };
$('png').onclick = () => {
  if (!cy) return;
  const link = document.createElement('a');
  link.href = cy.png({full: true, scale: 2, bg: css('--canvas')});
  link.download = assembly().id + '_rin.png';
  link.click();
};
$('theme').onclick = () => {
  const root = document.documentElement;
  const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
  root.dataset.theme = dark ? 'light' : 'dark';
  if (cy) cy.style().selector('node').style('color', css('--label')).update();
};
start();
})();
</script>
</body>
</html>
"""
