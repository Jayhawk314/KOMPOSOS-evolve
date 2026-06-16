#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Export a self-contained, zoomable mammal tree viewer with a CONVERGENCE OVERLAY.

SLICE 1 (tree)
--------------
A zoomable Order -> Family -> Genus -> Species browser over the REAL dated mammal
tips (PHYLACINE / Upham 2019) joined with PanTHERIA taxonomy + traits.

SLICE 2 (convergence overlay)  <-- new
------------------------------
For a handful of binary traits, this script runs the thesis convergence pipeline
on the real dated chronogram and bakes the verdict into the viewer:

  * Fitch parsimony -> ancestral state + number of INDEPENDENT origins of the
    derived state, and a per-tip "origin cluster" id (which independent gain a
    species descends from);
  * PTP randomization (is the homoplasy structured or random scatter?);
  * neutral-Mk excess test on the REAL dated branch lengths (does it exceed
    drift?);
  * a tier verdict (T0..T3), exactly as in `phenoscape_excess_dated.py`.

In the viewer, selecting a convergence trait:
  * colors each derived species by its independent origin,
  * heats internal nodes by how many independent origins they contain,
  * draws cross-lineage BRIDGES between frontier carriers of different origins
    (the "convergent pair" signature), and
  * stamps the TIER VERDICT badge (origins, PTP p, Mk p_excess).

PROVENANCE (honest)
-------------------
  * Every leaf is a species that is BOTH on the dated chronogram AND has a
    PanTHERIA row, so all tips are real dated taxa with real trait data.
  * Navigation is taxonomic (PanTHERIA Order/Family/Genus); the convergence
    analysis is run on the REAL phylogeny, then mapped onto those tips. So a
    trait that arose independently in distant orders literally shows up as the
    same colour family scattered across the taxonomy, bridged by arcs.

Run:  python evolve/viewer/export_mammal_tree.py
Out:  evolve/viewer/mammal_tree.html  (open in any browser)
"""

from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import quantiles

HERE = Path(__file__).resolve().parent
EVOLVE = HERE.parent
sys.path.insert(0, str(EVOLVE.parents[0] / "src"))
sys.path.insert(0, str(EVOLVE))

from ingest import pantheria                                                  # noqa: E402
from chronogram import (load_mammal_chronogram, all_tip_labels, prune_to,     # noqa: E402
                        normalize_unit_depth, branch_length_dict)
from polarity import analyze, randomization_test, count_origins, _leaves      # noqa: E402
from mk import excess_homoplasy_test                                          # noqa: E402

OUT = HERE / "mammal_tree.html"
MISSING = -999.0

COL_ORDER = "MSW05_Order"
COL_FAMILY = "MSW05_Family"
COL_GENUS = "MSW05_Genus"
COL_BODY = "5-1_AdultBodyMass_g"
COL_TEMP = "28-2_Temp_Mean_01degC"
COL_TROPHIC = "6-2_TrophicLevel"
COL_ACTIVITY = "1-1_ActivityCycle"

# Monte-Carlo budgets for the overlay (kept modest so the build stays snappy;
# bump for a final figure).
N_PERM = 499
N_SIM = 299


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f == MISSING else f


# ---------------------------------------------------------------------------
# Hierarchy (slice 1)
# ---------------------------------------------------------------------------
def build_hierarchy(df, chron):
    tree: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    species = set()
    for binom, row in df.iterrows():
        if not isinstance(binom, str):
            continue
        if binom.replace(" ", "_") not in chron:
            continue
        order = str(row.get(COL_ORDER) or "Incertae_sedis")
        family = str(row.get(COL_FAMILY) or "Incertae_sedis")
        genus = str(row.get(COL_GENUS) or binom.split(" ")[0])
        b = _num(row.get(COL_BODY))
        leaf = {
            "name": binom,
            "rank": "species",
            "body_log": round(math.log10(b), 3) if (b and b > 0) else None,
            "temp_c": (round(t / 10.0, 1) if (t := _num(row.get(COL_TEMP))) is not None else None),
            "trophic": (int(x) if (x := _num(row.get(COL_TROPHIC))) is not None else None),
            "activity": (int(x) if (x := _num(row.get(COL_ACTIVITY))) is not None else None),
        }
        tree[order][family][genus].append(leaf)
        species.add(binom)

    def genus_node(g, leaves):
        leaves.sort(key=lambda x: x["name"])
        return {"name": g, "rank": "genus", "children": leaves}

    def family_node(f, genera):
        return {"name": f, "rank": "family",
                "children": [genus_node(g, lv) for g, lv in sorted(genera.items())]}

    def order_node(o, families):
        return {"name": o, "rank": "order",
                "children": [family_node(f, gen) for f, gen in sorted(families.items())]}

    roots = [order_node(o, fam) for o, fam in sorted(tree.items())]
    root = {"name": "Mammalia", "rank": "root", "children": roots}
    meta = {
        "n_species": len(species),
        "n_orders": len(tree),
        "n_families": sum(len(f) for f in tree.values()),
        "n_genera": sum(len(g) for f in tree.values() for g in f.values()),
        "chronogram_tips": len(chron),
        "pantheria_rows": int(len(df)),
    }
    return root, meta, species


# ---------------------------------------------------------------------------
# Convergence overlay (slice 2)
# ---------------------------------------------------------------------------
def _label_origins(tree, derived: str):
    """Assign each derived tip the id of the independent gain it descends from.

    Walks the Fitch-reconstructed tree (node.state already set by analyze()).
    A new origin id is minted at every edge parent!=derived -> child==derived;
    descendants inherit it until the state reverts. Returns
    (tip_origin {tip_label: origin_id}, n_origins)."""
    counter = [0]
    tip_origin: dict = {}

    def walk(n, parent_state, cur):
        if n.state == derived and parent_state != derived:
            counter[0] += 1
            cur = counter[0]
        elif n.state != derived:
            cur = None
        if n.is_leaf and n.state == derived and cur is not None:
            tip_origin[n.label] = cur
        for c in n.children:
            walk(c, n.state, cur)

    walk(tree, None, None)
    return tip_origin, counter[0]


def _tier(origins, p_fewer, p_excess):
    if origins >= 2 and p_fewer < 0.05 and p_excess < 0.05:
        return 3, "excess beyond neutral drift"
    if origins >= 2 and p_fewer < 0.05:
        return 2, "structured homoplasy"
    if origins >= 2:
        return 1, "homoplasy present but labile"
    return 0, "single origin"


def _binary_states(df, species, value_fn):
    """value_fn(row) -> 'yes' | 'no' | None. Keyed by underscore label."""
    out = {}
    for binom, row in df.iterrows():
        if not isinstance(binom, str) or binom not in species:
            continue
        v = value_fn(row)
        if v is not None:
            out[binom.replace(" ", "_")] = v
    return out


def analyze_trait(full, df, species, key, label, derived_label, value_fn):
    tip = _binary_states(df, species, value_fn)
    keep = set(tip)
    if len(keep) < 8:
        print(f"  [{key}] too few scored tips ({len(keep)}) -- skipped")
        return None

    pruned = prune_to(full, keep)
    normalize_unit_depth(pruned)
    brlen = branch_length_dict(pruned)
    # restrict tip states to what survived pruning
    placed = {n.label: tip[n.label] for n in _leaves(pruned) if n.label in tip}

    if len(set(placed.values())) < 2 or sum(1 for v in placed.values() if v == "yes") < 2:
        print(f"  [{key}] not enough variation (need >=2 'yes') -- skipped")
        return None

    res = analyze(label, pruned, placed)
    derived = "yes" if res.root_state == "no" else ("no" if res.root_state == "yes" else
                                                     next(s for s in res.states if s != res.root_state))
    origins = count_origins(pruned, derived, None)
    tip_origin_us, n_org = _label_origins(pruned, derived)

    ptp = randomization_test(pruned, placed, n_perm=N_PERM, seed=0)
    mk = excess_homoplasy_test(pruned, placed, n_sim=N_SIM, seed=0, brlen=brlen)
    tn, verdict = _tier(origins, ptp["p_fewer"], mk["p_excess"])

    # map underscore tip keys -> viewer leaf names (space form)
    tip_origin = {k.replace("_", " "): v for k, v in tip_origin_us.items()}

    print(f"  [{key}] scored {len(placed)} | derived '{derived}' | origins {origins} "
          f"| PTP p_fewer {ptp['p_fewer']:.3f} | Mk p_excess {mk['p_excess']:.3f} -> T{tn}")

    return {
        "key": key, "label": label,
        "derived_label": derived_label if derived == "yes" else f"not-{derived_label}",
        "root_state": res.root_state, "derived": derived,
        "n_scored": len(placed), "n_derived": len(tip_origin),
        "origins": origins, "steps": res.steps, "ci": round(res.consistency_index, 3),
        "ptp_p_fewer": round(ptp["p_fewer"], 4),
        "mk": {"obs": mk["observed_steps"], "null": round(mk["null_mean"], 2),
               "p_excess": round(mk["p_excess"], 4), "q_hat": round(mk["q_hat"], 3)},
        "tier": tn, "tier_label": verdict,
        "tip_origin": tip_origin,
    }


# Landscape-defining mammalian convergences, coded by MSW05 family / order.
# Family membership is external taxonomy (independent of the trait being tested),
# so counting independent origins on the dated tree is a legitimate convergence
# test -- not fitting. Lists are the clean, unambiguous core of each syndrome.
MARINE_FAMS = {
    # cetaceans
    "Balaenidae", "Balaenopteridae", "Eschrichtiidae", "Neobalaenidae",
    "Delphinidae", "Monodontidae", "Phocoenidae", "Physeteridae", "Kogiidae",
    "Platanistidae", "Iniidae", "Pontoporiidae", "Lipotidae", "Ziphiidae",
    # sirenians
    "Dugongidae", "Trichechidae",
    # pinnipeds
    "Otariidae", "Phocidae", "Odobenidae",
}
GLIDER_FAMS = {  # gliding membranes (flight is coded separately via Chiroptera)
    "Cynocephalidae", "Anomaluridae", "Petauridae", "Acrobatidae", "Pseudocheiridae",
}
FOSSORIAL_FAMS = {
    "Talpidae", "Chrysochloridae", "Notoryctidae",           # "moles" x3 lineages
    "Bathyergidae", "Geomyidae", "Spalacidae", "Ctenomyidae",  # subterranean rodents
}
MYRMECOPHAGY_FAMS = {
    "Myrmecophagidae", "Cyclopedidae",   # anteaters
    "Manidae",                            # pangolins
    "Orycteropodidae",                    # aardvark
    "Tachyglossidae",                     # echidnas
    "Myrmecobiidae",                      # numbat
}
HOPPING_FAMS = {
    "Macropodidae", "Potoroidae",   # kangaroos / rat-kangaroos
    "Dipodidae",                    # jerboas
    "Pedetidae",                    # springhares
}


def _clade_fn(families=frozenset(), orders=frozenset()):
    def fn(r):
        fam = str(r.get(COL_FAMILY) or "")
        ordr = str(r.get(COL_ORDER) or "")
        return "yes" if (fam in families or ordr in orders) else "no"
    return fn


def build_convergence(full, df, species):
    print("Analyzing convergence on the real dated chronogram ...")
    body_vals = [math.log10(b) for binom in species
                 if (b := _num(df.at[binom, COL_BODY] if binom in df.index else None)) and b > 0]
    giant_cut = quantiles(body_vals, n=10)[-1] if len(body_vals) > 10 else max(body_vals)

    defs = [
        # --- landscape-defining body-plan / locomotion convergences ---
        ("marine", "marine body plan", "marine",
         _clade_fn(families=MARINE_FAMS)),
        ("aerial", "aerial locomotion (flight + gliding)", "aerial",
         _clade_fn(families=GLIDER_FAMS, orders={"Chiroptera"})),
        ("fossorial", "fossorial digging", "fossorial",
         _clade_fn(families=FOSSORIAL_FAMS)),
        ("myrmecophagy", "myrmecophagy (anteating)", "anteater",
         _clade_fn(families=MYRMECOPHAGY_FAMS)),
        ("hopping", "bipedal hopping", "hopper",
         _clade_fn(families=HOPPING_FAMS)),
        # --- secondary trait-coded overlays (life-history columns) ---
        ("carnivory", "carnivory", "carnivore",
         lambda r: ("yes" if (x := _num(r.get(COL_TROPHIC))) == 3 else
                    ("no" if x is not None else None))),
        ("diurnality", "diurnality", "diurnal",
         lambda r: ("yes" if (x := _num(r.get(COL_ACTIVITY))) == 3 else
                    ("no" if x is not None else None))),
        ("giant", "giant body (top 10%)", "giant",
         lambda r: ("yes" if (b := _num(r.get(COL_BODY))) and b > 0 and math.log10(b) >= giant_cut
                    else ("no" if (b := _num(r.get(COL_BODY))) and b > 0 else None))),
    ]
    out = {}
    for key, label, derived_label, fn in defs:
        r = analyze_trait(full, df, species, key, label, derived_label, fn)
        if r:
            out[key] = r
    return out


def main():
    print("Loading dated mammal chronogram (cached) ...")
    full = load_mammal_chronogram()
    chron = all_tip_labels(full)
    print(f"  chronogram tips: {len(chron)}")
    df = pantheria()
    print(f"  PanTHERIA rows: {len(df)}")

    root, meta, species = build_hierarchy(df, chron)
    print(f"  joined dated tips with traits: {len(species)}")

    conv = build_convergence(full, df, species)
    meta["conv_traits"] = list(conv.keys())

    html = (HTML_TEMPLATE
            .replace("__DATA__", json.dumps(root, separators=(",", ":")))
            .replace("__META__", json.dumps(meta))
            .replace("__CONV__", json.dumps(conv, separators=(",", ":"))))
    OUT.write_text(html, encoding="utf-8")
    kb = OUT.stat().st_size / 1024
    print(f"\nWrote {OUT}  ({kb:.0f} KB)")
    print(f"  {meta['n_species']} species | {meta['n_genera']} genera | "
          f"{meta['n_families']} families | {meta['n_orders']} orders")
    print(f"  convergence overlays: {', '.join(conv.keys()) or '(none)'}")
    print("  Open it in any browser.")


# ---------------------------------------------------------------------------
# Self-contained HTML template. __DATA__ / __META__ / __CONV__ replaced above.
# ---------------------------------------------------------------------------
HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Mammal Tree Explorer — convergence overlay</title>
<script src="https://cdn.jsdelivr.net/npm/d3@7"></script>
<style>
  :root { --bg:#0e1116; --panel:#161b22; --ink:#e6edf3; --muted:#8b949e; --line:#30363d; --accent:#58a6ff; }
  * { box-sizing: border-box; }
  html,body { margin:0; height:100%; background:var(--bg); color:var(--ink);
    font:13px/1.4 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif; }
  #bar { display:flex; gap:14px; align-items:center; flex-wrap:wrap;
    padding:8px 14px; background:var(--panel); border-bottom:1px solid var(--line); }
  #bar h1 { font-size:15px; margin:0 8px 0 0; font-weight:600; }
  #bar .stat { color:var(--muted); font-size:12px; }
  label { color:var(--muted); margin-right:4px; }
  select, input, button { background:#0d1117; color:var(--ink); border:1px solid var(--line);
    border-radius:6px; padding:5px 8px; font:inherit; }
  button { cursor:pointer; } button:hover { border-color:var(--accent); }
  #wrap { display:flex; height:calc(100% - 47px); }
  #stage { flex:1; overflow:hidden; position:relative; }
  #side { width:300px; border-left:1px solid var(--line); background:var(--panel);
    padding:14px; overflow:auto; }
  #side h2 { font-size:14px; margin:0 0 2px; }
  #side .rank { color:var(--muted); font-size:11px; text-transform:uppercase; letter-spacing:.06em; }
  #side table { width:100%; border-collapse:collapse; margin-top:10px; }
  #side td { padding:3px 0; border-bottom:1px solid var(--line); vertical-align:top; }
  #side td.k { color:var(--muted); width:46%; }
  #pic { margin-top:12px; min-height:20px; }
  #pic img { max-width:100%; background:#fff; border-radius:8px; padding:6px; }
  #pic .credit { color:var(--muted); font-size:10px; margin-top:4px; }
  .node circle { stroke:#0d1117; stroke-width:1.2px; cursor:pointer; }
  .node text { fill:var(--ink); font-size:11px; }
  .node text.count { fill:var(--muted); }
  .link { fill:none; stroke:var(--line); stroke-width:1.4px; }
  .bridge { fill:none; stroke:#f0883e; stroke-width:1.4px; opacity:.55; pointer-events:none; }
  .node.hit circle { stroke:var(--accent); stroke-width:2.5px; }
  #legend { position:absolute; bottom:12px; left:12px; background:rgba(22,27,34,.92);
    border:1px solid var(--line); border-radius:8px; padding:8px 10px; max-width:240px; }
  #legend .row { display:flex; align-items:center; gap:6px; margin:2px 0; }
  #legend .sw { width:12px; height:12px; border-radius:3px; flex:0 0 auto; }
  #legend .title { color:var(--muted); font-size:11px; margin-bottom:4px; }
  #legend .note { color:var(--muted); font-size:10.5px; margin-top:6px; line-height:1.35; }
  #verdict { position:absolute; top:12px; left:12px; background:rgba(22,27,34,.95);
    border:1px solid var(--line); border-radius:10px; padding:10px 12px; max-width:280px; display:none; }
  #verdict .v-title { font-weight:600; margin-bottom:2px; }
  #verdict .badge { display:inline-block; padding:2px 8px; border-radius:999px;
    font-weight:600; font-size:12px; margin-bottom:6px; }
  #verdict table { width:100%; border-collapse:collapse; }
  #verdict td { padding:2px 0; font-size:12px; } #verdict td.k { color:var(--muted); }
  #verdict .sub { color:var(--muted); font-size:10.5px; margin-top:6px; line-height:1.35; }
  .T0{background:#484f58;color:#fff;} .T1{background:#9e6a03;color:#fff;}
  .T2{background:#1f6feb;color:#fff;} .T3{background:#238636;color:#fff;}
  #tip { position:absolute; pointer-events:none; background:#000d; color:#fff;
    padding:4px 8px; border-radius:6px; font-size:12px; opacity:0; transition:opacity .1s; }
</style>
</head>
<body>
<div id="bar">
  <h1>🐾 Mammal Tree Explorer</h1>
  <span class="stat" id="meta"></span>
  <span style="flex:1"></span>
  <label for="trait">View</label>
  <select id="trait"></select>
  <input id="search" placeholder="find genus / species…" size="18">
  <button id="collapse">Collapse all</button>
</div>
<div id="wrap">
  <div id="stage">
    <svg id="svg"><g id="view"><g id="bridges"></g></g></svg>
    <div id="verdict"></div>
    <div id="legend"></div>
    <div id="tip"></div>
  </div>
  <div id="side">
    <div class="rank" id="s-rank">click a node</div>
    <h2 id="s-name">—</h2>
    <div id="s-body"></div>
    <div id="pic"></div>
  </div>
</div>
<script>
const DATA = __DATA__;
const META = __META__;
const CONV = __CONV__;

document.getElementById('meta').textContent =
  `${META.n_species.toLocaleString()} dated species · ${META.n_genera} genera · ` +
  `${META.n_families} families · ${META.n_orders} orders`;

const TROPHIC = {1:'herbivore', 2:'omnivore', 3:'carnivore'};
const ACTIVITY = {1:'nocturnal', 2:'crepuscular / mixed', 3:'diurnal'};

// build the View selector: base traits + convergence overlays
const sel = document.getElementById('trait');
[['order','Taxonomic order'],['body_log','Body mass (log g)'],['temp_c','Mean temperature (°C)'],
 ['trophic','Trophic level'],['activity','Activity cycle']].forEach(([v,t])=>{
  const o=document.createElement('option'); o.value=v; o.textContent=t; sel.appendChild(o);
});
if (Object.keys(CONV).length){
  const og=document.createElement('optgroup'); og.label='⚡ Convergence overlay';
  for (const k in CONV){ const o=document.createElement('option');
    o.value='conv:'+k; o.textContent='origins of '+CONV[k].label; og.appendChild(o); }
  sel.appendChild(og);
}

const svg = d3.select('#svg'), view = d3.select('#view'), bridges = d3.select('#bridges');
const stage = document.getElementById('stage');
const tip = d3.select('#tip');
let W = stage.clientWidth, H = stage.clientHeight;
svg.attr('width', W).attr('height', H);

const root = d3.hierarchy(DATA);
let i = 0;
root.x0 = H / 2; root.y0 = 0;

root.eachAfter(d => {
  const num = {};
  if (!d.children && !d._children) {
    for (const k of ['body_log','temp_c','trophic','activity'])
      if (d.data[k] != null) num[k] = {sum:d.data[k], n:1};
    d.agg = {num, count:1};
  } else {
    const kids = d.children || d._children || [];
    let count = 0;
    for (const k of ['body_log','temp_c','trophic','activity']) num[k] = {sum:0, n:0};
    for (const c of kids) { count += c.agg.count;
      for (const k in c.agg.num){ num[k].sum+=c.agg.num[k].sum; num[k].n+=c.agg.num[k].n; } }
    d.agg = {num, count};
  }
});
root.each(d => { let p=d; while(p && p.data.rank!=='order') p=p.parent; d.orderName = p?p.data.name:null; });

function aggValue(d, key){ const a=d.agg.num[key]; return (a&&a.n)?a.sum/a.n:null; }

// ---- base color scales -----------------------------------------------------
const orders = Array.from(new Set(root.leaves().map(l => l.orderName))).sort();
const orderColor = d3.scaleOrdinal(orders, orders.map((_,k)=>d3.interpolateRainbow(k/orders.length)));
const bodyExt = d3.extent(root.leaves(), l => l.data.body_log);
const tempExt = d3.extent(root.leaves(), l => l.data.temp_c);
const bodyColor = d3.scaleSequential(bodyExt, d3.interpolateViridis);
const tempColor = d3.scaleSequential([tempExt[1], tempExt[0]], d3.interpolateRdBu);
const trophicColor = d3.scaleOrdinal([1,2,3], ['#7ec97e','#e3c34d','#e36d6d']);
const activityColor = d3.scaleOrdinal([1,2,3], ['#3b4a8c','#9a7bd0','#f2c14e']);

let TRAIT = 'order';     // base trait key, or null when overlay on
let CONVMODE = null;     // convergence trait key, or null
let originColor = null;  // id -> color (set when overlay active)

function computeConv(){
  const to = CONV[CONVMODE].tip_origin;
  // manual post-order over the FULL tree (d3.eachAfter only follows .children,
  // which is null for collapsed nodes whose real kids live in ._children).
  (function rec(d){
    const kids = d.children || d._children;
    if (!kids || !kids.length){
      const o = to[d.data.name];
      d.convO = new Set(o!=null ? [o] : []);
      d.convId = (o!=null) ? o : null;
      return;
    }
    const s = new Set();
    kids.forEach(c => { rec(c); c.convO.forEach(x=>s.add(x)); });
    d.convO = s; d.convId = null;
  })(root);
  originColor = id => d3.interpolateRainbow(((id*0.61803)%1));  // golden-ratio hue spread
}

function colorOf(d){
  if (CONVMODE){
    if (!d.children && !d._children)
      return d.convId!=null ? originColor(d.convId) : '#2a2f37';
    const k = d.convO ? d.convO.size : 0;
    return k===0 ? '#22272e' : d3.interpolateYlOrRd(Math.min(1, 0.25 + k/8));
  }
  if (TRAIT === 'order') return d.orderName ? orderColor(d.orderName) : '#888';
  if (TRAIT === 'body_log'){ const v=aggValue(d,'body_log'); return v==null?'#444':bodyColor(v); }
  if (TRAIT === 'temp_c'){ const v=aggValue(d,'temp_c'); return v==null?'#444':tempColor(v); }
  if (TRAIT === 'trophic'){ const v=aggValue(d,'trophic'); return v==null?'#444':trophicColor(Math.round(v)); }
  if (TRAIT === 'activity'){ const v=aggValue(d,'activity'); return v==null?'#444':activityColor(Math.round(v)); }
  return '#888';
}

const tree = d3.tree().nodeSize([16, 200]);
const zoom = d3.zoom().scaleExtent([0.05, 4]).on('zoom', e => view.attr('transform', e.transform));
svg.call(zoom);

function collapse(d){ if(d.children){ d._children=d.children; d._children.forEach(collapse); d.children=null; } }
root.children && root.children.forEach(collapse);

function toggle(d){ if(d.children){ d._children=d.children; d.children=null; } else { d.children=d._children; d._children=null; } update(d); }

function update(source){
  tree(root);
  const nodes = root.descendants(), links = root.links();

  const node = view.selectAll('g.node').data(nodes, d => d.id || (d.id = ++i));
  const enter = node.enter().append('g').attr('class','node')
    .attr('transform', d => `translate(${source.y0},${source.x0})`)
    .on('click', (e,d) => { select(d); if (d.children || d._children) toggle(d); })
    .on('mouseover', (e,d) => tip.style('opacity',1).html(label(d)))
    .on('mousemove', e => tip.style('left',(e.offsetX+14)+'px').style('top',(e.offsetY+10)+'px'))
    .on('mouseout', () => tip.style('opacity',0));
  enter.append('circle').attr('r', d => radius(d));
  enter.append('text').attr('dy','0.32em')
    .attr('x', d => (d._children||d.children) ? -8 : 8)
    .attr('text-anchor', d => (d._children||d.children) ? 'end' : 'start')
    .text(d => d.data.name);
  enter.append('text').attr('class','count').attr('dy','0.32em').attr('x',8).attr('text-anchor','start');

  const merge = enter.merge(node);
  merge.transition().duration(300).attr('transform', d => `translate(${d.y},${d.x})`);
  merge.select('circle').attr('r', d=>radius(d)).attr('fill', colorOf)
    .attr('stroke', d => d._children ? 'var(--accent)' : '#0d1117');
  merge.select('text:not(.count)')
    .attr('x', d => (d._children||d.children) ? -8 : 8)
    .attr('text-anchor', d => (d._children||d.children) ? 'end' : 'start');
  merge.select('text.count').text(d => d._children ? '+'+d.agg.count : '');

  node.exit().transition().duration(300)
    .attr('transform', d => `translate(${source.y},${source.x})`).remove();

  const link = view.selectAll('path.link').data(links, d => d.target.id);
  link.enter().insert('path','#bridges').attr('class','link')
      .attr('d', elbow({source:{x:source.x0,y:source.y0}, target:{x:source.x0,y:source.y0}}))
    .merge(link).transition().duration(300).attr('d', d=>elbow(d));
  link.exit().transition().duration(300)
      .attr('d', elbow({source:{x:source.x,y:source.y}, target:{x:source.x,y:source.y}})).remove();

  nodes.forEach(d => { d.x0=d.x; d.y0=d.y; });
  drawBridges();
}

function radius(d){ return d.data.rank==='species' ? 3.5 : Math.min(9, 4+Math.log10(d.agg.count+1)*2.2); }
function elbow(d){ return `M${d.source.y},${d.source.x}C${(d.source.y+d.target.y)/2},${d.source.x} ${(d.source.y+d.target.y)/2},${d.target.x} ${d.target.y},${d.target.x}`; }
function label(d){
  if (d.data.rank==='species') return `<i>${d.data.name}</i>`;
  return `${d.data.name} <span style="color:#8b949e">(${d.data.rank}, ${d.agg.count} sp.)</span>`;
}

// ---- convergence bridges ---------------------------------------------------
function drawBridges(){
  bridges.selectAll('*').remove();
  if (!CONVMODE) return;
  const rendered = root.descendants();
  const carriers = rendered.filter(d => (d.data.rank==='species' || d._children)
                                        && d.convO && d.convO.size>0);
  let pairs = [];
  if (carriers.length <= 10){
    for (let a=0;a<carriers.length;a++) for (let b=a+1;b<carriers.length;b++){
      const u=carriers[a], v=carriers[b];
      if (new Set([...u.convO, ...v.convO]).size >= 2) pairs.push([u,v]);
    }
  } else {
    const s = carriers.slice().sort((p,q)=>p.x-q.x);
    for (let k=0;k+1<s.length;k++)
      if (new Set([...s[k].convO, ...s[k+1].convO]).size >= 2) pairs.push([s[k], s[k+1]]);
  }
  pairs = pairs.slice(0, 90);
  bridges.selectAll('path').data(pairs).enter().append('path')
    .attr('class','bridge').attr('d', p => arc(p[0], p[1]));
}
function arc(u,v){
  const x1=u.x,y1=u.y,x2=v.x,y2=v.y;
  const span=Math.abs(x1-x2);
  const cy=Math.min(y1,y2) - span*0.35 - 30;   // bow left, away from the right-growing tree
  const cx=(x1+x2)/2;
  return `M${y1},${x1} Q${cy},${cx} ${y2},${x2}`;
}

// ---- selection / detail panel ---------------------------------------------
function select(d){
  d3.selectAll('g.node').classed('hit', n => n===d);
  document.getElementById('s-rank').textContent = d.data.rank;
  document.getElementById('s-name').innerHTML =
    d.data.rank==='species' ? `<i>${d.data.name}</i>` : d.data.name;
  const rows = [];
  if (d.data.rank!=='species') rows.push(['species below', d.agg.count]);
  if (CONVMODE){
    if (d.data.rank==='species'){
      rows.push([CONV[CONVMODE].derived_label+'?', d.convId!=null ? 'yes (origin #'+d.convId+')' : 'no']);
    } else if (d.convO){
      rows.push(['independent origins here', d.convO.size]);
    }
  }
  const bv=aggValue(d,'body_log'), tv=aggValue(d,'temp_c'), tr=aggValue(d,'trophic'), ac=aggValue(d,'activity');
  if (d.data.rank==='species'){
    if (d.data.body_log!=null) rows.push(['body mass', (10**d.data.body_log).toLocaleString(undefined,{maximumFractionDigits:0})+' g']);
    if (d.data.temp_c!=null) rows.push(['mean temp', d.data.temp_c+' °C']);
    if (d.data.trophic!=null) rows.push(['trophic', TROPHIC[d.data.trophic]||d.data.trophic]);
    if (d.data.activity!=null) rows.push(['activity', ACTIVITY[d.data.activity]||d.data.activity]);
  } else {
    if (bv!=null) rows.push(['mean body', (10**bv).toLocaleString(undefined,{maximumFractionDigits:0})+' g']);
    if (tv!=null) rows.push(['mean temp', tv.toFixed(1)+' °C']);
    if (tr!=null) rows.push(['mean trophic', tr.toFixed(2)]);
    if (ac!=null) rows.push(['mean activity', ac.toFixed(2)]);
  }
  document.getElementById('s-body').innerHTML =
    '<table>' + rows.map(r=>`<tr><td class="k">${r[0]}</td><td>${r[1]}</td></tr>`).join('') + '</table>';
  loadPic(d.data.name);
}

function loadPic(name){
  const pic = document.getElementById('pic');
  pic.innerHTML = '<span style="color:#8b949e">searching PhyloPic…</span>';
  const q = encodeURIComponent(name);
  fetch(`https://api.phylopic.org/autocomplete?query=${q}`)
    .then(r=>r.json())
    .then(j => { const term=(j.matches&&j.matches[0])||name.split(' ')[0];
      return fetch(`https://api.phylopic.org/images?filter_name=${encodeURIComponent(term)}&page=0&embed_items=true`); })
    .then(r=>r.json())
    .then(j => { const items=j&&j._embedded&&j._embedded.items;
      const img=items&&items[0]&&items[0]._links&&(items[0]._links.rasterFiles||[]).slice(-1)[0];
      pic.innerHTML = (img&&img.href)
        ? `<img src="${img.href}" alt="${name}"><div class="credit">silhouette: PhyloPic.org</div>`
        : '<span style="color:#8b949e">no PhyloPic silhouette found</span>'; })
    .catch(()=> pic.innerHTML='<span style="color:#8b949e">PhyloPic unavailable (offline?)</span>');
}

// ---- verdict + legend ------------------------------------------------------
function drawVerdict(){
  const box = document.getElementById('verdict');
  if (!CONVMODE){ box.style.display='none'; return; }
  const c = CONV[CONVMODE];
  box.style.display='block';
  box.innerHTML =
    `<div class="v-title">Convergence verdict — ${c.label}</div>` +
    `<span class="badge T${c.tier}">TIER ${c.tier} · ${c.tier_label}</span>` +
    `<table>` +
    `<tr><td class="k">independent origins</td><td>${c.origins}</td></tr>` +
    `<tr><td class="k">derived species</td><td>${c.n_derived} / ${c.n_scored} scored</td></tr>` +
    `<tr><td class="k">parsimony steps · CI</td><td>${c.steps} · ${c.ci}</td></tr>` +
    `<tr><td class="k">PTP p(fewer)</td><td>${c.ptp_p_fewer}${c.ptp_p_fewer<0.05?' ✓ structured':''}</td></tr>` +
    `<tr><td class="k">Mk obs / null</td><td>${c.mk.obs} / ${c.mk.null}</td></tr>` +
    `<tr><td class="k">Mk p(excess)</td><td>${c.mk.p_excess}${c.mk.p_excess<0.05?' ✓ excess':''}</td></tr>` +
    `</table>` +
    `<div class="sub">Ancestral state: ${c.root_state}. Run on the real dated chronogram ` +
    `(Fitch origins, PTP randomization, neutral-Mk on dated branches). ` +
    `Tier ≥ 3 = excess beyond drift; lower tiers are downgraded, not discarded.</div>`;
}

function drawLegend(){
  const L = d3.select('#legend'); L.html('');
  if (CONVMODE){
    L.append('div').attr('class','title').text('Convergence overlay');
    const items=[['#2a2f37','ancestral (no trait)'],
                 ['#f0c419','tip: independent origin (hue = which origin)'],
                 ['#d8632a','internal node: more independent origins = hotter']];
    items.forEach(([c,t])=>{ const r=L.append('div').attr('class','row');
      r.append('div').attr('class','sw').style('background',c); r.append('div').text(t); });
    L.append('div').attr('class','note')
     .text('Orange arcs = cross-lineage bridges between independent origins (the convergent-pair signature). Expand a clade to refine the bridges.');
    return;
  }
  const s=document.getElementById('trait');
  L.append('div').attr('class','title').text(s.options[s.selectedIndex].text);
  if (TRAIT==='order'){ orders.forEach(o=>{ const r=L.append('div').attr('class','row');
      r.append('div').attr('class','sw').style('background',orderColor(o)); r.append('div').text(o); }); }
  else if (TRAIT==='trophic'){ [1,2,3].forEach(k=>{ const r=L.append('div').attr('class','row');
      r.append('div').attr('class','sw').style('background',trophicColor(k)); r.append('div').text(TROPHIC[k]); }); }
  else if (TRAIT==='activity'){ [1,2,3].forEach(k=>{ const r=L.append('div').attr('class','row');
      r.append('div').attr('class','sw').style('background',activityColor(k)); r.append('div').text(ACTIVITY[k]); }); }
  else { const ext=TRAIT==='body_log'?bodyExt:tempExt, sc=TRAIT==='body_log'?bodyColor:tempColor;
    L.append('div').style('height','12px').style('border-radius','3px')
      .style('background',`linear-gradient(90deg, ${sc(ext[0])}, ${sc((ext[0]+ext[1])/2)}, ${sc(ext[1])})`);
    const fmt=TRAIT==='body_log'?(v=>(10**v).toLocaleString(undefined,{maximumFractionDigits:0})+'g'):(v=>v.toFixed(0)+'°C');
    const r=L.append('div').attr('class','row').style('justify-content','space-between');
    r.append('span').text(fmt(ext[0])); r.append('span').text(fmt(ext[1])); }
}

// ---- controls --------------------------------------------------------------
document.getElementById('trait').addEventListener('change', e => {
  const v = e.target.value;
  if (v.startsWith('conv:')){ CONVMODE = v.slice(5); computeConv(); }
  else { CONVMODE = null; TRAIT = v; }
  view.selectAll('g.node circle').attr('fill', colorOf);
  drawBridges(); drawVerdict(); drawLegend();
});
document.getElementById('collapse').addEventListener('click', () => {
  root.children && root.children.forEach(collapse); update(root); setTimeout(fit, 380);
});
document.getElementById('search').addEventListener('input', e => {
  const q = e.target.value.trim().toLowerCase();
  if (!q){ view.selectAll('g.node').classed('hit', false); return; }
  let first=null;
  root.each(d => { if (d.data.name.toLowerCase().includes(q)){
    if(!first) first=d; let p=d.parent; while(p){ if(p._children){p.children=p._children;p._children=null;} p=p.parent; } } });
  update(root);
  view.selectAll('g.node').classed('hit', d => d.data.name.toLowerCase().includes(q));
  if (first) centerOn(first);
});

function fit(){
  const b = view.node().getBBox(); if(!b.width||!b.height) return;
  const s = Math.min(2, 0.9*Math.min(W/b.width, H/b.height));
  svg.transition().duration(400).call(zoom.transform,
    d3.zoomIdentity.translate(W/2 - s*(b.x+b.width/2), H/2 - s*(b.y+b.height/2)).scale(s));
}
function centerOn(d){ const t=d3.zoomTransform(svg.node());
  svg.transition().duration(400).call(zoom.transform,
    d3.zoomIdentity.translate(W/2 - t.k*d.y, H/2 - t.k*d.x).scale(t.k)); }
window.addEventListener('resize', () => { W=stage.clientWidth; H=stage.clientHeight; svg.attr('width',W).attr('height',H); });

update(root); drawLegend(); setTimeout(fit, 420);
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
