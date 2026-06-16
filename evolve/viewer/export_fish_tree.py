#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Export a self-contained, zoomable FISH tree viewer with a CONVERGENCE OVERLAY.

Phase-2/A sibling of the mammal + bird viewers, over the REAL dated ray-finned
fish tree (Rabosky et al. 2018, 11,638 species) joined with Phenoscape anatomical
annotations (fins, barbels). Same exporter pattern, same honesty contract, same
Fitch -> PTP -> dated-Mk -> tier pipeline.

PROVENANCE (honest)
-------------------
  * Every leaf is a species on the Rabosky et al. 2018 dated chronogram.
  * Order/Family for each genus is resolved via Open Tree of Life (cached).
  * Traits are Phenoscape KB presence/absence annotations on real anatomy terms.
  * Polarity is NOT assumed: Fitch reconstructs the ancestral state and we count
    independent origins of the DERIVED state (a loss for fins, a gain for barbels).
  * For LOSS characters we Dollo-recode first -- a taxon coding 'absent' OUTSIDE
    the clade that ever had the structure never had it (its absence is primitive,
    not a loss), so it is dropped, never counted as a convergent loss.
  * The convergence statistics are computed on the FULL 11.6k-tip tree; the
    interactive view shows a representative subsample (for browser speed) that
    always includes every derived-state carrier, so no origin is hidden.

Run:  python evolve/viewer/export_fish_tree.py
Out:  evolve/viewer/fish_tree.html
"""

from __future__ import annotations

import json
import random
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVOLVE = HERE.parent
sys.path.insert(0, str(EVOLVE.parents[0] / "src"))
sys.path.insert(0, str(EVOLVE))

from chronogram import (load_chronogram, all_tip_labels, prune_to,             # noqa: E402
                        normalize_unit_depth, branch_length_dict)
from polarity import (analyze, randomization_test, count_origins, _leaves,     # noqa: E402
                      dollo_recode)
from mk import excess_homoplasy_test                                           # noqa: E402
from ingest import reconcile, lineage                                          # noqa: E402
from phenoscape_ingest import resolve_entity, taxa_with_entity                 # noqa: E402
from phenoscape_polarity import clean_name, binary_state                       # noqa: E402
from export_bird_tree import HTML_TEMPLATE                                      # noqa: E402

OUT = HERE / "fish_tree.html"
CACHE = EVOLVE / "cache" / "fish_taxonomy.json"

N_PERM = 199
N_SIM = 99
VIEW_CAP = 3000   # max species rendered in the browser (stats use the full tree)

CHARACTERS = ["pelvic fin", "pectoral fin", "dorsal fin", "caudal fin",
              "adipose fin", "barbel"]


# ---------------------------------------------------------------------------
# Taxonomy (Open Tree, cached) -- unchanged, works well (97% of genera resolve)
# ---------------------------------------------------------------------------
def load_taxonomy(genera: set[str]) -> dict[str, dict]:
    tax = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    todo = [g for g in genera if g not in tax]
    if todo:
        print(f"Resolving taxonomy for {len(todo)} new genera via Open Tree ...")
        ids = reconcile(todo)
        for g in todo:
            order = family = "Incertae_sedis"
            if g in ids:
                try:
                    for rank_name in lineage(ids[g]):
                        if rank_name.endswith("iformes"):
                            order = rank_name
                        elif rank_name.endswith("idae"):
                            family = rank_name
                except Exception:
                    pass
            tax[g] = {"order": order, "family": family}
            if len(tax) % 25 == 0:
                CACHE.write_text(json.dumps(tax, indent=2))
        CACHE.write_text(json.dumps(tax, indent=2))
    return tax


# ---------------------------------------------------------------------------
# Phenoscape traits
# ---------------------------------------------------------------------------
def get_fish_traits():
    print("Fetching Phenoscape presence/absence annotations ...")
    traits: dict = defaultdict(dict)
    for ch in CHARACTERS:
        ent = resolve_entity(ch)
        if not ent:
            print(f"  [{ch}] no entity id -- skipped")
            continue
        per = defaultdict(list)
        for _t, tl, ph in taxa_with_entity(ent, limit=500):
            per[clean_name(tl)].append(ph)
        n = 0
        for name, phs in per.items():
            s = binary_state(phs)
            if s:
                traits[name][ch] = s
                n += 1
        print(f"  [{ch}] {n} species scored")
    return traits


# ---------------------------------------------------------------------------
# Hierarchy (slice 1)
# ---------------------------------------------------------------------------
def build_hierarchy(tips: set[str], tax: dict[str, dict]):
    tree: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for tip in tips:
        genus = tip.split("_")[0]
        info = tax.get(genus, {"order": "Incertae_sedis", "family": "Incertae_sedis"})
        tree[info["order"]][info["family"]][genus].append(
            {"name": tip.replace("_", " "), "rank": "species"})

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
    root = {"name": "Actinopterygii", "rank": "root", "children": roots}
    meta = {
        "n_species": len(tips),
        "n_orders": len(tree),
        "n_families": sum(len(f) for f in tree.values()),
        "n_genera": sum(len(g) for f in tree.values() for g in f.values()),
    }
    return root, meta


# ---------------------------------------------------------------------------
# Convergence overlay (slice 2) -- Fitch-derived polarity + Dollo for losses
# ---------------------------------------------------------------------------
def _origin_ids(tree, derived):
    counter = [0]
    tip_origin: dict = {}

    def walk(n, ps, cur):
        if n.state == derived and ps != derived:
            counter[0] += 1
            cur = counter[0]
        elif n.state != derived:
            cur = None
        if n.is_leaf and n.state == derived and cur is not None:
            tip_origin[n.label.replace("_", " ")] = cur
        for c in n.children:
            walk(c, n.state, cur)

    walk(tree, None, None)
    return tip_origin


def _tier(o, pf, pe):
    if o >= 2 and pf < 0.05 and pe < 0.05:
        return 3, "excess beyond neutral drift"
    if o >= 2 and pf < 0.05:
        return 2, "structured homoplasy"
    if o >= 2:
        return 1, "homoplasy present but labile"
    return 0, "single origin"


def analyze_convergence(full, traits):
    print("Analyzing convergence on the real dated fish chronogram ...")
    out = {}
    labels = all_tip_labels(full)
    for ch in CHARACTERS:
        key = ch.replace(" ", "_")
        scored = {t.replace(" ", "_"): traits[t][ch] for t in traits if ch in traits[t]}
        keep = set(scored) & labels
        if len(keep) < 10:
            print(f"  [{ch}] too few scored tips ({len(keep)}) -- skipped")
            continue
        pruned = prune_to(full, keep)
        if not pruned:
            continue
        normalize_unit_depth(pruned)
        placed = {n.label: scored[n.label] for n in _leaves(pruned) if n.label in scored}
        if len(set(placed.values())) < 2:
            print(f"  [{ch}] invariant -- skipped")
            continue

        # provisional Fitch -> which state is derived (non-ancestral)?
        res0 = analyze(ch, pruned, placed)
        derived = next((s for s in res0.states if s != res0.root_state), None)
        coded = "measured (Phenoscape)"

        # honesty: a LOSS coded outside the clade that ever had the structure is
        # primitive absence, not a convergent loss -> Dollo-recode it away.
        if derived == "absent":
            placed, inappl = dollo_recode(pruned, placed, present_state="present")
            coded = "measured (Phenoscape, Dollo-recoded losses)"
            if inappl:
                print(f"  [{ch}] Dollo dropped {len(inappl)} primitive-absent taxa")

        keep2 = set(placed)
        if len(keep2) < 10 or len(set(placed.values())) < 2:
            print(f"  [{ch}] too few after recoding -- skipped")
            continue
        pruned2 = prune_to(full, keep2)
        normalize_unit_depth(pruned2)
        brlen = branch_length_dict(pruned2)
        placed2 = {n.label: placed[n.label] for n in _leaves(pruned2) if n.label in placed}

        res = analyze(ch, pruned2, placed2)
        derived = next((s for s in res.states if s != res.root_state), derived)
        if sum(1 for v in placed2.values() if v == derived) < 2:
            print(f"  [{ch}] <2 derived carriers -- skipped")
            continue
        origins = count_origins(pruned2, derived, None)
        tip_origin = _origin_ids(pruned2, derived)

        ptp = randomization_test(pruned2, placed2, n_perm=N_PERM, seed=0)
        mk = excess_homoplasy_test(pruned2, placed2, n_sim=N_SIM, seed=0, brlen=brlen)
        tn, verdict = _tier(origins, ptp["p_fewer"], mk["p_excess"])
        derived_label = "loss" if derived == "absent" else "gain"

        print(f"  [{ch}] derived '{derived}' ({derived_label}) | origins {origins} "
              f"| PTP {ptp['p_fewer']:.3f} | Mk {mk['p_excess']:.3f} -> T{tn}")

        out[key] = {
            "key": key, "label": ch, "coded": coded,
            "derived_label": derived_label, "root_state": res.root_state, "derived": derived,
            "n_scored": len(placed2), "n_derived": len(tip_origin),
            "origins": origins, "steps": res.steps, "ci": round(res.consistency_index, 3),
            "ptp_p_fewer": round(ptp["p_fewer"], 4),
            "mk": {"obs": mk["observed_steps"], "null": round(mk["null_mean"], 2),
                   "p_excess": round(mk["p_excess"], 4), "q_hat": round(mk["q_hat"], 3)},
            "tier": tn, "tier_label": verdict, "tip_origin": tip_origin,
        }
    return out


def main():
    print("Loading dated fish chronogram (cached) ...")
    full = load_chronogram()
    tips = all_tip_labels(full)
    print(f"  {len(tips)} species in chronogram")

    traits = get_fish_traits()
    conv = analyze_convergence(full, traits)

    # subsample for the browser, but ALWAYS keep every derived-state carrier so the
    # overlays show all origins (stats above already used the full tree).
    carriers = set()
    for c in conv.values():
        carriers |= {n.replace(" ", "_") for n in c["tip_origin"]}
    carriers &= tips
    rng = random.Random(0)
    rest = list(tips - carriers)
    rng.shuffle(rest)
    sample = carriers | set(rest[:max(0, VIEW_CAP - len(carriers))])
    print(f"  rendering {len(sample)} species ({len(carriers)} trait carriers kept)")

    tax = load_taxonomy({t.split("_")[0] for t in sample})
    root, meta = build_hierarchy(sample, tax)
    meta["n_chronogram"] = len(tips)
    meta["subsampled"] = len(sample) < len(tips)
    meta["conv_traits"] = list(conv.keys())

    # substitute into the REAL template (with __DATA__/__META__/__CONV__ markers),
    # not a rendered HTML. Trim the bird-specific base-trait dropdown to Order.
    base_block = ("[['order','Taxonomic order'],['body_log','Body mass (log g)'],"
                  "['hwi','Hand-Wing index (dispersal)'],\n"
                  " ['trophic','Trophic level'],['lifestyle','Primary lifestyle']]")
    tmpl = HTML_TEMPLATE
    assert base_block in tmpl, "template base-trait block not found (template changed?)"
    tmpl = tmpl.replace(base_block, "[['order','Taxonomic order']]")

    html = (tmpl
            .replace("Bird Tree Explorer", "Fish Tree Explorer")
            .replace("🐦", "🐟")
            .replace("__DATA__", json.dumps(root, separators=(",", ":")))
            .replace("__META__", json.dumps(meta))
            .replace("__CONV__", json.dumps(conv, separators=(",", ":"))))
    assert "__DATA__" not in html and '"Actinopterygii"' in html, "substitution failed"
    OUT.write_text(html, encoding="utf-8")
    print(f"\nWrote {OUT}  ({OUT.stat().st_size/1024:.0f} KB)")
    print(f"  {meta['n_species']} species shown | {meta['n_orders']} orders | "
          f"overlays: {', '.join(conv.keys()) or '(none)'}")


if __name__ == "__main__":
    main()
