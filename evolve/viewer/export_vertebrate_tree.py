#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Export a unified VERTEBRATE backbone viewer.

Splices the three clade viewers -- Actinopterygii (fish), Mammalia, and Aves --
under a single dated backbone:

    Vertebrata -> Gnathostomata
                    |-- Actinopterygii            (ray-finned fish)
                    `-- Sarcopterygii -> Tetrapoda -> Amniota
                                                        |-- Mammalia
                                                        `-- Sauropsida -> Aves

It reads the `const DATA = {...}` hierarchy already baked into each clade's HTML
(so it inherits exactly the trees those exporters built and verified) and mounts
them under the backbone. Convergence overlays stay clade-local (each clade viewer
owns its own); this view is the navigable skeleton + a Major-clade colouring,
which is the precondition for cross-class convergence work (Phase 3+).

Run (after the three clade viewers exist):
    python evolve/viewer/export_mammal_tree.py
    python evolve/viewer/export_bird_tree.py
    python evolve/viewer/export_fish_tree.py
    python evolve/viewer/export_vertebrate_tree.py
Out: evolve/viewer/vertebrate_tree.html
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from export_bird_tree import HTML_TEMPLATE                                      # noqa: E402

OUT = HERE / "vertebrate_tree.html"


def get_data(filename: str):
    """Extract the embedded `const DATA = {...};` hierarchy from a clade viewer."""
    f = HERE / filename
    if not f.exists():
        return None
    txt = f.read_text(encoding="utf-8")
    a = txt.index("const DATA = ") + len("const DATA = ")
    b = txt.index(";\nconst META", a)
    return json.loads(txt[a:b])


def count_rank(node, rank) -> int:
    n = 1 if node.get("rank") == rank else 0
    return n + sum(count_rank(c, rank) for c in node.get("children", []))


def tag_clade(node, clade):
    node["clade"] = clade
    for c in node.get("children", []):
        tag_clade(c, clade)


def main():
    print("Splicing clade viewers into the vertebrate backbone ...")
    fish = get_data("fish_tree.html")
    mammal = get_data("mammal_tree.html")
    bird = get_data("bird_tree.html")

    present = [(n, d) for n, d in [("Actinopterygii", fish), ("Mammalia", mammal),
                                   ("Aves", bird)] if d]
    if not present:
        raise SystemExit("no clade viewers found -- build the clade HTMLs first")
    for name, d in present:
        tag_clade(d, name)
        print(f"  + {name:<16} {count_rank(d, 'species'):>5} species "
              f"(from {d.get('name')})")

    children = []
    if fish:
        children.append(fish)
    amniota = []
    if mammal:
        amniota.append(mammal)
    if bird:
        amniota.append({"name": "Sauropsida", "rank": "clade", "clade": "Aves",
                        "children": [bird]})
    if amniota:
        children.append({"name": "Sarcopterygii", "rank": "clade",
            "children": [{"name": "Tetrapoda", "rank": "clade",
                "children": [{"name": "Amniota", "rank": "clade", "children": amniota}]}]})

    root = {"name": "Vertebrata", "rank": "root", "children": children}

    meta = {
        "n_species": count_rank(root, "species"),
        "n_genera": count_rank(root, "genus"),
        "n_families": count_rank(root, "family"),
        "n_orders": count_rank(root, "order"),
        "clades": [n for n, _ in present],
    }
    print(f"  unified: {meta['n_species']} species across {len(present)} classes")

    # --- substitute into the REAL template; add a 'Major clade' colouring -----
    tmpl = HTML_TEMPLATE
    base_block = ("[['order','Taxonomic order'],['body_log','Body mass (log g)'],"
                  "['hwi','Hand-Wing index (dispersal)'],\n"
                  " ['trophic','Trophic level'],['lifestyle','Primary lifestyle']]")
    assert base_block in tmpl
    tmpl = tmpl.replace(base_block, "[['clade','Major clade'],['order','Taxonomic order']]")

    # clade colour scale (define alongside orderColor)
    tmpl = tmpl.replace(
        "const orderColor = d3.scaleOrdinal(",
        "const cladeColor = d3.scaleOrdinal(['Actinopterygii','Mammalia','Aves'],"
        "['#7ec97e','#e36d6d','#58a6ff']);\nconst orderColor = d3.scaleOrdinal(")
    # colour-by branch
    tmpl = tmpl.replace(
        "  if (TRAIT === 'order') return d.orderName ? orderColor(d.orderName) : '#888';",
        "  if (TRAIT === 'clade') return d.data.clade ? cladeColor(d.data.clade) : '#888';\n"
        "  if (TRAIT === 'order') return d.orderName ? orderColor(d.orderName) : '#888';")
    # legend: handle 'clade' (default branch assumes a numeric gradient and would throw)
    tmpl = tmpl.replace(
        "  if (TRAIT==='order'){ orders.forEach(o=>{ const r=L.append('div').attr('class','row');",
        "  if (TRAIT==='clade'){ ['Actinopterygii','Mammalia','Aves'].forEach(o=>{ "
        "const r=L.append('div').attr('class','row');\n"
        "      r.append('div').attr('class','sw').style('background',cladeColor(o)); "
        "r.append('div').text(o); }); }\n"
        "  else if (TRAIT==='order'){ orders.forEach(o=>{ const r=L.append('div').attr('class','row');")

    html = (tmpl
            .replace("Bird Tree Explorer", "Vertebrate Tree Explorer")
            .replace("🐦", "🦴")
            .replace("__DATA__", json.dumps(root, separators=(",", ":")))
            .replace("__META__", json.dumps(meta))
            .replace("__CONV__", "{}"))
    # start on the clade colouring so the three classes are obvious on open
    html = html.replace("let TRAIT = 'order';", "let TRAIT = 'clade';")
    assert "__DATA__" not in html and '"Vertebrata"' in html, "substitution failed"
    OUT.write_text(html, encoding="utf-8")
    print(f"\nWrote {OUT}  ({OUT.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
