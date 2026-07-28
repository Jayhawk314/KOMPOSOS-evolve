#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Phenoscape KB ingester -- REAL ontology-annotated anatomy, no curation.

This is the data source that actually measures the convergent axis (body plan /
anatomy), unlike PanTHERIA (life-history). It pulls discrete, ontology-annotated
evolutionary characters from the Phenoscape Knowledgebase and builds a
ConvergenceModel with two lenses:

    morphology : taxon --has--> phenotype-state   (e.g. "pelvic fin absent")
    ancestry   : taxon --member--> ancestor        (walked up the VTO taxonomy)

Live API (probed, working 2026):
    /term/search?text=...           name  -> ontology IRI (Uberon/ZFA/...)
    /taxon/annotations?entity=<iri> entity-> TSV (taxon IRI, label, phenotype label)
    /term?iri=<vto>                 taxon -> classification.subClassOf (parent)

Everything is cached under evolve/cache/phenoscape/ so re-runs are offline.

HONEST CAVEAT: Phenoscape gives character STATES, not POLARITY. Two taxa sharing
"pelvic fin present" share a likely-ANCESTRAL trait (plesiomorphy, not
convergence); two sharing a DERIVED loss ("pelvic fin absent") across distant
lineages is the convergence signal. The ancestry-subtracted score surfaces the
latter, but without an ancestral-state reconstruction this is a heuristic, not
proof. Flagged, not hidden.
"""

from __future__ import annotations

import csv
import io
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from convergence_engine import ConvergenceModel  # noqa: E402

BASE = "https://kb.phenoscape.org/api/v2"
CACHE = Path(__file__).resolve().parent / "cache" / "phenoscape"
CACHE.mkdir(parents=True, exist_ok=True)


def _strip(label: str) -> str:
    """Phenoscape serialises Option[String] as 'Some(x)' / 'None'."""
    label = label.strip()
    if label.startswith("Some(") and label.endswith(")"):
        label = label[5:-1]
    return label


def _cached_get(path: str, params: dict, cache_key: str, parse="json"):
    f = CACHE / cache_key
    if f.exists():
        txt = f.read_text(encoding="utf-8")
    else:
        r = requests.get(BASE + path, params=params, timeout=60)
        r.raise_for_status()
        txt = r.text
        f.write_text(txt, encoding="utf-8")
        time.sleep(0.2)  # be polite
    return json.loads(txt) if parse == "json" else txt


# ── API calls ───────────────────────────────────────────────────────────────
def resolve_entity(name: str) -> Optional[str]:
    """Anatomy name -> ontology IRI (prefer an exact Uberon match)."""
    key = "term_search__" + name.replace(" ", "_") + ".json"
    data = _cached_get("/term/search", {"text": name, "limit": 10}, key)
    results = data.get("results", [])
    if not results:
        return None
    # prefer exact Uberon, else first exact, else first
    for r in results:
        if r.get("matchType") == "exact" and "UBERON" in r["@id"]:
            return r["@id"]
    for r in results:
        if r.get("matchType") == "exact":
            return r["@id"]
    return results[0]["@id"]


def taxa_with_entity(entity_iri: str, limit: int = 200) -> List[Tuple[str, str, str]]:
    """(taxon_iri, taxon_label, phenotype_state) for taxa annotated on this entity."""
    key = "ann__" + entity_iri.rsplit("/", 1)[-1] + f"__{limit}.tsv"
    txt = _cached_get("/taxon/annotations",
                      {"entity": entity_iri, "limit": limit}, key, parse="tsv")
    rows = csv.DictReader(io.StringIO(txt), delimiter="\t")
    out = []
    for row in rows:
        tx = _strip(row.get("taxon IRI", ""))
        tl = _strip(row.get("taxon label", ""))
        ph = _strip(row.get("phenotype label", ""))
        if tx and tl and ph:
            out.append((tx, tl, ph))
    return out


def lineage(taxon_iri: str, max_depth: int = 14) -> List[str]:
    """Walk VTO subClassOf up; return ancestor labels (genus..order..class)."""
    key = "lin__" + taxon_iri.rsplit("/", 1)[-1] + ".json"
    cache = CACHE / key
    if cache.exists():
        return json.loads(cache.read_text())
    anc: List[str] = []
    cur = taxon_iri
    seen = set()
    for _ in range(max_depth):
        if cur in seen:
            break
        seen.add(cur)
        try:
            data = _cached_get("/term", {"iri": cur},
                               "term__" + cur.rsplit("/", 1)[-1] + ".json")
        except Exception:
            break
        parents = (data.get("classification", {}) or {}).get("subClassOf", []) or []
        if not parents:
            break
        p = parents[0]
        lbl = _strip(p.get("label", ""))
        if lbl:
            anc.append(lbl)
        cur = p["@id"]
        if lbl in ("Vertebrata", "Gnathostomata", "Chordata"):
            break
    cache.write_text(json.dumps(anc))
    return anc


# ── Build a ConvergenceModel from a character set ────────────────────────────
def build_model(character_names: List[str], per_char_limit: int = 200,
                min_states_per_taxon: int = 2, verbose: bool = True) -> ConvergenceModel:
    """Pull the given anatomical characters and assemble a ConvergenceModel.

    morphology lens : taxon --has--> "<entity>: <state>"
    ancestry lens   : taxon --member--> ancestor (VTO walk)
    Only taxa with >= min_states_per_taxon recorded states are kept (so the
    similarity has something to compare).
    """
    m = ConvergenceModel()
    taxon_states: Dict[str, set] = {}
    taxon_label: Dict[str, str] = {}

    for name in character_names:
        ent = resolve_entity(name)
        if not ent:
            if verbose:
                print(f"  [skip] no ontology term for '{name}'")
            continue
        rows = taxa_with_entity(ent, limit=per_char_limit)
        for tx, tl, ph in rows:
            state = f"{name}: {ph}"
            taxon_states.setdefault(tx, set()).add(state)
            taxon_label[tx] = tl
        if verbose:
            print(f"  {name:<16} -> {ent.rsplit('/',1)[-1]:<16} {len(rows)} annotations")

    kept = [tx for tx, st in taxon_states.items() if len(st) >= min_states_per_taxon]
    if verbose:
        print(f"\nTaxa with >= {min_states_per_taxon} states: {len(kept)} "
              f"(of {len(taxon_states)} annotated)")

    for tx in kept:
        node = taxon_label[tx]
        m.add_taxon(node)
        for state in taxon_states[tx]:
            m.relate("morphology", node, state, rel="has")
        for ancestor in lineage(tx):
            m.relate("ancestry", node, ancestor, rel="member")

    return m


if __name__ == "__main__":
    # smoke test: a few fin/limb characters where loss is famously convergent
    chars = ["pelvic fin", "pectoral fin", "dorsal fin", "caudal fin"]
    print("Phenoscape ingest smoke test:")
    model = build_model(chars, per_char_limit=120)
    print(f"\nConvergenceModel: {len(model.taxa)} taxa, "
          f"lenses={list(model.lenses)}")
    if model.taxa:
        a = model.taxa[0]
        print(f"example taxon: {a}")
        print(f"  morphology states: "
              f"{[m.target for m in model.lens('morphology').morphisms_from(a)][:5]}")
        print(f"  ancestry: "
              f"{[m.target for m in model.lens('ancestry').morphisms_from(a)][:6]}")
