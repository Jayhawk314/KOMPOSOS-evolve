#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Polarity-aware convergence on REAL Phenoscape fish + a REAL Open Tree topology.

Replaces the "absent keyword = derived" heuristic (phenoscape_convergence.py
View 1) with an actual ancestral-state reconstruction. Pipeline, all sources
INDEPENDENT of each other so nothing can be fitted:

  Phenoscape  -> tip presence/absence states for each fin character
  Open Tree   -> the phylogeny (induced subtree; trait data never touches it)
  Fitch       -> parameter-free count of independent origins + consistency index

We run EVERY fin character and report whatever comes out, CI=1 (single origin,
NOT convergence) included. No character is dropped or selected to win; no
threshold is tuned; polarity (which state is ancestral) is inferred, not assumed.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from phenoscape_ingest import resolve_entity, taxa_with_entity  # noqa: E402
from ingest import reconcile                                     # noqa: E402
from polarity import parse_newick, analyze                       # noqa: E402

OTT = "https://api.opentreeoflife.org/v3"
CHARS = ["pelvic fin", "pectoral fin", "dorsal fin", "caudal fin"]


def clean_name(label: str) -> str:
    """Phenoscape label -> a name Open Tree can resolve (drop authorship / 'sp.')."""
    label = re.sub(r"\(.*?\)", "", label).strip()
    label = re.sub(r"\bsp\.?\b.*$", "", label).strip()
    return label


def binary_state(phenotype_labels: list[str]) -> str | None:
    """Reduce a taxon's annotations on one entity to absent/present (else None)."""
    joined = " ".join(p.lower() for p in phenotype_labels)
    if "absent" in joined or "loss" in joined:
        return "absent"
    if "present" in joined:
        return "present"
    return None  # 'position'/'size'/etc. -> not a presence/absence call


def induced_subtree(ott_ids: list[int]) -> str | None:
    """Open Tree induced subtree Newick; retries after pruning ids not in the tree."""
    ids = list(dict.fromkeys(ott_ids))
    for _ in range(4):
        r = requests.post(f"{OTT}/tree_of_life/induced_subtree",
                          json={"ott_ids": ids}, timeout=60)
        if r.status_code == 200:
            return r.json().get("newick")
        try:
            data = r.json()
        except Exception:
            return None
        # Open Tree returns pruned/unknown ids under "unknown" {("ott<id>": reason)}
        bad = set()
        unknown = data.get("unknown", {})
        if isinstance(unknown, dict):
            for k in unknown:
                m = re.search(r"\d+", str(k))
                if m:
                    bad.add(int(m.group()))
        # also the message may name one offending node
        m = re.search(r"ott(\d+)", data.get("message", ""))
        if m:
            bad.add(int(m.group(1)))
        if not bad:
            return None
        ids = [i for i in ids if i not in bad]
        if len(ids) < 3:
            return None
    return None


def main():
    # 1. Phenoscape: per character, taxon-label -> presence/absence
    char_taxon_state: dict[str, dict[str, str]] = {}
    all_names: set[str] = set()
    for ch in CHARS:
        ent = resolve_entity(ch)
        if not ent:
            continue
        rows = taxa_with_entity(ent, limit=200)
        per_taxon: dict[str, list[str]] = defaultdict(list)
        for _tx, tl, ph in rows:
            per_taxon[clean_name(tl)].append(ph)
        states = {}
        for name, phs in per_taxon.items():
            st = binary_state(phs)
            if st and name:
                states[name] = st
        char_taxon_state[ch] = states
        all_names |= set(states)

    print(f"Phenoscape: {len(all_names)} cleaned taxa across {len(CHARS)} fin characters")

    # 2. Open Tree: reconcile names -> ott ids (INDEPENDENT of the trait data)
    name_to_ott = reconcile(sorted(all_names))
    print(f"Open Tree TNRS resolved: {len(name_to_ott)}/{len(all_names)} names")

    # 3. Induced subtree over all resolved taxa
    newick = induced_subtree(list(name_to_ott.values()))
    if not newick:
        print("Could not build an induced subtree.")
        return
    tree = parse_newick(newick)

    # 4. Fitch polarity per character (tip states keyed by OTT id string)
    print("\n" + "=" * 78)
    print("  POLARITY-AWARE CONVERGENCE  (Fitch parsimony, parameter-free)")
    print("=" * 78)
    print(f"{'character':<14}{'tips':>5}{'steps':>7}{'CI':>6}{'root(ancestral)':>18}"
          f"{'origins(derived)':>20}")
    print("-" * 78)
    for ch in CHARS:
        states = char_taxon_state.get(ch, {})
        tip_state = {str(name_to_ott[n]): s for n, s in states.items() if n in name_to_ott}
        if len(set(tip_state.values())) < 2 or len(tip_state) < 4:
            print(f"{ch:<14}{len(tip_state):>5}   (too few/one-state on tree -- skipped)")
            continue
        res = analyze(ch, tree, tip_state)
        derived = [s for s in res.states if s != res.root_state]
        d_orig = max((res.origins.get(s, 0) for s in derived), default=0)
        verdict = "CONVERGENT" if res.is_convergent else "single origin (NOT conv.)"
        print(f"{ch:<14}{res.n_tips_scored:>5}{res.steps:>7}{res.consistency_index:>6.2f}"
              f"{res.root_state:>18}{d_orig:>12}   {verdict}")

    print("\nHonest notes:")
    print("  - steps & CI are ROOT-INDEPENDENT (the robust convergence signal);")
    print("    steps>=2 == homoplasy == convergence. Root/origins are the polarity reading.")
    print("  - Tree is Open Tree's, built from names only -- it never saw the trait states.")
    print("  - Fitch has no parameters; CI=1 outcomes are reported, not hidden.")
    print("  - Coverage caveat: only taxa Open Tree could place are scored.")


if __name__ == "__main__":
    main()
