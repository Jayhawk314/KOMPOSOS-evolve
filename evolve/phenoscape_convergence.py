#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Convergence detection on REAL Phenoscape anatomy (fish fins).

Builds a ConvergenceModel from Phenoscape (phenoscape_ingest) and looks for the
convergence signature: distant lineages sharing a DERIVED state. Fin LOSS is a
textbook convergent event in fishes (independently in eels, some catfishes,
gymnotiforms, etc.), so we focus there rather than on shared ancestral presence.

Two views:
  1. Shared-loss view (polarity-aware): for each fin-LOSS state, which taxa have
     it, and do they span different orders? Cross-order shared loss = convergence.
  2. Naive convergence ranking (trait_sim - ancestry) -- shown WITH its caveat:
     it is inflated by plesiomorphy (shared 'present'), which is why view 1 exists.
"""

import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenoscape_ingest import build_model  # noqa: E402

CHARS = ["pelvic fin", "pectoral fin", "dorsal fin", "caudal fin"]
ORDER_SUFFIX = ("iformes", "ida", "i")   # crude order-level clade pick


def order_of(model, taxon: str) -> str:
    """Pick an order-ish clade from the ancestry lineage (first '-iformes')."""
    anc = [m.target for m in model.lens("ancestry").morphisms_from(taxon)]
    for a in anc:
        if a.endswith("iformes"):
            return a
    return anc[-1] if anc else "?"


def main():
    print("Building model from Phenoscape (cached)...")
    model = build_model(CHARS, per_char_limit=120, verbose=False)
    taxa = model.taxa
    print(f"{len(taxa)} taxa, lenses {list(model.lenses)}\n")

    # ---- View 1: shared DERIVED loss across lineages -----------------------
    state_taxa: dict[str, list[str]] = defaultdict(list)
    for t in taxa:
        for m in model.lens("morphology").morphisms_from(t):
            state_taxa[m.target].append(t)

    loss_states = {s: ts for s, ts in state_taxa.items()
                   if ("absent" in s.lower() or "loss" in s.lower()) and len(ts) >= 2}

    print("=" * 76)
    print("VIEW 1  shared DERIVED loss across orders  (the convergence signal)")
    print("=" * 76)
    if not loss_states:
        print("  no shared loss states found in this character set")
    for s, ts in sorted(loss_states.items(), key=lambda x: -len(x[1])):
        orders = {order_of(model, t) for t in ts}
        cross = len(orders) > 1
        flag = "  <== CONVERGENT (multiple orders)" if cross else ""
        print(f"\n  '{s}'  ({len(ts)} taxa, {len(orders)} order(s)){flag}")
        for t in ts[:6]:
            print(f"       {t[:40]:<42} [{order_of(model, t)}]")

    # ---- View 2: naive convergence ranking (with caveat) -------------------
    print("\n" + "=" * 76)
    print("VIEW 2  naive convergence ranking  trait_sim - ancestry  (top cross-order)")
    print("=" * 76)
    scored = []
    for a, b in combinations(taxa, 2):
        if order_of(model, a) == order_of(model, b):
            continue
        conv = model.convergence("morphology", a, b)   # trait_sim - ancestry_sim
        if conv > 0:
            scored.append((conv, a, b))
    scored.sort(reverse=True)
    for conv, a, b in scored[:10]:
        print(f"  {conv:+.2f}  {a[:30]:<32} ~ {b[:30]:<32}")
    print("\n  CAVEAT: this is inflated by PLESIOMORPHY -- two distant fish sharing")
    print("  'fin present' (an ancestral trait) score as 'convergent' here. Without")
    print("  ancestral-state polarity, View 1 (shared LOSS) is the trustworthy signal.")


if __name__ == "__main__":
    main()
