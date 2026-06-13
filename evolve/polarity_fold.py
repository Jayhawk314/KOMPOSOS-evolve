#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Fold ancestral-state polarity into the MAIN detectors, and show it works.

Pipeline (all sources independent -> un-fittable; see avoid-fitting-queries):
  1. polarity (Fitch on the real Open Tree topology) -> root state per character
  2. state_weight: ancestral(root) trait-states -> 0, derived states -> 1
  3. model.set_polarity(state_weight)  -> Yoneda/Ricci/spectral all become
     plesiomorphy-aware at once (they route through similarity)

Demonstration: the cross-lineage convergence ranking BEFORE polarity is inflated by
shared 'fin present' (ancestral); AFTER, it is dominated by shared DERIVED loss --
the View-2 plesiomorphy bug from phenoscape_convergence.py, fixed in the engine.

We do NOT hand-pick which states to down-weight: the parameter-free reconstruction
on a tree that never saw the traits decides which state is ancestral.
"""

import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from phenoscape_ingest import resolve_entity, taxa_with_entity, build_model  # noqa: E402
from phenoscape_polarity import clean_name, binary_state, induced_subtree, CHARS  # noqa: E402
from phenoscape_convergence import order_of                                  # noqa: E402
from ingest import reconcile                                                 # noqa: E402
from polarity import parse_newick, analyze                                   # noqa: E402


def character_roots() -> dict[str, str]:
    """Inferred ancestral (root) state per fin character, via Fitch on Open Tree."""
    char_states, names = {}, set()
    for ch in CHARS:
        ent = resolve_entity(ch)
        if not ent:
            continue
        per = defaultdict(list)
        for _tx, tl, ph in taxa_with_entity(ent, limit=200):
            per[clean_name(tl)].append(ph)
        states = {n: binary_state(p) for n, p in per.items()}
        states = {n: s for n, s in states.items() if s and n}
        char_states[ch] = states
        names |= set(states)
    name_to_ott = reconcile(sorted(names))
    newick = induced_subtree(list(name_to_ott.values()))
    tree = parse_newick(newick)
    roots = {}
    for ch, states in char_states.items():
        tip = {str(name_to_ott[n]): s for n, s in states.items() if n in name_to_ott}
        if len(set(tip.values())) >= 2 and len(tip) >= 4:
            roots[ch] = analyze(ch, tree, tip).root_state
    return roots


def polarity_weights(model, char_root: dict[str, str]) -> dict[str, float]:
    """Map every morphology trait-state -> 0 (reconstructed ancestral) or 1 (derived)."""
    weights = {}
    for t in model.taxa:
        for m in model.lens("morphology").morphisms_from(t):
            s = m.target
            if s in weights:
                continue
            char = s.split(":", 1)[0]
            low = s.lower()
            val = "absent" if ("absent" in low or "loss" in low) else (
                  "present" if "present" in low else None)
            root = char_root.get(char)
            # ancestral (matches reconstructed root) -> 0 ; derived/other -> 1
            weights[s] = 0.0 if (val is not None and val == root) else 1.0
    return weights


def top_convergent(model, k=10):
    scored = []
    for a, b in combinations(model.taxa, 2):
        if order_of(model, a) == order_of(model, b):
            continue
        c = model.convergence("morphology", a, b)
        if c > 0:
            scored.append((c, a, b))
    scored.sort(reverse=True)
    return scored[:k]


def main():
    print("Computing polarity (Fitch on Open Tree) ...")
    char_root = character_roots()
    print("  inferred ancestral state per character:",
          {c: r for c, r in char_root.items()})

    model = build_model(CHARS, per_char_limit=120, verbose=False)
    weights = polarity_weights(model, char_root)
    n_anc = sum(1 for w in weights.values() if w == 0.0)
    print(f"  trait-states: {len(weights)}  (ancestral/down-weighted: {n_anc}, "
          f"derived: {len(weights)-n_anc})\n")

    before = top_convergent(model)
    model.set_polarity(weights)
    after = top_convergent(model)

    def show(title, rows):
        print(title)
        for c, a, b in rows:
            print(f"   {c:+.2f}  {a[:26]:<28}[{order_of(model,a)[:14]:<14}] ~ "
                  f"{b[:22]:<24}[{order_of(model,b)[:14]}]")

    print("=" * 78)
    show("BEFORE polarity (top cross-order 'convergent' pairs) -- plesiomorphy-inflated:",
         before)
    print()
    show("AFTER polarity (ancestral states down-weighted) -- derived-sharing rises:",
         after)

    # aggregate evidence the fold did what it should
    def shares_only_present(a, b):
        sa = {m.target for m in model.lens("morphology").morphisms_from(a)}
        sb = {m.target for m in model.lens("morphology").morphisms_from(b)}
        shared = sa & sb
        return shared and all("present" in s.lower() for s in shared)

    plesio = [(a, b) for a, b in combinations(model.taxa, 2)
              if order_of(model, a) != order_of(model, b) and shares_only_present(a, b)]
    if plesio:
        model.state_weight = {}
        b_mean = sum(model.similarity("morphology", a, b) for a, b in plesio[:200]) / min(200, len(plesio))
        model.set_polarity(weights)
        a_mean = sum(model.similarity("morphology", a, b) for a, b in plesio[:200]) / min(200, len(plesio))
        print("\n" + "-" * 78)
        print(f"Cross-order pairs sharing ONLY 'present' (pure plesiomorphy), n={len(plesio)}:")
        print(f"   mean morphology similarity BEFORE polarity: {b_mean:.3f}")
        print(f"   mean morphology similarity AFTER  polarity: {a_mean:.3f}  "
              f"(should collapse toward 0 -- shared ancestry is no longer 'convergence')")
    print("\nPolarity from parameter-free Fitch on a tree that never saw the traits.")


if __name__ == "__main__":
    main()
