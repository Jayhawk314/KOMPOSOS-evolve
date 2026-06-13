#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Fold HOMOLOGY into the similarity metric (on top of polarity) and show it works.

Polarity removed plesiomorphy (ancestral shared states). Homology removes the next
false signal: a DERIVED state shared BY DESCENT (two taxa inherited one loss from a
common ancestor) is not convergence. Unlike polarity (a per-state weight), homology
is PAIRWISE -- it depends on each pair's MRCA -- so it enters the metric as a mask
predicate (convergence_engine yoneda_sim `homology=`), backed by Fitch reconstruction
on the real Open Tree topology (polarity.TreeHomology). Nothing is tuned.

Demonstration: same-genus fin-loss pairs (one inherited loss) should DROP after the
homology fold; independent cross-lineage losses should hold.
"""

import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from phenoscape_ingest import resolve_entity, taxa_with_entity, build_model  # noqa: E402
from phenoscape_polarity import clean_name, binary_state, induced_subtree, CHARS  # noqa: E402
from ingest import reconcile                                                # noqa: E402
from polarity import parse_newick, analyze, TreeHomology                    # noqa: E402
from polarity_fold import polarity_weights                                  # noqa: E402


def build_all():
    """Return (tree, name_to_ott, char_tip_states, char_root)."""
    char_states, names = {}, set()
    for ch in CHARS:
        ent = resolve_entity(ch)
        if not ent:
            continue
        per = defaultdict(list)
        for _tx, tl, ph in taxa_with_entity(ent, limit=200):
            per[clean_name(tl)].append(ph)
        st = {n: binary_state(p) for n, p in per.items()}
        char_states[ch] = {n: s for n, s in st.items() if s and n}
        names |= set(char_states[ch])
    name_to_ott = reconcile(sorted(names))
    tree = parse_newick(induced_subtree(list(name_to_ott.values())))
    char_tip = {ch: {str(name_to_ott[n]): s for n, s in st.items() if n in name_to_ott}
                for ch, st in char_states.items()}
    char_root = {}
    for ch, tip in char_tip.items():
        if len(set(tip.values())) >= 2 and len(tip) >= 4:
            char_root[ch] = analyze(ch, tree, tip).root_state
    return tree, name_to_ott, char_tip, char_root


def parse_state(state: str):
    char = state.split(":", 1)[0].strip()
    low = state.lower()
    val = "absent" if ("absent" in low or "loss" in low) else (
          "present" if "present" in low else None)
    return char, val


def main():
    print("Building tree + reconstruction (cached) ...")
    tree, name_to_ott, char_tip, char_root = build_all()
    ctx = TreeHomology(tree, char_tip)

    model = build_model(CHARS, per_char_limit=120, verbose=False)
    weights = polarity_weights(model, char_root)
    # keys must be STRINGS to match TreeHomology.leaves / char_tip keys
    node_to_ott = {n: str(name_to_ott[clean_name(n)]) for n in model.taxa
                   if clean_name(n) in name_to_ott}

    def homology(a, b, state):
        char, val = parse_state(state)
        if val is None:
            return True                       # non-binary state: keep
        oa, ob = node_to_ott.get(a), node_to_ott.get(b)
        if oa is None or ob is None:
            return True                       # can't place: keep (conservative)
        ms = ctx.mrca_state(oa, ob, char)
        if ms is None:
            return True
        return ms != val                      # independent iff MRCA lacked this state

    model.set_polarity(weights)               # polarity always on (plesiomorphy)

    # taxa with pelvic-fin loss, grouped to build same-genus vs cross-genus pairs
    loss_taxa = [t for t in model.taxa
                 if any("pelvic fin" in m.target and "absent" in m.target.lower()
                        for m in model.lens("morphology").morphisms_from(t))]
    genus = lambda t: t.split()[0]
    same_genus = [(a, b) for a, b in combinations(loss_taxa, 2) if genus(a) == genus(b)]
    cross_genus = [(a, b) for a, b in combinations(loss_taxa, 2) if genus(a) != genus(b)]

    def mean_sim(pairs, n=300):
        pairs = pairs[:n]
        if not pairs:
            return float("nan"), 0
        return sum(model.similarity("morphology", a, b) for a, b in pairs) / len(pairs), len(pairs)

    print(f"\npelvic-fin-loss taxa in model: {len(loss_taxa)}  "
          f"(same-genus pairs {len(same_genus)}, cross-genus {len(cross_genus)})")

    model.set_homology(None)
    sg_before, n_sg = mean_sim(same_genus)
    xg_before, n_xg = mean_sim(cross_genus)
    model.set_homology(homology)
    sg_after, _ = mean_sim(same_genus)
    xg_after, _ = mean_sim(cross_genus)

    print("\n" + "=" * 74)
    print("  HOMOLOGY FOLDED INTO SIMILARITY  (polarity already on)")
    print("=" * 74)
    print(f"{'pair group':<32}{'sim BEFORE':>12}{'sim AFTER':>12}{'change':>9}")
    print("-" * 74)
    print(f"{'same-genus loss (inherited)':<32}{sg_before:>12.3f}{sg_after:>12.3f}"
          f"{sg_after-sg_before:>+9.3f}")
    print(f"{'cross-genus loss (independent)':<32}{xg_before:>12.3f}{xg_after:>12.3f}"
          f"{xg_after-xg_before:>+9.3f}")
    print("\nExpect: same-genus DROPS (their shared loss is one inherited event,")
    print("excluded by the MRCA test); cross-genus holds (independent losses count).")
    print("Un-fitting: MRCA states are Fitch reconstructions on a tree that never")
    print("saw the trait data; the predicate has no tunable parameter.")

    # one concrete same-genus pair with the largest homology effect (both must map)
    best = None
    for a, b in same_genus:
        if a not in node_to_ott or b not in node_to_ott:
            continue
        model.set_homology(None); s0 = model.similarity("morphology", a, b)
        model.set_homology(homology); s1 = model.similarity("morphology", a, b)
        if best is None or (s0 - s1) > best[0]:
            best = (s0 - s1, a, b, s0, s1)
    if best:
        _, a, b, s0, s1 = best
        print(f"\nexample same-genus (inherited loss): {a[:26]} ~ {b[:26]}")
        print(f"   similarity {s0:.3f} -> {s1:.3f}  (MRCA already fin-less -> shared loss excluded)")


if __name__ == "__main__":
    main()
