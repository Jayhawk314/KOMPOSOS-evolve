#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Tiered convergence verdict with REAL DATED branch lengths (Fish Tree of Life).

Same COG-style tiered test as phenoscape_excess.py, but the neutral-Mk null now
runs on the time-calibrated fish chronogram (real divergence times, MY) instead of
Grafen topology branch lengths. Divergence times are independent of the trait data,
so this strengthens the null without fitting.

Coverage caveat: only taxa present in the fish chronogram are scored (extinct /
non-actinopterygian / higher taxa drop out). Reported honestly.
"""

import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from phenoscape_ingest import resolve_entity, taxa_with_entity   # noqa: E402
from phenoscape_polarity import clean_name, binary_state          # noqa: E402
from phenoscape_inapplicability import CHARACTERS                  # noqa: E402
from polarity import dollo_recode, randomization_test, analyze     # noqa: E402
from mk import excess_homoplasy_test                               # noqa: E402
from chronogram import (load_chronogram, all_tip_labels, prune_to,  # noqa: E402
                        normalize_unit_depth, branch_length_dict)


def char_states(ch):
    ent = resolve_entity(ch)
    if not ent:
        return {}
    per = defaultdict(list)
    for _t, tl, ph in taxa_with_entity(ent, limit=200):
        per[clean_name(tl)].append(ph)
    st = {n: binary_state(p) for n, p in per.items()}
    return {n: s for n, s in st.items() if s and n}


def tier(origins, p_fewer, p_excess):
    if origins >= 2 and p_fewer < 0.05 and p_excess < 0.05:
        return 3, "EXCESS homoplasy beyond neutral drift (strong)"
    if origins >= 2 and p_fewer < 0.05:
        return 2, "structured homoplasy (real, not excess)"
    if origins >= 2:
        return 1, "homoplasy present but labile"
    return 0, "single origin"


def main():
    print("Loading fish chronogram (cached) ...")
    full = load_chronogram()
    labels = all_tip_labels(full)
    print(f"chronogram: {len(labels)} dated tips\n")

    print("=" * 92)
    print("  TIERED CONVERGENCE VERDICT with REAL DATED branch lengths (Fish Tree of Life)")
    print("=" * 92)
    print(f"{'character':<13}{'matched':>8}{'orig':>5}{'PTP p_few':>10}{'q_hat':>7}"
          f"{'Mk obs/null':>13}{'p_excess':>9}  tier")
    print("-" * 92)
    for ch in CHARACTERS:
        states = char_states(ch)
        if not states:
            continue
        keep = {n.replace(" ", "_") for n in states} & labels
        if len(keep) < 6:
            print(f"{ch:<13}{len(keep):>8}   (too few on chronogram -- skipped)")
            continue
        pruned = prune_to(full, keep)
        normalize_unit_depth(pruned)
        brlen = branch_length_dict(pruned)
        tip = {us: states[us.replace('_', ' ')] for us in keep}

        recoded, _ = dollo_recode(pruned, tip, present_state="present")
        if len(set(recoded.values())) < 2 or len(recoded) < 6:
            print(f"{ch:<13}{len(keep):>8}   (one-state after recoding -- skipped)")
            continue
        res = analyze(ch, pruned, recoded)
        derived = [s for s in res.states if s != res.root_state]
        origins = max((res.origins.get(s, 0) for s in derived), default=0)
        ptp = randomization_test(pruned, recoded, n_perm=999, seed=0)
        mk = excess_homoplasy_test(pruned, recoded, n_sim=499, seed=0, brlen=brlen)
        tn, verdict = tier(origins, ptp["p_fewer"], mk["p_excess"])
        print(f"{ch:<13}{len(keep):>8}{origins:>5}{ptp['p_fewer']:>10}{mk['q_hat']:>7}"
              f"{str(mk['observed_steps'])+'/'+str(mk['null_mean']):>13}{mk['p_excess']:>9}   T{tn}")

    print("\nHONEST NOTES:")
    print("  - Null now uses REAL divergence times (treePL chronogram), independent of")
    print("    the trait data -> the strongest, least-assumption neutral null we can build.")
    print("  - Coverage drops vs Open Tree: only ray-finned fish SPECIES in the chronogram")
    print("    are scored (extinct/non-teleost/higher taxa excluded). Reported, not hidden.")
    print("  - Mk rate still ML-estimated under the no-convergence null (conservative).")


if __name__ == "__main__":
    main()
