#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Null-model significance test for convergence on REAL fish.

Counting independent origins is not enough: a state could change many times just
because change is easy on a big tree. This applies the phylogenetic randomization
(PTP) test (polarity.randomization_test) to ask whether the observed number of
independent changes is what RANDOM placement of the same states on the same tree
would give.

Interpretation (stated honestly, both directions reported):
  observed << null  (p_fewer small) : changes are CLUSTERED = phylogenetic signal;
                    the independent origins are constrained/structured, not random
                    scatter -- loss is NOT just 'easy', so where it recurs it is
                    meaningful convergence.
  observed ~= null  : the state is as scattered as chance = labile; the apparent
                    convergence is no more than random lability.

Inapplicable taxa (Dollo) are removed first, so we test TRUE losses. Run on ALL
characters, nothing selected. LIMIT: the randomization null is the no-structure
extreme; testing excess homoplasy beyond a neutral Mk rate needs branch lengths,
which the Open Tree topology lacks -- not claimed.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from polarity import dollo_recode, randomization_test, analyze  # noqa: E402
from phenoscape_inapplicability import build_character, CHARACTERS  # noqa: E402


def main():
    print("=" * 82)
    print("  CONVERGENCE SIGNIFICANCE  (phylogenetic randomization / PTP, 999 perms)")
    print("=" * 82)
    print(f"{'character':<14}{'origins':>8}{'obs.steps':>10}{'null mean':>11}"
          f"{'p_fewer':>9}   verdict")
    print("-" * 82)
    for ch in CHARACTERS:
        b = build_character(ch)
        if not b:
            print(f"{ch:<14}   (too sparse / unplaced -- skipped)")
            continue
        tree, o2n, tip = b
        recoded, _inappl = dollo_recode(tree, tip, present_state="present")
        if len(set(recoded.values())) < 2 or len(recoded) < 6:
            print(f"{ch:<14}   (after inapplicable recoding: too few -- skipped)")
            continue
        # independent origins of the derived state (homology-correct count)
        res = analyze(ch, tree, recoded)
        derived = [s for s in res.states if s != res.root_state]
        origins = max((res.origins.get(s, 0) for s in derived), default=0)
        rt = randomization_test(tree, recoded, n_perm=999, seed=0)
        sig = rt["p_fewer"] < 0.05
        verdict = ("clustered/structured (real signal)" if sig
                   else "as random as chance (labile)")
        print(f"{ch:<14}{origins:>8}{rt['observed_steps']:>10}{rt['null_mean']:>11}"
              f"{rt['p_fewer']:>9}   {verdict}")

    print("\nHONEST READING:")
    print("  - p_fewer < 0.05 => the independent origins are FEWER than random placement")
    print("    would give = the losses cluster phylogenetically = constrained, not 'easy")
    print("    random' change. The convergence events are real signal above lability.")
    print("  - p_fewer ~ large => indistinguishable from random scatter; the apparent")
    print("    convergence could be mere lability. Reported as-is.")
    print("  - The null is the NO-STRUCTURE extreme. A neutral-Mk 'excess homoplasy' test")
    print("    needs branch lengths (Open Tree topology has none) -- deliberately not claimed.")


if __name__ == "__main__":
    main()
