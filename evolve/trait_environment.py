#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Trait-ENVIRONMENT convergence on real mammals, with REAL DATED branch lengths.

Closes the loop to the original framing (environment as the driver) AND applies the
same dated-chronogram refinement we did for fish. The tree + divergence times are the
PHYLACINE/Upham 2019 dated mammal phylogeny (independent of the trait data); the test
(mk.trait_environment_test) asks whether the trait co-occurs with the environment MORE
than a NEUTRAL Mk process on that dated tree would produce.

PRE-SPECIFIED HYPOTHESIS (no scan -> no p-hacking): BERGMANN'S RULE -- larger body in
colder climates. trait = body mass above median ('large'); env = mean temp below median
('cold'); median splits (not tuned). Bergmann => positive phi(large, cold).

Controls bracket the test: permuted-env (should be non-significant) and env-tracks-trait
(should be significant). Whatever the verdict on Bergmann, it is reported.
"""

import sys
import random
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ingest import pantheria                                   # noqa: E402
from chronogram import load_mammal_chronogram, all_tip_labels, prune_dated  # noqa: E402
from polarity import _leaves                                    # noqa: E402
from mk import trait_environment_test                           # noqa: E402

BODY = "5-1_AdultBodyMass_g"
TEMP = "28-2_Temp_Mean_01degC"
MISSING = -999.0


def main():
    df = pantheria()
    have = df[(df[BODY] != MISSING) & (df[TEMP] != MISSING)]

    print("Loading dated mammal chronogram (Upham 2019, cached) ...")
    full = load_mammal_chronogram()
    chron = all_tip_labels(full)

    matched = [n for n in have.index if n.replace(" ", "_") in chron]
    rng = random.Random(0)
    sample = matched if len(matched) <= 600 else rng.sample(matched, 600)
    pruned, brlen, keep = prune_dated(full, sample)
    placed = [n for n in sample if n.replace(" ", "_") in keep]
    print(f"PanTHERIA body+temp: {len(have)}; on dated chronogram: {len(matched)}; "
          f"used: {len(placed)}")

    bodies = {n: have.at[n, BODY] for n in placed}
    temps = {n: have.at[n, TEMP] for n in placed}
    bmed, tmed = median(bodies.values()), median(temps.values())
    trait_state = {n.replace(" ", "_"): ("large" if bodies[n] > bmed else "small")
                   for n in placed}
    env_cold = {n.replace(" ", "_"): (temps[n] < tmed) for n in placed}
    print(f"median body {bmed:.0f} g, median temp {tmed/10:.1f} C; "
          f"dated tips used: {len(list(_leaves(pruned)))}\n")

    print("=" * 80)
    print("  BERGMANN'S RULE on the DATED mammal tree  (large body <-> cold climate)")
    print("=" * 80)
    res = trait_environment_test(pruned, trait_state, env_cold, derived_state="large",
                                 n_sim=999, seed=1, brlen=brlen)
    n11, n10, n01, n00 = res["counts_n11_n10_n01_n00"]
    print(f"  large&cold={n11}  large&warm={n10}  small&cold={n01}  small&warm={n00}")
    print(f"  observed phi(large,cold) = {res['observed_phi']:+.3f}")
    print(f"  neutral-Mk null phi       = {res['null_phi_mean']:+.3f}  (q_hat={res['q_hat']})")
    print(f"  p(assoc beyond phylogeny) = {res['p_assoc']}")
    sig = res["p_assoc"] < 0.05
    print(f"  => Bergmann {'SUPPORTED' if sig else 'NOT supported'} (p {'<' if sig else '>='} 0.05)")

    keys = list(env_cold)
    vals = list(env_cold.values())
    random.Random(7).shuffle(vals)
    ctrl = trait_environment_test(pruned, trait_state, dict(zip(keys, vals)),
                                  derived_state="large", n_sim=999, seed=1, brlen=brlen)
    rngp = random.Random(3)
    env_pos = {k: (trait_state[k] == "large") != (rngp.random() < 0.20) for k in keys}
    pos = trait_environment_test(pruned, trait_state, env_pos, derived_state="large",
                                 n_sim=999, seed=1, brlen=brlen)
    print(f"\n  NEGATIVE CONTROL (permuted env): phi={ctrl['observed_phi']:+.3f} "
          f"p={ctrl['p_assoc']}  (expect non-sig)")
    print(f"  POSITIVE CONTROL (env~trait+noise): phi={pos['observed_phi']:+.3f} "
          f"p={pos['p_assoc']}  (expect sig)")

    print("\nHONEST NOTES:")
    print("  - Tree + divergence times: dated Upham 2019 mammal phylogeny, independent of")
    print("    the trait data. Null = trait under neutral Mk on the DATED branches.")
    print("  - Earlier Grafen/Open-Tree run also found Bergmann NOT supported; the dated")
    print("    chronogram tests robustness of that to real branch lengths.")
    print("  - Pre-specified hypothesis + median splits; controls bracket the test.")


if __name__ == "__main__":
    main()
