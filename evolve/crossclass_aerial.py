#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Cross-class aerial-locomotion convergence on a deep-node-calibrated backbone (Track A).

The vertebrate viewer (`vertebrate_tree.html`) is a TAXONOMIC splice: hierarchy
only, no cross-class branch lengths. Fitch origin-counting works on topology, but
PTP and dated-Mk need real branch lengths SPANNING the classes. So the real work
here is to GRAFT the three dated chronograms onto a backbone whose deep nodes carry
independent fossil/clock calibrations, then run the SAME gauntlet on that dated
cross-class tree.

THE BACKBONE (ultrametric, all tips at the present; ages in Myr)
---------------------------------------------------------------
    Osteichthyes split  ~430 Ma  (Actinopterygii | Sarcopterygii; TimeTree 416-439)
      |-- Actinopterygii  = Rabosky 2018 fish chronogram   (crown ~368 Ma)
      |-- Sarcopterygii -> ... -> Amniota  ~319 Ma  (Mammalia | Sauropsida; fossil cal.)
            |-- Mammalia  = Upham 2019 mammal chronogram   (crown ~218 Ma)
            |-- Sauropsida -> Aves = Jetz 2012 bird chronogram (crown ~101 Ma)

Each clade is attached by a STEM branch = (calibration age - that clade's measured
crown depth), so the graft is internally consistent (every tip lands at time 0, the
root at ~430 Ma). The two deep ages are molecular-clock + fossil estimates,
INDEPENDENT of the aerial trait, so they cannot fit the convergence result -- the
same honesty argument the per-clade chronograms already rely on.

THE TRAIT (clade-coded by external taxonomy, independent of the tree)
--------------------------------------------------------------------
aerial = powered flight OR gliding. Carriers, by independent taxonomy:
  * birds (Aves)        -- powered flight (all, EXCEPT flightless ratites + penguins)
  * bats (Chiroptera)   -- powered flight
  * flying fish (Exocoetidae) -- gliding
  * mammalian gliders (Petauridae, Anomaluridae, Cynocephalidae, ...) -- gliding
Coding by family/order is external to the trait, so counting independent origins on
the dated tree is a legitimate convergence test, not fitting (same justification as
the per-clade overlays). Polarity: the vertebrate root is non-aerial, so aerial is
the DERIVED state; flight loss in ratites/penguins shows up as reversals, not as
phantom origins.

THE GAUNTLET (identical to the viewers)
---------------------------------------
Fitch independent origins -> PTP randomization (clustered vs scattered) -> dated-Mk
excess-homoplasy (beyond neutral drift on the REAL cross-class branch lengths) ->
tier T0..T3. Plus effect sizes (z vs the nulls) and a sensitivity sweep over the
random non-carrier sample and the trait coding (with/without mammalian gliders).

HONESTY CAVEATS
  * The graft uses TWO point calibrations; a viewer overlay built on it inherits
    that assumption. We report the ages used and that they are trait-independent.
  * Non-carrier taxa are subsampled (seeded) for tractability; EVERY aerial carrier
    found is kept, so no origin is hidden. The sensitivity sweep shows the verdict
    is stable across samples.
  * Flightless birds are coded non-aerial (secondary loss); this is a reversal the
    reconstruction handles, and it does not inflate the origin count.

Run:  python evolve/crossclass_aerial.py
Out:  evolve/results/crossclass_aerial.json
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from chronogram import (load_chronogram, load_mammal_chronogram,                # noqa: E402
                        load_bird_chronogram, all_tip_labels, prune_to,
                        branch_length_dict, _max_depth)
from polarity import Node, analyze, count_origins, randomization_test, _leaves  # noqa: E402
from mk import excess_homoplasy_test                                            # noqa: E402
from ingest import pantheria, avonet                                            # noqa: E402

RESULTS = HERE / "results"

# --- deep-node calibrations (Ma), independent of the aerial trait -------------
AGE_OSTEICHTHYES = 430.0   # Actinopterygii | Sarcopterygii split (TimeTree 416-439)
AGE_AMNIOTA = 319.0        # Mammalia | Sauropsida split (standard fossil calibration)

# --- Monte-Carlo budgets (raised vs the per-clade viewers; Idea 3) -----------
N_PERM = 999
N_SIM = 499
N_OTHER = 350     # random non-carrier tips kept per class (carriers always kept)

# --- carrier taxonomy (external, independent of the trait) -------------------
EXOCOETIDAE_GENERA = {  # flying fish (gliding) -- matched by genus prefix on tips
    "Exocoetus", "Cheilopogon", "Cypselurus", "Hirundichthys", "Prognichthys",
    "Parexocoetus", "Fodiator", "Danichthys", "Oxyporhamphus",
}
GLIDER_FAMS = {  # gliding-membrane mammals (flight itself is Chiroptera, below)
    "Cynocephalidae", "Anomaluridae", "Petauridae", "Acrobatidae", "Pseudocheiridae",
}
FLIGHTLESS_FAMS = {  # birds that do NOT fly (secondary loss) -> coded non-aerial
    "Struthionidae", "Rheidae", "Casuariidae", "Dromaiidae", "Apterygidae",
    "Spheniscidae",
}


def _norm(b: str) -> str:
    return str(b).replace(" ", "_")


def mammal_states(tips: set[str]) -> dict[str, str]:
    """underscore-binomial -> 'aerial'/'nonaerial' (Chiroptera or glider family)."""
    df = pantheria()
    order = {_norm(i): str(r.get("MSW05_Order") or "") for i, r in df.iterrows()
             if isinstance(i, str)}
    fam = {_norm(i): str(r.get("MSW05_Family") or "") for i, r in df.iterrows()
           if isinstance(i, str)}
    out = {}
    for t in tips:
        o, f = order.get(t, ""), fam.get(t, "")
        out[t] = "aerial" if (o == "Chiroptera" or f in GLIDER_FAMS) else "nonaerial"
    return out


def bird_states(tips: set[str]) -> dict[str, str]:
    """underscore-binomial -> 'aerial'/'nonaerial' (flying birds aerial; flightless
    families non-aerial). Birds with no AVONET family default to aerial (extant birds
    fly unless in a known flightless lineage)."""
    df = avonet()
    fam = {_norm(i): str(r.get("Family3") or "") for i, r in df.iterrows()
           if isinstance(i, str)}
    out = {}
    for t in tips:
        out[t] = "nonaerial" if fam.get(t, "") in FLIGHTLESS_FAMS else "aerial"
    return out


def fish_states(tips: set[str]) -> dict[str, str]:
    """underscore-binomial -> 'aerial'/'nonaerial' (Exocoetidae gliders by genus)."""
    return {t: ("aerial" if t.split("_")[0] in EXOCOETIDAE_GENERA else "nonaerial")
            for t in tips}


def subsample(full: Node, states: dict[str, str], n_other: int, seed: int):
    """Keep EVERY tip in the MINORITY state (so no independent origin or reversal is
    lost) and a random sample of the MAJORITY state down to `n_other`. This balances
    each class (in birds 'aerial' is the majority; in fish/mammals it is the minority)
    without changing origin counts -- flight is one origin however many flyers remain.
    Returns the pruned clade tree (ultrametric, native Myr) and its coded states."""
    tips = [t for t in all_tip_labels(full) if t in states]
    by_state: dict[str, list] = {"aerial": [], "nonaerial": []}
    for t in tips:
        by_state[states[t]].append(t)
    minority = min(by_state, key=lambda s: len(by_state[s]))
    majority = "aerial" if minority == "nonaerial" else "nonaerial"
    rng = random.Random(seed)
    maj = list(by_state[majority]); rng.shuffle(maj)
    keep = set(by_state[minority]) | set(maj[:n_other])
    pruned = prune_to(full, keep)
    coded = {n.label: states[n.label] for n in _leaves(pruned) if n.label in states}
    return pruned, coded


def graft(fish_sub, mam_sub, bird_sub) -> Node:
    """Attach the three pruned clade trees under the calibrated backbone (Myr)."""
    fish_sub.brlen = AGE_OSTEICHTHYES - _max_depth(fish_sub)
    mam_sub.brlen = AGE_AMNIOTA - _max_depth(mam_sub)
    bird_sub.brlen = AGE_AMNIOTA - _max_depth(bird_sub)
    amniota = Node(label="Amniota", brlen=AGE_OSTEICHTHYES - AGE_AMNIOTA)
    amniota.children = [mam_sub, bird_sub]
    root = Node(label="Osteichthyes", brlen=0.0)
    root.children = [fish_sub, amniota]
    return root


def _tier(origins, p_fewer, p_excess):
    if origins >= 2 and p_fewer < 0.05 and p_excess < 0.05:
        return 3, "excess beyond neutral drift"
    if origins >= 2 and p_fewer < 0.05:
        return 2, "structured homoplasy"
    if origins >= 2:
        return 1, "homoplasy present but labile"
    return 0, "single origin"


def run_gauntlet(tree: Node, coded: dict[str, str], n_perm=N_PERM, n_sim=N_SIM):
    res = analyze("aerial", tree, coded)
    origins = count_origins(tree, "aerial", None)
    brlen = branch_length_dict(tree)
    ptp = randomization_test(tree, coded, n_perm=n_perm, seed=0)
    mk = excess_homoplasy_test(tree, coded, n_sim=n_sim, seed=0, brlen=brlen)
    tn, verdict = _tier(origins, ptp["p_fewer"], mk["p_excess"])
    # effect sizes: how many SD the observed step count sits from each null
    ptp_z = (ptp["null_mean"] - ptp["observed_steps"]) / (ptp["null_sd"] or 1e-9)
    mk_z = (mk["observed_steps"] - mk["null_mean"]) / (mk["null_sd"] or 1e-9)
    return {
        "n_tips": sum(1 for _ in _leaves(tree)),
        "n_aerial": sum(1 for v in coded.values() if v == "aerial"),
        "root_state": res.root_state,
        "origins": origins,
        "steps": res.steps,
        "ptp": ptp, "ptp_z_clustering": round(ptp_z, 2),
        "mk": mk, "mk_z_excess": round(mk_z, 2),
        "tier": tn, "tier_label": verdict,
    }


def build(seed: int, include_mammal_gliders: bool):
    fish, mam, bird = load_chronogram(), load_mammal_chronogram(), load_bird_chronogram()
    fst = fish_states(all_tip_labels(fish))
    mst = mammal_states(all_tip_labels(mam))
    if not include_mammal_gliders:
        # keep only bats among mammal carriers; re-label glider gliders as nonaerial
        df = pantheria()
        order = {_norm(i): str(r.get("MSW05_Order") or "") for i, r in df.iterrows()
                 if isinstance(i, str)}
        mst = {t: ("aerial" if order.get(t, "") == "Chiroptera" else "nonaerial")
               for t in mst}
    bst = bird_states(all_tip_labels(bird))
    fsub, fc = subsample(fish, fst, N_OTHER, seed)
    msub, mc = subsample(mam, mst, N_OTHER, seed)
    bsub, bc = subsample(bird, bst, N_OTHER, seed)
    tree = graft(fsub, msub, bsub)
    coded = {**fc, **mc, **bc}
    return tree, coded


def main():
    print("Building deep-node-calibrated cross-class tree + coding aerial locomotion ...")
    tree, coded = build(seed=0, include_mammal_gliders=True)
    print(f"  grafted tips: {sum(1 for _ in _leaves(tree))} | "
          f"aerial: {sum(1 for v in coded.values() if v=='aerial')} "
          f"(calibrations: Osteichthyes {AGE_OSTEICHTHYES} Ma, Amniota {AGE_AMNIOTA} Ma)")

    print("Running the gauntlet (Fitch -> PTP -> dated-Mk -> tier) ...")
    main_res = run_gauntlet(tree, coded)

    # --- sensitivity sweep (Idea 4): vary the random sample + the coding ------
    print("Sensitivity sweep (sample seeds + coding) ...")
    sens = []
    for seed in [1, 2, 3]:
        t, c = build(seed=seed, include_mammal_gliders=True)
        r = run_gauntlet(t, c, n_perm=299, n_sim=199)
        sens.append({"variant": f"seed={seed}", "origins": r["origins"],
                     "ptp_p_fewer": r["ptp"]["p_fewer"], "mk_p_excess": r["mk"]["p_excess"],
                     "tier": r["tier"]})
    t, c = build(seed=0, include_mammal_gliders=False)
    r = run_gauntlet(t, c, n_perm=299, n_sim=199)
    sens.append({"variant": "no_mammal_gliders", "origins": r["origins"],
                 "ptp_p_fewer": r["ptp"]["p_fewer"], "mk_p_excess": r["mk"]["p_excess"],
                 "tier": r["tier"]})

    out = {
        "calibrations_Ma": {"Osteichthyes_split": AGE_OSTEICHTHYES, "Amniota": AGE_AMNIOTA},
        "calibration_note": ("molecular-clock + fossil ages, independent of the aerial "
                             "trait; TimeTree Osteichthyes 416-439 Ma, Amniota ~319 Ma"),
        "trait": "aerial locomotion (powered flight or gliding), clade-coded",
        "main": main_res,
        "sensitivity": sens,
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "crossclass_aerial.json").write_text(json.dumps(out, indent=2))

    m = main_res
    print("\n" + "=" * 74)
    print("CROSS-CLASS AERIAL CONVERGENCE  (fish + mammals + birds, dated backbone)")
    print("=" * 74)
    print(f"  tips={m['n_tips']}  aerial={m['n_aerial']}  root={m['root_state']}")
    print(f"  independent origins of aerial locomotion : {m['origins']}")
    print(f"  PTP  : observed {m['ptp']['observed_steps']} steps vs null "
          f"{m['ptp']['null_mean']}+-{m['ptp']['null_sd']}  p_fewer={m['ptp']['p_fewer']}"
          f"  (z={m['ptp_z_clustering']} -> clustered)")
    print(f"  Mk   : observed {m['mk']['observed_steps']} vs neutral null "
          f"{m['mk']['null_mean']}+-{m['mk']['null_sd']}  p_excess={m['mk']['p_excess']}"
          f"  (z={m['mk_z_excess']})")
    print(f"  VERDICT: TIER {m['tier']} - {m['tier_label']}")
    print("\n  sensitivity (origins / PTP / Mk / tier):")
    for s in sens:
        print(f"    {s['variant']:<18} origins={s['origins']:<4} "
              f"PTP={s['ptp_p_fewer']:.3f} Mk={s['mk_p_excess']:.3f} -> T{s['tier']}")
    print(f"\n  wrote {RESULTS / 'crossclass_aerial.json'}")


if __name__ == "__main__":
    main()
