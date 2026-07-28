#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Posterior integration (Idea 2) -- verdicts as DISTRIBUTIONS, not single-tree points.

Every bird verdict in the viewers rests on ONE posterior sample of the Jetz et al.
2012 tree. That hides phylogenetic uncertainty: a different posterior tree shifts
node ages and some relationships, which can move origin counts and p-values. This
runs the SAME gauntlet (Fitch origins -> PTP -> dated-Mk -> tier) across the first
N posterior trees from the VertLife Stage2 archive and reports each statistic as a
distribution. It is the single biggest honesty upgrade: a claim that holds on one
tree but not across the posterior is downgraded here, automatically.

DESIGN
  * The TIP SET is fixed once per overlay (every derived-state carrier + a seeded
    sample of non-carriers); each posterior tree is pruned to that same set, so the
    ONLY thing varying across runs is topology + branch lengths (the posterior),
    not the taxon sample. Divergence times remain trait-independent.
  * For each tree: independent origins, PTP p_fewer, dated-Mk p_excess, tier.
  * Report: origins (min/median/max), PTP and Mk p (median + fraction significant),
    and the tier distribution -- so "Tier 2 on 100/100 trees" vs "Tier 3 on 12/100"
    is visible.

HONESTY CONTRACT
  * Real posterior trees (VertLife) + real AVONET codings; nothing imputed.
  * Same fixed tip set across trees, so the distribution reflects tree uncertainty,
    not resampling noise. Subsampling keeps every carrier (no origin hidden).

Run:  python evolve/posterior_integration.py [n_trees]
Out:  evolve/results/posterior_integration.json
"""

from __future__ import annotations

import json
import random
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "viewer"))

from chronogram import (load_bird_posterior, all_tip_labels, prune_to,           # noqa: E402
                        normalize_unit_depth, branch_length_dict)
from polarity import analyze, count_origins, randomization_test, _leaves         # noqa: E402
from mk import excess_homoplasy_test                                             # noqa: E402
from ingest import avonet                                                         # noqa: E402
from export_bird_tree import (_niche_fn, _lifestyle_fn, _clade_fn, FLIGHTLESS_FAMS,  # noqa: E402
                              _str)

RESULTS = HERE / "results"

N_TREES = 100
N_OTHER = 1000   # non-carrier tips kept (carriers always kept), fixed across trees
N_PERM = 199
N_SIM = 99

OVERLAYS = {
    "aquatic": _lifestyle_fn("Aquatic"),
    "nectarivory": _niche_fn("Nectarivore"),
    "flightless": _clade_fn(families=FLIGHTLESS_FAMS),
}


def _norm(b):
    return str(b).replace(" ", "_")


def tip_states(df, fn) -> dict[str, str]:
    out = {}
    for b, r in df.iterrows():
        if not isinstance(b, str):
            continue
        v = fn(r)
        if v is not None:
            out[_norm(b)] = v
    return out


def _tier(origins, p_fewer, p_excess):
    if origins >= 2 and p_fewer < 0.05 and p_excess < 0.05:
        return 3
    if origins >= 2 and p_fewer < 0.05:
        return 2
    if origins >= 2:
        return 1
    return 0


def fixed_tipset(states: dict[str, str], on_any_tree: set[str], seed=0) -> set[str]:
    carriers = {t for t in states if states[t] == "yes" and t in on_any_tree}
    others = [t for t in states if states[t] == "no" and t in on_any_tree]
    rng = random.Random(seed)
    rng.shuffle(others)
    return carriers | set(others[:N_OTHER])


def gauntlet_on(tree, states, tipset):
    keep = {t for t in tipset if t in all_tip_labels(tree)}
    pruned = prune_to(tree, keep)
    if pruned is None:
        return None
    normalize_unit_depth(pruned)
    brlen = branch_length_dict(pruned)
    placed = {n.label: states[n.label] for n in _leaves(pruned) if n.label in states}
    if len(set(placed.values())) < 2 or sum(1 for v in placed.values() if v == "yes") < 2:
        return None
    analyze("t", pruned, placed)
    origins = count_origins(pruned, "yes", None)
    ptp = randomization_test(pruned, placed, n_perm=N_PERM, seed=0)
    mk = excess_homoplasy_test(pruned, placed, n_sim=N_SIM, seed=0, brlen=brlen)
    return origins, ptp["p_fewer"], mk["p_excess"], _tier(origins, ptp["p_fewer"], mk["p_excess"])


def summarize(name, rows):
    origins = [r[0] for r in rows]
    pf = [r[1] for r in rows]
    pe = [r[2] for r in rows]
    tiers = [r[3] for r in rows]
    from collections import Counter
    tier_dist = dict(sorted(Counter(tiers).items()))
    return {
        "overlay": name, "n_trees": len(rows),
        "origins_min": min(origins), "origins_median": int(statistics.median(origins)),
        "origins_max": max(origins),
        "ptp_p_fewer_median": round(statistics.median(pf), 4),
        "ptp_frac_significant": round(sum(1 for p in pf if p < 0.05) / len(pf), 3),
        "mk_p_excess_median": round(statistics.median(pe), 4),
        "mk_frac_significant": round(sum(1 for p in pe if p < 0.05) / len(pe), 3),
        "tier_distribution": tier_dist,
        "tier_mode": Counter(tiers).most_common(1)[0][0],
    }


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else N_TREES
    print(f"Loading {n} posterior bird trees (VertLife Stage2; cached after first run) ...")
    trees = load_bird_posterior(n)
    print(f"  {len(trees)} trees loaded")
    on_any = all_tip_labels(trees[0])
    df = avonet()

    out = {"n_trees": len(trees), "overlays": {}}
    for name, fn in OVERLAYS.items():
        states = tip_states(df, fn)
        tipset = fixed_tipset(states, on_any)
        print(f"\n[{name}] fixed tip set: {len(tipset)} "
              f"({sum(1 for t in tipset if states.get(t)=='yes')} carriers) "
              f"-- running gauntlet across {len(trees)} posterior trees ...")
        rows = []
        for i, tr in enumerate(trees):
            r = gauntlet_on(tr, states, tipset)
            if r is not None:
                rows.append(r)
        s = summarize(name, rows)
        out["overlays"][name] = s
        print(f"  origins {s['origins_min']}-{s['origins_max']} (median {s['origins_median']})"
              f" | PTP sig {s['ptp_frac_significant']:.0%}"
              f" | Mk p_excess median {s['mk_p_excess_median']} (sig {s['mk_frac_significant']:.0%})"
              f" | tiers {s['tier_distribution']}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "posterior_integration.json").write_text(json.dumps(out, indent=2))

    print("\n" + "=" * 74)
    print(f"POSTERIOR INTEGRATION  ({len(trees)} Jetz posterior trees; verdict distributions)")
    print("=" * 74)
    for name, s in out["overlays"].items():
        print(f"  {name:<12} origins {s['origins_min']}-{s['origins_max']} "
              f"(med {s['origins_median']}) | tier {s['tier_distribution']} "
              f"(mode T{s['tier_mode']})")
    print(f"\n  wrote {RESULTS / 'posterior_integration.json'}")


if __name__ == "__main__":
    main()
