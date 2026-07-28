#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Environment-as-driver SCAN (Track B) -- a GENERAL test, not one cherry-picked case.

The thesis question for the environment layer is: do a trait's independent origins
coincide with the SAME environmental shift -- convergence WITH A SHARED CAUSE?
Rather than hand-pick one trait + one climate axis, this scans the ENTIRE coded
convergence battery against EVERY available environmental axis, in both directions,
and reports the whole matrix under multiple-testing control. It is the general
"environment lens" applied across all of mammalian convergence at once.

WHAT IT DOES
  * Traits: every clade-/trait-coded mammal convergence the viewer already defines
    (marine, aerial, fossorial, myrmecophagy, hopping, carnivory, diurnality, giant).
  * Environment: every per-species environmental axis in PanTHERIA that has real
    coverage -- mean temperature, mean precipitation, actual evapotranspiration
    (AET, a water-energy / productivity proxy). These are REAL measured per-species
    values, not imputed; species lacking a value are dropped from that test.
  * Test: for each (trait x axis x direction) it runs `mk.trait_environment_test` --
    is the derived trait associated with that environment MORE than if the trait had
    evolved NEUTRALLY (Mk) on the real dated tree, independent of the environment?
    The statistic is a phylogenetically-controlled phi correlation; the null is the
    trait simulated under Mk with no environmental coupling. This is the existing
    "shared cause" hook, run generally.
  * Multiple testing: Benjamini-Hochberg FDR across the whole battery (q-values), so
    "many traits tested" cannot manufacture a hit. Effect size = observed phi.
  * Sensitivity: re-run with a tertile (extremes-only) environment split to show the
    significant couplings are not an artifact of the median threshold.

HONESTY CONTRACT
  * Real PanTHERIA environmental values only; missing dropped, never imputed.
  * The trait's phylogeny is controlled for (Mk null on the real dated chronogram),
    so a hit means association BEYOND shared ancestry -- environment as a driver.
  * Every test is reported, including the null results, with FDR q-values; nothing
    is promoted past what the correction allows.
  * This is mammals-first because PanTHERIA carries the environment locally; birds
    and fish need the GBIF -> WorldClim occurrence bridge (Layer C), the documented
    next data step, before the SAME scan generalizes to them.

Run:  python evolve/env_driver.py
Out:  evolve/results/env_driver.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "viewer"))

from chronogram import (load_mammal_chronogram, all_tip_labels, prune_to,        # noqa: E402
                        normalize_unit_depth, branch_length_dict)
from polarity import _leaves                                                      # noqa: E402
from mk import trait_environment_test                                            # noqa: E402
from ingest import pantheria                                                      # noqa: E402
# reuse the EXACT trait codings + column names the mammal viewer defines
from export_mammal_tree import (MARINE_FAMS, GLIDER_FAMS, FOSSORIAL_FAMS,         # noqa: E402
                                MYRMECOPHAGY_FAMS, HOPPING_FAMS,
                                COL_ORDER, COL_FAMILY, COL_BODY, COL_TROPHIC,
                                COL_ACTIVITY, _num)

RESULTS = HERE / "results"

# environmental axes: column, label, unit note (PanTHERIA, real per-species means)
ENV_AXES = [
    ("28-2_Temp_Mean_01degC", "mean_temperature", "0.1 degC"),
    ("28-1_Precip_Mean_mm", "mean_precipitation", "mm"),
    ("30-1_AET_Mean_mm", "AET_water_energy", "mm"),
]

MIN_YES = 5      # need at least this many derived-state species for a meaningful test
N_SIM = 499
N_SIM_SENS = 299


def _norm(b):
    return str(b).replace(" ", "_")


def trait_codings(df) -> dict[str, dict[str, str]]:
    """{trait_key: {underscore_binom: 'yes'/'no'}} using the viewer's definitions.
    Clade-coded traits score every species; column-coded ones skip missing values."""
    body_log = []
    for b, r in df.iterrows():
        v = _num(r.get(COL_BODY))
        if v and v > 0:
            import math
            body_log.append(math.log10(v))
    body_log.sort()
    giant_cut = body_log[int(0.9 * len(body_log))] if body_log else float("inf")

    def fam(r):
        return str(r.get(COL_FAMILY) or "")

    def order(r):
        return str(r.get(COL_ORDER) or "")

    import math
    defs = {
        "marine": lambda r: "yes" if fam(r) in MARINE_FAMS else "no",
        "aerial": lambda r: "yes" if (fam(r) in GLIDER_FAMS or order(r) == "Chiroptera") else "no",
        "fossorial": lambda r: "yes" if fam(r) in FOSSORIAL_FAMS else "no",
        "myrmecophagy": lambda r: "yes" if fam(r) in MYRMECOPHAGY_FAMS else "no",
        "hopping": lambda r: "yes" if fam(r) in HOPPING_FAMS else "no",
        "carnivory": lambda r: ("yes" if (x := _num(r.get(COL_TROPHIC))) == 3
                                else ("no" if x is not None else None)),
        "diurnality": lambda r: ("yes" if (x := _num(r.get(COL_ACTIVITY))) == 3
                                 else ("no" if x is not None else None)),
        "giant": lambda r: ("yes" if (b := _num(r.get(COL_BODY))) and b > 0 and math.log10(b) >= giant_cut
                            else ("no" if (b := _num(r.get(COL_BODY))) and b > 0 else None)),
    }
    out = {}
    for key, fn in defs.items():
        d = {}
        for b, r in df.iterrows():
            if not isinstance(b, str):
                continue
            v = fn(r)
            if v is not None:
                d[_norm(b)] = v
        out[key] = d
    return out


def env_values(df) -> dict[str, dict[str, float]]:
    """{axis_name: {underscore_binom: value}} for real (non -999) PanTHERIA values."""
    out = {}
    for col, name, _unit in ENV_AXES:
        d = {}
        for b, r in df.iterrows():
            if not isinstance(b, str):
                continue
            v = _num(r.get(col))
            if v is not None and v != -999:
                d[_norm(b)] = float(v)
        out[name] = d
    return out


def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg q-values for a list of p-values (order preserved)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        val = min(prev, pvals[i] * m / k)
        q[i] = val
        prev = val
    return q


def run_one(tree_full, trait_state, env_val, direction, threshold_mode, n_sim):
    """One (trait x axis x direction) test. direction 'high'/'low' picks which side of
    the environment split is the 'environment present'. threshold_mode 'median' or
    'tertile' (extremes only). Returns dict or None if not enough data."""
    import statistics
    shared = [t for t in trait_state if t in env_val]
    if len(shared) < 2 * MIN_YES:
        return None
    vals = sorted(env_val[t] for t in shared)
    if threshold_mode == "median":
        cut = statistics.median(vals)
        env_true = {t: (env_val[t] >= cut) if direction == "high" else (env_val[t] < cut)
                    for t in shared}
        keep = set(shared)
    else:  # tertile: keep only the extreme thirds, drop the middle
        lo = vals[len(vals) // 3]
        hi = vals[2 * len(vals) // 3]
        keep = {t for t in shared if env_val[t] <= lo or env_val[t] >= hi}
        env_true = {t: (env_val[t] >= hi) if direction == "high" else (env_val[t] <= lo)
                    for t in keep}
    # need variation in both trait and env among kept tips
    ts = {t: trait_state[t] for t in keep}
    n_yes = sum(1 for v in ts.values() if v == "yes")
    if n_yes < MIN_YES or n_yes == len(ts):
        return None
    if not (any(env_true.values()) and not all(env_true.values())):
        return None
    pruned = prune_to(tree_full, set(keep))
    if pruned is None:
        return None
    normalize_unit_depth(pruned)
    brlen = branch_length_dict(pruned)
    placed = {n.label: ts[n.label] for n in _leaves(pruned) if n.label in ts}
    et = {n.label: env_true[n.label] for n in _leaves(pruned) if n.label in env_true}
    if len(set(placed.values())) < 2:
        return None
    res = trait_environment_test(pruned, placed, et, "yes", n_sim=n_sim, seed=0, brlen=brlen)
    if "error" in res:
        return None
    res["n_tips"] = len(placed)
    res["n_derived"] = n_yes
    return res


def main():
    print("Loading dated mammal chronogram + PanTHERIA (traits + environment) ...")
    tree = load_mammal_chronogram()
    on_tree = all_tip_labels(tree)
    df = pantheria()
    traits = trait_codings(df)
    envs = env_values(df)
    # restrict trait states to species actually on the dated tree
    traits = {k: {t: v for t, v in d.items() if t in on_tree} for k, d in traits.items()}
    print(f"  {len(on_tree)} dated tips | traits: {list(traits)} | "
          f"env axes: {[n for _, n, _ in ENV_AXES]}")

    print("Scanning trait x environment (one two-sided test per pair, median split) ...")
    tests = []
    for tk, tstate in traits.items():
        for col, axis, unit in ENV_AXES:
            rh = run_one(tree, tstate, envs[axis], "high", "median", N_SIM)
            rl = run_one(tree, tstate, envs[axis], "low", "median", N_SIM)
            if rh is None and rl is None:
                continue
            # one hypothesis per (trait, axis): two-sided p ~ 2*min(one-sided high, low),
            # direction = the sign of association (which side gave the smaller p).
            p_high = rh["p_assoc"] if rh else 1.0
            p_low = rl["p_assoc"] if rl else 1.0
            if p_high <= p_low:
                direction, r = "high", rh
            else:
                direction, r = "low", rl
            p_two = min(1.0, 2.0 * min(p_high, p_low))
            tests.append({
                "trait": tk, "env_axis": axis, "direction": direction,
                "observed_phi": r["observed_phi"], "p_one_sided": r["p_assoc"],
                "p_two_sided": round(p_two, 4), "null_phi_mean": r["null_phi_mean"],
                "n_tips": r["n_tips"], "n_derived": r["n_derived"],
            })
            print(f"  {tk:<12} ~ {axis:<18} dir={direction:<4} "
                  f"phi={r['observed_phi']:+.3f} p2={p_two:.4f} "
                  f"(n={r['n_tips']}, derived={r['n_derived']})")

    # Benjamini-Hochberg FDR across the whole battery (Idea 3)
    qs = bh_fdr([t["p_two_sided"] for t in tests])
    for t, q in zip(tests, qs):
        t["q_fdr"] = round(q, 4)

    sig = [t for t in tests if t["q_fdr"] < 0.05]
    sig.sort(key=lambda t: t["q_fdr"])

    # Sensitivity (Idea 4): re-test the FDR-significant hits with the tertile split
    print("\nSensitivity: re-testing FDR-significant couplings on tertile (extremes) split ...")
    sens = []
    for t in sig:
        tstate = traits[t["trait"]]
        r = run_one(tree, tstate, envs[t["env_axis"]], t["direction"], "tertile", N_SIM_SENS)
        ok = (r is not None and r["p_assoc"] < 0.05 and
              (r["observed_phi"] > 0) == (t["observed_phi"] > 0))
        sens.append({"trait": t["trait"], "env_axis": t["env_axis"], "direction": t["direction"],
                     "tertile_phi": (r["observed_phi"] if r else None),
                     "tertile_p": (r["p_assoc"] if r else None), "holds": bool(ok)})
        print(f"  {t['trait']:<12} ~ {t['env_axis']:<18} ({t['direction']}): "
              f"tertile phi={(r['observed_phi'] if r else float('nan')):+.3f} "
              f"p={(r['p_assoc'] if r else float('nan')):.4f} -> {'HOLDS' if ok else 'weakens'}")

    out = {
        "clade": "Mammalia (Upham 2019 dated tree)",
        "environment_source": "PanTHERIA per-species means (real, not imputed)",
        "n_tests": len(tests),
        "fdr_alpha": 0.05,
        "tests": sorted(tests, key=lambda t: t["q_fdr"]),
        "significant_after_fdr": sig,
        "sensitivity_tertile": sens,
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "env_driver.json").write_text(json.dumps(out, indent=2))

    print("\n" + "=" * 74)
    print("ENVIRONMENT-AS-DRIVER SCAN  (mammals; FDR-controlled across the battery)")
    print("=" * 74)
    print(f"  {len(tests)} tests run; {len(sig)} survive BH-FDR q<0.05:")
    for t in sig:
        print(f"   {t['trait']:<12} ~ {t['env_axis']:<18} ({t['direction']:<4}) "
              f"phi={t['observed_phi']:+.3f}  p2={t['p_two_sided']:.4f}  q={t['q_fdr']:.4f}")
    if not sig:
        print("   (none) -- no trait-environment coupling exceeds the neutral-Mk null"
              " after FDR. Reported honestly, not inflated.")
    n_hold = sum(1 for s in sens if s["holds"])
    print(f"  sensitivity: {n_hold}/{len(sens)} significant couplings hold on the"
          f" extremes-only (tertile) split.")
    print(f"\n  wrote {RESULTS / 'env_driver.json'}")


if __name__ == "__main__":
    main()
