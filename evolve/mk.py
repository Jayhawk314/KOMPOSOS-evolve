#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Neutral-Mk excess-homoplasy test (parametric bootstrap).

Stronger than the randomization PTP test (which uses the no-structure extreme as
its null): here the null is a NEUTRAL continuous-time Markov (Mk2) process running
on the tree WITH branch lengths. We ask whether the observed homoplasy (parsimony
steps) EXCEEDS what plain drift at the best-fit rate produces -- i.e. whether the
state re-evolved MORE than neutral evolution alone would explain.

Careful rate selection (the discernment this needs):
  The transition rate q is a NUISANCE parameter. We estimate it by MAXIMUM
  LIKELIHOOD under the null Mk2 model (Felsenstein pruning) -- the rate that best
  explains the observed tips assuming NO convergence drive. It is not tuned toward
  a result. Then we simulate under that rate and compare. (This is statistical
  estimation; no categorical tool estimates a continuous rate -- see the project
  notes on not dressing decoration as method.)

Branch lengths: Grafen's method (parameter-free, from topology) since the Open Tree
topology has none. An empirical chronogram (TimeTree/DateLife) would be the gold
standard and is the stated next refinement.
"""

from __future__ import annotations

import math
import random
import statistics
from typing import Dict

from polarity import Node, fitch_up, ott_of


# ── Grafen branch lengths (parameter-free, topology only) ────────────────────
def grafen_lengths(tree: Node) -> Dict[int, float]:
    """Branch length of each node = parent_height - height, height = (#desc leaves - 1),
    normalised by root height. Standard substitute when no dated tree exists."""
    nleaf: Dict[int, int] = {}

    def count(n: Node) -> int:
        if n.is_leaf:
            nleaf[id(n)] = 1
            return 1
        c = sum(count(ch) for ch in n.children)
        nleaf[id(n)] = c
        return c
    count(tree)
    height = {k: (v - 1) for k, v in nleaf.items()}
    brlen: Dict[int, float] = {id(tree): 0.0}

    def setbl(n: Node):
        for ch in n.children:
            brlen[id(ch)] = max(1e-6, height[id(n)] - height[id(ch)])
            setbl(ch)
    setbl(tree)
    R = height[id(tree)] or 1
    return {k: v / R for k, v in brlen.items()}


# ── Mk2 likelihood (Felsenstein pruning) + ML rate ───────────────────────────
def _states(tip_state):
    return sorted(set(tip_state.values()))


def mk2_loglik(tree: Node, tip_state: Dict[str, str], q: float,
               brlen: Dict[int, float], s0: str) -> float:
    def L(n: Node):
        if n.is_leaf:
            st = tip_state.get(ott_of(n.label) or n.label)
            if st is None:
                return [1.0, 1.0]
            return [1.0, 0.0] if st == s0 else [0.0, 1.0]
        out = [1.0, 1.0]
        for ch in n.children:
            lc = L(ch)
            t = brlen.get(id(ch), 1.0)
            e = math.exp(-2 * q * t)
            ps, pd = 0.5 + 0.5 * e, 0.5 - 0.5 * e
            out[0] *= ps * lc[0] + pd * lc[1]
            out[1] *= pd * lc[0] + ps * lc[1]
        return out
    root = L(tree)
    tot = 0.5 * root[0] + 0.5 * root[1]
    return math.log(tot) if tot > 0 else -1e9


def mk2_mle(tree: Node, tip_state: Dict[str, str], brlen: Dict[int, float], s0: str) -> float:
    """ML estimate of the Mk2 rate q by grid + golden-section refinement."""
    grid = [10 ** e for e in [x / 4 for x in range(-8, 9)]]  # 1e-2 .. 1e2
    best = max(grid, key=lambda q: mk2_loglik(tree, tip_state, q, brlen, s0))
    lo, hi = best / 3, best * 3
    gr = (math.sqrt(5) - 1) / 2
    for _ in range(40):
        x1, x2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
        if mk2_loglik(tree, tip_state, x1, brlen, s0) > mk2_loglik(tree, tip_state, x2, brlen, s0):
            hi = x2
        else:
            lo = x1
    return (lo + hi) / 2


# ── Simulation + excess-homoplasy test ───────────────────────────────────────
def simulate_mk2(tree: Node, q: float, brlen: Dict[int, float], s0: str, s1: str,
                 rng: random.Random) -> Dict[str, str]:
    tip: Dict[str, str] = {}

    def rec(n: Node, state: str):
        if n.is_leaf:
            tip[ott_of(n.label) or n.label] = state
            return
        for ch in n.children:
            t = brlen.get(id(ch), 1.0)
            pd = 0.5 - 0.5 * math.exp(-2 * q * t)
            cs = (s1 if state == s0 else s0) if rng.random() < pd else state
            rec(ch, cs)
    rec(tree, s0 if rng.random() < 0.5 else s1)
    return tip


def excess_homoplasy_test(tree: Node, tip_state: Dict[str, str],
                          n_sim: int = 999, seed: int = 0,
                          brlen: Dict[int, float] | None = None) -> dict:
    """Parametric-bootstrap test: is observed homoplasy > neutral Mk expectation?

    `brlen` may be supplied (e.g. real dated branch lengths from a chronogram);
    if None, parameter-free Grafen lengths are computed from the topology.
    """
    sts = _states(tip_state)
    if len(sts) != 2:
        return {"error": "needs a binary character"}
    s0, s1 = sts
    if brlen is None:
        brlen = grafen_lengths(tree)
    q_hat = mk2_mle(tree, tip_state, brlen, s0)
    observed = fitch_up(tree, tip_state, set(sts))
    rng = random.Random(seed)
    null = []
    for _ in range(n_sim):
        sim = simulate_mk2(tree, q_hat, brlen, s0, s1, rng)
        null.append(fitch_up(tree, sim, set(sts)))
    ge = sum(1 for s in null if s >= observed)
    le = sum(1 for s in null if s <= observed)
    return {
        "observed_steps": observed,
        "q_hat": round(q_hat, 4),
        "null_mean": round(statistics.mean(null), 2),
        "null_sd": round(statistics.pstdev(null), 2),
        "p_excess": round((ge + 1) / (n_sim + 1), 4),    # observed MORE homoplasy than neutral
        "p_deficit": round((le + 1) / (n_sim + 1), 4),   # observed LESS (conserved)
        "n_sim": n_sim,
    }


def phi_coefficient(trait_true: Dict[str, bool], env_true: Dict[str, bool], tips):
    """2x2 phi correlation of two binary variables over shared tips (+counts)."""
    n11 = n10 = n01 = n00 = 0
    for t in tips:
        a, b = trait_true.get(t), env_true.get(t)
        if a is None or b is None:
            continue
        if a and b:
            n11 += 1
        elif a and not b:
            n10 += 1
        elif (not a) and b:
            n01 += 1
        else:
            n00 += 1
    den = math.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    phi = (n11 * n00 - n10 * n01) / den if den > 0 else 0.0
    return phi, (n11, n10, n01, n00)


def trait_environment_test(tree: Node, trait_state: Dict[str, str],
                           env_true: Dict[str, bool], derived_state: str,
                           n_sim: int = 999, seed: int = 0,
                           brlen: Dict[int, float] | None = None) -> dict:
    """Phylogenetically-controlled trait-environment association (convergence-to-env).

    Tests whether the derived trait co-occurs with the environment MORE than expected
    if the trait evolved NEUTRALLY on the tree, independent of the environment. The
    null: simulate the trait under Mk (rate ML-estimated under the no-association null)
    and recompute the trait-env correlation each time. Observed phi above the null =
    the trait converges WITH the environment beyond phylogeny (environment as driver).

    Un-fittable: Mk rate estimated under the null; environment fixed; statistic is a
    plain phi correlation. p_assoc = P(null phi >= observed) (one-sided, positive).
    """
    sts = _states(trait_state)
    if len(sts) != 2:
        return {"error": "trait must be binary"}
    s0, s1 = sts
    if brlen is None:
        brlen = grafen_lengths(tree)
    tips = [ott_of(n.label) or n.label for n in _leaves_local(tree)]
    trait_true = {t: (trait_state.get(t) == derived_state) for t in tips
                  if t in trait_state}
    obs_phi, counts = phi_coefficient(trait_true, env_true, tips)
    q_hat = mk2_mle(tree, trait_state, brlen, s0)
    rng = random.Random(seed)
    null = []
    for _ in range(n_sim):
        sim = simulate_mk2(tree, q_hat, brlen, s0, s1, rng)
        sim_true = {t: (v == derived_state) for t, v in sim.items()}
        null.append(phi_coefficient(sim_true, env_true, tips)[0])
    ge = sum(1 for p in null if p >= obs_phi)
    return {
        "observed_phi": round(obs_phi, 3),
        "counts_n11_n10_n01_n00": counts,
        "q_hat": round(q_hat, 4),
        "null_phi_mean": round(statistics.mean(null), 3),
        "p_assoc": round((ge + 1) / (n_sim + 1), 4),
        "n_sim": n_sim,
    }


def _leaves_local(node: Node):
    if node.is_leaf:
        yield node
    for c in node.children:
        yield from _leaves_local(c)


if __name__ == "__main__":
    from polarity import parse_newick
    # sanity: a perfectly clustered trait should NOT show EXCESS homoplasy
    t = parse_newick("((((A,B),C),(D,E)),((F,G),(H,I)));")
    clustered = {k: ("absent" if k in ("A", "B", "C") else "present") for k in "ABCDEFGHI"}
    r = excess_homoplasy_test(t, clustered, n_sim=499, seed=1)
    print("clustered trait:", r)
    print(f"  -> p_excess={r['p_excess']} (expect NOT significant: no excess homoplasy)")
