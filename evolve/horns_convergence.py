#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Simplicial horns over the CONVERGENCE category  +  leave-one-out retrodiction.

Ported from KOMPOSOS-IV-PHARM/oracle/horns.py and horns_retrodiction.py, with the
pharma-specific loaders removed and the construction re-pointed at convergent
evolution. Same mathematics, new domain.

The one-line point (pharma -> evolution):

    pharma:  Drug    --mech-->     Protein --assoc--> Disease   fill Drug--treats-->Disease
    ours:    Species --inhabits--> Niche   --selects--> Trait   fill Species--has-->Trait

An unfilled inner horn  Lambda^2_1 : Species --inhabits--> Niche --selects--> Trait
with NO direct Species --has--> Trait edge IS a convergence hypothesis, and FILLING
it = predicting the species evolves that trait (because its niche selects for it).
Filler confidence = product along the spine (multiplicative quantale).

Why this is a real test (no leakage): the Niche --selects--> Trait edges are built
from the OTHER species in that niche (cross-island members), excluding the held-out
species. So recovering a held-out trait is genuine convergent inference: "niche-mates
on other islands have this trait, therefore predict this species has it too."

Inner vs outer (Kan): inner-horn filling = composition (expected). Outer horns need
an edge INVERTED = an invertibility/equivalence claim; evolution is directional
(Dollo's law: lineages do not un-evolve), so the nerve should be a quasi-category
(inner-fillable) and NOT a Kan complex -- we check that below.

Run:
    python evolve/horns_convergence.py
    python evolve/horns_convergence.py --top 20
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict, namedtuple
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from komposos_core.core.category import Category  # noqa: E402
from anolis_validation import SPECIES, SYNDROME  # noqa: E402

Pair = Tuple[str, str]
Horn = namedtuple(
    "Horn",
    "a b c f_name f_conf g_name g_conf composite filled_any filled_has",
)


# ── Anolis data: species -> niche (ecomorph), species -> trait-state set ────────
# Mirrors anolis_validation's jitter so the two scripts agree.
def species_data() -> Tuple[Dict[str, str], Dict[str, Set[str]], Dict[str, Set[str]]]:
    """Returns (species_niche, species_traits, species_syndrome_traits).

    species_traits     : all trait-states a species has (jittered syndrome + an
                         ANCESTRAL island feature that the niche does NOT select).
    species_syndrome   : just the niche-relevant (convergent) trait-states -- these
                         are the positives we hold out and try to recover.
    """
    import random
    rng = random.Random(7)
    island_trait = {
        "Cuba": "scale=keeled", "Hispaniola": "scale=smooth",
        "Jamaica": "scale=granular", "Puerto Rico": "scale=mixed",
    }
    dims = list(next(iter(SYNDROME.values())).keys())

    niche: Dict[str, str] = {}
    traits: Dict[str, Set[str]] = {}
    syndrome: Dict[str, Set[str]] = {}
    for sp, (island, eco) in SPECIES.items():
        tr = dict(SYNDROME[eco])
        jdim = rng.choice(dims)                      # within-ecomorph variation
        states = sorted({SYNDROME[e][jdim] for e in SYNDROME})
        tr[jdim] = rng.choice(states)
        syn = {f"{d}={s}" for d, s in tr.items()}    # convergent (niche) traits
        niche[sp] = f"niche:{eco}"
        syndrome[sp] = syn
        traits[sp] = syn | {f"anc:{island_trait[island]}"}   # + ancestral feature
    return niche, traits, syndrome


# ── Build the nerve (typed convergence category) ────────────────────────────────
def build_nerve(
    species_niche: Dict[str, str],
    species_traits: Dict[str, Set[str]],
    skip_pair: Optional[Pair] = None,
) -> Tuple[Category, Dict[str, str]]:
    """Construct the Species/Niche/Trait category.

    Edges:
      species --inhabits--> niche                (conf 1.0)
      niche   --selects-->  trait                (conf = fraction of niche-members
                                                  that have the trait, computed from
                                                  species_traits EXCLUDING skip_pair)
      species --has-->      trait                (conf 1.0; the skip_pair edge omitted)
    Returns (category, type_by) where type_by[name] in {species, niche, trait}.
    """
    cat = Category(db_path=":memory:")
    type_by: Dict[str, str] = {}

    def add(name: str, typ: str):
        cat.add(name)
        type_by[name] = typ

    # niche membership
    members: Dict[str, List[str]] = defaultdict(list)
    for sp, nb in species_niche.items():
        add(sp, "species"); add(nb, "niche")
        cat.connect(sp, nb, name="inhabits", confidence=1.0)
        members[nb].append(sp)

    # niche --selects--> trait, support = fraction of members having it
    # (the held-out has-edge is removed from the counts -> no leakage)
    skip_sp, skip_tr = skip_pair if skip_pair else (None, None)
    niche_trait_support: Dict[Tuple[str, str], float] = {}
    for nb, mem in members.items():
        counts: Dict[str, int] = defaultdict(int)
        n = 0
        for sp in mem:
            n += 1
            for tr in species_traits[sp]:
                if sp == skip_sp and tr == skip_tr:
                    continue                          # hold out this membership
                counts[tr] += 1
        for tr, c in counts.items():
            niche_trait_support[(nb, tr)] = c / n if n else 0.0
    for (nb, tr), supp in niche_trait_support.items():
        if supp <= 0:
            continue
        add(tr, "trait")
        cat.connect(nb, tr, name="selects", confidence=float(supp))

    # species --has--> trait (ground truth), minus the held-out pair
    for sp, trs in species_traits.items():
        for tr in trs:
            if (sp, tr) == skip_pair:
                continue
            add(tr, "trait")
            cat.connect(sp, tr, name="has", confidence=1.0)

    return cat, type_by


# ── Nerve / horn core (ported, generic) ─────────────────────────────────────────
def _index(cat) -> Tuple[Dict[Pair, float], Set[Pair]]:
    """direct[(s,t)] = best confidence; has_edges = {(s,t) : a 'has' edge exists}."""
    direct: Dict[Pair, float] = {}
    has_edges: Set[Pair] = set()
    for m in cat.morphisms():
        key = (m.source, m.target)
        direct[key] = max(direct.get(key, 0.0), float(m.confidence))
        if m.name == "has":
            has_edges.add(key)
    return direct, has_edges


def inner_horns(cat, type_by: Dict[str, str],
                a_type: Optional[str] = None, c_type: Optional[str] = None) -> List[Horn]:
    """Every inner 2-horn (composable spine A->B->C with distinct A,B,C)."""
    direct, has_edges = _index(cat)
    horns: List[Horn] = []
    for B in cat.objects():
        b = B.name
        ins, outs = cat.morphisms_to(b), cat.morphisms_from(b)
        if not ins or not outs:
            continue
        for f in ins:
            a = f.source
            if a_type and type_by.get(a) != a_type:
                continue
            for g in outs:
                c = g.target
                if c_type and type_by.get(c) != c_type:
                    continue
                if len({a, b, c}) != 3:
                    continue
                horns.append(Horn(
                    a, b, c, f.name, float(f.confidence), g.name, float(g.confidence),
                    round(float(f.confidence) * float(g.confidence), 4),
                    (a, c) in direct, (a, c) in has_edges,
                ))
    return horns


def best_fillers(horns: List[Horn]) -> Dict[Pair, Horn]:
    best: Dict[Pair, Horn] = {}
    for h in horns:
        key = (h.a, h.c)
        if key not in best or h.composite > best[key].composite:
            best[key] = h
    return best


def invertible_pairs(cat) -> List[Pair]:
    """Pairs (s,t) with BOTH s->t and t->s -> outer horns fillable (Kan)."""
    direct, _ = _index(cat)
    out, seen = [], set()
    for (s, t) in direct:
        if s < t and (t, s) in direct and (s, t) not in seen:
            seen.add((s, t)); out.append((s, t))
    return out


def coherence_conflicts(horns: List[Horn], hi: float = 0.6, lo: float = 0.25) -> List[dict]:
    """Endpoint pairs reached by several niches whose fillers DISAGREE strongly."""
    by_pair: Dict[Pair, List[Horn]] = defaultdict(list)
    for h in horns:
        by_pair[(h.a, h.c)].append(h)
    conflicts = []
    for (a, c), hs in by_pair.items():
        comps = [h.composite for h in hs]
        if len(hs) >= 2 and max(comps) >= hi and min(comps) <= lo:
            strong, weak = max(hs, key=lambda h: h.composite), min(hs, key=lambda h: h.composite)
            conflicts.append({"a": a, "c": c, "n": len({h.b for h in hs}),
                              "spread": round(max(comps) - min(comps), 4),
                              "strong": (strong.b, strong.composite),
                              "weak": (weak.b, weak.composite)})
    conflicts.sort(key=lambda d: d["spread"], reverse=True)
    return conflicts


# ── Scoring: horn-max + corroboration (noisy-OR) ────────────────────────────────
def horn_pair_scores(cat, type_by) -> Tuple[Dict[Pair, float], Dict[Pair, float], Dict[Pair, int]]:
    """(max, noisy_or, support) over Species->Trait pairs, intermediates = niches."""
    by_pair_b: Dict[Pair, Dict[str, float]] = defaultdict(dict)
    for h in inner_horns(cat, type_by, a_type="species", c_type="trait"):
        key = (h.a, h.c)
        if h.composite > by_pair_b[key].get(h.b, 0.0):
            by_pair_b[key][h.b] = h.composite
    mx, nor, support = {}, {}, {}
    for key, b2c in by_pair_b.items():
        mx[key] = max(b2c.values())
        prod = 1.0
        for comp in b2c.values():
            prod *= (1.0 - comp)
        nor[key] = 1.0 - prod
        support[key] = len(b2c)
    return mx, nor, support


# ── Metrics (standalone; no pharma benchmark import) ─────────────────────────────
def pairwise_auroc(scores: List[float], labels: List[int]) -> float:
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return float("nan")
    conc = tie = 0
    for p in pos:
        for n in neg:
            if p > n:
                conc += 1
            elif p == n:
                tie += 1
    return (conc + 0.5 * tie) / (len(pos) * len(neg))


def hits_at_k(scores: List[float], labels: List[int], k: int) -> float:
    order = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    topk = set(order[:k])
    npos = sum(labels)
    return sum(1 for i in topk if labels[i] == 1) / npos if npos else 0.0


def mrr(scores: List[float], labels: List[int]) -> float:
    order = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    rr = [1.0 / (rank + 1) for rank, i in enumerate(order) if labels[i] == 1]
    return sum(rr) / len(rr) if rr else 0.0


def auprc(scores: List[float], labels: List[int]) -> float:
    order = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    tp = fp = 0
    npos = sum(labels)
    if not npos:
        return float("nan")
    ap = 0.0
    prev_recall = 0.0
    for i in order:
        if labels[i] == 1:
            tp += 1
        else:
            fp += 1
        recall = tp / npos
        precision = tp / (tp + fp)
        ap += precision * (recall - prev_recall)
        prev_recall = recall
    return ap


# ── Leave-one-(species,trait)-out retrodiction ──────────────────────────────────
def retrodict() -> dict:
    niche, traits, syndrome = species_data()
    all_species = list(niche)
    all_traits = sorted({t for ts in traits.values() for t in ts if not t.startswith("anc:")})

    # positives = each (species, syndrome-trait); negatives = trait-states a species lacks
    positives = [(sp, tr) for sp in all_species for tr in syndrome[sp]]
    negatives = [(sp, tr) for sp in all_species for tr in all_traits
                 if tr not in traits[sp]]

    n_folds = len(positives)
    held_max, held_nor, support_hist = [], [], []
    neg_max: Dict[Pair, float] = {n: 0.0 for n in negatives}
    neg_nor: Dict[Pair, float] = {n: 0.0 for n in negatives}

    for held in positives:
        cat, type_by = build_nerve(niche, traits, skip_pair=held)
        mx, nor, support = horn_pair_scores(cat, type_by)
        held_max.append(mx.get(held, 0.0))
        held_nor.append(nor.get(held, 0.0))
        support_hist.append(support.get(held, 0))
        for ng in negatives:
            neg_max[ng] += mx.get(ng, 0.0)
            neg_nor[ng] += nor.get(ng, 0.0)

    def metrics(held_scores, neg_sums):
        scores = list(held_scores) + [v / n_folds for v in neg_sums.values()]
        labels = [1] * len(held_scores) + [0] * len(neg_sums)
        return {
            "auroc": round(pairwise_auroc(scores, labels), 4),
            "auprc": round(auprc(scores, labels), 4),
            "hits_at_10": round(hits_at_k(scores, labels, 10), 3),
            "mrr": round(mrr(scores, labels), 4),
            "held_scored": sum(1 for s in held_scores if s > 0),
        }

    # baseline: global trait frequency (how often any species has the trait)
    freq: Dict[str, float] = {}
    for tr in all_traits:
        freq[tr] = sum(1 for sp in all_species if tr in traits[sp]) / len(all_species)
    base_scores = [freq.get(tr, 0.0) for (_, tr) in positives] + \
                  [freq.get(tr, 0.0) for (_, tr) in negatives]
    base_labels = [1] * len(positives) + [0] * len(negatives)

    return {
        "n_species": len(all_species), "n_trait_states": len(all_traits),
        "n_positives": n_folds, "n_negatives": len(negatives),
        "avg_intermediates_per_held": round(sum(support_hist) / len(support_hist), 2),
        "horn_max": metrics(held_max, neg_max),
        "horn_noisy_or": metrics(held_nor, neg_nor),
        "baseline_freq_auroc": round(pairwise_auroc(base_scores, base_labels), 4),
    }


# ── Structural diagnostics on the full nerve ────────────────────────────────────
def diagnostics(top: int) -> None:
    niche, traits, _ = species_data()
    cat, type_by = build_nerve(niche, traits)
    n_obj, n_mor = len(list(cat.objects())), len(list(cat.morphisms()))
    horns = inner_horns(cat, type_by, a_type="species", c_type="trait")
    best = best_fillers(horns)
    filled = [h for h in best.values() if h.filled_has]
    unfilled = sorted([h for h in best.values() if not h.filled_has],
                      key=lambda h: h.composite, reverse=True)

    print("=" * 76)
    print("  HORNS over the CONVERGENCE category  (nerve, dim <= 2)")
    print("=" * 76)
    print(f"0-simplices (objects)   : {n_obj}")
    print(f"1-simplices (morphisms) : {n_mor}")
    print(f"Inner 2-horns Species->Niche->Trait : {len(horns)}")
    print(f"  distinct Species->Trait pairs spanned : {len(best)}")
    print(f"  FILLED (species really has the trait) : {len(filled)}")
    print(f"  UNFILLED (predicted-but-absent)       : {len(unfilled)}")

    print("\n[1] FILLED inner horn  (species HAS the trait its niche selects)")
    for h in sorted(filled, key=lambda h: h.composite, reverse=True)[:3]:
        print(f"    {h.a} --has--> {h.c.split('=')[-1]:<10} "
              f"via {h.b.replace('niche:','')}  (composite {h.composite:.2f})")

    print(f"\n[2] UNFILLED inner horns -> FILLING = CONVERGENCE PREDICTION (top {top})")
    print("    niche selects the trait, but this species lacks it (gap / prediction):")
    for i, h in enumerate(unfilled[:top], 1):
        print(f"   {i:>2}. {h.a} --has?--> {h.c:<18} conf~{h.composite:.2f} "
              f"(niche {h.b.replace('niche:','')})")

    inv = invertible_pairs(cat)
    print("\n[3] OUTER horns / Kan condition")
    print(f"    Invertible edge pairs (s<->t): {len(inv)}")
    print("    Evolution is directional (Dollo): the nerve fills INNER horns")
    print("    (quasi-category) but not outer ones -> NOT a Kan complex. As it should be.")

    conflicts = coherence_conflicts(horns)
    print("\n[4] COHERENCE check (one trait reached via niches that DISAGREE)")
    if not conflicts:
        print("    No contradictory fillers (each trait routed through one niche).")
    for d in conflicts[:5]:
        print(f"    {d['a']} -> {d['c']}: {d['n']} niches, spread {d['spread']:.2f}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args(argv)

    diagnostics(args.top)

    print("\n" + "=" * 76)
    print("  RETRODICTION  (leave-one-(species,trait)-out; recover deleted truth)")
    print("=" * 76)
    r = retrodict()
    print(f"{r['n_species']} species x {r['n_trait_states']} trait-states")
    print(f"Held-out positives (folds): {r['n_positives']}   negatives: {r['n_negatives']}")
    print(f"Avg niches per held positive: {r['avg_intermediates_per_held']} "
          f"(=1 -> corroboration degenerate, expected)")

    def row(name, m):
        return (f"  {name:16s} AUROC {m['auroc']:.4f}  AUPRC {m['auprc']:.4f}  "
                f"Hits@10 {m['hits_at_10']:.3f}  MRR {m['mrr']:.4f}  "
                f"(scored {m['held_scored']}/{r['n_positives']})")
    print("\nHorn scorers:")
    print(row("horn-max", r["horn_max"]))
    print(row("horn-noisy_or", r["horn_noisy_or"]))
    print(f"\n  baseline (global trait frequency) AUROC : {r['baseline_freq_auroc']:.4f}")
    base = r["horn_max"]["auroc"]
    gain = base - r["baseline_freq_auroc"]
    print(f"  horn-max beats frequency baseline by    : {gain:+.4f}")
    print("\nReading: high horn-max AUROC means a held-out trait is recovered from its")
    print("niche-mates ON OTHER ISLANDS -- convergent inference, not memorization.")
    print("noisy_or == max here because each species sits in ONE niche (one intermediate);")
    print("the pharma finding that corroboration needs MANY intermediates transfers exactly.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
