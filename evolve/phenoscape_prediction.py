"""Real-data morphology-state prediction on Phenoscape fish characters.

This is a Phase 3 bridge beyond the Anolis horn validation. Phenoscape does not
provide explicit niche nodes, so this is not the same species->niche->trait horn.
Instead, it asks a narrower relation-completion question:

Given a taxon's other fin states, can we recover a held-out fin state from other
taxa that carry that state, especially cross-lineage donors?
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

sys.path.insert(0, str(Path(__file__).resolve().parent))

from horns_convergence import pairwise_auroc  # noqa: E402
from phenoscape_convergence import CHARS, order_of  # noqa: E402
from phenoscape_ingest import build_model  # noqa: E402


def _char_of(state: str) -> str:
    return state.split(":", 1)[0].strip()


def _jaccard(a: set[str], b: set[str]) -> float:
    den = len(a | b)
    return len(a & b) / den if den else 0.0


def _profiles(model) -> dict[str, set[str]]:
    return {
        taxon: {m.target for m in model.lens("morphology").morphisms_from(taxon)}
        for taxon in model.taxa
    }


def _candidate_pairs(profiles: dict[str, set[str]]) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    states_by_char: dict[str, set[str]] = defaultdict(set)
    for states in profiles.values():
        for state in states:
            states_by_char[_char_of(state)].add(state)

    positives: list[tuple[str, str]] = []
    negatives: list[tuple[str, str]] = []
    for taxon, states in profiles.items():
        observed_chars = {_char_of(state) for state in states}
        for char in sorted(observed_chars):
            for state in sorted(states_by_char[char]):
                pair = (taxon, state)
                if state in states:
                    positives.append(pair)
                else:
                    negatives.append(pair)
    return positives, negatives


def _frequency_scores(
    profiles: dict[str, set[str]],
    pairs: list[tuple[str, str]],
) -> list[float]:
    by_char = defaultdict(list)
    for taxon, states in profiles.items():
        for char in {_char_of(state) for state in states}:
            by_char[char].append(taxon)
    return [
        sum(1 for taxon in by_char[_char_of(state)] if state in profiles[taxon])
        / len(by_char[_char_of(state)])
        for _candidate_taxon, state in pairs
    ]


def _donor_score(model, profiles, candidate: tuple[str, str], mode: str) -> float:
    taxon, state = candidate
    taxon_order = order_of(model, taxon)
    candidate_profile = set(profiles[taxon])
    candidate_profile.discard(state)

    scores: list[float] = []
    for donor, donor_states in profiles.items():
        if donor == taxon or state not in donor_states:
            continue
        if mode in {"cross_order", "cross_order_minus_ancestry"} and order_of(model, donor) == taxon_order:
            continue
        donor_profile = set(donor_states)
        donor_profile.discard(state)
        sim = _jaccard(candidate_profile, donor_profile)
        if mode in {"minus_ancestry", "cross_order_minus_ancestry"}:
            sim -= model.similarity("ancestry", taxon, donor)
        scores.append(sim)
    return max(scores) if scores else 0.0


def _score_set(model, profiles, pairs: list[tuple[str, str]], mode: str) -> list[float]:
    return [_donor_score(model, profiles, pair, mode) for pair in pairs]


def _metrics(name: str, scores: list[float], labels: list[int]) -> tuple[str, float]:
    return name, round(pairwise_auroc(scores, labels), 4)


def _subset(
    positives: list[tuple[str, str]],
    negatives: list[tuple[str, str]],
    predicate,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    pos = [pair for pair in positives if predicate(pair)]
    neg = [pair for pair in negatives if predicate(pair)]
    return pos, neg


def _evaluate(model, profiles, positives, negatives) -> list[tuple[str, float]]:
    pairs = positives + negatives
    labels = [1] * len(positives) + [0] * len(negatives)
    results = [
        _metrics("frequency", _frequency_scores(profiles, pairs), labels),
        _metrics("nearest_profile", _score_set(model, profiles, pairs, "nearest"), labels),
        _metrics("profile_minus_ancestry", _score_set(model, profiles, pairs, "minus_ancestry"), labels),
        _metrics("cross_order_profile", _score_set(model, profiles, pairs, "cross_order"), labels),
        _metrics(
            "cross_order_minus_ancestry",
            _score_set(model, profiles, pairs, "cross_order_minus_ancestry"),
            labels,
        ),
    ]
    return results


def main() -> int:
    print("Building Phenoscape model from cached fin characters...")
    model = build_model(CHARS, per_char_limit=120, verbose=False)
    profiles = _profiles(model)
    positives, negatives = _candidate_pairs(profiles)
    absent_pos, absent_neg = _subset(
        positives,
        negatives,
        lambda pair: "absent" in pair[1].lower() or "loss" in pair[1].lower(),
    )

    print("=" * 86)
    print("PHENOSCAPE MORPHOLOGY-STATE PREDICTION")
    print("=" * 86)
    print(f"taxa: {len(model.taxa)}")
    print(f"states: {len({state for states in profiles.values() for state in states})}")
    print(f"all-state positives: {len(positives)}  negatives: {len(negatives)}")
    print(f"absence positives  : {len(absent_pos)}  negatives: {len(absent_neg)}")

    print("\n[1] All observed fin states")
    all_results = _evaluate(model, profiles, positives, negatives)
    for name, auroc in all_results:
        print(f"  {name:<28} AUROC {auroc:.4f}")

    print("\n[2] Absence/loss states only")
    absence_results = _evaluate(model, profiles, absent_pos, absent_neg)
    for name, auroc in absence_results:
        print(f"  {name:<28} AUROC {auroc:.4f}")

    best_all = max(all_results, key=lambda item: item[1])
    best_abs = max(absence_results, key=lambda item: item[1])
    print("\nHONEST READING")
    print(f"  Best all-state predictor    : {best_all[0]} ({best_all[1]:.4f})")
    print(f"  Best absence-state predictor: {best_abs[0]} ({best_abs[1]:.4f})")
    print("  This is relation completion over real Phenoscape morphology, not a")
    print("  niche-mediated horn test. It is weaker than Anolis because explicit")
    print("  environment/niche intermediates are absent.")
    print("  A useful result here is any improvement over frequency; a negative")
    print("  result would mean morphology context alone is not enough.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

