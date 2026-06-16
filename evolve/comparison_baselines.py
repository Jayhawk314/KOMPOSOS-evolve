"""Conventional baselines for the thesis convergence claims.

This script does not replace the categorical/relational analyses. It asks a
smaller question: do the thesis methods beat simple baselines that a skeptical
reader would expect to see?
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path
from statistics import mean

sys.path.insert(0, str(Path(__file__).resolve().parent))

from anolis_validation import SPECIES, build_lenses, yoneda  # noqa: E402
from horns_convergence import pairwise_auroc, retrodict, species_data  # noqa: E402


def _jaccard(a: set[str], b: set[str]) -> float:
    den = len(a | b)
    return len(a & b) / den if den else 0.0


def anolis_detection_baselines() -> dict[str, float]:
    """Compare pairwise convergence scores to simple Anolis baselines."""
    morph, ancestry = build_lenses()
    rows = []
    for a, b in combinations(SPECIES, 2):
        same_eco = SPECIES[a][1] == SPECIES[b][1]
        same_island = SPECIES[a][0] == SPECIES[b][0]
        morph_sim = yoneda(morph, a, b)
        ancestry_sim = yoneda(ancestry, a, b)
        rows.append(
            {
                "label": 1 if same_eco and not same_island else 0,
                "morph": morph_sim,
                "ancestry": ancestry_sim,
                "convergence": morph_sim - ancestry_sim,
                "not_same_island": 1.0 if not same_island else 0.0,
            }
        )

    labels = [r["label"] for r in rows]
    return {
        "n_pairs": len(rows),
        "n_positive_convergent_pairs": sum(labels),
        "morph_auroc": round(pairwise_auroc([r["morph"] for r in rows], labels), 4),
        "ancestry_auroc": round(pairwise_auroc([r["ancestry"] for r in rows], labels), 4),
        "not_same_island_auroc": round(
            pairwise_auroc([r["not_same_island"] for r in rows], labels), 4
        ),
        "convergence_auroc": round(
            pairwise_auroc([r["convergence"] for r in rows], labels), 4
        ),
        "mean_convergence_positive": round(
            mean(r["convergence"] for r in rows if r["label"] == 1), 4
        ),
        "mean_convergence_negative": round(
            mean(r["convergence"] for r in rows if r["label"] == 0), 4
        ),
    }


def _nearest_neighbor_score(
    candidate_species: str,
    candidate_trait: str,
    profiles: dict[str, set[str]],
    true_traits: dict[str, set[str]],
) -> float:
    """Best profile similarity to another species that has candidate_trait."""
    donors = [
        sp for sp, traits in true_traits.items()
        if sp != candidate_species and candidate_trait in traits
    ]
    if not donors:
        return 0.0
    profile = profiles[candidate_species]
    return max(_jaccard(profile, profiles[donor]) for donor in donors)


def horn_nearest_neighbor_baseline() -> dict[str, float]:
    """Leave-one-out nearest-neighbor baseline over the same folds as horns."""
    _, traits, syndrome = species_data()
    all_species = list(traits)
    all_traits = sorted({t for ts in traits.values() for t in ts if not t.startswith("anc:")})
    positives = [(sp, tr) for sp in all_species for tr in syndrome[sp]]
    negatives = [
        (sp, tr) for sp in all_species for tr in all_traits
        if tr not in traits[sp]
    ]

    positive_scores: list[float] = []
    negative_sums = {pair: 0.0 for pair in negatives}

    for held_sp, held_trait in positives:
        profiles = {
            sp: {t for t in ts if not t.startswith("anc:")}
            for sp, ts in traits.items()
        }
        profiles[held_sp] = set(profiles[held_sp])
        profiles[held_sp].discard(held_trait)
        positive_scores.append(
            _nearest_neighbor_score(held_sp, held_trait, profiles, syndrome)
        )
        for pair in negatives:
            negative_sums[pair] += _nearest_neighbor_score(pair[0], pair[1], profiles, syndrome)

    negative_scores = [value / len(positives) for value in negative_sums.values()]
    labels = [1] * len(positive_scores) + [0] * len(negative_scores)
    scores = positive_scores + negative_scores
    return {
        "n_positives": len(positive_scores),
        "n_negatives": len(negative_scores),
        "nearest_neighbor_auroc": round(pairwise_auroc(scores, labels), 4),
        "mean_positive_score": round(mean(positive_scores), 4),
        "mean_negative_score": round(mean(negative_scores), 4),
    }


def main() -> int:
    anolis = anolis_detection_baselines()
    horn_nn = horn_nearest_neighbor_baseline()
    horn = retrodict()

    print("=" * 78)
    print("COMPARISON BASELINES FOR THE THESIS CLAIMS")
    print("=" * 78)

    print("\n[1] Anolis detection: known convergent pair label")
    print(f"  pairs: {anolis['n_pairs']}  positives: {anolis['n_positive_convergent_pairs']}")
    print(f"  morph-only AUROC              : {anolis['morph_auroc']:.4f}")
    print(f"  ancestry-only AUROC           : {anolis['ancestry_auroc']:.4f}")
    print(f"  not-same-island AUROC         : {anolis['not_same_island_auroc']:.4f}")
    print(f"  morph-minus-ancestry AUROC    : {anolis['convergence_auroc']:.4f}")
    print(f"  mean convergence, positives   : {anolis['mean_convergence_positive']:.4f}")
    print(f"  mean convergence, negatives   : {anolis['mean_convergence_negative']:.4f}")

    print("\n[2] Horn retrodiction: nearest-neighbor baseline")
    print(f"  positives: {horn_nn['n_positives']}  negatives: {horn_nn['n_negatives']}")
    print(f"  nearest-neighbor AUROC        : {horn_nn['nearest_neighbor_auroc']:.4f}")
    print(f"  trait-frequency AUROC         : {horn['baseline_freq_auroc']:.4f}")
    print(f"  horn-max AUROC                : {horn['horn_max']['auroc']:.4f}")
    print(
        "  horn gain vs nearest-neighbor : "
        f"{horn['horn_max']['auroc'] - horn_nn['nearest_neighbor_auroc']:+.4f}"
    )
    print(
        "  horn gain vs trait-frequency  : "
        f"{horn['horn_max']['auroc'] - horn['baseline_freq_auroc']:+.4f}"
    )

    print("\nHONEST READING")
    print("  - The Anolis detection baseline asks whether subtracting ancestry improves")
    print("    a simple morphometric similarity score for known cross-island ecomorphs.")
    print("  - The nearest-neighbor retrodiction baseline is intentionally strong: it")
    print("    uses held-out trait profiles, but no explicit species->niche->trait horn.")
    print("  - If the horn score wins, the niche-mediated relational path is adding")
    print("    information beyond raw profile proximity.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

