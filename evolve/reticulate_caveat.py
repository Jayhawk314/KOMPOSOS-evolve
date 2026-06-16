"""Reticulate-history caveat for tree-only convergence inference.

The convergence gauntlet uses tree-based tools: Fitch polarity, MRCA homology,
randomization, and Mk nulls. Those are appropriate when the history is well
approximated by a bifurcating tree. They are not sufficient when the relevant
trait history includes hybridization, introgression, symbiosis, or horizontal
transfer.

This script demonstrates the limit with two histories that produce the same
observed tree and tip states:

1. true convergence: two independent origins of the derived state;
2. reticulate transfer: one origin plus a transfer edge to a distant lineage.

The tree-only pipeline must return the same answer for both, because it is given
the same tree and the same tips. That is the caveat the dissertation should make
explicit.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mk import excess_homoplasy_test  # noqa: E402
from polarity import analyze, homology_filter, parse_newick, randomization_test  # noqa: E402


TREE = "(((A1,A2),(B1,B2)),((C1,C2),(D1,D2)));"


@dataclass(frozen=True)
class Scenario:
    name: str
    biological_history: str
    tip_state: dict[str, str]
    innovations: int
    transfers: int


SCENARIOS = [
    Scenario(
        name="true_convergence",
        biological_history="derived arises independently on A1 and C1",
        tip_state={
            "A1": "derived", "A2": "ancestral",
            "B1": "ancestral", "B2": "ancestral",
            "C1": "derived", "C2": "ancestral",
            "D1": "ancestral", "D2": "ancestral",
        },
        innovations=2,
        transfers=0,
    ),
    Scenario(
        name="reticulate_transfer",
        biological_history="derived arises on A1, then transfers/introgresses to C1",
        tip_state={
            "A1": "derived", "A2": "ancestral",
            "B1": "ancestral", "B2": "ancestral",
            "C1": "derived", "C2": "ancestral",
            "D1": "ancestral", "D2": "ancestral",
        },
        innovations=1,
        transfers=1,
    ),
    Scenario(
        name="introgressed_cluster",
        biological_history="derived is inherited within A, then transfers/introgresses to C1",
        tip_state={
            "A1": "derived", "A2": "derived",
            "B1": "ancestral", "B2": "ancestral",
            "C1": "derived", "C2": "ancestral",
            "D1": "ancestral", "D2": "ancestral",
        },
        innovations=1,
        transfers=1,
    ),
]


def _tree_only_verdict(scenario: Scenario) -> dict[str, object]:
    tree = parse_newick(TREE)
    pol = analyze("reticulate_trait", tree, scenario.tip_state)
    independent, shared, carriers = homology_filter(tree, scenario.tip_state, "derived")
    ptp = randomization_test(parse_newick(TREE), scenario.tip_state, n_perm=499, seed=11)
    mk = excess_homoplasy_test(parse_newick(TREE), scenario.tip_state, n_sim=299, seed=11)
    return {
        "steps": pol.steps,
        "root": pol.root_state,
        "origins_derived": pol.origins.get("derived", 0),
        "tree_calls_convergent": pol.is_convergent,
        "derived_carriers": carriers,
        "mrca_independent_pairs": len(independent),
        "mrca_shared_pairs": len(shared),
        "ptp_p_fewer": ptp["p_fewer"],
        "mk_p_excess": mk["p_excess"],
        "mk_null_mean": mk["null_mean"],
    }


def main() -> int:
    print("=" * 92)
    print("RETICULATE-HISTORY CAVEAT: TREE-ONLY CONVERGENCE CAN OVERSTATE INDEPENDENCE")
    print("=" * 92)
    print(f"observed tree: {TREE}")
    print("\nThe first two scenarios have identical tip states. They differ only in the")
    print("unobserved history: independent origin versus transfer/introgression.\n")

    rows = []
    for scenario in SCENARIOS:
        verdict = _tree_only_verdict(scenario)
        rows.append((scenario, verdict))

    print(
        "scenario              bio innov transfer  tree steps origins conv "
        "indepPairs sharedPairs PTP p_fewer Mk p_excess"
    )
    print("-" * 92)
    for scenario, verdict in rows:
        print(
            f"{scenario.name:<21}"
            f"{scenario.innovations:>5}"
            f"{scenario.transfers:>9}"
            f"{verdict['steps']:>11}"
            f"{verdict['origins_derived']:>8}"
            f"{str(verdict['tree_calls_convergent']):>6}"
            f"{verdict['mrca_independent_pairs']:>11}"
            f"{verdict['mrca_shared_pairs']:>12}"
            f"{verdict['ptp_p_fewer']:>12}"
            f"{verdict['mk_p_excess']:>12}"
        )

    print("\nDETAIL")
    for scenario, verdict in rows:
        print(f"  {scenario.name}:")
        print(f"    biological history : {scenario.biological_history}")
        print(f"    derived carriers   : {', '.join(verdict['derived_carriers'])}")
        print(
            "    tree-only reading  : "
            f"{verdict['origins_derived']} derived origins, "
            f"convergent={verdict['tree_calls_convergent']}"
        )

    first = rows[0][1]
    second = rows[1][1]
    identical_tree_output = first == second

    print("\nHONEST READING")
    print(
        "  True convergence and reticulate transfer produce identical tree-only "
        f"outputs here: {identical_tree_output}."
    )
    print("  The MRCA/Fitch/Mk layer is doing what it can with the information it is")
    print("  given. The missing information is the non-tree edge.")
    print("  Therefore the dissertation should state: tree-based convergence claims are")
    print("  conditional on the assumed species tree being the relevant inheritance")
    print("  history for the trait. Where transfer, hybridization, introgression, or")
    print("  symbiosis are plausible, the future extension is a phylogenetic network or")
    print("  an explicit reticulate event model.")

    return 0 if identical_tree_output else 1


if __name__ == "__main__":
    raise SystemExit(main())

