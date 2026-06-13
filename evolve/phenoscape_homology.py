#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Homology-aware convergence on REAL fish: which shared fin losses are INDEPENDENT?

Polarity (phenoscape_polarity) showed pelvic-fin loss has multiple origins. But not
every pair of fin-less taxa is a convergence: two species in the SAME clade share
ONE inherited loss (homology by descent), not independent re-evolution. This script
applies the MRCA test (polarity.homology_filter) to partition the fin-loss pairs:

    INDEPENDENT  : the pair's MRCA still had the fin -> each lost it separately = convergence
    SHARED       : the pair's MRCA had already lost it -> one event, inherited = NOT convergence

No external assumption is injected (we do NOT hand-code "pelvic fins are gnathostome").
The tree is Open Tree's; Fitch is parameter-free. We report what the topology says,
including the residual limit it CANNOT resolve.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from phenoscape_ingest import resolve_entity, taxa_with_entity  # noqa: E402
from phenoscape_polarity import clean_name, binary_state, induced_subtree  # noqa: E402
from ingest import reconcile                                    # noqa: E402
from polarity import parse_newick, homology_filter             # noqa: E402

CHAR = "pelvic fin"


def genus(name: str) -> str:
    return name.split()[0] if name else "?"


def main():
    ent = resolve_entity(CHAR)
    per = defaultdict(list)
    for _tx, tl, ph in taxa_with_entity(ent, limit=200):
        per[clean_name(tl)].append(ph)
    states = {n: binary_state(p) for n, p in per.items()}
    states = {n: s for n, s in states.items() if s and n}

    name_to_ott = reconcile(sorted(states))
    ott_to_name = {str(o): n for n, o in name_to_ott.items()}
    newick = induced_subtree(list(name_to_ott.values()))
    tree = parse_newick(newick)
    tip_state = {str(name_to_ott[n]): s for n, s in states.items() if n in name_to_ott}

    indep, shared, carriers = homology_filter(tree, tip_state, "absent")
    nm = lambda k: ott_to_name.get(k, k)

    print("=" * 78)
    print(f"  HOMOLOGY-AWARE CONVERGENCE  '{CHAR} absent'  (MRCA independence test)")
    print("=" * 78)
    print(f"taxa placed with 'absent': {len(carriers)}")
    print(f"pairs sharing the loss     : {len(indep) + len(shared)}")
    print(f"  INDEPENDENT (convergent) : {len(indep)}")
    print(f"  SHARED-by-descent (not)  : {len(shared)}")

    print("\nSample INDEPENDENT losses (true convergence -- different lineages):")
    cross = [(a, b) for a, b in indep if genus(nm(a)) != genus(nm(b))]
    for a, b in cross[:8]:
        print(f"   {nm(a)[:32]:<34} ~ {nm(b)[:32]}")

    print("\nSample SHARED-by-descent (correctly NOT convergence):")
    for a, b in shared[:8]:
        tag = "  <-- same genus" if genus(nm(a)) == genus(nm(b)) else ""
        print(f"   {nm(a)[:32]:<34} ~ {nm(b)[:32]}{tag}")

    # un-fitting check: are SHARED pairs enriched for same-genus (as they should be)?
    def same_genus_frac(pairs):
        if not pairs:
            return float("nan")
        return sum(1 for a, b in pairs if genus(nm(a)) == genus(nm(b))) / len(pairs)
    print("\n" + "-" * 78)
    print("Cross-check (not tuned -- just observed):")
    print(f"  same-genus fraction among SHARED pairs      : {same_genus_frac(shared):.2f}")
    print(f"  same-genus fraction among INDEPENDENT pairs : {same_genus_frac(indep):.2f}")
    print("  (SHARED should be richer in same-genus pairs if the MRCA test is sound.)")

    print("\nHONEST RESIDUAL LIMIT (not faked):")
    print("  The MRCA test fixes shared-BY-DESCENT rigorously. It does NOT, on")
    print("  present/absent coding alone, distinguish a true LOSS from a lineage that")
    print("  NEVER HAD the fin (primitive absence, e.g. a jawless fish). That needs a")
    print("  homology statement / proper outgroup coding (mark the character INAPPLICABLE")
    print("  basal to the structure's origin) -- external knowledge we deliberately did")
    print("  NOT inject, to avoid fitting. It is a data-coding fix, not an algorithm tweak.")


if __name__ == "__main__":
    main()
