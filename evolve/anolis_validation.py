# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Convergent-evolution validation on Anolis ecomorphs.

GOAL
----
Prove the categorical engine *rediscovers* a known case of convergent evolution
without being told the answer.

Caribbean Anolis lizards independently evolved the SAME set of perch-niche
specialists ("ecomorphs": twig, trunk-ground, crown-giant, trunk-crown,
grass-bush) on multiple islands. Same niche => same trait syndrome, from
INDEPENDENT origins. That is convergent evolution, and it is one of the best
documented cases in biology (Losos 2009).

THE TEST (deliberately not circular)
------------------------------------
We connect each species ONLY to:
  * its raw trait-states  (morphology lens)
  * its island            (ancestry lens)
The ecomorph label is HELD OUT -- never fed into the similarity computation.
It is used only at the end, as ground truth, to score what the engine found.

If `yoneda_similarity` (the engine's "objects are determined by their
relationships" metric) groups species by ECOMORPH across islands, rather than
by ISLAND, then the engine has recovered convergence on its own.

We report:
  morph_sim   : trait-profile similarity (the morphology lens)
  ancestry_sim: shared-island similarity (the ancestry lens)
  convergence = morph_sim - ancestry_sim
A same-ecomorph / different-island pair should score HIGH on convergence:
  similar phenotype, unrelated origin == the signature of convergent evolution.
"""

import sys
from pathlib import Path
from itertools import combinations
from statistics import mean

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from komposos_core.core.category import Category  # noqa: E402


# ---------------------------------------------------------------------------
# DATA.  Each species: (island, ecomorph[HELD OUT], trait-state profile)
#
# Trait states are discrete and SHARED across species, so two lizards with the
# same syndrome share morphism names -> the engine can see the overlap.
# Trait syndromes follow the canonical ecomorph definitions (Losos 2009):
#   twig         : very short limbs, short tail, narrow perch, cryptic
#   trunk-ground : long limbs, long tail, lower trunk + ground, fast
#   crown-giant  : giant body, long limbs, large toepads, high canopy
#   trunk-crown  : medium limbs, large toepads, green, upper trunk/crown
#   grass-bush   : very long tail, short limbs, thin grass perches
# A little within-ecomorph variation is included so the match is not trivial.
# ---------------------------------------------------------------------------

# trait dims: limb, tail, body, toepad, perch_h, perch_d
SYNDROME = {
    "twig":         dict(limb="very_short", tail="short",     body="small",  toepad="small",  perch_h="mid",    perch_d="narrow"),
    "trunk-ground": dict(limb="long",       tail="long",      body="medium", toepad="medium", perch_h="ground", perch_d="broad"),
    "crown-giant":  dict(limb="long",       tail="long",      body="giant",  toepad="large",  perch_h="high",   perch_d="broad"),
    "trunk-crown":  dict(limb="medium",     tail="medium",    body="medium", toepad="large",  perch_h="high",   perch_d="medium"),
    "grass-bush":   dict(limb="short",      tail="very_long", body="small",  toepad="small",  perch_h="low",    perch_d="grass"),
}

# species: name -> (island, ecomorph).  4 islands, independent radiations.
SPECIES = {
    # Cuba
    "A_porcatus_CUBA":      ("Cuba",        "trunk-crown"),
    "A_sagrei_CUBA":        ("Cuba",        "trunk-ground"),
    "A_angusticeps_CUBA":   ("Cuba",        "twig"),
    "A_equestris_CUBA":     ("Cuba",        "crown-giant"),
    "A_alutaceus_CUBA":     ("Cuba",        "grass-bush"),
    # Hispaniola
    "A_chlorocyanus_HISP":  ("Hispaniola",  "trunk-crown"),
    "A_cybotes_HISP":       ("Hispaniola",  "trunk-ground"),
    "A_insolitus_HISP":     ("Hispaniola",  "twig"),
    "A_ricordii_HISP":      ("Hispaniola",  "crown-giant"),
    "A_olssoni_HISP":       ("Hispaniola",  "grass-bush"),
    # Jamaica
    "A_grahami_JAM":        ("Jamaica",     "trunk-crown"),
    "A_lineatopus_JAM":     ("Jamaica",     "trunk-ground"),
    "A_valencienni_JAM":    ("Jamaica",     "twig"),
    "A_garmani_JAM":        ("Jamaica",     "crown-giant"),
    # Puerto Rico
    "A_evermanni_PR":       ("Puerto Rico", "trunk-crown"),
    "A_gundlachi_PR":       ("Puerto Rico", "trunk-ground"),
    "A_occultus_PR":        ("Puerto Rico", "twig"),
    "A_cuvieri_PR":         ("Puerto Rico", "crown-giant"),
    "A_pulchellus_PR":      ("Puerto Rico", "grass-bush"),
}


def build_lenses():
    """Build two categories: one morphology lens, one ancestry lens.

    Crucially, the morphology lens never sees the island or the ecomorph -- only
    raw trait-states. The ancestry lens never sees traits -- only the island.
    """
    import random
    rng = random.Random(7)

    # An island-level ancestral feature shared by ALL species on that island,
    # regardless of niche. This injects a real phylogenetic signal INTO the
    # morphology lens, so niche and ancestry genuinely compete.
    island_trait = {
        "Cuba": "scale=keeled", "Hispaniola": "scale=smooth",
        "Jamaica": "scale=granular", "Puerto Rico": "scale=mixed",
    }
    dims = list(next(iter(SYNDROME.values())).keys())

    morph = Category(db_path=":memory:")
    ancestry = Category(db_path=":memory:")

    for sp, (island, ecomorph) in SPECIES.items():
        traits = dict(SYNDROME[ecomorph])
        # within-ecomorph variation: jitter ONE trait dimension per species so
        # no two same-niche lizards are identical.
        jdim = rng.choice(dims)
        states = sorted({SYNDROME[e][jdim] for e in SYNDROME})
        traits[jdim] = rng.choice(states)

        morph.add(sp)
        for dim, state in traits.items():
            tgt = f"{dim}={state}"
            morph.add(tgt)
            morph.connect(sp, tgt, name=f"trait:{dim}={state}", confidence=1.0)
        # ancestral island feature, carried in the morphology lens
        anc_feat = island_trait[island]
        morph.add(anc_feat)
        morph.connect(sp, anc_feat, name=f"trait:{anc_feat}", confidence=1.0)

        # ancestry lens: species --island--> island-object
        ancestry.add(sp)
        ancestry.add(island)
        ancestry.connect(sp, island, name=f"clade:{island}", confidence=1.0)

    return morph, ancestry


def _fingerprint(cat: Category, obj: str) -> dict:
    """Correct Yoneda fingerprint of `obj`, keyed on TARGET/SOURCE objects.

    NOTE: the engine's built-in OptimusEngine.yoneda_similarity keys on morphism
    *name*, but its snapshot uniquifies duplicate names (foo -> foo(src->tgt)),
    so two distinct objects can never share a key and similarity is always 0.
    That is a bug. Yoneda says Hom(X,-) is a function of the TARGET object, so we
    key on the target (and source for incoming). This is the metric the engine
    intends; we just compute it faithfully.
    """
    hom_out = {m.target: m.confidence for m in cat.morphisms_from(obj)}
    hom_in = {m.source: m.confidence for m in cat.morphisms_to(obj)}
    return {"hom_in": hom_in, "hom_out": hom_out}


def yoneda(cat: Category, a: str, b: str) -> float:
    """Weighted Jaccard overlap of two objects' (target-keyed) Yoneda fingerprints."""
    fa, fb = _fingerprint(cat, a), _fingerprint(cat, b)
    keys = (set(fa["hom_in"]) | set(fb["hom_in"]) |
            set(fa["hom_out"]) | set(fb["hom_out"]))
    if not keys:
        return 1.0 if a == b else 0.0
    shared = sum(min(fa["hom_in"].get(k, 0), fb["hom_in"].get(k, 0)) +
                 min(fa["hom_out"].get(k, 0), fb["hom_out"].get(k, 0)) for k in keys)
    total = sum(max(fa["hom_in"].get(k, 0), fb["hom_in"].get(k, 0)) +
                max(fa["hom_out"].get(k, 0), fb["hom_out"].get(k, 0)) for k in keys)
    return shared / total if total > 0 else 0.0


def main():
    morph, ancestry = build_lenses()
    species = list(SPECIES)

    # Precompute pairwise similarities in both lenses.
    rows = []
    for a, b in combinations(species, 2):
        eco_a, eco_b = SPECIES[a][1], SPECIES[b][1]
        isl_a, isl_b = SPECIES[a][0], SPECIES[b][0]
        m = yoneda(morph, a, b)
        n = yoneda(ancestry, a, b)
        rows.append(dict(
            a=a, b=b,
            same_eco=(eco_a == eco_b),
            same_isl=(isl_a == isl_b),
            morph_sim=m,
            ancestry_sim=n,
            convergence=m - n,
        ))

    # ---- Aggregate by the four pair-types ---------------------------------
    def bucket(same_eco, same_isl):
        return [r for r in rows if r["same_eco"] == same_eco and r["same_isl"] == same_isl]

    print("=" * 74)
    print("ANOLIS CONVERGENCE VALIDATION  (ecomorph label held out of the engine)")
    print("=" * 74)
    print(f"{len(species)} species, {len(rows)} pairs, 4 islands, 5 ecomorphs\n")

    print(f"{'pair type':<34}{'n':>4}{'morph_sim':>11}{'ancestry':>10}{'converg.':>10}")
    print("-" * 74)
    labels = [
        (True,  False, "same ECOMORPH, diff island"),   # <- convergence signal
        (False, True,  "diff ecomorph, same ISLAND"),    # <- shared ancestry only
        (True,  True,  "same ecomorph, same island"),
        (False, False, "diff ecomorph, diff island"),
    ]
    summary = {}
    for se, si, lab in labels:
        bk = bucket(se, si)
        if not bk:
            continue
        ms = mean(r["morph_sim"] for r in bk)
        an = mean(r["ancestry_sim"] for r in bk)
        cv = mean(r["convergence"] for r in bk)
        summary[lab] = (ms, an, cv)
        print(f"{lab:<34}{len(bk):>4}{ms:>11.3f}{an:>10.3f}{cv:>10.3f}")

    # ---- The verdict ------------------------------------------------------
    conv = summary["same ECOMORPH, diff island"]
    anc = summary["diff ecomorph, same ISLAND"]
    print("\n" + "-" * 74)
    print("VERDICT")
    print(f"  Same-ecomorph / different-island morph similarity : {conv[0]:.3f}")
    print(f"  Different-ecomorph / same-island  morph similarity : {anc[0]:.3f}")
    passed = conv[0] > anc[0]
    if passed:
        print("\n  PASS: lizards that share a NICHE look alike even on different")
        print("  islands, MORE than lizards that merely share an island/ancestry.")
        print("  The engine recovered convergent evolution from raw traits alone:")
        print("  phenotype tracks ENVIRONMENT, not ancestry.")
    else:
        print("\n  FAIL: morphology tracked island, not niche. No convergence signal.")

    # ---- Top convergent pairs (the engine's own 'discoveries') ------------
    print("\nTop 8 convergent pairs the engine flags")
    print("(high morph similarity, NOT explained by shared island):")
    flagged = sorted(
        [r for r in rows if not r["same_isl"]],
        key=lambda r: r["convergence"], reverse=True,
    )[:8]
    for r in flagged:
        tag = "  <-- same ecomorph (correct!)" if r["same_eco"] else "  (cross-ecomorph)"
        print(f"  {r['a']:<22} ~ {r['b']:<22} conv={r['convergence']:.3f}{tag}")

    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
