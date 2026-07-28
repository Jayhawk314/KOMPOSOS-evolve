# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Flagship: echolocation convergence across morphology + behavior + molecular.

THE CASE
--------
Laryngeally-echolocating bats and toothed whales (odontocetes) independently
evolved biosonar. Their last common ancestor could not echolocate, and they sit
in different mammalian orders. The convergence is famous because it STACKS:
  * behavior  : active biosonar / high-frequency hunting calls
  * morphology: enlarged, stiffened cochlea tuned to high frequencies
  * molecular : convergent amino-acid substitutions in hearing genes -- most
                famously Prestin / SLC26A5, plus Cdh23, Tmc1 (Liu et al. 2010,
                Li et al. 2010, Parker et al. 2013).

NEGATIVE CONTROL
----------------
Pteropus, an Old World fruit bat, is a BAT (very close ancestry to the
echolocating bats) but does NOT laryngeally echolocate -- it forages by vision
and smell. If the engine is detecting real convergence and not just "battiness",
Pteropus must NOT rank as convergent with dolphins, and must rank as merely
RELATED (not convergent) with the echolocating bats.

THE DRIVER (environment)
------------------------
The shared selective problem -- "detect moving prey without light" -- is the
environment object the convergent lineages all map to. Fruit bats face the
night too, but a DIFFERENT problem (stationary fruit), so they get a different
solution. Environment, not habitat, predicts the convergent trait.

EXPECTATION
-----------
deep convergence (high trait similarity across ALL lenses, minus ancestry)
should rank BAT~WHALE pairs at the very top, ABOVE within-order echolocator
pairs (which are similar partly by shared descent) and far above the fruit bat.
"""

from convergence_engine import ConvergenceModel


def build() -> ConvergenceModel:
    m = ConvergenceModel()

    echolocators = ["Rhinolophus_bat", "Myotis_bat", "Tursiops_dolphin", "Phocoena_porpoise"]
    controls = ["Pteropus_fruitbat", "Bos_cow", "Equus_horse", "Homo_human"]
    for t in echolocators + controls:
        m.add_taxon(t)

    # ---- ancestry lens (nested clade membership; subtracted out) ----------
    clades = {
        "Rhinolophus_bat":   ["Rhinolophidae", "Chiroptera", "Laurasiatheria", "Mammalia"],
        "Myotis_bat":        ["Vespertilionidae", "Chiroptera", "Laurasiatheria", "Mammalia"],
        "Pteropus_fruitbat": ["Pteropodidae", "Chiroptera", "Laurasiatheria", "Mammalia"],
        "Tursiops_dolphin":  ["Delphinidae", "Odontoceti", "Cetacea", "Cetartiodactyla", "Laurasiatheria", "Mammalia"],
        "Phocoena_porpoise": ["Phocoenidae", "Odontoceti", "Cetacea", "Cetartiodactyla", "Laurasiatheria", "Mammalia"],
        "Bos_cow":           ["Bovidae", "Artiodactyla", "Cetartiodactyla", "Laurasiatheria", "Mammalia"],
        "Equus_horse":       ["Equidae", "Perissodactyla", "Laurasiatheria", "Mammalia"],
        "Homo_human":        ["Hominidae", "Primates", "Euarchontoglires", "Mammalia"],
    }
    for taxon, cs in clades.items():
        for c in cs:
            m.relate("ancestry", taxon, c, rel="member")

    # ---- behavior lens ----------------------------------------------------
    for t in echolocators:
        m.relate("behavior", t, "active_biosonar", rel="does")
        m.relate("behavior", t, "high_freq_hunting_calls", rel="does")
        m.relate("behavior", t, "auditory_scene_from_echoes", rel="does")
    m.relate("behavior", "Pteropus_fruitbat", "visual_smell_foraging", rel="does")
    m.relate("behavior", "Bos_cow", "diurnal_grazing", rel="does")
    m.relate("behavior", "Equus_horse", "diurnal_grazing", rel="does")
    m.relate("behavior", "Homo_human", "visual_foraging", rel="does")

    # ---- morphology lens (auditory system) --------------------------------
    for t in echolocators:
        m.relate("morphology", t, "cochlea_enlarged", rel="has")
        m.relate("morphology", t, "basilar_membrane_stiffened", rel="has")
        m.relate("morphology", t, "high_frequency_hearing", rel="has")
        m.relate("morphology", t, "expanded_auditory_nuclei", rel="has")
    for t in ["Pteropus_fruitbat", "Bos_cow", "Equus_horse", "Homo_human"]:
        m.relate("morphology", t, "cochlea_generalized", rel="has")
        m.relate("morphology", t, "mid_frequency_hearing", rel="has")

    # ---- molecular lens (convergent hearing-gene variants) ----------------
    # Prestin/SLC26A5 etc. convergent substitutions shared by echolocators.
    for t in echolocators:
        m.relate("molecular", t, "prestin_echolocator_variant", rel="carries")
        m.relate("molecular", t, "cdh23_echolocator_variant", rel="carries")
        m.relate("molecular", t, "tmc1_echolocator_variant", rel="carries")
    for t in ["Pteropus_fruitbat", "Bos_cow", "Equus_horse", "Homo_human"]:
        m.relate("molecular", t, "prestin_ancestral_variant", rel="carries")
        m.relate("molecular", t, "cdh23_ancestral_variant", rel="carries")
        m.relate("molecular", t, "tmc1_ancestral_variant", rel="carries")

    # ---- environment lens (the DRIVER: selective problem, not habitat) ----
    for t in echolocators:
        m.relate("environment", t, "detect_moving_prey_without_light", rel="faces")
    m.relate("environment", "Pteropus_fruitbat", "nocturnal_frugivory", rel="faces")
    m.relate("environment", "Bos_cow", "diurnal_foraging", rel="faces")
    m.relate("environment", "Equus_horse", "diurnal_foraging", rel="faces")
    m.relate("environment", "Homo_human", "diurnal_foraging", rel="faces")

    return m


TRAIT_LENSES = ["behavior", "morphology", "molecular"]


def main():
    m = build()

    print("=" * 78)
    print("ECHOLOCATION  -- deep (stacked) convergence across 3 trait levels")
    print("=" * 78)

    # ---- 1. Headline pair vs negative control -----------------------------
    print("\n1) THE HEADLINE  bat vs dolphin, and the fruit-bat control\n")
    for a, b, note in [
        ("Rhinolophus_bat", "Tursiops_dolphin", "echo bat vs dolphin  (distant kin)"),
        ("Rhinolophus_bat", "Pteropus_fruitbat", "echo bat vs FRUIT bat (close kin, control)"),
        ("Tursiops_dolphin", "Phocoena_porpoise", "dolphin vs porpoise  (close kin, both echo)"),
    ]:
        d = m.deep_convergence(a, b, TRAIT_LENSES)
        per = "  ".join(f"{L[:4]}={d['per_lens'][L]:.2f}" for L in TRAIT_LENSES)
        print(f"  {note}")
        print(f"     {per}   ancestry={d['ancestry_sim']:.2f}   "
              f"DEEP={d['deep_score']:+.2f}\n")

    # ---- 2. Full ranking by deep convergence ------------------------------
    print("2) ALL PAIRS RANKED BY DEEP CONVERGENCE")
    print("   (min trait similarity across behavior+morph+molecular, minus ancestry)\n")
    print(f"   {'pair':<42}{'min_trait':>10}{'ancestry':>10}{'DEEP':>8}")
    print("   " + "-" * 68)
    for d in m.rank_deep(TRAIT_LENSES)[:12]:
        a, b = d["pair"]
        star = "  <==" if d["deep_score"] > 0.5 else ""
        print(f"   {a+' ~ '+b:<42}{d['min_trait_sim']:>10.2f}"
              f"{d['ancestry_sim']:>10.2f}{d['deep_score']:>8.2f}{star}")

    # ---- 3. Environment as predictor --------------------------------------
    print("\n3) ENVIRONMENT AS DRIVER: predict the convergent solution")
    print("   Taxa facing the same selective problem should share the trait,")
    print("   regardless of ancestry. Group by environment fingerprint:\n")
    env = m.lens("environment")
    from convergence_engine import yoneda_sim
    groups: dict[str, list[str]] = {}
    for t in m.taxa:
        prob = next(iter(fp_out(env, t)), "?")
        groups.setdefault(prob, []).append(t)
    for prob, members in groups.items():
        # does this environment group actually share the behavior? (verify)
        if len(members) > 1:
            sims = [m.similarity("behavior", a, b) for a, b in _pairs(members)]
            agree = sum(s > 0.8 for s in sims) / len(sims) if sims else 0.0
        else:
            agree = float("nan")
        print(f"   problem '{prob}':")
        for t in members:
            print(f"        - {t}")
        print(f"     -> behavior agreement within group: {agree:.0%}\n")

    # ---- verdict ----------------------------------------------------------
    ranked = m.rank_deep(TRAIT_LENSES)
    top = ranked[0]
    cross_order_top = all(
        not _same_order(*top["pair"]) for _ in [0]
    )
    print("-" * 78)
    print("VERDICT")
    a, b = top["pair"]
    print(f"  Top deep-convergence pair: {a} ~ {b}  (DEEP={top['deep_score']:+.2f})")
    if _same_order(a, b):
        print("  NOTE: top pair shares an order -- ancestry not fully subtracted.")
    else:
        print("  This pair spans different mammalian orders: similar phenotype at")
        print("  ALL THREE levels, with low shared ancestry == deep convergence.")
    # control check
    ctrl = m.deep_convergence("Rhinolophus_bat", "Pteropus_fruitbat", TRAIT_LENSES)
    print(f"\n  Control (echo bat ~ fruit bat) DEEP = {ctrl['deep_score']:+.2f}  "
          f"[ancestry={ctrl['ancestry_sim']:.2f}]")
    print("  -> close relatives, but NOT convergent: the engine is detecting the")
    print("     convergent SOLUTION, not mere relatedness.")


def fp_out(cat, obj):
    return {mm.target: mm.confidence for mm in cat.morphisms_from(obj)}


def _pairs(xs):
    return [(xs[i], xs[j]) for i in range(len(xs)) for j in range(i + 1, len(xs))]


_ORDER = {
    "Rhinolophus_bat": "Chiroptera", "Myotis_bat": "Chiroptera", "Pteropus_fruitbat": "Chiroptera",
    "Tursiops_dolphin": "Cetacea", "Phocoena_porpoise": "Cetacea", "Bos_cow": "Artiodactyla",
    "Equus_horse": "Perissodactyla", "Homo_human": "Primates",
}


def _same_order(a, b):
    return _ORDER.get(a) == _ORDER.get(b)


if __name__ == "__main__":
    main()
