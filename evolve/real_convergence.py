# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Convergence on REAL data: secondarily-aquatic mammals across 3 orders.

Sample (the only thing we supply -- traits/tree/environment come from the DBs):
  aquatic, independently evolved in:
    Cetacea     : dolphin, orca, beluga, porpoise, fin whale
    Carnivora   : harbor seal, sea lion, walrus, elephant seal, SEA OTTER
    Sirenia     : manatee, dugong
  terrestrial relatives (controls):
    Artiodactyla: cow, pig, hippo(semi-aquatic)
    Carnivora   : wolf, brown bear, lion
    Afrotheria  : African elephant, hyrax
    others      : horse, human

Hypothesis: aquatic mammals converge (body size, climate niche) ACROSS orders,
so deep-convergence should surface cross-order aquatic pairs above their
terrestrial relatives -- detected from PanTHERIA + Open Tree, no curation.
"""

from ingest import build_model

AQUATIC = {
    "Tursiops truncatus", "Orcinus orca", "Delphinapterus leucas",
    "Phocoena phocoena", "Balaenoptera physalus",
    "Phoca vitulina", "Zalophus californianus", "Odobenus rosmarus",
    "Mirounga angustirostris", "Enhydra lutris",
    "Trichechus manatus", "Dugong dugon",
}
TERRESTRIAL = {
    "Bos taurus", "Sus scrofa", "Hippopotamus amphibius",
    "Canis lupus", "Ursus arctos", "Panthera leo",
    "Loxodonta africana", "Procavia capensis",
    "Equus caballus", "Homo sapiens",
}
ORDER = {  # for reading results
    "Tursiops truncatus": "Cetacea", "Orcinus orca": "Cetacea",
    "Delphinapterus leucas": "Cetacea", "Phocoena phocoena": "Cetacea",
    "Balaenoptera physalus": "Cetacea", "Hippopotamus amphibius": "Artiodactyla",
    "Bos taurus": "Artiodactyla", "Sus scrofa": "Artiodactyla",
    "Phoca vitulina": "Carnivora", "Zalophus californianus": "Carnivora",
    "Odobenus rosmarus": "Carnivora", "Mirounga angustirostris": "Carnivora",
    "Enhydra lutris": "Carnivora", "Canis lupus": "Carnivora",
    "Ursus arctos": "Carnivora", "Panthera leo": "Carnivora",
    "Trichechus manatus": "Sirenia", "Dugong dugon": "Sirenia",
    "Loxodonta africana": "Proboscidea", "Procavia capensis": "Hyracoidea",
    "Equus caballus": "Perissodactyla", "Homo sapiens": "Primates",
}
TRAIT_LENSES = ["morphology", "ecology", "environment"]


def life(t):  # aquatic / terrestrial tag
    return "AQUA" if t in AQUATIC else "land"


def main():
    names = sorted(AQUATIC | TERRESTRIAL)
    m, report = build_model(names)
    taxa = report["present"]

    print("\n" + "=" * 80)
    print("DEEP CONVERGENCE on real data (morphology+ecology+environment, minus ancestry)")
    print("=" * 80)
    print(f"{'pair':<50}{'minT':>6}{'anc':>6}{'DEEP':>7}  type")
    print("-" * 80)
    ranked = m.rank_deep(TRAIT_LENSES)
    for d in ranked[:18]:
        a, b = d["pair"]
        cross = ORDER.get(a) != ORDER.get(b)
        kind = f"{life(a)}~{life(b)}" + ("  X-order" if cross else "  same-order")
        print(f"{a+' ~ '+b:<50}{d['min_trait_sim']:>6.2f}"
              f"{d['ancestry_sim']:>6.2f}{d['deep_score']:>7.2f}  {kind}")

    # ---- aggregate: do cross-order aquatic pairs out-converge the rest? ----
    def avg(pred):
        xs = [d["deep_score"] for d in ranked if pred(d)]
        return (sum(xs) / len(xs), len(xs)) if xs else (float("nan"), 0)

    aa_x = avg(lambda d: life(d["pair"][0]) == life(d["pair"][1]) == "AQUA"
               and ORDER.get(d["pair"][0]) != ORDER.get(d["pair"][1]))
    ll_x = avg(lambda d: life(d["pair"][0]) == life(d["pair"][1]) == "land"
               and ORDER.get(d["pair"][0]) != ORDER.get(d["pair"][1]))
    mix = avg(lambda d: life(d["pair"][0]) != life(d["pair"][1]))

    print("\nMean DEEP convergence by pair type (cross-order only):")
    print(f"  aquatic ~ aquatic (diff order) : {aa_x[0]:+.3f}  (n={aa_x[1]})")
    print(f"  land    ~ land    (diff order) : {ll_x[0]:+.3f}  (n={ll_x[1]})")
    print(f"  aquatic ~ land                 : {mix[0]:+.3f}  (n={mix[1]})")
    print("\nVERDICT")
    if aa_x[0] > ll_x[0] and aa_x[0] > mix[0]:
        print("  Aquatic mammals from DIFFERENT orders converge more than expected:")
        print("  real traits track lifestyle/environment across ancestry == convergence.")
    else:
        print("  No clean aquatic-convergence signal in these PanTHERIA lenses.")
        print("  (Likely missing data and/or plesiomorphy -- see notes.)")


if __name__ == "__main__":
    main()
