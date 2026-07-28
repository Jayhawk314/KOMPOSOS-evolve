#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
ESM-C molecular-convergence pilot: Prestin (SLC26A5) and echolocation.

THE QUESTION
------------
`echolocation_flagship.py` shows the bat + toothed-whale echolocation convergence
phenotypically -- but its "molecular lens" is HAND-ASSERTED (a literal
`prestin_echolocator_variant` label). This pilot replaces that asserted label with
a REAL, MEASURED molecular signal and asks the honest question:

    Does a function-aware protein embedding of Prestin INDEPENDENTLY recover the
    echolocation convergence -- i.e. do echolocating bats and toothed whales land
    CLOSER in Prestin embedding space than their phylogenetic distance predicts?

This is the molecular analogue of the morphology lens, and it adds a "molecular
tier": does a phenotypic convergence have an independent molecular signature?
Prestin is the textbook case (Liu et al. 2010, Li et al. 2010: convergent
amino-acid substitutions shared by echolocating bats and toothed whales).

THE LENS
--------
ESM-C (EvolutionaryScale; the lineage after Meta's ESM-2/ESMFold). Input = the
amino-acid sequence; output = a function-aware per-residue embedding. We mean-pool
to one 960-d vector per species -- a molecular relational fingerprint (the Yoneda
lens at the molecular level).

THE GAUNTLET (same philosophy as the rest of the project)
---------------------------------------------------------
Embedding similarity alone is NOT convergence -- close embeddings can also be
homology (it is the same protein in every mammal, >90% identical). So we SUBTRACT
ANCESTRY before claiming anything, exactly as the morphology pipeline does:

  1. embedding distance matrix  d_emb(i,j) = 1 - cosine(v_i, v_j)
  2. patristic distance matrix  d_phylo(i,j) from the REAL dated mammal chronogram
     (Upham 2019) -- independent of the sequences, so it cannot fit the result;
  3. regress d_emb on d_phylo over ALL pairs; the RESIDUAL is the part of molecular
     similarity NOT explained by shared ancestry. A negative residual = closer in
     function space than divergence-time predicts = candidate convergence;
  4. permutation null: shuffle which taxa are "echolocators" and recompute the
     cross-order echolocator residual. p = P(null <= observed).

NEGATIVE CONTROLS (the honesty teeth)
-------------------------------------
  * Pteropus (Old World fruit bat): a BAT that does NOT laryngeally echolocate.
    Must NOT read convergent with toothed whales.
  * Balaenoptera (minke, a baleen whale): a WHALE that does NOT echolocate.
    Must NOT read convergent with echo bats.
If either control reads convergent, the lens is detecting "battiness/whaleyness",
not the echolocation solution -- and the pilot says so.

HONESTY CAVEATS (stated, per the project contract)
--------------------------------------------------
  * Real sequences only (NCBI RefSeq/Swiss-Prot). A taxon with no full-length
    Prestin is dropped and reported, never imputed.
  * ESM distances are used ONLY as a function lens -- never to build a tree
    (alignment-free embedding phylogeny is confounded). The tree is external.
  * One protein, ~14 taxa: this is a PILOT that proves the lens, not a survey.

Run:  python evolve/esmc_prestin_pilot.py
Out:  evolve/results/prestin_esmc.json   (+ cached FASTA / embeddings under cache/)
"""

from __future__ import annotations

import json
import sys
import time
from io import StringIO
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from chronogram import load_mammal_chronogram, all_tip_labels                  # noqa: E402
from polarity import _index_parents, _mrca                                     # noqa: E402

CACHE = HERE / "cache" / "prestin"
RESULTS = HERE / "results"
EMAIL = "jhawk314@gmail.com"   # NCBI Entrez requires a contact address

# group: "echo"   = laryngeal echolocator (the convergent trait)
#        "control"= bat/whale that does NOT echolocate (negative control)
#        "out"    = outgroup mammal
# chronogram = Upham 2019 tip binomial (underscore); genus fallback if absent.
TAXA = [
    # organism (NCBI)             label              group      chronogram
    ("Rhinolophus ferrumequinum", "Rhinolophus_bat",  "echo",    "Rhinolophus_ferrumequinum"),
    ("Myotis lucifugus",          "Myotis_bat",       "echo",    "Myotis_lucifugus"),
    ("Eptesicus fuscus",          "Eptesicus_bat",    "echo",    "Eptesicus_fuscus"),
    ("Tursiops truncatus",        "Tursiops_dolphin", "echo",    "Tursiops_truncatus"),
    ("Phocoena phocoena",         "Phocoena_porpoise","echo",    "Phocoena_phocoena"),
    ("Physeter catodon",          "Physeter_whale",   "echo",    "Physeter_macrocephalus"),
    ("Orcinus orca",              "Orcinus_orca",     "echo",    "Orcinus_orca"),
    ("Pteropus vampyrus",         "Pteropus_fruitbat","control", "Pteropus_vampyrus"),
    ("Balaenoptera acutorostrata","Balaenoptera_baleen","control","Balaenoptera_acutorostrata"),
    ("Bos taurus",                "Bos_cow",          "out",     "Bos_taurus"),
    ("Equus caballus",            "Equus_horse",      "out",     "Equus_caballus"),
    ("Loxodonta africana",        "Loxodonta_elephant","out",    "Loxodonta_africana"),
    ("Mus musculus",              "Mus_mouse",        "out",     "Mus_musculus"),
    ("Homo sapiens",              "Homo_human",       "out",     "Homo_sapiens"),
]

# order, for the cross-order permutation statistic (echolocation spans Chiroptera
# x Cetacea; same-order echo pairs are similar partly by descent and are excluded).
ORDER = {
    "Rhinolophus_bat": "Chiroptera", "Myotis_bat": "Chiroptera",
    "Eptesicus_bat": "Chiroptera", "Pteropus_fruitbat": "Chiroptera",
    "Tursiops_dolphin": "Cetacea", "Phocoena_porpoise": "Cetacea",
    "Physeter_whale": "Cetacea", "Orcinus_orca": "Cetacea",
    "Balaenoptera_baleen": "Cetacea",
    "Bos_cow": "Artiodactyla", "Equus_horse": "Perissodactyla",
    "Loxodonta_elephant": "Proboscidea", "Mus_mouse": "Rodentia",
    "Homo_human": "Primates",
}


# ---------------------------------------------------------------------------
# 1. Sequences (NCBI Entrez -> cached FASTA)
# ---------------------------------------------------------------------------
def _retry(fn, *a, tries=4, **k):
    """Entrez over a flaky link drops connections; retry with backoff."""
    for i in range(tries):
        try:
            return fn(*a, **k)
        except Exception as e:
            if i == tries - 1:
                raise
            print(f"    (retry {i+1}/{tries} after {type(e).__name__})")
            time.sleep(1.5 * (i + 1))


def fetch_sequences() -> dict[str, dict]:
    from Bio import Entrez, SeqIO
    Entrez.email = EMAIL
    CACHE.mkdir(parents=True, exist_ok=True)
    out: dict[str, dict] = {}
    for organism, label, group, chron in TAXA:
        f = CACHE / f"{label}.fasta"
        if f.exists():
            rec = next(SeqIO.parse(str(f), "fasta"))
            out[label] = {"seq": str(rec.seq), "acc": rec.id, "organism": organism,
                          "group": group, "chron": chron}
            continue
        term = (f'(SLC26A5[Gene Name] OR prestin[Protein Name]) '
                f'AND "{organism}"[Organism] AND 650:820[Sequence Length]')
        h = _retry(Entrez.esearch, db="protein", term=term, retmax=25)
        ids = Entrez.read(h)["IdList"]; h.close()
        time.sleep(0.34)
        if not ids:
            print(f"  [{label:<20}] NO full-length Prestin found -- dropped (honest)")
            continue
        h = _retry(Entrez.efetch, db="protein", id=ids, rettype="fasta", retmode="text")
        recs = [r for r in SeqIO.parse(StringIO(h.read()), "fasta") if 650 <= len(r.seq) <= 820]
        h.close(); time.sleep(0.34)
        if not recs:
            print(f"  [{label:<20}] only fragments -- dropped (honest)")
            continue
        best = max(recs, key=lambda r: len(r.seq))
        f.write_text(f">{best.id} {organism}\n{best.seq}\n")
        out[label] = {"seq": str(best.seq), "acc": best.id, "organism": organism,
                      "group": group, "chron": chron}
        print(f"  [{label:<20}] {best.id}  ({len(best.seq)} aa)")
    return out


# ---------------------------------------------------------------------------
# 2. ESM-C embeddings (mean-pooled, cached)
# ---------------------------------------------------------------------------
def embed(seqs: dict[str, dict]):
    import numpy as np
    npz = CACHE / "embeddings.npz"
    labels = list(seqs)
    key = {l: seqs[l]["acc"] for l in labels}
    if npz.exists():
        z = np.load(npz, allow_pickle=True)
        if list(z["labels"]) == labels and dict(z["accs"].item()) == key:
            return {l: z["V"][i] for i, l in enumerate(labels)}
    import torch
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig
    print("  loading ESM-C (esmc_300m) ...")
    client = ESMC.from_pretrained("esmc_300m").to("cpu").eval()
    V = []
    for l in labels:
        p = client.encode(ESMProtein(sequence=seqs[l]["seq"]))
        with torch.no_grad():
            out = client.logits(p, LogitsConfig(sequence=True, return_embeddings=True))
        v = out.embeddings[0, 1:-1, :].mean(0).float().cpu().numpy()  # drop BOS/EOS
        V.append(v / (np.linalg.norm(v) + 1e-9))
        print(f"    embedded {l:<20} dim={v.shape[0]}")
    V = np.vstack(V)
    np.savez(npz, labels=np.array(labels), V=V, accs=np.array(key, dtype=object))
    return {l: V[i] for i, l in enumerate(labels)}


# ---------------------------------------------------------------------------
# 3. Patristic distances from the REAL dated mammal chronogram
# ---------------------------------------------------------------------------
def patristic(labels_to_chron: dict[str, str]) -> dict[tuple[str, str], float]:
    tree = load_mammal_chronogram()
    tips = all_tip_labels(tree)
    parent, leaves = _index_parents(tree)

    def resolve(binom: str) -> str | None:
        if binom in tips:
            return binom
        genus = binom.split("_")[0] + "_"
        cands = sorted(t for t in tips if t.startswith(genus))
        return cands[0] if cands else None

    placed = {l: resolve(c) for l, c in labels_to_chron.items()}
    for l, t in placed.items():
        if t is None:
            print(f"  [phylo] {l}: not on chronogram -- excluded from phylo control")
        elif t != labels_to_chron[l]:
            print(f"  [phylo] {l}: {labels_to_chron[l]} absent -> using congener {t}")

    def dist_to_root(node) -> float:
        d, n = 0.0, node
        while id(n) in parent:
            d += n.brlen or 0.0
            n = parent[id(n)]
        return d

    droot = {l: dist_to_root(leaves[t]) for l, t in placed.items() if t}
    D: dict[tuple[str, str], float] = {}
    have = [l for l, t in placed.items() if t]
    for i in range(len(have)):
        for j in range(i + 1, len(have)):
            a, b = have[i], have[j]
            m = _mrca(leaves[placed[a]], leaves[placed[b]], parent)
            D[(a, b)] = D[(b, a)] = droot[a] + droot[b] - 2 * dist_to_root(m)
    return D


# ---------------------------------------------------------------------------
# 4. Subtract ancestry (regress out phylogeny) + permutation test
# ---------------------------------------------------------------------------
def analyze(vecs, phylo, groups):
    import numpy as np
    labels = [l for l in vecs if any((l, o) in phylo for o in vecs)]
    pairs = [(a, b) for i, a in enumerate(labels) for b in labels[labels.index(a) + 1:]
             if (a, b) in phylo]

    d_emb = np.array([1.0 - float(np.dot(vecs[a], vecs[b])) for a, b in pairs])
    d_phy = np.array([phylo[(a, b)] for a, b in pairs])
    d_phy = d_phy / d_phy.max()                      # normalize timescale to [0,1]

    # OLS d_emb = a + b*d_phylo ; residual = ancestry-subtracted molecular distance
    A = np.vstack([np.ones_like(d_phy), d_phy]).T
    coef, *_ = np.linalg.lstsq(A, d_emb, rcond=None)
    resid = d_emb - A @ coef
    res = {p: float(r) for p, r in zip(pairs, resid)}

    echo = {l for l in labels if groups[l] == "echo"}

    def cross_order_echo_mean(echo_set):
        rs = [res[p] for p in pairs
              if p[0] in echo_set and p[1] in echo_set and ORDER[p[0]] != ORDER[p[1]]]
        return float(np.mean(rs)) if rs else float("nan")

    s_obs = cross_order_echo_mean(echo)

    # permutation null: reassign the 'echo' label to a random same-size taxon set,
    # recompute the cross-order echo-pair mean residual. Phylogeny already removed.
    rng = np.random.default_rng(0)
    pool = list(labels); k = len(echo)
    null = []
    for _ in range(5000):
        e2 = set(rng.choice(pool, size=k, replace=False))
        if sum(1 for o in {ORDER[x] for x in e2}) < 2:
            continue
        v = cross_order_echo_mean(e2)
        if v == v:
            null.append(v)
    null = np.array(null)
    p = float((np.sum(null <= s_obs) + 1) / (len(null) + 1))
    return res, pairs, s_obs, p, sorted(res.items(), key=lambda kv: kv[1])


# ---------------------------------------------------------------------------
def main():
    print("Fetching real Prestin (SLC26A5) sequences from NCBI ...")
    seqs = fetch_sequences()
    print(f"  got {len(seqs)} sequences")

    print("Embedding with ESM-C ...")
    vecs = embed(seqs)

    print("Patristic distances on the dated mammal chronogram (Upham 2019) ...")
    phylo = patristic({l: seqs[l]["chron"] for l in seqs})

    print("Subtracting ancestry + permutation test ...")
    groups = {l: seqs[l]["group"] for l in seqs}
    res, pairs, s_obs, p, ranked = analyze(vecs, phylo, groups)

    def named(a, b):
        return res.get((a, b), res.get((b, a)))

    echo_bats = [l for l in seqs if groups[l] == "echo" and ORDER[l] == "Chiroptera"]
    odonto = [l for l in seqs if groups[l] == "echo" and ORDER[l] == "Cetacea"]

    print("\n" + "=" * 74)
    print("ESM-C PRESTIN PILOT  --  molecular signature of echolocation convergence")
    print("=" * 74)

    print("\nMost convergent pairs (closest in Prestin function-space AFTER removing")
    print("phylogeny; negative = closer than divergence-time predicts):\n")
    for (a, b), r in ranked[:8]:
        tag = "  bat x toothed-whale" if {ORDER[a], ORDER[b]} == {"Chiroptera", "Cetacea"} \
              and groups[a] == groups[b] == "echo" else ""
        print(f"   {a+' ~ '+b:<42} resid={r:+.4f}{tag}")

    print(f"\nHEADLINE statistic  (mean residual, echo bat x toothed whale pairs):")
    print(f"   observed = {s_obs:+.4f}   permutation p = {p:.4f}"
          f"   ({'SIGNIFICANT' if p < 0.05 else 'not significant'})")

    print("\nNEGATIVE CONTROLS (must NOT read convergent with the other clade):")
    fb = "Pteropus_fruitbat"
    bw = "Balaenoptera_baleen"
    if fb in seqs:
        vals = [named(fb, o) for o in odonto if named(fb, o) is not None]
        m = sum(vals) / len(vals) if vals else float("nan")
        print(f"   fruit bat  (Pteropus)  vs toothed whales : mean resid {m:+.4f}"
              f"   {'OK (not convergent)' if m >= 0 else 'WARN: reads convergent'}")
    if bw in seqs:
        vals = [named(bw, o) for o in echo_bats if named(bw, o) is not None]
        m = sum(vals) / len(vals) if vals else float("nan")
        print(f"   baleen whale (Balaenoptera) vs echo bats : mean resid {m:+.4f}"
              f"   {'OK (not convergent)' if m >= 0 else 'WARN: reads convergent'}")

    # echo bat x toothed whale, the convergent contrast
    conv_vals = [named(a, b) for a in echo_bats for b in odonto if named(a, b) is not None]
    conv_mean = sum(conv_vals) / len(conv_vals) if conv_vals else float("nan")

    controls_ok = True
    for ctrl, others in [(fb, odonto), (bw, echo_bats)]:
        if ctrl in seqs:
            vals = [named(ctrl, o) for o in others if named(ctrl, o) is not None]
            if vals and sum(vals) / len(vals) < conv_mean:
                controls_ok = False

    print("\n" + "-" * 74)
    print("VERDICT")
    verdict = (p < 0.05 and conv_mean < 0 and controls_ok)
    print(f"  Echo bat x toothed-whale Prestin is {'CLOSER' if conv_mean < 0 else 'NOT closer'} "
          f"in ESM-C space than phylogeny predicts (mean resid {conv_mean:+.4f}),")
    print(f"  permutation p = {p:.4f}; negative controls "
          f"{'behave (fruit bat + baleen whale NOT convergent)' if controls_ok else 'FAIL'}.")
    if verdict:
        print("  => The MOLECULAR lens INDEPENDENTLY recovers the echolocation convergence")
        print("     that echolocation_flagship.py shows phenotypically. Convergence-of-")
        print("     convergence: a phenotypic convergence with an independent molecular tier.")
    else:
        print("  => Molecular tier NOT cleanly recovered at this sample/threshold "
              "(reported honestly, not inflated).")

    RESULTS.mkdir(parents=True, exist_ok=True)
    outp = RESULTS / "prestin_esmc.json"
    outp.write_text(json.dumps({
        "taxa": {l: {"acc": seqs[l]["acc"], "organism": seqs[l]["organism"],
                     "group": seqs[l]["group"]} for l in seqs},
        "n_seqs": len(seqs),
        "headline_mean_resid_bat_x_toothedwhale": conv_mean,
        "permutation_p": p,
        "controls_ok": controls_ok,
        "verdict_molecular_convergence": verdict,
        "ranked_residuals": [{"pair": list(k), "resid": v} for k, v in ranked],
    }, indent=2))
    print(f"\n  wrote {outp}")


if __name__ == "__main__":
    main()
