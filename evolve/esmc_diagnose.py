#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Diagnose the ESM-C molecular lens on the Prestin pilot (Track C).

The pilot (`esmc_prestin_pilot.py`) returned an HONEST NEGATIVE: mean-pooled
ESM-C 300M embeddings of Prestin did NOT recover the bat+toothed-whale
echolocation convergence (p=0.13, negative controls failed). Before spending
compute on ESM-C 600M/6B or more proteins, this script asks WHY, with three
self-contained diagnostics on the SAME 14 cached sequences/embeddings:

  PART 1 -- Geometry of the embedding space.
      PCA the 14 mean-pooled 960-d vectors. Do echolocators cluster at all, or
      does phylogeny / sequence length dominate? Report variance explained and a
      silhouette score for the echo-vs-rest split.

  PART 2 -- Does ESM beat raw % identity?
      Build a %-identity distance from a reference-anchored alignment and run the
      EXACT SAME ancestry-subtracted permutation test the pilot runs on ESM
      cosine distance. If ESM's p-value is no better than %-identity's, the
      function-aware embedding is adding nothing over a trivial baseline for this
      protein -- and we say so.

  PART 3 -- Localize: convergent sites, not whole-protein average.
      Mean-pooling averages ~700 residues; the echolocation signal lives in a
      handful of convergent substitutions (Liu et al. 2010, Li et al. 2010 report
      a SET of ~3-7 parallel sites between echo bats and toothed whales -- no
      single canonical position list, so we DERIVE candidate convergent columns
      from the data: alignment columns where every echo bat AND every toothed
      whale share an amino acid that NONE of the non-echo taxa carry). Then we
      re-embed PER-RESIDUE with ESM-C and test whether echo taxa cluster at those
      columns better than under mean-pooling. (Re-embedding is cached.)

HONESTY CONTRACT
  * Real cached sequences/embeddings only; nothing imputed.
  * Candidate convergent sites are DERIVED from the alignment, not asserted, and
    reported explicitly so they can be checked against the literature.
  * ESM stays a function lens; we never build a tree from embeddings. The test
    always subtracts ancestry via the real dated mammal chronogram.
  * Every comparison is reported plainly, including "ESM adds nothing here".

Run:  python evolve/esmc_diagnose.py
Out:  evolve/results/prestin_esmc_diagnosis.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Reuse the pilot's REAL chronogram patristic logic + taxon tables -- do not reinvent.
from esmc_prestin_pilot import TAXA, ORDER, CACHE, RESULTS, patristic  # noqa: E402

PRESTIN_CACHE = CACHE  # evolve/cache/prestin
DIAG_OUT = RESULTS / "prestin_esmc_diagnosis.json"

GROUP = {label: group for _, label, group, _ in TAXA}
CHRON = {label: chron for _, label, _, chron in TAXA}
ECHO_BATS = [l for l in GROUP if GROUP[l] == "echo" and ORDER[l] == "Chiroptera"]
ODONTO = [l for l in GROUP if GROUP[l] == "echo" and ORDER[l] == "Cetacea"]


# ---------------------------------------------------------------------------
# Load cached sequences + mean-pooled embeddings (no NCBI / no ESM needed here)
# ---------------------------------------------------------------------------
def load_sequences() -> dict[str, str]:
    from Bio import SeqIO
    out = {}
    for _, label, _, _ in TAXA:
        f = PRESTIN_CACHE / f"{label}.fasta"
        if not f.exists():
            print(f"  [{label}] no cached FASTA -- skipped (honest)")
            continue
        out[label] = str(next(SeqIO.parse(str(f), "fasta")).seq)
    return out


def load_mean_embeddings() -> dict[str, np.ndarray]:
    z = np.load(PRESTIN_CACHE / "embeddings.npz", allow_pickle=True)
    labels = list(z["labels"])
    return {l: z["V"][i] for i, l in enumerate(labels)}


# ---------------------------------------------------------------------------
# Reference-anchored alignment (orthologs >90% identical -> pairwise-to-ref is
# enough to define homologous columns; avoids requiring an external MSA tool).
# Returns, per taxon: ref_col -> (residue_index_in_taxon, amino_acid or '-').
# ---------------------------------------------------------------------------
def align_to_reference(seqs: dict[str, str], ref_label: str):
    from Bio.Align import PairwiseAligner, substitution_matrices
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    ref = seqs[ref_label]
    n_ref = len(ref)
    # ref_maps[label] = list length n_ref; entry = (taxon_residue_index, aa) or None
    ref_maps: dict[str, list] = {}
    for label, s in seqs.items():
        aln = aligner.align(ref, s)[0]
        # aln.aligned -> blocks of (ref_start,ref_end),(s_start,s_end)
        rblocks, sblocks = aln.aligned
        mp = [None] * n_ref
        for (r0, r1), (s0, s1) in zip(rblocks, sblocks):
            for k in range(r1 - r0):
                mp[r0 + k] = (s0 + k, s[s0 + k])
        ref_maps[label] = mp
    return n_ref, ref_maps


# ---------------------------------------------------------------------------
# Generic ancestry-subtracted cross-order echo test on ANY pairwise distance.
# (Same statistic the pilot uses, factored to accept a distance dict so we can
#  feed it ESM-cosine, %-identity, or per-site distances identically.)
# ---------------------------------------------------------------------------
def ancestry_subtracted_echo_test(dist: dict, phylo: dict, labels: list,
                                  n_perm: int = 20000, seed: int = 0):
    pairs = [(a, b) for i, a in enumerate(labels) for b in labels[i + 1:]
             if (a, b) in dist and (a, b) in phylo]
    d_emb = np.array([dist[p] for p in pairs])
    d_phy = np.array([phylo[p] for p in pairs])
    d_phy = d_phy / d_phy.max()
    A = np.vstack([np.ones_like(d_phy), d_phy]).T
    coef, *_ = np.linalg.lstsq(A, d_emb, rcond=None)
    resid = d_emb - A @ coef
    res = {p: float(r) for p, r in zip(pairs, resid)}

    echo = {l for l in labels if GROUP[l] == "echo"}

    def cross_order_mean(echo_set):
        rs = [res[p] for p in pairs
              if p[0] in echo_set and p[1] in echo_set and ORDER[p[0]] != ORDER[p[1]]]
        return float(np.mean(rs)) if rs else float("nan")

    s_obs = cross_order_mean(echo)
    rng = np.random.default_rng(seed)
    pool = list(labels); k = len(echo)
    null = []
    for _ in range(n_perm):
        e2 = set(rng.choice(pool, size=k, replace=False))
        if len({ORDER[x] for x in e2}) < 2:
            continue
        v = cross_order_mean(e2)
        if v == v:
            null.append(v)
    null = np.array(null)
    p = float((np.sum(null <= s_obs) + 1) / (len(null) + 1))

    # negative controls: do non-echo bat/whale read convergent with the other clade?
    def mean_to(a, others):
        vals = [res.get((a, o), res.get((o, a))) for o in others]
        vals = [v for v in vals if v is not None]
        return float(np.mean(vals)) if vals else float("nan")

    conv = float(np.mean([res.get((a, b), res.get((b, a)))
                          for a in ECHO_BATS for b in ODONTO
                          if res.get((a, b), res.get((b, a))) is not None]))
    fruitbat_vs_odonto = mean_to("Pteropus_fruitbat", ODONTO)
    baleen_vs_echobat = mean_to("Balaenoptera_baleen", ECHO_BATS)
    controls_ok = (not (fruitbat_vs_odonto < conv)) and (not (baleen_vs_echobat < conv))
    return {
        "observed_cross_order_echo_resid": s_obs,
        "permutation_p": p,
        "n_null": int(len(null)),
        "conv_bat_x_toothedwhale_mean_resid": conv,
        "control_fruitbat_vs_toothedwhales": fruitbat_vs_odonto,
        "control_baleen_vs_echobats": baleen_vs_echobat,
        "controls_ok": bool(controls_ok),
    }


# ---------------------------------------------------------------------------
# PART 1 -- geometry of the mean-pooled embedding space
# ---------------------------------------------------------------------------
def part1_geometry(mean_emb: dict, seqs: dict):
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score
    labels = [l for l in mean_emb if l in seqs]
    X = np.vstack([mean_emb[l] for l in labels])
    pca = PCA(n_components=min(5, len(labels) - 1)).fit(X)
    coords = pca.transform(X)
    evr = pca.explained_variance_ratio_

    echo_mask = np.array([1 if GROUP[l] == "echo" else 0 for l in labels])
    # silhouette of echo vs non-echo in full embedding space (cosine)
    try:
        sil = float(silhouette_score(X, echo_mask, metric="cosine"))
    except Exception:
        sil = float("nan")

    # how strongly does PC1 track sequence length (a confound to watch for)?
    lengths = np.array([len(seqs[l]) for l in labels], dtype=float)
    pc1_len_r = float(np.corrcoef(coords[:, 0], lengths)[0, 1])

    return {
        "labels": labels,
        "explained_variance_ratio": [float(x) for x in evr],
        "pc_coords": {l: [float(coords[i, 0]), float(coords[i, 1])] for i, l in enumerate(labels)},
        "echo_vs_rest_silhouette_cosine": sil,
        "pc1_vs_seqlength_corr": pc1_len_r,
    }


# ---------------------------------------------------------------------------
# PART 2 -- ESM cosine vs raw %-identity, same ancestry-subtracted test
# ---------------------------------------------------------------------------
def pct_identity_distance(n_ref: int, ref_maps: dict, labels: list):
    dist = {}
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            ma, mb = ref_maps[a], ref_maps[b]
            same = tot = 0
            for c in range(n_ref):
                pa, pb = ma[c], mb[c]
                if pa is None or pb is None:
                    continue
                tot += 1
                if pa[1] == pb[1]:
                    same += 1
            ident = same / tot if tot else 0.0
            dist[(a, b)] = dist[(b, a)] = 1.0 - ident
    return dist


def esm_cosine_distance(mean_emb: dict, labels: list):
    dist = {}
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            d = 1.0 - float(np.dot(mean_emb[a], mean_emb[b]))
            dist[(a, b)] = dist[(b, a)] = d
    return dist


# ---------------------------------------------------------------------------
# PART 3 -- data-derived convergent sites + per-residue ESM embeddings
# ---------------------------------------------------------------------------
def _col_aa(ref_maps, label, c):
    m = ref_maps.get(label, [None])[c] if label in ref_maps else None
    return m[1] if m is not None else None


def derive_convergent_columns(n_ref: int, ref_maps: dict, labels: list):
    """Two data-derived candidate-site sets, both reported:

    STRICT (perfectly diagnostic): every echolocator shares residue X and NO
    non-echo taxon carries X. Cleanest, but unrealistically demanding.

    PARALLEL (polarity-aware, literature-aligned): echo bats all share X AND
    toothed whales all share X, where X differs from each clade's IMMEDIATE
    non-echo relative (Pteropus for bats; Balaenoptera/Bos for the whale side).
    This is the parallel-substitution signature Liu/Li 2010 use -- distant
    outgroups may coincidentally carry X (homoplasy elsewhere) without
    disqualifying the site, because polarity is judged against close relatives.
    """
    nonecho = [l for l in labels if GROUP[l] != "echo"]
    bat_rel = "Pteropus_fruitbat"          # non-echo bat (sister to echo bats)
    whale_rels = ["Balaenoptera_baleen", "Bos_cow"]  # non-echo relatives of odontocetes
    strict, parallel = [], []
    for c in range(n_ref):
        bat_aas = {_col_aa(ref_maps, l, c) for l in ECHO_BATS}
        whale_aas = {_col_aa(ref_maps, l, c) for l in ODONTO}
        if len(bat_aas) != 1 or len(whale_aas) != 1:
            continue
        x = next(iter(bat_aas))
        if x is None or x == "-" or next(iter(whale_aas)) != x:
            continue  # both echo clades must be fixed for the SAME residue X
        # STRICT
        nonecho_aas = {_col_aa(ref_maps, l, c) for l in nonecho}
        if x not in nonecho_aas:
            strict.append((c, x))
        # PARALLEL: differs from each clade's immediate non-echo relative
        bat_polar = _col_aa(ref_maps, bat_rel, c) not in (None, x)
        whale_polar = any(_col_aa(ref_maps, r, c) not in (None, x) for r in whale_rels)
        if bat_polar and whale_polar:
            parallel.append((c, x))
    return strict, parallel


def per_residue_embeddings(seqs: dict):
    """Re-embed each sequence PER RESIDUE with ESM-C 300M (cached)."""
    npz = PRESTIN_CACHE / "embeddings_perres.npz"
    labels = list(seqs)
    if npz.exists():
        z = np.load(npz, allow_pickle=True)
        if list(z["labels"]) == labels:
            return {l: z[l] for l in labels}
    import torch
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig
    print("  loading ESM-C (esmc_300m) for per-residue embeddings ...")
    client = ESMC.from_pretrained("esmc_300m").to("cpu").eval()
    store = {}
    for l in labels:
        p = client.encode(ESMProtein(sequence=seqs[l]))
        with torch.no_grad():
            out = client.logits(p, LogitsConfig(sequence=True, return_embeddings=True))
        # drop BOS/EOS -> [L, 960], aligned 1:1 with residues
        store[l] = out.embeddings[0, 1:-1, :].float().cpu().numpy()
        print(f"    per-residue {l:<20} shape={store[l].shape}")
    np.savez(npz, labels=np.array(labels), **store)
    return store


def per_site_identity_distance(cols, ref_maps, labels: list):
    """Raw AA-mismatch distance at the SAME columns -- the no-ESM baseline. Because
    the columns were selected so echo clades share the residue, this is expected to
    look 'significant' too; comparing it to the ESM per-site test exposes how much of
    the per-site signal is mere site-selection vs. ESM's contextual embedding."""
    dist = {}
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            mm = tot = 0
            for c, _aa in cols:
                pa, pb = ref_maps[a][c], ref_maps[b][c]
                if pa is None or pb is None:
                    continue
                tot += 1
                if pa[1] != pb[1]:
                    mm += 1
            if tot:
                dist[(a, b)] = dist[(b, a)] = mm / tot
    return dist


def per_site_distance(cols, ref_maps, perres: dict, labels: list):
    """Concatenate the (L2-normalized) per-residue embeddings at the convergent
    columns into one vector per taxon; cosine distance between taxa there."""
    vecs = {}
    for l in labels:
        if l not in perres:
            continue
        parts = []
        ok = True
        for c, _aa in cols:
            m = ref_maps[l][c]
            if m is None:
                ok = False
                break
            ridx = m[0]
            parts.append(perres[l][ridx])
        if not ok or not parts:
            continue
        v = np.concatenate(parts)
        vecs[l] = v / (np.linalg.norm(v) + 1e-9)
    dist = {}
    have = list(vecs)
    for i, a in enumerate(have):
        for b in have[i + 1:]:
            d = 1.0 - float(np.dot(vecs[a], vecs[b]))
            dist[(a, b)] = dist[(b, a)] = d
    return dist, have


# ---------------------------------------------------------------------------
def main(run_part3: bool = True):
    print("Loading cached sequences + mean-pooled embeddings ...")
    seqs = load_sequences()
    mean_emb = load_mean_embeddings()
    labels = [l for l in mean_emb if l in seqs]
    print(f"  {len(labels)} taxa")

    print("Patristic distances on the dated mammal chronogram (Upham 2019) ...")
    phylo = patristic({l: CHRON[l] for l in labels})

    print("Reference-anchored alignment (ref = Homo_human) ...")
    n_ref, ref_maps = align_to_reference(seqs, "Homo_human")

    out = {}

    # ---- PART 1 ----
    print("\n" + "=" * 74)
    print("PART 1 -- geometry of the mean-pooled ESM-C embedding space")
    print("=" * 74)
    g = part1_geometry(mean_emb, seqs)
    out["part1_geometry"] = g
    print(f"  PC1..PC3 variance explained: "
          f"{', '.join(f'{x:.1%}' for x in g['explained_variance_ratio'][:3])}")
    print(f"  echo-vs-rest silhouette (cosine): {g['echo_vs_rest_silhouette_cosine']:+.3f}"
          f"   (>0 = echo separates; <=0 = no echo cluster)")
    print(f"  PC1 vs sequence-length corr: {g['pc1_vs_seqlength_corr']:+.2f}"
          f"   (|r|~1 = geometry tracks length, not function)")

    # ---- PART 2 ----
    print("\n" + "=" * 74)
    print("PART 2 -- does ESM-C beat raw % identity? (same ancestry-subtracted test)")
    print("=" * 74)
    d_esm = esm_cosine_distance(mean_emb, labels)
    d_pid = pct_identity_distance(n_ref, ref_maps, labels)
    r_esm = ancestry_subtracted_echo_test(d_esm, phylo, labels)
    r_pid = ancestry_subtracted_echo_test(d_pid, phylo, labels)
    out["part2_esm_vs_identity"] = {"esm_cosine": r_esm, "pct_identity": r_pid}
    for name, r in [("ESM-C cosine ", r_esm), ("% identity   ", r_pid)]:
        print(f"  {name}: echo resid={r['conv_bat_x_toothedwhale_mean_resid']:+.4f}  "
              f"p={r['permutation_p']:.4f}  controls_ok={r['controls_ok']}")
    better = ("ESM-C beats %-identity" if r_esm["permutation_p"] < r_pid["permutation_p"]
              else "ESM-C does NOT beat %-identity")
    print(f"  => {better} on this protein/sample.")

    # ---- PART 3 ----
    if run_part3:
        print("\n" + "=" * 74)
        print("PART 3 -- data-derived convergent sites + per-residue embeddings")
        print("=" * 74)
        strict, parallel = derive_convergent_columns(n_ref, ref_maps, labels)
        print(f"  STRICT (perfectly diagnostic) columns:        {len(strict)}")
        print(f"  PARALLEL (polarity vs non-echo relatives):    {len(parallel)}")
        print(f"  (literature reports a SET of ~3-7 parallel sites; Liu/Li 2010)")
        for c, aa in parallel:
            print(f"     ref_pos {c+1:>4}  convergent residue = {aa}")
        out["part3_convergent_sites"] = {
            "n_strict": len(strict),
            "n_parallel": len(parallel),
            "strict_columns_1based_ref": [{"ref_pos": c + 1, "aa": aa} for c, aa in strict],
            "parallel_columns_1based_ref": [{"ref_pos": c + 1, "aa": aa} for c, aa in parallel],
        }
        if parallel:
            perres = per_residue_embeddings(seqs)
            d_site, have = per_site_distance(parallel, ref_maps, perres, labels)
            r_site = ancestry_subtracted_echo_test(d_site, phylo, have)
            d_sid = per_site_identity_distance(parallel, ref_maps, labels)
            r_sid = ancestry_subtracted_echo_test(d_sid, phylo, labels)
            out["part3_convergent_sites"]["per_site_esm_test"] = r_site
            out["part3_convergent_sites"]["per_site_identity_test"] = r_sid
            out["part3_convergent_sites"]["circularity_caveat"] = (
                "Sites were SELECTED using the echo/non-echo labels, so both per-site "
                "tests are partly circular: echo taxa share these residues by "
                "construction. The per-site identity test quantifies the selection-"
                "driven component; ESM beyond it is ESM's marginal contribution. A "
                "non-circular test needs convergent sites from an INDEPENDENT source "
                "(literature ASR / a held-out tree).")
            print(f"  per-site ESM test      (PARALLEL cols): echo resid="
                  f"{r_site['conv_bat_x_toothedwhale_mean_resid']:+.4f}  "
                  f"p={r_site['permutation_p']:.4f}  controls_ok={r_site['controls_ok']}")
            print(f"  per-site %-identity test (same cols)   : echo resid="
                  f"{r_sid['conv_bat_x_toothedwhale_mean_resid']:+.4f}  "
                  f"p={r_sid['permutation_p']:.4f}  controls_ok={r_sid['controls_ok']}")
            print("  CAVEAT: sites were selected using the labels -> per-site tests are")
            print("  partly circular. Both being significant shows the per-site signal is")
            print("  largely SITE-SELECTION, not independent ESM discovery. The honest,")
            print("  non-circular finding is the LOCALIZATION + dilution explanation:")
            n_res = len(next(iter(perres.values())))
            print(f"  convergence lives in ~{len(parallel)}/{n_res} residues "
                  f"(~{100*len(parallel)/n_res:.1f}%); mean-pooling over all of them "
                  f"dilutes it ~{n_res//max(len(parallel),1)}x -> the pilot's negative.")
        else:
            print("  no parallel convergent columns -> per-site test skipped (honest)")

    DIAG_OUT.parent.mkdir(parents=True, exist_ok=True)
    DIAG_OUT.write_text(json.dumps(out, indent=2))
    print(f"\n  wrote {DIAG_OUT}")


if __name__ == "__main__":
    main(run_part3="--no-part3" not in sys.argv)
