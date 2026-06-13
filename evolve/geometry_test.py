"""
Do the repo's Ricci / spectral modules (a) run, and (b) add anything for
convergence?

We build a TAXON-TO-TAXON similarity graph from the Anolis morphology lens
(edge = trait-similarity above a threshold), then ask:

  SPECTRAL : does the graph Laplacian's spectral embedding cluster lizards by
             ECOMORPH (the convergent classes) rather than by island?
  RICCI    : Ollivier-Ricci curvature flags "bridge" edges between communities.
             A convergent pair = distant ancestry linked by shared traits = a
             bridge. Are the most NEGATIVELY curved edges the cross-island
             same-ecomorph (convergent) ones?
"""

import sys
from pathlib import Path
from itertools import combinations

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from anolis_validation import build_lenses, yoneda, SPECIES  # noqa: E402
from komposos_core.core.category import Category  # noqa: E402
from komposos_wesys.geometry.grid_spectral import SpectralGraphAnalyzer  # noqa: E402
from komposos_wesys.geometry.grid_ricci import OllivierRicciCurvature  # noqa: E402

THRESH = 0.30


def taxon_graph() -> Category:
    """Taxon-taxon similarity graph from the morphology lens."""
    morph, _ = build_lenses()
    g = Category(db_path=":memory:")
    species = list(SPECIES)
    for s in species:
        g.add(s)
    for a, b in combinations(species, 2):
        sim = yoneda(morph, a, b)
        if sim >= THRESH:
            g.connect(a, b, name=f"sim:{a}:{b}", confidence=float(sim))
    return g


class _RicciAdapter:
    """Expose a komposos Category through the interface OllivierRicciCurvature
    expects: list_morphisms(limit) -> objs with .source_name/.target_name/.confidence."""
    def __init__(self, cat): self.cat = cat
    def list_morphisms(self, limit=100000):
        class M:
            def __init__(s, m): s.source_name, s.target_name, s.confidence = m.source, m.target, m.confidence
        return [M(m) for m in self.cat.morphisms()]


def test_spectral(g: Category):
    print("\n" + "=" * 70)
    print("SPECTRAL  — does the Laplacian embedding cluster by ecomorph?")
    print("=" * 70)
    sa = SpectralGraphAnalyzer(category=g)
    sa.build_laplacian()
    evals, evecs = sa.compute_spectrum()
    print(f"  ran OK: {len(evals)} eigenvalues, "
          f"Fiedler value (algebraic connectivity) = {evals[1]:.4f}")

    # spectral embedding = first k non-trivial eigenvectors
    k = 4
    names = sa.node_names
    emb = {names[i]: evecs[i, 1:1 + k] for i in range(len(names))}

    def d(a, b):
        return float(np.linalg.norm(emb[a] - emb[b]))

    same_eco, diff_eco, same_isl = [], [], []
    for a, b in combinations(names, 2):
        ea, eb = SPECIES[a][1], SPECIES[b][1]
        ia, ib = SPECIES[a][0], SPECIES[b][0]
        if ea == eb:
            same_eco.append(d(a, b))
        else:
            diff_eco.append(d(a, b))
        if ia == ib and ea != eb:
            same_isl.append(d(a, b))
    print(f"  mean embedding distance:")
    print(f"     same ecomorph      : {np.mean(same_eco):.3f}")
    print(f"     diff ecomorph      : {np.mean(diff_eco):.3f}")
    print(f"     same island/diff eco: {np.mean(same_isl):.3f}")
    ok = np.mean(same_eco) < np.mean(diff_eco)
    print(f"  -> ecomorph-mates are spectrally {'CLOSER' if ok else 'NOT closer'} "
          f"than non-mates : {'PASS' if ok else 'FAIL'}")


def test_ricci(g: Category):
    print("\n" + "=" * 70)
    print("RICCI  — are the most negatively-curved edges convergent bridges?")
    print("=" * 70)
    orc = OllivierRicciCurvature(_RicciAdapter(g), alpha=0.5)
    edges = []
    for m in g.morphisms():
        kappa = orc.compute_edge_curvature(m.source, m.target)
        edges.append((kappa, m.source, m.target))
    edges.sort()  # most negative first

    def tag(a, b):
        ea, eb = SPECIES[a][1], SPECIES[b][1]
        ia, ib = SPECIES[a][0], SPECIES[b][0]
        same_eco = ea == eb
        cross_isl = ia != ib
        if same_eco and cross_isl:
            return "CONVERGENT bridge (same ecomorph, diff island)"
        if same_eco:
            return "same ecomorph, same island"
        return "cross-ecomorph"

    print("  Most negatively curved edges (candidate convergent bridges):")
    conv_in_top = 0
    topn = edges[:8]
    for kappa, a, b in topn:
        t = tag(a, b)
        if t.startswith("CONVERGENT"):
            conv_in_top += 1
        print(f"    kappa={kappa:+.3f}  {a:<20} ~ {b:<20} {t}")
    print(f"  -> {conv_in_top}/{len(topn)} most-negative edges are convergent bridges")


def lineage_community_graph() -> Category:
    """Graph where LINEAGES (islands) are the communities.

    - ancestry edges: every pair of island-mates connected (forms the community)
    - trait edges: cross-island pairs connected IF trait-similar (the bridges)
    A convergent pair = trait-similar across islands = a bridge between two
    lineage communities, which Ollivier-Ricci should flag as negatively curved.
    """
    morph, _ = build_lenses()
    g = Category(db_path=":memory:")
    species = list(SPECIES)
    for s in species:
        g.add(s)
    for a, b in combinations(species, 2):
        same_island = SPECIES[a][0] == SPECIES[b][0]
        sim = yoneda(morph, a, b)
        if same_island:
            g.connect(a, b, name=f"anc:{a}:{b}", confidence=1.0)        # community tie
        elif sim >= THRESH:
            g.connect(a, b, name=f"conv:{a}:{b}", confidence=float(sim))  # bridge candidate
    return g


def test_ricci_lineage(g: Category):
    print("\n" + "=" * 70)
    print("RICCI v2 — lineages as communities; convergent = cross-island bridge")
    print("=" * 70)
    orc = OllivierRicciCurvature(_RicciAdapter(g), alpha=0.5)
    edges = []
    for m in g.morphisms():
        kappa = orc.compute_edge_curvature(m.source, m.target)
        cross_island = SPECIES[m.source][0] != SPECIES[m.target][0]
        same_eco = SPECIES[m.source][1] == SPECIES[m.target][1]
        edges.append((kappa, m.source, m.target, cross_island, same_eco))
    edges.sort()
    print("  Most negatively curved edges:")
    hits = 0
    topn = edges[:8]
    for kappa, a, b, cross, eco in topn:
        kind = ("CONVERGENT (cross-island, same ecomorph)" if cross and eco
                else "cross-island, diff ecomorph" if cross else "within-island")
        if cross and eco:
            hits += 1
        print(f"    kappa={kappa:+.3f}  {a:<20} ~ {b:<20} {kind}")
    # are cross-island bridges more negative than within-island ties?
    cross_k = [e[0] for e in edges if e[3]]
    within_k = [e[0] for e in edges if not e[3]]
    print(f"\n  mean curvature  cross-island bridges : {np.mean(cross_k):+.3f}")
    print(f"  mean curvature  within-island ties   : {np.mean(within_k):+.3f}")
    ok = np.mean(cross_k) < np.mean(within_k)
    print(f"  -> convergent bridges are {'MORE negatively curved' if ok else 'NOT more curved'}"
          f" : {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    g = taxon_graph()
    print(f"Taxon graph: {len(g.objects())} taxa, {len(g.morphisms())} similarity edges "
          f"(threshold {THRESH})")
    test_spectral(g)
    test_ricci(g)
    g2 = lineage_community_graph()
    print(f"\nLineage-community graph: {len(g2.morphisms())} edges "
          f"(island ties + cross-island convergence bridges)")
    test_ricci_lineage(g2)
