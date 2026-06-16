"""Threshold sensitivity for the Anolis spectral/Ricci geometry claim."""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path
from statistics import mean

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from anolis_validation import SPECIES, build_lenses, yoneda  # noqa: E402
from komposos_core.core.category import Category  # noqa: E402
from komposos_wesys.geometry.grid_ricci import OllivierRicciCurvature  # noqa: E402
from komposos_wesys.geometry.grid_spectral import SpectralGraphAnalyzer  # noqa: E402


class _RicciAdapter:
    def __init__(self, cat: Category):
        self.cat = cat

    def list_morphisms(self, limit: int = 100000):
        class M:
            def __init__(self, m):
                self.source_name = m.source
                self.target_name = m.target
                self.confidence = m.confidence

        return [M(m) for m in self.cat.morphisms()]


def _taxon_graph(threshold: float) -> Category:
    morph, _ = build_lenses()
    g = Category(db_path=":memory:")
    for species in SPECIES:
        g.add(species)
    for a, b in combinations(SPECIES, 2):
        sim = yoneda(morph, a, b)
        if sim >= threshold:
            g.connect(a, b, name=f"sim:{a}:{b}", confidence=float(sim))
    return g


def _lineage_graph(threshold: float) -> Category:
    morph, _ = build_lenses()
    g = Category(db_path=":memory:")
    for species in SPECIES:
        g.add(species)
    for a, b in combinations(SPECIES, 2):
        same_island = SPECIES[a][0] == SPECIES[b][0]
        sim = yoneda(morph, a, b)
        if same_island:
            g.connect(a, b, name=f"anc:{a}:{b}", confidence=1.0)
        elif sim >= threshold:
            g.connect(a, b, name=f"conv:{a}:{b}", confidence=float(sim))
    return g


def _spectral_score(g: Category) -> dict[str, float | bool]:
    if len(list(g.morphisms())) == 0:
        return {"same_eco": float("nan"), "diff_eco": float("nan"), "pass": False}

    analyzer = SpectralGraphAnalyzer(category=g)
    analyzer.build_laplacian()
    evals, evecs = analyzer.compute_spectrum()
    names = analyzer.node_names
    k = min(4, max(1, len(evals) - 1))
    emb = {names[i]: evecs[i, 1:1 + k] for i in range(len(names))}

    same_eco: list[float] = []
    diff_eco: list[float] = []
    same_island_diff_eco: list[float] = []
    for a, b in combinations(names, 2):
        dist = float(np.linalg.norm(emb[a] - emb[b]))
        same_ecomorph = SPECIES[a][1] == SPECIES[b][1]
        same_island = SPECIES[a][0] == SPECIES[b][0]
        if same_ecomorph:
            same_eco.append(dist)
        else:
            diff_eco.append(dist)
        if same_island and not same_ecomorph:
            same_island_diff_eco.append(dist)

    se = mean(same_eco)
    de = mean(diff_eco)
    si = mean(same_island_diff_eco)
    return {
        "same_eco": round(se, 4),
        "diff_eco": round(de, 4),
        "same_island_diff_eco": round(si, 4),
        "fiedler": round(float(evals[1]), 4) if len(evals) > 1 else float("nan"),
        "pass": se < de and se < si,
    }


def _ricci_score(g: Category) -> dict[str, float | int | bool]:
    edges = list(g.morphisms())
    if not edges:
        return {
            "cross_mean": float("nan"),
            "within_mean": float("nan"),
            "top8_convergent": 0,
            "pass": False,
        }

    orc = OllivierRicciCurvature(_RicciAdapter(g), alpha=0.5)
    rows = []
    for edge in edges:
        cross = SPECIES[edge.source][0] != SPECIES[edge.target][0]
        same_eco = SPECIES[edge.source][1] == SPECIES[edge.target][1]
        rows.append(
            {
                "kappa": orc.compute_edge_curvature(edge.source, edge.target),
                "cross": cross,
                "same_eco": same_eco,
            }
        )

    cross_k = [r["kappa"] for r in rows if r["cross"]]
    within_k = [r["kappa"] for r in rows if not r["cross"]]
    top8 = sorted(rows, key=lambda r: r["kappa"])[:8]
    top8_conv = sum(1 for r in top8 if r["cross"] and r["same_eco"])
    ok = bool(cross_k and within_k and mean(cross_k) < mean(within_k))
    return {
        "cross_mean": round(mean(cross_k), 4) if cross_k else float("nan"),
        "within_mean": round(mean(within_k), 4) if within_k else float("nan"),
        "top8_convergent": top8_conv,
        "pass": ok,
    }


def main() -> int:
    thresholds = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    print("=" * 88)
    print("ANOLIS GEOMETRY THRESHOLD SENSITIVITY")
    print("=" * 88)
    print(
        "threshold  sim_edges  spectral_pass  sameEco  diffEco  "
        "ricci_pass  crossK  withinK  top8conv"
    )
    print("-" * 88)

    spectral_passes = 0
    ricci_passes = 0
    usable = 0
    for threshold in thresholds:
        sim_graph = _taxon_graph(threshold)
        lineage_graph = _lineage_graph(threshold)
        spectral = _spectral_score(sim_graph)
        ricci = _ricci_score(lineage_graph)
        usable += 1
        spectral_passes += int(bool(spectral["pass"]))
        ricci_passes += int(bool(ricci["pass"]))
        print(
            f"{threshold:>9.2f}"
            f"{len(list(sim_graph.morphisms())):>11}"
            f"{str(spectral['pass']):>15}"
            f"{spectral['same_eco']:>9}"
            f"{spectral['diff_eco']:>9}"
            f"{str(ricci['pass']):>12}"
            f"{ricci['cross_mean']:>8}"
            f"{ricci['within_mean']:>9}"
            f"{ricci['top8_convergent']:>10}"
        )

    print("\nHONEST READING")
    print(
        f"  Spectral ecomorph-closeness passed {spectral_passes}/{usable} "
        "thresholds."
    )
    print(
        f"  Ricci cross-island bridge negativity passed {ricci_passes}/{usable} "
        "thresholds."
    )
    print("  A robust geometry claim should not depend on exactly threshold 0.30.")
    return 0 if spectral_passes and ricci_passes else 1


if __name__ == "__main__":
    raise SystemExit(main())

