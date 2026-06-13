#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
CategoricalVerifier -- CAT structural meta-verifier (ported from PHARM).

Ported from KOMPOSOS-IV-PHARM/oracle/categorical_verifier.py, with:
  * the legacy `from core.category import Category` dropped (we duck-type),
  * Ricci wired to komposos_wesys.geometry via the _RicciAdapter (the wesys
    OllivierRicciCurvature takes a store, not a Category),
  * a bug fix: the original never set verdict.mean_curvature, so the geometric
    class was always UNKNOWN -- we set it from the path curvature here.

CAT verifies a (source, target, relation) claim structurally, independently of
ZFC, across up to four dimensions that DISCRIMINATE among claims (unlike the
dual-engine's flat path-count fallback):

  1. curvature    : Ricci along the connecting path. Spherical (kappa>0) = tight
                    cluster = strong; hyperbolic (kappa<0) = fragile bridge.
  2. neighborhood : Jaccard overlap of source/target neighborhoods (Kan-like).
  3. stability    : edge-connectivity (min edges to disconnect) = persistence.
  4. enriched     : morphism confidence weight.

Returns a StructuralVerdict with structural_confidence in [0,1] + geometric_class.
Inject into the dual engine via `bridge._verifier = CategoricalVerifier(cat)`.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from convergence_engine import _RicciAdapter  # noqa: E402

try:
    from komposos_wesys.geometry.grid_ricci import OllivierRicciCurvature
    RICCI_AVAILABLE = True
except Exception:
    RICCI_AVAILABLE = False


class GeometricClass(Enum):
    SPHERICAL = auto()    # kappa > 0: tight cluster, high confidence
    EUCLIDEAN = auto()    # kappa ~ 0: chain/pathway, moderate confidence
    HYPERBOLIC = auto()   # kappa < 0: bridge between clusters, fragile
    UNKNOWN = auto()


@dataclass
class StructuralVerdict:
    source: str
    target: str
    relation: str
    structural_confidence: float = 0.0
    geometric_class: GeometricClass = GeometricClass.UNKNOWN
    curvature_score: float = 0.0
    neighborhood_score: float = 0.0
    stability_score: float = 0.0
    enriched_score: float = 0.0
    path_count: int = 0
    path_lengths: List[int] = field(default_factory=list)
    mean_curvature: float = 0.0
    explanation: str = ""
    dimensions_used: List[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return (f"StructuralVerdict({self.source}->{self.target} "
                f"conf={self.structural_confidence:.3f} geo={self.geometric_class.name})")


class CategoricalVerifier:
    """CAT structural verifier over a komposos Category (duck-typed)."""

    def __init__(self, category):
        self.category = category
        self._graph: Optional[nx.Graph] = None
        self._curvature_computed = False
        self._edge_curvatures: Dict[Tuple[str, str], float] = {}

    # -- structural view ----------------------------------------------------
    def _ensure_graph(self) -> nx.Graph:
        if self._graph is not None:
            return self._graph
        G = nx.Graph()
        for obj in self.category.objects():
            G.add_node(obj.name, type_name=getattr(obj, "type_name", None),
                       metadata=getattr(obj, "metadata", {}))
        for m in self.category.morphisms():
            G.add_edge(m.source, m.target, relation=m.name,
                       confidence=m.confidence, metadata=getattr(m, "metadata", {}))
        self._graph = G
        return G

    def _ensure_curvature(self) -> None:
        if self._curvature_computed or not RICCI_AVAILABLE:
            self._curvature_computed = True
            return
        G = self._ensure_graph()
        if G.number_of_edges() == 0:
            self._curvature_computed = True
            return
        try:
            ricci = OllivierRicciCurvature(_RicciAdapter(self.category))
            result = ricci.compute_all_curvatures()
            for (src, tgt), kappa in result.edge_curvatures.items():
                self._edge_curvatures[(src, tgt)] = kappa
                self._edge_curvatures[(tgt, src)] = kappa
        except Exception:
            pass
        self._curvature_computed = True

    # -- core verification --------------------------------------------------
    def verify(self, source: str, target: str, relation: str) -> StructuralVerdict:
        v = StructuralVerdict(source=source, target=target, relation=relation)
        G = self._ensure_graph()
        if source not in G or target not in G:
            v.explanation = "source/target not in structural graph"
            return v

        paths, lengths = self._find_paths(G, source, target)
        v.path_count = len(paths)
        v.path_lengths = lengths

        cscore = self._curvature_check(G, source, target, paths, v)
        if cscore is not None:
            v.curvature_score = cscore
            v.dimensions_used.append("curvature")

        v.neighborhood_score = self._neighborhood_check(G, source, target)
        v.dimensions_used.append("neighborhood")

        sscore = self._stability_check(G, source, target)
        if sscore is not None:
            v.stability_score = sscore
            v.dimensions_used.append("stability")

        escore = self._enriched_check(source, target)
        if escore is not None:
            v.enriched_score = escore
            v.dimensions_used.append("enriched")

        scores = []
        if "curvature" in v.dimensions_used:
            scores.append(v.curvature_score)
        scores.append(v.neighborhood_score)
        if "stability" in v.dimensions_used:
            scores.append(v.stability_score)
        if "enriched" in v.dimensions_used:
            scores.append(v.enriched_score)

        path_bonus = min(0.3, len(paths) * 0.1) if paths else 0.0
        if scores:
            v.structural_confidence = min(1.0, sum(scores) / len(scores) + path_bonus)
        elif paths:
            v.structural_confidence = path_bonus

        v.geometric_class = self._classify_geometry(v.mean_curvature)
        v.explanation = self._build_explanation(v)
        return v

    # -- dimensions ---------------------------------------------------------
    def _curvature_check(self, G, source, target, paths, verdict) -> Optional[float]:
        if not RICCI_AVAILABLE:
            return None
        self._ensure_curvature()
        if not self._edge_curvatures:
            return None
        curvatures = []
        for path in paths:
            for i in range(len(path) - 1):
                e = (path[i], path[i + 1])
                kappa = self._edge_curvatures.get(e, self._edge_curvatures.get((e[1], e[0]), 0.0))
                curvatures.append(kappa)
        if not curvatures:
            kappa = self._edge_curvatures.get((source, target),
                                              self._edge_curvatures.get((target, source), None))
            if kappa is not None:
                curvatures.append(kappa)
        if not curvatures:
            return None
        mean_kappa = sum(curvatures) / len(curvatures)
        verdict.mean_curvature = mean_kappa          # <-- bug fix vs original
        score = 0.5 + mean_kappa * 0.5               # [-1,1] -> [0,1]
        return max(0.0, min(1.0, score))

    def _neighborhood_check(self, G, source, target) -> float:
        sn = set(G.neighbors(source)) if source in G else set()
        tn = set(G.neighbors(target)) if target in G else set()
        if not sn and not tn:
            return 0.0
        union = sn | tn
        if not union:
            return 0.0
        jaccard = len(sn & tn) / len(union)
        if G.has_edge(source, target):
            return min(1.0, 0.5 + jaccard * 0.5)
        return jaccard * 0.7

    def _stability_check(self, G, source, target) -> Optional[float]:
        if source not in G or target not in G:
            return None
        try:
            connectivity = nx.edge_connectivity(G, source, target)
            return min(1.0, connectivity * 0.2)
        except (nx.NetworkXError, nx.NetworkXUnfeasible):
            return None

    def _enriched_check(self, source, target) -> Optional[float]:
        ms = [m for m in self.category.morphisms_from(source) if m.target == target]
        if not ms:
            ms = [m for m in self.category.morphisms_from(target) if m.target == source]
        if not ms:
            return None
        return sum(m.confidence for m in ms) / len(ms)

    # -- helpers ------------------------------------------------------------
    def _find_paths(self, G, source, target, max_length: int = 5):
        try:
            paths = list(nx.all_simple_paths(G, source, target, cutoff=max_length))
            return paths, [len(p) - 1 for p in paths]
        except (nx.NetworkXError, nx.NodeNotFound):
            return [], []

    def _classify_geometry(self, mean_curvature: float) -> GeometricClass:
        if mean_curvature > 0.1:
            return GeometricClass.SPHERICAL
        elif mean_curvature < -0.1:
            return GeometricClass.HYPERBOLIC
        elif mean_curvature != 0.0:
            return GeometricClass.EUCLIDEAN
        return GeometricClass.UNKNOWN

    def _build_explanation(self, v: StructuralVerdict) -> str:
        parts = []
        parts.append("No structural paths found" if v.path_count == 0
                     else f"{v.path_count} path(s), lengths {v.path_lengths}")
        geo = v.geometric_class
        if geo == GeometricClass.SPHERICAL:
            parts.append("Spherical (tight cluster, high support)")
        elif geo == GeometricClass.HYPERBOLIC:
            parts.append("Hyperbolic (bridge, fragile)")
        elif geo == GeometricClass.EUCLIDEAN:
            parts.append("Euclidean (chain, moderate support)")
        if v.neighborhood_score > 0.5:
            parts.append(f"strong neighborhood overlap ({v.neighborhood_score:.2f})")
        if v.stability_score > 0.6:
            parts.append(f"topologically stable ({v.stability_score:.2f})")
        return "; ".join(parts)


class ConvergenceCATVerifier:
    """CAT verifier specialised for convergence, on the TAXON-TAXON geometry.

    The generic CategoricalVerifier fails to filter on the hub-dense Species/Niche/
    Trait nerve (everything is 2-3 hops apart). The right structural question for a
    prediction "species X should have trait T" is convergence-shaped:

        is there a CROSS-LINEAGE, TRAIT-SIMILAR donor species that has T?

    Real prediction  -> a different-lineage look-alike has T (that IS convergence)
                        -> high structural confidence.
    Spurious (T not in X's syndrome) -> only DISSIMILAR donors have T
                        -> low structural confidence -> REJECT.

    structural_confidence = max over donors Y-of-T of  conv(X,Y),
        conv(X,Y) = trait_jaccard(X,Y) - ancestry(X,Y)   (cross-lineage similarity)
    This is the Yoneda/Ricci taxon-taxon signal we showed separates convergent pairs
    from non-convergent ones (ricci_bridges: -0.54 vs +0.30).
    """

    def __init__(self, species_traits: dict, species_lineage: dict):
        self.traits = species_traits            # sp -> set of trait-states
        self.lineage = species_lineage          # sp -> lineage tag (island/clade)

    def _syndrome(self, sp: str) -> set:
        return {t for t in self.traits.get(sp, ()) if not t.startswith("anc:")}

    def _trait_jaccard(self, a: str, b: str) -> float:
        A, B = self._syndrome(a), self._syndrome(b)
        if not (A | B):
            return 0.0
        return len(A & B) / len(A | B)

    def verify(self, source: str, target: str, relation: str) -> StructuralVerdict:
        v = StructuralVerdict(source=source, target=target, relation=relation)
        if source not in self.traits:
            v.explanation = "unknown taxon"
            return v
        donors = [y for y in self.traits
                  if y != source and target in self.traits[y]]
        v.path_count = len(donors)
        if not donors:
            v.explanation = "no donor taxon carries this trait"
            return v
        best, best_y, best_cross = 0.0, None, False
        for y in donors:
            sim = self._trait_jaccard(source, y)
            anc = 1.0 if self.lineage.get(source) == self.lineage.get(y) else 0.0
            conv = max(0.0, sim - anc)          # cross-lineage trait similarity
            if conv > best:
                best, best_y, best_cross = conv, y, (anc == 0.0)
        v.structural_confidence = best
        v.neighborhood_score = best
        # convergent support across lineages reads as a tight cross-lineage cluster
        v.geometric_class = (GeometricClass.SPHERICAL if best > 0.3 and best_cross
                             else GeometricClass.HYPERBOLIC if best > 0
                             else GeometricClass.UNKNOWN)
        v.explanation = (f"{len(donors)} donor(s); best cross-lineage convergence "
                         f"{best:.2f} via {best_y}")
        return v
