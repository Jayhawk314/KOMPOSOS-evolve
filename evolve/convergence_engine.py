# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
A small convergent-evolution engine built on the komposos categorical core.

Idea (honest version)
---------------------
Every taxon, trait, behavior, molecular state, environment and clade is an
OBJECT. Relationships are MORPHISMS, grouped into named LENSES:

    morphology : taxon --has--> morph-trait-state
    behavior   : taxon --does--> behavior-state
    molecular  : taxon --carries--> substitution/gene-state
    environment: taxon --inhabits--> niche          (the DRIVER)
    ancestry   : taxon --member-of--> clade         (to be subtracted out)

A taxon's Yoneda fingerprint in a lens = its profile of relationships in that
lens. Two taxa are similar-in-a-lens when they map to the same targets.

    convergence(A,B | lens) = yoneda_sim(A,B | lens) - yoneda_sim(A,B | ancestry)

high similarity in a trait lens, NOT explained by shared ancestry.

DEEP convergence = convergence is high in morphology AND behavior AND molecular
at once, for a pair whose ancestry similarity is low. That stacking is the
signature the categorical/relational view makes legible.

We compute the Yoneda fingerprint keyed on the TARGET object (Hom(X,-) as a
function of where X maps), because the repo's built-in metric keys on morphism
name and its snapshot mangles duplicate names to 0. See komposos-repo-reality.
"""

from __future__ import annotations

import sys
from pathlib import Path
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from komposos_core.core.category import Category  # noqa: E402


class _RicciAdapter:
    """Expose a komposos Category through the interface OllivierRicciCurvature
    expects: list_morphisms(limit) -> objs with .source_name/.target_name/.confidence.
    (komposos_wesys.geometry lives in a sibling package; see komposos-repo-reality.)
    """

    def __init__(self, cat: Category):
        self.cat = cat

    def list_morphisms(self, limit: int = 100000):
        class M:
            def __init__(s, m):
                s.source_name, s.target_name, s.confidence = m.source, m.target, m.confidence
        return [M(m) for m in self.cat.morphisms()]


def fingerprint(cat: Category, obj: str) -> dict:
    """Correct Yoneda fingerprint keyed on connected objects (not morphism name)."""
    return {
        "out": {m.target: m.confidence for m in cat.morphisms_from(obj)},
        "in": {m.source: m.confidence for m in cat.morphisms_to(obj)},
    }


def yoneda_sim(cat: Category, a: str, b: str,
               weight: "dict[str, float] | None" = None,
               homology: "callable | None" = None) -> float:
    """Weighted Jaccard overlap of two objects' Yoneda fingerprints in `cat`.

    POLARITY (`weight`): each connected object k contributes a factor
    weight.get(k, 1.0). Ancestral/plesiomorphic states get ~0 so sharing them is
    not convergence evidence; derived states keep 1.

    HOMOLOGY (`homology`): an optional predicate homology(a, b, k) -> bool. For a
    state k that BOTH a and b have (a shared state), if homology(a,b,k) is False
    the state is SHARED BY DESCENT (their MRCA already had it = one inherited
    event, not convergence) and is EXCLUDED entirely (neutral: out of both
    numerator and denominator). States only one taxon has are unaffected.

    Both default to None -> original equal-weight similarity (backward-compatible).
    """
    fa, fb = fingerprint(cat, a), fingerprint(cat, b)
    keys = set(fa["out"]) | set(fb["out"]) | set(fa["in"]) | set(fb["in"])
    if not keys:
        return 1.0 if a == b else 0.0
    wf = (lambda k: 1.0) if weight is None else (lambda k: weight.get(k, 1.0))
    shared = total = 0.0
    for k in keys:
        out_min = min(fa["out"].get(k, 0), fb["out"].get(k, 0))
        in_min = min(fa["in"].get(k, 0), fb["in"].get(k, 0))
        # a shared-by-descent state (fails the homology test) is excluded as neutral
        if homology is not None and (out_min > 0 or in_min > 0) and not homology(a, b, k):
            continue
        w = wf(k)
        shared += w * (out_min + in_min)
        total += w * (max(fa["out"].get(k, 0), fb["out"].get(k, 0)) +
                      max(fa["in"].get(k, 0), fb["in"].get(k, 0)))
    return shared / total if total else 0.0


@dataclass
class ConvergenceModel:
    """Multi-lens categorical model of a set of taxa.

    Each lens is an independent Category. `ancestry` is the lens we subtract to
    isolate genuine convergence (similarity not due to shared descent).
    """

    taxa: list[str] = field(default_factory=list)
    lenses: dict[str, Category] = field(default_factory=dict)
    # Optional POLARITY weights: trait-state -> [0,1]. ~0 = ancestral (sharing it is
    # not convergence evidence), 1 = derived. Set via set_polarity(); when empty all
    # states weigh 1 (original behaviour). Applies to trait lenses, never ancestry.
    state_weight: dict = field(default_factory=dict)
    # Optional HOMOLOGY mask: homology(taxonA, taxonB, state) -> bool. False means
    # the shared state is inherited (MRCA had it) and is excluded from similarity.
    # Set via set_homology(); pairwise, so it is a callable not a dict. See polarity.TreeHomology.
    homology: object = None

    def set_homology(self, predicate) -> None:
        """Fold MRCA-based homology into trait similarity: a shared state counts
        only if independently derived in both taxa (their MRCA lacked it). The
        predicate should be backed by Fitch reconstruction on an independent tree
        (polarity.TreeHomology) so it cannot be fitted."""
        self.homology = predicate

    def set_polarity(self, state_weight: dict) -> None:
        """Fold ancestral-state polarity into every trait-similarity computation.

        `state_weight` maps a trait-state object -> weight in [0,1]. It should come
        from a parameter-free reconstruction on an INDEPENDENT tree (e.g.
        polarity.analyze via Fitch) so this cannot be tuned toward a wanted answer.
        Down-weighting ancestral states makes the Yoneda/Ricci/spectral detectors
        (which all route through `similarity`) plesiomorphy-aware at once.
        """
        self.state_weight = dict(state_weight)

    def lens(self, name: str) -> Category:
        if name not in self.lenses:
            self.lenses[name] = Category(db_path=":memory:")
        return self.lenses[name]

    def add_taxon(self, name: str) -> None:
        if name not in self.taxa:
            self.taxa.append(name)
        for cat in self.lenses.values():
            cat.add(name)

    def relate(self, lens: str, taxon: str, target: str, rel: str, conf: float = 1.0) -> None:
        """taxon --rel--> target, recorded in the named lens."""
        cat = self.lens(lens)
        cat.add(taxon)
        cat.add(target)
        cat.connect(taxon, target, name=f"{rel}:{target}", confidence=conf)

    # ---- queries ----------------------------------------------------------
    def similarity(self, lens: str, a: str, b: str) -> float:
        # polarity weights + homology mask apply to trait lenses; ancestry stays raw
        if lens == "ancestry":
            return yoneda_sim(self.lens(lens), a, b)
        weight = self.state_weight if self.state_weight else None
        return yoneda_sim(self.lens(lens), a, b, weight=weight, homology=self.homology)

    def convergence(self, lens: str, a: str, b: str) -> float:
        """Similarity in `lens` minus similarity in the ancestry lens."""
        anc = self.similarity("ancestry", a, b) if "ancestry" in self.lenses else 0.0
        return self.similarity(lens, a, b) - anc

    def deep_convergence(self, a: str, b: str, trait_lenses: list[str]) -> dict:
        """Convergence across every trait lens at once, plus the ancestry it is
        measured against. 'deep' = min across lenses is still high while
        ancestry is low.
        """
        anc = self.similarity("ancestry", a, b) if "ancestry" in self.lenses else 0.0
        per_lens = {L: self.similarity(L, a, b) for L in trait_lenses}
        return {
            "pair": (a, b),
            "ancestry_sim": anc,
            "per_lens": per_lens,
            "min_trait_sim": min(per_lens.values()) if per_lens else 0.0,
            "deep_score": (min(per_lens.values()) if per_lens else 0.0) - anc,
        }

    def rank_deep(self, trait_lenses: list[str]) -> list[dict]:
        """All taxon pairs ranked by deep convergence (stacked across lenses)."""
        out = [self.deep_convergence(a, b, trait_lenses)
               for a, b in combinations(self.taxa, 2)]
        out.sort(key=lambda d: d["deep_score"], reverse=True)
        return out

    # ======================================================================
    # Geometric detectors (Ricci + spectral), wired from komposos_wesys.
    # These give INDEPENDENT signatures of convergence that triangulate with
    # the Yoneda/relational score above:
    #   * spectral : global partition of morphospace into convergent classes
    #   * ricci    : per-edge scalar; convergent pairs are negatively-curved
    #                bridges between lineage communities
    # ======================================================================
    def _trait_sim(self, a: str, b: str, trait_lenses: list[str]) -> float:
        """Mean Yoneda similarity across the trait lenses."""
        sims = [self.similarity(L, a, b) for L in trait_lenses]
        return sum(sims) / len(sims) if sims else 0.0

    def trait_similarity_graph(self, trait_lenses: list[str],
                               threshold: float = 0.3) -> Category:
        """Taxon-taxon graph; edge = mean trait similarity above `threshold`.
        Communities here are the convergent CLASSES (use for spectral)."""
        g = Category(db_path=":memory:")
        for t in self.taxa:
            g.add(t)
        for a, b in combinations(self.taxa, 2):
            s = self._trait_sim(a, b, trait_lenses)
            if s >= threshold:
                g.connect(a, b, name=f"sim:{a}:{b}", confidence=float(s))
        return g

    def lineage_community_graph(self, trait_lenses: list[str],
                                anc_threshold: float = 0.5,
                                trait_threshold: float = 0.3) -> Category:
        """Graph where LINEAGES are the communities (use for Ricci).
        - ancestry tie if ancestry similarity >= anc_threshold (binds a community)
        - else trait BRIDGE if mean trait similarity >= trait_threshold
        A convergent pair = trait-similar across lineages = a bridge."""
        g = Category(db_path=":memory:")
        for t in self.taxa:
            g.add(t)
        for a, b in combinations(self.taxa, 2):
            anc = self.similarity("ancestry", a, b) if "ancestry" in self.lenses else 0.0
            if anc >= anc_threshold:
                g.connect(a, b, name=f"anc:{a}:{b}", confidence=float(anc))
            else:
                s = self._trait_sim(a, b, trait_lenses)
                if s >= trait_threshold:
                    g.connect(a, b, name=f"bridge:{a}:{b}", confidence=float(s))
        return g

    def spectral_embedding(self, trait_lenses: list[str], k: int = 4,
                           threshold: float = 0.3) -> dict[str, np.ndarray]:
        """Laplacian spectral embedding of taxa from the trait-similarity graph.
        Returns taxon -> coordinates in the first k non-trivial eigenvectors."""
        from komposos_wesys.geometry.grid_spectral import SpectralGraphAnalyzer
        g = self.trait_similarity_graph(trait_lenses, threshold)
        sa = SpectralGraphAnalyzer(category=g)
        sa.build_laplacian()
        _, evecs = sa.compute_spectrum()
        names = sa.node_names
        kk = min(k, max(1, len(names) - 1))
        return {names[i]: evecs[i, 1:1 + kk] for i in range(len(names))}

    def spectral_convergence(self, a: str, b: str, trait_lenses: list[str],
                             k: int = 4, threshold: float = 0.3) -> float:
        """Spectral closeness of a,b (1 / (1+embedding distance)) minus ancestry.
        High => spectrally co-clustered despite low shared ancestry."""
        emb = self.spectral_embedding(trait_lenses, k, threshold)
        if a not in emb or b not in emb:
            return 0.0
        dist = float(np.linalg.norm(emb[a] - emb[b]))
        closeness = 1.0 / (1.0 + dist)
        anc = self.similarity("ancestry", a, b) if "ancestry" in self.lenses else 0.0
        return closeness - anc

    def ricci_bridges(self, trait_lenses: list[str], alpha: float = 0.5,
                      anc_threshold: float = 0.5, trait_threshold: float = 0.3) -> list[dict]:
        """Ollivier-Ricci curvature on the lineage-community graph.
        Returns edges sorted by curvature ASCENDING (most negative first); the
        most-negative are the strongest convergent bridges."""
        from komposos_wesys.geometry.grid_ricci import OllivierRicciCurvature
        g = self.lineage_community_graph(trait_lenses, anc_threshold, trait_threshold)
        orc = OllivierRicciCurvature(_RicciAdapter(g), alpha=alpha)
        out = []
        for m in g.morphisms():
            kappa = orc.compute_edge_curvature(m.source, m.target)
            is_bridge = m.name.startswith("bridge:")
            out.append({
                "pair": (m.source, m.target),
                "kappa": float(kappa),
                "is_bridge": is_bridge,
                "ancestry_sim": self.similarity("ancestry", m.source, m.target)
                if "ancestry" in self.lenses else 0.0,
                "trait_sim": self._trait_sim(m.source, m.target, trait_lenses),
            })
        out.sort(key=lambda d: d["kappa"])
        return out

    @staticmethod
    def classify_routes(paths: list[list[str]]) -> dict:
        """HoTT-flavored route analysis: do trajectories to a shared endpoint
        take the SAME route or INDEPENDENT routes?

        An evolutionary trajectory is a path (ancestor ... -> phenotype). Two
        lineages reaching the same phenotype form a loop (down path1, up path2).
        - loop CONTRACTIBLE (shared interior spine)  => PARALLEL evolution
        - loop NON-CONTRACTIBLE (disjoint interiors) => CONVERGENT evolution

        This is the distinction the other detectors (yoneda/ricci/spectral)
        CANNOT make -- they only see that the endpoints are similar, not whether
        the routes were independent. NOTE: the repo's PathHomotopyChecker calls
        everything homotopic (its DISTINCT branch is unreachable), so we compute
        route independence directly from interior overlap.

        Returns classification + the homotopy reading.
        """
        if len(paths) < 2:
            return {"classification": "trivial", "reason": "need >= 2 trajectories"}
        endpoints = {(p[0], p[-1]) for p in paths if p}
        interiors = [set(p[1:-1]) for p in paths]
        shared = set.intersection(*interiors) if interiors else set()
        pairwise_disjoint = all(
            not (interiors[i] & interiors[j])
            for i in range(len(interiors)) for j in range(i + 1, len(interiors))
        )
        same_endpoint = len({e[1] for e in endpoints}) == 1

        if not same_endpoint:
            cls, reason = "no convergence", "trajectories do not share an endpoint"
        elif all(iv == interiors[0] for iv in interiors):
            cls, reason = "PARALLEL", "identical interior route (contractible loop)"
        elif pairwise_disjoint:
            cls, reason = "CONVERGENT", "disjoint interior routes (non-contractible loop)"
        else:
            cls, reason = "PARTIALLY PARALLEL", f"shared spine {sorted(shared)} then diverge"
        return {
            "classification": cls,
            "reason": reason,
            "shared_endpoint": same_endpoint,
            "shared_spine": sorted(shared),
            "homotopic": cls in ("PARALLEL",),
        }

    def detect(self, a: str, b: str, trait_lenses: list[str]) -> dict:
        """All three convergence signatures for one pair, side by side."""
        ricci = {tuple(sorted(d["pair"])): d["kappa"] for d in self.ricci_bridges(trait_lenses)}
        return {
            "pair": (a, b),
            "yoneda_deep": self.deep_convergence(a, b, trait_lenses)["deep_score"],
            "spectral": self.spectral_convergence(a, b, trait_lenses),
            "ricci_kappa": ricci.get(tuple(sorted((a, b)))),  # None if no edge
            "ancestry_sim": self.similarity("ancestry", a, b) if "ancestry" in self.lenses else 0.0,
        }
