# Methods at a glance

Every method, what it actually computes, where it lives, and why it cannot be
tuned toward a wanted answer. Read alongside `BOOK_OF_FINDINGS.md` (the plain-
language version).

## Data provenance (read this first)

The analyses fall into two classes, and conclusions must be weighted accordingly.

- **Constructed control — Anolis (`anolis_validation.py`, `horns_convergence.py`,
  `geometry_test.py`, and everything importing `SPECIES`/`SYNDROME`).** Trait
  profiles are *not* measured morphometrics. Each species' states are coded from
  the canonical ecomorph syndromes (Losos 2009), with one jittered trait dimension
  per species and a shared per-island "scale" feature added as a competing ancestry
  signal. The ecomorph label is held out of the similarity computation but the
  traits are derived from it, so high Anolis AUROC values are a **sanity check that
  the pipeline recovers convergence when it is present by construction**, not a
  discovery from independent data.
- **Real data — Phenoscape fish + dated mammals (`phenoscape_*.py`,
  `trait_environment.py`).** States come from Phenoscape ontology queries; trees
  from Open Tree (TNRS) and the dated Fish Tree of Life / Upham 2019 chronograms.
  These carry the dissertation's empirical weight, and their results are
  appropriately tiered or negative.

## The detectors (sound the alarm)

| Method | What it computes | Precise technique | File | Key result |
|---|---|---|---|---|
| Yoneda similarity | how much two taxa share relationships, minus shared ancestry | weighted Jaccard of target-keyed relational fingerprints (computed in the analysis layer; see foundational fixes below) | `convergence_engine.py` | Anolis *control*: niche 0.54 vs ancestry-only 0.18 |
| Spectral clustering | global partition of taxa into convergent classes | graph Laplacian eigen-embedding | `convergence_engine.py`, `geometry_test.py` | clusters anoles by niche not island |
| Ollivier-Ricci curvature | which edges are cross-lineage bridges | Wasserstein-1 transport between neighbor distributions | `convergence_engine.py` (via `komposos_wesys`) | convergent edges −0.54 vs within-lineage +0.30 |
| Homotopy route class | convergent (different roads) vs parallel (same road) | shared-spine / interior-overlap of trajectories | `convergence_engine.classify_routes` | bat~whale = independent routes |

## Predict & verify

| Method | What | Technique | File | Result |
|---|---|---|---|---|
| Horn-filling | predict environment-driven traits | inner-horn (Λ²₁) completion of the category nerve; confidence = product along spine | `horns_convergence.py` | leave-one-out AUROC 0.94 vs 0.65 baseline |
| Retrodiction harness | recover deleted truth | leave-one-(taxon,trait)-out | `horns_convergence.py` | the AUROC above |
| Dual ZFC/CAT engine | two independent judges | ZFC logical entailment + CAT structural verdict → AGREE/ORPHAN/HOLLOW/REJECT | `horns_verified.py` (+ `src/.../zfc/bridge.py`) | known → AGREE, predicted → HOLLOW |
| CategoricalVerifier | discriminating structural confidence | curvature + neighborhood + connectivity + enriched weight | `categorical_verifier.py` | real preds 0.66 vs spurious 0.14 |

## The gauntlet (rule out the four impostors)

| Impostor | Test | Technique | File | Un-fittable because |
|---|---|---|---|---|
| Old news (plesiomorphy) | polarity weighting | Fitch parsimony root-state; ancestral states weighted ~0 | `polarity.py`, `polarity_fold.py` | parsimony has zero parameters |
| Inherited loss (homology) | MRCA independence | shared state counts only if the pair's MRCA lacked it | `polarity.homology_filter`, folded into `yoneda_sim` | pure tree topology + Fitch |
| Never had it (inapplicability) | Dollo recoding | applicable clade = MRCA of all taxa with the structure; outsiders → inapplicable | `polarity.dollo_recode`, `phenoscape_inapplicability.py` | MRCA-of-present is pure topology; no taxonomy injected |
| Dumb luck (chance) | randomization / PTP | permute which tips carry each state; recompute parsimony steps | `polarity.randomization_test`, `phenoscape_significance.py` | only parameter is #permutations |
| Dumb luck (drift) | neutral-Mk excess homoplasy | ML rate under null + parametric bootstrap on real dated branch lengths | `mk.py`, `phenoscape_excess*.py` | rate estimated under the no-convergence null |

## Environment as driver

| Method | What | Technique | File | Result |
|---|---|---|---|---|
| Trait-environment test | does the trait track the environment beyond phylogeny? | simulate trait under neutral Mk; compare observed phi(trait,env) to null | `mk.trait_environment_test`, `trait_environment.py` | Bergmann NOT supported (phi −0.07, p 0.78) |

## Real dated trees (the gold-standard null)

| Source | What | File |
|---|---|---|
| Fish Tree of Life (Rabosky 2018, treePL) | 11,638 dated ray-finned-fish tips | `chronogram.load_chronogram` |
| Upham 2019 mammal phylogeny (PHYLACINE) | 4,253 dated mammal tips | `chronogram.load_mammal_chronogram` |

Both chronograms are built from molecular clocks + fossils, **independent of the
trait data**, so using them as the evolutionary null cannot rig the convergence
result. Branch lengths are normalized to unit root-to-tip depth (the absolute
timescale is absorbed into the Mk rate; relative divergence structure is what informs
the null).

## The two foundational fixes (inherited code was wrong)

| Issue | Fix |
|---|---|
| Built-in `OptimusEngine.yoneda_similarity` returns 0 for all cross-object comparisons (its snapshot mangles duplicate morphism names, so no two objects share a key) | the analysis scripts compute a faithful Yoneda overlap keyed on the **target/source object** (what Yoneda actually says). The reported numbers use this faithful version, not the engine's built-in method, which remains unfixed in `core/optimus.py`. |
| Dead geometry bridge / unreachable HoTT branch / legacy `zfc` import paths | wired the real implementations; documented in `src/.../core` and `komposos_wesys` |
