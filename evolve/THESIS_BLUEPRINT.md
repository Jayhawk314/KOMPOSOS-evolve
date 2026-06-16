# Doctoral Thesis Blueprint

## Working title

**Convergence as a Tiered Inference Problem in Adaptive Relational Systems**

## Core thesis

Convergent evolution should not be treated as a single similarity score. It is a
tiered inference problem: a candidate similarity becomes a defensible convergence
claim only after it survives increasingly severe tests for ancestry, polarity,
homology, applicability, chance, neutral drift, and environmental association.

The categorical and complex-systems machinery is useful when it represents the
adaptive structure of the problem: organisms occupy niches, niches impose
constraints, traits recur as solutions, and missing trait relations can be
predicted by constraint propagation. It is not a substitute for phylogenetic
null models. The dissertation should defend that division of labor.

## One-sentence contribution

This work builds and tests a convergence engine that detects repeated biological
solutions, predicts missing adaptive traits, and then downgrades or rejects its
own claims under independent evolutionary nulls.

## Dissertation shape

### Chapter 1: The problem of repeated solutions

Question: When two lineages resemble each other, what would make that resemblance
evidence of repeated adaptation rather than inherited similarity?

Argument:

- Convergence is independent similarity, not mere similarity.
- Homoplasy can be adaptive, neutral, labile, or mis-coded.
- A serious convergence claim must be staged, not binary.

Repo evidence:

- `evolve/echolocation_flagship.py`
- `evolve/anolis_validation.py`

External positioning:

- Convergent evolution is commonly defined as independently evolved similarity
  in distinct lineages.
- Phylogenetic comparative methods exist because species are not independent
  samples.

### Chapter 2: Relational detection of candidate convergence

Question: On a constructed positive control, can a relational model recover
convergence without being handed the convergence label at scoring time?

Argument:

- Taxa can be represented by relational fingerprints over traits, niches, and
  ancestry.
- Similarity must be corrected by lineage structure.
- In the Anolis control, ecomorph similarity is recovered from syndrome-coded
  trait relations (the traits are derived from the ecomorph label, so this is a
  sanity check, not a discovery — see Risks and limits).

Key result:

- Same-ecomorph, different-island Anolis pairs: morphology similarity `0.540`.
- Different-ecomorph, same-island pairs: morphology similarity `0.180`.

Repo evidence:

- `evolve/convergence_engine.py`
- `evolve/anolis_validation.py`
- `evolve/polarity_fold.py`
- `evolve/homology_fold.py`

Claim strength:

- Strong as a positive validation on a known biological system.
- Not yet sufficient as a general convergence test without the later nulls.

### Chapter 3: Adaptive geometry and emergent structure

Question: Does the relational system show complex adaptive structure beyond
pairwise scores?

Argument:

- Convergence should appear as emergent organization in a graph of lineages,
  niches, and traits.
- Spectral structure should cluster organisms by adaptive role when niche is the
  organizing pressure.
- Ricci curvature should expose cross-lineage bridges between otherwise separate
  communities.

Key results:

- Anolis spectral embedding places same-ecomorph taxa closer than different
  ecomorph taxa.
- Lineage-community Ricci test: cross-island convergence bridges have mean
  curvature about `-0.537`; within-island ties about `+0.296`.

Repo evidence:

- `evolve/geometry_test.py`
- `evolve/categorical_verifier.py`
- `src/komposos_wesys/geometry/grid_ricci.py`
- `src/komposos_wesys/geometry/grid_spectral.py`

Claim strength:

- This is the best complex adaptive systems chapter.
- The language should be "emergent relational structure," not "geometry proves
  adaptation."

### Chapter 4: Horn filling as adaptive prediction

Question: If niches impose constraints, can the system predict missing trait
relations?

Argument:

- A species-niche-trait triangle can be treated as a partially filled adaptive
  relation.
- Inner horn filling proposes trait states selected by a niche and observed in
  cross-lineage niche-mates.
- Retrodiction tests whether the prediction recovers hidden known traits.

Key results:

- Leave-one-(species, trait)-out AUROC: `0.9391`.
- Frequency baseline AUROC: `0.6491`.
- Horn prediction beats the baseline by about `0.2900`.

Repo evidence:

- `evolve/horns_convergence.py`
- `evolve/horns_verified.py`

Claim strength:

- Strong within the Anolis validation setting.
- Needs a real-data extension before being claimed as a general biological
  predictor.

### Chapter 5: Verification as agreement and disagreement

Question: Can independent structural checks separate genuine adaptive support
from spurious relation completion?

Argument:

- Truth and convergence support are not identical.
- A trait can be true for a species but unsupported as a convergence claim.
- Dual verification should distinguish known, predicted, orphan, hollow, and
  rejected relations.

Key results:

- Known filled horns split into `AGREE` and `ORPHAN`.
- Unfilled predictions appear as structurally supported `HOLLOW` candidates.
- Real niche-routed predictions average CAT confidence about `0.659`.
- Spurious candidates average CAT confidence about `0.139`.

Repo evidence:

- `evolve/horns_verified.py`
- `evolve/categorical_verifier.py`
- `src/komposos_core/zfc/bridge.py`

Claim strength:

- Useful as ranking and filtering.
- Should not be presented as a hard proof system for biology.

### Chapter 6: The fish gauntlet

Question: Which apparent fin-loss convergence claims survive real biological
controls?

Argument:

- Shared absence is not automatically convergence.
- The pipeline must rule out ancestral retention, inherited loss, never-had-it
  coding, random scatter, and neutral drift.
- A convergence claim should be tiered by the strongest test it survives.

Tier ladder:

- Tier 0: no convergence; single origin or no meaningful signal.
- Tier 1: homoplasy present, but random or labile.
- Tier 2: structured homoplasy; survives topology and randomization.
- Tier 3: excess beyond neutral drift or strong phylogeny-controlled
  environment association.

Key results:

- Phenoscape candidate view flags pelvic, pectoral, dorsal, and caudal fin
  absences as apparent cross-order losses.
- Polarity test: pelvic fin, pectoral fin, and caudal fin show multiple origins;
  dorsal fin is single-origin in the tested set.
- Pelvic fin homology test: `152` independent pairs versus `148` shared-by-descent
  pairs among `300` pairs sharing absence.
- PTP/randomization: pelvic fin and pectoral fin are structured signals; dorsal
  and caudal are weak or chance-like in this run.
- Dated neutral-Mk: pelvic fin remains Tier 2, not Tier 3; structured homoplasy
  is real but not excess beyond neutral drift.

Repo evidence:

- `evolve/phenoscape_convergence.py`
- `evolve/phenoscape_polarity.py`
- `evolve/phenoscape_homology.py`
- `evolve/phenoscape_inapplicability.py`
- `evolve/phenoscape_significance.py`
- `evolve/phenoscape_excess.py`
- `evolve/phenoscape_excess_dated.py`
- `evolve/chronogram.py`

Claim strength:

- This is the honesty chapter.
- It should emphasize downgrading claims, not maximizing discoveries.

### Chapter 7: Environment as driver, with a famous negative result

Question: Does a trait track environment beyond phylogeny?

Argument:

- Process-based convergence needs an environmental driver, not only repeated
  state transitions.
- The test must control for phylogeny using a dated tree.
- Negative results are evidence of calibration, not failure.

Key result:

- Bergmann's rule across sampled mammals is not supported in this run:
  observed phi about `-0.067`, neutral-Mk p about `0.777`.
- Negative and positive controls behave as expected.

Repo evidence:

- `evolve/trait_environment.py`
- `evolve/mk.py`
- `evolve/chronogram.py`
- cached PanTHERIA and mammal chronogram inputs

Claim strength:

- Strong as a demonstration that the framework can reject a famous hypothesis.
- The dissertation should avoid claiming Bergmann's rule is false in all scopes;
  the claim is narrower: this operational cross-mammal test does not support it.

### Chapter 8: What the framework is, and what it is not

Question: What role do category theory and complex adaptive systems actually play?

Argument:

- Category-style representation is a disciplined way to encode relations,
  compositions, constraints, and missing completions.
- Complex adaptive systems language is justified when the system shows local
  constraints, modularity, emergent structure, and adaptive recurrence.
- Statistical phylogenetics remains necessary for inference.

The defensible boundary:

- Use relational/categorical tools for representation, detection, prediction,
  and structural filtering.
- Use phylogenetic methods for ancestry, nulls, rates, branch lengths,
  applicability, and significance.

## Research questions and executable answers

| Research question | Script evidence | Current answer |
|---|---|---|
| Can the engine recover convergence on a positive control? | `anolis_validation.py` | Yes, on the syndrome-coded Anolis control (sanity check, not discovery). |
| Can it separate convergence from kinship? | `echolocation_flagship.py` | Yes, in the demonstration case. |
| Does adaptive structure appear geometrically? | `geometry_test.py` | Yes, spectral and Ricci signals align with ecomorph bridges. |
| Is the geometry claim threshold-fragile? | `geometry_sensitivity.py` | No in the tested range: spectral and Ricci checks pass `7/7` thresholds. |
| Can niche constraints predict held-out traits? | `horns_convergence.py` | Yes, AUROC about `0.94` in Anolis. |
| Can structural verification reject spurious completions? | `horns_verified.py` | Partly; strong separation, not a hard gate. |
| Do thesis methods beat simple baselines? | `comparison_baselines.py` | Yes in this run: ancestry-corrected Anolis AUROC `0.9942`; horn AUROC beats nearest-neighbor by `0.0257`. |
| Do fish fin-loss claims survive polarity and homology? | `phenoscape_polarity.py`, `phenoscape_homology.py` | Pelvic and pectoral are strongest. |
| Can real Phenoscape states be relation-completed beyond Anolis? | `phenoscape_prediction.py` | Partly: morphology context predicts states, but nearest-profile beats cross-order/ancestry-subtracted variants. |
| Do they exceed random placement? | `phenoscape_significance.py` | Pelvic and pectoral do. |
| Do they exceed neutral drift on dated branches? | `phenoscape_excess_dated.py` | No; strongest claims remain Tier 2. |
| Does Bergmann's rule hold in this cross-mammal test? | `trait_environment.py` | No support in this operational test. |
| Where can tree-only inference fail? | `reticulate_caveat.py` | True convergence and reticulate transfer can produce identical tree-only outputs. |

## Literature positioning

The dissertation should explicitly compare against these families of methods:

- Pattern-based convergence measures: phenotypic distance corrected by
  phylogenetic distance, convergence strength, and frequency/surprise measures.
- Process-based methods: Ornstein-Uhlenbeck adaptive-regime methods such as
  SURFACE.
- Parsimony and ancestral-state methods: Fitch parsimony, consistency index,
  MRCA independence, Dollo-style applicability recoding.
- Phylogenetic null models: randomization/PTP, Brownian or Mk models, dated-tree
  simulations.
- Complex adaptive systems: networks of interacting adaptive agents, emergent
  organization, modularity, robustness, and evolvability.

The comparative claim should be modest and useful:

This framework complements existing phylogenetic comparative methods by making
the convergence claim explicit as a staged evidential object. It does not replace
OU models, Brownian models, Mk models, or biological judgment about character
coding.

## Risks and limits

- Tree assumption: most analyses assume bifurcating trees. Hybridization,
  introgression, and reticulate histories can mislead tree-based comparative
  methods.
- Character coding: absence can mean loss, ancestral absence, inapplicability,
  or missing data. The inapplicability chapter must stay central.
- Coverage: dated chronograms exclude many Phenoscape taxa, extinct taxa, higher
  taxa, or unmatched names.
- Constructed control, not discovery: the Anolis trait profiles are coded from
  the canonical ecomorph syndromes (Losos 2009), so the morphology features are a
  near-deterministic function of the ecomorph label. The echolocation demo is
  similarly illustrative. The thesis must mark these as constructed positive
  controls that verify the pipeline recovers convergence when it is present, not
  as discoveries from independent data or broad empirical surveys. High Anolis
  AUROC largely reflects reconstruction of the syndrome table, not biological
  signal. The empirical claims rest on the real-data Phenoscape and mammal work.
- Statistical power: Tier 3 negative results can reflect genuine lack of excess,
  conservative nulls, sparse data, or imperfect matching.
- Categorical overreach: the language is powerful but must not outrun the
  evidence. Biological claims require biological nulls.

## Execution plan

### Phase 1: Freeze the evidential baseline

Deliverable: a reproducible table of all current script outputs.

Actions:

- Add a single runner that executes the core thesis scripts.
- Capture outputs under `evolve/results/`.
- Record cache state and input data provenance.

Core scripts:

- `polarity.py`
- `mk.py`
- `anolis_validation.py`
- `echolocation_flagship.py`
- `geometry_test.py`
- `geometry_sensitivity.py`
- `horns_convergence.py`
- `horns_verified.py`
- `comparison_baselines.py`
- `phenoscape_convergence.py`
- `phenoscape_polarity.py`
- `phenoscape_homology.py`
- `phenoscape_inapplicability.py`
- `phenoscape_prediction.py`
- `phenoscape_significance.py`
- `phenoscape_excess.py`
- `phenoscape_excess_dated.py`
- `trait_environment.py`
- `reticulate_caveat.py`

### Phase 2: Add comparison baselines

Deliverable: a chapter table comparing this framework to conventional measures.

Actions:

- Add simple distance-versus-phylogeny convergence baselines where possible.
- Compare horn ranking against trait-frequency and nearest-neighbor baselines.
- Add sensitivity to threshold choices in graph construction.

Status:

- First baseline script added: `comparison_baselines.py`.
- Geometry threshold sensitivity added: `geometry_sensitivity.py`.
- Current Anolis detection comparison:
  - morphology-only AUROC: `0.9879`
  - ancestry-only AUROC: `0.3750`
  - not-same-island AUROC: `0.6250`
  - morphology-minus-ancestry AUROC: `0.9942`
- Current horn retrodiction comparison:
  - nearest-neighbor AUROC: `0.9134`
  - trait-frequency AUROC: `0.6491`
  - horn-max AUROC: `0.9391`
  - horn gain versus nearest-neighbor: `+0.0257`
  - horn gain versus trait-frequency: `+0.2900`
- Current geometry sensitivity:
  - thresholds tested: `0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50`
  - spectral ecomorph-closeness passed: `7/7`
  - Ricci cross-island bridge negativity passed: `7/7`

### Phase 3: Strengthen real-data prediction

Deliverable: a real-data horn or relation-completion experiment beyond Anolis.

Actions:

- Use Phenoscape characters as held-out morphology states.
- Predict held-out fin state from phylogenetically independent relatives and
  ontology neighborhood.
- Score recovery against frequency and phylogenetic-nearest baselines.

Status:

- First real-data relation-completion script added: `phenoscape_prediction.py`.
- Current all-state prediction:
  - frequency AUROC: `0.7661`
  - nearest-profile AUROC: `0.8867`
  - profile-minus-ancestry AUROC: `0.8110`
  - cross-order profile AUROC: `0.5949`
  - cross-order-minus-ancestry AUROC: `0.7499`
- Current absence/loss-only prediction:
  - frequency AUROC: `0.7023`
  - nearest-profile AUROC: `0.8526`
  - profile-minus-ancestry AUROC: `0.7912`
  - cross-order profile AUROC: `0.7559`
  - cross-order-minus-ancestry AUROC: `0.7909`
- Honest interpretation: Phenoscape morphology context predicts real held-out
  state relations, but this is not yet a general cross-lineage adaptive horn
  result. Without explicit niche/environment intermediates, the simple nearest
  profile remains the strongest predictor.

### Phase 4: Reticulate-history caveat

Deliverable: a limitations or simulation appendix.

Actions:

- Simulate a small tree plus reticulate transfer events.
- Show where the current MRCA/Mk pipeline can misclassify signal.
- State the future extension: phylogenetic networks instead of trees.

Status:

- Caveat script added: `reticulate_caveat.py`.
- Simulation design:
  - same observed bifurcating tree,
  - same observed tip states for the first two scenarios,
  - different unobserved histories.
- Current result:
  - `true_convergence`: two biological innovations, no transfer; tree-only
    reading reports `2` derived origins and `convergent=True`.
  - `reticulate_transfer`: one biological innovation plus one transfer; tree-only
    reading still reports `2` derived origins and `convergent=True`.
  - the two tree-only outputs are identical because the pipeline receives the
    same tree and same tips.
- Honest interpretation: MRCA/Fitch/Mk inference is conditional on the species
  tree being the relevant inheritance history for the trait. If hybridization,
  introgression, symbiosis, or horizontal transfer are plausible, convergence
  claims require a phylogenetic network or an explicit reticulate event model.

### Phase 5: Dissertation manuscript

Deliverable: full thesis draft.

Order:

1. Introduction and definitions.
2. Related work and comparative methods.
3. Relational convergence engine.
4. Adaptive geometry and horn prediction.
5. Verification and tiered inference.
6. Fish real-data gauntlet.
7. Mammal environment test and negative controls.
8. Limits, reticulation, and future work.

## Minimum defensible dissertation claim

The framework provides a reproducible, tiered method for turning apparent
similarity into a graded convergence claim. It succeeds on known positive
systems, rejects or downgrades weaker real-data claims, and clarifies where
relational/categorical structure is useful versus where statistical phylogenetic
nulls are indispensable.

## Stronger claim, if Phase 3 succeeds

Relational horn filling can predict missing adaptive trait states from
cross-lineage niche or ontology structure, and those predictions can be ranked
and filtered by independent structural and phylogenetic tests.

## Sources to cite in the thesis

- Convergent evolution and homoplasy definitions:
  https://en.wikipedia.org/wiki/Convergent_evolution and
  https://en.wikipedia.org/wiki/Homoplasy
- Phylogenetic comparative methods overview:
  https://en.wikipedia.org/wiki/Phylogenetic_comparative_methods
- Stayton, C. T. 2015. The definition, recognition, and interpretation of
  convergent evolution, and two new measures for quantifying and assessing the
  significance of convergence.
- Ingram, T. and Mahler, D. L. 2013. SURFACE: detecting convergent evolution
  from comparative data by fitting Ornstein-Uhlenbeck models.
- Arbuckle, K., Bennett, C. M., and Speed, M. P. 2014. A simple measure of the
  strength of convergent evolution.
- Adaptive radiation background:
  https://en.wikipedia.org/wiki/Adaptive_radiation
- Complex adaptive systems background:
  https://en.wikipedia.org/wiki/Complex_adaptive_system
- Evolvability background:
  https://en.wikipedia.org/wiki/Evolvability
- Reticulate-history caution:
  https://arxiv.org/abs/2603.25986
