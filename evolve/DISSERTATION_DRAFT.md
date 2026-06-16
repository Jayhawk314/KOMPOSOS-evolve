# Convergence as a Tiered Inference Problem in Adaptive Relational Systems

Working dissertation draft, version 0.1.

Source of truth for numeric results:
`evolve/results/20260613T225636Z_phase4_network/`

That archive was produced by `evolve/run_thesis_baseline.py` and contains 19
passing scripts, with no failed or timed-out analyses.

## Abstract

Convergent evolution is often introduced as independent lineages arriving at the
same biological solution. The practical difficulty is that similarity alone is
not convergence. Similarity can be inherited from a common ancestor, retained
from an ancestral state, produced by random lability, created by neutral drift,
or misread from a character state that was never applicable to one of the
lineages being compared. This dissertation treats convergence as a tiered
inference problem rather than a single score.

The system developed here represents taxa, traits, niches, environments, and
ancestry as typed relations. This relational layer detects candidate repeated
solutions, exposes adaptive structure through graph geometry, and predicts
missing trait relations through niche-mediated horn filling. The biological
claim is then passed through a gauntlet of independent checks: polarity-aware
Fitch parsimony, MRCA-based homology filtering, Dollo-style applicability
recoding, phylogenetic randomization, neutral Mk simulation on dated trees, and
phylogenetically controlled trait-environment tests.

The framework is first exercised on a constructed positive control. In a
Caribbean Anolis lizard model whose trait profiles are coded from the canonical
ecomorph syndromes (Losos 2009), same-ecomorph pairs on different islands have
mean morphology similarity 0.540, while different-ecomorph pairs on the same
island have mean morphology similarity 0.180. A morphology-minus-ancestry score
separates same-ecomorph cross-island pairs with AUROC 0.9942, and niche-mediated
horn filling recovers held-out trait states with AUROC 0.9391, beating trait
frequency by 0.2900 and a nearest-profile baseline by 0.0257. Because the trait
profiles are generated from the ecomorph definitions, these numbers verify that
the pipeline behaves correctly when convergence is genuinely present; they are a
sanity check on the machinery, not a discovery from measured morphology. The
empirical weight of the dissertation therefore rests on the real-data chapters
that follow.

On those real data, the framework refuses stronger claims when the evidence does
not support them. On Phenoscape fish anatomy, pelvic fin loss survives polarity,
homology,
and randomization filters, but the dated neutral-Mk test leaves it at Tier 2:
structured homoplasy, not excess convergence beyond drift. A mammal-wide
Bergmann test on a dated chronogram is not supported in this operational test
(phi = -0.067, p = 0.777), while positive and negative controls behave as
expected. A reticulate-history simulation shows the core tree-based limitation:
true convergence and one innovation plus transfer can produce identical
tree-only outputs.

The contribution is therefore not a universal convergence detector. It is a
reproducible evidential framework for turning apparent similarity into a graded
convergence claim, while preserving the ability to downgrade or reject that
claim.

## Central Claim

A convergence claim is not a single number. It is a staged evidential object.

The stages are:

1. Candidate similarity: two lineages share a trait relation.
2. Ancestry correction: the similarity is not explained by close relatedness.
3. Polarity: the shared state is derived, not ancestral retention.
4. Homology: the shared state was not inherited from the pair's MRCA.
5. Applicability: the state is a true loss or gain, not a never-applicable code.
6. Significance: the pattern is not random scatter over the tree.
7. Neutral null: the pattern exceeds, or fails to exceed, neutral drift.
8. Environmental/process support: the trait tracks an external adaptive driver.
9. Historical caveat: the tree is the relevant inheritance history for the trait.

The relational and categorical machinery helps with stages 1, 2, 3 in folded
form, 4 in folded form, and prediction. The phylogenetic and statistical
machinery is responsible for stages 3 through 9. The dissertation should keep
that division clean.

## Tier Ladder

Tier 0: no convergence claim.

The trait has a single origin, lacks meaningful signal, or cannot be scored with
enough coverage.

Tier 1: homoplasy present, but weak.

The state appears more than once on the tree, but the pattern may be random,
labile, or too sparsely sampled.

Tier 2: structured homoplasy.

The state survives polarity, homology, applicability, and randomization checks.
It is real structure on the phylogeny, but not necessarily excess beyond neutral
drift.

Tier 3: excess or process-supported convergence.

The pattern exceeds a neutral model on dated branches, or a trait-environment
test supports an adaptive driver beyond phylogeny.

## Evidence Table

| Question | Script | Result | Interpretation |
|---|---|---:|---|
| Does the engine recover convergence when it is present? | `anolis_validation.py` | same-ecomorph cross-island morphology 0.540 vs same-island different-ecomorph 0.180 | Constructed control (syndrome-coded traits): the engine groups by niche, not ancestry. |
| Does ancestry correction help (control)? | `comparison_baselines.py` | morph-only AUROC 0.9879; morph-minus-ancestry AUROC 0.9942 | On the control, subtracting ancestry sharpens an already strong signal. |
| Does geometry depend on one threshold (control)? | `geometry_sensitivity.py` | spectral 7/7; Ricci 7/7 | The control geometry claim is not fragile over tested thresholds. |
| Can niche relations predict traits (control)? | `horns_convergence.py` | horn-max AUROC 0.9391 | On the control, niche-mediated horn filling recovers held-out traits. |
| Does horn prediction beat baselines (control)? | `comparison_baselines.py` | +0.0257 over nearest-profile; +0.2900 over frequency | Horns add modest information over a strong profile baseline. |
| Does structural verification filter spurious predictions (control)? | `horns_verified.py` | CAT prediction mean 0.659 vs spurious 0.139 | Verification ranks and filters; it is not a hard proof. |
| Do fish fin losses survive polarity? | `phenoscape_polarity.py` | pelvic 4 origins; pectoral 2; dorsal 1; caudal 2 | Pelvic and pectoral are strongest among fin characters. |
| Is pelvic fin loss independent or inherited? | `phenoscape_homology.py` | 152 independent pairs; 148 shared-by-descent pairs | MRCA filtering splits real recurrence from inherited loss. |
| Does dated neutral drift explain fish losses? | `phenoscape_excess_dated.py` | pelvic Tier 2; adipose Tier 2; barbel Tier 1 | Strongest fish claims are structured, not Tier 3. |
| Does morphology context predict Phenoscape states? | `phenoscape_prediction.py` | nearest-profile AUROC 0.8867 all states; 0.8526 absence states | Real-data prediction works, but simple profile similarity wins. |
| Does the mammal Bergmann test support association? | `trait_environment.py` | phi -0.067; p 0.777 | Not supported in this cross-mammal operational test. |
| Where can tree-only inference fail? | `reticulate_caveat.py` | transfer and convergence give identical tree outputs | Tree-based claims are conditional on tree-like inheritance. |

## Chapter 1: Why Similarity Is Not Enough

The biological idea is simple: lineages sometimes arrive independently at the
same solution. The methodological problem is that the data rarely say
"independently." The data say that two or more taxa share a state. The rest is
inference.

This distinction is the starting point for the dissertation. A shark and a
dolphin are streamlined for reasons that are not reducible to shared descent.
Two close relatives may also resemble each other, but that resemblance is not
convergence unless the relevant state arose independently. A fish that lacks a
structure may have lost it, inherited the absence, never had the structure in
the first place, or simply lack a usable observation. The same surface pattern
can therefore arise from several different histories.

This work treats convergence as a claim that must survive a sequence of
questions. Does the similarity remain after subtracting ancestry? Is the shared
state derived? Did the pair's MRCA already have it? Is the character applicable
to both taxa? Is the pattern more structured than random placement? Is it more
than neutral drift? Does it track an environmental driver? Is the species tree
the right history for this trait?

The value of the system is that these questions are executable. They are not
only cautions in prose. Each is represented by a script, an output, and a
reported boundary.

## Chapter 2: Relational Detection of Candidate Convergence

The relational model represents taxa by what they map to: traits, niches,
environments, molecular states, and ancestry. A taxon's relational fingerprint is
compared to another taxon's fingerprint. Candidate convergence is high trait
similarity after subtracting ancestry similarity.

The Anolis system is used here as a constructed positive control, and its status
must be stated plainly. The trait profiles are not measured morphometrics: each
species' trait states are coded from the canonical ecomorph syndromes (Losos
2009), perturbed by a single jittered trait dimension per species and a shared
per-island "scale" feature that injects a competing ancestry signal. The
morphology lens sees only those trait states; the ancestry lens sees only island
membership; the ecomorph label is held out of the similarity computation and used
only at scoring time. Because the trait states are themselves derived from the
ecomorph label, this is not a discovery from independent data. It is a check that
the relational layer recovers convergence when convergence is, by construction,
present — and that an injected ancestry signal does not overwhelm it.

Read that way, the control behaves as it should. Same-ecomorph pairs on different
islands have mean morphology similarity 0.540; different-ecomorph pairs on the
same island have mean morphology similarity 0.180. The morphology-only AUROC for
the same-ecomorph cross-island label is 0.9879, and the morphology-minus-ancestry
score is 0.9942: ancestry correction sharpens the signal rather than weakening
it.

Interpretation: the relational layer groups taxa by ecological role rather than
by the injected island/ancestry feature, which is the behaviour a convergence
detector should show on a positive control.

Boundary: this is a sanity check on synthetic, syndrome-coded traits, not
evidence about nature. High AUROC here largely reflects that the morphology lens
reconstructs the syndrome table it was given. No biological convergence claim
rests on this chapter; the empirical claims live in the real-data fish and mammal
chapters, which supply the filters that decide how far a claim can be carried.

## Chapter 3: Adaptive Geometry

If convergence is more than pairwise resemblance, it should appear as structure
in a graph. The Anolis system supplies that graph: taxa are nodes, and edges
encode trait similarity. Spectral structure asks whether taxa occupy nearby
positions in the graph embedding when they share an ecomorph. Ricci curvature
asks whether cross-island same-ecomorph edges behave like bridges between
lineage communities.

At the default threshold, spectral embedding places same-ecomorph taxa much
closer than different-ecomorph taxa. The sensitivity test then sweeps thresholds
from 0.20 to 0.50. Spectral ecomorph-closeness passes 7 out of 7 thresholds.
Ricci cross-island bridge negativity also passes 7 out of 7 thresholds.

Interpretation: the geometry claim is not an artifact of selecting exactly one
edge threshold. The same adaptive organization is visible across the tested
range.

Boundary: graph geometry is a structural signature, not a biological null
model. It helps identify candidate organization. It does not replace polarity,
homology, or drift tests. This chapter still runs on the constructed Anolis
control of Chapter 2, so it shows that the geometry is recoverable when
convergence is present, not that it was discovered in measured data.

## Chapter 4: Horn Filling as Adaptive Prediction

The prediction chapter asks whether niche constraints can recover held-out
traits. The construction is a triangle:

`Species -> Niche -> Trait`

If a species inhabits a niche, and that niche selects a trait, but the direct
species-to-trait relation is missing, the unfilled triangle becomes a prediction.
In the Anolis validation, the system contains 50 objects and 212 morphisms. It
finds 229 inner horns from species to niche to trait, with 133 filled and 96
unfilled.

The retrodiction harness removes known species-trait relations and asks whether
the horn score recovers them. There are 114 held-out positives and 304
negatives. Horn-max AUROC is 0.9391. The global trait-frequency baseline is
0.6491. The nearest-profile baseline is much stronger at 0.9134, but horn-max
still beats it by 0.0257.

Interpretation: niche-mediated relation completion adds information beyond
frequency and slightly beyond raw trait-profile proximity in the Anolis setting.

Boundary: the gain over nearest-profile is modest, and it is measured on the
constructed Anolis control, where niche, trait, and ecomorph are coupled by
construction. The defensible claim is that explicit niche-mediated structure
improves a strong profile baseline in this controlled setting. The real-data
analogue is Chapter 8, where no niche intermediates exist and the gain does not
hold.

## Chapter 5: Verification and Filtering

The verification layer separates truth from convergence support. A trait may be
true for a species but not supported as a convergence claim. Conversely, an
unfilled horn may be structurally plausible but not yet known as a true trait.

In the archived run, filled known horns split into 107 AGREE and 26 ORPHAN.
Unfilled predictions are 96 HOLLOW. The filter test is more informative:
niche-routed predictions have mean CAT confidence 0.659, while spurious
candidates average 0.139. Of 57 spurious candidates, 48 are rejected.

Interpretation: the verifier is useful as a ranking and filtering layer. It
distinguishes niche-routed convergence candidates from many structurally weak
candidates.

Boundary: it is not a hard proof system for biology, and it is demonstrated on
the constructed Anolis control. Some spurious candidates survive, and some true
traits read as ORPHAN. The correct role is filtering, not final biological
adjudication.

## Chapter 6: The Fish Gauntlet

The Phenoscape fish analyses are the honesty chapter. They begin with apparent
fin-loss convergence and then ask how far that claim can survive.

The polarity analysis scores 123 cleaned taxa across four fin characters, with
108 names resolved by Open Tree. Pelvic fin has 65 scored tips, 8 steps, CI
0.12, ancestral state present, and 4 derived origins. Pectoral fin has 2
derived origins. Dorsal fin has only 1 derived origin in this scoring and is not
called convergent. Caudal fin has 2 derived origins.

The homology filter for pelvic fin absence places 25 absent taxa. Among 300
pairs sharing absence, 152 are classified as independent and 148 as
shared-by-descent. This split matters: without the MRCA test, inherited loss and
independent recurrence would be mixed together.

The dated neutral-Mk test is the strongest downgrading step. It uses a real fish
chronogram with 11,638 dated tips. Pelvic fin remains Tier 2: 30 matched tips,
4 origins, PTP p_fewer 0.015, and Mk observed/null 4/4.58 with p_excess 0.738.
Adipose fin is also Tier 2. Barbel is Tier 1. Caudal fin is Tier 0 in the dated
run.

Interpretation: pelvic fin loss is structured homoplasy. It survives several
filters and is not dismissed as random scatter.

Boundary: it is not Tier 3. The dated neutral-Mk null does not support excess
homoplasy beyond drift. The system therefore downgrades the claim instead of
overstating it.

## Chapter 7: Environment as Driver

The mammal environment test asks whether a trait tracks environment beyond
phylogeny. The operational test is Bergmann's rule across a sampled mammal set:
large body size versus cold climate, both median-split, evaluated on a dated
mammal chronogram.

The archived run uses 600 dated tips. Counts are balanced: large-cold 140,
large-warm 160, small-cold 160, small-warm 140. The observed phi is -0.067. The
neutral-Mk null mean is -0.003, and p_assoc is 0.777. The negative control is
non-significant. The positive control is significant, with phi 0.617 and p
0.001.

Interpretation: this operational mammal-wide Bergmann test is not supported.
The controls show that the test can reject noise and detect a strong injected
association.

Boundary: this does not prove that Bergmann's rule is false in every taxonomic
or ecological scope. It says that this pre-specified broad test does not support
the association beyond phylogeny.

## Chapter 8: Real-Data Prediction Boundary

The Phenoscape prediction experiment asks whether morphology states can be
relation-completed beyond the Anolis validation. The result is useful because it
does not give the strongest hoped-for answer.

For all observed fin states, nearest-profile AUROC is 0.8867. Frequency is
0.7661. Profile-minus-ancestry is 0.8110. Cross-order profile is 0.5949.
Cross-order-minus-ancestry is 0.7499. For absence/loss states only,
nearest-profile remains best at 0.8526. Cross-order-minus-ancestry reaches
0.7909.

Interpretation: morphology context predicts held-out Phenoscape state relations
well.

Boundary: without explicit niche or environment intermediates, the simple
nearest-profile baseline remains strongest. The real-data relation-completion
result should be presented as a boundary on the current representation, not as a
general horn-prediction success.

## Chapter 9: Reticulate Histories

The tree-based gauntlet assumes that the species tree is the relevant inheritance
history for the trait. That assumption can fail. Hybridization, introgression,
symbiosis, and horizontal transfer can move traits across the tree.

The reticulate caveat makes the failure mode executable. Two scenarios are given
the same observed tree and the same tip states. In the first, a derived state
arises independently on A1 and C1. In the second, it arises once on A1 and is
transferred or introgressed to C1. The tree-only pipeline reports identical
outputs for both: 2 derived origins and convergent=True.

Interpretation: the tree-only layer is doing what it can with the information it
receives. The missing information is the non-tree edge.

Boundary: any convergence claim from this system is conditional on a tree-like
inheritance model. Where reticulate histories are plausible, the future method
must use a phylogenetic network or explicit reticulate event model.

## Chapter 10: Conclusion

The dissertation's contribution is a disciplined convergence workflow. The
system detects repeated trait structure, predicts missing adaptive relations,
and then submits its own claims to increasingly severe tests. On the constructed
Anolis control, the pipeline behaves as a convergence detector should — it
recovers the planted signal relationally, geometrically, and predictively — which
establishes that the machinery works when convergence is present. The substantive
findings appear in the real-data fish and mammal analyses, where apparent
convergence is honestly limited to Tier 2 or rejected as unsupported.

The final claim should stay precise:

This framework turns similarity into a graded convergence claim. It behaves
correctly on a constructed positive control, finds structured homoplasy in real
fish anatomy, rejects an unsupported broad Bergmann test, and exposes its own
tree-based limits under reticulate history.

The next extension is not more confidence. It is richer history: explicit
environment nodes in real data, stronger external trait coding, and phylogenetic
networks where the inheritance history is not tree-like.

## Manuscript To-Do

- Add a formal related-work chapter with full citations and DOI checks.
- Convert the evidence table into manuscript tables.
- Add methods details for each script without duplicating the code.
- Decide whether to include all result archives or only the official Phase 4
  archive.
- Add figures:
  - tier ladder,
  - Anolis relational diagram,
  - fish gauntlet flow,
  - reticulate caveat diagram.
- Tighten the prose for chapter transitions.
- Add a reproducibility appendix pointing to `run_thesis_baseline.py`.

