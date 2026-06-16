# Glossary — the jargon, in plain words

Terms are grouped by where they appear. Each definition is meant to be correct *and*
readable.

## The biology

- **Convergent evolution** — the same trait arising *independently* in lineages that
  did not inherit it from a shared ancestor. Sharks and dolphins both became torpedoes.
- **Parallel evolution** — reaching the same trait by the *same* evolutionary route
  (same intermediate steps). Cousins of convergence; we distinguish them with the
  "route" lens.
- **Homoplasy** — the umbrella term for similarity that is *not* due to common
  descent. Convergence and parallelism are both kinds of homoplasy.
- **Homology** — similarity that *is* due to common descent (inherited from a shared
  ancestor). The opposite of homoplasy. Your arm and a bat's wing are homologous.
- **Plesiomorphy** — an *ancestral* trait, inherited and old (e.g. "has a backbone"
  among vertebrates). Sharing one is no evidence of convergence. ("Old news.")
- **Derived trait (apomorphy)** — a *new* trait, evolved within the group. Sharing a
  derived trait independently *is* convergence. ("News.")
- **Lineage / clade** — a branch of the tree of life; a group of organisms and all
  descendants of their common ancestor.
- **Lability** — how easily a trait flips on and off over evolutionary time. Highly
  labile traits appear in scattered lineages by chance, mimicking convergence.
- **Dollo's law** — complex structures, once lost, are not re-evolved from scratch;
  evolution is effectively one-directional for them.
- **Bergmann's rule** — the hypothesis that within warm-blooded groups, body size
  increases in colder climates (bigger bodies conserve heat). We test it; across all
  mammals it is not supported.
- **Ecomorph** — a body-plan associated with a particular niche (e.g. "twig anole"),
  repeated across regions. Convergence made into a category.

## Trees & reconstruction

- **Phylogeny** — the family tree of species, showing who is related to whom.
- **Chronogram** — a phylogeny whose branch lengths are *time* (millions of years),
  calibrated by molecular clocks and fossils.
- **Branch length** — how much change (or time, in a chronogram) is on an edge of
  the tree.
- **MRCA** — Most Recent Common Ancestor; the node where two species' lineages meet.
- **Ancestral-state reconstruction** — inferring what trait an ancestor had, from
  the traits of its descendants and the tree.
- **Fitch parsimony** — the simplest reconstruction: assign ancestral states so as to
  minimize the total number of changes. Has **no tunable parameters** — which is why
  we trust it not to be rigged.
- **Consistency Index (CI)** — minimum possible changes ÷ observed changes on the
  tree. For a binary trait the minimum is 1, so CI = 1 ÷ (observed changes). CI = 1
  means a single origin (no convergence); CI < 1 means the trait arose more than once
  (homoplasy).
- **Inapplicable** — a character that doesn't apply to a lineage because the relevant
  structure never existed there (a jellyfish's "leg length" is inapplicable, not zero).

## The categorical / mathematical machinery

- **Category theory** — a branch of mathematics about objects and the relationships
  (arrows) between them, and how relationships compose.
- **Yoneda lemma** — the principle that an object is fully determined by its
  relationships to everything else. Basis of our "relational fingerprint."
- **Relational fingerprint** — the full list of what a species connects to (traits,
  environment, lineage). We compare fingerprints to measure convergence.
- **Ollivier-Ricci curvature** — a notion of curvature for networks. Negative
  curvature marks an edge that *bridges* two otherwise-separate communities, which we
  *interpret* as a candidate convergent pair across distant lineages (a structural
  signature, not proof of convergence).
- **Spectral clustering** — grouping points using the eigenvectors of a graph's
  Laplacian; finds natural clusters (here, convergent classes).
- **Homotopy** — in topology, whether two paths can be continuously deformed into
  each other. We use it to ask whether two lineages took the *same* road (parallel)
  or *different* roads (convergent) to the same trait.
- **Simplicial nerve / inner horn** — a structure built from composable chains
  A→B→C. An unfilled "inner horn" is a pattern missing its conclusion; *filling* it
  is our prediction step (predict A→C from A→B→C).
- **Dual ZFC/CAT engine** — two independent verifiers: ZFC (logic — is the claim
  entailed?) and CAT (structure — does the shape support it?). Their agreement
  classifies each claim (AGREE / ORPHAN / HOLLOW / REJECT).

## The statistics

- **AUROC** — Area Under the ROC Curve; a score from 0.5 (coin flip) to 1.0 (perfect)
  for how well a method ranks true cases above false ones. Our predictor scored 0.94
  on the constructed Anolis control; on real Phenoscape data the best predictor is
  lower (~0.89) and a simple profile baseline wins.
- **phi coefficient** — a correlation between two yes/no variables (e.g. "large body"
  vs "cold climate"). Positive = they co-occur.
- **PTP (Permutation Tail Probability)** — a randomization test: shuffle the data
  many times to build a "by chance" distribution, then see where reality falls.
- **Mk model** — a simple continuous-time model of a trait changing along a tree at
  some rate. Our null for "blind, goalless evolution."
- **Parametric bootstrap** — simulate data under a fitted null model many times to
  ask whether the real data is unusual.
- **Maximum likelihood (ML) estimate** — the parameter value that makes the observed
  data most probable under a model. We use it to set the Mk rate *under the no-
  convergence null* — estimated, never hand-tuned.
- **Null model** — a model of "what would happen with no special cause." A result is
  significant only if it stands out against the null.
- **Grafen branch lengths** — branch lengths assigned from tree *shape* alone, used
  when no dated tree is available. A parameter-free stand-in for real divergence times.
