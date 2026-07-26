# Future Directions — the convergence engine beyond biology

> Companion to [`viewer/EXPANSION_PLAN.md`](viewer/EXPANSION_PLAN.md). That document
> grows the **biological** tree-of-life coverage (more clades, environment, scale).
> This one asks the bigger question: the engine is not a biology tool, it is a
> **general detector of independent rediscovery** — "the same solution reached by
> disjoint paths, proven to exceed chance after subtracting shared history." This is
> the brainstorm of where that generality can go, kept to the same honesty contract
> (real data, ancestry subtracted, claims downgraded not inflated, nothing "done"
> until verified). It is a menu and a sequencing argument, not a promise of dates.

---

## The one reframing

The whole gauntlet needs exactly three ingredients:

| Ingredient | In biology | The general role |
|---|---|---|
| **Tree** | dated chronogram | the shared-history scaffold to subtract |
| **Traits** | morphology / ecology / sequence | states on the tips |
| **Null** | neutral Mk / PTP randomization | what drift / chance produces |

Anywhere those three exist, `polarity.py` (Fitch origins, PTP) + `mk.py`
(excess-homoplasy, trait-environment) + `convergence_engine.py` (Yoneda ancestry
subtraction, Ricci bridges, spectral morphospace) run unchanged. Every direction
below is "what else has a tree + traits + a null?"

---

## The menu (all directions at a glance)

| # | Direction | Tier | Effort | Risk | Reuses | Payoff |
|---|---|---|---|---|---|---|
| 1 | **Neural-net representational convergence** | adjacent | medium | low | ESM lens, ancestry subtraction, gauntlet | novel + current + fundable |
| 2 | **Real-time pathogen convergence early-warning** | adjacent | medium | medium | posterior integration, molecular lens | high public-health impact |
| 3 | **The Convergence Atlas of Life** | adjacent | high | medium | the entire biological pipeline | a public good; the product of Phase 5/6 |
| 4 | **Linguistic / cultural convergence** | new universe | medium | medium | Fitch/PTP/Mk, env-driver | a whole new domain, same code |
| 5 | **Invert it: convergence as design target** | new universe | high | high | ESM lens (run backwards) | detection → engineering |
| 6 | **The Inevitability Index (replay-the-tape)** | outside | medium | medium | origin counts, dated-Mk null | quantifies a famous open question |
| 7 | **Convergence-based biosignature / technosignature prior** | outside | low | high | the Atlas (#3) | a principled astrobiology prior |
| 8 | **Topological convergence (make the HoTT real)** | outside | high | high | ricci_bridges, topology_bridge, hott_bridge | convergence as an invariant |
| 9 | **Genotype↔phenotype convergence routing** | outside | high | medium | morphology + molecular lenses | the deepest scientific payoff |

Cross-cutting infrastructure (needed by several) is listed at the end.

---

## Tier 1 — adjacent, the machinery already fits

### 1. Neural-network representational convergence
**Claim under test:** independently-trained models converge to the *same* internal
representations (the "platonic representation hypothesis") — but that claim is made
*without subtracting shared history* (shared pretraining data, shared architecture
lineage). That is exactly the gap this engine closes.

- **Mapping:** models = taxa; the pretrain → fine-tune lineage = phylogeny; learned
  features / probe-recoverable concepts = traits; representational-similarity (CKA /
  SVCCA) = the lens distance.
- **First experiment:** take ~10–20 open checkpoints with a known training lineage
  (e.g. a base model and its descendants vs an independently-trained family). Build
  the "phylogeny" from documented lineage. For a battery of probes (does feature X
  exist?), run the gauntlet: is a shared feature present *beyond* shared training
  ancestry? Reuse `esmc_diagnose.py`'s structure verbatim (it already does
  embedding-distance − ancestry-distance + permutation, and already exposes the
  circularity trap).
- **Honesty caveats:** "ancestry" for models is messy (data overlap is the real
  confounder, not just lineage) — this must be stated and, where possible, measured.
  Representation similarity metrics are themselves contested; report more than one.
- **Exit criteria:** a verdict on whether a named feature is *convergent* across two
  model families after ancestry subtraction, with the same T0–T3 tiering.

### 2. Real-time pathogen convergence early-warning
**Claim under test:** resistance / immune-escape mutations recur convergently across
viral and bacterial lineages; the convergent ones are the dangerous ones.

- **Mapping:** strains = taxa (a growing posterior from GISAID/Nextstrain); mutations
  = traits; the dated phylogeny is rebuilt as data streams in.
- **First experiment:** one protein (e.g. flu HA or SARS-CoV-2 spike), a fixed set of
  candidate sites, the gauntlet across the posterior (`posterior_integration.py`
  already runs the whole battery across N trees). Flag sites whose independent origins
  exceed the neutral-Mk null — candidate convergent-escape positions — and
  retrospectively check against known escape mutations.
- **Honesty caveats:** sampling bias is severe (surveillance is non-random); the null
  must account for it. This is forecasting, so calibrate the false-positive rate
  end-to-end before any claim.
- **Exit criteria:** on held-out historical data, the engine flags known convergent
  escape sites earlier than they fixed, with a stated FP rate.

### 3. The Convergence Atlas of Life
**Claim under test:** none — this is the *product*. A queryable, FDR-controlled,
ancestry-subtracted, posterior-integrated verdict for **every** (clade × trait) the
data supports.

- **Mapping:** it is the biological pipeline (Phases 3–6 of `EXPANSION_PLAN.md`) run
  exhaustively and cached, not hand-picked.
- **First experiment:** ship the *vertebrate* slice first (mammals + birds + fish are
  already dated + coded): a static, browsable atlas of every coded trait's verdict
  distribution. `verify_all.py` is the standing guard; the data layer is Phase 3.
- **Honesty caveats:** trait deserts must be shown, not hidden; reticulation caveat
  for microbes; every cell carries its tier + q-value, never a bare "yes".
- **Exit criteria:** a public, reproducible atlas where every claim links to its
  gauntlet output.

---

## Tier 2 — same engine, different universe

### 4. Linguistic / cultural / technological convergence
**Claim under test:** features like tone, click consonants, SOV word order, or
inventions like agriculture / writing / metallurgy arose independently more than
drift predicts ("phylomemetics").

- **Mapping:** Glottolog language tree + WALS typological features is *structurally
  identical* to chronogram + traits. Cultures/technologies need a phylogeny proxy
  (descent-with-modification of traditions).
- **First experiment:** WALS feature × Glottolog tree, one feature (e.g. tone), the
  exact Fitch → PTP → (dated-)Mk gauntlet. `env_driver.py` ports directly if you
  treat geography/contact as the "environment" axis (areal diffusion is the
  linguistic analogue of HGT — and the same reticulation caveat applies).
- **Honesty caveats:** contact/borrowing breaks the tree model exactly as HGT does —
  this is the reticulation chapter in a new domain; refuse or switch to a network.
- **Exit criteria:** one typological convergence with a defensible tier, contact
  explicitly controlled for.

### 5. Invert the engine — convergence as a design target
**Claim under test (flipped):** instead of *detecting* that lineages fell into a
functional attractor, *aim* at the attractor on purpose.

- **Mapping:** the ESM functional-neighborhood lens, run backwards: locate the basin
  a convergent function occupies, then design sequences toward it.
- **First experiment:** for a convergence the engine already finds (or the Prestin
  per-site analysis), characterize the embedding-space attractor, then score *de
  novo* candidate sequences by distance to it; validate against known functional
  variants before any wet-lab claim.
- **Honesty caveats:** embedding proximity ≠ function (the whole Track-C lesson);
  in-silico only until validated. Do **not** overclaim design success from geometry.
- **Exit criteria:** designed candidates rank known-functional variants above
  known-nonfunctional ones, held out.

---

## Tier 3 — genuinely outside the box

### 6. The Inevitability Index — a number for "replay the tape of life"
**Claim under test:** how repeatable is each evolutionary solution? Gould said replay
the tape and you'd get a different world; Conway Morris said you'd get the same
attractors. This makes it measurable.

- **Mapping:** origins-per-opportunity, drift-controlled, per trait — the engine
  already computes independent origins and the neutral-Mk expectation; the index is
  observed-origins normalized by the null's expected origins, across the Atlas.
- **First experiment:** compute the index for every coded vertebrate trait already in
  the repo; rank from "deep attractor" (eyes-like: many origins beyond drift) to
  "frozen accident" (singular). The data exists *today* — this is mostly a new
  summary statistic over existing gauntlet output.
- **Honesty caveats:** "opportunity" is hard to define (denominator problem); report
  several denominators and the sensitivity across them (reuse the Idea-4 pattern).
- **Exit criteria:** a ranked, drift-controlled inevitability table with explicit
  denominator sensitivity.

### 7. Convergence-based biosignature / technosignature prior
**Claim under test:** traits with many independent origins on Earth are the safe bets
for life (or technology) *anywhere*.

- **Mapping:** the Inevitability Index (#6) over the Atlas (#3) *is* the prior — high
  origin-count traits are convergence-robust and least contingent on Earth history.
- **First experiment:** publish the top-of-the-index traits as an explicit
  "expect-elsewhere" list with the origin counts behind each; same idea for
  convergent technologies (#4) as a technosignature prior.
- **Honesty caveats:** Earth is n=1 for the *origin* of life; this is a prior over
  *outcomes given life*, not over life arising. Say so loudly.
- **Exit criteria:** a defensible, origin-count-grounded prior, clearly scoped.

### 8. Topological convergence — make the HoTT metaphor real
**Claim under test:** two lineages reaching the same phenotype by disjoint interior
paths form a *non-contractible loop* in evolutionary space. Right now (`CLAUDE.md`,
`hott_bridge.py`) that is poetry; this makes it an invariant.

- **Mapping:** build a morphospace from the trait-similarity graph; convergent pairs
  are the most negatively-curved bridges (already in `convergence_engine.ricci_bridges`).
  Persistent homology (`topology_bridge.py`) over morphospace-through-time asks
  whether convergent "holes" appear and persist.
- **First experiment:** for one well-understood convergence, compute the Ricci
  curvature field + a persistence diagram of morphospace; show convergent pairs as
  persistent topological features, not just low-distance pairs.
- **Honesty caveats:** this is the most speculative mathematically; an invariant is
  only worth it if it *predicts* something a plain distance does not. Hold it to that.
- **Exit criteria:** a topological signature that distinguishes true convergence from
  homology on held-out cases better than distance alone.

### 9. Genotype↔phenotype convergence routing
**Claim under test:** same phenotype — same genes (true molecular convergence) or
different genes (developmental system drift)? Classify *how* each convergence is
built.

- **Mapping:** triangulate the morphology lens (phenotype convergent?) against the
  molecular lens (`esmc_diagnose.py`: same protein route?). Echolocation/Prestin is
  the proof-of-concept already in the repo.
- **First experiment:** for a handful of classic convergences with sequenced genomes,
  run both lenses and tabulate the route: convergent-phenotype × {same-gene,
  different-gene}. This is the scientifically deepest and most publishable thread.
- **Honesty caveats:** carries the Track-C circularity lesson — molecular-site
  selection must be independent of the phenotype labels.
- **Exit criteria:** a routing table for ≥5 convergences, each lens independently
  verified.

---

## Cross-cutting infrastructure (several directions need these)

- **A data layer** (Phase 3 of `EXPANSION_PLAN.md`): the Atlas (#3), pathogen (#2),
  and any scale needs taxa-as-objects + typed-morphism storage — which is exactly
  where KOMPOSOS's `Category` belongs, tying this back into the main project.
- **A pluggable "domain adapter"**: formalize the (tree, traits, null) contract so
  languages (#4), models (#1), and pathogens (#2) drop in without touching the core
  gauntlet. The honesty contract becomes a shared interface, not per-domain code.
- **Calibration harness**: an adversarial "manufacture fake convergence and confirm
  the gauntlet rejects it" suite — turns the honesty contract into a measured
  false-positive rate, and is reusable across every domain.

---

## Recommended sequencing

1. **Inevitability Index (#6)** first — it needs almost no new data (it is a new
   statistic over gauntlet output you already produce), and it reframes the whole
   project from "a pile of T2 verdicts" into "a quantitative answer to replay-the-
   tape." It also produces the Atlas (#3) as a by-product and the biosignature prior
   (#7) for free.
2. **Neural-net convergence (#1)** in parallel — most novel, most current, and
   `esmc_diagnose.py` ports almost 1:1; it proves the engine is genuinely general.
3. **Genotype↔phenotype routing (#9)** as the deep-science spine — it is the most
   publishable and the molecular lens already exists.
4. Everything else (pathogen, linguistic, design, topology, technosignatures) becomes
   a domain adapter once the cross-cutting infra and the Atlas exist.

**The throughline:** this was never a biology tool. It is a general epistemics of
*independent rediscovery* — and that is the thing worth being grandiose about.
