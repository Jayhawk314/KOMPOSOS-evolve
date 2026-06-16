# Tree-of-Life Expansion Plan

How the mammal viewer + convergence overlay grows toward the full animal kingdom,
eventually all life, and an integrated non-living environment layer — without
trying to do it all at once.

> Status: Phase 0 complete (mammals); Phase 1 complete (birds — Jetz tree + AVONET,
> `export_bird_tree.py`); Phase 2 in progress — fish viewer (`export_fish_tree.py`,
> Rabosky tree + Phenoscape) and a unified Vertebrata backbone
> (`export_vertebrate_tree.py`, splices fish+mammals+birds, 16.7k species) both
> built and headless-verified. Molecular pilot (Phase M): ESM-C Prestin run is an
> honest negative at 14 taxa, now **diagnosed** (`esmc_diagnose.py`): the negative
> is a mean-pooling artifact — the convergence signal is localized to ~2/741
> residues, so averaging dilutes it ~370× (see below). This is the roadmap, not a
> promise of dates.

---

## North star

An interactive, zoomable tree where you can drop into **any clade**, color it by a
trait, and ask the thesis question — *did this trait converge independently, and
how far does the claim survive?* — with the same tiered verdict (Fitch origins →
PTP → dated-Mk → tier) the dissertation defends. Eventually the same view can
overlay the **non-living environment** as the driver, turning "independent
origins" into "independent origins that share a cause."

## Guiding principles

1. **Convergence is clade-local.** You never need all of life loaded at once; you
   need a navigable backbone plus the ability to analyze whatever clade is open.
   This is what makes the moonshot decomposable.
2. **Traits are the bottleneck, not trees.** Trees exist at scale; comparable
   trait codings do not. Coverage is rich for vertebrates + plants and thins to
   near-empty for most invertebrates, fungi, and microbes. Plan around that.
3. **Honesty over coverage.** Missing data is greyed, never imputed. Synthetic /
   clade-coded traits are labeled as such. Verdicts are downgraded, not inflated.
4. **Each phase ships something useful on its own.** No phase is wasted scaffolding
   for a later one.

---

## The three layers

### Layer A — the tree (scales well)
- **Backbone:** Open Tree of Life synthetic tree (~2.3M named-species tips; API +
  bulk download). Topology only; soft spots in places.
- **Dated subtrees, grafted per clade where they exist:**
  | Clade | Chronogram | Status |
  |---|---|---|
  | Ray-finned fish | Rabosky et al. 2018 (Fish Tree of Life, 11.6k) | in repo |
  | Mammals | Upham et al. 2019 (PHYLACINE, ~4–5k) | in repo |
  | Birds | Jetz et al. 2012 (birdtree.org, ~9993) | in repo (streamed + cached) |
  | Squamates | Tonini et al. 2016 | future |
  | Amphibians | Jetz & Pyron 2018 | future |
  | Seed plants | Smith & Brown 2018 (~356k) | future |
  | Bacteria/Archaea | GTDB | far future (see reticulation caveat) |
- **Pattern:** Open Tree topology as the global skeleton; splice dated clade trees
  for the analysis (Mk needs real branch lengths).

### Layer B — traits (the bottleneck)
| Clade | Trait sources |
|---|---|
| Mammals | PanTHERIA (in repo), COMBINE, EltonTraits |
| Birds | **AVONET** (morphometrics — excellent), EltonTraits |
| Fish | FishBase, Phenoscape (in repo) |
| Amphibians | AmphiBIO |
| Plants | **TRY** (huge), BIEN |
| Inverts / fungi / microbes | sparse — expect data deserts |

Body-plan/locomotion convergences (marine, fossorial, aerial, anteating, hopping)
are codable by **taxonomic family/order** even when no trait column exists —
external taxonomy is independent of the trait tested, so counting origins stays a
legitimate test, not fitting. This is the technique already used in the mammal
overlay and it generalizes to every clade.

### Layer C — the non-living environment (different shape)
Not a tree — **geospatial/temporal layers**, linked to organisms via occurrences.
- **Climate:** WorldClim, CHELSA; paleoclimate via PaleoClim.
- **Biomes / ecoregions:** Olson et al.
- **Bridge to organisms:** GBIF occurrences → sample environment at each point →
  per-species environmental profile.
- **Payoff:** the thesis's environment-as-driver lens at scale. Pick a trait, see
  whether its independent origins coincide with the **same** environmental shift —
  convergence *with a shared cause*. This is the most novel and most
  thesis-aligned endpoint.

### Layer D — molecular convergence (ESM-C protein embeddings)
A new **detection lens at the molecular level**, not a tree source and not a fix
for the body-plan overlays. Convergence often happens in proteins: independent
lineages reaching the *same* molecular solution (Prestin → echolocation in bats
and toothed whales; opsin spectral tuning; antifreeze glycoproteins; foregut
lysozymes; toxin-resistant ion channels; high-altitude hemoglobins).

- **Tool:** ESM-C (EvolutionaryScale, the lineage after Meta's ESM-2/ESMFold).
  Input = amino-acid sequence; output = function-aware embeddings. 300M/600M run
  locally on a GPU; 6B via the Forge API.
- **What it buys:** a similarity that reflects structure/function, not raw sequence
  identity — so two distant lineages whose protein lands in the same functional
  neighborhood (despite low identity / distant ancestry) read as a candidate
  *molecular* convergence. This is the molecular analogue of the morphology lens.
- **Data:** per-taxon orthologous protein sequences — UniProt, NCBI RefSeq,
  **OrthoDB**/Ensembl. Rich for genome-sequenced organisms (a few hundred mammals,
  model species); a desert for most taxa. Same rich-for-some constraint as traits.
- **How it slots in:** a species' set of protein embeddings is a molecular
  relational fingerprint (the Yoneda lens, at the molecular level). It reuses the
  **same gauntlet** — subtract ancestry, then Fitch origins → PTP → dated-Mk →
  tier. ESM does **not** replace the nulls; it is a new lens feeding the same
  staged verdict, and it adds a **molecular tier**: does a phenotypic convergence
  have an independent molecular signature? (Prestin/echolocation says yes — a
  convergence-of-convergence corroboration.)

**Honesty caveats:**
- Embedding similarity ≠ convergence — similar embeddings can also be homology, so
  the ancestry correction is mandatory, exactly as for morphology.
- Do **not** use ESM distances for tree-building (alignment-free embedding
  phylogeny is confounded). Keep it as a trait/function lens only.

**Pilot (entry point): the Prestin experiment.** Pull Prestin (SLC26A5) across the
mammals already in the viewer, embed with ESM-C, and test whether bats + toothed
whales cluster in embedding space against a phylogenetic null — i.e. whether the
molecular lens *independently recovers* the echolocation convergence that
`echolocation_flagship.py` already shows phenotypically. Tiny data, high payoff;
proves the lens before committing to large-scale sequence acquisition.

---

## Architecture evolution

### Current (good ≤ ~10k tips)
Python exporter → one self-contained HTML with all data embedded as JSON.
Everything precomputed at build time. Perfect for a single clade; cannot embed or
render millions of nodes.

### Scaled (needed beyond ~tens of thousands of tips)
1. **Data layer** instead of an embedded blob — a store of taxonomy + dated
   subtrees + trait tables + environment/occurrences. Start SQLite/DuckDB, grow to
   Postgres. **This is where KOMPOSOS's `Category` fits:** taxa as objects;
   taxonomy, traits, and environment as typed morphisms — the tree-of-life store
   *is* a Category, tying the viewer back into the main project instead of being a
   side tool.
2. **On-demand compute** — the Fitch → PTP → Mk pipeline runs per requested
   `(clade, trait)` and caches the verdict. Never precompute everything.
3. **Streaming renderer** — the client requests only the subtree under the node
   being expanded; Canvas/WebGL for large fans. It holds only what is on screen.

---

## Honesty cautions (sharper at scale)

- **Reticulation breaks the tree model for microbes.** Horizontal gene transfer
  makes "convergence vs. transfer" indistinguishable on a tree. The reticulate
  caveat (already an executable chapter) becomes central: refuse tree-based
  convergence claims for bacteria/archaea, or switch to a phylogenetic network.
- **Trait data deserts** must be shown, not hidden — large branches will have no
  scorable traits.
- **Compute budgets** — Mk/PTP on huge trees need larger Monte-Carlo budgets to
  keep verdicts honest; do not ship under-powered tiers.
- **Dating gaps** — many clades lack chronograms; analyses there fall back to
  topology-only nulls (PTP), with no dated-Mk tier (report T2-max, not T3).

---

## Staged roadmap

| Phase | Deliverable | Data | Effort | Risk | Exit criteria |
|---|---|---|---|---|---|
| **0 ✅** | Mammal viewer + 8 convergence overlays | Upham + PanTHERIA (in repo) | done | — | done |
| **1 ✅** | **Birds** viewer, same pattern | Jetz tree + AVONET | done | — | done — nectarivory (48 origins), aerial (71), aquatic (35), vertivory (34), flightlessness (ratites+penguins) all recovered |
| **2 ◑** | Vertebrate multi-class view (fish + mammals + birds folded in; amphibians/squamates next) | Open Tree vertebrate subtree + dated clade trees | medium | medium | ✅ one view spanning 3 classes (`vertebrate_tree.html`, 16.7k species); amphibians/squamates remain |
| **A ◑** | **Cross-class convergence on a dated backbone** (`crossclass_aerial.py`): graft fish+mammal+bird chronograms via deep-node calibrations, run the full gauntlet | dated clade trees + TimeTree/fossil ages | medium | medium | ✅ aerial locomotion = 6 origins, PTP p=0.001, Mk p_excess=1.0 → **T2** (stable across sensitivity sweep); interactive viewer overlay remains |
| **3** | **Architecture shift**: data layer + on-demand compute + streaming renderer | (re-platform) | high | medium | render & analyze a >50k-tip clade smoothly |
| **4** | Plants | Smith & Brown tree + TRY | high | medium | a plant convergence (e.g. C4, succulence) overlay |
| **5** | **Environment layer** + occurrences → convergence-with-driver | WorldClim/CHELSA + GBIF | high | medium | a trait's origins testable against a shared environmental shift |
| **6** | All-life backbone | Open Tree (2.3M) + GTDB | very high | high | navigable to any clade; microbe branches flagged for reticulation |
| **M (parallel) ◑** | **Molecular convergence lens (ESM-C)** — Prestin pilot (negative) + diagnosis done; next is non-circular site selection at 600M (Layer D) | OrthoDB/UniProt sequences + ESM-C | low (pilot) | low | pilot negative at 14 taxa; diagnosed as a mean-pool dilution artifact (signal in ~2/741 residues); independent recovery still pending |

**Recommended immediate step:** Phase 1 (birds). It reuses the exact exporter
pattern, has superb trait data (AVONET), and famous convergences — so it cheaply
proves the whole approach generalizes beyond mammals *before* committing to the
Phase 3 re-platform.

**Parallel low-cost option:** the Phase M **Prestin pilot** is independent of the
tree-scaling track — it needs only the mammals already in the viewer plus one
protein's sequences, so it can run anytime as a self-contained molecular-tier
proof and extends `echolocation_flagship.py` directly.

## What stays the same across all phases
- The verdict pipeline: Fitch origins → PTP randomization → neutral-Mk on dated
  branches → tier (T0–T3).
- The honesty contract: real tips only, missing data greyed, clade-coded traits
  labeled, claims downgraded not inflated.
- The interaction: zoom, click-to-drill by rank, color-by-trait, convergence
  overlay with origin colors + cross-lineage bridges + verdict badge.
