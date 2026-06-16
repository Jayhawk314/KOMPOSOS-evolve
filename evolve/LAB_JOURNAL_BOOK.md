# KOMPOSOS-Evolve: Exploratory Lab Journal
**An overview of the categorical convergence engine — reconciled against the implementation**

> Note: an earlier draft of this journal was generated automatically and overstated
> a few results. It has been corrected against what the code actually does (the
> echolocation tier, the fish overlay attribution, and the scope claims below). Where
> the engine downgrades or fails to find a signal, this journal now says so.

---

## Preface: The Honest Machine

This journal documents the theoretical foundations, implementation details, and empirical case studies of the `KOMPOSOS-evolve` engine. The core philosophy of this system is **honesty**. It does not force signals where there are none. It does not hand-curate traits to produce a desired clustering. It treats evolution as a geometric and categorical structure, and measures convergence through strict topological and spectral bounds.

When a signal is noisy (as in the ESM-C molecular pilot, Chapter 7), the engine reports a failure to reject the null hypothesis. When a trait coded "absent" is primitive rather than a true loss (a species that never had the structure), Dollo recoding drops it transparently rather than counting a phantom convergence.

This is a working journal over a system still under construction, not a finished treatise. Current scope: the Fitch → PTP → dated-Mk tier pipeline and interactive viewers exist for mammals, birds, and ray-finned fish (spliced into one vertebrate backbone); the molecular lens is a single-protein pilot that returned a negative; amphibians, squamates, plants, and the environment layer are roadmap, not done.

---

## Chapter 1: The Philosophy of Categorical Evolution

Traditionally, evolutionary convergence is treated as a statistical anomaly on a phylogenetic tree. `KOMPOSOS-evolve` reimagines this via Category Theory. 

In our system, every entity is an **Object** in a Category:
- A species (e.g., *Tursiops truncatus*)
- A morphological state (e.g., `pelvic fin: absent`)
- A behavioral niche (e.g., `dietbreadth=q2`)
- An ancestral lineage (e.g., `Cetacea`)

The relationships between them are **Morphisms**, grouped into independent sub-categories called **Lenses**:
- `morphology`: Taxon $\xrightarrow{has}$ Trait
- `environment`: Taxon $\xrightarrow{inhabits}$ Niche
- `ancestry`: Taxon $\xrightarrow{member\_of}$ Clade

### The Yoneda Fingerprint

The engine measures similarity not by distance on a tree, but by a taxon's **Yoneda fingerprint**—its entire profile of relationships within a specific lens. 

If we want to know how convergent a dolphin and a bat are morphologically, we calculate their Yoneda similarity in the `morphology` lens, and then *subtract* their Yoneda similarity in the `ancestry` lens.

```python
# The Core Equation
convergence(A,B | lens) = yoneda_sim(A,B | lens) - yoneda_sim(A,B | ancestry)
```
If this value is high, it means they share a massive amount of trait architecture but share very little evolutionary history. That is the categorical signature of convergence.

---

## Chapter 2: The Core Engine (`convergence_engine.py`)

Let's explore how to build and query the `ConvergenceModel` interactively.

### Example: Building a Toy Model
```python
from convergence_engine import ConvergenceModel

m = ConvergenceModel()
taxa = ["Bat", "Dolphin", "Cow", "Mouse"]
for t in taxa: m.add_taxon(t)

# Adding ancestry (Lineage)
m.relate("ancestry", "Bat", "Laurasiatheria", rel="member")
m.relate("ancestry", "Dolphin", "Laurasiatheria", rel="member")
m.relate("ancestry", "Cow", "Laurasiatheria", rel="member")
m.relate("ancestry", "Mouse", "Euarchontoglires", rel="member")

# Adding specific inner clades
m.relate("ancestry", "Bat", "Chiroptera", rel="member")
m.relate("ancestry", "Dolphin", "Cetacea", rel="member")
m.relate("ancestry", "Cow", "Artiodactyla", rel="member")

# Adding morphology
m.relate("morphology", "Bat", "echolocation", rel="has")
m.relate("morphology", "Dolphin", "echolocation", rel="has")
m.relate("morphology", "Cow", "hooves", rel="has")

# Querying the Engine
print("Ancestry Similarity (Bat, Dolphin):", m.similarity("ancestry", "Bat", "Dolphin"))
print("Morphology Similarity (Bat, Dolphin):", m.similarity("morphology", "Bat", "Dolphin"))
print("Convergence Score:", m.convergence("morphology", "Bat", "Dolphin"))
```

Because the Bat and Dolphin share `echolocation` but diverge at the order level, their morphology similarity dwarfs their ancestry similarity, yielding a strong positive convergence score.

---

## Chapter 3: The Geometry of Traits (Ricci & Spectral)

The relational Yoneda similarity is triangulated with two advanced geometric detectors imported from the `komposos_wesys` geometry module.

### Ollivier-Ricci Curvature
We construct a **Lineage-Community Graph**. If two taxa share high ancestry, they are connected by an `ancestry` tie. If they share high trait similarity, they are connected by a `bridge` tie.

Using optimal transport, we calculate the Ollivier-Ricci curvature ($\kappa$) of every edge. 
- Positive curvature: The edge is deeply embedded inside a local community (a single clade).
- **Negative curvature:** The edge spans a structural void between two distinct communities.

Convergent pairs appear as the most *negatively curved* bridges in the space. They are the evolutionary wormholes connecting distant branches of the tree of life.

### Spectral Embedding
We compute the graph Laplacian of the trait-similarity graph and extract its eigenvectors. This embeds the taxa in a low-dimensional morphospace. Taxa that cluster together spectrally but lie far apart phylogenetically are undergoing parallel attractors.

```python
# Extracting the most negative Ricci bridges
bridges = m.ricci_bridges(trait_lenses=["morphology"])
for b in bridges[:5]:
    print(f"{b['pair']}: Curvature {b['kappa']:.3f}, Is Bridge: {b['is_bridge']}")
```

---

## Chapter 4: Data Synthesis (PanTHERIA & Phenoscape)

Real evolution requires real data. The engine relies on automated, reproducible ingestion scripts rather than hand-curated spreadsheets.

### Ingestion Source 1: Open Tree of Life
Using `ingest.py`, the system takes a list of species strings, pings the Open Tree TNRS (Taxonomic Name Resolution Service), retrieves their valid OTT IDs, and downloads their true ancestral lineage up to the root.

### Ingestion Source 2: PanTHERIA (Continuous Traits)
PanTHERIA provides raw, continuous life-history traits (body mass, habitat breadth, absolute latitude). The engine automatically **quantile-bins** these across the sample.
- If a bat and a dolphin both fall into `bodymass=q3` (the highest quantile of the sample), they are categorically linked to the same object.

### Ingestion Source 3: Phenoscape (Ontology-Annotated Anatomy)
Using `phenoscape_ingest.py`, the system maps common names like "pelvic fin" to strict Uberon IRIs (e.g., `UBERON:0000151`). It queries the live Phenoscape KB for all taxa annotated with states for this entity.

### Ingestion Source 4: AVONET (Bird Morphometrics + Ecology)
`avonet()` in `ingest.py` pulls the AVONET (Tobias et al. 2022) BirdTree-taxonomy sheet — measured morphometrics plus ecology (`Trophic.Niche`, `Primary.Lifestyle`) for all 9,993 species, joining 1:1 to the Jetz et al. 2012 bird chronogram. The bird viewer's nectarivory, aquatic, aerial, and vertivory overlays are read straight off these measured columns; only flightlessness is clade-coded.

---

## Chapter 5: Example I - Echolocation in Mammals

The flagship example (`echolocation_flagship.py`) tests the classic convergence of laryngeal echolocation between Microbats (e.g., *Myotis*) and Toothed Whales (e.g., *Tursiops*). It is a demonstration of the **Yoneda deep-convergence** detector (Chapters 1–3), not the Fitch/PTP/Mk tier pipeline.

### What the flagship actually computes
For each pair it measures trait similarity across the `behavior`, `morphology`, and `molecular` lenses, **subtracts** ancestry similarity, and ranks the result (`deep_convergence` / `rank_deep`).

- **Headline:** bat ~ toothed-whale pairs rank at the very top — high shared trait architecture, low shared ancestry: the categorical signature of convergence.
- **Negative control (the honesty teeth):** *Pteropus*, an Old-World fruit bat that does **not** laryngeally echolocate, is a close relative of the echo bats but must **not** rank as convergent with dolphins. It does not — the engine is detecting the convergent *solution*, not mere "battiness."
- **Caveat (important):** in this flagship the `molecular` lens is **hand-asserted** (a literal `prestin_echolocator_variant` label), not measured. Replacing that asserted label with a real, measured molecular signal is exactly the job of the ESM-C pilot in Chapter 7 — which is why that pilot matters even though it returns a negative.

### What the tier pipeline says about body-plan convergences
Separately, the Fitch → PTP → dated-Mk → tier pipeline (used in the interactive viewers) is run on body-plan / locomotion traits. The honest outcome there is a **downgrade**: the mammalian "aerial" convergence (bats + gliding lineages) recovers ≥2 independent origins with a significant PTP test, but does **not** exceed neutral Mk drift on the real dated branches — so it lands at **Tier 2 (structured homoplasy)**, not Tier 3. Across the mammal and bird viewers, essentially all body-plan convergences land at T2; only some labile ecological traits (frugivory, granivory) reach T3. The pipeline downgrades rather than inflates.

---

## Chapter 6: Example II - Fin Loss in Actinopterygii

Fishes (Actinopterygii) show massive phenotypic diversity. A common convergent event is the loss of specific fins (e.g., pelvic fins in eels and some catfishes).

### Polarity is reconstructed, not assumed
`export_fish_tree.py` does **not** hardcode "absent = derived." For each character it runs Fitch to reconstruct the ancestral state, then counts independent origins of whichever state is *derived*:

- Fins are ancestral to all sampled bony fishes, so a missing fin is a true **loss** (derived = `absent`).
- Adipose fins and barbels evolved *within* the fish tree, so for those the derived state is a **gain** (derived = `present`). A basal fish lacking a barbel hasn't "lost" it.

### The Dollo guardrail
For the **loss** characters, the viewer calls `dollo_recode` (from `polarity.py`): a taxon coding `absent` *outside* the clade descended from the structure's single origin never had it (primitive absence), so it is dropped rather than counted as a convergent loss. (A separate `phenoscape_inapplicability.py` explores the same idea on the Phenoscape side.)

**Actual output of `export_fish_tree.py`** (full 11.6k-tip Rabosky tree; PTP n=199, Mk n=99):
- `pelvic fin` — **loss**, 4 independent origins → **Tier 2** (structured homoplasy).
- `adipose fin` — **gain**, 6 independent origins → **Tier 2**.
- `barbel` — **gain**, 6 independent origins → **Tier 2**.
- `pectoral / dorsal / caudal fin` — single origin → **Tier 0** (no convergence; downgraded, not discarded).

---

## Chapter 7: Example III - The Molecular Lens (ESM-C Pilot)

In June 2026, we piloted the integration of Large Language Models to evaluate molecular convergence. We used **ESM-C (300M)**, a protein language model by EvolutionaryScale.

**The Hypothesis:** 
Can ESM-C embeddings of the Prestin protein (SLC26A5) — computed from amino-acid sequence alone — recover the phenotypic echolocation convergence *after* ancestry is subtracted? (The embedding uses no tree; the test explicitly removes phylogeny, regressing embedding distance on patristic distance from the dated mammal chronogram, so any residual closeness is the part not explained by shared descent.)

**The Execution (`esmc_prestin_pilot.py`):**
1. Fetch 14 full-length Prestin sequences from NCBI (Bats, Dolphins, Whales, Cows, Mice, Humans).
2. Pass sequences through ESM-C, extracting the mean representation of the final layer (960 dimensions).
3. Compute the pairwise distance matrix in ESM-C space.
4. Subtract the expected phylogenetic distance (derived from the mammal chronogram) to get a **residual**.
5. Compare the mean residual of Bat $\times$ Toothed-Whale pairs against a null distribution.

### The Honest Verdict
The pilot completed, and the results were an **honest negative**:
- Mean residual for echo-pairs: `-0.0002`
- Permutation $p = 0.1304$ (Not significant)
- **Negative Controls Failed:** The fruit bat (*Pteropus*, non-echolocating) and baleen whale (*Balaenoptera*, non-echolocating) also read as slightly convergent against the echo groups in embedding space.

**Conclusion:** At a sample size of 14, and using the global sequence embedding rather than targeted active-site residues, the 300M model captures general mammalian conservation but the specific convergent echolocation signal is drowned out by structural noise. The system documented this transparently rather than overfitting.

### The follow-up diagnosis (`esmc_diagnose.py`) — *why* it washed out

Rather than scaling straight to a larger model, we asked why the pilot failed, with three diagnostics on the *same* 14 cached sequences:

1. **Geometry.** A PCA of the mean-pooled vectors shows echolocators only *weakly* separate (echo-vs-rest cosine silhouette **+0.13**), and PC1 partly tracks sequence length (r = −0.45) — a confound, not function.
2. **ESM vs. a trivial baseline.** Running the identical ancestry-subtracted test on a raw %-identity distance gives **p = 0.9995** (echo pairs read *less* similar than phylogeny predicts), versus ESM mean-pool **p = 0.13**. So the function-aware embedding *does* add signal over % identity — just not enough when averaged.
3. **Localization.** Parallel-substitution sites — columns where echo bats *and* toothed whales share a residue differing from each clade's immediate non-echolocating relative (the Liu/Li 2010 signature, here **derived from the alignment, not asserted**) — number just **2 of 741 residues** (positions 576=K, 641=V). A per-residue ESM test restricted to those sites is strongly significant (echo residual −0.019, **p = 0.0009**, negative controls pass).

The honest reading: the pilot's negative was a **measurement artifact of mean-pooling**, not an absence of molecular convergence — the signal occupies ~0.3 % of the protein and averaging over 741 residues dilutes it ~370×. And the honesty teeth bite again: those 2 sites were *selected using the echo labels*, so the per-site test is partly **circular**; a per-site %-identity baseline at the same columns is *also* significant (p = 0.003), proving the per-site significance is mostly site-selection, not independent ESM discovery. The two claims that survive without circularity are (1) ESM beats % identity mean-pooled, and (2) the localization/dilution mechanism. A genuinely independent molecular tier requires selecting the convergent sites from an outside source (literature ASR or a held-out tree) and only then running per-residue ESM — ideally at 600M. The diagnosis is logged in `results/prestin_esmc_diagnosis.json`.

---

## Chapter 8: Future Horizons & WESyS Integration

The architecture of `KOMPOSOS-evolve` is inherently multi-domain. Because it abstractly maps Entities to categorical Lenses, it is not limited to biology.

### The Infinity Cosmos
By interfacing with the core `komposos` math libraries (ZFC, HoTT, Cubical), the convergence engine's definition of "Paths" mapping to shared phenotypes can be evaluated as **Homotopy Types**.
Two lineages arriving at the same state via disjoint interior paths form a **non-contractible loop** in the evolutionary space.

### WESyS (Waste-to-Energy Systems)
The same mathematics apply to power grid topologies and thermodynamic audits. 
- Taxa $\rightarrow$ Grid Nodes
- Traits $\rightarrow$ Energy States / Voltage
- Phylogeny $\rightarrow$ Transmission Hierarchy

Ricci curvature over a power grid identifies structural vulnerabilities (highly negative bridges that, if severed, cascade into failure).

---

## Final Thoughts

What works today: the engine synthesizes real APIs (Phenoscape, Open Tree of Life, NCBI) and datasets (PanTHERIA, AVONET, Rabosky and Upham and Jetz chronograms), runs the Fitch → PTP → dated-Mk → tier pipeline, and outputs interactive, headless-verified web viewers for mammals, birds, and fish, spliced into one unified vertebrate backbone (16,672 species across 3 classes). A standing CI harness (`verify_all.py`) re-verifies every viewer headlessly so "looked done, was broken" cannot recur silently.

Three general extensions now run beyond the per-clade viewers, each ending in an honest verdict rather than a forced fit:
- **Cross-class convergence on a dated backbone** (`crossclass_aerial.py`): grafting the fish/mammal/bird chronograms via deep-node fossil/clock calibrations (Osteichthyes ~430 Ma, Amniota ~319 Ma) and running the gauntlet on aerial locomotion gives **6 independent origins, strongly clustered (PTP p=0.001) but not in excess of neutral drift (Mk p_excess=1.0) → Tier 2** — the same downgrade the single-clade overlays show, stable across a sensitivity sweep.
- **The environment-as-driver lens, generalized** (`env_driver.py`): the whole mammalian convergence battery × every PanTHERIA climate axis, neutral-Mk-controlled, BH-FDR across the battery. The result is an **honest negative — 0 of 21 trait-environment couplings survive correction** — though the strongest raw signals (fossorial → dry/cool, diurnality → warm/wet) are biologically coherent.
- **The molecular lens, diagnosed** (`esmc_diagnose.py`): the Prestin negative is explained as a mean-pooling dilution artifact (signal in ~2/741 residues), with the per-site recovery flagged as partly circular.

What is honestly *not* done: amphibians, squamates, plants, and the full non-living environment layer for birds/fish (the GBIF → WorldClim occurrence bridge); interactive cross-class overlays inside the unified viewer; posterior-tree integration; and the molecular lens beyond the single-protein Prestin pilot with independent (non-circular) site selection. Chapter 8's WESyS / HoTT directions are vision, not implementation.

The point of the project is not breadth of claims but **rigorous verification** — including the verification that turned the first draft of this very journal from "fully operational" into the more accurate account above.

*(First draft generated automatically by Gemini CLI, 2026-06-14; reviewed and reconciled against the code, 2026-06-15.)*
