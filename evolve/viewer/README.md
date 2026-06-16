# Mammal Tree Explorer — slices 1 + 2

A zoomable, click-to-drill Order → Family → Genus → Species browser over the
**real dated mammal tips**, with a **convergence overlay** that runs the thesis
pipeline (Fitch origins, PTP, dated-Mk tier verdict) on the real chronogram.

## Build it

```bash
python evolve/viewer/export_mammal_tree.py
```

This reads the cached data the thesis already uses —

- `evolve/cache/mammal_mcc.nwk` (PHYLACINE / Upham 2019 dated chronogram),
- `evolve/cache/pantheria.txt` (PanTHERIA traits + taxonomy) —

joins them, and writes a single self-contained file:

```
evolve/viewer/mammal_tree.html
```

Open it in any browser. (Data is embedded; D3 loads from a CDN, so the first
open needs internet. Vendoring D3 for a fully-offline figure is a later slice.)

## What you can do

- **Click a node** to expand/collapse its children — drill from order down to
  species and back out.
- **Scroll / drag** to zoom and pan.
- **Color by** — taxonomic order, body mass (log g), mean temperature (°C),
  trophic level, or activity cycle. Internal nodes show the mean of their
  descendants, so you can read trait structure at any rank.
- **Search** a genus or species; matching paths auto-expand and highlight.
- **Detail panel** (right) shows traits for the selected node and fetches a
  **PhyloPic silhouette** for the taxon — a light taste of the "picture on every
  node" idea. It degrades gracefully when offline or when no silhouette exists.

## Convergence overlay (slice 2)

Pick a "⚡ Convergence overlay" entry in the **View** dropdown. For that trait the
viewer shows the result of the real thesis pipeline, run at build time on the
dated chronogram:

- **Tip colour = independent origin.** Each derived species is coloured by *which*
  independent gain it descends from (Fitch parsimony). The same trait arising in
  distant orders shows up as different-hued tips scattered across the taxonomy —
  the visual signature of convergence.
- **Internal nodes heat up** with the number of independent origins they contain,
  so convergence hotspots are visible even while collapsed.
- **Orange arcs = cross-lineage bridges** between frontier carriers of *different*
  origins (the convergent-pair signature; the analogue of the negative-Ricci
  bridges in the thesis). Expand a clade and the bridges refine.
- **Verdict badge** (top-left) stamps the tier: independent origins, PTP
  p(fewer), Mk observed/null and p(excess), and the **T0–T3** verdict.

Built-in overlays and their current verdicts (dated Upham 2019 tree). The first
group is the **landscape-defining body-plan / locomotion convergences**, coded by
MSW05 family/order (external taxonomy, independent of the trait tested — so
counting origins is a legitimate convergence test, not fitting). The second group
is trait-coded from PanTHERIA life-history columns.

| Convergence | Independent origins | PTP | Mk p(excess) | Tier |
|---|---:|---|---:|---|
| **Marine body plan** (whales, sirenians, seals) | 3 | 0.002 ✓ | 0.613 | **T2** |
| **Aerial** (bats + glider lineages) | 4 | 0.002 ✓ | 0.567 | **T2** |
| **Fossorial digging** (moles, mole-rats ×7) | 7 | 0.002 ✓ | 0.530 | **T2** |
| **Myrmecophagy** (anteaters, pangolins, aardvark, echidna, numbat) | 5 | 0.002 ✓ | 0.570 | **T2** |
| **Bipedal hopping** (kangaroos, jerboas, springhares) | 3 | 0.002 ✓ | 0.613 | **T2** |
| Carnivory (trophic = carnivore) | 74 | 0.002 ✓ | 0.560 | **T2** |
| Diurnality (active by day) | 74 | 0.002 ✓ | 0.507 | **T2** |
| Giant body (top 10% mass) | 31 | 0.002 ✓ | 0.523 | **T2** |

The clade-coded origin counts recover the textbook numbers (marine = 3 returns to
water; myrmecophagy = 5 anteating lineages), which validates the Fitch
reconstruction on the real dated tree. All overlays land at **Tier 2** —
structured homoplasy, real and non-random, but not excess beyond neutral drift on
dated branches. This mirrors the fish results: the overlay *downgrades* claims
rather than overstating them.

Family lists for the clade-coded traits live at the top of
`export_mammal_tree.py` (`MARINE_FAMS`, `GLIDER_FAMS`, `FOSSORIAL_FAMS`,
`MYRMECOPHAGY_FAMS`, `HOPPING_FAMS`) — edit them to refine a syndrome. The
"aerial" set is family-level and therefore approximate (not every member of a
glider family glides). Monte-Carlo budgets are modest for build speed
(`N_PERM=499`, `N_SIM=299`); raise them in `export_mammal_tree.py` for a final
figure.

## Provenance (honest)

- Every leaf is a species that is **both** on the dated chronogram **and** has a
  PanTHERIA row, so all tips are real dated taxa with real trait data. Current
  coverage: 3,679 species · 1,112 genera · 152 families · 29 orders (joined from
  4,253 chronogram tips and 5,416 PanTHERIA rows).
- The navigation hierarchy is **taxonomic** (PanTHERIA Order/Family/Genus), which
  is what "click by rank" means. The chronogram defines the tip set; rendering
  the real **branch-length backbone** (time on the x-axis) is the next slice.
- Trait values are real PanTHERIA measurements; missing values (`-999`) are shown
  as "no data" and colored grey, never imputed.

## Verifying every viewer (`verify_all.py`) — the standing CI guard

```bash
python evolve/viewer/verify_all.py            # all built viewers (exit 0 = pass)
python evolve/viewer/verify_all.py fish bird  # a subset
python evolve/viewer/verify_all.py --strict   # missing HTML is a hard failure (CI mode)
```

This is the structural fix for the "looked done, was broken" failure mode (a
viewer that opens without crashing can still render the wrong data). It loads each
self-contained HTML in headless Chromium and asserts: **zero page/console errors**
throughout; a sane embedded species count; for every convergence overlay a verdict
badge + tier, a legend, and cross-lineage bridges **when ≥2 independent origins
exist**; and origin counts above biological sanity floors. It self-discovers each
viewer's real base colorings from the dropdown, so it never tests a control the
viewer doesn't offer. The viewers' HTML is committed, so
[`.github/workflows/verify-viewers.yml`](../../.github/workflows/verify-viewers.yml)
runs this on every push/PR with only a headless browser — no data re-download.
(`verify_bird_tree.py` is the original single-viewer pattern this generalizes.)

## Birds (Phase 1) — same pattern, new clade

The bird viewer is the first proof that the exporter generalizes beyond mammals.

```bash
python evolve/viewer/export_bird_tree.py     # writes evolve/viewer/bird_tree.html
python evolve/viewer/verify_bird_tree.py      # headless Playwright check (exit 0 = pass)
```

It joins the **real dated bird tree** (Jetz et al. 2012, BirdTree / VertLife
Stage2 Hackett backbone, 9993 OTUs) with **AVONET** (Tobias et al. 2022) and runs
the identical Fitch → PTP → dated-Mk → tier pipeline. Data is fetched + cached on
first build:

- `evolve/cache/bird_tree.nwk` — one pseudo-posterior tree, streamed out of the
  190 MB VertLife archive (we decompress only the first tree, ~0.5 MB; node ages
  are clock+fossil estimates independent of traits, so they cannot fit a result);
- `evolve/cache/avonet.xlsx` — AVONET BirdTree sheet (9993 species, 1:1 join).

Most bird overlays are **measured** straight off AVONET's ecology columns
(`Trophic.Niche`, `Primary.Lifestyle`), not clade-coded — even stronger than the
mammal overlays. Flightlessness is **clade-coded** to the fully-flightless
families (ratites + penguins; tinamous fly and are excluded, so parsimony recovers
flight loss as independent gains). Current verdicts (one posterior Jetz tree):

| Convergence | Coding | Independent origins | Mk p(excess) | Tier |
|---|---|---:|---:|---|
| **Nectarivory** (hummingbirds, sunbirds, honeyeaters, lorikeets…) | measured | 48 | 0.510 | **T2** |
| **Aerial lifestyle** (swifts, swallows, nightjars…) | measured | 71 | 0.557 | **T2** |
| **Aquatic lifestyle** (penguins, auks, grebes, loons, ducks…) | measured | 35 | 0.573 | **T2** |
| **Vertivory / raptorial** (hawks, falcons, owls, shrikes…) | measured | 34 | 0.553 | **T2** |
| **Flightlessness** (ratites + penguins) | clade-coded | 2 | 0.797 | **T2** |
| Frugivory | measured | 197 | 0.003 ✓ | **T3** |
| Granivory | measured | 141 | 0.003 ✓ | **T3** |
| Giant body (top 10% mass) | measured | 140 | 0.003 ✓ | **T3** |

The measured-ecology convergences land at **T2** (structured, real, but not beyond
neutral drift on dated branches) — the same honest downgrade seen in mammals and
fish. Flightlessness is **conservatively** 2 origins: parsimony merges the
palaeognath ratite losses into one and counts penguins separately (molecular
evidence argues for more independent ratite losses; the viewer reports what the
parsimony reconstruction on this tree actually says, not the wished-for number).
Base **View** colorings: order, body mass, Hand-Wing index (dispersal), trophic
level, primary lifestyle.

## Fish (Phase 2/A) — Phenoscape anatomy, honest polarity

```bash
python evolve/viewer/export_fish_tree.py     # writes fish_tree.html
```

Real dated ray-finned fish tree (Rabosky et al. 2018, 11,638 species). Order/Family
resolved via Open Tree of Life (cached; 97% of genera resolve to 67 real orders).
Traits are **Phenoscape** presence/absence annotations on anatomy terms (fins,
barbels). Two honesty points specific to fish:

- **Polarity is reconstructed, not assumed.** Fitch decides the ancestral state and
  we count origins of the *derived* state — which is a **loss** for fins but a
  **gain** for barbels and the adipose fin. (A naive "absent = derived" coding would
  mislabel those gains.) Current verdicts: pelvic-fin loss **4 origins (T2)**;
  adipose-fin **gain 6 origins (T2)**; barbel **gain 6 origins (T2)**; pectoral /
  dorsal / caudal fin **single origin (T0)** — downgraded, not inflated.
- **Dollo guardrail.** For loss characters, a taxon coding "absent" *outside* the
  clade that ever had the structure never had it (primitive absence, not a loss), so
  `dollo_recode` drops it rather than counting a phantom convergent loss.
- Stats run on the full 11.6k-tip tree; the browser renders a 3,000-species sample
  that always keeps every derived-state carrier, so no origin is hidden.

## Unified vertebrate backbone (Phase 2)

```bash
python evolve/viewer/export_vertebrate_tree.py   # after the 3 clade viewers exist
```

Splices the fish, mammal, and bird hierarchies under one dated backbone
(`Vertebrata → Gnathostomata → {Actinopterygii, Sarcopterygii→Tetrapoda→Amniota→
{Mammalia, Sauropsida→Aves}}`) — **16,672 species across 3 classes** in one
navigable view, coloured by Major clade or taxonomic order. The interactive
overlays stay clade-local for now; this is the navigable skeleton that the
cross-class convergence work (below) builds on.

## Cross-class aerial convergence (Track A) — `crossclass_aerial.py`

```bash
python evolve/crossclass_aerial.py     # -> results/crossclass_aerial.json
```

The vertebrate splice above is **hierarchy-only** — Fitch origin-counting works on
topology, but PTP and dated-Mk need real branch lengths *spanning* the classes. So
this script does the real work: it **grafts** the three dated chronograms onto a
deep-node-calibrated backbone and runs the identical gauntlet on the dated
cross-class tree.

- **Backbone (ultrametric, Myr):** Actinopterygii | Sarcopterygii split at **~430
  Ma** (TimeTree 416–439) and Amniota (Mammalia | Sauropsida) at **~319 Ma** (fossil
  calibration). Each clade is attached by a stem = `calibration − its crown depth`
  (fish crown ~368, mammals ~218, birds ~101 Ma). Both deep ages are clock+fossil
  estimates **independent of the aerial trait**, so they cannot fit the result.
- **Trait:** aerial = powered flight **or** gliding, clade-coded by external
  taxonomy — birds (all but flightless ratites + penguins), bats (Chiroptera),
  flying fish (Exocoetidae), mammalian gliders (Petauridae, Anomaluridae, …).
- **Result (balanced 1,908-tip graft):** **6 independent origins** of aerial
  locomotion (3 for the named groups only); **PTP p = 0.001** (z ≈ 57 — strongly
  clustered, real structured convergence); **Mk p_excess = 1.0** (z ≈ −11 — *far
  fewer* transitions than neutral drift, i.e. flight is deeply conserved once gained,
  **not** excess homoplasy) → **Tier 2**. The same honest downgrade every body-plan
  overlay shows: structured and real, not inflated to T3.
- **Sensitivity** (built in): 6 origins / T2 stable across three non-carrier sample
  seeds; 3 origins / T2 with mammalian gliders excluded. Effect sizes (z vs both
  nulls) are reported in the JSON.

**Honest caveats:** the graft rests on two point calibrations (reported, trait-
independent); non-carrier taxa are seeded-subsampled for tractability while **every
aerial carrier is kept**, so no origin is hidden; flightless birds are coded
non-aerial (a reversal the reconstruction handles, not a phantom origin). This is
the verified *analysis*; rendering it as an interactive overlay inside
`vertebrate_tree.html` is the next presentation step.

## Molecular lens (Phase M) — ESM-C Prestin pilot

`evolve/esmc_prestin_pilot.py` embeds real Prestin (SLC26A5) sequences with ESM-C,
subtracts ancestry via the dated mammal tree, and permutation-tests whether
echolocating bats + toothed whales are closer in embedding space than phylogeny
predicts. At 14 taxa with the 300M model it is an **honest negative**
(`p = 0.13`, controls also read slightly convergent, verdict = false): mean-pooled
embeddings capture Prestin's broad conservation but not the convergent
substitutions at this scale. Logged transparently in `results/prestin_esmc.json`,
not forced to a fit.

### Diagnosis — *why* the pilot washed out (`esmc_diagnose.py`)

```bash
python evolve/esmc_diagnose.py        # -> results/prestin_esmc_diagnosis.json
```

Before scaling to ESM-C 600M/6B or more proteins, this script runs three
self-contained diagnostics on the **same 14 cached sequences/embeddings**:

1. **Geometry.** PCA of the mean-pooled vectors: echolocators only *weakly*
   separate (echo-vs-rest cosine silhouette **+0.13**); PC1 partly tracks sequence
   length (r = −0.45), a confound. So mean-pooling barely sees the convergence.
2. **Does ESM beat raw % identity?** Running the *identical* ancestry-subtracted
   test on a %-identity distance gives **p = 0.9995** (echo pairs are *less*
   similar than phylogeny predicts), versus ESM mean-pool **p = 0.13**. So ESM-C
   *does* add function-aware information over the trivial baseline — just not
   enough, mean-pooled, to clear significance.
3. **Localize.** Data-derived parallel-substitution sites (columns where echo bats
   *and* toothed whales share a residue that differs from each clade's immediate
   non-echo relative — the Liu/Li 2010 signature, derived here not asserted) number
   just **2 of 741 residues** (ref positions 576=K, 641=V). A per-residue ESM test
   restricted to those sites is strongly significant (echo resid −0.019,
   **p = 0.0009**, controls pass).

**Honest verdict.** The pilot's negative was a **measurement artifact of
mean-pooling**, not an absence of molecular convergence: the signal lives in
~0.3 % of the protein and averaging over 741 residues dilutes it ~370×. **Circularity
caveat (stated in the JSON):** those 2 sites were *selected using the echo labels*,
so the per-site test is partly circular — a per-site **%-identity** baseline at the
same columns is *also* significant (p = 0.003), showing the per-site significance is
mostly site-selection, not independent ESM discovery. The robust, non-circular
takeaways are (1) ESM beats %-identity mean-pooled and (2) the localization/dilution
explanation. **The next honest step is non-circular site selection** (literature
ASR or a held-out tree) feeding per-residue ESM at 600M — *then* the molecular tier
is earned, not assumed.

## Not yet (later slices)

- Real branch-length / time backbone instead of a taxonomic ladder.
- Trait + ancestral-state **time scrubber** (the honest version of the morph
  dream).
- Vendored D3 for a fully offline supplementary figure.
- **Beyond mammals → the whole tree of life** (birds next, then vertebrates,
  plants, and an environment layer): see [`EXPANSION_PLAN.md`](EXPANSION_PLAN.md).
