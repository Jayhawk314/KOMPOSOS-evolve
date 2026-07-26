# Convergent-Evolution Engine — Continuation Prompt

Paste the block below into a fresh Claude Code session. It is self-contained:
what the project is, what to read, what to reuse, the honesty contract, the next
tracks, and new ideas that are NOT yet in the roadmap. Scope is strictly the
convergent-evolution program (the WESyS / power-grid / HoTT material lives in a
different repo and is out of scope here).

---

You're continuing the **convergent-evolution engine** in
`C:\Users\JAMES\github\komposos-evolve` (branch `evolve-clean`, pushed to the
`evolve` remote = github.com/Jayhawk314/KOMPOSOS-evolve).

## READ FIRST (in order)
- `evolve/viewer/EXPANSION_PLAN.md` — the Phase 0–6 + M roadmap and the 3-layer
  (tree / traits / environment) + Layer D (molecular) design.
- `evolve/viewer/README.md` — what each viewer is and the honesty contract.
- `evolve/LAB_JOURNAL_BOOK.md` — the reconciled narrative + worked examples
  (treat it as accurate now; an earlier auto-generated draft was corrected).
- `evolve/esmc_prestin_pilot.py` + `evolve/results/prestin_esmc.json` — the
  molecular pilot and its honest-negative result.

## WHAT EXISTS (built + headless-verified this far)
- Viewers over REAL dated trees, each running the full gauntlet:
  `export_mammal_tree.py` (Upham + PanTHERIA), `export_bird_tree.py` (Jetz +
  AVONET, 9993 sp.), `export_fish_tree.py` (Rabosky + Phenoscape),
  `export_vertebrate_tree.py` (taxonomic splice of all three, 16,672 sp.).
  `verify_bird_tree.py` is the Playwright verification pattern — copy it.
- Molecular lens pilot (`esmc_prestin_pilot.py`): Prestin/SLC26A5 → ESM-C 300M
  (EvolutionaryScale, installed locally; weights cached) → ancestry subtracted via
  the dated mammal tree → permutation test. Honest negative at 14 taxa (p=0.13).

## KEY MODULES — reuse, do not reinvent
- `polarity.py` — Newick parse, Fitch `analyze`/`count_origins`,
  `randomization_test` (PTP), and **`dollo_recode`** (drop primitive-absent taxa).
- `mk.py` — `excess_homoplasy_test` (neutral-Mk on dated branches) AND
  **`trait_environment_test`** / `phi_coefficient` (phylogenetically-controlled
  trait↔environment association — the existing hook for "shared cause").
- `chronogram.py` — `load_chronogram` (fish), `load_mammal_chronogram`,
  `load_bird_chronogram`, `prune_to`, `normalize_unit_depth`, `branch_length_dict`.
- `ingest.py` — `pantheria()`, `avonet()`, `reconcile()`/`lineage()` (Open Tree).
- `phenoscape_ingest.py` / `phenoscape_polarity.py` — fish anatomy traits.
- `convergence_engine.py` — `ConvergenceModel`: Yoneda `deep_convergence`,
  `ricci_bridges` (Ollivier-Ricci), `spectral_embedding`, `detect`.
- `cache/` is gitignored (re-fetchable). ESM-C/torch/Biopython/Playwright are
  installed; no GPU (CPU is fine for tens of sequences). NCBI Entrez email is set.

## NON-NEGOTIABLE HONESTY CONTRACT
- Real data only; missing = greyed/dropped, never imputed.
- Polarity is reconstructed (Fitch decides ancestral state), not assumed; losses
  get Dollo-recoded so a "never had it" baseline can't fake convergence.
- Every claim goes through Fitch origins → PTP → dated-Mk → tier. Tier ≥3 = excess
  beyond drift; lower tiers are downgraded, not discarded or inflated.
- NOTHING is "done" until headless-verified (the last session shipped two broken
  viewers that "ran fine" but rendered the wrong data). Verify, then claim.

## TRACKS (pick with the user; A is the cheapest first win)

### A. Cross-class convergence overlay on the vertebrate tree
Code a trait that recurs across classes and run the gauntlet on a CROSS-CLASS
dated tree. The cleanest fully-in-data example is **aerial / gliding locomotion**:
bats (Chiroptera), birds (powered flight), and flying fish (Exocoetidae, gliding)
— all three classes are already in the loaded trees.
- HONEST CAVEAT (important): `vertebrate_tree.html` is a *taxonomic* splice
  (hierarchy only). Fitch origin-counting works on topology, but PTP/dated-Mk need
  real branch lengths SPANNING the classes. So this track requires grafting the
  three chronograms onto a deep-node-calibrated backbone (use TimeTree / DateLife
  ages for the Actinopterygii–Tetrapoda–Mammalia/Aves splits), then running the
  gauntlet on that dated cross-class tree. That backbone graft IS the real work;
  budget for it rather than running Mk on a branch-length-less splice.

### B. Environment-as-driver layer (Phase 5) — the most thesis-aligned endpoint
Turn "independent origins" into "independent origins WITH A SHARED CAUSE."
- GBIF occurrence API → per-species lat/long points → sample WorldClim/CHELSA (and
  PaleoClim for deep time) rasters → per-species climate niche.
- Then use the EXISTING `mk.trait_environment_test`: does the derived trait
  co-occur with the same environmental shift MORE than a neutral-Mk null predicts?
- Start small: one trait already overlaid (e.g. fossorial mammals, or marine) vs a
  single climate axis, on the clade where you already have the dated tree.

### C. Understand and improve the ESM-C molecular lens
The pilot returned a clean negative; the goal now is to UNDERSTAND why before
scaling, then improve it.
- Diagnose: load `results/prestin_esmc.json` + the cached embeddings; visualize the
  embedding space (PCA/UMAP), and compare against a RAW sequence-identity baseline —
  does ESM add anything over % identity for this protein? (If not, say so.)
- Localize: mean-pooling washes out convergent sites. Switch to PER-RESIDUE
  embeddings at the KNOWN convergent positions (Prestin "N7" group; Liu/Li 2010),
  i.e. test the active-site residues, not the whole-protein average.
- Scale only if justified: ESM-C 600M/6B (Forge API), and MULTIPLE convergence
  proteins (opsins / spectral tuning, antifreeze glycoproteins, toxin-resistant Na
  channels) — one protein is an anecdote.
- Keep the contract: ESM is a trait/function lens feeding the SAME gauntlet; never
  use embedding distance to build trees, and always subtract ancestry.

## NEW IDEAS — not yet in EXPANSION_PLAN.md (raise with the user)
1. **Posterior integration.** Every current verdict rests on ONE tree (the bird
   tree is a single pseudo-posterior sample). Run the gauntlet across N≈100 sampled
   posterior trees and report origins/PTP/Mk as DISTRIBUTIONS. Biggest single
   honesty upgrade; the VertLife/posterior archives are already the download source.
2. **A standing verification harness** (`verify_all.py` in CI) that headless-loads
   every viewer, asserts zero errors + biologically sane origin counts. This is the
   structural fix for the "looked done, was broken" failure mode.
3. **Multiple-testing discipline.** You test many traits; report effect sizes and
   FDR-controlled q-values, and raise the Monte-Carlo budgets (PTP/Mk are currently
   modest, ~199–499) so the "everything lands at T2" result is properly powered.
4. **Sensitivity reporting.** Show how each verdict moves with tree choice, trait
   threshold, and clade-coding decisions — make the robustness visible, not implicit.

## DEFINITION OF DONE
For any track: real data in, full gauntlet run, output headless-verified (zero page
errors, sane numbers), honesty caveats written into the viewer/README, and the
result stated plainly — including negatives. Then update `EXPANSION_PLAN.md` status.

Start by reading the four files above and confirming with the user which track
(A/B/C) and which of the new ideas to fold in.
