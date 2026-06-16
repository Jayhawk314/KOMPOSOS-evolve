# Working Session Summary

## Current objective

Execute the doctoral thesis plan for the `evolve/` convergence system, starting
with a reproducible evidence baseline.

## Completed so far

- Explored the repository structure and identified `evolve/` as the active
  evolutionary biology subsystem.
- Ran the core self-tests and demonstrations manually:
  - `polarity.py`
  - `mk.py`
  - `anolis_validation.py`
  - `echolocation_flagship.py`
  - `geometry_test.py`
  - `horns_convergence.py`
  - `horns_verified.py`
  - `trait_environment.py`
  - Phenoscape convergence, polarity, homology, inapplicability, significance,
    excess, and dated-excess analyses.
- Used web research to refresh the thesis positioning around convergent
  evolution, phylogenetic comparative methods, SURFACE/OU approaches, homoplasy,
  adaptive radiation, complex adaptive systems, evolvability, and reticulate
  evolution cautions.
- Added `evolve/THESIS_BLUEPRINT.md`, a doctoral-style prospectus that defines:
  - the core thesis,
  - chapter structure,
  - the tiered convergence ladder,
  - script-backed evidence,
  - limitations,
  - and a phased execution plan.
- Added `evolve/run_thesis_baseline.py`, a runner that archives outputs from the
  thesis-critical scripts into timestamped directories under `evolve/results/`.
- Ran the baseline runner once inside the normal sandbox:
  - results directory: `evolve/results/20260613T211322Z_initial/`
  - outcome: 10 passed, 5 failed
  - failure cause: Open Tree network/proxy failures in Phenoscape scripts that
    require taxon reconciliation.
- Reran the baseline runner with network approval:
  - results directory: `evolve/results/20260613T211444Z_network/`
  - outcome: 15 passed, 0 failed
  - archive size: about 36 KB
  - `SUMMARY.md` and `manifest.json` were written with per-script status,
    duration, stdout file, and stderr file.
- Added `evolve/comparison_baselines.py`, which compares the thesis methods to
  simpler baselines:
  - Anolis morphology-only AUROC: `0.9879`
  - Anolis morphology-minus-ancestry AUROC: `0.9942`
  - horn nearest-neighbor AUROC: `0.9134`
  - horn-max AUROC: `0.9391`
  - horn gain over nearest-neighbor: `+0.0257`
  - horn gain over trait frequency: `+0.2900`
- Updated `evolve/run_thesis_baseline.py` so future archives include
  `comparison_baselines.py`.
- Added `evolve/geometry_sensitivity.py`, which sweeps Anolis similarity
  thresholds and tests whether the geometry claim depends on the fixed `0.30`
  threshold:
  - thresholds tested: `0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50`
  - spectral ecomorph-closeness passed: `7/7`
  - Ricci cross-island bridge negativity passed: `7/7`
- Updated `evolve/run_thesis_baseline.py` so future archives include
  `geometry_sensitivity.py`.
- A run that was interrupted at the conversation layer had already completed:
  - results directory: `evolve/results/20260613T212713Z_baselines/`
  - outcome: 16 passed, 0 failed
  - includes the first Phase 2 comparison baseline.
- Updated `evolve/THESIS_BLUEPRINT.md` to record the new comparison baseline
  evidence, geometry sensitivity evidence, and mark the first Phase 2 baselines
  as added.
- Reran the full thesis baseline after adding both Phase 2 scripts:
  - sandbox/offline attempt: `evolve/results/20260613T213749Z_phase2/`
  - outcome: 12 passed, 5 failed
  - failure cause: same Open Tree network/proxy restriction for uncached
    Phenoscape reconciliation.
- Reran the full thesis baseline with network approval:
  - results directory: `evolve/results/20260613T214829Z_phase2_network/`
  - outcome: 17 passed, 0 failed
- Added `evolve/phenoscape_prediction.py`, a real-data relation-completion
  experiment beyond Anolis:
  - all-state nearest-profile AUROC: `0.8867`
  - all-state cross-order-minus-ancestry AUROC: `0.7499`
  - absence/loss nearest-profile AUROC: `0.8526`
  - absence/loss cross-order-minus-ancestry AUROC: `0.7909`
  - interpretation: morphology context predicts held-out Phenoscape states, but
    without explicit niche/environment intermediates the simple nearest-profile
    baseline remains strongest.
- Updated `evolve/run_thesis_baseline.py` so future archives include
  `phenoscape_prediction.py`.
- Reran the full thesis baseline with Phase 3 included:
  - results directory: `evolve/results/20260613T220554Z_phase3_network/`
  - outcome: 18 passed, 0 failed
  - archive size: about 41 KB
  - this was the official manuscript baseline until Phase 4 superseded it (see below).
- Added `evolve/results/README.md` to index evidence archives and identify the
  current official baseline.
- Added `evolve/reticulate_caveat.py` for Phase 4:
  - simulates true convergence versus one innovation plus reticulate transfer;
  - first two scenarios have the same tree and same tip states;
  - tree-only pipeline reports identical outputs for those two histories:
    `2` derived origins and `convergent=True`;
  - interpretation: MRCA/Fitch/Mk claims are conditional on the species tree
    being the relevant trait-inheritance history.
- Updated `evolve/run_thesis_baseline.py` so future archives include
  `reticulate_caveat.py`.
- Updated `evolve/THESIS_BLUEPRINT.md` with Phase 4 status and interpretation.
- Reran the full thesis baseline with Phase 4 included:
  - results directory: `evolve/results/20260613T225636Z_phase4_network/`
  - outcome: 19 passed, 0 failed
  - archive size: about 44 KB
  - this is the current official manuscript baseline.
- Updated `evolve/results/README.md` to mark the Phase 4 run as official.
- Started Phase 5 manuscript drafting:
  - added `evolve/DISSERTATION_DRAFT.md`
  - draft includes abstract, central claim, tier ladder, evidence table, chapter
    scaffold, chapter-level interpretations, boundaries, conclusion, and
    manuscript to-do list.
  - numeric claims are drawn from
    `evolve/results/20260613T225636Z_phase4_network/`.

## In progress

- Phase 1 is complete.
- Phase 2 has comparison and geometry-threshold baselines implemented and
  archived.
- Phase 3 has started with `phenoscape_prediction.py` and is archived in the
  current official baseline.
- Phase 4 has an executable reticulate-history caveat and is archived in the
  current official baseline.
- Phase 5 has started with `DISSERTATION_DRAFT.md`.

## Next steps

- Expand the related-work chapter with full citations and checked DOI metadata.
- Convert the evidence table into polished manuscript tables.
- Consider whether to keep all earlier result archives or prune failed/partial
  runs before manuscript packaging.
