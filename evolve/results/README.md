# Thesis Evidence Archives

Each timestamped directory was produced by `evolve/run_thesis_baseline.py`.
Every run contains:

- `SUMMARY.md`: human-readable pass/fail table.
- `manifest.json`: machine-readable run metadata.
- `*.out.txt`: captured stdout for each script.
- `*.err.txt`: captured stderr for each script.

## Current Official Baseline

`20260613T225636Z_phase4_network/`

- scripts: 19
- passed: 19
- failed: 0
- includes:
  - Phase 1 evidence freeze,
  - Phase 2 comparison baselines,
  - Phase 2 geometry threshold sensitivity,
  - Phase 3 Phenoscape relation-completion experiment,
  - Phase 4 reticulate-history caveat.

## Earlier Runs

- `20260613T211322Z_initial/`: sandbox run; 10 passed, 5 Open Tree network failures.
- `20260613T211444Z_network/`: network-approved run; 15 passed, 0 failed.
- `20260613T212713Z_baselines/`: included `comparison_baselines.py`; 16 passed, 0 failed.
- `20260613T213749Z_phase2/`: sandbox run; 12 passed, 5 Open Tree network failures.
- `20260613T214829Z_phase2_network/`: included comparison and geometry sensitivity; 17 passed, 0 failed.
- `20260613T220554Z_phase3_network/`: included Phenoscape prediction; 18 passed, 0 failed.
