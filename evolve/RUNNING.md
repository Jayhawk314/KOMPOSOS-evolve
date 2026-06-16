# Running the engine

Every script is standalone and prints a self-explaining report. On first run each
fetches and **caches** its data under `evolve/cache/` (git-ignored); re-runs are
offline. Expect the first run of the dated-tree scripts to download large files
(fish chronogram ~0.2 MB compressed; mammal chronogram ~96 MB).

## Requirements

- Python 3.10+
- `requests`, `pandas`, `numpy`, `scipy`, `networkx`
- Internet access on first run (Open Tree, Phenoscape, PanTHERIA, Fish Tree of Life,
  PHYLACINE). All cached afterward.

Run from the repository root, e.g. `python evolve/polarity.py`.

## Self-tests first (no network, instant)

```
python evolve/polarity.py        # Fitch, homology, Dollo, randomization — hand-checkable toys
python evolve/mk.py              # neutral-Mk sanity on a clustered toy
```
These assert known answers on tiny trees; if they pass, the math core is sound.

## Synthetic validation (proves the method recovers a known answer)

```
python evolve/anolis_validation.py      # recovers Anolis ecomorph convergence from raw traits
python evolve/echolocation_flagship.py  # bat~dolphin convergence stacks across 3 levels; fruit-bat control
python evolve/geometry_test.py          # spectral + Ricci recover the same convergence
```

## Prediction & verification

```
python evolve/horns_convergence.py   # predict traits + leave-one-out AUROC 0.94
python evolve/horns_verified.py      # dual ZFC/CAT verdict on every prediction
```

## Real-data pipeline (the honest tests)

```
python evolve/phenoscape_ingest.py          # pull ontology-annotated fish anatomy
python evolve/phenoscape_convergence.py      # naive view (shows the plesiomorphy trap)
python evolve/phenoscape_polarity.py         # ancestral vs derived (parameter-free)
python evolve/phenoscape_homology.py         # inherited vs independent loss
python evolve/phenoscape_inapplicability.py  # lost vs never-had-it (Dollo)
python evolve/phenoscape_significance.py     # is it more than chance? (PTP)
python evolve/phenoscape_excess.py           # is it more than drift? (Mk, Grafen branch lengths)
python evolve/phenoscape_excess_dated.py     # same, with REAL fish chronogram branch lengths
python evolve/polarity_fold.py               # polarity folded into the similarity metric
python evolve/homology_fold.py               # homology folded into the similarity metric
```

## Environment as driver

```
python evolve/trait_environment.py   # Bergmann's rule on the dated mammal tree + controls
```

## Reading the output

Each script prints, in order: what it loaded, the per-character / per-pair numbers,
a verdict, and an "HONEST NOTES" / "HONEST READING" block stating the limits and any
place a negative result is reported. The negatives are not errors — they are the
engine doing its job.

## Reproducibility notes

- All randomness is seeded (`seed=`/`random.Random(n)`), so reruns are identical.
- Network results are cached, so a rerun reproduces the exact numbers offline.
- The dated-tree scripts set `sys.setrecursionlimit(200000)` to parse the large
  Newick trees.
