# The convergent-evolution engine

**A pipeline that detects convergent evolution on real biological data — and rules
out the four things that masquerade as it — without fooling itself.**

Convergence is when nature solves the same problem twice, independently: sharks and
dolphins both became torpedoes; bats and dolphins both "see" with sound. The hard
part isn't spotting the resemblance — it's *proving it was earned twice, not
inherited once*. This engine does the proving, end to end, and is built so it
**cannot be rigged** to find what we hoped.

> **Start here:** [`BOOK_OF_FINDINGS.md`](BOOK_OF_FINDINGS.md) — the whole thing
> explained in plain language, no math required. Then
> [`METHODS.md`](METHODS.md), [`RUNNING.md`](RUNNING.md),
> [`GLOSSARY.md`](GLOSSARY.md).

## The honest headline

The engine recovers real convergence from raw data (Anolis ecomorphs, echolocation,
fish fin-loss). But when pushed to the strongest claims, **it says no** — and that
is the point:

- Fin-loss convergence is real and structured, but **not beyond neutral drift**
  (confirmed with real dated fish trees).
- **Bergmann's rule is not supported** across mammals (confirmed with the real dated
  mammal tree), validated by positive *and* negative controls.

A tool that only ever confirms is measuring nothing. This one is willing to reject —
including the answers its author hoped for.

## The pipeline

```
ingest        Open Tree · PanTHERIA · Phenoscape · dated fish & mammal chronograms
detect        Yoneda relational similarity · spectral · Ricci bridges · homotopy routes
predict       simplicial horn-filling          leave-one-out AUROC 0.94
verify        dual ZFC/CAT engine              filters spurious 0.66 vs 0.14
─── the gauntlet: rule out the four impostors ───────────────────────────────
polarity      ancestral vs derived            (old news)        Fitch parsimony
homology      independent vs inherited        (inherited loss)  MRCA test
inapplicable  lost vs never-had-it            (category error)  Dollo recoding
significance  signal vs chance                (dumb luck)       randomization (PTP)
excess        convergence vs drift            (dumb luck)       neutral-Mk + dated trees
─── the real driver ─────────────────────────────────────────────────────────
environment   does the trait track the world beyond ancestry?  phylo-controlled test
```

Every method is parameter-free or estimated-under-the-null; every tree is built from
names and DNA, independent of the traits being tested; every test is bracketed by
controls that prove it can say both yes and no.

## Quick start

```bash
pip install requests pandas numpy scipy networkx

python evolve/polarity.py            # self-tests (no network, instant)
python evolve/anolis_validation.py   # recovers a known convergence from raw traits
python evolve/trait_environment.py   # Bergmann's rule on the dated mammal tree + controls
```

First run of each script downloads and caches its data under `evolve/cache/`
(git-ignored); re-runs are offline and reproduce the exact numbers.

## What's inside

- **Engine**: `convergence_engine.py` (the four detectors + folded-in polarity &
  homology), `polarity.py` (Fitch, MRCA homology, Dollo, randomization, dated-tree
  homoplasy), `mk.py` (neutral-Mk + trait-environment), `chronogram.py` (real dated
  trees), `categorical_verifier.py`, `horns_convergence.py` / `horns_verified.py`.
- **Data ingest**: `ingest.py` (Open Tree + PanTHERIA), `phenoscape_ingest.py`.
- **Demonstrations**: the `phenoscape_*` and `*_fold` scripts, each a self-explaining
  report.

## A note on provenance

This grew on a copy of a categorical-AI codebase whose documentation oversold its
abilities. Several inherited functions were broken (the core similarity returned zero
for everything; a geometry bridge was dead; a logic branch was unreachable). They
were found, fixed, or honestly flagged. Where the categorical framing did real work
(relational fingerprints, Ricci bridges, horn-prediction, tiered verification) it is
used; where it would have been decoration (e.g. using a category-theory tool to
estimate a statistical rate) it was refused. Skepticism was load-bearing throughout.

License: Apache-2.0 (see repository).
