# KOMPOSOS-evolve

**A pipeline that detects convergent evolution on real biological data — and rules
out the four things that masquerade as it — without fooling itself.**

Convergence is when nature solves the same problem twice, independently: sharks
and dolphins both became torpedoes; bats and dolphins both "see" with sound. The
hard part isn't spotting the resemblance — it's *proving it was earned twice,
not inherited once*. This engine does the proving, end to end, and is built so
it **cannot be rigged** to find what it hoped to find.

> **Start here:** [`evolve/README.md`](evolve/README.md) for the full pipeline
> (ingest → detect → predict → verify → the four-impostor gauntlet), then
> [`evolve/BOOK_OF_FINDINGS.md`](evolve/BOOK_OF_FINDINGS.md) for the plain-language
> version, no math required.

## The honest headline

The engine recovers real convergence from raw data (Anolis ecomorphs,
echolocation, fish fin-loss). But when pushed to the strongest claims, **it says
no** — and that is the point:

- Fin-loss convergence is real and structured, but **not beyond neutral drift**
  (confirmed with real dated fish trees).
- **Bergmann's rule is not supported** across mammals (confirmed with the real
  dated mammal tree), validated by positive *and* negative controls.

A tool that only ever confirms is measuring nothing. This one is willing to
reject — including the answers its author hoped for.

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

## Quick start

```bash
pip install requests pandas numpy scipy networkx

python evolve/polarity.py            # self-tests (no network, instant)
python evolve/anolis_validation.py   # recovers a known convergence from raw traits
python evolve/trait_environment.py   # Bergmann's rule on the dated mammal tree + controls
```

The `evolve/` package also shares a categorical substrate (`src/komposos_core/`)
with the rest of the KOMPOSOS family — see [`CLAUDE.md`](CLAUDE.md) for that
architecture; this README covers what's specific to convergent-evolution
analysis.

## License

Apache License 2.0. See `LICENSE`.

Author: James Ray Hawkins
