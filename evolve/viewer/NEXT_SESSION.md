# Next Session — Handoff Prompt

Paste the block below into a fresh Claude Code session to continue work on the
convergence viewer. It is self-contained: it tells the new session what the
project is, which modules to reuse, the honesty contract, and the next task.

> **Phase 1 (Birds) is DONE.** `evolve/viewer/export_bird_tree.py` builds
> `bird_tree.html` (Jetz 2012 tree + AVONET, 9993 species, 8 overlays;
> nectarivory 48 origins, aerial 71, aquatic 35, vertivory 34, flightlessness 2).
> Verified headlessly by `verify_bird_tree.py`. Loaders: `avonet()` in
> `ingest.py`, `load_bird_chronogram()` in `chronogram.py`.
>
> The **default task below is now Phase M — the Prestin molecular pilot**. The
> original Phase 1 (Birds) block is kept at the bottom for reference.

---

You're working in the KOMPOSOS-evolve repo (`C:\Users\JAMES\github\komposos-evolve`,
branch `evolve-clean`). Read `evolve/viewer/EXPANSION_PLAN.md` and
`evolve/viewer/README.md` first — they describe the system and roadmap.

## CONTEXT

This project turns convergent evolution into a tiered, reproducible inference
problem. There's an interactive tree viewer in `evolve/viewer/` that renders the
real dated mammal tree (Upham 2019) joined with PanTHERIA traits, with a
"convergence overlay" that runs the real thesis pipeline (Fitch independent
origins -> PTP randomization -> neutral-Mk on dated branches -> tier verdict
T0-T3) and draws origin colors + cross-lineage bridges + a verdict badge.
Build it with: `python evolve/viewer/export_mammal_tree.py` (outputs a
self-contained `mammal_tree.html`). It currently has 8 overlays; the clade-coded
ones (marine=3 origins, myrmecophagy=5, fossorial=7, aerial=4, hopping=3)
recover the textbook numbers.

## KEY MODULES (reuse, don't reinvent)

- `evolve/polarity.py`   : Newick parse, Fitch (analyze, count_origins), randomization_test
- `evolve/mk.py`         : excess_homoplasy_test (neutral-Mk on dated branches)
- `evolve/chronogram.py` : load/prune/normalize dated trees, branch_length_dict
- `evolve/ingest.py`     : trait loaders (e.g. pantheria())
- `evolve/cache/`        : cached data (gitignored, re-fetchable). `mammal_mcc.nwk` is
  the mammal tree used; `pantheria.txt` the traits; `actinopt_*.tre` the fish tree.

## NON-NEGOTIABLE CONVENTIONS (this project is about honesty)

- Real data only. Missing values are shown as "no data" / greyed, never imputed.
- Clade-coded traits (coded by MSW05 family/order) are legitimate because family
  membership is external taxonomy, independent of the trait tested — label them
  as clade-coded, not measured.
- Every convergence claim goes through the full gauntlet (Fitch -> PTP ->
  dated-Mk -> tier). Don't ship a verdict that skips the nulls. Tier >= 3 means
  excess beyond drift; lower tiers are downgraded, not discarded.

## TASK: Phase 1 — add a BIRDS viewer using the same exporter pattern.

- Tree: Jetz et al. 2012 (birdtree.org) dated phylogeny.
- Traits: AVONET (morphometrics) + EltonTraits (diet/foraging).
- Reproduce the mammal viewer for birds: Order->Family->Genus->Species nav,
  trait coloring, and convergence overlays for known bird convergences
  (e.g. flightlessness, nectarivory, raptorial bill) — clade-coded where needed.
- Verify headlessly before declaring done: Python Playwright + Chromium are
  installed (`C:\Users\JAMES\AppData\Local\ms-playwright`). Load the HTML, select
  each overlay, assert zero page errors, bridges drawn, verdict shown.
- Confirm origin counts are biologically sane (flightlessness should show many
  independent origins across ratites/penguins/rails/etc.).

Start by reading `EXPANSION_PLAN.md` (Phase 1 row + Layer A/B sources) and
`export_mammal_tree.py` to mirror its structure. Keep the same honesty contract.

---

## ALTERNATIVE TASK: Phase M — Prestin molecular pilot (Layer D)

Pull SLC26A5 (Prestin) across the viewer's mammals, embed with ESM-C, test
whether bats + toothed whales cluster against a phylogenetic null — does the
molecular lens recover the echolocation convergence in
`echolocation_flagship.py`?
