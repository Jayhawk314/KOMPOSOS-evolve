#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Standing headless verification of EVERY viewer (Idea 1 / CI harness).

This is the structural fix for the "looked done, was broken" failure mode: a
viewer that opens without crashing can still render the WRONG data. Here, for each
built viewer, we load the self-contained HTML in Chromium (Playwright) and assert:

  * zero page errors / zero console errors throughout,
  * META embedded with a sane species count (per-viewer floor),
  * every convergence overlay: verdict badge shown + tier stamped, a legend drawn,
    cross-lineage bridges present WHEN there are >=2 independent origins, and the
    origin count clears a biological sanity floor,
  * a few base (non-overlay) colorings switch without error.

Generalizes `verify_bird_tree.py` to all viewers; same DOM/JS contract
(`META`, `CONV`, `#trait`, `#verdict`, `#bridges`, `#legend`, `conv:` prefix).

Run:  python evolve/viewer/verify_all.py            # all viewers
      python evolve/viewer/verify_all.py mammal fish  # a subset
Exit code 0 = all present viewers pass. Missing HTML is reported, not a hard fail
unless --strict is given (CI passes --strict so an un-exported viewer fails loudly).
"""

from __future__ import annotations

import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent

# Per-viewer config: HTML file, species floor, base colorings to exercise, and
# biological origin floors for named overlays (real numbers are >= these).
VIEWERS = {
    "mammal": {
        "html": "mammal_tree.html",
        "min_species": 3000,
        "origin_floors": {"marine": 3, "aerial": 2, "fossorial": 5, "myrmecophagy": 4,
                          "hopping": 2, "carnivory": 10, "diurnality": 10, "giant": 10},
    },
    "bird": {
        "html": "bird_tree.html",
        "min_species": 9000,
        "origin_floors": {"flightless": 2, "nectarivory": 10, "aquatic": 8,
                          "aerial": 8, "vertivory": 5},
    },
    "fish": {
        "html": "fish_tree.html",
        "min_species": 2000,
        # single-origin fins (pectoral/dorsal/caudal) are T0 -> floor 1, no bridges
        "origin_floors": {"pelvic_fin": 3, "adipose_fin": 4, "barbel": 4,
                          "pectoral_fin": 1, "dorsal_fin": 1, "caudal_fin": 1},
    },
    "vertebrate": {
        "html": "vertebrate_tree.html",
        "min_species": 16000,
        "origin_floors": {},   # cross-class overlays not added yet (Track A)
    },
}


def verify_one(name: str, cfg: dict) -> tuple[bool, bool]:
    """Returns (present, ok). present=False if the HTML is missing."""
    html = HERE / cfg["html"]
    if not html.exists():
        print(f"[{name}] MISSING: {cfg['html']} (run its exporter)")
        return False, False

    errors: list[str] = []
    ok = True
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.on("pageerror", lambda e: errors.append(f"pageerror: {e}"))
        page.on("console", lambda m: errors.append(f"console.{m.type}: {m.text}")
                if m.type == "error" else None)
        page.goto(html.as_uri())
        page.wait_for_function("typeof META!=='undefined' && typeof CONV!=='undefined'",
                               timeout=20000)

        meta = page.evaluate("META")
        conv = page.evaluate("CONV")
        n_sp = meta.get("n_species", 0)
        print(f"[{name}] {n_sp} species | overlays: {list(conv.keys()) or '(none)'}")
        if n_sp < cfg["min_species"]:
            errors.append(f"too few species: {n_sp} < {cfg['min_species']}")
            ok = False

        for key, c in conv.items():
            page.evaluate(
                """k => { const s=document.querySelector('#trait');
                    s.value='conv:'+k; s.dispatchEvent(new Event('change')); }""", key)
            page.wait_for_timeout(150)
            verdict_visible = page.evaluate(
                "getComputedStyle(document.querySelector('#verdict')).display!=='none'")
            badge = page.evaluate(
                "(document.querySelector('#verdict .badge')||{}).textContent||''")
            n_bridges = page.evaluate(
                "document.querySelectorAll('#bridges path.bridge').length")
            n_legend = page.evaluate("document.querySelectorAll('#legend .row').length")

            origins = c["origins"]
            floor = cfg["origin_floors"].get(key)
            sane = (floor is None) or (origins >= floor)
            # bridges only expected when there are >=2 origins to bridge between
            bridges_ok = (n_bridges >= 1) if origins >= 2 else True
            row_ok = (verdict_visible and badge.upper().startswith("TIER")
                      and bridges_ok and n_legend >= 1 and sane)
            ok = ok and row_ok
            flag = "OK " if row_ok else "FAIL"
            extra = "" if sane else f"  <-- origins {origins} < floor {floor}"
            print(f"   [{flag}] {key:<14} origins={origins:<4} tier=T{c['tier']} "
                  f"badge='{badge}' bridges={n_bridges} legend={n_legend}{extra}")

        # exercise the viewer's ACTUAL base colorings (options without a conv: prefix),
        # read from the DOM so we never feed a value the dropdown doesn't offer.
        bases = page.evaluate(
            "[...document.querySelector('#trait').options].map(o=>o.value)"
            ".filter(v=>!v.startsWith('conv:'))")
        for base in bases:
            page.evaluate("""v => { const s=document.querySelector('#trait');
                s.value=v; s.dispatchEvent(new Event('change')); }""", base)
            page.wait_for_timeout(60)
        print(f"   base colorings exercised: {bases}")
        browser.close()

    if errors:
        print(f"   PAGE/CONSOLE ERRORS ({name}):")
        for e in errors[:20]:
            print("     " + e)
        ok = False
    print(f"   => {name}: {'PASS' if ok else 'FAIL'}")
    return True, ok


def main(argv: list[str]) -> int:
    strict = "--strict" in argv
    wanted = [a for a in argv if not a.startswith("--")] or list(VIEWERS)
    any_fail = False
    missing = []
    for name in wanted:
        if name not in VIEWERS:
            print(f"unknown viewer '{name}' (known: {list(VIEWERS)})")
            any_fail = True
            continue
        present, ok = verify_one(name, VIEWERS[name])
        if not present:
            missing.append(name)
        elif not ok:
            any_fail = True

    print("\n" + "=" * 60)
    if missing:
        print(f"MISSING (not built): {', '.join(missing)}"
              + ("  -> FAIL (--strict)" if strict else "  (skipped)"))
    if any_fail or (strict and missing):
        print("RESULT: FAILURES PRESENT")
        return 1
    print("RESULT: ALL PRESENT VIEWERS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
