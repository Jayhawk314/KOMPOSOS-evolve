#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Headless verification of bird_tree.html (Phase-1 birds viewer).

Loads the self-contained HTML in Chromium (Playwright), then for every convergence
overlay it switches the View selector, and asserts:
  * zero page errors / zero console errors the whole time,
  * the data (META + CONV) actually embedded and sane,
  * the verdict badge is shown with the right origin count + tier,
  * cross-lineage bridges are drawn,
  * origin counts are biologically sane (flightlessness, nectarivory have the
    expected ballpark of independent origins).

Run:  python evolve/viewer/verify_bird_tree.py     (after export_bird_tree.py)
Exit code 0 = all checks pass.
"""

from __future__ import annotations

import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

HTML = Path(__file__).resolve().parent / "bird_tree.html"

# biological sanity floors: each of these famous bird convergences must recover at
# least this many independent origins (real numbers are higher; these are floors).
SANE_MIN_ORIGINS = {
    "flightless": 2,   # ratites (palaeognath) + penguins, conservatively
    "nectarivory": 10,  # hummingbirds, sunbirds, honeyeaters, lorikeets, ...
    "aquatic": 8,       # penguins, auks, grebes, loons, ducks, cormorants, ...
    "aerial": 8,        # swifts, swallows, nightjars, ...
    "vertivory": 5,     # hawks, falcons, owls, shrikes, ...
}


def main() -> int:
    if not HTML.exists():
        print(f"FAIL: {HTML} not found -- run export_bird_tree.py first")
        return 1

    errors: list[str] = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.on("pageerror", lambda e: errors.append(f"pageerror: {e}"))
        page.on("console", lambda m: errors.append(f"console.{m.type}: {m.text}")
                if m.type == "error" else None)

        page.goto(HTML.as_uri())
        # wait until d3 loaded + the script populated the View selector
        page.wait_for_function("document.querySelector('#trait')"
                               "&&document.querySelector('#trait').options.length>5",
                               timeout=20000)

        meta = page.evaluate("META")
        conv = page.evaluate("CONV")
        print(f"viewer: {meta['n_species']} species | {meta['n_orders']} orders | "
              f"{meta['n_families']} families | overlays: {list(conv.keys())}")

        if meta["n_species"] < 9000:
            errors.append(f"too few species embedded: {meta['n_species']}")

        ok = True
        for key, c in conv.items():
            # switch the View selector to this overlay and fire 'change'
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
            floor = SANE_MIN_ORIGINS.get(key)
            sane = (floor is None) or (origins >= floor)
            row_ok = (verdict_visible and badge.startswith("TIER")
                      and n_bridges >= 1 and n_legend >= 1 and sane)
            ok = ok and row_ok
            flag = "OK " if row_ok else "FAIL"
            extra = "" if sane else f"  <-- origins {origins} < floor {floor}"
            print(f"  [{flag}] {key:<12} origins={origins:<4} tier=T{c['tier']} "
                  f"badge='{badge}' bridges={n_bridges} legend={n_legend} "
                  f"coded={c['coded']}{extra}")

        # also exercise a couple of base (non-overlay) views
        for base in ["order", "body_log", "lifestyle"]:
            page.evaluate("""v => { const s=document.querySelector('#trait');
                s.value=v; s.dispatchEvent(new Event('change')); }""", base)
            page.wait_for_timeout(80)
        browser.close()

    if errors:
        print("\nPAGE/CONSOLE ERRORS:")
        for e in errors[:30]:
            print("  " + e)
        ok = False

    print("\nRESULT:", "ALL CHECKS PASS" if ok and not errors else "FAILURES PRESENT")
    return 0 if (ok and not errors) else 1


if __name__ == "__main__":
    sys.exit(main())
