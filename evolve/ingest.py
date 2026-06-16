"""
Real-data ingestion for the convergence engine. NO hand-curation of traits.

Two public sources, joined on taxon name inside the Category (which is itself
SQLite-backed -- the join you asked for):

  * Open Tree of Life  (REST API)      -> ancestry lens (real lineage / clades)
  * PanTHERIA          (one bulk file) -> morphology / ecology / environment lenses

The only thing a user supplies is the SAMPLE: a list of species names to study.
That is not curating the answer -- the traits, the tree, and the environment all
come from the databases. Continuous traits are quantile-binned across the sample
so that "shared trait state" becomes a shared object in the Category.

Caches everything under evolve/cache/ so re-runs are offline.
"""

from __future__ import annotations

import json
import pathlib
import time

import pandas as pd
import requests

from convergence_engine import ConvergenceModel

CACHE = pathlib.Path(__file__).resolve().parent / "cache"
CACHE.mkdir(exist_ok=True)
OTT = "https://api.opentreeoflife.org/v3"
PANTHERIA_URL = "https://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR05_Aug2008.txt"
MISSING = -999.0  # PanTHERIA's missing-value sentinel


# ----------------------------------------------------------------------------
# Open Tree of Life
# ----------------------------------------------------------------------------
def reconcile(names: list[str]) -> dict[str, int]:
    """Map species names -> OTT ids via Open Tree's name-resolution service."""
    cache = CACHE / "ott_ids.json"
    known = json.loads(cache.read_text()) if cache.exists() else {}
    todo = [n for n in names if n not in known]
    if todo:
        r = requests.post(f"{OTT}/tnrs/match_names", json={"names": todo}, timeout=30)
        r.raise_for_status()
        for res in r.json()["results"]:
            if res["matches"]:
                known[res["name"]] = res["matches"][0]["taxon"]["ott_id"]
        cache.write_text(json.dumps(known, indent=2))
    return {n: known[n] for n in names if n in known}


def lineage(ott_id: int) -> list[str]:
    """Real ancestral clade names for an OTT id (order/family/.. up to root)."""
    cache = CACHE / "lineages.json"
    known = json.loads(cache.read_text()) if cache.exists() else {}
    key = str(ott_id)
    if key not in known:
        r = requests.post(f"{OTT}/taxonomy/taxon_info",
                          json={"ott_id": ott_id, "include_lineage": True}, timeout=30)
        r.raise_for_status()
        anc = r.json().get("lineage", [])
        # keep only major ranks so the ancestry lens is meaningful, not noisy
        keep = {"genus", "family", "superfamily", "order", "superorder",
                "infraclass", "class"}
        known[key] = [a["name"] for a in anc if a.get("rank") in keep]
        cache.write_text(json.dumps(known, indent=2))
        time.sleep(0.2)  # be polite to the API
    return known[key]


# ----------------------------------------------------------------------------
# PanTHERIA
# ----------------------------------------------------------------------------
def pantheria() -> pd.DataFrame:
    f = CACHE / "pantheria.txt"
    if not f.exists():
        f.write_bytes(requests.get(PANTHERIA_URL, timeout=60).content)
    df = pd.read_csv(f, sep="\t")
    return df.set_index("MSW05_Binomial")


# ----------------------------------------------------------------------------
# AVONET (Tobias et al. 2022, Ecology Letters) -- bird traits
# ----------------------------------------------------------------------------
AVONET_URL = "https://ndownloader.figshare.com/files/34480856"  # figshare suppl. 1


def avonet() -> pd.DataFrame:
    """AVONET trait table in the BirdTree (Jetz et al. 2012) taxonomy, indexed by
    binomial. Real measured morphometrics (Mass, beak, wing, Hand-Wing.Index, ...)
    plus ecology (Habitat, Migration, Trophic.Level, Trophic.Niche,
    Primary.Lifestyle) for all 9993 BirdTree species -- a 1:1 join to the bird
    chronogram tips (which use the same Genus_species names)."""
    x = CACHE / "avonet.xlsx"
    if not x.exists():
        x.write_bytes(requests.get(AVONET_URL, timeout=300).content)
    df = pd.read_excel(x, sheet_name="AVONET3_BirdTree")
    return df.set_index("Species3")


# Which PanTHERIA columns feed which lens, and a short tag for each.
LENS_COLUMNS = {
    "morphology": {
        "5-1_AdultBodyMass_g": "bodymass",
        "8-1_AdultForearmLen_mm": "forearm",
        "13-1_AdultHeadBodyLen_mm": "headbody",
    },
    "ecology": {
        "1-1_ActivityCycle": "activity",
        "6-2_TrophicLevel": "trophic",
        "6-1_DietBreadth": "dietbreadth",
        "12-1_HabitatBreadth": "habitat",
        "12-2_Terrestriality": "terrestrial",
    },
    "environment": {  # the DRIVER: real climate niche
        "28-2_Temp_Mean_01degC": "temp",
        "28-1_Precip_Mean_mm": "precip",
        "30-1_AET_Mean_mm": "aet",
        "26-4_GR_MidRangeLat_dd": "abslat",
    },
}
# columns that are already categorical (do not quantile-bin)
CATEGORICAL = {"activity", "trophic", "terrestrial"}


def _bin_series(tag: str, values: dict[str, float], n: int = 3) -> dict[str, str]:
    """Turn a column of numeric values into shared discrete state-labels.

    Categorical columns map their integer code directly; continuous columns are
    cut into n quantile bins (low/mid/high...) ACROSS THE SAMPLE so that two
    taxa in the same bin share a state object.
    """
    present = {k: v for k, v in values.items() if v != MISSING and pd.notna(v)}
    if not present:
        return {}
    if tag in CATEGORICAL:
        return {k: f"{tag}={int(v)}" for k, v in present.items()}
    s = pd.Series(present)
    if tag == "abslat":
        s = s.abs()
    try:
        labels = [f"{tag}={lab}" for lab in (["lo", "mid", "hi"][:n] if n == 3
                                             else list(range(n)))]
        binned = pd.qcut(s, q=min(n, s.nunique()), labels=False, duplicates="drop")
        return {k: f"{tag}=q{int(b)}" for k, b in binned.items()}
    except Exception:
        return {}


def build_model(names: list[str], verbose: bool = True) -> tuple[ConvergenceModel, dict]:
    """Build a multi-lens ConvergenceModel for `names` purely from real data."""
    ids = reconcile(names)
    df = pantheria()
    present = [n for n in names if n in ids and n in df.index]
    dropped = [n for n in names if n not in present]

    m = ConvergenceModel()
    for n in present:
        m.add_taxon(n)

    # ancestry lens from Open Tree lineage
    for n in present:
        for clade in lineage(ids[n]):
            m.relate("ancestry", n, clade, rel="member")

    # trait lenses from PanTHERIA, quantile-binned across the sample
    report = {"present": present, "dropped": dropped, "lens_states": {}}
    for lens, cols in LENS_COLUMNS.items():
        for col, tag in cols.items():
            raw = {n: df.at[n, col] for n in present if col in df.columns}
            states = _bin_series(tag, raw)
            for n, state in states.items():
                m.relate(lens, n, state, rel="has")
            if states:
                report["lens_states"].setdefault(lens, []).append(tag)

    if verbose:
        print(f"Ingested {len(present)}/{len(names)} taxa "
              f"({len(dropped)} not found: {dropped[:5]}{'...' if len(dropped)>5 else ''})")
        for lens in ["ancestry", "morphology", "ecology", "environment"]:
            tags = report["lens_states"].get(lens, []) if lens != "ancestry" else ["lineage"]
            print(f"  {lens:<12} <- {tags}")
    return m, report
