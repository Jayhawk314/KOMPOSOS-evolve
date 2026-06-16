"""Run the thesis-critical convergence scripts and archive their outputs.

This runner keeps the evidence trail separate from the analysis scripts. It
creates a timestamped directory under evolve/results/ with one stdout/stderr log
per script plus a machine-readable manifest.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EVOLVE = ROOT / "evolve"
RESULTS = EVOLVE / "results"

SCRIPTS = [
    "polarity.py",
    "mk.py",
    "anolis_validation.py",
    "echolocation_flagship.py",
    "geometry_test.py",
    "geometry_sensitivity.py",
    "horns_convergence.py",
    "horns_verified.py",
    "comparison_baselines.py",
    "phenoscape_convergence.py",
    "phenoscape_polarity.py",
    "phenoscape_homology.py",
    "phenoscape_inapplicability.py",
    "phenoscape_prediction.py",
    "phenoscape_significance.py",
    "phenoscape_excess.py",
    "phenoscape_excess_dated.py",
    "trait_environment.py",
    "reticulate_caveat.py",
]


def _run_script(script: str, out_dir: Path, timeout: int) -> dict[str, object]:
    script_path = EVOLVE / script
    started = datetime.now(timezone.utc)
    t0 = time.perf_counter()
    record: dict[str, object] = {
        "script": f"evolve/{script}",
        "started_utc": started.isoformat(),
        "timeout_seconds": timeout,
    }

    try:
        proc = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            env={
                **os.environ,
                "PYTHONUTF8": "1",
            },
        )
        record.update(
            {
                "returncode": proc.returncode,
                "duration_seconds": round(time.perf_counter() - t0, 3),
                "stdout_file": f"{script_path.stem}.out.txt",
                "stderr_file": f"{script_path.stem}.err.txt",
            }
        )
        (out_dir / f"{script_path.stem}.out.txt").write_text(
            proc.stdout, encoding="utf-8"
        )
        (out_dir / f"{script_path.stem}.err.txt").write_text(
            proc.stderr, encoding="utf-8"
        )
    except subprocess.TimeoutExpired as exc:
        record.update(
            {
                "returncode": None,
                "duration_seconds": round(time.perf_counter() - t0, 3),
                "timeout_expired": True,
                "stdout_file": f"{script_path.stem}.out.txt",
                "stderr_file": f"{script_path.stem}.err.txt",
            }
        )
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode("utf-8", errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode("utf-8", errors="replace")
        (out_dir / f"{script_path.stem}.out.txt").write_text(
            stdout, encoding="utf-8"
        )
        (out_dir / f"{script_path.stem}.err.txt").write_text(
            stderr, encoding="utf-8"
        )

    return record


def _write_summary(out_dir: Path, records: list[dict[str, object]]) -> None:
    passed = [r for r in records if r.get("returncode") == 0]
    failed = [r for r in records if r.get("returncode") != 0]
    lines = [
        "# Thesis Baseline Run",
        "",
        f"- Scripts: {len(records)}",
        f"- Passed: {len(passed)}",
        f"- Failed or timed out: {len(failed)}",
        "",
        "| Script | Status | Seconds | Output |",
        "|---|---:|---:|---|",
    ]
    for rec in records:
        status = rec.get("returncode")
        status_text = "PASS" if status == 0 else f"FAIL {status}"
        if rec.get("timeout_expired"):
            status_text = "TIMEOUT"
        lines.append(
            "| {script} | {status} | {seconds} | {stdout} |".format(
                script=rec["script"],
                status=status_text,
                seconds=rec.get("duration_seconds", ""),
                stdout=rec.get("stdout_file", ""),
            )
        )
    lines.append("")
    (out_dir / "SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--timeout",
        type=int,
        default=240,
        help="Per-script timeout in seconds.",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Optional suffix for the timestamped results directory.",
    )
    args = parser.parse_args()

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_name = stamp if not args.label else f"{stamp}_{args.label}"
    out_dir = RESULTS / run_name
    out_dir.mkdir(parents=True, exist_ok=False)

    manifest: dict[str, object] = {
        "run_name": run_name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "python": sys.version,
        "platform": platform.platform(),
        "scripts": [],
    }

    print(f"Writing thesis baseline to {out_dir}")
    records: list[dict[str, object]] = []
    for script in SCRIPTS:
        print(f"RUN evolve/{script}")
        rec = _run_script(script, out_dir, args.timeout)
        records.append(rec)
        print(
            "  -> returncode={returncode} seconds={duration_seconds}".format(
                returncode=rec.get("returncode"),
                duration_seconds=rec.get("duration_seconds"),
            )
        )

    manifest["scripts"] = records
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["passed"] = sum(1 for r in records if r.get("returncode") == 0)
    manifest["failed"] = sum(1 for r in records if r.get("returncode") != 0)
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    _write_summary(out_dir, records)

    print(f"Summary: {manifest['passed']} passed, {manifest['failed']} failed")
    return 0 if manifest["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
