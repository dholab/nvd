"""Freeze a full-size e2e run's detections into tests/data/e2e_parity.json.

The reduced e2e samplesheet keeps one SRA run per code path. The recording this
script writes is what proves that reduction did not change what the pipeline
finds: the e2e compares its own detections against it.

Usage:

    pixi run -e dev python tests/scripts/record_e2e_parity.py <results_root>

where ``<results_root>`` is the ``results/nvd`` directory of a full-size run,
i.e. the one holding ``12_experiment_summary/experiment_blast_results.tsv``.

Re-record only from a run of the *full* five-row samplesheet. Recording from a
reduced run would make the check compare the reduction against itself, which
passes unconditionally and proves nothing.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from tests.e2e_parity import PARITY_PATH, summarize_detections  # noqa: E402

BLAST_RESULTS_RELATIVE_PATH = Path("12_experiment_summary") / "experiment_blast_results.tsv"


def read_blast_rows(results_root: Path) -> list[dict[str, str]]:
    """Read the experiment-wide blast result rows from a run's results tree."""
    blast_results_path = results_root / BLAST_RESULTS_RELATIVE_PATH
    if not blast_results_path.is_file():
        message = f"No blast results at {blast_results_path}"
        raise SystemExit(message)
    with blast_results_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def current_git_commit() -> str:
    """Record which commit produced the run, so a stale recording is traceable."""
    try:
        completed = subprocess.run(  # noqa: S603
            ["git", "-C", str(REPOSITORY_ROOT), "rev-parse", "HEAD"],  # noqa: S607
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return completed.stdout.strip()


def build_recording(
    results_root: Path,
    rows: list[dict[str, str]],
    detections: dict[str, dict[str, list[str]]],
) -> dict[str, object]:
    """Assemble the JSON document, provenance block first."""
    return {
        "recorded_from": {
            "results_root": str(results_root),
            "git_commit": current_git_commit(),
            "recorded_at": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "blast_result_rows": len(rows),
            "samples": sorted(detections),
        },
        "detections": detections,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "results_root",
        type=Path,
        help="The results/nvd directory of a full-size e2e run.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PARITY_PATH,
        help=f"Where to write the recording (default: {PARITY_PATH}).",
    )
    args = parser.parse_args()

    results_root = args.results_root.resolve()
    rows = read_blast_rows(results_root)
    detections = summarize_detections(rows)
    recording = build_recording(results_root, rows, detections)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(recording, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Recorded {len(rows)} rows from {len(detections)} samples")  # noqa: T201
    for sample_id, detection in sorted(detections.items()):
        print(f"  {sample_id}: taxids={detection['taxids']}")  # noqa: T201
    print(f"Wrote {args.out}")  # noqa: T201


if __name__ == "__main__":
    main()
