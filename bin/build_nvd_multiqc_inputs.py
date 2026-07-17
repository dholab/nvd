#!/usr/bin/env python3
"""Thin CLI for building NVD-owned MultiQC inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from py_nvd.multiqc_report import build_multiqc_inputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved-reads", type=Path, required=True)
    parser.add_argument("--nvd-version", type=Path, required=True)
    parser.add_argument("--fastqc-root", type=Path, required=True)
    parser.add_argument(
        "--experimental-enabled",
        choices=("true", "false"),
        required=True,
    )
    parser.add_argument("--output-dir", type=Path, default=Path("nvd_inputs"))
    args = parser.parse_args()

    build_multiqc_inputs(
        roster_path=args.resolved_reads,
        version_path=args.nvd_version,
        fastqc_root=args.fastqc_root,
        output_dir=args.output_dir,
        experimental_enabled=args.experimental_enabled == "true",
    )


if __name__ == "__main__":
    main()
