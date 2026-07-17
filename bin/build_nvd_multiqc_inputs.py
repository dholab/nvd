#!/usr/bin/env python3
"""Thin CLI for building NVD-owned MultiQC inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from py_nvd.multiqc_report import (
    CompileRequest,
    ReportConfiguration,
    ReportRoots,
    build_multiqc_inputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved-reads", type=Path, required=True)
    parser.add_argument("--nvd-version", type=Path, required=True)
    parser.add_argument("--fastqc-root", type=Path, required=True)
    parser.add_argument("--target-enrichment-root", type=Path)
    parser.add_argument("--depletion-root", type=Path)
    parser.add_argument("--fastx-root", type=Path)
    parser.add_argument("--assembly-root", type=Path)
    parser.add_argument(
        "--experimental-enabled",
        choices=("true", "false"),
        required=True,
    )
    parser.add_argument("--output-dir", type=Path, default=Path("nvd_inputs"))
    parser.add_argument(
        "--target-enrichment-enabled",
        choices=("true", "false"),
        required=True,
    )
    parser.add_argument(
        "--depletion-enabled",
        choices=("true", "false"),
        required=True,
    )
    parser.add_argument(
        "--assembly-enabled",
        choices=("true", "false"),
        required=True,
    )
    args = parser.parse_args()

    build_multiqc_inputs(
        CompileRequest(
            roster_path=args.resolved_reads,
            version_path=args.nvd_version,
            fastqc_root=args.fastqc_root,
            output_dir=args.output_dir,
            configuration=ReportConfiguration(
                experimental_enabled=args.experimental_enabled == "true",
                target_enrichment_enabled=args.target_enrichment_enabled == "true",
                depletion_enabled=args.depletion_enabled == "true",
                assembly_enabled=args.assembly_enabled == "true",
            ),
            report_roots=ReportRoots(
                target_enrichment=args.target_enrichment_root,
                depletion=args.depletion_root,
                fastx=args.fastx_root,
                assembly=args.assembly_root,
            ),
        ),
    )


if __name__ == "__main__":
    main()
