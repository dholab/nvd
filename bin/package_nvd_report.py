#!/usr/bin/env python3
"""Materialize one internal NVD report input package."""

from __future__ import annotations

import argparse
from pathlib import Path

from py_nvd.multiqc_packages import (
    DeaconReceipt,
    Domain,
    EligibilityReceipt,
    LongReadEligibilityReceipt,
    PreparedQueryBatchesReceipt,
    ProfileReceipt,
    UnionReceipt,
    write_package,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    deacon = commands.add_parser("deacon")
    deacon.add_argument(
        "--domain",
        type=Domain,
        choices=(Domain.TARGET_ENRICHMENT, Domain.DEPLETION),
        required=True,
    )
    deacon.add_argument("--sample-id", required=True)
    deacon.add_argument("--query-class", default="")
    deacon.add_argument("--stats", type=Path, required=True)

    profile = commands.add_parser("profile")
    profile.add_argument(
        "--domain",
        type=Domain,
        choices=(Domain.FASTX, Domain.ASSEMBLY),
        required=True,
    )
    profile.add_argument("--sample-id", required=True)
    profile.add_argument("--stage", required=True)
    profile.add_argument("--query-class", default="")
    profile.add_argument("--producer", default="")
    profile.add_argument("--profile", type=Path, required=True)
    profile.add_argument("--length-histogram", type=Path, required=True)
    profile.add_argument("--quality-histogram", type=Path)

    eligibility = commands.add_parser("eligibility")
    eligibility.add_argument("--sample-id", required=True)
    eligibility.add_argument("--producer", required=True)
    eligibility.add_argument("--decision", choices=("run", "skip"), required=True)
    eligibility.add_argument("--sequence-count", type=int, required=True)
    eligibility.add_argument("--threshold", type=int, required=True)
    eligibility.add_argument("--reason", required=True)

    long_read_eligibility = commands.add_parser("long-read-eligibility")
    long_read_eligibility.add_argument("--sample-id", required=True)
    long_read_eligibility.add_argument("--summary", type=Path, required=True)

    union = commands.add_parser("union")
    union.add_argument("--sample-id", required=True)
    union.add_argument("--summary", type=Path, required=True)

    prepared_query_batches = commands.add_parser("prepared-query-batches")
    prepared_query_batches.add_argument("--sample-id", required=True)
    prepared_query_batches.add_argument("--summary", type=Path, required=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "deacon":
        receipt = DeaconReceipt(
            domain=args.domain,
            sample_id=args.sample_id,
            query_class=args.query_class or None,
        )
        payloads = (args.stats,)
    elif args.command == "profile":
        receipt = ProfileReceipt(
            domain=args.domain,
            sample_id=args.sample_id,
            stage=args.stage,
            query_class=args.query_class or None,
            producer=args.producer or None,
            quality_histogram=(
                "quality_histogram.payload"
                if args.quality_histogram is not None
                else None
            ),
        )
        payloads = (args.profile, args.length_histogram)
        if args.quality_histogram is not None:
            payloads = (*payloads, args.quality_histogram)
    elif args.command == "eligibility":
        receipt = EligibilityReceipt(
            sample_id=args.sample_id,
            producer=args.producer,
            decision=args.decision,
            sequence_count=args.sequence_count,
            threshold=args.threshold,
            reason=args.reason,
        )
        payloads = ()
    elif args.command == "long-read-eligibility":
        receipt = LongReadEligibilityReceipt(sample_id=args.sample_id)
        payloads = (args.summary,)
    elif args.command == "union":
        receipt = UnionReceipt(sample_id=args.sample_id)
        payloads = (args.summary,)
    else:
        receipt = PreparedQueryBatchesReceipt(sample_id=args.sample_id)
        payloads = (args.summary,)
    write_package(receipt, payloads, output_root=Path.cwd())


if __name__ == "__main__":
    main()
