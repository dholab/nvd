#!/usr/bin/env python3
"""Write the receipt for one successfully completed raw FastQC unit."""

from __future__ import annotations

import argparse
from pathlib import Path

from py_nvd.multiqc_report import FastqcReceipt, ReadEnd, ReadSource


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--source", type=ReadSource, required=True)
    parser.add_argument("--read-end", type=ReadEnd, required=True)
    parser.add_argument("--input-ordinal", type=int, required=True)
    parser.add_argument("--alias", required=True)
    parser.add_argument("--zip-alias", required=True)
    parser.add_argument("--html-alias", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    receipt = FastqcReceipt(
        schema_version="nvd.fastqc-raw-unit/v1",
        sample_id=args.sample_id,
        platform=args.platform,
        source=args.source,
        read_end=args.read_end,
        input_ordinal=args.input_ordinal,
        alias=args.alias,
        zip_alias=args.zip_alias,
        html_alias=args.html_alias,
    )
    args.output.write_text(
        receipt.model_dump_json() + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
