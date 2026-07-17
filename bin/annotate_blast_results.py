#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

import argparse
import csv
import logging
import os
import sys
from dataclasses import dataclass

from py_nvd import taxonomy

# Set up logging
logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AnnotationRequest:
    input_file: str
    output_file: str
    sample_name: str
    task: str
    taxonomy_dir: str | None
    sync: bool
    taxonomy_mode: str | None
    taxonomy_max_age_days: int | None


def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate BLAST results with taxonomy info",
    )
    parser.add_argument("--input_file", required=True, help="Path to BLAST hits file")
    parser.add_argument("--output_file", required=True, help="Output TSV file")
    parser.add_argument("--sample_name", required=True, help="Sample name")
    parser.add_argument("--task", required=True, help="Task name")
    parser.add_argument(
        "--taxonomy-dir",
        required=False,
        default=None,
        help="Explicit taxonomy directory",
    )
    parser.add_argument(
        "--sync",
        action="store_true",
        help="Require pre-cached taxonomy database (fail if unavailable)",
    )
    parser.add_argument(
        "--taxonomy-mode",
        choices=[mode.value for mode in taxonomy.TaxonomyMode],
        default=None,
        help="Pipeline taxonomy mode: read_only or missing",
    )
    parser.add_argument(
        "--taxonomy-max-age-days",
        type=int,
        default=None,
        help="Freshness warning threshold for existing taxonomy",
    )

    return parser.parse_args()


def annotate_blast_results(request: AnnotationRequest) -> None:
    try:
        if not os.path.exists(request.input_file):
            logger.error(f"Input file not found: {request.input_file}")
            sys.exit(1)

        if os.path.getsize(request.input_file) == 0:
            logger.warning("Input file is empty. Creating an empty output file.")
            open(request.output_file, "w").close()
            return

        with (
            taxonomy.open(
                taxonomy_dir=request.taxonomy_dir,
                sync=request.sync,
                taxonomy_mode=request.taxonomy_mode,
                max_age_days=request.taxonomy_max_age_days,
            ) as tax,
            open(request.input_file) as infile,
            open(request.output_file, "w", newline="") as outfile,
        ):
            reader = csv.reader(infile, delimiter="\t")
            # Skip header line
            next(reader)
            writer = csv.writer(outfile, delimiter="\t")

            # Write header
            header = [
                "task",
                "sample",
                "qseqid",
                "qlen",
                "sseqid",
                "stitle",
                "length",
                "pident",
                "evalue",
                "bitscore",
                "sscinames",
                "staxids",
                "rank",
            ]
            writer.writerow(header)

            for row in reader:
                (
                    _qseqid,
                    _qlen,
                    _sseqid,
                    _stitle,
                    _length,
                    _pident,
                    _evalue,
                    _bitscore,
                    _sscinames,
                    staxids,
                ) = row

                # Get lineages for all taxids
                taxids = staxids.split(";")
                lineages = [tax.get_lineage_string(int(taxid)) for taxid in taxids]
                full_lineage = "; ".join(lineages)

                # Create new row with additional information
                new_row = [request.task, request.sample_name, *row, full_lineage]
                writer.writerow(new_row)

        logger.info(f"Annotated results written to {request.output_file}")

    except Exception:
        logger.exception("An unexpected error occurred")
        sys.exit(1)


def main() -> None:
    args = parse_command_line_args()
    annotate_blast_results(
        AnnotationRequest(
            input_file=args.input_file,
            output_file=args.output_file,
            sample_name=args.sample_name,
            task=args.task,
            taxonomy_dir=args.taxonomy_dir,
            sync=args.sync,
            taxonomy_mode=args.taxonomy_mode,
            taxonomy_max_age_days=args.taxonomy_max_age_days,
        ),
    )


if __name__ == "__main__":
    main()
