#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "snakemake",
# ]
# ///

import argparse
import csv
import logging
import os
import sys
from typing import Literal

from py_nvd import taxonomy
from snakemake.logging import logger as snakemake_logger

# Set up logging
logger = logging.getLogger(__name__)

# snakemake setup
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"


def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate BLAST results with taxonomy info"
    )
    # Keep --sqlite_cache for backward compatibility with Nextflow,
    # but it's no longer used (taxonomy.open() handles DB location)
    parser.add_argument(
        "--sqlite_cache",
        required=False,
        default=None,
        help="(Deprecated) Path to SQLite taxonomy DB - now handled automatically",
    )
    parser.add_argument("--input_file", required=True, help="Path to BLAST hits file")
    parser.add_argument("--output_file", required=True, help="Output TSV file")
    parser.add_argument("--sample_name", required=True, help="Sample name")
    parser.add_argument("--task", required=True, help="Task name")
    parser.add_argument(
        "--state-dir",
        required=False,
        default=None,
        help="State directory containing taxonomy cache (default: NVD_STATE_DIR or ~/.nvd/)",
    )

    return parser.parse_args()


def main() -> None:
    state_dir = None
    if "snakemake" in globals() and MODE == "snakemake":
        # sqlite_cache from snakemake is ignored - taxonomy.open() handles it
        input_file = snakemake.input.blast_results
        output_file = snakemake.output[0]
        sample_name = snakemake.params.sample
        task = snakemake.params.task
        # snakemake mode doesn't support state_dir yet
    else:
        args = parse_command_line_args()
        # args.sqlite_cache is ignored - taxonomy.open() handles DB location
        input_file = args.input_file
        output_file = args.output_file
        sample_name = args.sample_name
        task = args.task
        state_dir = args.state_dir

    try:
        if not os.path.exists(input_file):
            logger.error(f"Input file not found: {input_file}")
            sys.exit(1)

        if os.path.getsize(input_file) == 0:
            logger.warning("Input file is empty. Creating an empty output file.")
            open(output_file, "w").close()
            return

        with (
            taxonomy.open(state_dir=state_dir) as tax,
            open(input_file) as infile,
            open(output_file, "w", newline="") as outfile,
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
                new_row = [task, sample_name, *row, full_lineage]
                writer.writerow(new_row)

        logger.info(f"Annotated results written to {output_file}")

    except Exception:
        logger.exception("An unexpected error occurred")
        sys.exit(1)


if __name__ == "__main__":
    # Use Snakemake's logger
    logger = snakemake_logger
    main()
