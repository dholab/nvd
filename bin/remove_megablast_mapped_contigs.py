#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "snakemake",
# ]
# ///

import argparse
import logging
import subprocess
from contextlib import redirect_stderr
from types import SimpleNamespace
from typing import Literal

# snakemake setup
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"


def setup_logger(log_file: str | None) -> logging.Logger:
    logger = logging.getLogger("remove_megablast_mapped_contigs")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    if not log_file:
        return logger

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger


def extract_unique_ids(input_file, output_file, logger):
    unique_ids = set()
    logger.info(f"Reading megablast results from {input_file}")
    with open(input_file) as f:
        for line in f:
            unique_ids.add(line.rstrip("\n").split("\t")[0] + "\n")

    logger.info(f"Writing classified contigs to {output_file}")
    with open(output_file, "w") as f:
        for contig_id in unique_ids:
            f.write(contig_id)
    logger.info(f"Wrote {len(unique_ids)} unique contig IDs")


def run_filterbyname(input_fasta, output_fasta, names_file, logger):
    command = [
        "seqkit",
        "grep",
        "-n",
        "-f",
        f"{names_file}",
        "-o",
        f"{output_fasta}",
        f"{input_fasta}",
    ]
    logger.info(f"Running command: {' '.join(command)}")
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        logger.info(result.stdout)
        if result.stderr:
            logger.warning(result.stderr)
    except subprocess.CalledProcessError:
        logger.exception("Error running filterbyname.sh")
        raise


def parse_args() -> SimpleNamespace:
    if MODE == "snakemake":
        return SimpleNamespace(input=snakemake.input, output=snakemake.output, log=snakemake.log)
    parser = argparse.ArgumentParser(description="Remove contigs mapped by megablast")
    parser.add_argument("--megablast_results", required=True, help="Path to megablast results")
    parser.add_argument("--contigs_fasta", required=True, help="Path to original contigs FASTA")
    parser.add_argument(
        "--classified_contigs",
        required=True,
        help="Path to output classified contigs list",
    )
    parser.add_argument("--pruned_contigs", required=True, help="Path to output pruned FASTA")
    parser.add_argument("--log", default=None, required=False, help="Path to log file")
    args = parser.parse_args()

    return SimpleNamespace(
        input=SimpleNamespace(
            megablast_results=args.megablast_results,
            contigs_fasta=args.contigs_fasta,
        ),
        output=SimpleNamespace(
            classified_contigs=args.classified_contigs,
            pruned_contigs=args.pruned_contigs,
        ),
        log=args.log,
    )


def main() -> None:
    args = parse_args()

    if args.log and len(args.log > 1):
        with open(args.log[0], "w") as log_file, redirect_stderr(log_file):
            logger = setup_logger(args.log[0])
            logger.info("Starting remove_megablast_mapped_contigs script")
    else:
        logger = setup_logger(args.log)
        logger.info("Starting remove_megablast_mapped_contigs script")

    try:
        extract_unique_ids(
            args.input.megablast_results,
            args.output.classified_contigs,
            logger,
        )

        run_filterbyname(
            args.input.contigs_fasta,
            args.output.pruned_contigs,
            args.output.classified_contigs,
            logger,
        )

        logger.info("Script completed successfully")
    except Exception:
        logger.exception("An error occurred.")
        raise


if __name__ == "__main__":
    main()
