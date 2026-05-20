#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

import argparse
import logging
import subprocess


def setup_logger(log_file: str | None) -> logging.Logger:
    logger = logging.getLogger("remove_megablast_mapped_contigs")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    if not log_file:
        return logger

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger


def extract_unique_ids(input_file: str, output_file: str, logger: logging.Logger) -> None:
    unique_ids = set()
    logger.info(f"Reading megablast results from {input_file}")
    with open(input_file) as f:
        for line in f:
            unique_ids.add(line.rstrip("\n").split("\t")[0] + "\n")

    logger.info(f"Writing classified contigs to {output_file}")
    with open(output_file, "w") as f:
        f.writelines(unique_ids)
    logger.info(f"Wrote {len(unique_ids)} unique contig IDs")


def run_filterbyname(
    input_fasta: str,
    output_fasta: str,
    names_file: str,
    logger: logging.Logger,
) -> None:
    command = [
        "seqkit",
        "grep",
        "-v",
        "-n",
        "-f",
        f"{names_file}",
        "-o",
        f"{output_fasta}",
        f"{input_fasta}",
    ]
    logger.info(f"Running command: {' '.join(command)}")
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)  # noqa: S603
        logger.info(result.stdout)
        if result.stderr:
            logger.warning(result.stderr)
    except subprocess.CalledProcessError:
        logger.exception(f"Error running {command}")
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Remove contigs mapped by megablast")
    parser.add_argument(
        "--megablast_results",
        required=True,
        help="Path to megablast results",
    )
    parser.add_argument(
        "--contigs_fasta",
        required=True,
        help="Path to original contigs FASTA",
    )
    parser.add_argument(
        "--classified_contigs",
        required=True,
        help="Path to output classified contigs list",
    )
    parser.add_argument(
        "--pruned_contigs",
        required=True,
        help="Path to output pruned FASTA",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    logger = setup_logger(None)
    logger.info("Starting remove_megablast_mapped_contigs script")

    try:
        extract_unique_ids(
            args.megablast_results,
            args.classified_contigs,
            logger,
        )

        run_filterbyname(
            args.contigs_fasta,
            args.pruned_contigs,
            args.classified_contigs,
            logger,
        )

        logger.info("Script completed successfully")
    except Exception:
        logger.exception("An error occurred.")
        raise


if __name__ == "__main__":
    main()
