#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

import argparse
import csv
import logging
import subprocess


def setup_logger(log_file: str | None) -> logging.Logger:
    logger = logging.getLogger("partition_megablast_queries")
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


def write_accounted_query_ids(
    megablast_results: str,
    accounted_query_ids: str,
    logger: logging.Logger,
) -> int:
    unique_ids = set()
    logger.info(f"Reading MEGABLAST results from {megablast_results}")
    with open(megablast_results) as f:
        for line in f:
            unique_ids.add(line.rstrip("\n").split("\t")[0] + "\n")

    logger.info(f"Writing accounted query IDs to {accounted_query_ids}")
    with open(accounted_query_ids, "w") as f:
        f.writelines(sorted(unique_ids))
    logger.info(f"Wrote {len(unique_ids)} unique query IDs")
    return len(unique_ids)


def write_blastn_candidate_fasta(
    query_fasta: str,
    blastn_candidate_fasta: str,
    accounted_query_ids: str,
    logger: logging.Logger,
) -> None:
    command = [
        "seqkit",
        "grep",
        "-v",
        "-f",
        f"{accounted_query_ids}",
        "-o",
        f"{blastn_candidate_fasta}",
        f"{query_fasta}",
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


def count_fasta_records(path: str) -> int:
    with open(path) as handle:
        return sum(line.startswith(">") for line in handle)


def write_partition_summary(
    *,
    output: str,
    sample_id: str,
    query_class: str,
    query_count: int,
    accounted_count: int,
    candidate_count: int,
) -> None:
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "sample_id",
                "query_class",
                "queries_in",
                "megablast_accounted",
                "blastn_candidates",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": sample_id,
                "query_class": query_class,
                "queries_in": query_count,
                "megablast_accounted": accounted_count,
                "blastn_candidates": candidate_count,
            },
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Partition MEGABLAST queries into accounted IDs and BLASTN candidates",
    )
    parser.add_argument(
        "--megablast-results",
        required=True,
        help="Path to raw MEGABLAST results",
    )
    parser.add_argument(
        "--query-fasta",
        required=True,
        help="Path to the original query FASTA",
    )
    parser.add_argument(
        "--accounted-query-ids",
        required=True,
        help="Path to query IDs accounted for by MEGABLAST",
    )
    parser.add_argument(
        "--blastn-candidate-fasta",
        required=True,
        help="Path to residual queries requiring BLASTN",
    )
    parser.add_argument("--sample-id")
    parser.add_argument("--query-class")
    parser.add_argument("--summary-tsv")
    args = parser.parse_args()
    summary_args = (args.sample_id, args.query_class, args.summary_tsv)
    if any(summary_args) and not all(summary_args):
        parser.error(
            "--sample-id, --query-class, and --summary-tsv must be supplied together",
        )
    return args


def main() -> None:
    args = parse_args()

    logger = setup_logger(None)
    logger.info("Starting partition_megablast_queries script")

    try:
        accounted_count = write_accounted_query_ids(
            args.megablast_results,
            args.accounted_query_ids,
            logger,
        )

        write_blastn_candidate_fasta(
            args.query_fasta,
            args.blastn_candidate_fasta,
            args.accounted_query_ids,
            logger,
        )
        if args.summary_tsv is not None:
            write_partition_summary(
                output=args.summary_tsv,
                sample_id=args.sample_id,
                query_class=args.query_class,
                query_count=count_fasta_records(args.query_fasta),
                accounted_count=accounted_count,
                candidate_count=count_fasta_records(args.blastn_candidate_fasta),
            )

        logger.info("Script completed successfully")
    except Exception:
        logger.exception("An error occurred.")
        raise


if __name__ == "__main__":
    main()
