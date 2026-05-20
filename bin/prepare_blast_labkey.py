#!/usr/bin/env python3
"""
Prepare BLAST data for LabKey upload.

The input TSV is already enriched with mapped_reads, total_reads,
blast_db_version, stat_kmer_db, and nextflow_run_id by ADD_READ_COUNTS_TO_BLAST.
This script adds the experiment_id column, renames columns to match
the LabKey schema, and converts TSV → CSV.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Column renames from BLAST TSV names → LabKey schema names
COLUMN_RENAMES = {
    "task": "blast_task",
    "sample": "sample_id",
    "rank": "tax_rank",
}


def rename_columns(row: dict[str, str]) -> dict[str, str]:
    """Rename TSV columns to match LabKey schema."""
    return {COLUMN_RENAMES.get(key, key): value for key, value in row.items()}


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Prepare BLAST data for LabKey upload")
    parser.add_argument(
        "--blast-csv",
        type=Path,
        required=True,
        help="Path to enriched BLAST TSV file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path to output CSV file",
    )
    parser.add_argument("--meta", type=str, required=True, help="Sample metadata/ID")
    parser.add_argument(
        "--experiment-id",
        type=int,
        required=True,
        help="Experiment ID",
    )
    args = parser.parse_args()

    try:
        records = []
        row_count = 0
        skipped_count = 0

        with open(args.blast_csv) as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                row_count += 1

                # Skip header rows accidentally included as data
                if row.get("qseqid", "").lower() == "qseqid":
                    skipped_count += 1
                    continue

                # Rename columns to match LabKey schema, then add experiment_id
                record = {"experiment": args.experiment_id}
                record.update(rename_columns(row))
                records.append(record)

        if records:
            fieldnames = list(records[0].keys())
            with open(args.output, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(records)
            logger.info(f"Wrote {len(records)} records to {args.output}")
        else:
            args.output.touch()
            logger.warning(f"No valid records. Created empty file: {args.output}")

        logger.info(
            f"Processed {row_count} rows: {len(records)} valid, {skipped_count} skipped",
        )

    except Exception as e:
        logger.error(f"Processing failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
