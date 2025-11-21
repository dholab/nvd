#!/usr/bin/env python3
"""
Prepare BLAST data for LabKey upload with data validation and type checking.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field, asdict
from enum import Enum

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class BlastTask(Enum):
    """Valid BLAST task types."""

    BLAST = "blast"
    MEGABLAST = "megablast"


class AdjustmentMethod(Enum):
    """Valid adjustment methods."""

    DOMINANT = "dominant"
    LCA = "lca"


@dataclass
class BlastRecord:
    """Validated BLAST record for LabKey upload."""

    experiment: int
    blast_task: str
    sample_id: str
    qseqid: str
    qlen: str
    sseqid: str
    stitle: str
    tax_rank: str
    length: int
    pident: float
    evalue: float
    bitscore: float
    sscinames: str
    staxids: int
    blast_db_version: str
    snakemake_run_id: str
    mapped_reads: int
    total_reads: int
    stat_db_version: str
    adjusted_taxid: Optional[int] = None
    adjustment_method: Optional[str] = None
    adjusted_taxid_name: Optional[str] = None
    adjusted_taxid_rank: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for CSV writing."""
        data = asdict(self)
        # Convert None values to empty strings for CSV
        for key, value in data.items():
            if value is None:
                data[key] = ""
        return data


class DataValidator:
    """Validates and converts BLAST data fields."""

    @staticmethod
    def validate_int(
        value: Any, field_name: str, allow_empty: bool = False
    ) -> Optional[int]:
        """Validate and convert to integer."""
        if value is None or value == "":
            if allow_empty:
                return None
            raise ValueError(f"{field_name} cannot be empty")

        try:
            return int(value)
        except (ValueError, TypeError) as e:
            raise ValueError(f"{field_name} must be an integer, got '{value}'") from e

    @staticmethod
    def validate_float(
        value: Any, field_name: str, allow_empty: bool = False
    ) -> Optional[float]:
        """Validate and convert to float."""
        if value is None or value == "":
            if allow_empty:
                return None
            raise ValueError(f"{field_name} cannot be empty")

        try:
            return float(value)
        except (ValueError, TypeError) as e:
            raise ValueError(f"{field_name} must be a float, got '{value}'") from e

    @staticmethod
    def validate_blast_task(value: str) -> str:
        """Validate BLAST task type."""
        value_lower = value.lower()
        try:
            task = BlastTask(value_lower)
            return task.value
        except ValueError:
            raise ValueError(
                f"blast_task must be 'blast' or 'megablast', got '{value}'"
            )

    @staticmethod
    def validate_adjustment_method(value: str) -> Optional[str]:
        """Validate adjustment method."""
        if value is None or value == "":
            return None

        value_lower = value.lower()
        try:
            method = AdjustmentMethod(value_lower)
            return method.value
        except ValueError:
            raise ValueError(
                f"adjustment_method must be 'dominant' or 'lca', got '{value}'"
            )

    @staticmethod
    def validate_string(value: Any, field_name: str, allow_empty: bool = True) -> str:
        """Validate string field."""
        if value is None:
            value = ""

        value_str = str(value).strip()

        if not allow_empty and not value_str:
            raise ValueError(f"{field_name} cannot be empty")

        return value_str


def read_contig_counts(filepath: Path) -> Dict[str, int]:
    """
    Read contig mapped read counts from TSV file.

    Args:
        filepath: Path to contig counts file

    Returns:
        Dictionary mapping contig ID to read count
    """
    contig_counts = {}

    try:
        with open(filepath, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) >= 2:
                    contig_id = parts[0]
                    try:
                        count = int(parts[1])
                        contig_counts[contig_id] = count
                    except ValueError:
                        logger.warning(
                            f"Line {line_num}: Invalid count value '{parts[1]}' for contig '{contig_id}'"
                        )
                else:
                    logger.warning(
                        f"Line {line_num}: Invalid format in contig counts file"
                    )

    except FileNotFoundError:
        logger.error(f"Contig counts file not found: {filepath}")
        raise
    except Exception as e:
        logger.error(f"Error reading contig counts file: {e}")
        raise

    logger.info(f"Loaded {len(contig_counts)} contig read counts")
    return contig_counts


def validate_and_convert_blast_row(
    row: Dict[str, str],
    contig_counts: Dict[str, int],
    experiment_id: int,
    run_id: str,
    blast_db_version: str,
    stat_db_version: str,
    total_reads: int,
    sample_id: str,
    validator: DataValidator,
) -> Optional[BlastRecord]:
    """
    Validate and convert a BLAST result row.

    Returns:
        BlastRecord if valid, None if row should be skipped
    """
    try:
        # Skip header rows that might be accidentally included
        if row.get("qseqid", "").lower() == "qseqid":
            logger.debug("Skipping header row")
            return None

        # Validate blast_task - skip if invalid
        try:
            blast_task = validator.validate_blast_task(row.get("task", "megablast"))
        except ValueError as e:
            logger.warning(f"Skipping row due to invalid blast_task: {e}")
            return None

        # Get mapped reads count
        qseqid = row.get("qseqid", "")
        mapped_reads = contig_counts.get(qseqid, 0)

        # Create validated record
        record = BlastRecord(
            experiment=experiment_id,
            blast_task=blast_task,
            sample_id=sample_id
            or validator.validate_string(
                row.get("sample", ""), "sample_id", allow_empty=False
            ),
            qseqid=validator.validate_string(qseqid, "qseqid", allow_empty=False),
            qlen=validator.validate_string(row.get("qlen", ""), "qlen"),
            sseqid=validator.validate_string(row.get("sseqid", ""), "sseqid"),
            stitle=validator.validate_string(row.get("stitle", ""), "stitle"),
            tax_rank=validator.validate_string(row.get("rank", ""), "tax_rank"),
            length=validator.validate_int(row.get("length", ""), "length"),
            pident=validator.validate_float(row.get("pident", ""), "pident"),
            evalue=validator.validate_float(row.get("evalue", ""), "evalue"),
            bitscore=validator.validate_float(row.get("bitscore", ""), "bitscore"),
            sscinames=validator.validate_string(row.get("sscinames", ""), "sscinames"),
            staxids=validator.validate_int(row.get("staxids", ""), "staxids"),
            blast_db_version=blast_db_version,
            snakemake_run_id=run_id,
            mapped_reads=mapped_reads,
            total_reads=total_reads,
            stat_db_version=stat_db_version,
            adjusted_taxid=validator.validate_int(
                row.get("adjusted_taxid", ""), "adjusted_taxid", allow_empty=True
            ),
            adjustment_method=validator.validate_adjustment_method(
                row.get("adjustment_method", "")
            ),
            adjusted_taxid_name=validator.validate_string(
                row.get("adjusted_taxid_name", ""), "adjusted_taxid_name"
            ),
            adjusted_taxid_rank=validator.validate_string(
                row.get("adjusted_taxid_rank", ""), "adjusted_taxid_rank"
            ),
        )

        return record

    except ValueError as e:
        logger.error(f"Validation error for row: {e}")
        logger.debug(f"Problematic row data: {row}")
        raise


def process_blast_data(
    blast_csv: Path,
    contig_counts: Dict[str, int],
    experiment_id: int,
    run_id: str,
    blast_db_version: str,
    stat_db_version: str,
    total_reads: int,
    meta: str,
) -> List[BlastRecord]:
    """
    Process BLAST CSV file and return validated records.
    """
    blast_records = []
    validator = DataValidator()
    sample_ids = set()
    row_count = 0
    skipped_count = 0
    error_count = 0

    try:
        with open(blast_csv, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row_num, row in enumerate(reader, 1):
                row_count += 1

                try:
                    record = validate_and_convert_blast_row(
                        row=row,
                        contig_counts=contig_counts,
                        experiment_id=experiment_id,
                        run_id=run_id,
                        blast_db_version=blast_db_version,
                        stat_db_version=stat_db_version,
                        total_reads=total_reads,
                        sample_id=meta,
                        validator=validator,
                    )

                    if record is None:
                        skipped_count += 1
                        continue

                    blast_records.append(record)
                    sample_ids.add(record.sample_id)

                except Exception as e:
                    error_count += 1
                    logger.error(f"Row {row_num}: Failed to process - {e}")
                    # Continue processing other rows instead of failing completely
                    continue

    except FileNotFoundError:
        logger.error(f"BLAST CSV file not found: {blast_csv}")
        raise
    except Exception as e:
        logger.error(f"Error reading BLAST CSV file: {e}")
        raise

    # Validate that all records have the same sample ID
    if len(sample_ids) > 1:
        logger.warning(
            f"Multiple sample IDs found: {sample_ids}. Expected only one sample per upload."
        )
        # You might want to make this an error instead:
        # raise ValueError(f"Multiple sample IDs found: {sample_ids}. Expected only one sample per upload.")

    logger.info(
        f"Processed {row_count} rows: {len(blast_records)} valid, {skipped_count} skipped, {error_count} errors"
    )

    return blast_records


def write_output(blast_records: List[BlastRecord], output_path: Path) -> None:
    """
    Write validated BLAST records to CSV file.
    """
    if blast_records:
        with open(output_path, "w", newline="") as f:
            # Get fieldnames from the first record
            fieldnames = list(blast_records[0].to_dict().keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for record in blast_records:
                writer.writerow(record.to_dict())

        logger.info(f"Wrote {len(blast_records)} records to {output_path}")
    else:
        # Create empty file if no data
        output_path.touch()
        logger.warning(f"No valid records to write. Created empty file: {output_path}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Prepare BLAST data for LabKey upload")

    parser.add_argument(
        "--blast-csv", type=Path, required=True, help="Path to BLAST results CSV file"
    )
    parser.add_argument(
        "--contig-counts",
        type=Path,
        required=True,
        help="Path to contig mapped read counts file",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Path to output CSV file"
    )
    parser.add_argument("--meta", type=str, required=True, help="Sample metadata/ID")
    parser.add_argument(
        "--experiment-id", type=int, required=True, help="Experiment ID"
    )
    parser.add_argument("--run-id", type=str, required=True, help="Snakemake run ID")
    parser.add_argument(
        "--total-reads", type=int, required=True, help="Total number of reads"
    )
    parser.add_argument(
        "--blast-db-version", type=str, required=True, help="BLAST database version"
    )
    parser.add_argument(
        "--stat-db-version", type=str, required=True, help="Statistics database version"
    )

    args = parser.parse_args()

    try:
        # Read contig counts
        contig_counts = read_contig_counts(args.contig_counts)

        # Process BLAST data
        blast_records = process_blast_data(
            blast_csv=args.blast_csv,
            contig_counts=contig_counts,
            experiment_id=args.experiment_id,
            run_id=args.run_id,
            blast_db_version=args.blast_db_version,
            stat_db_version=args.stat_db_version,
            total_reads=args.total_reads,
            meta=args.meta,
        )

        # Write output
        write_output(blast_records, args.output)

        logger.info("Processing completed successfully")

    except Exception as e:
        logger.error(f"Processing failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
