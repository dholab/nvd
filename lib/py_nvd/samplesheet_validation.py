"""Reusable samplesheet validation for NVD CLI and runtime preflight paths."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from py_nvd.read_inputs import ResolutionError, read_rows, resolve_rows

if TYPE_CHECKING:
    from pathlib import Path

REQUIRED_COLUMNS = ("sample_id", "srr", "platform", "fastq1", "fastq2")
GLOB_COLUMNS = ("fastq1_glob", "fastq2_glob")
OUTPUT_COLUMNS = REQUIRED_COLUMNS + GLOB_COLUMNS


@dataclass
class SamplesheetValidationResult:
    """
    Result of samplesheet validation.

    Attributes:
        valid: True if samplesheet passed all validation checks
        samples: List of valid sample IDs found
        errors: List of validation error messages
        warnings: List of warning messages (non-fatal)
        header_valid: True if header contains all required columns
    """

    valid: bool
    samples: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    header_valid: bool = False
    header_columns: list[str] = field(default_factory=list)


def validate_samplesheet(path: Path) -> SamplesheetValidationResult:
    """
    Validate a samplesheet CSV file.

    Checks:
    - File exists and is readable
    - Header contains all required columns
    - Each non-comment row can be resolved by the same read-input preflight used
      when NVD starts a run
    """
    result = SamplesheetValidationResult(valid=False)

    if not path.exists():
        result.errors.append(f"Samplesheet not found: {path}")
        return result

    try:
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)

            if not reader.fieldnames:
                result.errors.append("Samplesheet is empty or has no header")
                return result

            result.header_columns = list(reader.fieldnames)

            fieldnames: set[str] = set(reader.fieldnames)
            missing_cols = set(REQUIRED_COLUMNS) - fieldnames
            if missing_cols:
                result.errors.append(
                    f"Missing required columns: {', '.join(sorted(missing_cols))}",
                )
                return result

            result.header_valid = True

            resolution = resolve_rows(read_rows(path))
            result.samples = [sample.sample_id for sample in resolution.samples]
            result.warnings.extend(resolution.warnings)
            result.valid = True

    except csv.Error as e:
        result.errors.append(f"CSV parsing error: {e}")
    except OSError as e:
        result.errors.append(f"Failed to read samplesheet: {e}")
    except ResolutionError as e:
        result.errors.append(str(e))

    return result
