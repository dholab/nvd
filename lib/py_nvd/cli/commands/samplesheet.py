# ruff: noqa: B008
"""
Samplesheet commands for the NVD CLI.

Commands:
    nvd samplesheet generate  - Generate samplesheet from FASTQ directory or SRA list
    nvd samplesheet validate  - Validate samplesheet format and content

The generate command scans a directory for FASTQ files, pairs them by sample,
and creates a properly formatted samplesheet CSV. It validates its own output
using the same logic as the validate command.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import typer
from rich.table import Table

from py_nvd.cli.prompts import (
    PROMPT_TIMEOUT_SECONDS,
    PromptTimeout,
    confirm_with_timeout,
    is_interactive,
    prompt_with_timeout,
)
from py_nvd.cli.utils import console, error, info, success, warning
from py_nvd.state import get_sample_history, was_sample_ever_uploaded

REQUIRED_COLUMNS = ("sample_id", "srr", "platform", "fastq1", "fastq2")
VALID_PLATFORMS = {"illumina", "ont", ""}

# Expected number of FASTQ files for paired-end samples
PAIRED_END_FILE_COUNT = 2

# FASTQ file extensions to recognize
FASTQ_EXTENSIONS = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}

# Regex patterns for extracting sample stem and read number
# Matches: _R1, _R2, _1, _2, .R1, .R2, .1, .2 (optionally followed by _001 etc)
READ_INDICATOR_PATTERN = re.compile(
    r"[._](?:R)?([12])(?:_\d+)?(?:\.fastq|\.fq)?(?:\.gz)?$",
    re.IGNORECASE,
)

# Pattern to extract the sample stem (everything before the read indicator)
STEM_PATTERN = re.compile(
    r"^(.+?)[._](?:R)?[12](?:_\d+)?(?:\.fastq|\.fq)?(?:\.gz)?$",
    re.IGNORECASE,
)


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


def validate_samplesheet(path: Path) -> SamplesheetValidationResult:
    """
    Validate a samplesheet CSV file.

    Checks:
    - File exists and is readable
    - Header contains all required columns
    - Each row has a non-empty sample_id (unless commented with #)
    - Platform is valid (illumina, ont, or empty)
    - Either srr or fastq1 is provided
    - Warns on duplicate sample IDs

    Args:
        path: Path to the samplesheet CSV file

    Returns:
        SamplesheetValidationResult with validation status and details
    """
    result = SamplesheetValidationResult(valid=False)

    if not path.exists():
        result.errors.append(f"Samplesheet not found: {path}")
        return result

    try:
        with path.open() as f:
            reader = csv.DictReader(f)

            if not reader.fieldnames:
                result.errors.append("Samplesheet is empty or has no header")
                return result

            # Check for required columns
            fieldnames: set[str] = set(reader.fieldnames)
            required: set[str] = set(REQUIRED_COLUMNS)
            missing_cols = required - fieldnames
            if missing_cols:
                result.errors.append(
                    f"Missing required columns: {', '.join(sorted(missing_cols))}"
                )
                return result

            result.header_valid = True

            # Validate rows
            samples: list[str] = []
            for i, row in enumerate(reader, start=2):  # Line 1 is header
                sample_id = row.get("sample_id", "").strip()

                # Skip comment lines
                if sample_id.startswith("#"):
                    continue

                # Check sample_id not empty
                if not sample_id:
                    result.errors.append(f"Line {i}: sample_id is empty")
                    continue

                # Check platform is valid
                platform = row.get("platform", "").strip().lower()
                if platform and platform not in VALID_PLATFORMS:
                    result.errors.append(
                        f"Line {i} ({sample_id}): invalid platform '{platform}'"
                    )

                # Check that either SRR or FASTQ files are provided
                srr = row.get("srr", "").strip()
                fastq1 = row.get("fastq1", "").strip()

                if not srr and not fastq1:
                    result.errors.append(
                        f"Line {i} ({sample_id}): must provide either 'srr' or 'fastq1'"
                    )

                samples.append(sample_id)

            result.samples = samples

            # Check for duplicate sample IDs
            duplicates = [s for s, count in Counter(samples).items() if count > 1]
            if duplicates:
                result.warnings.append(
                    f"Duplicate sample_id found: {', '.join(duplicates)}"
                )

            # Valid if no errors
            result.valid = len(result.errors) == 0

    except csv.Error as e:
        result.errors.append(f"CSV parsing error: {e}")
    except OSError as e:
        result.errors.append(f"Failed to read samplesheet: {e}")

    return result


@dataclass
class SampleEntry:
    """
    Represents a sample entry for the samplesheet.

    Attributes:
        sample_id: Unique sample identifier (derived from filename stem)
        srr: SRA accession (empty for local files)
        platform: Sequencing platform (illumina or ont)
        fastq1: Path to read 1 FASTQ file
        fastq2: Path to read 2 FASTQ file (empty for single-end)
    """

    sample_id: str
    srr: str = ""
    platform: str = ""
    fastq1: str = ""
    fastq2: str = ""


def _get_fastq_files(directory: Path) -> list[Path]:
    """
    Find all FASTQ files in a directory.

    Searches for files with extensions: .fastq, .fq, .fastq.gz, .fq.gz

    Args:
        directory: Directory to search

    Returns:
        List of Path objects for FASTQ files, sorted alphabetically
    """
    fastq_files: list[Path] = []

    for ext in FASTQ_EXTENSIONS:
        # Handle double extensions like .fastq.gz
        if ext.startswith(".fastq") or ext.startswith(".fq"):
            fastq_files.extend(directory.glob(f"*{ext}"))

    # Deduplicate (glob patterns may overlap) and sort
    return sorted(set(fastq_files))


def _extract_stem(filename: str) -> str | None:
    """
    Extract the sample stem from a FASTQ filename.

    Removes read indicators (_R1, _R2, _1, _2, etc.) and extensions.

    Args:
        filename: FASTQ filename (not full path)

    Returns:
        Sample stem, or None if pattern doesn't match

    Examples:
        >>> _extract_stem("sample1_R1.fastq.gz")
        'sample1'
        >>> _extract_stem("sample1_R2_001.fastq.gz")
        'sample1'
        >>> _extract_stem("sample1.1.fq.gz")
        'sample1'
    """
    match = STEM_PATTERN.match(filename)
    if match:
        return match.group(1)

    # Fallback: strip known extensions for single-end files
    name = filename
    for ext in sorted(FASTQ_EXTENSIONS, key=len, reverse=True):
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break

    return name if name else None


def _get_read_number(filename: str) -> int | None:
    """
    Extract the read number (1 or 2) from a FASTQ filename.

    Args:
        filename: FASTQ filename

    Returns:
        1 or 2 if found, None otherwise
    """
    match = READ_INDICATOR_PATTERN.search(filename)
    if match:
        return int(match.group(1))
    return None


def scan_fastq_directory(directory: Path) -> tuple[list[SampleEntry], list[str]]:
    """
    Scan a directory for FASTQ files and pair them by sample.

    Groups files by sample stem, then pairs R1/R2 files. Files without
    read indicators are treated as single-end.

    Args:
        directory: Directory containing FASTQ files

    Returns:
        Tuple of (list of SampleEntry, list of warning messages)
    """
    fastq_files = _get_fastq_files(directory)
    warnings_list: list[str] = []

    if not fastq_files:
        return [], [f"No FASTQ files found in {directory}"]

    # Group files by stem
    stem_to_files: dict[str, list[Path]] = {}
    for fq in fastq_files:
        stem = _extract_stem(fq.name)
        if stem is None:
            warnings_list.append(f"Could not parse filename: {fq.name}")
            continue

        if stem not in stem_to_files:
            stem_to_files[stem] = []
        stem_to_files[stem].append(fq)

    # Build sample entries
    samples: list[SampleEntry] = []

    for stem, files in sorted(stem_to_files.items()):
        if len(files) == 1:
            # Single-end
            samples.append(
                SampleEntry(
                    sample_id=stem,
                    fastq1=str(files[0]),
                )
            )
        elif len(files) == PAIRED_END_FILE_COUNT:
            # Paired-end: determine which is R1 and R2
            r1, r2 = None, None
            for f in files:
                read_num = _get_read_number(f.name)
                if read_num == 1:
                    r1 = f
                elif read_num == 2:
                    r2 = f

            if r1 and r2:
                samples.append(
                    SampleEntry(
                        sample_id=stem,
                        fastq1=str(r1),
                        fastq2=str(r2),
                    )
                )
            else:
                # Couldn't determine pairing, use alphabetical order
                warnings_list.append(
                    f"Could not determine R1/R2 for {stem}, using alphabetical order"
                )
                sorted_files = sorted(files, key=lambda p: p.name)
                samples.append(
                    SampleEntry(
                        sample_id=stem,
                        fastq1=str(sorted_files[0]),
                        fastq2=str(sorted_files[1]),
                    )
                )
        else:
            # More than 2 files - ambiguous
            warnings_list.append(
                f"Sample {stem} has {len(files)} FASTQ files (expected 1-2), skipping"
            )

    return samples, warnings_list


def parse_sra_accessions(path: Path) -> tuple[list[SampleEntry], list[str]]:
    """
    Parse SRA accessions from a text file.

    Expects one accession per line. Lines starting with # are comments.
    Blank lines are skipped.

    Args:
        path: Path to text file with SRA accessions

    Returns:
        Tuple of (list of SampleEntry, list of warning messages)
    """
    warnings_list: list[str] = []
    samples: list[SampleEntry] = []

    try:
        with path.open() as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue

                # Basic SRA accession validation (SRR, ERR, DRR followed by digits)
                if not re.match(r"^[SED]RR\d+$", line, re.IGNORECASE):
                    warnings_list.append(
                        f"Line {line_num}: '{line}' doesn't look like an SRA accession"
                    )

                # Use accession as sample_id
                samples.append(
                    SampleEntry(
                        sample_id=line,
                        srr=line,
                    )
                )

    except OSError as e:
        warnings_list.append(f"Failed to read SRA file: {e}")

    return samples, warnings_list


def write_samplesheet(samples: list[SampleEntry], output: Path, platform: str) -> None:
    """
    Write samples to a samplesheet CSV file.

    Args:
        samples: List of SampleEntry objects
        output: Output file path
        platform: Platform to set for all samples
    """
    with output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(REQUIRED_COLUMNS)

        for sample in samples:
            writer.writerow(
                [
                    sample.sample_id,
                    sample.srr,
                    platform,
                    sample.fastq1,
                    sample.fastq2,
                ]
            )


def format_samples_table(samples: list[SampleEntry], platform: str) -> Table:
    """
    Create a Rich table for previewing samples.

    Args:
        samples: List of SampleEntry objects
        platform: Platform to display

    Returns:
        Rich Table object
    """
    table = Table(title="Samplesheet Preview", show_lines=False)
    table.add_column("Sample ID", style="cyan", no_wrap=True)
    table.add_column("SRR", style="dim")
    table.add_column("Platform", style="green")
    table.add_column("Read 1", style="white", overflow="ellipsis")
    table.add_column("Read 2", style="white", overflow="ellipsis")

    # Show up to 15 samples, then summarize
    max_display = 15
    for sample in samples[:max_display]:
        # Truncate long paths for display
        fastq1_display = Path(sample.fastq1).name if sample.fastq1 else "-"
        fastq2_display = Path(sample.fastq2).name if sample.fastq2 else "-"

        table.add_row(
            sample.sample_id,
            sample.srr or "-",
            platform,
            fastq1_display,
            fastq2_display,
        )

    if len(samples) > max_display:
        remaining = len(samples) - max_display
        table.add_row(
            f"[dim]... and {remaining} more[/dim]",
            "",
            "",
            "",
            "",
        )

    return table


@dataclass
class PreviousProcessingInfo:
    """Information about a sample's previous processing history."""

    sample_id: str
    run_ids: list[str]
    was_uploaded: bool


def check_samples_against_state(
    sample_ids: list[str],
    state_dir: Path | str | None = None,
) -> list[PreviousProcessingInfo]:
    """
    Check if any samples were previously processed (read-only state query).

    This is a lightweight check that queries the state database without
    modifying it. Useful for warning users about potential duplicates
    before they run the pipeline.

    Args:
        sample_ids: List of sample IDs to check
        state_dir: Optional state directory override

    Returns:
        List of PreviousProcessingInfo for samples that were previously processed
    """
    previously_processed: list[PreviousProcessingInfo] = []

    for sid in sample_ids:
        history = get_sample_history(sid, state_dir=state_dir)
        if history:
            # Get unique run IDs from history
            run_ids = list(dict.fromkeys(h.run_id for h in history))
            was_uploaded = was_sample_ever_uploaded(sid, state_dir=state_dir)
            previously_processed.append(
                PreviousProcessingInfo(
                    sample_id=sid,
                    run_ids=run_ids,
                    was_uploaded=was_uploaded,
                )
            )

    return previously_processed


def format_state_warnings(previously_processed: list[PreviousProcessingInfo]) -> None:
    """Display warnings about previously processed samples."""
    if not previously_processed:
        return

    console.print()
    warning(
        f"[bold]{len(previously_processed)} sample(s) were previously processed:[/bold]"
    )

    # Show details for each previously processed sample
    for prev in previously_processed[:10]:  # Show first 10
        run_str = ", ".join(prev.run_ids[:3])  # Show first 3 runs
        if len(prev.run_ids) > 3:
            run_str += f" (+{len(prev.run_ids) - 3} more)"

        upload_marker = " [dim](uploaded)[/dim]" if prev.was_uploaded else ""
        console.print(f"  • [cyan]{prev.sample_id}[/cyan] → {run_str}{upload_marker}")

    if len(previously_processed) > 10:
        console.print(f"  [dim]... and {len(previously_processed) - 10} more[/dim]")

    console.print()
    info(
        "These samples may be skipped or cause duplicate detection errors. "
        "Use --force to proceed anyway."
    )


samplesheet_app = typer.Typer(
    name="samplesheet",
    help="Samplesheet generation and validation",
    no_args_is_help=True,
)


@samplesheet_app.command("generate")
@samplesheet_app.command("gen", hidden=True)  # Alias
def generate(  # noqa: PLR0912, PLR0913, C901
    from_dir: Path | None = typer.Option(
        None,
        "--from-dir",
        "-d",
        help="Directory containing FASTQ files to scan",
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    from_sra: Path | None = typer.Option(
        None,
        "--from-sra",
        help="Text file with SRA accessions (one per line)",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        Path("samplesheet.csv"),
        "--output",
        "-o",
        help="Output samplesheet path",
    ),
    platform: str | None = typer.Option(
        None,
        "--platform",
        "-p",
        help="Sequencing platform (illumina or ont). Prompts if not provided.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        "-n",
        help="Preview samplesheet without writing file",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Overwrite output file without confirmation",
    ),
) -> None:
    """
    Generate a samplesheet CSV from FASTQ files or SRA accessions.

    Scans a directory for FASTQ files, pairs them by sample name, and creates
    a properly formatted samplesheet. The generated samplesheet is validated
    before writing.

    If a state database exists, automatically checks for previously processed
    samples and warns about potential duplicates.

    Examples:

        # Generate from FASTQ directory (will prompt for platform)
        nvd samplesheet generate --from-dir ./fastqs -o samples.csv

        # Generate with platform specified (no prompts)
        nvd samplesheet generate -d ./fastqs -o samples.csv --platform illumina

        # Generate from SRA accession list
        nvd samplesheet generate --from-sra accessions.txt -o samples.csv -p illumina

        # Preview without writing
        nvd samplesheet generate -d ./fastqs --dry-run
    """
    console.print()

    # Validate input source
    if from_dir is None and from_sra is None:
        error("Must specify either --from-dir or --from-sra")

    if from_dir is not None and from_sra is not None:
        error("Cannot specify both --from-dir and --from-sra")

    # Handle platform
    if platform is not None:
        # Validate provided platform
        platform = platform.lower()
        if platform not in ("illumina", "ont"):
            error(f"Invalid platform: {platform}. Must be 'illumina' or 'ont'.")
    else:
        # Need to prompt for platform
        if not is_interactive():
            error("--platform is required in non-interactive mode")

        try:
            platform = prompt_with_timeout(
                console,
                "Sequencing platform",
                choices=["illumina", "ont"],
                default="illumina",
                timeout_seconds=PROMPT_TIMEOUT_SECONDS,
            )
            console.print()
        except PromptTimeout:
            console.print()
            error(
                f"Prompt timed out after {PROMPT_TIMEOUT_SECONDS // 60} minutes. "
                "Use --platform to specify non-interactively."
            )
        except KeyboardInterrupt:
            console.print("\n[yellow]Cancelled[/yellow]")
            raise typer.Abort from None

    # At this point, platform is guaranteed to be a string
    # (either provided via CLI or obtained from prompt)
    assert platform is not None

    # Scan for samples
    if from_dir is not None:
        info(f"Scanning directory: {from_dir}")
        samples, scan_warnings = scan_fastq_directory(from_dir)
    else:
        assert from_sra is not None  # for type checker
        info(f"Reading SRA accessions: {from_sra}")
        samples, scan_warnings = parse_sra_accessions(from_sra)

    # Display warnings from scanning
    for warn in scan_warnings:
        warning(warn)

    if not samples:
        error("No samples found")

    # Check for previously processed samples (read-only state query)
    sample_ids = [s.sample_id for s in samples]
    try:
        previously_processed = check_samples_against_state(sample_ids)
        if previously_processed:
            format_state_warnings(previously_processed)
    except Exception:
        # State database may not exist yet - that's fine, just skip the check
        pass

    # Display preview
    console.print()
    table = format_samples_table(samples, platform)
    console.print(table)
    console.print()

    # Summary
    paired = sum(1 for s in samples if s.fastq2)
    single = sum(1 for s in samples if s.fastq1 and not s.fastq2)
    sra = sum(1 for s in samples if s.srr)

    info(f"Found {len(samples)} samples:")
    if paired:
        info(f"  • {paired} paired-end")
    if single:
        info(f"  • {single} single-end")
    if sra:
        info(f"  • {sra} SRA accessions")

    console.print()

    # Dry run stops here
    if dry_run:
        info("Dry run mode - no file written")
        return

    # Check output file
    if output.exists() and not force:
        warning(f"Output file already exists: {output}")
        try:
            if not confirm_with_timeout(
                console,
                "Overwrite?",
                default=False,
                timeout_seconds=PROMPT_TIMEOUT_SECONDS,
            ):
                info("Aborted")
                raise typer.Abort from None
        except PromptTimeout:
            console.print()
            error("Confirmation timed out. Use --force to overwrite.")
        except KeyboardInterrupt:
            console.print("\n[yellow]Cancelled[/yellow]")
            raise typer.Abort from None
        console.print()

    # Confirm before writing (unless force)
    if not force:
        try:
            if not confirm_with_timeout(
                console,
                f"Write samplesheet to [cyan]{output}[/cyan]?",
                default=True,
                timeout_seconds=PROMPT_TIMEOUT_SECONDS,
            ):
                info("Aborted")
                raise typer.Abort from None
        except PromptTimeout:
            console.print()
            error("Confirmation timed out. Use --force to skip confirmation.")
        except KeyboardInterrupt:
            console.print("\n[yellow]Cancelled[/yellow]")
            raise typer.Abort from None
        console.print()

    # Write samplesheet
    write_samplesheet(samples, output, platform)
    success(f"Wrote samplesheet to {output}")

    # Validate the output
    console.print()
    info("Validating generated samplesheet...")
    result = validate_samplesheet(output)

    if result.valid:
        success(f"Samplesheet is valid ({len(result.samples)} samples)")
    else:
        warning("Generated samplesheet has validation issues:")
        for err in result.errors:
            console.print(f"  [red]•[/red] {err}")
        # This shouldn't happen if our generation logic is correct
        sys.exit(1)

    for warn in result.warnings:
        warning(warn)


@samplesheet_app.command("validate")
@samplesheet_app.command("val", hidden=True)  # Alias
def validate_cmd(
    samplesheet: Path = typer.Argument(..., help="Path to samplesheet CSV"),
) -> None:
    """
    Validate samplesheet format and content.

    Checks that the samplesheet has the required columns, valid platforms,
    and that each row has either an SRR accession or FASTQ file paths.

    This is equivalent to 'nvd validate samplesheet'.
    """
    console.print(f"\n[bold]Validating Samplesheet[/bold]: {samplesheet}\n")

    result = validate_samplesheet(samplesheet)

    # Display header status
    if result.header_valid:
        success(f"Header valid: {', '.join(REQUIRED_COLUMNS)}")
    elif result.errors:
        # Header errors will be in the errors list
        pass

    # Display errors
    if result.errors:
        console.print("\n[red]Validation errors:[/red]")
        for err in result.errors:
            console.print(f"  • {err}")
        console.print()
        error(f"Found {len(result.errors)} validation error(s)")

    # Display warnings
    for warn in result.warnings:
        warning(warn)

    # Success
    success(f"Found {len(result.samples)} valid samples")
    console.print()
    success("Samplesheet is valid")
