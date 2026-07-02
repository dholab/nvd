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
import glob
import re
import sys
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import typer
from rich.table import Table

from py_nvd.cli.prompts import (
    PROMPT_TIMEOUT_SECONDS,
    PromptTimeoutError,
    confirm_with_timeout,
    is_interactive,
    prompt_with_timeout,
)
from py_nvd.cli.samplesheet_validation import render_samplesheet_validation_result
from py_nvd.cli.utils import console, error, info, success, warning
from py_nvd.read_filenames import (
    GENERATED_FASTQ_SUFFIXES,
    CasavaLaneName,
    parse_casava_lane_name,
)
from py_nvd.read_inputs import (
    ResolutionError,
    normalize_platform,
)
from py_nvd.samplesheet_validation import (
    OUTPUT_COLUMNS,
    REQUIRED_COLUMNS,
    validate_samplesheet,
)
from py_nvd.sra_accessions import looks_like_sra_accession

# Expected number of FASTQ files for paired-end samples
PAIRED_END_FILE_COUNT = 2

# Regex patterns for extracting sample stem and read number
# Matches: _R1, _R2, _1, _2, .R1, .R2, .1, .2 (optionally followed by _001 etc)
READ_INDICATOR_PATTERN = re.compile(
    r"[._](?:R)?([12])(?:_\d+)?(?:\.fastq|\.fq)?(?:\.gz)?$",
    re.IGNORECASE,
)

# Pattern to extract the sample stem (everything before the read indicator)
STEM_PATTERN = re.compile(
    r"^(.+)[._](?:R)?[12](?:_\d+)?(?:\.fastq|\.fq)?(?:\.gz)?$",
    re.IGNORECASE,
)

# Illumina/CASAVA sample/lane suffixes left after read indicators are removed.
# Example: sample_S1_L001_R1_001.fastq.gz -> sample_S1_L001 -> sample
CASAVA_SAMPLE_LANE_SUFFIX_PATTERN = re.compile(r"_S\d+(?:_L\d{3})?$", re.IGNORECASE)


class SamplesheetGenerationError(ValueError):
    """Raised when samplesheet generation cannot safely infer inputs."""


class FastqGrouping(StrEnum):
    """FASTQ grouping modes supported by samplesheet generation."""

    ILLUMINA_LANES = "illumina-lanes"
    ONT_BARCODES = "ont-barcodes"


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
        fastq1_glob: Glob pattern for read 1 FASTQ files
        fastq2_glob: Glob pattern for read 2 FASTQ files
    """

    sample_id: str
    srr: str = ""
    platform: str = ""
    fastq1: str = ""
    fastq2: str = ""
    fastq1_glob: str = ""
    fastq2_glob: str = ""


def _generated_fastq_extension(path: Path) -> str:
    name = path.name
    for suffix in sorted(GENERATED_FASTQ_SUFFIXES, key=len, reverse=True):
        if name.lower().endswith(suffix):
            return suffix
    msg = f"Unsupported FASTQ suffix for generated samplesheet input: {path}"
    raise SamplesheetGenerationError(msg)


def _get_fastq_files(directory: Path) -> list[Path]:
    """
    Find all FASTQ files in a directory.

    Searches for files with extensions: .fastq, .fq, .fastq.gz, .fq.gz

    Args:
        directory: Directory to search

    Returns:
        List of absolute Path objects for FASTQ files, sorted alphabetically
    """
    fastq_files: list[Path] = []
    directory = directory.absolute()

    for ext in GENERATED_FASTQ_SUFFIXES:
        # Handle double extensions like .fastq.gz
        if ext.startswith((".fastq", ".fq")):
            fastq_files.extend(directory.glob(f"*{ext}"))

    # Deduplicate (glob patterns may overlap) and sort
    return sorted(set(fastq_files))


def _ensure_unique_generated_sample_ids(samples: list[SampleEntry]) -> None:
    """Reject generated samplesheets that would contain duplicate sample IDs."""
    seen: set[str] = set()
    duplicates: set[str] = set()
    for sample in samples:
        if sample.sample_id in seen:
            duplicates.add(sample.sample_id)
        seen.add(sample.sample_id)

    if not duplicates:
        return

    duplicate_list = ", ".join(f"'{sample_id}'" for sample_id in sorted(duplicates))
    msg = (
        f"Generated sample_id values are not unique: {duplicate_list}. "
        "If these rows are separate Illumina lanes for the same biological sample, "
        "rerun with --group-by illumina-lanes. If they should remain separate samples, omit "
        "--sanitize or use more specific FASTQ filenames."
    )
    raise SamplesheetGenerationError(msg)


def _strip_terminal_casava_sample_lane_suffix(stem: str) -> str:
    """Strip Illumina/CASAVA sample and lane suffixes from a sample stem."""
    return CASAVA_SAMPLE_LANE_SUFFIX_PATTERN.sub("", stem)


def _extract_stem(filename: str, *, sanitize: bool = False) -> str | None:
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
        stem = match.group(1)
        return _strip_terminal_casava_sample_lane_suffix(stem) if sanitize else stem

    # Fallback: strip known extensions for single-end files
    name = filename
    for ext in sorted(GENERATED_FASTQ_SUFFIXES, key=len, reverse=True):
        if name.lower().endswith(ext):  # ty:ignore[invalid-argument-type]
            name = name[: -len(ext)]
            break

    if sanitize:
        name = _strip_terminal_casava_sample_lane_suffix(name)

    return name or None


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


def _glob_for_casava_group(
    directory: Path,
    prefix: str,
    read: int,
    extension: str,
) -> str:
    """Build a tight glob for all CASAVA lanes/chunks in one read direction."""
    escaped_prefix_path = glob.escape(str(directory / prefix))
    return f"{escaped_prefix_path}_L[0-9][0-9][0-9]_R{read}_[0-9][0-9][0-9]{extension}"


def _scan_fastq_directory_grouped_by_illumina_lanes(
    fastq_files: list[Path],
    *,
    sanitize: bool,
) -> tuple[list[SampleEntry], list[str]]:
    """Scan Illumina lane FASTQs and emit one glob-backed row per sample."""
    parsed_names: list[CasavaLaneName] = []
    unparseable: list[Path] = []
    for fq in fastq_files:
        parsed = parse_casava_lane_name(fq)
        if parsed is None:
            unparseable.append(fq)
        else:
            parsed_names.append(parsed)

    if unparseable:
        examples = ", ".join(path.name for path in unparseable[:3])
        msg = (
            "--group-by illumina-lanes requires CASAVA-style Illumina lane FASTQ names such "
            "as sample_S1_L001_R1_001.fastq.gz and "
            f"sample_S1_L001_R2_001.fastq.gz. Could not parse: {examples}"
        )
        raise SamplesheetGenerationError(msg)

    prefix_to_names: dict[str, list[CasavaLaneName]] = {}
    for parsed in parsed_names:
        prefix_to_names.setdefault(parsed.prefix, []).append(parsed)

    samples: list[SampleEntry] = []
    sample_id_to_prefix: dict[str, str] = {}
    for prefix, names in sorted(prefix_to_names.items()):
        extensions = {name.extension for name in names}
        if len(extensions) != 1:
            msg = (
                f"Sample {prefix} has mixed FASTQ extensions. Lane grouping requires "
                "one extension/compression type per sample."
            )
            raise SamplesheetGenerationError(msg)

        lanes: dict[tuple[str, str], set[int]] = {}
        for name in names:
            lanes.setdefault(name.lane_chunk_key, set()).add(name.read)

        for (lane, chunk), reads in sorted(lanes.items()):
            if reads == {1, 2}:
                continue
            missing = "R2" if reads == {1} else "R1"
            msg = (
                f"Sample {prefix} lane L{lane} chunk {chunk} is missing matching "
                f"{missing}. Fix the FASTQ directory before generating a grouped "
                "lane samplesheet."
            )
            raise SamplesheetGenerationError(msg)

        sample_id = (
            _strip_terminal_casava_sample_lane_suffix(prefix) if sanitize else prefix
        )
        previous_prefix = sample_id_to_prefix.get(sample_id)
        if previous_prefix is not None and previous_prefix != prefix:
            msg = (
                f"Grouping lanes would emit duplicate sample_id '{sample_id}' for "
                f"{previous_prefix} and {prefix}. Use more specific filenames or "
                "generate without --sanitize."
            )
            raise SamplesheetGenerationError(msg)
        sample_id_to_prefix[sample_id] = prefix

        directory = names[0].path.parent
        extension = extensions.pop()
        samples.append(
            SampleEntry(
                sample_id=sample_id,
                fastq1_glob=_glob_for_casava_group(directory, prefix, 1, extension),
                fastq2_glob=_glob_for_casava_group(directory, prefix, 2, extension),
            ),
        )

    return samples, []


def _scan_fastq_directory_grouped_by_ont_barcodes(
    directory: Path,
) -> tuple[list[SampleEntry], list[str]]:
    """Scan ONT barcode directories and emit one glob-backed row per barcode."""
    directory = directory.absolute()
    samples: list[SampleEntry] = []
    warnings_list: list[str] = []

    for barcode_dir in sorted(directory.iterdir()):
        if not barcode_dir.is_dir():
            continue
        if not re.fullmatch(r"barcode\d+", barcode_dir.name, flags=re.IGNORECASE):
            warnings_list.append(f"Skipping non-barcode directory: {barcode_dir.name}")
            continue

        fastq_files = _get_fastq_files(barcode_dir)
        if not fastq_files:
            warnings_list.append(
                f"Skipping empty barcode directory: {barcode_dir.name}",
            )
            continue

        extensions = {_generated_fastq_extension(path) for path in fastq_files}
        if len(extensions) != 1:
            msg = (
                f"Barcode {barcode_dir.name} has mixed FASTQ extensions. ONT barcode "
                "grouping requires one extension/compression type per barcode."
            )
            raise SamplesheetGenerationError(msg)

        extension = extensions.pop()
        samples.append(
            SampleEntry(
                sample_id=barcode_dir.name,
                fastq1_glob=f"{glob.escape(str(barcode_dir))}/*{extension}",
            ),
        )

    return samples, warnings_list


def scan_fastq_directory(
    directory: Path,
    *,
    sanitize: bool = False,
    group_by: FastqGrouping | None = None,
) -> tuple[list[SampleEntry], list[str]]:
    """
    Scan a directory for FASTQ files and pair them by sample.

    Groups files by sample stem, then pairs R1/R2 files. Files without
    read indicators are treated as single-end.

    Args:
        directory: Directory containing FASTQ files
        sanitize: Strip Illumina/CASAVA suffixes from generated sample IDs
        group_by: Optional FASTQ grouping mode before writing samplesheet rows

    Returns:
        Tuple of (list of SampleEntry, list of warning messages)
    """
    warnings_list: list[str] = []

    if group_by == FastqGrouping.ONT_BARCODES:
        return _scan_fastq_directory_grouped_by_ont_barcodes(directory)

    fastq_files = _get_fastq_files(directory)

    if not fastq_files:
        return [], [f"No FASTQ files found in {directory}"]

    if group_by == FastqGrouping.ILLUMINA_LANES:
        return _scan_fastq_directory_grouped_by_illumina_lanes(
            fastq_files,
            sanitize=sanitize,
        )

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
        sample_id = (
            _strip_terminal_casava_sample_lane_suffix(stem) if sanitize else stem
        )
        if len(files) == 1:
            # Single-end
            samples.append(
                SampleEntry(
                    sample_id=sample_id,
                    fastq1=str(files[0]),
                ),
            )
        elif len(files) == PAIRED_END_FILE_COUNT:
            # Paired-end: determine which is R1 and R2
            r1, r2 = None, None
            for f in files:
                read_num = _get_read_number(f.name)
                if read_num == 1:
                    r1 = f
                elif read_num == 2:  # noqa: PLR2004
                    r2 = f

            if r1 and r2:
                samples.append(
                    SampleEntry(
                        sample_id=sample_id,
                        fastq1=str(r1),
                        fastq2=str(r2),
                    ),
                )
            else:
                # Couldn't determine pairing, use alphabetical order
                warnings_list.append(
                    f"Could not determine R1/R2 for {stem}, using alphabetical order",
                )
                sorted_files = sorted(files, key=lambda p: p.name)
                samples.append(
                    SampleEntry(
                        sample_id=sample_id,
                        fastq1=str(sorted_files[0]),
                        fastq2=str(sorted_files[1]),
                    ),
                )
        else:
            # More than 2 files - ambiguous
            warnings_list.append(
                f"Sample {stem} has {len(files)} FASTQ files (expected 1-2), skipping",
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
                line = line.strip()  # noqa: PLW2901

                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue

                # Basic SRA accession validation (SRR, ERR, DRR followed by digits)
                if not looks_like_sra_accession(line):
                    warnings_list.append(
                        f"Line {line_num}: '{line}' doesn't look like an SRA accession",
                    )

                # Use accession as sample_id
                samples.append(
                    SampleEntry(
                        sample_id=line,
                        srr=line,
                    ),
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
    columns = (
        OUTPUT_COLUMNS
        if any(sample.fastq1_glob or sample.fastq2_glob for sample in samples)
        else REQUIRED_COLUMNS
    )

    with output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)

        for sample in samples:
            row = [
                sample.sample_id,
                sample.srr,
                platform,
                sample.fastq1,
                sample.fastq2,
            ]
            if columns == OUTPUT_COLUMNS:
                row.extend([sample.fastq1_glob, sample.fastq2_glob])
            writer.writerow(row)


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
        fastq1 = sample.fastq1 or sample.fastq1_glob
        fastq2 = sample.fastq2 or sample.fastq2_glob
        fastq1_display = Path(fastq1).name if fastq1 else "-"
        fastq2_display = Path(fastq2).name if fastq2 else "-"

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


samplesheet_app = typer.Typer(
    name="samplesheet",
    help="Samplesheet generation and validation",
    no_args_is_help=True,
)


@samplesheet_app.command("generate")
@samplesheet_app.command("gen", hidden=True)  # Alias
def generate(
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
    sanitize: bool = typer.Option(
        False,
        "--sanitize",
        help="Strip Illumina/CASAVA suffixes from generated sample IDs",
    ),
    group_by: FastqGrouping | None = typer.Option(
        None,
        "--group-by",
        help=(
            "Group discovered FASTQs before writing samplesheet rows. "
            "Choices: illumina-lanes, ont-barcodes."
        ),
    ),
) -> None:
    """
    Generate a samplesheet CSV from FASTQ files or SRA accessions.

    Scans a directory for FASTQ files, pairs them by sample name, and creates
    a properly formatted samplesheet. The generated samplesheet is validated
    before writing.

    Examples:

        # Generate from FASTQ directory (will prompt for platform)
        nvd samplesheet generate --from-dir ./fastqs -o samples.csv

        # Generate with platform specified (no prompts)
        nvd samplesheet generate -d ./fastqs -o samples.csv --platform illumina

        # Generate from SRA accession list
        nvd samplesheet generate --from-sra accessions.txt -o samples.csv -p illumina

        # Preview without writing
        nvd samplesheet generate -d ./fastqs --dry-run

        # Strip Illumina/CASAVA suffixes from sample IDs
        nvd samplesheet generate -d ./fastqs -p illumina --sanitize

        # Group Illumina lanes into glob-backed samplesheet rows
        nvd samplesheet generate -d ./fastqs -p illumina --group-by illumina-lanes --sanitize

        # Group ONT barcode directories into glob-backed samplesheet rows
        nvd samplesheet generate -d ./fastq_pass -p ont --group-by ont-barcodes
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
        try:
            platform, platform_warnings = normalize_platform(
                "generated samplesheet",
                platform,
            )
        except ResolutionError as exc:
            error(str(exc))
        for warn in platform_warnings:
            warning(warn)
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
        except PromptTimeoutError:
            console.print()
            error(
                f"Prompt timed out after {PROMPT_TIMEOUT_SECONDS // 60} minutes. "
                "Use --platform to specify non-interactively.",
            )
        except KeyboardInterrupt:
            console.print("\n[yellow]Cancelled[/yellow]")
            raise typer.Abort from None

    # At this point, platform is guaranteed to be a string
    # (either provided via CLI or obtained from prompt)
    assert platform is not None

    if group_by is not None:
        if from_dir is None:
            error("--group-by requires --from-dir")
        if group_by == FastqGrouping.ILLUMINA_LANES and platform != "illumina":
            error("--group-by illumina-lanes requires --platform illumina")
        if group_by == FastqGrouping.ONT_BARCODES and platform != "ont":
            error("--group-by ont-barcodes requires --platform ont")

    # Scan for samples
    if from_dir is not None:
        info(f"Scanning directory: {from_dir}")
        try:
            samples, scan_warnings = scan_fastq_directory(
                from_dir,
                sanitize=sanitize,
                group_by=group_by,
            )
        except SamplesheetGenerationError as exc:
            error(str(exc))
    else:
        assert from_sra is not None  # for type checker
        info(f"Reading SRA accessions: {from_sra}")
        samples, scan_warnings = parse_sra_accessions(from_sra)

    # Display warnings from scanning
    for warn in scan_warnings:
        warning(warn)

    if not samples:
        error("No samples found")

    try:
        _ensure_unique_generated_sample_ids(samples)
    except SamplesheetGenerationError as exc:
        error(str(exc))

    # Display preview
    console.print()
    table = format_samples_table(samples, platform)
    console.print(table)
    console.print()

    # Summary
    paired = sum(1 for s in samples if s.fastq2 or s.fastq2_glob)
    single = sum(
        1
        for s in samples
        if (s.fastq1 or s.fastq1_glob) and not (s.fastq2 or s.fastq2_glob)
    )
    sra = sum(1 for s in samples if s.srr)
    grouped = sum(1 for s in samples if s.fastq1_glob)

    info(f"Found {len(samples)} samples:")
    if grouped:
        info(f"  • {grouped} lane-grouped")
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
        except PromptTimeoutError:
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
                f"Write samplesheet to {output}?",
                default=True,
                timeout_seconds=PROMPT_TIMEOUT_SECONDS,
            ):
                info("Aborted")
                raise typer.Abort from None
        except PromptTimeoutError:
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
    render_samplesheet_validation_result(result)
