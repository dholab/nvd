#!/usr/bin/env python3
"""Resolve NVD samplesheet read declarations into a canonical JSONL manifest."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

Source = Literal["single_file", "paired_files", "single_glob", "paired_globs", "sra"]

SUPPORTED_COMPRESSIONS = {"none", "gz", "zst", "xz"}
SUPPORTED_FASTQ_SUFFIXES = (
    ".fastq",
    ".fq",
    ".fastq.gz",
    ".fq.gz",
    ".fastq.zst",
    ".fq.zst",
    ".fastq.xz",
    ".fq.xz",
)
PLATFORM_ALIASES = {
    "illumina": "illumina",
    "ont": "ont",
    "nanopore": "ont",
    "oxford_nanopore": "ont",
    "oxford-nanopore": "ont",
}
CASAVA_NAME = re.compile(r"^(?P<prefix>.+)_L(?P<lane>\d{3})_R(?P<read>[12])_(?P<chunk>\d{3})$")
NATURAL_PART = re.compile(r"(\d+)")


class ResolutionError(ValueError):
    """Raised when a samplesheet declaration cannot be safely resolved."""


@dataclass(frozen=True)
class SamplesheetRow:
    sample_id: str
    srr: str
    platform: str
    fastq1: str
    fastq2: str
    fastq1_glob: str
    fastq2_glob: str


@dataclass(frozen=True)
class PairKey:
    prefix: str
    lane: str
    chunk: str


@dataclass(frozen=True)
class ResolvedSample:
    sample_id: str
    platform: str
    source: Source
    r1: tuple[Path, ...] = ()
    r2: tuple[Path, ...] = ()
    reads: tuple[Path, ...] = ()
    srr: str = ""
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class Resolution:
    samples: tuple[ResolvedSample, ...]
    warnings: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve NVD samplesheet read declarations into JSONL.",
    )
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--output-jsonl", required=True, type=Path)
    return parser.parse_args()


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def read_rows(samplesheet: Path) -> tuple[SamplesheetRow, ...]:
    with samplesheet.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return tuple(
            SamplesheetRow(
                sample_id=clean(row.get("sample_id")),
                srr=clean(row.get("srr")),
                platform=clean(row.get("platform")),
                fastq1=clean(row.get("fastq1")),
                fastq2=clean(row.get("fastq2")),
                fastq1_glob=clean(row.get("fastq1_glob")),
                fastq2_glob=clean(row.get("fastq2_glob")),
            )
            for row in reader
            if not clean(row.get("sample_id")).startswith("#")
        )


def absolute_preserving_symlink(path: str, cwd: Path) -> Path:
    expanded = Path(path).expanduser()
    if not expanded.is_absolute():
        expanded = cwd / expanded
    return Path(os.path.abspath(os.fspath(expanded)))


def normalize_platform(sample_id: str, platform: str) -> tuple[str, tuple[str, ...]]:
    if not platform:
        message = (
            f'Invalid platform for sample "{sample_id}"\n\n'
            "The platform column is required. Supported values are:\n\n"
            "  illumina\n"
            "  ont"
        )
        raise ResolutionError(message)

    key = platform.lower()
    canonical = PLATFORM_ALIASES.get(key)
    if canonical is None:
        message = (
            f'Invalid platform for sample "{sample_id}"\n\n'
            "NVD supports these platform values:\n\n"
            "  illumina\n"
            "  ont\n\n"
            "Aliases are also accepted:\n\n"
            "  nanopore -> ont\n"
            "  oxford_nanopore -> ont\n"
            "  oxford-nanopore -> ont\n\n"
            f"Found:\n\n  {platform}"
        )
        raise ResolutionError(message)

    if canonical == key:
        return canonical, ()
    return canonical, (f"platform alias '{platform}' was normalized to '{canonical}'",)


def has_glob_magic(value: str) -> bool:
    return glob.has_magic(value)


def reject_magic_in_exact_path(sample_id: str, column: str, value: str) -> None:
    if not has_glob_magic(value):
        return
    glob_column = f"{column}_glob"
    message = (
        f'Invalid FASTQ declaration for sample "{sample_id}"\n\n'
        f"The {column} column contains glob characters:\n\n"
        f"  {value}\n\n"
        f"The fastq1 and fastq2 columns are for exact file paths only. "
        f"Move this pattern to {glob_column} if it should be expanded."
    )
    raise ResolutionError(message)


def supported_fastq_suffix(path: Path) -> bool:
    return path.name.endswith(SUPPORTED_FASTQ_SUFFIXES)


def compression_kind(path: Path) -> str:
    name = path.name
    if name.endswith(".gz"):
        return "gz"
    if name.endswith(".zst"):
        return "zst"
    if name.endswith(".xz"):
        return "xz"
    return "none"


def validate_supported_fastqs(sample_id: str, paths: tuple[Path, ...]) -> None:
    invalid = [path for path in paths if not supported_fastq_suffix(path)]
    if not invalid:
        return
    supported = ", ".join(SUPPORTED_FASTQ_SUFFIXES)
    message = (
        f'Invalid FASTQ declaration for sample "{sample_id}"\n\n'
        "These files do not use a supported FASTQ suffix:\n\n"
        + "\n".join(f"  {path}" for path in invalid)
        + f"\n\nSupported suffixes are: {supported}"
    )
    raise ResolutionError(message)


def validate_one_compression(sample_id: str, paths: tuple[Path, ...]) -> None:
    compressions = {compression_kind(path) for path in paths}
    if compressions <= SUPPORTED_COMPRESSIONS and len(compressions) <= 1:
        return
    message = (
        f'Invalid FASTQ declaration for sample "{sample_id}"\n\n'
        "Glob inputs for one sample must use one compression type. Found:\n\n"
        + "\n".join(f"  {kind}" for kind in sorted(compressions))
    )
    raise ResolutionError(message)


def validate_existing_file(sample_id: str, column: str, path: Path) -> None:
    if path.is_file():
        return
    message = (
        f'Invalid FASTQ declaration for sample "{sample_id}"\n\n'
        f"The {column} column expects an exact path to an existing file. "
        "This file was not found:\n\n"
        f"  {path}\n\n"
        "If this sample has multiple lane files, use fastq1_glob and fastq2_glob instead."
    )
    raise ResolutionError(message)


def natural_key(value: str) -> tuple[object, ...]:
    parts = NATURAL_PART.split(value)
    return tuple(int(part) if part.isdigit() else part for part in parts)


def single_sort_key(path: Path) -> tuple[tuple[object, ...], str]:
    return natural_key(path.name), str(path)


def strip_fastq_suffix(path: Path) -> str:
    name = path.name
    for suffix in sorted(SUPPORTED_FASTQ_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def parse_casava_key(sample_id: str, path: Path, expected_read: str) -> PairKey:
    stem = strip_fastq_suffix(path)
    match = CASAVA_NAME.fullmatch(stem)
    if match is None:
        message = (
            f'Invalid paired glob declaration for sample "{sample_id}"\n\n'
            "Paired glob files must use CASAVA-style Illumina lane names:\n\n"
            "  <prefix>_L001_R1_001.fastq.gz\n"
            "  <prefix>_L001_R2_001.fastq.gz\n\n"
            f"This file did not match that convention:\n\n  {path}"
        )
        raise ResolutionError(message)
    read = match.group("read")
    if read == expected_read:
        return PairKey(match.group("prefix"), match.group("lane"), match.group("chunk"))
    message = (
        f'Invalid paired glob declaration for sample "{sample_id}"\n\n'
        f"A fastq{expected_read}_glob pattern matched an R{read} file:\n\n"
        f"  {path}"
    )
    raise ResolutionError(message)


def pair_sort_key(pair_key: PairKey) -> tuple[tuple[object, ...], str, str]:
    return natural_key(pair_key.prefix), pair_key.lane, pair_key.chunk


def keyed_paths(
    sample_id: str,
    paths: tuple[Path, ...],
    expected_read: str,
) -> dict[PairKey, Path]:
    result: dict[PairKey, Path] = {}
    duplicates: list[PairKey] = []
    for path in paths:
        key = parse_casava_key(sample_id, path, expected_read)
        if key in result:
            duplicates.append(key)
        result[key] = path
    if not duplicates:
        return result
    message = (
        f'Invalid paired glob declaration for sample "{sample_id}"\n\n'
        "The glob matched duplicate lane/chunk keys:\n\n"
        + "\n".join(f"  {key.prefix}_L{key.lane}_{key.chunk}" for key in duplicates)
    )
    raise ResolutionError(message)


def resolve_glob(sample_id: str, column: str, pattern: str, cwd: Path) -> tuple[Path, ...]:
    absolute_pattern = absolute_preserving_symlink(pattern, cwd)
    matches = tuple(
        Path(match)
        for match in glob.glob(os.fspath(absolute_pattern), recursive=True)
        if Path(match).is_file()
    )
    if matches:
        return tuple(sorted(matches, key=single_sort_key))
    message = (
        f'Invalid FASTQ declaration for sample "{sample_id}"\n\n'
        f"The {column} pattern did not match any files:\n\n"
        f"  {pattern}\n\n"
        "Check that the directory exists, the pattern is quoted where needed, "
        "and the files are visible to the machine running NVD."
    )
    raise ResolutionError(message)


def shadowed_source_warnings(row: SamplesheetRow, source: Source) -> tuple[str, ...]:
    lower_sources: list[str] = []
    if source in {"single_file", "paired_files"}:
        if row.fastq1_glob:
            lower_sources.append("fastq1_glob")
        if row.fastq2_glob:
            lower_sources.append("fastq2_glob")
        if row.srr:
            lower_sources.append("srr")
        used = "fastq1/fastq2" if source == "paired_files" else "fastq1"
    elif source in {"single_glob", "paired_globs"}:
        if row.srr:
            lower_sources.append("srr")
        used = "fastq1_glob/fastq2_glob" if source == "paired_globs" else "fastq1_glob"
    else:
        lower_sources = []
        used = "srr"

    if not lower_sources:
        return ()
    ignored = join_human(lower_sources)
    return (f"using {used}; ignoring {ignored}",)


def join_human(values: list[str]) -> str:
    if len(values) <= 1:
        return "".join(values)
    return ", ".join(values[:-1]) + f", and {values[-1]}"


def resolve_exact_files(
    row: SamplesheetRow,
    platform: str,
    platform_warnings: tuple[str, ...],
    cwd: Path,
) -> ResolvedSample:
    reject_magic_in_exact_path(row.sample_id, "fastq1", row.fastq1)
    fastq1 = absolute_preserving_symlink(row.fastq1, cwd)
    validate_existing_file(row.sample_id, "fastq1", fastq1)

    if row.fastq2:
        reject_magic_in_exact_path(row.sample_id, "fastq2", row.fastq2)
        fastq2 = absolute_preserving_symlink(row.fastq2, cwd)
        validate_existing_file(row.sample_id, "fastq2", fastq2)
        validate_supported_fastqs(row.sample_id, (fastq1, fastq2))
        warnings = platform_warnings + shadowed_source_warnings(row, "paired_files")
        return ResolvedSample(
            sample_id=row.sample_id,
            platform=platform,
            source="paired_files",
            r1=(fastq1,),
            r2=(fastq2,),
            warnings=warnings,
        )

    validate_supported_fastqs(row.sample_id, (fastq1,))
    warnings = platform_warnings + shadowed_source_warnings(row, "single_file")
    return ResolvedSample(
        sample_id=row.sample_id,
        platform=platform,
        source="single_file",
        reads=(fastq1,),
        warnings=warnings,
    )


def resolve_paired_globs(
    row: SamplesheetRow,
    platform: str,
    platform_warnings: tuple[str, ...],
    cwd: Path,
) -> ResolvedSample:
    r1_files = resolve_glob(row.sample_id, "fastq1_glob", row.fastq1_glob, cwd)
    r2_files = resolve_glob(row.sample_id, "fastq2_glob", row.fastq2_glob, cwd)
    all_files = r1_files + r2_files
    validate_supported_fastqs(row.sample_id, all_files)
    validate_one_compression(row.sample_id, all_files)

    r1_by_key = keyed_paths(row.sample_id, r1_files, "1")
    r2_by_key = keyed_paths(row.sample_id, r2_files, "2")
    if r1_by_key.keys() != r2_by_key.keys():
        message = (
            f'Invalid paired glob declaration for sample "{row.sample_id}"\n\n'
            "The fastq1_glob and fastq2_glob patterns do not describe the same lane set. "
            "Every R1 lane/chunk must have exactly one matching R2 lane/chunk."
        )
        raise ResolutionError(message)

    ordered_keys = tuple(sorted(r1_by_key, key=pair_sort_key))
    warnings = platform_warnings + shadowed_source_warnings(row, "paired_globs")
    return ResolvedSample(
        sample_id=row.sample_id,
        platform=platform,
        source="paired_globs",
        r1=tuple(r1_by_key[key] for key in ordered_keys),
        r2=tuple(r2_by_key[key] for key in ordered_keys),
        warnings=warnings,
    )


def resolve_single_glob(
    row: SamplesheetRow,
    platform: str,
    platform_warnings: tuple[str, ...],
    cwd: Path,
) -> ResolvedSample:
    reads = resolve_glob(row.sample_id, "fastq1_glob", row.fastq1_glob, cwd)
    validate_supported_fastqs(row.sample_id, reads)
    validate_one_compression(row.sample_id, reads)
    warnings = platform_warnings + shadowed_source_warnings(row, "single_glob")
    return ResolvedSample(
        sample_id=row.sample_id,
        platform=platform,
        source="single_glob",
        reads=reads,
        warnings=warnings,
    )


def resolve_sra(
    row: SamplesheetRow,
    platform: str,
    platform_warnings: tuple[str, ...],
) -> ResolvedSample:
    return ResolvedSample(
        sample_id=row.sample_id,
        platform=platform,
        source="sra",
        srr=row.srr,
        warnings=platform_warnings,
    )


def resolve_row(row: SamplesheetRow, cwd: Path) -> ResolvedSample:
    if not row.sample_id:
        message = "Invalid FASTQ declaration\n\nThe sample_id column is required for every row."
        raise ResolutionError(message)

    platform, platform_warnings = normalize_platform(row.sample_id, row.platform)

    if row.fastq1:
        return resolve_exact_files(row, platform, platform_warnings, cwd)
    if row.fastq2:
        message = (
            f'Invalid FASTQ declaration for sample "{row.sample_id}"\n\n'
            "fastq2 was provided without fastq1. For paired-end reads, provide both "
            "fastq1 and fastq2."
        )
        raise ResolutionError(message)
    if row.fastq1_glob and row.fastq2_glob:
        return resolve_paired_globs(row, platform, platform_warnings, cwd)
    if row.fastq1_glob:
        return resolve_single_glob(row, platform, platform_warnings, cwd)
    if row.fastq2_glob:
        message = (
            f'Invalid FASTQ declaration for sample "{row.sample_id}"\n\n'
            "fastq2_glob was provided without fastq1_glob. For paired-end lane files, "
            "provide both fastq1_glob and fastq2_glob."
        )
        raise ResolutionError(message)
    if row.srr:
        return resolve_sra(row, platform, platform_warnings)

    message = (
        f'Invalid FASTQ declaration for sample "{row.sample_id}"\n\n'
        "Each row must provide exact FASTQ path(s), FASTQ glob pattern(s), or an SRA accession."
    )
    raise ResolutionError(message)


def validate_unique_sample_ids(rows: tuple[SamplesheetRow, ...]) -> None:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for row in rows:
        if row.sample_id in seen:
            duplicates.add(row.sample_id)
        seen.add(row.sample_id)
    if not duplicates:
        return
    duplicate_list = ", ".join(f'"{sample_id}"' for sample_id in sorted(duplicates))
    message = (
        "Invalid samplesheet\n\n"
        f"The sample_id {duplicate_list} appears more than once.\n\n"
        "NVD expects one row per sample. If these rows describe separate sequencing "
        "lanes for the same sample, use fastq1_glob and fastq2_glob in a single row."
    )
    raise ResolutionError(message)


def resolve_rows(rows: tuple[SamplesheetRow, ...], cwd: Path | None = None) -> Resolution:
    working_dir = cwd or Path.cwd()
    validate_unique_sample_ids(rows)
    samples = tuple(resolve_row(row, working_dir) for row in rows)
    warnings = tuple(
        f'sample "{sample.sample_id}": {warning}'
        for sample in samples
        for warning in sample.warnings
    )
    return Resolution(samples=samples, warnings=warnings)


def sample_to_json(sample: ResolvedSample) -> dict[str, object]:
    record: dict[str, object] = {
        "sample_id": sample.sample_id,
        "platform": sample.platform,
        "source": sample.source,
    }
    if sample.r1:
        record["r1"] = [str(path) for path in sample.r1]
    if sample.r2:
        record["r2"] = [str(path) for path in sample.r2]
    if sample.reads:
        record["reads"] = [str(path) for path in sample.reads]
    if sample.srr:
        record["srr"] = sample.srr
    record["warnings"] = list(sample.warnings)
    return record


def write_jsonl(samples: tuple[ResolvedSample, ...], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for sample in samples:
            handle.write(json.dumps(sample_to_json(sample), sort_keys=True) + "\n")


def main() -> None:
    args = parse_args()
    try:
        resolution = resolve_rows(read_rows(args.samplesheet), cwd=Path.cwd())
    except ResolutionError as error:
        sys.stderr.write(f"{error}\n")
        raise SystemExit(1) from error

    write_jsonl(resolution.samples, args.output_jsonl)
    for warning in resolution.warnings:
        sys.stderr.write(f"Warning: {warning}\n")


if __name__ == "__main__":
    main()
