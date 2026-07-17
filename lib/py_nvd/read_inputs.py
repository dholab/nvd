"""Resolve NVD samplesheet read declarations into a canonical JSONL manifest."""

from __future__ import annotations

import csv
import glob
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from py_nvd.read_filenames import (
    SUPPORTED_COMPRESSIONS,
    SUPPORTED_FASTQ_SUFFIXES,
    ReadPairKey,
    compression_kind,
    natural_key,
    parse_read_marker,
    read_pair_sort_key,
    supported_fastq_suffix,
)
from py_nvd.sra_accessions import looks_like_sra_accession

Source = Literal["single_file", "paired_files", "single_glob", "paired_globs", "sra"]
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")

PLATFORM_ALIASES = {
    "illumina": "illumina",
    "ont": "ont",
    "nanopore": "ont",
    "oxford_nanopore": "ont",
    "oxford-nanopore": "ont",
}


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


def _clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def read_rows(samplesheet: Path) -> tuple[SamplesheetRow, ...]:
    with samplesheet.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return tuple(
            SamplesheetRow(
                sample_id=_clean(row.get("sample_id")),
                srr=_clean(row.get("srr")),
                platform=_clean(row.get("platform")),
                fastq1=_clean(row.get("fastq1")),
                fastq2=_clean(row.get("fastq2")),
                fastq1_glob=_clean(row.get("fastq1_glob")),
                fastq2_glob=_clean(row.get("fastq2_glob")),
            )
            for row in reader
            if not _clean(row.get("sample_id")).startswith("#")
        )


def _absolute_preserving_symlink(path: str) -> Path:
    expanded = Path(path).expanduser()
    if not expanded.is_absolute():
        message = (
            "FASTQ paths and glob patterns must be absolute. "
            "NVD resolves read inputs inside Nextflow task directories, where relative paths are ambiguous. "
            "Use an absolute path, or an absolute symlink path from the stable cluster namespace."
        )
        raise ResolutionError(message)
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


def _has_glob_magic(value: str) -> bool:
    return glob.has_magic(value)


def _reject_magic_in_exact_path(sample_id: str, column: str, value: str) -> None:
    if not _has_glob_magic(value):
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


def _validate_supported_fastqs(sample_id: str, paths: tuple[Path, ...]) -> None:
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


def _validate_one_compression(sample_id: str, paths: tuple[Path, ...]) -> None:
    compressions = {compression_kind(path) for path in paths}
    if compressions <= SUPPORTED_COMPRESSIONS and len(compressions) <= 1:
        return
    message = (
        f'Invalid FASTQ declaration for sample "{sample_id}"\n\n'
        "Glob inputs for one sample must use one compression type. Found:\n\n"
        + "\n".join(f"  {kind}" for kind in sorted(compressions))
    )
    raise ResolutionError(message)


def _validate_existing_file(sample_id: str, column: str, path: Path) -> None:
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


def _single_sort_key(path: Path) -> tuple[tuple[object, ...], str]:
    return natural_key(path.name), str(path)


def _parse_pair_key(sample_id: str, path: Path, expected_read: str) -> ReadPairKey:
    marker = parse_read_marker(path)
    if marker is None:
        message = (
            f'Invalid paired glob declaration for sample "{sample_id}"\n\n'
            "Paired glob files must use CASAVA-style Illumina lane names or "
            "simple terminal read markers:\n\n"
            "  <prefix>_L001_R1_001.fastq.gz\n"
            "  <prefix>_L001_R2_001.fastq.gz\n"
            "  <prefix>_R1.fastq.gz / <prefix>_R2.fastq.gz\n"
            "  <prefix>.1.fq.gz / <prefix>.2.fq.gz\n\n"
            f"This file did not match a supported convention:\n\n  {path}"
        )
        raise ResolutionError(message)
    if marker.read == expected_read:
        return marker.key
    message = (
        f'Invalid paired glob declaration for sample "{sample_id}"\n\n'
        f"A fastq{expected_read}_glob pattern matched an R{marker.read} file:\n\n"
        f"  {path}"
    )
    raise ResolutionError(message)


def _keyed_paths(
    sample_id: str,
    paths: tuple[Path, ...],
    expected_read: str,
) -> dict[ReadPairKey, Path]:
    result: dict[ReadPairKey, Path] = {}
    duplicates: list[ReadPairKey] = []
    for path in paths:
        key = _parse_pair_key(sample_id, path, expected_read)
        if key in result:
            duplicates.append(key)
        result[key] = path
    if not duplicates:
        return result
    message = (
        f'Invalid paired glob declaration for sample "{sample_id}"\n\n'
        "The glob matched duplicate lane/chunk keys:\n\n"
        + "\n".join(
            f"  {key.prefix} ({key.style}, lane={key.lane or '-'}, chunk={key.chunk or '-'})"
            for key in duplicates
        )
    )
    raise ResolutionError(message)


def _resolve_glob(sample_id: str, column: str, pattern: str) -> tuple[Path, ...]:
    absolute_pattern = _absolute_preserving_symlink(pattern)
    matches = tuple(
        Path(match)
        for match in glob.glob(os.fspath(absolute_pattern), recursive=True)
        if Path(match).is_file()
    )
    if matches:
        return tuple(sorted(matches, key=_single_sort_key))
    message = (
        f'Invalid FASTQ declaration for sample "{sample_id}"\n\n'
        f"The {column} pattern did not match any files:\n\n"
        f"  {pattern}\n\n"
        "Check that the directory exists, the pattern is quoted where needed, "
        "and the files are visible to the machine running NVD."
    )
    raise ResolutionError(message)


def _shadowed_source_warnings(row: SamplesheetRow, source: Source) -> tuple[str, ...]:
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
    ignored = _join_human(lower_sources)
    return (f"using {used}; ignoring {ignored}",)


def _join_human(values: list[str]) -> str:
    if len(values) <= 1:
        return "".join(values)
    return ", ".join(values[:-1]) + f", and {values[-1]}"


def _resolve_exact_files(
    row: SamplesheetRow,
    platform: str,
    platform_warnings: tuple[str, ...],
) -> ResolvedSample:
    _reject_magic_in_exact_path(row.sample_id, "fastq1", row.fastq1)
    fastq1 = _absolute_preserving_symlink(row.fastq1)
    _validate_existing_file(row.sample_id, "fastq1", fastq1)

    if row.fastq2:
        _reject_magic_in_exact_path(row.sample_id, "fastq2", row.fastq2)
        fastq2 = _absolute_preserving_symlink(row.fastq2)
        _validate_existing_file(row.sample_id, "fastq2", fastq2)
        _validate_supported_fastqs(row.sample_id, (fastq1, fastq2))
        warnings = platform_warnings + _shadowed_source_warnings(row, "paired_files")
        return ResolvedSample(
            sample_id=row.sample_id,
            platform=platform,
            source="paired_files",
            r1=(fastq1,),
            r2=(fastq2,),
            warnings=warnings,
        )

    _validate_supported_fastqs(row.sample_id, (fastq1,))
    warnings = platform_warnings + _shadowed_source_warnings(row, "single_file")
    return ResolvedSample(
        sample_id=row.sample_id,
        platform=platform,
        source="single_file",
        reads=(fastq1,),
        warnings=warnings,
    )


def _resolve_paired_globs(
    row: SamplesheetRow,
    platform: str,
    platform_warnings: tuple[str, ...],
) -> ResolvedSample:
    r1_files = _resolve_glob(row.sample_id, "fastq1_glob", row.fastq1_glob)
    r2_files = _resolve_glob(row.sample_id, "fastq2_glob", row.fastq2_glob)
    all_files = r1_files + r2_files
    _validate_supported_fastqs(row.sample_id, all_files)
    _validate_one_compression(row.sample_id, all_files)

    r1_by_key = _keyed_paths(row.sample_id, r1_files, "1")
    r2_by_key = _keyed_paths(row.sample_id, r2_files, "2")
    if r1_by_key.keys() != r2_by_key.keys():
        message = (
            f'Invalid paired glob declaration for sample "{row.sample_id}"\n\n'
            "The fastq1_glob and fastq2_glob patterns do not describe the same lane set. "
            "Every R1 lane/chunk must have exactly one matching R2 lane/chunk."
        )
        raise ResolutionError(message)

    ordered_keys = tuple(sorted(r1_by_key, key=read_pair_sort_key))
    warnings = platform_warnings + _shadowed_source_warnings(row, "paired_globs")
    return ResolvedSample(
        sample_id=row.sample_id,
        platform=platform,
        source="paired_globs",
        r1=tuple(r1_by_key[key] for key in ordered_keys),
        r2=tuple(r2_by_key[key] for key in ordered_keys),
        warnings=warnings,
    )


def _resolve_single_glob(
    row: SamplesheetRow,
    platform: str,
    platform_warnings: tuple[str, ...],
) -> ResolvedSample:
    reads = _resolve_glob(row.sample_id, "fastq1_glob", row.fastq1_glob)
    _validate_supported_fastqs(row.sample_id, reads)
    _validate_one_compression(row.sample_id, reads)
    warnings = platform_warnings + _shadowed_source_warnings(row, "single_glob")
    return ResolvedSample(
        sample_id=row.sample_id,
        platform=platform,
        source="single_glob",
        reads=reads,
        warnings=warnings,
    )


def _resolve_sra(
    row: SamplesheetRow,
    platform: str,
    platform_warnings: tuple[str, ...],
) -> ResolvedSample:
    sra_warnings = ()
    if not looks_like_sra_accession(row.srr):
        sra_warnings = (
            f"SRA accession '{row.srr}' does not look like an SRR/ERR/DRR accession",
        )
    return ResolvedSample(
        sample_id=row.sample_id,
        platform=platform,
        source="sra",
        srr=row.srr,
        warnings=platform_warnings + sra_warnings,
    )


def _resolve_row(row: SamplesheetRow) -> ResolvedSample:
    if not row.sample_id:
        message = "Invalid FASTQ declaration\n\nThe sample_id column is required for every row."
        raise ResolutionError(message)
    if not SAFE_SAMPLE_ID.fullmatch(row.sample_id):
        message = (
            f"Invalid sample_id {row.sample_id!r}.\n\n"
            "Sample IDs may contain only letters, numbers, dots, underscores, "
            "and hyphens."
        )
        raise ResolutionError(message)

    platform, platform_warnings = normalize_platform(row.sample_id, row.platform)

    if row.fastq1:
        return _resolve_exact_files(row, platform, platform_warnings)
    if row.fastq2:
        message = (
            f'Invalid FASTQ declaration for sample "{row.sample_id}"\n\n'
            "fastq2 was provided without fastq1. For paired-end reads, provide both "
            "fastq1 and fastq2."
        )
        raise ResolutionError(message)
    if row.fastq1_glob and row.fastq2_glob:
        return _resolve_paired_globs(row, platform, platform_warnings)
    if row.fastq1_glob:
        return _resolve_single_glob(row, platform, platform_warnings)
    if row.fastq2_glob:
        message = (
            f'Invalid FASTQ declaration for sample "{row.sample_id}"\n\n'
            "fastq2_glob was provided without fastq1_glob. For paired-end lane files, "
            "provide both fastq1_glob and fastq2_glob."
        )
        raise ResolutionError(message)
    if row.srr:
        return _resolve_sra(row, platform, platform_warnings)

    message = (
        f'Invalid FASTQ declaration for sample "{row.sample_id}"\n\n'
        "Each row must provide exact FASTQ path(s), FASTQ glob pattern(s), or an SRA accession."
    )
    raise ResolutionError(message)


def _validate_unique_sample_ids(rows: tuple[SamplesheetRow, ...]) -> None:
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


def resolve_rows(rows: tuple[SamplesheetRow, ...]) -> Resolution:
    if not rows:
        message = "Invalid samplesheet\n\nThe samplesheet must contain at least one non-comment sample row."
        raise ResolutionError(message)
    _validate_unique_sample_ids(rows)
    samples = tuple(_resolve_row(row) for row in rows)
    warnings = tuple(
        f'sample "{sample.sample_id}": {warning}'
        for sample in samples
        for warning in sample.warnings
    )
    return Resolution(samples=samples, warnings=warnings)


def _sample_to_json(sample: ResolvedSample) -> dict[str, object]:
    if sample.r1:
        read_structure = "paired"
        read_counts: dict[str, int] | None = {
            "R1": len(sample.r1),
            "R2": len(sample.r2),
        }
    elif sample.reads:
        read_structure = "single"
        read_counts = {"single": len(sample.reads)}
    else:
        read_structure = "unknown"
        read_counts = None
    record: dict[str, object] = {
        "sample_id": sample.sample_id,
        "platform": sample.platform,
        "source": sample.source,
        "read_structure": read_structure,
        "read_counts": read_counts,
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
            handle.write(json.dumps(_sample_to_json(sample), sort_keys=True) + "\n")
