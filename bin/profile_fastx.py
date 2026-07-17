#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = []
# ///
"""Profile FASTA/FASTQ artifacts once and emit routing/count sidecars."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Any


BBTOOLS_HISTOGRAM_COUNT_COLUMNS = {
    ("Length", "Count"): (1,),
    ("Quality", "count1", "fraction1"): (1,),
    ("Quality", "count1", "fraction1", "count2", "fraction2"): (1, 3),
}


@dataclass(frozen=True, slots=True)
class Threshold:
    name: str
    axis: str
    value: float


@dataclass(frozen=True, slots=True)
class ThresholdSummary:
    threshold: Threshold
    sequence_count_at_or_above: int
    bases_at_or_above: int | None


@dataclass(frozen=True, slots=True)
class RawProfile:
    length_counts: Counter[int]
    quality_counts: Counter[int]


@dataclass(frozen=True, slots=True)
class ProfileRequest:
    sample_id: str
    stage: str
    input_files: tuple[Path, ...]
    input_format: str
    length_bin_size: int
    quality_bin_size: int
    thresholds: tuple[Threshold, ...]


@dataclass(frozen=True, slots=True)
class FastxProfile:
    sample_id: str
    stage: str
    input_files: tuple[str, ...]
    sequence_count: int
    total_bases: int
    min_length: int | None
    max_length: int | None
    mean_length: float | None
    n50: int | None
    l50: int | None
    mean_quality: float | None
    format: str
    length_histogram: Counter[int]
    quality_histogram: Counter[int]
    length_bin_size: int
    quality_bin_size: int
    thresholds: tuple[ThresholdSummary, ...]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Profile FASTA/FASTQ sequence counts, lengths, and qualities.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--input-format", choices=("fastq", "fasta"), required=True)
    parser.add_argument("--output-prefix")
    parser.add_argument("--length-bin-size", type=int, default=50)
    parser.add_argument("--quality-bin-size", type=int, default=1)
    parser.add_argument(
        "--threshold",
        action="append",
        default=[],
        metavar="NAME:AXIS:VALUE",
        help="Threshold metadata to embed in the profile JSON.",
    )
    return parser.parse_args(argv)


def parse_threshold(raw: str) -> Threshold:
    name, axis, value = raw.split(":", maxsplit=2)
    if axis not in {"length", "quality"}:
        message = f"Unsupported threshold axis {axis!r}; expected length or quality."
        raise ValueError(message)
    return Threshold(name=name, axis=axis, value=float(value))


def bin_start(value: float, bin_size: int) -> int:
    return int(value) // bin_size * bin_size


def parse_bbtools_histogram(path: Path, *, expected_axis: str) -> Counter[int]:
    counts: Counter[int] = Counter()
    header: tuple[str, ...] | None = None
    count_columns: tuple[int, ...] | None = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                header = tuple(stripped.removeprefix("#").split())
                count_columns = BBTOOLS_HISTOGRAM_COUNT_COLUMNS.get(header)
                if count_columns is None:
                    message = (
                        f"Unsupported BBTools histogram header in {path}: {header}"
                    )
                    raise ValueError(message)
                if header[0] != expected_axis:
                    message = (
                        f"BBTools histogram axis is {header[0]}; "
                        f"expected {expected_axis} in {path}"
                    )
                    raise ValueError(message)
                continue
            if header is None or count_columns is None:
                message = f"BBTools histogram data precedes its header in {path}"
                raise ValueError(message)
            fields = stripped.split()
            if len(fields) != len(header):
                message = (
                    f"BBTools histogram row has {len(fields)} columns; "
                    f"expected {len(header)} in {path}: {stripped}"
                )
                raise ValueError(message)
            value = int(fields[0])
            count = sum(int(fields[column]) for column in count_columns)
            if count:
                counts[value] += count
    if header is None:
        message = f"BBTools histogram is missing its header in {path}"
        raise ValueError(message)
    return counts


def run_reformat(input_file: Path, *, tmpdir: Path) -> RawProfile:
    reformat = shutil.which("reformat.sh")
    if reformat is None:
        message = "Could not find reformat.sh on PATH; BBTools is required."
        raise FileNotFoundError(message)
    length_histogram = tmpdir / f"{input_file.name}.length_histogram.tsv"
    quality_histogram = tmpdir / f"{input_file.name}.quality_histogram.tsv"
    subprocess.run(  # noqa: S603
        [
            reformat,
            f"in={input_file}",
            "out=/dev/null",
            "int=f",
            f"lhist={length_histogram}",
            f"aqhist={quality_histogram}",
            "overwrite=true",
        ],
        check=True,
    )
    return RawProfile(
        length_counts=parse_bbtools_histogram(
            length_histogram,
            expected_axis="Length",
        ),
        quality_counts=parse_bbtools_histogram(
            quality_histogram,
            expected_axis="Quality",
        ),
    )


def merge_profiles(profiles: Iterator[RawProfile]) -> RawProfile:
    length_counts: Counter[int] = Counter()
    quality_counts: Counter[int] = Counter()
    for profile in profiles:
        length_counts.update(profile.length_counts)
        quality_counts.update(profile.quality_counts)
    return RawProfile(length_counts=length_counts, quality_counts=quality_counts)


def rebin_counts(counts: Counter[int], bin_size: int) -> Counter[int]:
    rebinned: Counter[int] = Counter()
    for value, count in counts.items():
        rebinned[bin_start(value, bin_size)] += count
    return rebinned


def weighted_mean(counts: Counter[int]) -> float | None:
    total_count = sum(counts.values())
    if not total_count:
        return None
    return sum(value * count for value, count in counts.items()) / total_count


def n50_l50(length_counts: Counter[int]) -> tuple[int | None, int | None]:
    total_bases = sum(length * count for length, count in length_counts.items())
    if total_bases == 0:
        return None, None
    target = (total_bases + 1) // 2
    cumulative_bases = 0
    cumulative_sequences = 0
    for length in sorted(length_counts, reverse=True):
        count = length_counts[length]
        bases = length * count
        if cumulative_bases + bases >= target:
            needed = (target - cumulative_bases + length - 1) // length
            return length, cumulative_sequences + max(needed, 1)
        cumulative_bases += bases
        cumulative_sequences += count
    raise AssertionError("N50/L50 scan exhausted nonempty length counts")


def summarize_threshold(
    threshold: Threshold,
    *,
    length_counts: Counter[int],
    quality_counts: Counter[int],
) -> ThresholdSummary:
    counts = length_counts if threshold.axis == "length" else quality_counts
    sequence_count = 0
    bases = 0
    for value, count in counts.items():
        if value < threshold.value:
            continue
        sequence_count += count
        if threshold.axis == "length":
            bases += value * count
    return ThresholdSummary(
        threshold=threshold,
        sequence_count_at_or_above=sequence_count,
        bases_at_or_above=bases if threshold.axis == "length" else None,
    )


def profile_fastx(request: ProfileRequest) -> FastxProfile:
    with tempfile.TemporaryDirectory(prefix="nvd-fastx-profile-") as tmpdir_name:
        raw = merge_profiles(
            run_reformat(input_file, tmpdir=Path(tmpdir_name))
            for input_file in request.input_files
        )
    sequence_count = sum(raw.length_counts.values())
    total_bases = sum(length * count for length, count in raw.length_counts.items())
    quality_count = sum(raw.quality_counts.values())
    if request.input_format == "fastq" and quality_count != sequence_count:
        message = "FASTQ quality histogram count must equal sequence count"
        raise ValueError(message)
    if request.input_format == "fasta" and quality_count != 0:
        message = "FASTA input must not produce quality histogram counts"
        raise ValueError(message)
    n50, l50 = n50_l50(raw.length_counts)
    return FastxProfile(
        sample_id=request.sample_id,
        stage=request.stage,
        input_files=tuple(str(path) for path in request.input_files),
        sequence_count=sequence_count,
        total_bases=total_bases,
        min_length=min(raw.length_counts) if raw.length_counts else None,
        max_length=max(raw.length_counts) if raw.length_counts else None,
        mean_length=(total_bases / sequence_count) if sequence_count else None,
        n50=n50,
        l50=l50,
        mean_quality=weighted_mean(raw.quality_counts),
        format="fastq" if request.input_format == "fastq" else "fasta_or_empty",
        length_histogram=rebin_counts(raw.length_counts, request.length_bin_size),
        quality_histogram=rebin_counts(raw.quality_counts, request.quality_bin_size),
        length_bin_size=request.length_bin_size,
        quality_bin_size=request.quality_bin_size,
        thresholds=tuple(
            summarize_threshold(
                threshold,
                length_counts=raw.length_counts,
                quality_counts=raw.quality_counts,
            )
            for threshold in request.thresholds
        ),
    )


def histogram_rows(
    *,
    profile: FastxProfile,
    counts: Counter[int],
    bin_size: int,
) -> Iterator[dict[str, str | int]]:
    for start, count in sorted(counts.items()):
        yield {
            "sample_id": profile.sample_id,
            "stage": profile.stage,
            "bin_start": start,
            "bin_end": start + bin_size - 1,
            "sequence_count": count,
        }


def write_histogram(
    path: Path,
    *,
    profile: FastxProfile,
    counts: Counter[int],
    bin_size: int,
) -> None:
    fieldnames = ["sample_id", "stage", "bin_start", "bin_end", "sequence_count"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            histogram_rows(profile=profile, counts=counts, bin_size=bin_size),
        )


def profile_json(profile: FastxProfile) -> dict[str, Any]:
    return {
        "sample_id": profile.sample_id,
        "stage": profile.stage,
        "input_files": list(profile.input_files),
        "format": profile.format,
        "sequence_count": profile.sequence_count,
        "total_bases": profile.total_bases,
        "length": {
            "min": profile.min_length,
            "max": profile.max_length,
            "mean": profile.mean_length,
            "n50": profile.n50,
            "l50": profile.l50,
            "bin_size": profile.length_bin_size,
        },
        "quality": {
            "mean": profile.mean_quality,
            "bin_size": profile.quality_bin_size,
            "encoding": "phred33_integer_average"
            if profile.mean_quality is not None
            else None,
        },
        "thresholds": [
            {
                "name": summary.threshold.name,
                "axis": summary.threshold.axis,
                "value": summary.threshold.value,
                "sequence_count_at_or_above": summary.sequence_count_at_or_above,
                "bases_at_or_above": summary.bases_at_or_above,
            }
            for summary in profile.thresholds
        ],
    }


def write_profile(profile: FastxProfile, output_prefix: str) -> None:
    Path(f"{output_prefix}.fastx_profile.json").write_text(
        json.dumps(profile_json(profile), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    Path(f"{output_prefix}.sequence_count.txt").write_text(
        f"{profile.sequence_count}\n",
        encoding="utf-8",
    )
    write_histogram(
        Path(f"{output_prefix}.length_histogram.tsv"),
        profile=profile,
        counts=profile.length_histogram,
        bin_size=profile.length_bin_size,
    )
    if profile.quality_histogram:
        write_histogram(
            Path(f"{output_prefix}.quality_histogram.tsv"),
            profile=profile,
            counts=profile.quality_histogram,
            bin_size=profile.quality_bin_size,
        )


def strip_fastx_suffix(path: Path) -> Path:
    return next(
        (
            path.with_name(path.name[: -len(suffix)])
            for suffix in (
                ".fastq.gz",
                ".fq.gz",
                ".fasta.gz",
                ".fa.gz",
                ".fastq",
                ".fq",
                ".fasta",
                ".fa",
            )
            if path.name.endswith(suffix)
        ),
        path.with_suffix(""),
    )


def output_prefix(
    *,
    sample_id: str,
    stage: str,
    input_files: Sequence[Path],
    explicit_prefix: str | None,
) -> str:
    if explicit_prefix:
        return explicit_prefix
    if len(input_files) == 1:
        return str(strip_fastx_suffix(input_files[0]))
    return f"{sample_id}.{stage}"


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    write_profile(
        profile_fastx(
            ProfileRequest(
                sample_id=args.sample_id,
                stage=args.stage,
                input_files=tuple(args.input),
                input_format=args.input_format,
                length_bin_size=args.length_bin_size,
                quality_bin_size=args.quality_bin_size,
                thresholds=tuple(
                    parse_threshold(threshold) for threshold in args.threshold
                ),
            ),
        ),
        output_prefix(
            sample_id=args.sample_id,
            stage=args.stage,
            input_files=args.input,
            explicit_prefix=args.output_prefix,
        ),
    )


if __name__ == "__main__":
    main()
