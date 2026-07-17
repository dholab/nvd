"""Typed FASTX profile parsing and report-specific presentation."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, TypeAlias

from pydantic import (
    BaseModel,
    ConfigDict,
    NonNegativeInt,
    PositiveInt,
    ValidationError,
)

from py_nvd.multiqc_packages import PackageEnvelopeError, ProfileReceipt, ReportPackage

if TYPE_CHECKING:
    from pathlib import Path

FASTX_MAX_POINTS = 1000
Number: TypeAlias = int | float


class ProducerModel(BaseModel):
    """Immutable projection of a trusted producer artifact."""

    model_config = ConfigDict(extra="ignore", frozen=True, strict=True)


class LengthSummary(ProducerModel):
    min: NonNegativeInt | None
    max: NonNegativeInt | None
    mean: Number | None
    n50: PositiveInt | None
    l50: PositiveInt | None


class QualitySummary(ProducerModel):
    mean: Number | None


class ThresholdSummary(ProducerModel):
    name: str
    sequence_count_at_or_above: NonNegativeInt


class FastxProfile(ProducerModel):
    sample_id: str
    stage: str
    sequence_count: NonNegativeInt
    total_bases: NonNegativeInt
    length: LengthSummary
    quality: QualitySummary
    thresholds: tuple[ThresholdSummary, ...]


class HistogramRow(ProducerModel):
    sample_id: str
    stage: str
    bin_start: NonNegativeInt
    sequence_count: PositiveInt


class FastxRow(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    sample_id: str
    stage: str
    query_class: str | None = None
    producer: str | None = None
    status: Literal["observed", "complete_empty", "invalid"]
    sequence_count: NonNegativeInt | None = None
    total_bases: NonNegativeInt | None = None
    min_length: NonNegativeInt | None = None
    max_length: NonNegativeInt | None = None
    mean_length: Number | None = None
    n50: PositiveInt | None = None
    l50: PositiveInt | None = None
    mean_quality: Number | None = None
    thresholds: str | None = None
    validation_message: str | None = None


@dataclass(frozen=True, slots=True)
class ParsedFastxProfile:
    receipt: ProfileReceipt
    profile: FastxProfile
    length_histogram: tuple[HistogramRow, ...]
    quality_histogram: tuple[HistogramRow, ...]


FastxParseError: TypeAlias = (
    ValidationError | UnicodeDecodeError | csv.Error | PackageEnvelopeError
)


@dataclass(frozen=True, slots=True)
class InvalidFastxProfile:
    receipt: ProfileReceipt
    error: FastxParseError


def parse_histogram(path: Path) -> tuple[HistogramRow, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(HistogramRow.model_fields).issubset(
            reader.fieldnames,
        ):
            message = "histogram is missing required columns"
            raise PackageEnvelopeError(message)
        return tuple(
            HistogramRow.model_validate_strings(row, strict=True) for row in reader
        )


def parse_fastx_package(package: ReportPackage) -> ParsedFastxProfile:
    receipt = package.receipt
    if not isinstance(receipt, ProfileReceipt):
        message = "FASTX parser received a non-profile package"
        raise TypeError(message)
    profile = FastxProfile.model_validate_json(
        package.payload(receipt.profile).read_text(encoding="utf-8"),
        strict=True,
    )
    if profile.sample_id != receipt.sample_id or profile.stage != receipt.stage:
        message = "FASTX profile identity does not match its receipt"
        raise PackageEnvelopeError(message)
    length = parse_histogram(package.payload(receipt.length_histogram))
    quality = (
        parse_histogram(package.payload(receipt.quality_histogram))
        if receipt.quality_histogram is not None
        else ()
    )
    if any(
        row.sample_id != receipt.sample_id or row.stage != receipt.stage
        for rows in (length, quality)
        for row in rows
    ):
        message = "FASTX histogram identity does not match its receipt"
        raise PackageEnvelopeError(message)
    return ParsedFastxProfile(
        receipt=receipt,
        profile=profile,
        length_histogram=length,
        quality_histogram=quality,
    )


def collect_fastx_profiles(
    packages: tuple[ReportPackage, ...],
    *,
    domain: Literal["fastx", "assembly"],
) -> tuple[ParsedFastxProfile | InvalidFastxProfile, ...]:
    profiles: list[ParsedFastxProfile | InvalidFastxProfile] = []
    for package in packages:
        receipt = package.receipt
        if not isinstance(receipt, ProfileReceipt) or receipt.domain != domain:
            continue
        try:
            profiles.append(parse_fastx_package(package))
        except (
            ValidationError,
            UnicodeDecodeError,
            csv.Error,
            PackageEnvelopeError,
        ) as error:
            profiles.append(InvalidFastxProfile(receipt=receipt, error=error))
    return tuple(profiles)


def fastx_row(profile_input: ParsedFastxProfile | InvalidFastxProfile) -> FastxRow:
    receipt = profile_input.receipt
    if isinstance(profile_input, InvalidFastxProfile):
        return FastxRow(
            sample_id=receipt.sample_id,
            stage=receipt.stage,
            query_class=receipt.query_class,
            producer=receipt.producer,
            status="invalid",
            validation_message=str(profile_input.error),
        )
    profile = profile_input.profile
    return FastxRow(
        sample_id=receipt.sample_id,
        stage=receipt.stage,
        query_class=receipt.query_class,
        producer=receipt.producer,
        status="complete_empty" if profile.sequence_count == 0 else "observed",
        sequence_count=profile.sequence_count,
        total_bases=profile.total_bases,
        min_length=profile.length.min,
        max_length=profile.length.max,
        mean_length=profile.length.mean,
        n50=profile.length.n50,
        l50=profile.length.l50,
        mean_quality=profile.quality.mean,
        thresholds="; ".join(
            f"{item.name}={item.sequence_count_at_or_above}"
            for item in profile.thresholds
        ),
    )


def histogram_series(
    profiles: tuple[ParsedFastxProfile | InvalidFastxProfile, ...],
    *,
    quality: bool,
) -> dict[str, dict[int, int]]:
    series: dict[str, dict[int, int]] = {}
    for item in profiles:
        if isinstance(item, InvalidFastxProfile):
            continue
        rows = item.quality_histogram if quality else item.length_histogram
        if not rows:
            continue
        points: dict[int, int] = {}
        for row in rows:
            points[row.bin_start] = points.get(row.bin_start, 0) + row.sequence_count
        series[series_key(item.receipt)] = rebin_series(points)
    return series


def series_key(receipt: ProfileReceipt) -> str:
    parts = [receipt.sample_id, receipt.stage]
    for part in (receipt.query_class, receipt.producer):
        if part is not None and part not in parts:
            parts.append(part)
    return " | ".join(parts)


def rebin_series(points: dict[int, int]) -> dict[int, int]:
    """Bound presentation volume without changing total counts."""
    if len(points) <= FASTX_MAX_POINTS:
        return dict(sorted(points.items()))
    starts = sorted(points)
    span = starts[-1] - starts[0] + 1
    bin_width = max(1, math.ceil(span / FASTX_MAX_POINTS))
    origin = starts[0]
    rebinned: dict[int, int] = {}
    for start in starts:
        presentation_start = origin + ((start - origin) // bin_width) * bin_width
        rebinned[presentation_start] = (
            rebinned.get(presentation_start, 0) + points[start]
        )
    return dict(sorted(rebinned.items()))
