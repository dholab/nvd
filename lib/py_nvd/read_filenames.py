"""Shared FASTQ filename, compression, and pairing primitives."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from pathlib import Path

Compression = Literal["none", "gz", "xz", "zst"]

SUPPORTED_COMPRESSIONS: frozenset[Compression] = frozenset({"none", "gz", "xz", "zst"})
SUPPORTED_FASTQ_SUFFIXES = (
    ".fastq",
    ".fq",
    ".fastq.gz",
    ".fq.gz",
    ".fastq.xz",
    ".fq.xz",
    ".fastq.zst",
    ".fq.zst",
)
GENERATED_FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

CASAVA_LANE_NAME = re.compile(
    r"^(?P<prefix>.+)_L(?P<lane>\d{3})_R(?P<read>[12])_"
    r"(?P<chunk>\d{3})(?P<extension>\.(?:fastq|fq)(?:\.(?:gz|xz|zst))?)$",
    re.IGNORECASE,
)
TERMINAL_READ_MARKER = re.compile(
    r"^(?P<prefix>.+)(?P<separator>[_.])(?P<r>R?)(?P<read>[12])(?P<chunk>_\d+)?$",
    re.IGNORECASE,
)
NATURAL_PART = re.compile(r"(\d+)")


@dataclass(frozen=True)
class CasavaLaneName:
    """Parsed Illumina/CASAVA lane FASTQ filename components."""

    path: Path
    prefix: str
    lane: str
    read: int
    chunk: str
    extension: str

    @property
    def lane_chunk_key(self) -> tuple[str, str]:
        """Key that must have both R1 and R2 for paired lane grouping."""
        return (self.lane, self.chunk)


@dataclass(frozen=True)
class ReadPairKey:
    """Stable key for pairing R1/R2 FASTQs without relying on sorted order."""

    prefix: str
    style: str
    lane: str = ""
    chunk: str = ""


@dataclass(frozen=True)
class ReadMarker:
    """Read number and pair key parsed from a FASTQ filename."""

    read: str
    key: ReadPairKey


def natural_key(value: str) -> tuple[object, ...]:
    parts = NATURAL_PART.split(value)
    return tuple(int(part) if part.isdigit() else part for part in parts)


def supported_fastq_suffix(path: Path) -> bool:
    return path.name.endswith(SUPPORTED_FASTQ_SUFFIXES)


def compression_kind(path: Path) -> Compression:
    name = path.name
    if name.endswith(".gz"):
        return "gz"
    if name.endswith(".xz"):
        return "xz"
    if name.endswith(".zst"):
        return "zst"
    return "none"


def strip_fastq_suffix(path: Path) -> str:
    name = path.name
    for suffix in sorted(SUPPORTED_FASTQ_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def parse_casava_lane_name(path: Path) -> CasavaLaneName | None:
    match = CASAVA_LANE_NAME.fullmatch(path.name)
    if match is None:
        return None

    return CasavaLaneName(
        path=path,
        prefix=match.group("prefix"),
        lane=match.group("lane"),
        read=int(match.group("read")),
        chunk=match.group("chunk"),
        extension=match.group("extension"),
    )


def parse_read_marker(path: Path) -> ReadMarker | None:
    casava = parse_casava_lane_name(path)
    if casava is not None:
        return ReadMarker(
            read=str(casava.read),
            key=ReadPairKey(
                prefix=casava.prefix,
                style="casava",
                lane=casava.lane,
                chunk=casava.chunk,
            ),
        )

    stem = strip_fastq_suffix(path)
    match = TERMINAL_READ_MARKER.fullmatch(stem)
    if match is None:
        return None

    marker = match.group("r").upper()
    chunk = match.group("chunk") or ""
    style = f"terminal:{match.group('separator')}{marker}:{chunk}"
    return ReadMarker(
        read=match.group("read"),
        key=ReadPairKey(
            prefix=match.group("prefix"),
            style=style,
            chunk=chunk,
        ),
    )


def read_pair_sort_key(
    pair_key: ReadPairKey,
) -> tuple[tuple[object, ...], str, str, str]:
    return natural_key(pair_key.prefix), pair_key.style, pair_key.lane, pair_key.chunk
