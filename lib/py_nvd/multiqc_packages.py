"""Typed construction and parsing of internal report-package receipts."""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass
from enum import StrEnum
from typing import TYPE_CHECKING, Annotated, Literal, TypeAlias

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    StringConstraints,
    TypeAdapter,
    ValidationError,
)

if TYPE_CHECKING:
    from collections.abc import Collection, Mapping, Sequence
    from pathlib import Path

REPORT_PACKAGE_SCHEMA = "nvd.report-package/v1"

SafeName = Annotated[
    str,
    StringConstraints(min_length=1, pattern=r"^[A-Za-z0-9_.-]+$"),
]
NonEmptyString = Annotated[str, StringConstraints(min_length=1)]


class Domain(StrEnum):
    """Report package domains."""

    TARGET_ENRICHMENT = "target_enrichment"
    DEPLETION = "depletion"
    FASTX = "fastx"
    ASSEMBLY = "assembly"
    QUERY_PREPARATION = "query_preparation"
    BLAST = "blast"


class ArtifactType(StrEnum):
    """Report package artifact kinds."""

    DEACON_STATS = "deacon_stats"
    PROFILE = "profile"
    ELIGIBILITY = "eligibility"
    LONG_READ_ELIGIBILITY = "long_read_eligibility"
    UNION_SUMMARY = "union_summary"
    PREPARED_QUERY_BATCHES = "prepared_query_batches"
    MEGABLAST_QUERY_PARTITION = "megablast_query_partition"


DeaconDomain: TypeAlias = Literal[Domain.TARGET_ENRICHMENT, Domain.DEPLETION]
ProfileDomain: TypeAlias = Literal[Domain.FASTX, Domain.ASSEMBLY]


class ReceiptBase(BaseModel):
    """Fields shared by all report-package receipts."""

    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    schema_version: Literal["nvd.report-package/v1"] = REPORT_PACKAGE_SCHEMA
    domain: Domain
    artifact_type: ArtifactType
    sample_id: SafeName

    @property
    def logical_key(self) -> tuple[str, ...]:
        raise NotImplementedError

    @property
    def payload_names(self) -> tuple[str, ...]:
        return ()

    @property
    def package_name(self) -> str:
        encoded = "\0".join(self.logical_key).encode()
        digest = hashlib.sha256(encoded).hexdigest()[:16]
        return f"{self.domain}.{self.artifact_type}.{digest}.pkg"


class DeaconReceipt(ReceiptBase):
    """A target-enrichment or depletion DEACON summary."""

    domain: DeaconDomain
    artifact_type: Literal["deacon_stats"] = ArtifactType.DEACON_STATS.value
    query_class: SafeName | None = None
    stats: Literal["stats.payload"] = "stats.payload"

    @property
    def logical_key(self) -> tuple[str, ...]:
        return (
            self.domain,
            self.artifact_type,
            self.sample_id,
            self.query_class or "",
        )

    @property
    def payload_names(self) -> tuple[str, ...]:
        return (self.stats,)


class ProfileReceipt(ReceiptBase):
    """A producer-owned FASTX profile and presentation histograms."""

    domain: ProfileDomain
    artifact_type: Literal["profile"] = ArtifactType.PROFILE.value
    stage: SafeName
    query_class: SafeName | None = None
    producer: SafeName | None = None
    profile: Literal["profile.payload"] = "profile.payload"
    length_histogram: Literal["length_histogram.payload"] = "length_histogram.payload"
    quality_histogram: Literal["quality_histogram.payload"] | None = None

    @property
    def logical_key(self) -> tuple[str, ...]:
        return (
            self.domain,
            self.artifact_type,
            self.sample_id,
            self.stage,
            self.query_class or "",
            self.producer or "",
        )

    @property
    def payload_names(self) -> tuple[str, ...]:
        names = (self.profile, self.length_histogram)
        if self.quality_histogram is not None:
            return (*names, self.quality_histogram)
        return names


class EligibilityReceipt(ReceiptBase):
    """One explicit assembly routing decision."""

    domain: Literal["assembly"] = Domain.ASSEMBLY.value
    artifact_type: Literal["eligibility"] = ArtifactType.ELIGIBILITY.value
    producer: SafeName
    decision: Literal["run", "skip"]
    sequence_count: NonNegativeInt
    threshold: NonNegativeInt
    reason: NonEmptyString

    @property
    def logical_key(self) -> tuple[str, ...]:
        return (self.domain, self.artifact_type, self.sample_id, self.producer)


class LongReadEligibilityReceipt(ReceiptBase):
    """Producer summary covering long-read assembler eligibility."""

    domain: Literal["assembly"] = Domain.ASSEMBLY.value
    artifact_type: Literal["long_read_eligibility"] = (
        ArtifactType.LONG_READ_ELIGIBILITY.value
    )
    producer: Literal["experimental_long_read"] = "experimental_long_read"
    summary: Literal["summary.payload"] = "summary.payload"

    @property
    def logical_key(self) -> tuple[str, ...]:
        return (self.domain, self.artifact_type, self.sample_id, self.producer)

    @property
    def payload_names(self) -> tuple[str, ...]:
        return (self.summary,)


class UnionReceipt(ReceiptBase):
    """Producer summaries for the long-read contig union."""

    domain: Literal["assembly"] = Domain.ASSEMBLY.value
    artifact_type: Literal["union_summary"] = ArtifactType.UNION_SUMMARY.value
    producer: Literal["long_read_union"] = "long_read_union"
    summary: Literal["summary.payload"] = "summary.payload"

    @property
    def logical_key(self) -> tuple[str, ...]:
        return (self.domain, self.artifact_type, self.sample_id, self.producer)

    @property
    def payload_names(self) -> tuple[str, ...]:
        return (self.summary,)


class PreparedQueryBatchesReceipt(ReceiptBase):
    """One producer summary of prepared query batches for a sample."""

    domain: Literal["query_preparation"] = Domain.QUERY_PREPARATION.value
    artifact_type: Literal["prepared_query_batches"] = (
        ArtifactType.PREPARED_QUERY_BATCHES.value
    )
    summary: Literal["summary.payload"] = "summary.payload"

    @property
    def logical_key(self) -> tuple[str, ...]:
        return (self.domain, self.artifact_type, self.sample_id)

    @property
    def payload_names(self) -> tuple[str, ...]:
        return (self.summary,)


class MegablastQueryPartitionReceipt(ReceiptBase):
    """One MEGABLAST query partition summary for a query class."""

    domain: Literal["blast"] = Domain.BLAST.value
    artifact_type: Literal["megablast_query_partition"] = (
        ArtifactType.MEGABLAST_QUERY_PARTITION.value
    )
    query_class: SafeName
    summary: Literal["summary.payload"] = "summary.payload"

    @property
    def logical_key(self) -> tuple[str, ...]:
        return (self.domain, self.artifact_type, self.sample_id, self.query_class)

    @property
    def payload_names(self) -> tuple[str, ...]:
        return (self.summary,)


ReportReceipt: TypeAlias = Annotated[
    DeaconReceipt
    | ProfileReceipt
    | EligibilityReceipt
    | LongReadEligibilityReceipt
    | UnionReceipt
    | PreparedQueryBatchesReceipt
    | MegablastQueryPartitionReceipt,
    Field(discriminator="artifact_type"),
]
REPORT_RECEIPT = TypeAdapter(ReportReceipt)
DOMAIN_ORDER = {domain: index for index, domain in enumerate(Domain)}


class PackageEnvelopeError(ValueError):
    """A staged package does not match its typed receipt."""


@dataclass(frozen=True, slots=True)
class ReportPackage:
    """A parsed receipt paired with its confined payload directory."""

    root: Path
    receipt: ReportReceipt

    def payload(self, name: str) -> Path:
        return self.root / name


@dataclass(frozen=True, slots=True)
class InvalidPackage:
    """An optional package that failed at its parsing boundary."""

    domain: Domain
    package_name: str
    error: (
        ValidationError | UnicodeDecodeError | FileNotFoundError | PackageEnvelopeError
    )


ParsedPackage: TypeAlias = ReportPackage | InvalidPackage


def write_package(
    receipt: ReportReceipt,
    payloads: Sequence[Path],
    *,
    output_root: Path,
) -> Path:
    """Materialize one self-describing package."""
    if len(receipt.payload_names) != len(payloads):
        message = "receipt payload count does not match supplied files"
        raise ValueError(message)

    package_dir = output_root / receipt.package_name
    package_dir.mkdir()
    for target_name, source in zip(receipt.payload_names, payloads, strict=True):
        if not source.is_file():
            message = f"payload is not a regular file: {source.name}"
            raise FileNotFoundError(message)
        shutil.copyfile(source, package_dir / target_name, follow_symlinks=True)
    (package_dir / "receipt.json").write_text(
        REPORT_RECEIPT.dump_json(receipt).decode() + "\n",
        encoding="utf-8",
    )
    return package_dir


def parse_receipt(path: Path) -> ReportReceipt:
    """Parse one receipt without replacing Pydantic diagnostics."""
    return REPORT_RECEIPT.validate_json(path.read_text(encoding="utf-8"), strict=True)


def load_report_package(package_dir: Path, *, expected_domain: Domain) -> ReportPackage:
    """Parse and confine one staged package."""
    receipt_path = package_dir / "receipt.json"
    if not receipt_path.is_file():
        message = f"{package_dir.name} has no receipt.json"
        raise FileNotFoundError(message)
    receipt = parse_receipt(receipt_path)
    if receipt.domain != expected_domain:
        message = f"{package_dir.name} is staged under the wrong domain"
        raise PackageEnvelopeError(message)
    if package_dir.name != receipt.package_name:
        message = f"{package_dir.name} does not match its receipt identity"
        raise PackageEnvelopeError(message)
    for payload_name in receipt.payload_names:
        if not (package_dir / payload_name).is_file():
            message = f"{package_dir.name} is missing {payload_name}"
            raise FileNotFoundError(message)
    return ReportPackage(root=package_dir, receipt=receipt)


def parse_report_packages(
    roots: Mapping[Domain, Path | None],
    roster_ids: Collection[str],
) -> tuple[ParsedPackage, ...]:
    """Parse optional packages and reject cross-package wiring defects."""
    roster = set(roster_ids)
    seen: set[tuple[str, ...]] = set()
    parsed: list[ParsedPackage] = []
    for domain, root in roots.items():
        if root is None or not root.exists():
            continue
        if not root.is_dir():
            message = f"{domain} package root is not a directory"
            raise NotADirectoryError(message)
        for package_dir in sorted(path for path in root.iterdir() if path.is_dir()):
            try:
                package = load_report_package(package_dir, expected_domain=domain)
            except (
                ValidationError,
                UnicodeDecodeError,
                FileNotFoundError,
                PackageEnvelopeError,
            ) as error:
                parsed.append(
                    InvalidPackage(
                        domain=domain,
                        package_name=package_dir.name,
                        error=error,
                    ),
                )
                continue
            if package.receipt.sample_id not in roster:
                message = f"Unknown report package sample: {package.receipt.sample_id}"
                raise PackageEnvelopeError(message)
            if package.receipt.logical_key in seen:
                message = (
                    "Duplicate report package logical identity: "
                    f"{package.receipt.logical_key}"
                )
                raise PackageEnvelopeError(message)
            seen.add(package.receipt.logical_key)
            parsed.append(package)
    return tuple(
        sorted(
            parsed,
            key=lambda package: (
                DOMAIN_ORDER[package.domain]
                if isinstance(package, InvalidPackage)
                else DOMAIN_ORDER[Domain(package.receipt.domain)],
                package.package_name
                if isinstance(package, InvalidPackage)
                else package.receipt.package_name,
            ),
        ),
    )
