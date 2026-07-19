"""Build NVD-owned inputs for the terminal MultiQC report."""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from enum import StrEnum
from typing import TYPE_CHECKING, Annotated, Literal

import yaml
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PositiveInt,
    StringConstraints,
)

from py_nvd.multiqc_domains import (
    FASTX_MAX_POINTS,
    DomainConfiguration,
    LineGraphSection,
    ReportSection,
    TableSection,
    build_domain_sections,
    report_row_label,
)
from py_nvd.multiqc_packages import (
    Domain,
    InvalidPackage,
    ReportPackage,
    ReportReceipt,
    parse_report_packages,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from pathlib import Path

MANIFEST_SCHEMA = "nvd.multiqc-report/v1"
REPORT_PLAN_SCHEMA = "nvd.report-plan/v1"
FASTQC_RECEIPT_SCHEMA = "nvd.fastqc-raw-unit/v1"

SafeName = Annotated[
    str,
    StringConstraints(min_length=1, pattern=r"^[A-Za-z0-9_.-]+$"),
]
NonEmptyString = Annotated[str, StringConstraints(min_length=1)]


class NvdMultiqcInputError(ValueError):
    """Invalid NVD MultiQC compiler input."""


class FrozenModel(BaseModel):
    """Strict immutable model for parsed report inputs."""

    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)


class ReadSource(StrEnum):
    """Resolved read-input source variants."""

    SINGLE_FILE = "single_file"
    PAIRED_FILES = "paired_files"
    SINGLE_GLOB = "single_glob"
    PAIRED_GLOBS = "paired_globs"
    SRA = "sra"


class ReadEnd(StrEnum):
    """Physical read-unit positions understood by FastQC reporting."""

    SINGLE = "single"
    R1 = "R1"
    R2 = "R2"


class RosterSample(FrozenModel):
    """Report-facing projection of one resolved sample."""

    model_config = ConfigDict(extra="ignore", frozen=True, strict=True)

    sample_id: SafeName
    platform: NonEmptyString
    source: ReadSource
    read_structure: Literal["single", "paired", "unknown"]
    read_counts: dict[ReadEnd, PositiveInt] | None = Field(exclude=True)
    warnings: tuple[str, ...] = ()


class SourceIdentity(FrozenModel):
    """The redacted source identity retained in the report."""

    version: NonEmptyString
    revision: NonEmptyString


class FastqcReceipt(FrozenModel):
    """Self-describing receipt for one observed raw FastQC unit."""

    schema_version: Literal["nvd.fastqc-raw-unit/v1"]
    sample_id: SafeName
    platform: NonEmptyString
    source: ReadSource
    read_end: ReadEnd
    input_ordinal: PositiveInt
    alias: SafeName
    zip_alias: SafeName
    html_alias: SafeName

    @property
    def logical_key(self) -> tuple[str, ReadEnd, int]:
        return (self.sample_id, self.read_end, self.input_ordinal)


class ReportPlan(FrozenModel):
    """Run-level reporting configuration recorded in the manifest."""

    schema_version: Literal["nvd.report-plan/v1"] = REPORT_PLAN_SCHEMA
    experimental_enabled: bool
    target_enrichment_enabled: bool = True
    depletion_enabled: bool = True
    assembly_enabled: bool = True
    blast_enabled: bool = True


class PackageWarning(FrozenModel):
    domain: Domain
    package_name: str
    message: str


class ReportManifest(FrozenModel):
    """Machine-readable index of inputs presented to MultiQC."""

    schema_version: Literal["nvd.multiqc-report/v1"] = MANIFEST_SCHEMA
    report_plan: ReportPlan
    source_identity: SourceIdentity
    samples: tuple[RosterSample, ...]
    raw_fastqc: tuple[FastqcReceipt, ...]
    report_packages: tuple[ReportReceipt, ...] = ()
    report_package_warnings: tuple[PackageWarning, ...] = ()
    fastx_max_points: PositiveInt = FASTX_MAX_POINTS
    sections: tuple[SafeName, ...]


@dataclass(frozen=True, slots=True)
class ReportRoots:
    target_enrichment: Path | None = None
    depletion: Path | None = None
    fastx: Path | None = None
    assembly: Path | None = None
    query_preparation: Path | None = None
    blast: Path | None = None
    taxonomy: Path | None = None


@dataclass(frozen=True, slots=True)
class ReportConfiguration:
    experimental_enabled: bool
    target_enrichment_enabled: bool = True
    depletion_enabled: bool = True
    assembly_enabled: bool = True
    blast_enabled: bool = True


@dataclass(frozen=True, slots=True)
class CompileRequest:
    roster_path: Path
    version_path: Path
    fastqc_root: Path
    output_dir: Path
    configuration: ReportConfiguration
    report_roots: ReportRoots = field(default_factory=ReportRoots)


def build_multiqc_inputs(request: CompileRequest) -> Path:
    roster = parse_roster(request.roster_path)
    samples = tuple(roster.values())
    source_identity = parse_source_identity(request.version_path)
    raw_fastqc = parse_fastqc_packages(
        discover_fastqc_packages(request.fastqc_root),
        roster,
    )
    parsed_packages = parse_report_packages(
        {
            Domain.TARGET_ENRICHMENT: request.report_roots.target_enrichment,
            Domain.DEPLETION: request.report_roots.depletion,
            Domain.FASTX: request.report_roots.fastx,
            Domain.ASSEMBLY: request.report_roots.assembly,
            Domain.QUERY_PREPARATION: request.report_roots.query_preparation,
            Domain.BLAST: request.report_roots.blast,
            Domain.TAXONOMY: request.report_roots.taxonomy,
        },
        roster,
    )
    packages = tuple(
        package for package in parsed_packages if isinstance(package, ReportPackage)
    )
    invalid_packages = tuple(
        package for package in parsed_packages if isinstance(package, InvalidPackage)
    )
    domain_sections = build_domain_sections(
        packages=packages,
        sample_ids=tuple(roster),
        sample_platforms={
            sample_id: sample.platform for sample_id, sample in roster.items()
        },
        configuration=DomainConfiguration(
            experimental_enabled=request.configuration.experimental_enabled,
            target_enrichment_enabled=request.configuration.target_enrichment_enabled,
            depletion_enabled=request.configuration.depletion_enabled,
            assembly_enabled=request.configuration.assembly_enabled,
            blast_enabled=request.configuration.blast_enabled,
        ),
    )
    sections = ["nvd_sample_roster"]
    if raw_fastqc:
        sections.append("nvd_raw_fastqc_inventory")
    sections.extend(domain_sections)
    if not request.configuration.experimental_enabled:
        sections.append("nvd_experimental_capabilities")

    manifest = ReportManifest(
        report_plan=ReportPlan(
            experimental_enabled=request.configuration.experimental_enabled,
            target_enrichment_enabled=request.configuration.target_enrichment_enabled,
            depletion_enabled=request.configuration.depletion_enabled,
            assembly_enabled=request.configuration.assembly_enabled,
            blast_enabled=request.configuration.blast_enabled,
        ),
        source_identity=source_identity,
        samples=samples,
        raw_fastqc=raw_fastqc,
        report_packages=tuple(package.receipt for package in packages),
        report_package_warnings=tuple(
            PackageWarning(
                domain=package.domain,
                package_name=package.package_name,
                message=str(package.error),
            )
            for package in invalid_packages
        ),
        sections=tuple(sections),
    )

    if request.output_dir.exists():
        shutil.rmtree(request.output_dir)
    request.output_dir.mkdir(parents=True)
    (request.output_dir / "nvd_report_manifest.json").write_text(
        manifest.model_dump_json(indent=2) + "\n",
        encoding="utf-8",
    )
    write_software_versions(
        request.output_dir / "nvd_mqc_versions.yaml",
        source_identity,
    )
    write_roster(
        request.output_dir / "nvd_sample_roster_mqc.yaml",
        samples,
        raw_fastqc,
    )
    if raw_fastqc:
        write_fastqc_inventory(
            request.output_dir / "nvd_raw_fastqc_inventory_mqc.yaml",
            raw_fastqc,
        )
    write_domain_sections(request.output_dir, domain_sections)
    if not request.configuration.experimental_enabled:
        write_experimental_invitation(
            request.output_dir / "nvd_experimental_capabilities_mqc.yaml",
        )
    return request.output_dir


def parse_roster(path: Path) -> dict[str, RosterSample]:
    samples: dict[str, RosterSample] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        sample = RosterSample.model_validate_json(line, strict=True)
        if sample.sample_id in samples:
            message = f"Duplicate roster sample_id: {sample.sample_id}"
            raise NvdMultiqcInputError(message)
        samples[sample.sample_id] = sample
    if not samples:
        message = "Roster is empty"
        raise NvdMultiqcInputError(message)
    return samples


def discover_fastqc_packages(fastqc_root: Path) -> tuple[Path, ...]:
    if not fastqc_root.is_dir():
        message = f"FastQC root is not a directory: {fastqc_root}"
        raise NvdMultiqcInputError(message)
    return tuple(
        sorted(
            path
            for path in fastqc_root.iterdir()
            if path.is_dir() and path.name.endswith(".fastqc")
        ),
    )


def parse_source_identity(path: Path) -> SourceIdentity:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in {"version", "revision"}:
            values[key] = value
    return SourceIdentity.model_validate(values, strict=True)


def parse_fastqc_packages(
    packages: Sequence[Path],
    roster: Mapping[str, RosterSample],
) -> tuple[FastqcReceipt, ...]:
    sample_order = {sample_id: index for index, sample_id in enumerate(roster)}
    logical_keys: set[tuple[str, ReadEnd, int]] = set()
    observed: list[FastqcReceipt] = []
    for package in sorted(packages, key=lambda item: item.name):
        unit = parse_fastqc_package(package)
        if unit.logical_key in logical_keys:
            message = f"Duplicate FastQC logical key: {unit.logical_key}"
            raise NvdMultiqcInputError(message)
        logical_keys.add(unit.logical_key)
        sample = roster.get(unit.sample_id)
        if sample is None:
            message = f"Unknown FastQC sample: {unit.sample_id}"
            raise NvdMultiqcInputError(message)
        if unit.platform != sample.platform:
            message = f"FastQC receipt {unit.alias} platform mismatch for sample {unit.sample_id}"
            raise NvdMultiqcInputError(message)
        if unit.source != sample.source:
            message = f"FastQC receipt {unit.alias} source mismatch for sample {unit.sample_id}"
            raise NvdMultiqcInputError(message)
        require_observed_unit_bound(unit, sample)
        observed.append(unit)
    return tuple(
        sorted(
            observed,
            key=lambda unit: (
                sample_order[unit.sample_id],
                unit.read_end.value,
                unit.input_ordinal,
                unit.alias,
            ),
        ),
    )


def require_observed_unit_bound(
    unit: FastqcReceipt,
    sample: RosterSample,
) -> None:
    if sample.read_counts is None:
        return
    maximum = sample.read_counts.get(unit.read_end)
    if maximum is None:
        raise incompatible_read_end(unit)
    if unit.input_ordinal > maximum:
        message = (
            f"FastQC receipt {unit.alias} input_ordinal out of range "
            f"for sample {unit.sample_id}"
        )
        raise NvdMultiqcInputError(message)


def incompatible_read_end(unit: FastqcReceipt) -> NvdMultiqcInputError:
    message = (
        f"FastQC receipt {unit.alias} incompatible read_end for sample {unit.sample_id}"
    )
    return NvdMultiqcInputError(message)


def parse_fastqc_package(package: Path) -> FastqcReceipt:
    receipt_path = package / "unit.json"
    if not receipt_path.is_file():
        message = f"Malformed FastQC package {package.name}: missing unit.json"
        raise NvdMultiqcInputError(message)
    receipt = FastqcReceipt.model_validate_json(
        receipt_path.read_text(encoding="utf-8"),
        strict=True,
    )
    if package.name != f"{receipt.alias}.fastqc":
        message = (
            f"Malformed FastQC package {package.name}: directory does not match alias"
        )
        raise NvdMultiqcInputError(message)
    if (
        not (package / receipt.zip_alias).is_file()
        or not (package / receipt.html_alias).is_file()
    ):
        message = f"Malformed FastQC package {package.name}: missing payload"
        raise NvdMultiqcInputError(message)
    return receipt


def write_software_versions(path: Path, source_identity: SourceIdentity) -> None:
    path.write_text(
        yaml.safe_dump(
            {"NVD": f"{source_identity.version} ({source_identity.revision})"},
            sort_keys=False,
        ),
        encoding="utf-8",
    )


def write_roster(
    path: Path,
    samples: Sequence[RosterSample],
    raw_fastqc: Sequence[FastqcReceipt],
) -> None:
    read_ends_by_sample: dict[str, set[ReadEnd]] = {}
    for unit in raw_fastqc:
        read_ends_by_sample.setdefault(unit.sample_id, set()).add(unit.read_end)

    data: dict[str, dict[str, object]] = {}
    for sample in samples:
        row = sample.model_dump(mode="json", exclude={"sample_id"})
        read_ends = read_ends_by_sample.get(sample.sample_id, set())
        may_refine_sra_structure = (
            sample.source is ReadSource.SRA and sample.read_structure == "unknown"
        )
        if may_refine_sra_structure and {
            ReadEnd.R1,
            ReadEnd.R2,
        }.issubset(read_ends):
            row["read_structure"] = "paired"
        elif may_refine_sra_structure and read_ends == {ReadEnd.SINGLE}:
            row["read_structure"] = "single"
        data[sample.sample_id] = row

    path.write_text(
        yaml.safe_dump(
            {
                "id": "nvd_sample_roster",
                "section_name": "Sample Roster",
                "description": "Samples and read inputs for this run.",
                "plot_type": "table",
                "data": data,
                "headers": {
                    "read_structure": {
                        "title": "Read structure",
                        "scale": False,
                    },
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )


def write_fastqc_inventory(path: Path, raw_fastqc: Sequence[FastqcReceipt]) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "id": "nvd_raw_fastqc_inventory",
                "section_name": "Raw FastQC Inventory",
                "description": "Observed raw-read FastQC units that entered native MultiQC.",
                "plot_type": "table",
                "data": {
                    unit.alias: unit.model_dump(
                        mode="json",
                        exclude={"schema_version"},
                    )
                    for unit in raw_fastqc
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )


def write_domain_sections(
    output_dir: Path,
    sections: dict[str, ReportSection],
) -> None:
    filenames = {
        "nvd_higher_risk_taxonomic_findings": "nvd_higher_risk_taxonomic_findings_mqc.yaml",
        "nvd_strongest_taxonomic_signals": "nvd_strongest_taxonomic_signals_mqc.yaml",
        "nvd_faintest_taxonomic_signals": "nvd_faintest_taxonomic_signals_mqc.yaml",
        "nvd_target_enrichment": "nvd_target_enrichment_mqc.yaml",
        "nvd_depletion": "nvd_depletion_mqc.yaml",
        "nvd_fastx_profiles": "nvd_fastx_profiles_mqc.yaml",
        "nvd_fastx_single_read_length_distribution": "nvd_fastx_single_read_length_distribution_mqc.yaml",
        "nvd_fastx_overlap_merged_pair_length_distribution": "nvd_fastx_overlap_merged_pair_length_distribution_mqc.yaml",
        "nvd_fastx_filtered_contigs_length_distribution": "nvd_fastx_filtered_contigs_length_distribution_mqc.yaml",
        "nvd_fastx_quality_distribution": "nvd_fastx_quality_distribution_mqc.yaml",
        "nvd_assembly": "nvd_assembly_mqc.yaml",
        "nvd_prepared_blast_query_batches": "nvd_prepared_blast_query_batches_mqc.yaml",
        "nvd_blast_task_breakdown": "nvd_blast_task_breakdown_mqc.yaml",
    }
    for section_id, section in sections.items():
        content = {
            "id": section_id,
            "section_name": section.section_name,
            "description": section.description,
            "plot_type": section.plot_type,
        }
        if isinstance(section, TableSection):
            rows = {}
            for row in section.rows:
                label = report_row_label(row)
                if label in rows:
                    message = f"Duplicate report row label in {section_id}: {label}"
                    raise NvdMultiqcInputError(message)
                rows[label] = row.model_dump(
                    mode="json",
                    exclude={"sample_id"},
                    exclude_none=True,
                )
            content["data"] = rows
            if section.headers:
                content["headers"] = {
                    column: header.model_dump(mode="json", exclude_none=True)
                    for column, header in section.headers.items()
                }
            if section.first_column_title is not None:
                content["pconfig"] = {"col1_header": section.first_column_title}
        elif isinstance(section, LineGraphSection):
            content["data"] = section.data
            content["pconfig"] = {
                "xlab": section.xlab,
                "ylab": section.ylab,
                "style": "lines",
                "smooth_points": section.smooth_points,
                "xlog": False,
                "ylog": False,
                "logswitch": False,
            }
        (output_dir / filenames[section_id]).write_text(
            yaml.safe_dump(content, sort_keys=False),
            encoding="utf-8",
        )


def write_experimental_invitation(path: Path) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "id": "nvd_experimental_capabilities",
                "section_name": "Experimental Capabilities",
                "description": "Experimental mode can add CRUMBS profiling, rapid-screening evaluation, sample-similarity QC, and related results when enabled.",
                "plot_type": "table",
                "data": {"experimental_enabled": {"enabled": "false"}},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
