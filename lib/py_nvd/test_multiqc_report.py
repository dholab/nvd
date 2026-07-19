"""Tests for NVD-owned MultiQC report input generation."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest
import yaml
from pydantic import ValidationError

from py_nvd.multiqc_packages import (
    MegablastQueryPartitionReceipt,
    PreparedQueryBatchesReceipt,
    write_package,
)
from py_nvd.multiqc_report import (
    CompileRequest,
    ReportConfiguration,
    ReportRoots,
    build_multiqc_inputs,
    parse_fastqc_packages,
    parse_roster,
)

if TYPE_CHECKING:
    from pathlib import Path


def write_roster(path: Path, *, shuffled: bool = False) -> Path:
    records = [
        {
            "sample_id": "sample_A",
            "platform": "illumina",
            "source": "paired_files",
            "read_structure": "paired",
            "read_counts": {"R1": 1, "R2": 1},
            "r1": ["/secret/source/A_R1.fastq.gz"],
            "r2": ["/secret/source/A_R2.fastq.gz"],
            "warnings": ["neutral warning"],
        },
        {
            "sample_id": "sample_B.label",
            "platform": "ont",
            "source": "single_file",
            "read_structure": "single",
            "read_counts": {"single": 1},
            "reads": ["/secret/private/B.fastq"],
        },
    ]
    if shuffled:
        records.reverse()
    path.write_text("\n".join(json.dumps(record) for record in records) + "\n")
    return path


def write_version(path: Path) -> Path:
    path.write_text(
        "version=3.3.0\nrevision=abc123\ndirty=true\ncommand=/must/not/copy\n",
        encoding="utf-8",
    )
    return path


def write_fastqc_package(
    root: Path,
    *,
    alias: str = "sample_A.raw.R1.001",
    package_alias: str | None = None,
    receipt_overrides: dict[str, object] | None = None,
) -> Path:
    package = root / f"{package_alias or alias}.fastqc"
    package.mkdir(parents=True)
    receipt = {
        "schema_version": "nvd.fastqc-raw-unit/v1",
        "sample_id": "sample_A",
        "platform": "illumina",
        "source": "paired_files",
        "read_end": "R1",
        "input_ordinal": 1,
        "alias": alias,
        "zip_alias": f"{alias}_fastqc.zip",
        "html_alias": f"{alias}_fastqc.html",
    }
    receipt.update(receipt_overrides or {})
    (package / "unit.json").write_text(json.dumps(receipt), encoding="utf-8")
    (package / receipt["zip_alias"]).write_text("zip", encoding="utf-8")
    (package / receipt["html_alias"]).write_text("html", encoding="utf-8")
    return package


def read_manifest(output_dir: Path) -> dict[str, object]:
    return json.loads((output_dir / "nvd_report_manifest.json").read_text())


def build_inputs(
    *,
    roster_path: Path,
    version_path: Path,
    fastqc_root: Path,
    output_dir: Path,
    experimental_enabled: bool,
) -> Path:
    return build_multiqc_inputs(
        CompileRequest(
            roster_path=roster_path,
            version_path=version_path,
            fastqc_root=fastqc_root,
            output_dir=output_dir,
            configuration=ReportConfiguration(
                experimental_enabled=experimental_enabled,
            ),
        ),
    )


def test_compiler_generates_minimal_manifest_and_custom_content(tmp_path: Path) -> None:
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")
    fastqc_root = tmp_path / "fastqc_packages"
    write_fastqc_package(fastqc_root)
    output_dir = tmp_path / "nvd_inputs"

    build_inputs(
        roster_path=roster,
        version_path=version,
        fastqc_root=fastqc_root,
        output_dir=output_dir,
        experimental_enabled=False,
    )

    manifest = read_manifest(output_dir)
    assert manifest == {
        "schema_version": "nvd.multiqc-report/v1",
        "report_plan": {
            "schema_version": "nvd.report-plan/v1",
            "experimental_enabled": False,
            "target_enrichment_enabled": True,
            "depletion_enabled": True,
            "assembly_enabled": True,
            "blast_enabled": True,
        },
        "source_identity": {"version": "3.3.0", "revision": "abc123"},
        "samples": [
            {
                "sample_id": "sample_A",
                "platform": "illumina",
                "source": "paired_files",
                "read_structure": "paired",
                "warnings": ["neutral warning"],
            },
            {
                "sample_id": "sample_B.label",
                "platform": "ont",
                "source": "single_file",
                "read_structure": "single",
                "warnings": [],
            },
        ],
        "raw_fastqc": [
            {
                "schema_version": "nvd.fastqc-raw-unit/v1",
                "sample_id": "sample_A",
                "platform": "illumina",
                "source": "paired_files",
                "read_end": "R1",
                "input_ordinal": 1,
                "alias": "sample_A.raw.R1.001",
                "zip_alias": "sample_A.raw.R1.001_fastqc.zip",
                "html_alias": "sample_A.raw.R1.001_fastqc.html",
            },
        ],
        "report_packages": [],
        "report_package_warnings": [],
        "fastx_max_points": 1000,
        "sections": [
            "nvd_sample_roster",
            "nvd_raw_fastqc_inventory",
            "nvd_experimental_capabilities",
        ],
    }
    assert not (output_dir / "nvd_source_identity_mqc.yaml").exists()
    assert yaml.safe_load((output_dir / "nvd_mqc_versions.yaml").read_text()) == {
        "NVD": "3.3.0 (abc123)",
    }
    assert (output_dir / "nvd_sample_roster_mqc.yaml").is_file()
    assert (output_dir / "nvd_raw_fastqc_inventory_mqc.yaml").is_file()
    assert (output_dir / "nvd_experimental_capabilities_mqc.yaml").is_file()

    roster_section = yaml.safe_load(
        (output_dir / "nvd_sample_roster_mqc.yaml").read_text(encoding="utf-8"),
    )
    assert all("sample_id" not in row for row in roster_section["data"].values())

    generated_text = "\n".join(path.read_text() for path in output_dir.iterdir())
    assert "/secret/source" not in generated_text
    assert "SRR_SECRET" not in generated_text
    assert "/must/not/copy" not in generated_text


def test_no_observed_fastqc_is_valid_and_invents_no_state(tmp_path: Path) -> None:
    output_dir = tmp_path / "nvd_inputs"

    fastqc_root = tmp_path / "fastqc_packages"
    fastqc_root.mkdir()

    build_inputs(
        roster_path=write_roster(tmp_path / "resolved_reads.jsonl"),
        version_path=write_version(tmp_path / "nvd_version.txt"),
        fastqc_root=fastqc_root,
        output_dir=output_dir,
        experimental_enabled=True,
    )

    manifest = read_manifest(output_dir)
    assert manifest["raw_fastqc"] == []
    assert "nvd_experimental_capabilities" not in manifest["sections"]
    assert not (output_dir / "nvd_raw_fastqc_inventory_mqc.yaml").exists()
    generated_text = "\n".join(path.read_text() for path in output_dir.iterdir())
    assert "missing" not in generated_text.lower()
    assert "failed" not in generated_text.lower()
    assert "skipped" not in generated_text.lower()


def test_compiler_renders_prepared_query_batches_and_localizes_invalid_payload(
    tmp_path: Path,
) -> None:
    fastqc_root = tmp_path / "fastqc_packages"
    fastqc_root.mkdir()
    package_root = tmp_path / "query_preparation_packages"
    package_root.mkdir()

    valid_summary = tmp_path / "sample_A.blast_query_batches.tsv"
    valid_summary.write_text(
        "sample_id\tplatform\tquery_class\tquery_source\tn_query_sequences\tquery_fasta_present\tquery_lookup_present\tproducer_note\n"
        "sample_A\tillumina\tshort_assembly_contig\tcontig\t2\ttrue\ttrue\tignored\n"
        "sample_A\tillumina\tlong_assembly_contig\tcontig\t0\ttrue\ttrue\tignored\n"
        "sample_A\tillumina\toverlap_merged_pair\tread_query\t0\tfalse\tfalse\tignored\n"
        "sample_A\tillumina\tsingle_read\tread_query\t0\tfalse\tfalse\tignored\n",
        encoding="utf-8",
    )
    valid_receipt = PreparedQueryBatchesReceipt(sample_id="sample_A")
    write_package(valid_receipt, (valid_summary,), output_root=package_root)

    invalid_summary = tmp_path / "sample_B.blast_query_batches.tsv"
    invalid_summary.write_text("wrong\nvalue\n", encoding="utf-8")
    invalid_receipt = PreparedQueryBatchesReceipt(sample_id="sample_B.label")
    write_package(invalid_receipt, (invalid_summary,), output_root=package_root)

    output_dir = build_multiqc_inputs(
        CompileRequest(
            roster_path=write_roster(tmp_path / "resolved_reads.jsonl"),
            version_path=write_version(tmp_path / "nvd_version.txt"),
            fastqc_root=fastqc_root,
            output_dir=tmp_path / "nvd_inputs",
            configuration=ReportConfiguration(experimental_enabled=False),
            report_roots=ReportRoots(query_preparation=package_root),
        ),
    )

    manifest = read_manifest(output_dir)
    assert {
        (receipt["domain"], receipt["artifact_type"], receipt["sample_id"])
        for receipt in manifest["report_packages"]
    } == {
        ("query_preparation", "prepared_query_batches", "sample_A"),
        ("query_preparation", "prepared_query_batches", "sample_B.label"),
    }
    assert "nvd_prepared_blast_query_batches" in manifest["sections"]

    section = yaml.safe_load(
        (output_dir / "nvd_prepared_blast_query_batches_mqc.yaml").read_text(
            encoding="utf-8",
        ),
    )
    assert section["section_name"] == "Prepared BLAST Query Batches"
    assert section["headers"]["platform"]["hidden"] is True
    assert section["headers"]["validation_details"]["hidden"] is True
    assert section["data"]["sample_A | short_assembly_contig"] == {
        "query_class": "Short assembly contig",
        "sequence_type": "Assembly contig",
        "availability": "enabled",
        "query_sequences": 2,
        "platform": "illumina",
    }
    assert section["data"]["sample_A | long_assembly_contig"]["query_sequences"] == 0
    assert section["data"]["sample_A | single_read"] == {
        "query_class": "Single read",
        "sequence_type": "Unassembled read",
        "availability": "disabled",
        "platform": "illumina",
    }
    assert section["data"]["sample_B.label"]["availability"] == "invalid"
    assert section["data"]["sample_B.label"]["problem"] == (
        "Prepared query summary could not be read"
    )
    assert section["data"]["sample_B.label"]["validation_details"]
    assert "SUMMARIZE_BLAST_QUERY_BATCHES" in section["description"]
    assert "producer_note" not in json.dumps(section)


def test_compiler_renders_blast_task_breakdown_and_localizes_invalid_payload(
    tmp_path: Path,
) -> None:
    fastqc_root = tmp_path / "fastqc_packages"
    fastqc_root.mkdir()
    package_root = tmp_path / "blast_packages"
    package_root.mkdir()

    valid_summary = tmp_path / "sample_A.megablast_query_partition.tsv"
    valid_summary.write_text(
        "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\tproducer_note\n"
        "sample_A\tshort_assembly_contig\t4\t2\t2\tignored\n",
        encoding="utf-8",
    )
    valid_receipt = MegablastQueryPartitionReceipt(
        sample_id="sample_A",
        query_class="short_assembly_contig",
    )
    write_package(valid_receipt, (valid_summary,), output_root=package_root)

    invalid_summary = tmp_path / "sample_B.megablast_query_partition.tsv"
    invalid_summary.write_text("wrong\nvalue\n", encoding="utf-8")
    invalid_receipt = MegablastQueryPartitionReceipt(
        sample_id="sample_B.label",
        query_class="long_assembly_contig",
    )
    write_package(invalid_receipt, (invalid_summary,), output_root=package_root)

    output_dir = build_multiqc_inputs(
        CompileRequest(
            roster_path=write_roster(tmp_path / "resolved_reads.jsonl"),
            version_path=write_version(tmp_path / "nvd_version.txt"),
            fastqc_root=fastqc_root,
            output_dir=tmp_path / "nvd_inputs",
            configuration=ReportConfiguration(
                experimental_enabled=False,
                blast_enabled=True,
            ),
            report_roots=ReportRoots(blast=package_root),
        ),
    )

    manifest = read_manifest(output_dir)
    assert {
        (
            receipt["domain"],
            receipt["artifact_type"],
            receipt["sample_id"],
            receipt["query_class"],
        )
        for receipt in manifest["report_packages"]
    } == {
        (
            "blast",
            "megablast_query_partition",
            "sample_A",
            "short_assembly_contig",
        ),
        (
            "blast",
            "megablast_query_partition",
            "sample_B.label",
            "long_assembly_contig",
        ),
    }
    assert "nvd_blast_task_breakdown" in manifest["sections"]

    section = yaml.safe_load(
        (output_dir / "nvd_blast_task_breakdown_mqc.yaml").read_text(
            encoding="utf-8",
        ),
    )
    assert section["section_name"] == "BLAST Task Breakdown"
    assert section["data"]["sample_A | short_assembly_contig"] == {
        "query_class": "Short assembly contig",
        "queries_searched": 4,
        "matched_by_megablast": 2,
        "forwarded_to_blastn": 2,
        "forwarded_to_blastn_percentage": 50.0,
    }
    assert section["data"]["sample_B.label | long_assembly_contig"]["problem"] == (
        "BLAST task breakdown could not be read"
    )
    assert section["headers"]["forwarded_to_blastn_percentage"]["title"] == (
        "Forwarded to BLASTN (%)"
    )
    assert section["headers"]["forwarded_to_blastn_percentage"]["format"] == ("{:,.1f}")
    assert section["headers"]["validation_details"]["hidden"] is True
    assert section["description"] == (
        "Queries are searched first with MEGABLAST. Queries without a MEGABLAST "
        "match are forwarded to BLASTN."
    )
    assert "producer_note" not in json.dumps(section)


def test_compiler_reports_blast_disabled_once_without_synthesizing_sample_rows(
    tmp_path: Path,
) -> None:
    fastqc_root = tmp_path / "fastqc_packages"
    fastqc_root.mkdir()

    output_dir = build_multiqc_inputs(
        CompileRequest(
            roster_path=write_roster(tmp_path / "resolved_reads.jsonl"),
            version_path=write_version(tmp_path / "nvd_version.txt"),
            fastqc_root=fastqc_root,
            output_dir=tmp_path / "nvd_inputs",
            configuration=ReportConfiguration(
                experimental_enabled=False,
                blast_enabled=False,
            ),
        ),
    )

    manifest = read_manifest(output_dir)
    assert manifest["report_plan"]["blast_enabled"] is False
    assert "nvd_blast_task_breakdown" in manifest["sections"]
    section = yaml.safe_load(
        (output_dir / "nvd_blast_task_breakdown_mqc.yaml").read_text(
            encoding="utf-8",
        ),
    )
    assert section["data"] == {
        "Run configuration": {
            "notice": "BLAST searching was disabled for this run.",
        },
    }
    assert "sample_A" not in section["data"]
    assert "sample_B.label" not in section["data"]


def test_compiler_rejects_duplicate_logical_key_unknown_sample_and_malformed_receipt(
    tmp_path: Path,
) -> None:
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")
    package = write_fastqc_package(tmp_path / "first")
    duplicate = write_fastqc_package(tmp_path / "second")

    with pytest.raises(ValueError, match="Duplicate FastQC logical key"):
        parse_fastqc_packages([package, duplicate], parse_roster(roster))

    unknown_root = tmp_path / "unknown_packages"
    unknown = write_fastqc_package(
        unknown_root,
        alias="unknown.raw.single.001",
        receipt_overrides={"sample_id": "unknown", "read_end": "single"},
    )
    with pytest.raises(ValueError, match="Unknown FastQC sample"):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=unknown.parent,
            output_dir=tmp_path / "unknown_out",
            experimental_enabled=True,
        )

    malformed_root = tmp_path / "malformed_packages"
    malformed = malformed_root / "bad.fastqc"
    malformed.mkdir(parents=True)
    (malformed / "unit.json").write_text("{}", encoding="utf-8")
    with pytest.raises(ValidationError):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=malformed_root,
            output_dir=tmp_path / "bad_out",
            experimental_enabled=True,
        )


def test_compiler_preserves_roster_order_and_sorts_fastqc_units(
    tmp_path: Path,
) -> None:
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")
    fastqc_root = tmp_path / "fastqc_packages"
    write_fastqc_package(
        fastqc_root,
        alias="sample_A.raw.R2.001",
        receipt_overrides={"read_end": "R2"},
    )
    write_fastqc_package(
        fastqc_root,
        alias="sample_A.raw.R1.001",
    )

    build_inputs(
        roster_path=roster,
        version_path=version,
        fastqc_root=fastqc_root,
        output_dir=tmp_path / "out_a",
        experimental_enabled=True,
    )
    shuffled_roster_path = tmp_path / "shuffled_resolved_reads.jsonl"
    shuffled_roster = write_roster(shuffled_roster_path, shuffled=True)
    build_inputs(
        roster_path=shuffled_roster,
        version_path=version,
        fastqc_root=fastqc_root,
        output_dir=tmp_path / "out_b",
        experimental_enabled=True,
    )

    manifest_a = read_manifest(tmp_path / "out_a")
    manifest_b = read_manifest(tmp_path / "out_b")

    assert [sample["sample_id"] for sample in manifest_a["samples"]] == [
        "sample_A",
        "sample_B.label",
    ]
    assert [sample["sample_id"] for sample in manifest_b["samples"]] == [
        "sample_B.label",
        "sample_A",
    ]
    assert [unit["alias"] for unit in manifest_a["raw_fastqc"]] == [
        "sample_A.raw.R1.001",
        "sample_A.raw.R2.001",
    ]
    assert [unit["alias"] for unit in manifest_b["raw_fastqc"]] == [
        "sample_A.raw.R1.001",
        "sample_A.raw.R2.001",
    ]


def test_compiler_rejects_unsafe_alias_identity(tmp_path: Path) -> None:
    package = tmp_path / "unsafe.fastqc"
    package.mkdir()
    (package / "unit.json").write_text(
        json.dumps(
            {
                "schema_version": "nvd.fastqc-raw-unit/v1",
                "sample_id": "sample_A",
                "platform": "illumina",
                "source": "paired_files",
                "read_end": "R1",
                "input_ordinal": 1,
                "alias": "sample_A/unsafe",
                "zip_alias": "sample_A/unsafe_fastqc.zip",
                "html_alias": "sample_A/unsafe_fastqc.html",
            },
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        build_inputs(
            roster_path=write_roster(tmp_path / "resolved_reads.jsonl"),
            version_path=write_version(tmp_path / "nvd_version.txt"),
            fastqc_root=tmp_path,
            output_dir=tmp_path / "out",
            experimental_enabled=True,
        )


def test_compiler_rejects_receipt_metadata_and_identity_mismatches(
    tmp_path: Path,
) -> None:
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")

    platform_root = tmp_path / "platform_mismatch"
    write_fastqc_package(platform_root, receipt_overrides={"platform": "ont"})
    with pytest.raises(ValueError, match="platform mismatch"):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=platform_root,
            output_dir=tmp_path / "platform_out",
            experimental_enabled=True,
        )

    source_root = tmp_path / "source_mismatch"
    write_fastqc_package(source_root, receipt_overrides={"source": "sra"})
    with pytest.raises(ValueError, match="source mismatch"):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=source_root,
            output_dir=tmp_path / "source_out",
            experimental_enabled=True,
        )

    package_root = tmp_path / "package_mismatch"
    write_fastqc_package(package_root, package_alias="sample_A.raw.R1.999")
    with pytest.raises(ValueError, match="directory does not match alias"):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=package_root,
            output_dir=tmp_path / "package_out",
            experimental_enabled=True,
        )


def test_compiler_rejects_invalid_read_end_ordinal_and_malformed_roster_json(
    tmp_path: Path,
) -> None:
    version = write_version(tmp_path / "nvd_version.txt")
    roster = write_roster(tmp_path / "resolved_reads.jsonl")

    read_end_root = tmp_path / "bad_read_end"
    write_fastqc_package(
        read_end_root,
        alias="sample_A.raw.bad.001",
        receipt_overrides={"read_end": "bad"},
    )
    with pytest.raises(ValidationError):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=read_end_root,
            output_dir=tmp_path / "read_end_out",
            experimental_enabled=True,
        )

    ordinal_root = tmp_path / "bad_ordinal"
    write_fastqc_package(
        ordinal_root,
        alias="sample_A.raw.R1.000",
        receipt_overrides={"input_ordinal": 0},
    )
    with pytest.raises(ValidationError):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=ordinal_root,
            output_dir=tmp_path / "ordinal_out",
            experimental_enabled=True,
        )

    malformed_roster = tmp_path / "malformed_roster.jsonl"
    malformed_roster.write_text('{"sample_id": "sample_A"\n', encoding="utf-8")
    empty_root = tmp_path / "empty_fastqc"
    empty_root.mkdir()
    with pytest.raises(ValidationError):
        build_inputs(
            roster_path=malformed_roster,
            version_path=version,
            fastqc_root=empty_root,
            output_dir=tmp_path / "malformed_roster_out",
            experimental_enabled=True,
        )


def test_compiler_rejects_roster_shape_edges(tmp_path: Path) -> None:
    version = write_version(tmp_path / "nvd_version.txt")
    empty_root = tmp_path / "empty_fastqc"
    empty_root.mkdir()

    empty_roster = tmp_path / "empty.jsonl"
    empty_roster.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="Roster is empty"):
        build_inputs(
            roster_path=empty_roster,
            version_path=version,
            fastqc_root=empty_root,
            output_dir=tmp_path / "empty_out",
            experimental_enabled=True,
        )

    duplicate_roster = tmp_path / "duplicate.jsonl"
    duplicate_record = {
        "sample_id": "sample_A",
        "platform": "illumina",
        "source": "paired_files",
        "read_structure": "paired",
        "read_counts": {"R1": 1, "R2": 1},
    }
    duplicate_roster.write_text(
        f"{json.dumps(duplicate_record)}\n{json.dumps(duplicate_record)}\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Duplicate roster sample_id"):
        build_inputs(
            roster_path=duplicate_roster,
            version_path=version,
            fastqc_root=empty_root,
            output_dir=tmp_path / "duplicate_out",
            experimental_enabled=True,
        )

    invalid_warnings = tmp_path / "invalid_warnings.jsonl"
    invalid_warnings.write_text(
        json.dumps(
            {
                "sample_id": "sample_A",
                "platform": "illumina",
                "source": "paired_files",
                "read_structure": "paired",
                "read_counts": {"R1": 1, "R2": 1},
                "warnings": [1],
            },
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValidationError):
        build_inputs(
            roster_path=invalid_warnings,
            version_path=version,
            fastqc_root=empty_root,
            output_dir=tmp_path / "warnings_out",
            experimental_enabled=True,
        )

    malformed_counts = tmp_path / "malformed_counts.jsonl"
    malformed_counts.write_text(
        json.dumps(
            {
                "sample_id": "sample_A",
                "platform": "illumina",
                "source": "paired_files",
                "read_structure": "paired",
                "read_counts": {"R1": 0, "R2": 1},
            },
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValidationError):
        build_inputs(
            roster_path=malformed_counts,
            version_path=version,
            fastqc_root=empty_root,
            output_dir=tmp_path / "shape_out",
            experimental_enabled=True,
        )


def test_compiler_rejects_fastqc_identity_edges(tmp_path: Path) -> None:
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")

    duplicate_first = write_fastqc_package(tmp_path / "duplicate_first")
    duplicate_second = write_fastqc_package(tmp_path / "duplicate_second")
    with pytest.raises(ValueError, match="Duplicate FastQC logical key"):
        parse_fastqc_packages(
            [duplicate_first, duplicate_second],
            parse_roster(roster),
        )

    bool_ordinal_root = tmp_path / "bool_ordinal"
    write_fastqc_package(
        bool_ordinal_root,
        receipt_overrides={"input_ordinal": True},
    )
    with pytest.raises(ValidationError):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=bool_ordinal_root,
            output_dir=tmp_path / "bool_ordinal_out",
            experimental_enabled=True,
        )

    missing_payload_root = tmp_path / "missing_payload"
    package = write_fastqc_package(missing_payload_root)
    (package / "sample_A.raw.R1.001_fastqc.html").unlink()
    with pytest.raises(ValueError, match="missing payload"):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=missing_payload_root,
            output_dir=tmp_path / "missing_payload_out",
            experimental_enabled=True,
        )

    single_roster = tmp_path / "single_resolved_reads.jsonl"
    single_roster.write_text(
        json.dumps(
            {
                "sample_id": "sample_A",
                "platform": "illumina",
                "source": "single_file",
                "read_structure": "single",
                "read_counts": {"single": 1},
                "reads": ["/private/single.fastq"],
            },
        )
        + "\n",
        encoding="utf-8",
    )
    incompatible_root = tmp_path / "incompatible_read_end"
    write_fastqc_package(
        incompatible_root,
        alias="sample_A.raw.R1.001",
        receipt_overrides={"source": "single_file"},
    )
    with pytest.raises(ValueError, match="incompatible read_end"):
        build_inputs(
            roster_path=single_roster,
            version_path=version,
            fastqc_root=incompatible_root,
            output_dir=tmp_path / "incompatible_out",
            experimental_enabled=True,
        )

    out_of_range_root = tmp_path / "out_of_range"
    write_fastqc_package(
        out_of_range_root,
        alias="sample_A.raw.R1.002",
        receipt_overrides={"input_ordinal": 2},
    )
    with pytest.raises(ValueError, match="input_ordinal out of range"):
        build_inputs(
            roster_path=roster,
            version_path=version,
            fastqc_root=out_of_range_root,
            output_dir=tmp_path / "out_of_range_out",
            experimental_enabled=True,
        )
