"""Smoke tests for the retained v3 CLI surface."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import TYPE_CHECKING

from typer.testing import CliRunner

from py_nvd.cli.app import app

if TYPE_CHECKING:
    import pytest


runner = CliRunner()


def test_retained_top_level_commands_show_help() -> None:
    for command in [
        "run",
        "params",
        "preset",
        "config",
        "taxonomy",
        "version",
        "setup",
    ]:
        result = runner.invoke(app, [command, "--help"])
        assert result.exit_code == 0, result.output


def test_removed_state_command_is_not_registered() -> None:
    result = runner.invoke(app, ["state", "--help"])
    assert result.exit_code != 0


def test_preset_lifecycle_uses_config_dir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("NVD_CONFIG_DIR", str(tmp_path))
    monkeypatch.delenv("NVD_PRESET_STORE", raising=False)

    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps({"samplesheet": "samples.csv", "dedup": True}))

    imported = runner.invoke(
        app,
        [
            "preset",
            "import",
            str(params_path),
            "--name",
            "smoke",
            "--description",
            "test",
        ],
    )
    assert imported.exit_code == 0, imported.output
    assert (tmp_path / "presets.sqlite").exists()
    assert not (tmp_path / "state.sqlite").exists()

    listed = runner.invoke(app, ["preset", "list"])
    assert listed.exit_code == 0, listed.output
    assert "smoke" in listed.output

    shown = runner.invoke(app, ["preset", "show", "smoke"])
    assert shown.exit_code == 0, shown.output
    assert "dedup" in shown.output

    export_path = tmp_path / "smoke.json"
    exported = runner.invoke(
        app,
        ["preset", "export", "smoke", "--output", str(export_path), "--format", "json"],
    )
    assert exported.exit_code == 0, exported.output
    assert '"dedup": true' in export_path.read_text()

    deleted = runner.invoke(app, ["preset", "delete", "smoke", "--force"])
    assert deleted.exit_code == 0, deleted.output


def test_preset_lifecycle_uses_explicit_preset_store(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_dir = tmp_path / "config"
    shared_store = tmp_path / "shared" / "presets.sqlite"
    monkeypatch.setenv("NVD_CONFIG_DIR", str(config_dir))
    monkeypatch.setenv("NVD_PRESET_STORE", str(shared_store))

    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps({"samplesheet": "samples.csv", "dedup": True}))

    imported = runner.invoke(
        app,
        [
            "preset",
            "import",
            str(params_path),
            "--name",
            "shared",
            "--description",
            "test",
        ],
    )
    assert imported.exit_code == 0, imported.output
    assert shared_store.exists()
    assert not (config_dir / "presets.sqlite").exists()
    assert not (tmp_path / "state.sqlite").exists()

    listed = runner.invoke(app, ["preset", "list"])
    assert listed.exit_code == 0, listed.output
    assert "shared" in listed.output


def test_taxonomy_status_uses_explicit_taxonomy_dir(tmp_path: Path) -> None:
    taxonomy_dir = tmp_path / "taxdump"
    result = runner.invoke(
        app,
        ["taxonomy", "status", "--taxonomy-dir", str(taxonomy_dir), "--json"],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["taxdump_dir"] == str(taxonomy_dir)


def test_run_uses_taxonomy_env_as_pipeline_param(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    samplesheet = tmp_path / "samples.csv"
    taxonomy_dir = Path("/shared/taxdump")
    samplesheet.write_text("sample_id,srr,platform,fastq1,fastq2\n", encoding="utf-8")
    monkeypatch.setenv("NVD_TAXONOMY_DB", str(taxonomy_dir))

    result = runner.invoke(
        app,
        ["run", "--samplesheet", str(samplesheet), "--dry-run"],
    )

    assert result.exit_code == 0, result.output
    assert "--taxonomy_dir" in result.output
    assert str(taxonomy_dir) in result.output


def test_run_accepts_skip_stage_flags(tmp_path: Path) -> None:
    """Hyphenated CLI flags map to underscore Nextflow params."""
    samplesheet = tmp_path / "samples.csv"
    samplesheet.write_text("sample_id,srr,platform,fastq1,fastq2\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "run",
            "--samplesheet",
            str(samplesheet),
            "--skip-assembly",
            "--skip-blast",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "--skip_assembly" in result.output
    assert "--skip_blast" in result.output
    assert "--skip-assembly" not in result.output
    assert "--skip-blast" not in result.output


def test_run_accepts_no_enrichment_flag(tmp_path: Path) -> None:
    """Hyphenated CLI flag maps to the underscore Nextflow param."""
    samplesheet = tmp_path / "samples.csv"
    samplesheet.write_text("sample_id,srr,platform,fastq1,fastq2\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "run",
            "--samplesheet",
            str(samplesheet),
            "--no-enrichment",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "--no_enrichment" in result.output
    assert "--no-enrichment" not in result.output


def test_run_accepts_low_complexity_read_filter_options(tmp_path: Path) -> None:
    """Read-complexity CLI options map to underscore Nextflow params."""
    samplesheet = tmp_path / "samples.csv"
    samplesheet.write_text("sample_id,srr,platform,fastq1,fastq2\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "run",
            "--samplesheet",
            str(samplesheet),
            "--filter-low-complexity-reads",
            "--min-read-entropy",
            "0.65",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "--filter_low_complexity_reads" in result.output
    assert "--min_read_entropy" in result.output
    assert "0.65" in result.output


def test_samplesheet_generate_sanitizes_illumina_ids(tmp_path: Path) -> None:
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    (fastq_dir / "patient-001_S7_L003_R1_001.fastq.gz").write_text("", encoding="utf-8")
    (fastq_dir / "patient-001_S7_L003_R2_001.fastq.gz").write_text("", encoding="utf-8")
    samplesheet = tmp_path / "samplesheet.csv"

    result = runner.invoke(
        app,
        [
            "samplesheet",
            "generate",
            "--from-dir",
            str(fastq_dir),
            "--platform",
            "illumina",
            "--output",
            str(samplesheet),
            "--force",
            "--sanitize",
        ],
    )

    assert result.exit_code == 0, result.output
    with samplesheet.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    assert reader.fieldnames == ["sample_id", "srr", "platform", "fastq1", "fastq2"]
    assert rows[0]["sample_id"] == "patient-001"


def test_samplesheet_generate_accepts_platform_aliases(tmp_path: Path) -> None:
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    (fastq_dir / "nanopore.fastq.gz").write_text("", encoding="utf-8")
    samplesheet = tmp_path / "samplesheet.csv"

    result = runner.invoke(
        app,
        [
            "samplesheet",
            "generate",
            "--from-dir",
            str(fastq_dir),
            "--platform",
            "nanopore",
            "--output",
            str(samplesheet),
            "--force",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "platform alias 'nanopore' was normalized to 'ont'" in result.output
    with samplesheet.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["platform"] == "ont"


def test_samplesheet_generate_groups_illumina_lanes_in_glob_columns(
    tmp_path: Path,
) -> None:
    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    for lane in ("L001", "L002"):
        (fastq_dir / f"patient-001_S7_{lane}_R1_001.fastq.gz").write_text(
            "",
            encoding="utf-8",
        )
        (fastq_dir / f"patient-001_S7_{lane}_R2_001.fastq.gz").write_text(
            "",
            encoding="utf-8",
        )
    samplesheet = tmp_path / "samplesheet.csv"

    result = runner.invoke(
        app,
        [
            "samplesheet",
            "generate",
            "--from-dir",
            str(fastq_dir),
            "--platform",
            "illumina",
            "--output",
            str(samplesheet),
            "--force",
            "--sanitize",
            "--group-by",
            "illumina-lanes",
        ],
    )

    assert result.exit_code == 0, result.output
    with samplesheet.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    assert reader.fieldnames == [
        "sample_id",
        "srr",
        "platform",
        "fastq1",
        "fastq2",
        "fastq1_glob",
        "fastq2_glob",
    ]
    assert rows == [
        {
            "sample_id": "patient-001",
            "srr": "",
            "platform": "illumina",
            "fastq1": "",
            "fastq2": "",
            "fastq1_glob": str(
                fastq_dir.absolute()
                / "patient-001_S7_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz",
            ),
            "fastq2_glob": str(
                fastq_dir.absolute()
                / "patient-001_S7_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9].fastq.gz",
            ),
        },
    ]


def test_samplesheet_validation_subcommands_use_same_read_preflight(
    tmp_path: Path,
) -> None:
    samplesheet = tmp_path / "samples.csv"
    samplesheet.write_text(
        "sample_id,srr,platform,fastq1,fastq2\nS1,,illumina,reads.fastq.gz,\n",
        encoding="utf-8",
    )

    for command in (
        ["samplesheet", "validate", str(samplesheet)],
        ["validate", "samplesheet", str(samplesheet)],
    ):
        result = runner.invoke(app, command)

        assert result.exit_code != 0
        assert "FASTQ paths and glob patterns must be absolute" in result.output
