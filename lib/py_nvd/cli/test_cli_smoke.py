"""Smoke tests for CLI commands - verify they don't crash."""

import pytest
from typer.testing import CliRunner

from py_nvd.cli.app import app

runner = CliRunner()


class TestCLIHelp:
    """Verify help commands work."""

    def test_main_help(self):
        """--help exits cleanly."""
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "Usage:" in result.stdout

    def test_state_help(self):
        """state --help exits cleanly."""
        result = runner.invoke(app, ["state", "--help"])
        assert result.exit_code == 0
        assert "state" in result.stdout.lower()

    def test_validate_help(self):
        """validate --help exits cleanly."""
        result = runner.invoke(app, ["validate", "--help"])
        assert result.exit_code == 0

    def test_preset_help(self):
        """preset --help exits cleanly."""
        result = runner.invoke(app, ["preset", "--help"])
        assert result.exit_code == 0

    def test_params_help(self):
        """params --help exits cleanly."""
        result = runner.invoke(app, ["params", "--help"])
        assert result.exit_code == 0

    def test_config_help(self):
        """config --help exits cleanly."""
        result = runner.invoke(app, ["config", "--help"])
        assert result.exit_code == 0


class TestVersionCommand:
    """Verify version command works."""

    def test_version_runs(self):
        """version command outputs version info."""
        result = runner.invoke(app, ["version"])
        assert result.exit_code == 0
        # Should contain version info
        assert "NVD" in result.stdout or "Version" in result.stdout


class TestStateCommands:
    """Verify state commands work with isolated state directory."""

    def test_state_path(self, tmp_path, monkeypatch):
        """state path runs without crash."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        result = runner.invoke(app, ["state", "path"])
        assert result.exit_code == 0
        # Path may be wrapped across lines in output, so check for key parts
        assert "State Directory" in result.stdout
        assert "NVD_STATE_DIR" in result.stdout

    def test_state_path_json(self, tmp_path, monkeypatch):
        """state path --json outputs valid structure."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        result = runner.invoke(app, ["state", "path", "--json"])
        assert result.exit_code == 0
        assert "state_dir" in result.stdout

    def test_state_init(self, tmp_path, monkeypatch):
        """state init creates database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        result = runner.invoke(app, ["state", "init"])
        assert result.exit_code == 0
        assert (tmp_path / "state.sqlite").exists()

    def test_state_info_after_init(self, tmp_path, monkeypatch):
        """state info runs after init."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        # Initialize first
        runner.invoke(app, ["state", "init"])
        # Then check info
        result = runner.invoke(app, ["state", "info"])
        assert result.exit_code == 0
        assert "runs" in result.stdout.lower() or "Runs" in result.stdout

    def test_state_info_json(self, tmp_path, monkeypatch):
        """state info --json outputs valid structure."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["state", "info", "--json"])
        assert result.exit_code == 0
        assert "counts" in result.stdout

    def test_state_runs_empty(self, tmp_path, monkeypatch):
        """state runs handles empty database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["state", "runs"])
        assert result.exit_code == 0
        # Should indicate no runs or show empty table
        assert "No runs" in result.stdout or result.exit_code == 0

    def test_state_runs_json_empty(self, tmp_path, monkeypatch):
        """state runs --json outputs empty array for empty db."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["state", "runs", "--json"])
        assert result.exit_code == 0
        assert "[]" in result.stdout

    def test_state_samples_empty(self, tmp_path, monkeypatch):
        """state samples handles empty database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["state", "samples"])
        assert result.exit_code == 0

    def test_state_uploads_empty(self, tmp_path, monkeypatch):
        """state uploads handles empty database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["state", "uploads"])
        assert result.exit_code == 0

    def test_state_databases_empty(self, tmp_path, monkeypatch):
        """state databases handles empty database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["state", "databases"])
        assert result.exit_code == 0

    def test_state_taxonomy(self, tmp_path, monkeypatch):
        """state taxonomy runs without crash."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        result = runner.invoke(app, ["state", "taxonomy"])
        assert result.exit_code == 0

    def test_state_info_includes_hits(self, tmp_path, monkeypatch):
        """state info shows hit counts."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGT", "2024-01-01T00:00:00Z", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["state", "info", "--json"])
        assert result.exit_code == 0
        assert '"hits": 1' in result.stdout
        assert '"hit_observations": 1' in result.stdout

    def test_state_prune_cascades_to_hit_observations(self, tmp_path, monkeypatch):
        """Pruning a run deletes its hit observations."""
        from datetime import datetime, timedelta

        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import (
            count_hit_observations,
            record_hit_observation,
            register_hit,
        )
        from py_nvd.state import complete_run, register_run

        # Create an old run and mark it completed
        old_date = (datetime.now() - timedelta(days=100)).strftime("%Y-%m-%d %H:%M:%S")
        register_run("old_run", "old_set", state_dir=tmp_path)
        complete_run("old_run", "completed", state_dir=tmp_path)

        # Manually set the started_at to be old
        from py_nvd.db import connect

        with connect(tmp_path) as conn:
            conn.execute(
                "UPDATE runs SET started_at = ? WHERE run_id = ?", (old_date, "old_run")
            )
            conn.commit()

        # Register hit and observation for this run
        hit, _ = register_hit("ACGTACGT", old_date, tmp_path)
        record_hit_observation(
            hit.hit_key, "old_set", "sample_a", old_date, state_dir=tmp_path
        )

        assert count_hit_observations(tmp_path) == 1

        # Prune runs older than 90 days
        result = runner.invoke(app, ["state", "prune", "--older-than", "90d", "--yes"])
        assert result.exit_code == 0
        assert "hit observation" in result.stdout

        # Observation should be deleted
        assert count_hit_observations(tmp_path) == 0

    def test_state_prune_garbage_collects_orphaned_hits(self, tmp_path, monkeypatch):
        """Hits with no observations are deleted after prune."""
        from datetime import datetime, timedelta

        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import count_hits, record_hit_observation, register_hit
        from py_nvd.state import complete_run, register_run

        old_date = (datetime.now() - timedelta(days=100)).strftime("%Y-%m-%d %H:%M:%S")
        register_run("old_run", "old_set", state_dir=tmp_path)
        complete_run("old_run", "completed", state_dir=tmp_path)

        from py_nvd.db import connect

        with connect(tmp_path) as conn:
            conn.execute(
                "UPDATE runs SET started_at = ? WHERE run_id = ?", (old_date, "old_run")
            )
            conn.commit()

        # Register hit with single observation
        hit, _ = register_hit("ACGTACGT", old_date, tmp_path)
        record_hit_observation(
            hit.hit_key, "old_set", "sample_a", old_date, state_dir=tmp_path
        )

        assert count_hits(tmp_path) == 1

        # Prune
        result = runner.invoke(app, ["state", "prune", "--older-than", "90d", "--yes"])
        assert result.exit_code == 0
        assert "orphaned hit" in result.stdout

        # Hit should be garbage collected
        assert count_hits(tmp_path) == 0

    def test_state_prune_preserves_shared_hits(self, tmp_path, monkeypatch):
        """Hits observed in multiple runs survive partial prune."""
        from datetime import datetime, timedelta

        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import count_hits, record_hit_observation, register_hit
        from py_nvd.state import complete_run, register_run

        old_date = (datetime.now() - timedelta(days=100)).strftime("%Y-%m-%d %H:%M:%S")
        new_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Create old run and new run, mark both completed
        register_run("old_run", "old_set", state_dir=tmp_path)
        complete_run("old_run", "completed", state_dir=tmp_path)
        register_run("new_run", "new_set", state_dir=tmp_path)
        complete_run("new_run", "completed", state_dir=tmp_path)

        from py_nvd.db import connect

        with connect(tmp_path) as conn:
            conn.execute(
                "UPDATE runs SET started_at = ? WHERE run_id = ?", (old_date, "old_run")
            )
            conn.commit()

        # Same hit observed in both runs
        hit, _ = register_hit("ACGTACGT", old_date, tmp_path)
        record_hit_observation(
            hit.hit_key, "old_set", "sample_a", old_date, state_dir=tmp_path
        )
        record_hit_observation(
            hit.hit_key, "new_set", "sample_b", new_date, state_dir=tmp_path
        )

        assert count_hits(tmp_path) == 1

        # Prune old run only
        result = runner.invoke(app, ["state", "prune", "--older-than", "90d", "--yes"])
        assert result.exit_code == 0

        # Hit should still exist (has observation from new run)
        assert count_hits(tmp_path) == 1

        # Prune old run only
        result = runner.invoke(app, ["state", "prune", "--older-than", "90d", "--yes"])
        assert result.exit_code == 0

        # Hit should still exist (has observation from new run)
        assert count_hits(tmp_path) == 1


class TestValidateCommands:
    """Verify validate commands run without crashing."""

    def test_validate_params(self):
        """validate params runs (may warn but shouldn't crash)."""
        result = runner.invoke(app, ["validate", "params"])
        # Exit code 0 or 1 acceptable (validation pass/fail)
        # Exit code 2 would indicate crash
        assert result.exit_code in (0, 1)

    def test_validate_deps(self):
        """validate deps runs (may fail if deps missing but shouldn't crash)."""
        result = runner.invoke(app, ["validate", "deps"])
        # May exit 1 if deps missing, but shouldn't crash
        assert result.exit_code in (0, 1)


class TestPresetCommands:
    """Verify preset commands work with isolated state directory."""

    def test_preset_list_empty(self, tmp_path, monkeypatch):
        """preset list handles empty database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        # Initialize state first
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["preset", "list"])
        assert result.exit_code == 0

    def test_preset_register_with_inline_params(self, tmp_path, monkeypatch):
        """preset register with inline params creates a preset."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        # Register a preset with inline parameters (valid schema params)
        result = runner.invoke(
            app,
            [
                "preset",
                "register",
                "test-preset",
                "--cutoff-percent",
                "0.01",
                "--description",
                "A test preset",
            ],
        )
        assert result.exit_code == 0

        # Verify it appears in list
        result = runner.invoke(app, ["preset", "list"])
        assert result.exit_code == 0
        assert "test-preset" in result.stdout


class TestRunCommandHelp:
    """Verify run command help works (don't actually run pipeline)."""

    def test_run_help(self):
        """run --help shows usage."""
        result = runner.invoke(app, ["run", "--help"])
        assert result.exit_code == 0
        assert "samplesheet" in result.stdout.lower()


class TestSamplesheetCommands:
    """Verify samplesheet commands work."""

    def test_samplesheet_help(self):
        """samplesheet --help exits cleanly."""
        result = runner.invoke(app, ["samplesheet", "--help"])
        assert result.exit_code == 0
        assert "generate" in result.stdout.lower()
        assert "validate" in result.stdout.lower()

    def test_samplesheet_generate_help(self):
        """samplesheet generate --help shows options."""
        result = runner.invoke(app, ["samplesheet", "generate", "--help"])
        assert result.exit_code == 0
        assert "--from-dir" in result.stdout
        assert "--from-sra" in result.stdout
        assert "--platform" in result.stdout
        assert "--dry-run" in result.stdout

    def test_samplesheet_validate_help(self):
        """samplesheet validate --help shows usage."""
        result = runner.invoke(app, ["samplesheet", "validate", "--help"])
        assert result.exit_code == 0

    def test_samplesheet_generate_requires_input(self):
        """samplesheet generate fails without --from-dir or --from-sra."""
        result = runner.invoke(
            app, ["samplesheet", "generate", "--platform", "illumina"]
        )
        assert result.exit_code == 1
        assert "Must specify" in result.stdout

    def test_samplesheet_generate_dry_run(self, tmp_path):
        """samplesheet generate --dry-run works with FASTQ directory."""
        # Create test FASTQ files
        (tmp_path / "sample1_R1.fastq.gz").touch()
        (tmp_path / "sample1_R2.fastq.gz").touch()

        result = runner.invoke(
            app,
            [
                "samplesheet",
                "generate",
                "--from-dir",
                str(tmp_path),
                "--platform",
                "illumina",
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "sample1" in result.stdout
        assert "Dry run" in result.stdout

    def test_samplesheet_generate_writes_file(self, tmp_path):
        """samplesheet generate creates valid CSV file."""
        # Create test FASTQ files
        fastq_dir = tmp_path / "fastqs"
        fastq_dir.mkdir()
        (fastq_dir / "sample1_R1.fastq.gz").touch()
        (fastq_dir / "sample1_R2.fastq.gz").touch()

        output_file = tmp_path / "samplesheet.csv"

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
                str(output_file),
                "--force",
            ],
        )
        assert result.exit_code == 0
        assert output_file.exists()
        assert "Wrote samplesheet" in result.stdout

        # Verify CSV content
        content = output_file.read_text()
        assert "sample_id" in content
        assert "sample1" in content
        assert "illumina" in content

    def test_samplesheet_generate_from_sra(self, tmp_path):
        """samplesheet generate --from-sra works with accession list."""
        # Create test accession file
        accession_file = tmp_path / "accessions.txt"
        accession_file.write_text("SRR123456\nSRR789012\n")

        result = runner.invoke(
            app,
            [
                "samplesheet",
                "generate",
                "--from-sra",
                str(accession_file),
                "--platform",
                "ont",
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "SRR123456" in result.stdout
        assert "SRR789012" in result.stdout

    def test_samplesheet_validate_valid_file(self, tmp_path):
        """samplesheet validate passes for valid samplesheet."""
        # Create valid samplesheet
        samplesheet = tmp_path / "samples.csv"
        samplesheet.write_text(
            "sample_id,srr,platform,fastq1,fastq2\nsample1,SRR123456,illumina,,\n"
        )

        result = runner.invoke(app, ["samplesheet", "validate", str(samplesheet)])
        assert result.exit_code == 0
        assert "valid" in result.stdout.lower()

    def test_samplesheet_validate_invalid_file(self, tmp_path):
        """samplesheet validate fails for invalid samplesheet."""
        # Create invalid samplesheet (missing required column)
        samplesheet = tmp_path / "samples.csv"
        samplesheet.write_text("sample_id,platform\nsample1,illumina\n")

        result = runner.invoke(app, ["samplesheet", "validate", str(samplesheet)])
        assert result.exit_code == 1
        assert "Missing required columns" in result.stdout

    def test_samplesheet_validate_missing_file(self, tmp_path):
        """samplesheet validate fails for missing file."""
        result = runner.invoke(
            app, ["samplesheet", "validate", str(tmp_path / "nonexistent.csv")]
        )
        assert result.exit_code == 1

    def test_samplesheet_alias_ss(self):
        """ss alias works for samplesheet."""
        result = runner.invoke(app, ["ss", "--help"])
        assert result.exit_code == 0
        assert "generate" in result.stdout.lower()

    def test_samplesheet_alias_sheet(self):
        """sheet alias works for samplesheet."""
        result = runner.invoke(app, ["sheet", "--help"])
        assert result.exit_code == 0
        assert "generate" in result.stdout.lower()

    def test_samplesheet_generate_warns_previously_processed(
        self, tmp_path, monkeypatch
    ):
        """samplesheet generate warns about previously processed samples."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))

        # Initialize state and register a previous run with sample1
        from py_nvd.state import (
            compute_sample_set_id,
            record_upload,
            register_processed_sample,
            register_run,
        )

        sample_ids = ["sample1", "other_sample"]
        sample_set_id = compute_sample_set_id(sample_ids)
        register_run("previous_run", sample_set_id)
        register_processed_sample("sample1", sample_set_id, "previous_run")
        record_upload("sample1", sample_set_id, "blast", "labkey", "hash123")

        # Create test FASTQ files including sample1
        fastq_dir = tmp_path / "fastqs"
        fastq_dir.mkdir()
        (fastq_dir / "sample1_R1.fastq.gz").touch()
        (fastq_dir / "sample1_R2.fastq.gz").touch()

        # State check happens automatically - no flag needed
        result = runner.invoke(
            app,
            [
                "samplesheet",
                "generate",
                "--from-dir",
                str(fastq_dir),
                "--platform",
                "illumina",
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "previously processed" in result.stdout
        assert "sample1" in result.stdout
        assert "previous_run" in result.stdout
        assert "uploaded" in result.stdout


class TestSecretsCommands:
    """Verify secrets commands work (requires nextflow in PATH)."""

    def test_secrets_help(self):
        """secrets --help exits cleanly."""
        result = runner.invoke(app, ["secrets", "--help"])
        assert result.exit_code == 0
        assert "list" in result.stdout.lower()
        assert "set" in result.stdout.lower()

    def test_secrets_list_help(self):
        """secrets list --help shows options."""
        result = runner.invoke(app, ["secrets", "list", "--help"])
        assert result.exit_code == 0
        assert "--show-values" in result.stdout

    def test_secrets_set_help(self):
        """secrets set --help shows usage."""
        result = runner.invoke(app, ["secrets", "set", "--help"])
        assert result.exit_code == 0
        assert "LABKEY_API_KEY" in result.stdout

    def test_secrets_check_help(self):
        """secrets check --help shows usage."""
        result = runner.invoke(app, ["secrets", "check", "--help"])
        assert result.exit_code == 0


class TestResumeCommand:
    """Verify resume command works."""

    def test_resume_help(self):
        """resume --help exits cleanly."""
        result = runner.invoke(app, ["resume", "--help"])
        assert result.exit_code == 0
        assert "--interactive" in result.stdout
        assert "-i" in result.stdout
        assert ".nfresume" in result.stdout

    def test_resume_no_file(self, tmp_path, monkeypatch):
        """resume fails gracefully when no .nfresume exists."""
        monkeypatch.chdir(tmp_path)
        result = runner.invoke(app, ["resume"])
        assert result.exit_code == 1
        assert "No previous run detected" in result.stdout
        assert ".nfresume not found" in result.stdout

    def test_resume_interactive_no_file(self, tmp_path, monkeypatch):
        """resume -i fails gracefully when no .nfresume exists."""
        monkeypatch.chdir(tmp_path)
        result = runner.invoke(app, ["resume", "-i"])
        assert result.exit_code == 1
        assert "No previous run detected" in result.stdout


class TestPipelineRoot:
    """Verify PIPELINE_ROOT detection works correctly."""

    def test_pipeline_root_exists(self):
        """PIPELINE_ROOT points to a valid pipeline directory."""
        from py_nvd.cli.utils import PIPELINE_ROOT

        assert PIPELINE_ROOT.exists()
        assert (PIPELINE_ROOT / "main.nf").exists()
        assert (PIPELINE_ROOT / "nextflow.config").exists()

    def test_pipeline_root_is_absolute(self):
        """PIPELINE_ROOT is an absolute path."""
        from py_nvd.cli.utils import PIPELINE_ROOT

        assert PIPELINE_ROOT.is_absolute()


class TestHitsCommands:
    """Verify hits commands work with isolated state directory."""

    def test_hits_help(self):
        """hits --help exits cleanly."""
        result = runner.invoke(app, ["hits", "--help"])
        assert result.exit_code == 0
        assert "export" in result.stdout.lower()

    def test_hits_export_help(self):
        """hits export --help shows options."""
        result = runner.invoke(app, ["hits", "export", "--help"])
        assert result.exit_code == 0
        assert "--format" in result.stdout
        assert "--output" in result.stdout
        assert "--include-sequence" in result.stdout

    def test_hits_export_empty_database(self, tmp_path, monkeypatch):
        """hits export handles empty database gracefully."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["hits", "export"])
        assert result.exit_code == 0
        assert "No hits found" in result.stdout

    def test_hits_export_tsv(self, tmp_path, monkeypatch):
        """hits export outputs TSV format."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        # Register a hit
        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01T00:00:00Z", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            contig_id="NODE_1",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "export", "--format", "tsv"])
        assert result.exit_code == 0
        assert "hit_key" in result.stdout
        assert "sample_a" in result.stdout
        assert hit.hit_key in result.stdout

    def test_hits_export_csv(self, tmp_path, monkeypatch):
        """hits export outputs CSV format."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01T00:00:00Z", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "export", "--format", "csv"])
        assert result.exit_code == 0
        assert "hit_key" in result.stdout

    def test_hits_export_fasta(self, tmp_path, monkeypatch):
        """hits export outputs FASTA format with just hit_key header."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01T00:00:00Z", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "export", "--format", "fasta"])
        assert result.exit_code == 0
        assert f">{hit.hit_key}" in result.stdout
        assert "ACGTACGTACGT" in result.stdout

    def test_hits_export_parquet_requires_output(self, tmp_path, monkeypatch):
        """hits export --format parquet requires --output."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        result = runner.invoke(app, ["hits", "export", "--format", "parquet"])
        assert result.exit_code == 1
        assert "requires --output" in result.stdout

    def test_hits_export_parquet_to_file(self, tmp_path, monkeypatch):
        """hits export --format parquet writes to file."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01T00:00:00Z", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=tmp_path,
        )

        output_file = tmp_path / "hits.parquet"
        result = runner.invoke(
            app, ["hits", "export", "--format", "parquet", "--output", str(output_file)]
        )
        assert result.exit_code == 0
        assert output_file.exists()

    def test_hits_export_include_sequence(self, tmp_path, monkeypatch):
        """hits export --include-sequence includes full sequence."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01T00:00:00Z", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "export", "--include-sequence"])
        assert result.exit_code == 0
        assert "sequence" in result.stdout
        assert "ACGTACGTACGT" in result.stdout

    def test_hits_export_to_file(self, tmp_path, monkeypatch):
        """hits export --output writes to file."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01T00:00:00Z", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=tmp_path,
        )

        output_file = tmp_path / "hits.tsv"
        result = runner.invoke(app, ["hits", "export", "--output", str(output_file)])
        assert result.exit_code == 0
        assert output_file.exists()
        assert "Wrote" in result.stdout

        content = output_file.read_text()
        assert "hit_key" in content
        assert hit.hit_key in content

    def test_hits_stats_help(self):
        """hits stats --help shows options."""
        result = runner.invoke(app, ["hits", "stats", "--help"])
        assert result.exit_code == 0
        assert "--json" in result.stdout
        assert "--top" in result.stdout

    def test_hits_stats_empty_database(self, tmp_path, monkeypatch):
        """hits stats handles empty database gracefully."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["hits", "stats"])
        assert result.exit_code == 0
        assert "No hits in database" in result.stdout

    def test_hits_stats_empty_database_json(self, tmp_path, monkeypatch):
        """hits stats --json handles empty database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])
        result = runner.invoke(app, ["hits", "stats", "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        assert data["total_hits"] == 0
        assert data["top_recurring"] == []

    def test_hits_stats_with_data(self, tmp_path, monkeypatch):
        """hits stats shows statistics for hits."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "stats"])
        assert result.exit_code == 0
        assert "Total unique hits:" in result.stdout
        assert "1" in result.stdout

    def test_hits_stats_json_with_data(self, tmp_path, monkeypatch):
        """hits stats --json outputs valid JSON with data."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "stats", "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        assert data["total_hits"] == 1
        assert data["total_observations"] == 1
        assert data["unique_samples"] == 1

    def test_hits_stats_top_recurring(self, tmp_path, monkeypatch):
        """hits stats shows recurring hits."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        # Same hit in two samples
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "stats"])
        assert result.exit_code == 0
        assert hit.hit_key[:16] in result.stdout  # Truncated key shown

    def test_hits_stats_no_recurring(self, tmp_path, monkeypatch):
        """hits stats shows message when no recurring hits."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        # Only one sample
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "stats"])
        assert result.exit_code == 0
        assert "No recurring hits found" in result.stdout

    def test_hits_lookup_help(self):
        """hits lookup --help shows options."""
        result = runner.invoke(app, ["hits", "lookup", "--help"])
        assert result.exit_code == 0
        assert "--json" in result.stdout
        assert "QUERY" in result.stdout

    def test_hits_lookup_not_found(self, tmp_path, monkeypatch):
        """hits lookup shows message when not found."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        result = runner.invoke(app, ["hits", "lookup", "ACGTACGT"])
        assert result.exit_code == 0
        assert "not found" in result.stdout.lower()

    def test_hits_lookup_not_found_json(self, tmp_path, monkeypatch):
        """hits lookup --json shows found=false when not found."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        result = runner.invoke(app, ["hits", "lookup", "ACGTACGT", "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        assert data["found"] is False

    def test_hits_lookup_by_sequence(self, tmp_path, monkeypatch):
        """hits lookup finds hit by sequence."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "lookup", seq])
        assert result.exit_code == 0
        assert "Hit found" in result.stdout
        assert hit.hit_key in result.stdout

    def test_hits_lookup_by_key(self, tmp_path, monkeypatch):
        """hits lookup finds hit by key."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "lookup", hit.hit_key])
        assert result.exit_code == 0
        assert "Hit found" in result.stdout

    def test_hits_lookup_json(self, tmp_path, monkeypatch):
        """hits lookup --json outputs valid JSON."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "lookup", hit.hit_key, "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        assert data["found"] is True
        assert data["hit"]["hit_key"] == hit.hit_key
        assert len(data["observations"]) == 1

    def test_hits_lookup_shows_observations(self, tmp_path, monkeypatch):
        """hits lookup shows all observations."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-02",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "lookup", hit.hit_key])
        assert result.exit_code == 0
        assert "sample_a" in result.stdout
        assert "sample_b" in result.stdout
        assert "2 total" in result.stdout

    def test_hits_trace_help(self):
        """hits trace --help shows options."""
        result = runner.invoke(app, ["hits", "trace", "--help"])
        assert result.exit_code == 0
        assert "--json" in result.stdout
        assert "HIT_KEY" in result.stdout

    def test_hits_trace_invalid_key(self, tmp_path, monkeypatch):
        """hits trace rejects invalid key format."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        result = runner.invoke(app, ["hits", "trace", "not-a-valid-key"])
        assert result.exit_code == 1
        assert "Invalid hit key" in result.stdout

    def test_hits_trace_not_found(self, tmp_path, monkeypatch):
        """hits trace shows message when not found."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        # Valid format but doesn't exist
        fake_key = "a" * 32
        result = runner.invoke(app, ["hits", "trace", fake_key])
        assert result.exit_code == 0
        assert "not found" in result.stdout.lower()

    def test_hits_trace_not_found_json(self, tmp_path, monkeypatch):
        """hits trace --json shows found=false when not found."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        fake_key = "b" * 32
        result = runner.invoke(app, ["hits", "trace", fake_key, "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        assert data["found"] is False

    def test_hits_trace_shows_sequence(self, tmp_path, monkeypatch):
        """hits trace shows full sequence."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "trace", hit.hit_key])
        assert result.exit_code == 0
        assert "Hit found" in result.stdout
        assert "Sequence:" in result.stdout
        assert seq in result.stdout

    def test_hits_trace_json_includes_sequence(self, tmp_path, monkeypatch):
        """hits trace --json includes full sequence."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "trace", hit.hit_key, "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        assert data["found"] is True
        assert data["sequence"] == seq
        assert "sequence_compressed" not in data["hit"]

    def test_hits_trace_shows_observations(self, tmp_path, monkeypatch):
        """hits trace shows all observations."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-02",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "trace", hit.hit_key])
        assert result.exit_code == 0
        assert "sample_a" in result.stdout
        assert "sample_b" in result.stdout
        assert "2 total" in result.stdout

    def test_hits_timeline_help(self):
        """hits timeline --help shows options."""
        result = runner.invoke(app, ["hits", "timeline", "--help"])
        assert result.exit_code == 0
        assert "--granularity" in result.stdout
        assert "--json" in result.stdout

    def test_hits_timeline_empty(self, tmp_path, monkeypatch):
        """hits timeline handles empty database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        result = runner.invoke(app, ["hits", "timeline"])
        assert result.exit_code == 0
        assert "No hits in database" in result.stdout

    def test_hits_timeline_invalid_granularity(self, tmp_path, monkeypatch):
        """hits timeline rejects invalid granularity."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        result = runner.invoke(app, ["hits", "timeline", "--granularity", "hourly"])
        assert result.exit_code == 1
        assert "Invalid granularity" in result.stdout

    def test_hits_timeline_shows_histogram(self, tmp_path, monkeypatch):
        """hits timeline shows ASCII histogram."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-15", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-15",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "timeline"])
        assert result.exit_code == 0
        assert "Discovery Timeline" in result.stdout
        assert "2024-01" in result.stdout
        assert "█" in result.stdout
        assert "1 unique hits discovered" in result.stdout

    def test_hits_timeline_fills_gaps(self, tmp_path, monkeypatch):
        """hits timeline fills gaps between periods."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        # Hit in January
        hit1, _ = register_hit("ACGTACGTACGT", "2024-01-15", tmp_path)
        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_a",
            "2024-01-15",
            state_dir=tmp_path,
        )
        # Hit in March (skip February)
        hit2, _ = register_hit("GGGGCCCCAAAA", "2024-03-15", tmp_path)
        record_hit_observation(
            hit2.hit_key,
            "set_001",
            "sample_b",
            "2024-03-15",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "timeline"])
        assert result.exit_code == 0
        # All three months should appear
        assert "2024-01" in result.stdout
        assert "2024-02" in result.stdout
        assert "2024-03" in result.stdout
        assert "2 unique hits discovered" in result.stdout

    def test_hits_timeline_json(self, tmp_path, monkeypatch):
        """hits timeline --json outputs valid JSON."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-15", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-15",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "timeline", "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]["period"] == "2024-01"
        assert data[0]["new_hits"] == 1

    def test_hits_timeline_json_no_gap_fill(self, tmp_path, monkeypatch):
        """hits timeline --json does not fill gaps (sparse data)."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        # Hit in January and March (skip February)
        hit1, _ = register_hit("ACGTACGTACGT", "2024-01-15", tmp_path)
        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_a",
            "2024-01-15",
            state_dir=tmp_path,
        )
        hit2, _ = register_hit("GGGGCCCCAAAA", "2024-03-15", tmp_path)
        record_hit_observation(
            hit2.hit_key,
            "set_001",
            "sample_b",
            "2024-03-15",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "timeline", "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        # JSON output should be sparse (no February)
        assert len(data) == 2
        periods = [d["period"] for d in data]
        assert "2024-01" in periods
        assert "2024-02" not in periods
        assert "2024-03" in periods

    def test_hits_recur_help(self):
        """hits recur --help shows options."""
        result = runner.invoke(app, ["hits", "recur", "--help"])
        assert result.exit_code == 0
        assert "--min-samples" in result.stdout
        assert "--min-runs" in result.stdout
        assert "--limit" in result.stdout
        assert "--json" in result.stdout

    def test_hits_recur_empty(self, tmp_path, monkeypatch):
        """hits recur handles empty database."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        result = runner.invoke(app, ["hits", "recur"])
        assert result.exit_code == 0
        assert "No recurring hits found" in result.stdout

    def test_hits_recur_no_recurring(self, tmp_path, monkeypatch):
        """hits recur shows message when no hits recur."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "recur"])
        assert result.exit_code == 0
        assert "No recurring hits found" in result.stdout

    def test_hits_recur_finds_recurring(self, tmp_path, monkeypatch):
        """hits recur finds hits in multiple samples."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "recur"])
        assert result.exit_code == 0
        assert "1 recurring hit" in result.stdout

    def test_hits_recur_json(self, tmp_path, monkeypatch):
        """hits recur --json outputs valid JSON."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01",
            state_dir=tmp_path,
        )

        result = runner.invoke(app, ["hits", "recur", "--json"])
        assert result.exit_code == 0
        import json

        data = json.loads(result.stdout)
        assert len(data) == 1
        assert data[0]["hit_key"] == hit.hit_key
        assert data[0]["sample_count"] == 2

    def test_hits_recur_min_samples_filter(self, tmp_path, monkeypatch):
        """hits recur --min-samples filters correctly."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        runner.invoke(app, ["state", "init"])

        from py_nvd.hits import record_hit_observation, register_hit

        hit, _ = register_hit("ACGTACGTACGT", "2024-01-01", tmp_path)
        # Only 2 samples
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01",
            state_dir=tmp_path,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01",
            state_dir=tmp_path,
        )

        # Should find with min-samples=2
        result = runner.invoke(app, ["hits", "recur", "--min-samples", "2"])
        assert result.exit_code == 0
        assert "1 recurring hit" in result.stdout

        # Should not find with min-samples=3
        result = runner.invoke(app, ["hits", "recur", "--min-samples", "3"])
        assert result.exit_code == 0
        assert "No recurring hits found" in result.stdout


class TestResumeFile:
    """Verify .nfresume file handling."""

    def test_resume_file_constant(self):
        """RESUME_FILE is correctly defined."""
        from py_nvd.cli.utils import RESUME_FILE

        assert RESUME_FILE.name == ".nfresume"

    def test_get_editor_fallback(self, monkeypatch):
        """get_editor falls back to vi when env vars not set."""
        from py_nvd.cli.utils import get_editor

        monkeypatch.delenv("VISUAL", raising=False)
        monkeypatch.delenv("EDITOR", raising=False)
        assert get_editor() == "vi"

    def test_get_editor_visual(self, monkeypatch):
        """get_editor prefers $VISUAL."""
        from py_nvd.cli.utils import get_editor

        monkeypatch.setenv("VISUAL", "code")
        monkeypatch.setenv("EDITOR", "nano")
        assert get_editor() == "code"

    def test_get_editor_editor(self, monkeypatch):
        """get_editor uses $EDITOR when $VISUAL not set."""
        from py_nvd.cli.utils import get_editor

        monkeypatch.delenv("VISUAL", raising=False)
        monkeypatch.setenv("EDITOR", "nano")
        assert get_editor() == "nano"

    def test_editor_min_duration_constant(self):
        """EDITOR_MIN_DURATION_SECONDS is defined for fast-exit detection."""
        from py_nvd.cli.commands.resume import EDITOR_MIN_DURATION_SECONDS

        assert EDITOR_MIN_DURATION_SECONDS > 0
        assert EDITOR_MIN_DURATION_SECONDS <= 2  # Reasonable threshold
