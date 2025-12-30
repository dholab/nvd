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
