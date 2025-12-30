"""Smoke tests for CLI commands - verify they don't crash."""

import pytest
from typer.testing import CliRunner

from py_nvd.cli_legacy import app

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
