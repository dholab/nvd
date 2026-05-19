"""Smoke tests for the retained v3 CLI surface."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

from typer.testing import CliRunner

from py_nvd.cli.app import app

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


runner = CliRunner()


def test_retained_top_level_commands_show_help() -> None:
    for command in ["run", "params", "preset", "config", "taxonomy", "version", "setup"]:
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

    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps({"samplesheet": "samples.csv", "dedup": True}))

    imported = runner.invoke(
        app,
        ["preset", "import", str(params_path), "--name", "smoke", "--description", "test"],
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


def test_taxonomy_status_uses_explicit_taxonomy_dir(tmp_path: Path) -> None:
    taxonomy_dir = tmp_path / "taxdump"
    result = runner.invoke(
        app,
        ["taxonomy", "status", "--taxonomy-dir", str(taxonomy_dir), "--json"],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["taxdump_dir"] == str(taxonomy_dir)
