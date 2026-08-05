"""Tests for setup helpers after the v3 config-dir refactor."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from py_nvd.cli.commands.setup import (
    CHTC_DEFAULT_CONFIG_DIR,
    CHTC_DEFAULT_PRESET_STORE,
    CHTC_DEFAULT_TAXONOMY_DIR,
    _detect_shell,
    _generate_setup_conf,
    _generate_shell_hook_block,
    _generate_wrapper_script,
    _guard_bash_completion_script,
    _install_shell_hook,
    _is_oconnor_chtc,
    _migrate_legacy_wrapper,
    _path_contains_directory,
    _validate_config_dir,
    _write_setup_conf,
    setup_app,
)
from py_nvd.cli.utils import get_pipeline_root

runner = CliRunner()


def test_detect_shell_bash(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SHELL", "/bin/bash")
    assert _detect_shell() == "bash"


def test_shell_hook_block_bootstraps_wrapper_without_path() -> None:
    block = _generate_shell_hook_block("zsh")

    assert '$HOME/.local/bin/nvd" setup shell-hook --shell zsh' in block
    assert "$(nvd setup shell-hook)" not in block


def test_install_shell_hook_replaces_legacy_path_dependent_hook(
    tmp_path: Path,
) -> None:
    rc_file = tmp_path / ".zshrc"
    rc_file.write_text(
        'export EDITOR=vim\n\n# NVD shell integration\neval "$(nvd setup shell-hook)"\n',
    )

    assert _install_shell_hook(rc_file, "zsh")

    updated = rc_file.read_text()
    assert "export EDITOR=vim" in updated
    assert updated.count("# NVD shell integration") == 1
    assert "$(nvd setup shell-hook)" not in updated
    assert '$HOME/.local/bin/nvd" setup shell-hook --shell zsh' in updated


def test_path_directory_check_uses_complete_path_entries(tmp_path: Path) -> None:
    wrapper_dir = tmp_path / ".local" / "bin"

    assert _path_contains_directory(
        wrapper_dir,
        f"/usr/bin:{wrapper_dir}:/bin",
    )
    assert not _path_contains_directory(
        wrapper_dir,
        f"/usr/bin:{wrapper_dir}-old:/bin",
    )


def test_is_oconnor_chtc_true_on_matching_host() -> None:
    with (
        patch("socket.gethostname", return_value="oconnor-ap2001"),
        patch("socket.getfqdn", return_value="oconnor-ap2001.chtc.wisc.edu"),
    ):
        assert _is_oconnor_chtc()


def test_generate_setup_conf_writes_config_dir() -> None:
    conf = _generate_setup_conf(
        Path("/example/nvd-config"),
        default_profile="chtc_htc",
    )

    assert "NVD_REPO=" not in conf
    assert "NVD_CONFIG_DIR=/example/nvd-config" in conf
    assert "NVD_DEFAULT_PROFILE=chtc_htc" in conf
    assert "NVD_STATE_DIR" not in conf


def test_generate_setup_conf_writes_chtc_preset_store() -> None:
    conf = _generate_setup_conf(
        CHTC_DEFAULT_CONFIG_DIR,
        default_profile="chtc_htc",
        preset_store=CHTC_DEFAULT_PRESET_STORE,
    )

    assert f"NVD_CONFIG_DIR={Path.home() / '.nvd'}" in conf
    assert "NVD_PRESET_STORE=/staging/groups/oconnor_group/nvd/presets.sqlite" in conf
    assert "NVD_DEFAULT_PROFILE=chtc_htc" in conf
    assert "NVD_STATE_DIR" not in conf


def test_generate_setup_conf_writes_chtc_taxonomy_dir() -> None:
    conf = _generate_setup_conf(
        CHTC_DEFAULT_CONFIG_DIR,
        default_profile="chtc_htc",
        taxonomy_dir=CHTC_DEFAULT_TAXONOMY_DIR,
    )

    assert "NVD_TAXONOMY_DB=/staging/groups/oconnor_group/nvd/taxdump" in conf
    assert "NVD_DEFAULT_PROFILE=chtc_htc" in conf
    assert "NVD_STATE_DIR" not in conf


def test_shell_hook_exports_preset_store_from_setup_conf(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_dir = tmp_path / "config"
    preset_store = tmp_path / "shared" / "presets.sqlite"
    config_dir.mkdir()
    (config_dir / "setup.conf").write_text(
        _generate_setup_conf(
            config_dir,
            default_profile="chtc_htc",
            preset_store=preset_store,
        ),
    )
    monkeypatch.setenv("NVD_CONFIG_DIR", str(config_dir))

    result = runner.invoke(setup_app, ["shell-hook", "--shell", "bash"])

    assert result.exit_code == 0, result.output
    assert f'export NVD_CONFIG_DIR="{config_dir}"' in result.output
    assert f'export NVD_PRESET_STORE="{preset_store}"' in result.output
    assert f"{config_dir / 'completions.bash'}" in result.output


def test_shell_hook_exports_taxonomy_dir_from_setup_conf(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_dir = tmp_path / "config"
    taxonomy_dir = tmp_path / "shared" / "taxdump"
    config_dir.mkdir()
    (config_dir / "setup.conf").write_text(
        _generate_setup_conf(
            config_dir,
            default_profile="chtc_htc",
            taxonomy_dir=taxonomy_dir,
        ),
    )
    monkeypatch.setenv("NVD_CONFIG_DIR", str(config_dir))

    result = runner.invoke(setup_app, ["shell-hook", "--shell", "bash"])

    assert result.exit_code == 0, result.output
    assert f'export NVD_CONFIG_DIR="{config_dir}"' in result.output
    assert f'export NVD_TAXONOMY_DB="{taxonomy_dir}"' in result.output
    assert f"{config_dir / 'completions.bash'}" in result.output


def test_generated_wrapper_selects_one_repository(
    tmp_path: Path,
) -> None:
    repos = {
        "configured": tmp_path / "configured $repo",
        "canonical": tmp_path / "canonical",
        "legacy": tmp_path / "legacy",
    }
    for label, repo in repos.items():
        nvd_bin = repo / ".pixi" / "envs" / "default" / "bin" / "nvd"
        nvd_bin.parent.mkdir(parents=True)
        nvd_bin.write_text(f"#!/bin/bash\nprintf '%s\\n' '{label}'\n")
        nvd_bin.chmod(0o755)

    wrapper = tmp_path / "nvd"
    wrapper.write_text(_generate_wrapper_script(repos["configured"]))
    wrapper.chmod(0o755)

    cases = (
        ({}, "configured"),
        ({"NVD_REPO": str(repos["canonical"])}, "canonical"),
        ({"NVD_PIPELINE_ROOT": str(repos["legacy"])}, "legacy"),
    )
    for overrides, expected in cases:
        env = os.environ.copy()
        env.pop("NVD_REPO", None)
        env.pop("NVD_PIPELINE_ROOT", None)
        env.update(overrides)

        result = subprocess.run(  # noqa: S603
            [wrapper],
            check=True,
            capture_output=True,
            env=env,
            text=True,
        )

        assert result.stdout.strip() == expected


def test_pipeline_root_rejects_conflicting_repository_overrides(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("NVD_REPO", str(tmp_path / "canonical"))
    monkeypatch.setenv("NVD_PIPELINE_ROOT", str(tmp_path / "legacy"))

    with pytest.raises(
        RuntimeError,
        match="NVD_REPO and NVD_PIPELINE_ROOT select different repositories",
    ):
        get_pipeline_root()


def test_setup_migrates_only_managed_wrapper_and_preserves_setup_conf(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    managed_wrapper = tmp_path / "managed" / "nvd"
    managed_wrapper.parent.mkdir()
    configured_repo = Path("/home/user/.nvd/latest")
    managed_wrapper.write_text(
        """#!/bin/bash
# NVD wrapper - runs nvd directly from pixi environment
# Generated by: nvd setup v3.4.0
# Date: 2026-07-31
#
# This wrapper sets PATH to include pixi environment binaries (for nextflow, etc.)
# but only for the nvd process - it doesn't pollute the user's shell PATH.

NVD_REPO="/home/user/.nvd/latest"
NVD_BIN="$NVD_REPO/.pixi/envs/default/bin"

# Verify the environment exists
if [[ ! -x "$NVD_BIN/nvd" ]]; then
    echo "Error: nvd not found at $NVD_BIN/nvd" >&2
    echo "Run 'pixi install' in $NVD_REPO to set up the environment." >&2
    exit 1
fi

# Run nvd with pixi bin on PATH (for subprocess access to nextflow, etc.)
PATH="$NVD_BIN:$PATH" exec "$NVD_BIN/nvd" "$@"
""",
    )
    unknown_wrapper = tmp_path / "unknown" / "nvd"
    unknown_wrapper.parent.mkdir()
    unknown_wrapper.write_text("#!/bin/bash\nprintf 'custom wrapper\\n'\n")

    config_dir = tmp_path / "config"
    config_dir.mkdir()
    setup_conf = config_dir / "setup.conf"
    original_setup_conf = "NVD_REPO=/home/user/.nvd/latest\nCUSTOM=value\n"
    setup_conf.write_text(original_setup_conf)
    monkeypatch.setenv("NVD_CONFIG_DIR", str(config_dir))

    assert _migrate_legacy_wrapper(configured_repo, managed_wrapper)
    assert "export NVD_REPO=" in managed_wrapper.read_text()
    assert not _migrate_legacy_wrapper(
        configured_repo,
        unknown_wrapper,
    )
    assert unknown_wrapper.read_text() == "#!/bin/bash\nprintf 'custom wrapper\\n'\n"

    _write_setup_conf(config_dir)
    assert setup_conf.read_text() == original_setup_conf


def test_validate_config_dir_creates_missing_directory(tmp_path: Path) -> None:
    target = tmp_path / "nvd-config"
    valid, message = _validate_config_dir(target)
    assert valid, message
    assert target.is_dir()


def test_bash_completion_registration_is_guarded() -> None:
    script = """_nvd_completion() {
    return 0
}

complete -o default -F _nvd_completion nvd
"""

    guarded = _guard_bash_completion_script(script)

    assert "if type complete >/dev/null 2>&1; then" in guarded
    assert "    complete -o default -F _nvd_completion nvd" in guarded
    assert guarded.rstrip().endswith("fi")
