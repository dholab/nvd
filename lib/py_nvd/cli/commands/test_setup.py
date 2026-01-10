"""Tests for the setup command."""

import os
from pathlib import Path
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from py_nvd import __version__
from py_nvd.cli.app import app
from py_nvd.cli.commands.setup import (
    CONTAINER_IMAGE_URI,
    SHELL_HOOK_MARKER,
    _check_pixi_installed,
    _detect_shell,
    _generate_setup_conf,
    _generate_wrapper_script,
    _get_rc_file,
    _is_hook_installed,
    _is_oconnor_chtc,
    _load_and_substitute_template,
    _validate_state_dir,
)

runner = CliRunner()


class TestEnvironmentDetection:
    """Tests for environment detection functions."""

    def test_is_oconnor_chtc_false_on_generic_system(self):
        """_is_oconnor_chtc returns False on non-CHTC systems."""
        with patch("socket.gethostname", return_value="generic-workstation"):
            with patch("socket.getfqdn", return_value="generic-workstation.local"):
                assert _is_oconnor_chtc() is False

    def test_is_oconnor_chtc_true_on_chtc(self):
        """_is_oconnor_chtc returns True when hostname matches pattern."""
        with patch("socket.gethostname", return_value="oconnor-ap2001"):
            with patch("socket.getfqdn", return_value="oconnor-ap2001.chtc.wisc.edu"):
                assert _is_oconnor_chtc() is True

    def test_is_oconnor_chtc_false_wrong_hostname(self):
        """_is_oconnor_chtc returns False when hostname doesn't match."""
        with patch("socket.gethostname", return_value="some-other-host"):
            with patch("socket.getfqdn", return_value="some-other-host.chtc.wisc.edu"):
                assert _is_oconnor_chtc() is False

    def test_is_oconnor_chtc_false_wrong_domain(self):
        """_is_oconnor_chtc returns False when domain doesn't match."""
        with patch("socket.gethostname", return_value="oconnor-ap2001"):
            with patch("socket.getfqdn", return_value="oconnor-ap2001.example.com"):
                assert _is_oconnor_chtc() is False

    def test_detect_shell_bash(self, monkeypatch):
        """_detect_shell returns 'bash' for bash shell."""
        monkeypatch.setenv("SHELL", "/bin/bash")
        assert _detect_shell() == "bash"

    def test_detect_shell_zsh(self, monkeypatch):
        """_detect_shell returns 'zsh' for zsh shell."""
        monkeypatch.setenv("SHELL", "/bin/zsh")
        assert _detect_shell() == "zsh"

    def test_detect_shell_zsh_usr_local(self, monkeypatch):
        """_detect_shell handles non-standard zsh paths."""
        monkeypatch.setenv("SHELL", "/usr/local/bin/zsh")
        assert _detect_shell() == "zsh"

    def test_detect_shell_unknown(self, monkeypatch):
        """_detect_shell returns None for unknown shells."""
        monkeypatch.setenv("SHELL", "/bin/fish")
        assert _detect_shell() is None

    def test_detect_shell_empty(self, monkeypatch):
        """_detect_shell returns None when SHELL is empty."""
        monkeypatch.delenv("SHELL", raising=False)
        assert _detect_shell() is None

    def test_check_pixi_installed_true(self):
        """_check_pixi_installed returns True when pixi is available."""
        with patch("shutil.which", return_value="/usr/local/bin/pixi"):
            assert _check_pixi_installed() is True

    def test_check_pixi_installed_false(self):
        """_check_pixi_installed returns False when pixi is not found."""
        with patch("shutil.which", return_value=None):
            assert _check_pixi_installed() is False


class TestTemplateSubstitution:
    """Tests for template loading and variable substitution."""

    def test_load_and_substitute_template(self, tmp_path):
        """_load_and_substitute_template substitutes variables correctly."""
        template = tmp_path / "test.config"
        template.write_text(
            "// User: ${USER}\n"
            "// Hostname: ${HOSTNAME}\n"
            "// Date: ${DATE}\n"
            "// Version: ${VERSION}\n"
            "// Home: ${HOME}\n"  # Should NOT be substituted
        )

        result = _load_and_substitute_template(template)

        assert "${USER}" not in result
        assert "${HOSTNAME}" not in result
        assert "${DATE}" not in result
        assert "${VERSION}" not in result
        assert "${HOME}" in result  # Should remain for Nextflow
        assert __version__ in result

    def test_load_and_substitute_template_missing_file(self, tmp_path):
        """_load_and_substitute_template raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            _load_and_substitute_template(tmp_path / "nonexistent.config")


class TestWrapperScriptGeneration:
    """Tests for wrapper script generation."""

    def test_generate_wrapper_script(self, tmp_path):
        """_generate_wrapper_script generates valid bash script."""
        nvd_repo = tmp_path / "nvd"
        script = _generate_wrapper_script(nvd_repo)

        assert script.startswith("#!/bin/bash")
        assert str(nvd_repo) in script
        assert "pixi run" in script
        assert "nvd" in script
        assert __version__ in script

    def test_generate_setup_conf(self, tmp_path):
        """_generate_setup_conf generates valid config."""
        nvd_repo = tmp_path / "nvd"
        state_dir = tmp_path / "state"

        conf = _generate_setup_conf(nvd_repo, state_dir)

        assert f"NVD_REPO={nvd_repo}" in conf
        assert f"NVD_STATE_DIR={state_dir}" in conf
        assert __version__ in conf


class TestRCFileIntegration:
    """Tests for RC file detection and hook installation."""

    def test_get_rc_file_bash(self):
        """_get_rc_file returns .bashrc for bash."""
        rc = _get_rc_file("bash")
        assert rc.name == ".bashrc"
        assert rc.parent == Path.home()

    def test_get_rc_file_zsh(self):
        """_get_rc_file returns .zshrc for zsh."""
        rc = _get_rc_file("zsh")
        assert rc.name == ".zshrc"
        assert rc.parent == Path.home()

    def test_is_hook_installed_true(self, tmp_path):
        """_is_hook_installed returns True when hook is present."""
        rc_file = tmp_path / ".bashrc"
        rc_file.write_text(f'eval "$({SHELL_HOOK_MARKER})"\n')

        assert _is_hook_installed(rc_file) is True

    def test_is_hook_installed_false(self, tmp_path):
        """_is_hook_installed returns False when hook is absent."""
        rc_file = tmp_path / ".bashrc"
        rc_file.write_text("# some other content\n")

        assert _is_hook_installed(rc_file) is False

    def test_is_hook_installed_missing_file(self, tmp_path):
        """_is_hook_installed returns False for missing file."""
        rc_file = tmp_path / ".bashrc"
        assert _is_hook_installed(rc_file) is False


class TestStateDirValidation:
    """Tests for state directory validation."""

    def test_validate_state_dir_exists_writable(self, tmp_path):
        """_validate_state_dir succeeds for existing writable directory."""
        is_valid, message = _validate_state_dir(tmp_path)
        assert is_valid is True
        assert "writable" in message.lower()

    def test_validate_state_dir_creates_new(self, tmp_path):
        """_validate_state_dir creates new directory if needed."""
        new_dir = tmp_path / "new_state"
        is_valid, message = _validate_state_dir(new_dir)

        assert is_valid is True
        assert new_dir.exists()
        assert "Created" in message

    def test_validate_state_dir_parent_missing(self, tmp_path):
        """_validate_state_dir fails when parent doesn't exist."""
        bad_path = tmp_path / "nonexistent" / "state"
        is_valid, message = _validate_state_dir(bad_path)

        assert is_valid is False
        assert "Parent directory does not exist" in message

    def test_validate_state_dir_not_a_directory(self, tmp_path):
        """_validate_state_dir fails when path is a file."""
        file_path = tmp_path / "not_a_dir"
        file_path.touch()

        is_valid, message = _validate_state_dir(file_path)

        assert is_valid is False
        assert "not a directory" in message


class TestShellHookOutput:
    """Tests for shell-hook command output."""

    def test_shell_hook_bash(self, monkeypatch, tmp_path):
        """shell-hook outputs bash-compatible code."""
        monkeypatch.setenv("SHELL", "/bin/bash")
        # Create setup.conf so state_dir is read
        nvd_home = tmp_path / ".nvd"
        nvd_home.mkdir()
        (nvd_home / "setup.conf").write_text(f"NVD_STATE_DIR={tmp_path}\n")

        with patch(
            "py_nvd.cli.commands.setup.NVD_HOME",
            nvd_home,
        ):
            result = runner.invoke(app, ["setup", "shell-hook", "--shell", "bash"])

        assert result.exit_code == 0
        assert "_NVD_COMPLETE=bash_source" in result.stdout
        assert "export PATH" in result.stdout
        assert "export NVD_STATE_DIR" in result.stdout

    def test_shell_hook_zsh(self, monkeypatch, tmp_path):
        """shell-hook outputs zsh-compatible code."""
        monkeypatch.setenv("SHELL", "/bin/zsh")
        nvd_home = tmp_path / ".nvd"
        nvd_home.mkdir()
        (nvd_home / "setup.conf").write_text(f"NVD_STATE_DIR={tmp_path}\n")

        with patch(
            "py_nvd.cli.commands.setup.NVD_HOME",
            nvd_home,
        ):
            result = runner.invoke(app, ["setup", "shell-hook", "--shell", "zsh"])

        assert result.exit_code == 0
        assert "_NVD_COMPLETE=zsh_source" in result.stdout

    def test_shell_hook_unsupported_shell(self):
        """shell-hook rejects unsupported shells."""
        result = runner.invoke(app, ["setup", "shell-hook", "--shell", "fish"])
        assert result.exit_code == 1
        assert "Unsupported shell" in result.stdout


class TestSetupCommandHelp:
    """Tests for setup command help."""

    def test_setup_help(self):
        """setup --help exits cleanly."""
        result = runner.invoke(app, ["setup", "--help"])
        assert result.exit_code == 0
        assert "--non-interactive" in result.stdout
        assert "--shell" in result.stdout
        assert "--state-dir" in result.stdout
        assert "--skip-shell-hook" in result.stdout
        assert "--skip-container" in result.stdout
        assert "--force" in result.stdout

    def test_shell_hook_help(self):
        """setup shell-hook --help exits cleanly."""
        result = runner.invoke(app, ["setup", "shell-hook", "--help"])
        assert result.exit_code == 0
        assert "--shell" in result.stdout


class TestSetupCommandNonInteractive:
    """Tests for non-interactive setup mode."""

    def test_setup_non_interactive_creates_wrapper(self, tmp_path, monkeypatch):
        """setup --non-interactive creates wrapper script."""
        # Set up isolated environment
        nvd_home = tmp_path / ".nvd"
        local_bin = tmp_path / ".local" / "bin"
        monkeypatch.setenv("HOME", str(tmp_path))

        with patch("py_nvd.cli.commands.setup.NVD_HOME", nvd_home):
            with patch("py_nvd.cli.commands.setup.DEFAULT_WRAPPER_DIR", local_bin):
                with patch("py_nvd.cli.commands.setup._get_rc_file") as mock_rc:
                    # Point RC file to temp location
                    mock_rc.return_value = tmp_path / ".zshrc"
                    (tmp_path / ".zshrc").touch()

                    result = runner.invoke(
                        app,
                        [
                            "setup",
                            "--non-interactive",
                            "--skip-shell-hook",
                            "--state-dir",
                            str(tmp_path / "state"),
                        ],
                    )

        assert result.exit_code == 0
        assert "Wrapper script installed" in result.stdout
        assert (local_bin / "nvd").exists()

    def test_setup_non_interactive_creates_setup_conf(self, tmp_path, monkeypatch):
        """setup --non-interactive creates setup.conf."""
        nvd_home = tmp_path / ".nvd"
        local_bin = tmp_path / ".local" / "bin"
        monkeypatch.setenv("HOME", str(tmp_path))

        with patch("py_nvd.cli.commands.setup.NVD_HOME", nvd_home):
            with patch("py_nvd.cli.commands.setup.DEFAULT_WRAPPER_DIR", local_bin):
                with patch("py_nvd.cli.commands.setup._get_rc_file") as mock_rc:
                    mock_rc.return_value = tmp_path / ".zshrc"
                    (tmp_path / ".zshrc").touch()

                    result = runner.invoke(
                        app,
                        [
                            "setup",
                            "--non-interactive",
                            "--skip-shell-hook",
                            "--state-dir",
                            str(tmp_path / "state"),
                        ],
                    )

        assert result.exit_code == 0
        assert (nvd_home / "setup.conf").exists()

        conf_content = (nvd_home / "setup.conf").read_text()
        assert "NVD_REPO=" in conf_content
        assert f"NVD_STATE_DIR={tmp_path / 'state'}" in conf_content

    def test_setup_non_interactive_invalid_state_dir(self, tmp_path, monkeypatch):
        """setup --non-interactive fails with invalid state dir."""
        nvd_home = tmp_path / ".nvd"
        monkeypatch.setenv("HOME", str(tmp_path))

        with patch("py_nvd.cli.commands.setup.NVD_HOME", nvd_home):
            result = runner.invoke(
                app,
                [
                    "setup",
                    "--non-interactive",
                    "--state-dir",
                    "/nonexistent/parent/path",
                ],
            )

        assert result.exit_code == 1
        assert "Invalid state directory" in result.stdout
