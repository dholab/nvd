from __future__ import annotations

import base64
from pathlib import Path

from write_run_metadata import (
    describe_version_skew,
    find_installed_py_nvd_version,
    write_run_metadata,
)


def test_write_run_metadata_preserves_command_and_source_identity(
    tmp_path: Path,
) -> None:
    output_dir = Path(tmp_path)
    command = "nextflow run './pipeline dir' --label '$HOME `literal`'"
    encoded = base64.b64encode(command.encode("utf-8")).decode("ascii")

    write_run_metadata(
        command_base64=encoded,
        version="3.0.0",
        revision="abc123",
        dirty="false",
        output_dir=output_dir,
    )

    assert (output_dir / "run_command.sh").read_text(encoding="utf-8") == f"{command}\n"
    assert (output_dir / "nvd_version.txt").read_text(encoding="utf-8") == (
        "version=3.0.0\nrevision=abc123\ndirty=false\n"
    )


def test_version_skew_is_silent_when_the_runtime_matches() -> None:
    """The normal case must not produce a message."""
    assert describe_version_skew("3.6.0", "3.6.0") is None


def test_version_skew_names_both_versions_and_how_to_fix_it() -> None:
    """The message has to be actionable from the Nextflow log alone.

    This is the failure that silently cost the MultiQC report on the scheduled
    runs: bin/ was newer than the installed py_nvd, so a helper script called
    ReportConfiguration with a keyword it did not have.
    """
    message = describe_version_skew("3.6.0", "3.3.2")

    assert message is not None
    assert "3.6.0" in message
    assert "3.3.2" in message
    # Names the mechanism, so the reader knows why two versions exist at all.
    assert "stages bin/" in message
    # And says what to do about it.
    assert "pixi install" in message


def test_version_skew_reports_a_missing_py_nvd_differently() -> None:
    """An absent package is a different problem from a stale one."""
    message = describe_version_skew("3.6.0", None)

    assert message is not None
    assert "not importable at all" in message


def test_the_installed_version_is_discoverable_in_this_environment() -> None:
    """In the repo's own environment the two must already agree.

    If this fails, the dev environment itself is skewed and every helper script
    is a coin flip.
    """
    from py_nvd import __version__

    assert find_installed_py_nvd_version() == __version__
