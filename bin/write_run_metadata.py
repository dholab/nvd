#!/usr/bin/env python3
"""Write top-level NVD invocation metadata without shell reinterpretation."""

from __future__ import annotations

import argparse
import base64
from pathlib import Path


def write_run_metadata(
    *,
    command_base64: str,
    version: str,
    revision: str,
    dirty: str,
    output_dir: Path,
) -> None:
    """Decode and write the run command and source identity artifacts."""
    command = base64.b64decode(command_base64, validate=True).decode("utf-8")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "run_command.sh").write_text(f"{command}\n", encoding="utf-8")
    (output_dir / "nvd_version.txt").write_text(
        f"version={version}\nrevision={revision}\ndirty={dirty}\n",
        encoding="utf-8",
    )


def find_installed_py_nvd_version() -> str | None:
    """Return the py_nvd version importable here, or None if it is missing.

    Nextflow stages ``bin/`` from the pipeline checkout but does NOT stage
    ``lib/``, so py_nvd comes from whatever is installed in the runtime
    environment. The two are versioned independently, which is why they can
    drift apart without any repository test noticing.
    """
    # Imported here rather than at module scope on purpose. A missing py_nvd is
    # one of the cases this function reports, and a top-level import would make
    # it a traceback before argparse ever runs, losing the actionable message.
    try:
        from py_nvd import __version__  # noqa: PLC0415
    except ImportError:
        return None
    return str(__version__)


def describe_version_skew(
    pipeline_version: str,
    installed_version: str | None,
) -> str | None:
    """Return an actionable message if the runtime does not match the pipeline.

    Returns None when they agree, which is the normal case.
    """
    if installed_version == pipeline_version:
        return None

    found = (
        "not importable at all"
        if installed_version is None
        else f"version {installed_version}"
    )
    return (
        f"NVD runtime mismatch: the pipeline is version {pipeline_version}, but the "
        f"py_nvd package in this environment is {found}.\n"
        "\n"
        "Nextflow stages bin/ from the pipeline checkout but not lib/, so the helper "
        "scripts come from the checkout while py_nvd comes from the environment. When "
        "they disagree, scripts call py_nvd with arguments it does not have and fail "
        "mid-run, which has previously cost the MultiQC report silently because that "
        "step ignores its own errors.\n"
        "\n"
        "Fix it by installing py_nvd from this checkout into the environment the run "
        "uses, for example `pixi install -e dev`, or by using a container image built "
        "from the same revision. Failing the run here is deliberate: it costs seconds "
        "now instead of a whole run with missing outputs."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--command-base64", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--revision", default="unknown")
    parser.add_argument(
        "--dirty",
        choices=("true", "false", "unknown"),
        default="unknown",
    )
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    parser.add_argument(
        "--skip-runtime-version-check",
        action="store_true",
        help=(
            "Write the metadata even if py_nvd does not match the pipeline version. "
            "For debugging a mismatch, not for routine runs."
        ),
    )
    args = parser.parse_args()

    if not args.skip_runtime_version_check:
        skew = describe_version_skew(args.version, find_installed_py_nvd_version())
        if skew is not None:
            raise SystemExit(skew)

    write_run_metadata(
        command_base64=args.command_base64,
        version=args.version,
        revision=args.revision,
        dirty=args.dirty,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
