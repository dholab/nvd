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
    args = parser.parse_args()
    write_run_metadata(
        command_base64=args.command_base64,
        version=args.version,
        revision=args.revision,
        dirty=args.dirty,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
