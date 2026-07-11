from __future__ import annotations

import base64
from pathlib import Path

from write_run_metadata import write_run_metadata


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
