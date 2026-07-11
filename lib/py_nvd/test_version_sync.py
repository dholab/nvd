from __future__ import annotations

import re
import tomllib
from pathlib import Path

from py_nvd import __version__

ROOT = Path(__file__).resolve().parents[2]


def test_python_and_nextflow_versions_match_project_version() -> None:
    project_version = tomllib.loads(
        (ROOT / "pyproject.toml").read_text(encoding="utf-8"),
    )["project"]["version"]
    nextflow_config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    manifest_match = re.search(
        r"^\s*version\s*=\s*['\"]([^'\"]+)['\"]",
        nextflow_config,
        re.MULTILINE,
    )

    assert manifest_match is not None, "nextflow.config manifest.version is missing"
    assert __version__ == project_version
    assert manifest_match.group(1) == project_version
