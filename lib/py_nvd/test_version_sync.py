from __future__ import annotations

import json
import re
import tomllib
from pathlib import Path

from py_nvd import __version__
from py_nvd.params import SCHEMA_URL

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


def test_latest_params_schema_points_to_v3_2() -> None:
    """The rolling schema link should expose the v3.2 parameter contract."""
    latest_schema = ROOT / "schemas" / "nvd-params.latest.schema.json"

    assert latest_schema.is_symlink()
    assert latest_schema.readlink() == Path("nvd-params.v3.2.0.schema.json")
    assert SCHEMA_URL.endswith("/nvd-params.v3.2.0.schema.json")


def test_v3_2_schema_accepts_disabled_optional_read_limits() -> None:
    """The public schema should accept Nextflow's null filter defaults."""
    schema_path = ROOT / "schemas" / "nvd-params.v3.2.0.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))

    assert schema["properties"]["filter_reads"]["type"] == ["boolean", "null"]
    assert schema["properties"]["max_read_length"]["type"] == ["integer", "null"]
