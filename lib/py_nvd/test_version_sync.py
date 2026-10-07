from __future__ import annotations

import json
import re
import tomllib
from pathlib import Path

from py_nvd import __version__
from py_nvd.models import NvdParams
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


def test_read_entropy_defaults_match_runtime_and_schema() -> None:
    """Native Nextflow and the Python wrapper should apply the same default."""
    nextflow_config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    config_match = re.search(
        r"^\s*min_read_entropy\s*=\s*([0-9.]+)",
        nextflow_config,
        re.MULTILINE,
    )
    latest_schema = json.loads(
        (ROOT / "schemas" / "nvd-params.latest.schema.json").read_text(
            encoding="utf-8",
        ),
    )

    assert config_match is not None, "nextflow.config min_read_entropy is missing"
    assert float(config_match.group(1)) == 0.5
    assert NvdParams().min_read_entropy == 0.5
    assert latest_schema["properties"]["min_read_entropy"]["default"] == 0.5


def test_merge_pairs_default_matches_runtime_and_schema() -> None:
    """Pair merging is on by default in Nextflow, the model, and the schema."""
    nextflow_config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    config_match = re.search(
        r"^\s*merge_pairs\s*=\s*(\w+)",
        nextflow_config,
        re.MULTILINE,
    )
    latest_schema = json.loads(
        (ROOT / "schemas" / "nvd-params.latest.schema.json").read_text(
            encoding="utf-8",
        ),
    )

    assert config_match is not None, "nextflow.config merge_pairs is missing"
    assert config_match.group(1) == "true"
    assert NvdParams().merge_pairs is True
    assert latest_schema["properties"]["merge_pairs"]["default"] is True


def test_skip_unassembled_read_queries_is_declared_everywhere() -> None:
    """The read-query skip exists in Nextflow, the model, and the schema."""
    nextflow_config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    config_match = re.search(
        r"^\s*skip_unassembled_read_queries\s*=\s*(\S+)",
        nextflow_config,
        re.MULTILINE,
    )
    latest_schema = json.loads(
        (ROOT / "schemas" / "nvd-params.latest.schema.json").read_text(
            encoding="utf-8",
        ),
    )

    assert config_match is not None, (
        "nextflow.config skip_unassembled_read_queries is missing"
    )
    # null keeps bare `--skip_unassembled_read_queries` usable as a Nextflow flag,
    # matching the other skip_* params.
    assert config_match.group(1) == "null"
    assert NvdParams().skip_unassembled_read_queries is False
    assert (
        latest_schema["properties"]["skip_unassembled_read_queries"]["default"] is False
    )


def test_preprocess_param_is_gone() -> None:
    """preprocess was never read by the pipeline and is removed in v3.4.0."""
    nextflow_config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    latest_schema = json.loads(
        (ROOT / "schemas" / "nvd-params.latest.schema.json").read_text(
            encoding="utf-8",
        ),
    )

    assert not re.search(
        r"^\s*preprocess\s*=",
        nextflow_config,
        re.MULTILINE,
    ), "nextflow.config still declares the removed preprocess param"
    assert "preprocess" not in latest_schema["properties"]
    assert "preprocess" not in NvdParams.model_fields


def test_latest_params_schema_points_to_v3_6() -> None:
    """The rolling schema link should expose the v3.6 parameter contract."""
    latest_schema = ROOT / "schemas" / "nvd-params.latest.schema.json"

    assert latest_schema.is_symlink()
    assert latest_schema.readlink() == Path("nvd-params.v3.6.0.schema.json")
    assert SCHEMA_URL.endswith("/nvd-params.v3.6.0.schema.json")


def test_v3_6_schema_starts_as_the_v3_5_contract() -> None:
    """Until the release bump, v3.6.0 only adds properties on top of v3.5.0."""
    previous = json.loads(
        (ROOT / "schemas" / "nvd-params.v3.5.0.schema.json").read_text(
            encoding="utf-8",
        ),
    )
    current = json.loads(
        (ROOT / "schemas" / "nvd-params.v3.6.0.schema.json").read_text(
            encoding="utf-8",
        ),
    )

    assert current["$id"].endswith("/nvd-params.v3.6.0.schema.json")
    for name, definition in previous["properties"].items():
        assert current["properties"][name] == definition, name
    current_without_id = {key: value for key, value in current.items() if key != "$id"}
    previous_without_id = {
        key: value for key, value in previous.items() if key != "$id"
    }
    current_without_id["properties"] = previous["properties"]
    assert current_without_id == previous_without_id


def test_v3_3_2_schema_corrects_only_the_read_entropy_default() -> None:
    """The patch schema preserves the published v3.3.0 default."""
    original_path = ROOT / "schemas" / "nvd-params.v3.3.0.schema.json"
    corrected_path = ROOT / "schemas" / "nvd-params.v3.3.2.schema.json"
    original = json.loads(original_path.read_text(encoding="utf-8"))
    corrected = json.loads(corrected_path.read_text(encoding="utf-8"))

    assert original["properties"]["min_read_entropy"]["default"] == 0.9
    assert corrected["properties"]["min_read_entropy"]["default"] == 0.5

    corrected["$id"] = original["$id"]
    corrected["properties"]["min_read_entropy"]["default"] = 0.9
    assert corrected == original


def test_v3_2_schema_accepts_disabled_optional_read_limits() -> None:
    """The public schema should accept Nextflow's null filter defaults."""
    schema_path = ROOT / "schemas" / "nvd-params.v3.2.0.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))

    assert schema["properties"]["filter_reads"]["type"] == ["boolean", "null"]
    assert schema["properties"]["max_read_length"]["type"] == ["integer", "null"]


def test_skip_big_tables_is_declared_everywhere() -> None:
    """The big-table skip exists in Nextflow, the model, and the schema."""
    nextflow_config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    config_match = re.search(
        r"^\s*skip_big_tables\s*=\s*(\S+)",
        nextflow_config,
        re.MULTILINE,
    )
    latest_schema = json.loads(
        (ROOT / "schemas" / "nvd-params.latest.schema.json").read_text(
            encoding="utf-8",
        ),
    )

    assert config_match is not None, "nextflow.config skip_big_tables is missing"
    # null keeps bare `--skip_big_tables` usable as a Nextflow flag, matching
    # the other skip_* params.
    assert config_match.group(1) == "null"
    assert NvdParams().skip_big_tables is False
    assert latest_schema["properties"]["skip_big_tables"]["default"] is False


def test_background_depletion_params_are_declared_everywhere() -> None:
    """background_index and its thresholds exist in Nextflow, the model, and the schema."""
    nextflow_config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    latest_schema = json.loads(
        (ROOT / "schemas" / "nvd-params.latest.schema.json").read_text(
            encoding="utf-8",
        ),
    )
    expected = {
        "background_index": "null",
        "background_abs_threshold": "1",
        "background_rel_threshold": "0.0",
    }
    for name, value in expected.items():
        match = re.search(rf"^\s*{name}\s*=\s*(\S+)", nextflow_config, re.MULTILINE)
        assert match is not None, f"nextflow.config {name} is missing"
        assert match.group(1) == value, name

    defaults = NvdParams()
    assert defaults.background_index is None
    assert defaults.background_abs_threshold == 1
    assert defaults.background_rel_threshold == 0.0
    assert latest_schema["properties"]["background_index"]["default"] is None
    assert latest_schema["properties"]["background_abs_threshold"]["default"] == 1
    assert latest_schema["properties"]["background_rel_threshold"]["default"] == 0.0


def test_skip_contig_filter_is_declared_everywhere() -> None:
    """The contig-screen skip exists in Nextflow, the model, and the schema."""
    nextflow_config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    config_match = re.search(
        r"^\s*skip_contig_filter\s*=\s*(\S+)",
        nextflow_config,
        re.MULTILINE,
    )
    latest_schema = json.loads(
        (ROOT / "schemas" / "nvd-params.latest.schema.json").read_text(
            encoding="utf-8",
        ),
    )

    assert config_match is not None, "nextflow.config skip_contig_filter is missing"
    assert config_match.group(1) == "null"
    assert NvdParams().skip_contig_filter is False
    assert latest_schema["properties"]["skip_contig_filter"]["default"] is False
