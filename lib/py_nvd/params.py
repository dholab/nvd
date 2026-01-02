"""
Parameter file handling and validation.

Provides functions for loading, validating, and generating parameter files
using the NVD JSON Schema as the source of truth.

The schema is loaded from the repository's schemas/ directory. This module
handles locating the schema whether running from a development checkout or
an installed package.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema
import yaml

# Schema filename (symlink to current version)
SCHEMA_FILENAME = "nvd-params.latest.schema.json"

# GitHub raw URL for schema (fallback and for generated templates)
SCHEMA_URL = "https://raw.githubusercontent.com/dhoconno/nvd/state-management/schemas/nvd-params.v2.3.0.schema.json"


def _find_schema_path() -> Path:
    """
    Locate the schema file.

    Searches for the schema in the following order:
    1. Relative to this module (development: lib/py_nvd/ -> ../../schemas/)
    2. Relative to the package installation

    Returns:
        Path to the schema file

    Raises:
        FileNotFoundError: If schema cannot be located
    """
    # Try relative to this module (works in development)
    module_dir = Path(__file__).parent
    repo_root = module_dir.parent.parent  # lib/py_nvd/ -> lib/ -> repo root
    schema_path = repo_root / "schemas" / SCHEMA_FILENAME

    if schema_path.exists():
        return schema_path

    # If we get here, schema wasn't found
    msg = (
        f"Could not locate schema file '{SCHEMA_FILENAME}'.\n"
        f"Searched: {schema_path}\n"
        f"This may indicate an incomplete installation."
    )
    raise FileNotFoundError(msg)


def get_schema() -> dict[str, Any]:
    """
    Load the NVD params schema.

    Returns:
        The parsed JSON schema as a dictionary

    Raises:
        FileNotFoundError: If schema file cannot be located
        json.JSONDecodeError: If schema is not valid JSON
    """
    schema_path = _find_schema_path()
    with open(schema_path, encoding="utf-8") as f:
        return json.load(f)


def get_schema_url() -> str:
    """
    Get the URL for the schema (for use in generated files).

    Returns the $id from the schema if available, otherwise falls back
    to the default GitHub raw URL.
    """
    try:
        schema = get_schema()
        return schema.get("$id", SCHEMA_URL)
    except FileNotFoundError:
        return SCHEMA_URL


def validate_params(params: dict[str, Any]) -> list[str]:
    """
    Validate parameters against the NVD schema.

    Args:
        params: Dictionary of parameter names to values

    Returns:
        List of error messages (empty if valid)
    """
    schema = get_schema()
    validator = jsonschema.Draft202012Validator(schema)

    errors = []
    for error in validator.iter_errors(params):
        # Format error message nicely
        if error.path:
            path = ".".join(str(p) for p in error.path)
            errors.append(f"{path}: {error.message}")
        else:
            errors.append(error.message)

    return errors


def load_params_file(path: Path | str) -> dict[str, Any]:
    """
    Load parameters from a YAML or JSON file.

    Strips schema reference ($schema) and returns just the params.

    Args:
        path: Path to the params file (.yaml, .yml, or .json)

    Returns:
        Dictionary of parameter names to values

    Raises:
        FileNotFoundError: If file doesn't exist
        yaml.YAMLError: If YAML parsing fails
        json.JSONDecodeError: If JSON parsing fails
    """
    path = Path(path)

    with open(path, encoding="utf-8") as f:
        data = yaml.safe_load(f) if path.suffix in (".yaml", ".yml") else json.load(f)

    # Handle empty files
    if data is None:
        return {}

    # Remove schema reference if present
    if isinstance(data, dict):
        data.pop("$schema", None)

    return data


def validate_params_file(path: Path | str) -> list[str]:
    """
    Load and validate a params file against the NVD schema.

    Args:
        path: Path to the params file

    Returns:
        List of error messages (empty if valid)
    """
    params = load_params_file(path)
    return validate_params(params)


def generate_template(path: Path | str, output_format: str = "yaml") -> None:
    """
    Generate a template params file with schema reference.

    The generated file includes common parameters with defaults and
    comments explaining each option. Edit the file in your IDE to get
    autocomplete and validation.

    Args:
        path: Output file path
        output_format: Output format ("yaml" or "json")
    """
    path = Path(path)
    schema = get_schema()
    schema_url = get_schema_url()

    if output_format == "yaml":
        _generate_yaml_template(path, schema, schema_url)
    else:
        _generate_json_template(path, schema, schema_url)


def _format_yaml_value(value: Any) -> str:
    """Format a value for YAML output."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, str):
        # Quote strings that might be misinterpreted
        if value in ("true", "false", "null", "yes", "no", "on", "off"):
            return f'"{value}"'
        return value
    return str(value)


def _generate_yaml_template(path: Path, schema: dict, schema_url: str) -> None:
    """Generate YAML template with comments."""
    lines = [
        f"# yaml-language-server: $schema={schema_url}",
        "#",
        "# NVD Pipeline Parameters",
        "# Generated by: nvd params init",
        "#",
        "# Edit this file in your IDE for autocomplete and validation.",
        "#",
        "# Use with:",
        f"#   nextflow run dhoconno/nvd -params-file {path.name}",
        "#",
        "# Or register as a preset:",
        f"#   nvd preset register my-preset --from-file {path.name}",
        "",
    ]

    properties = schema.get("properties", {})

    # Required/common params first
    lines.append("# === Required Parameters ===")
    common = ["samplesheet", "experiment_id", "tools", "results"]
    for name in common:
        if name in properties:
            prop = properties[name]
            desc = prop.get("description", "")
            # These are required, so show as null (user must fill in)
            lines.append(f"{name}: null  # {desc}")
    lines.append("")

    # Analysis settings (with defaults)
    lines.append("# === Analysis Settings ===")
    analysis = ["cutoff_percent", "entropy", "tax_stringency"]
    for name in analysis:
        if name in properties:
            prop = properties[name]
            default = prop.get("default")
            desc = prop.get("description", "")
            lines.append(f"{name}: {_format_yaml_value(default)}  # {desc}")
    lines.append("")

    # Preprocessing
    lines.append("# === Preprocessing ===")
    preprocess_params = [
        "preprocess",
        "merge_pairs",
        "dedup",
        "trim_adapters",
        "scrub_host_reads",
        "filter_reads",
    ]
    for name in preprocess_params:
        if name in properties:
            prop = properties[name]
            default = prop.get("default")
            desc = prop.get("description", "")
            value = _format_yaml_value(default)
            if name == "preprocess":
                # Show preprocess uncommented since it's the main toggle
                lines.append(f"{name}: {value}  # {desc}")
            else:
                # Comment out individual preprocessing options
                lines.append(f"# {name}: {value}  # {desc}")
    lines.append("")

    # LabKey integration (commented out - optional)
    lines.append("# === LabKey Integration ===")
    lines.append("# Set labkey: true and configure these for LabKey uploads.")
    labkey_params = [
        "labkey",
        "labkey_server",
        "labkey_project_name",
        "labkey_webdav",
        "labkey_schema",
        "labkey_gottcha_fasta_list",
        "labkey_gottcha_full_list",
        "labkey_gottcha_blast_verified_full_list",
        "labkey_blast_meta_hits_list",
        "labkey_blast_fasta_list",
        "labkey_exp_id_guard_list",
    ]
    for name in labkey_params:
        if name in properties:
            prop = properties[name]
            default = prop.get("default")
            desc = prop.get("description", "")
            value = _format_yaml_value(default)
            lines.append(f"# {name}: {value}  # {desc}")
    lines.append("")

    # Database paths (commented out - environment specific)
    lines.append("# === Database Paths ===")
    lines.append(
        "# These paths are environment-specific. Set them in your config or here.",
    )
    databases = [
        "gottcha2_db",
        "blast_db",
        "blast_db_prefix",
        "stat_index",
        "stat_dbss",
        "stat_annotation",
        "human_virus_taxlist",
        "nvd_files",
    ]
    for name in databases:
        if name in properties:
            prop = properties[name]
            desc = prop.get("description", "")
            lines.append(f"# {name}: null  # {desc}")
    lines.append("")

    # Read quality/length filtering (commented out)
    lines.append("# === Read Filtering ===")
    filtering = [
        "min_read_quality_illumina",
        "min_read_quality_nanopore",
        "min_read_length",
        "max_read_length",
    ]
    for name in filtering:
        if name in properties:
            prop = properties[name]
            default = prop.get("default")
            desc = prop.get("description", "")
            value = _format_yaml_value(default)
            lines.append(f"# {name}: {value}  # {desc}")
    lines.append("")

    # Tool-specific settings (commented out)
    lines.append("# === Tool-Specific Settings ===")
    tool_specific = [
        "min_gottcha_reads",
        "max_blast_targets",
        "blast_retention_count",
        "min_consecutive_bases",
        "qtrim",
        "include_children",
    ]
    for name in tool_specific:
        if name in properties:
            prop = properties[name]
            default = prop.get("default")
            desc = prop.get("description", "")
            value = _format_yaml_value(default)
            lines.append(f"# {name}: {value}  # {desc}")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _generate_json_template(path: Path, schema: dict, schema_url: str) -> None:
    """Generate JSON template."""
    # Include commonly-used parameters with sensible defaults
    template = {
        "$schema": schema_url,
        "samplesheet": None,
        "experiment_id": None,
        "tools": None,
        "results": None,
        "cutoff_percent": 0.001,
        "entropy": 0.9,
        "tax_stringency": 0.7,
        "preprocess": False,
    }

    with open(path, "w", encoding="utf-8") as f:
        json.dump(template, f, indent=2)
        f.write("\n")
