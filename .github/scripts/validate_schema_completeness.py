#!/usr/bin/env python3
"""
Validate that the NVD params JSON schema covers all configurable parameters
from nextflow.config AND that the Pydantic NvdParams model stays in sync.

This script is run in CI to ensure:
1. JSON schema covers all configurable params from nextflow.config
2. NvdParams model fields match JSON schema properties
3. CLI run options match NvdParams fields

Usage:
    python .github/scripts/validate_schema_completeness.py

Exit codes:
    0: All checks pass
    1: Schema is incomplete or NvdParams is out of sync
"""

from __future__ import annotations

import ast
import json
import re
import sys
from pathlib import Path

# Parameters that are internal implementation details and should NOT be in schema
# (from nextflow.config perspective)
INTERNAL_PARAMS_NEXTFLOW = {
    # Runtime/build-time values
    "date",
    "resources",
    # Version tracking (managed separately)
    "gottcha2_db_version",
    "blast_db_version",
    # Complex types not suitable for simple schema validation
    "human_virus_families",
    # Internal file paths
    "nvd_files",
}

# Parameters in NvdParams that are intentionally NOT in the JSON schema
# (internal runtime params that users shouldn't set via params file)
PYDANTIC_ONLY_PARAMS = {
    "date",
    "resources",
}

# Parameters in JSON schema that are intentionally NOT in NvdParams
# (currently none, but placeholder for future exceptions)
SCHEMA_ONLY_PARAMS: set[str] = set()

# CLI options that are NOT NvdParams fields (Nextflow-native or CLI-only options)
CLI_ONLY_OPTIONS = {
    # Typer internals
    "ctx",  # Typer context object
    # Nextflow-native options (not pipeline params)
    "profile",
    "config",
    "resume",
    # CLI-only options
    "params_file",  # Handled specially by CLI, merged into NvdParams
    "preset",  # Resolved to params before passing to Nextflow
    "dry_run",  # CLI behavior flag
    "no_slack",  # Negated form of slack_enabled
}

# NvdParams fields that don't need CLI options
# (set via params-file or preset only, or internal)
PYDANTIC_NO_CLI_NEEDED = {
    # Internal params
    "date",
    "resources",
    "refman_registry",
    "monoimage",
    # Version tracking (set in config, not CLI)
    "gottcha2_db_version",
    "blast_db_version",
    # Complex types (set via params-file)
    "human_virus_families",
    # Typically set via params-file or preset
    "nvd_files",
    "state_dir",
    # Exposed via negated CLI flag --no-slack
    "slack_enabled",
}


def extract_nextflow_params(config_path: Path) -> dict[str, dict]:
    """
    Parse nextflow.config and extract all params from params { } block.

    Returns:
        dict mapping parameter name to info dict with keys:
        - value: the default value as a string
        - line: line number in config file
        - category: "configurable" or "internal"
    """
    params = {}

    with open(config_path) as f:
        lines = f.readlines()

    in_params_block = False
    brace_depth = 0
    preceding_comment = ""

    for line_num, line in enumerate(lines, start=1):
        stripped = line.strip()

        # Track comments
        if stripped.startswith("//"):
            preceding_comment += stripped + "\n"
            continue

        # Detect params block start
        if not in_params_block and re.match(r"params\s*\{", stripped):
            in_params_block = True
            brace_depth = 1
            continue

        if not in_params_block:
            continue

        # Track brace depth
        brace_depth += stripped.count("{") - stripped.count("}")

        # Exit params block
        if brace_depth == 0:
            break

        # Parse parameter line: name = value
        match = re.match(r"(\w+)\s*=\s*(.+)", stripped)
        if match:
            param_name = match.group(1)
            param_value = match.group(2).strip()

            # Remove trailing comments
            param_value = re.sub(r"\s*//.*$", "", param_value)

            # Categorize
            category = categorize_param(param_name, param_value, preceding_comment)

            params[param_name] = {
                "value": param_value,
                "line": line_num,
                "category": category,
            }

            # Reset comment accumulator
            preceding_comment = ""

    return params


def categorize_param(name: str, value: str, preceding_comment: str) -> str:
    """
    Determine if a parameter should be in the schema.

    Returns:
        "internal" or "configurable"
    """
    # Explicit internal list
    if name in INTERNAL_PARAMS_NEXTFLOW:
        return "internal"

    # Version tracking params
    if name.endswith("_version"):
        return "internal"

    # Complex types (lists)
    if value.strip().startswith("["):
        return "internal"

    # Java expressions
    if "new java" in value:
        return "internal"

    # Project/launch directory paths (internal paths)
    if "${projectDir}" in value:
        return "internal"

    # Check for "SHOULD NOT BE CHANGED" in comments
    if "SHOULD NOT BE CHANGED" in preceding_comment:
        return "internal"

    return "configurable"


def load_schema(schema_path: Path) -> dict:
    """Load and parse the JSON schema."""
    with open(schema_path) as f:
        return json.load(f)


def get_pydantic_fields() -> set[str]:
    """
    Get all field names from NvdParams model.

    Returns:
        Set of field names defined in NvdParams
    """
    # Import here to avoid import errors if py_nvd isn't installed
    from py_nvd.models import NvdParams

    return set(NvdParams.model_fields.keys())


def compare_pydantic_to_schema(
    pydantic_fields: set[str],
    schema: dict,
) -> tuple[list[str], list[str]]:
    """
    Compare NvdParams fields against schema properties.

    Returns:
        (missing_from_schema, missing_from_pydantic)
        - missing_from_schema: Pydantic fields not in JSON schema
        - missing_from_pydantic: Schema properties not in Pydantic model
    """
    schema_properties = set(schema.get("properties", {}).keys())

    # Fields in Pydantic but not in schema (excluding known exceptions)
    pydantic_only = pydantic_fields - schema_properties - PYDANTIC_ONLY_PARAMS
    missing_from_schema = sorted(pydantic_only)

    # Fields in schema but not in Pydantic (excluding known exceptions)
    schema_only = schema_properties - pydantic_fields - SCHEMA_ONLY_PARAMS
    missing_from_pydantic = sorted(schema_only)

    return missing_from_schema, missing_from_pydantic


def extract_cli_options(run_py_path: Path) -> set[str]:
    """
    Extract CLI option names from run.py using AST parsing.

    Returns:
        Set of parameter names from the run() function signature
    """
    with open(run_py_path) as f:
        source = f.read()

    tree = ast.parse(source)
    options = set()

    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        if node.name != "run":
            continue

        # Get all argument names from the run() function
        for arg in node.args.args:
            options.add(arg.arg)

        break  # Found run(), done

    return options


def compare_cli_to_pydantic(
    cli_options: set[str],
    pydantic_fields: set[str],
) -> tuple[list[str], list[str]]:
    """
    Compare CLI options against NvdParams fields.

    Returns:
        (cli_missing_from_pydantic, pydantic_missing_from_cli)
        - cli_missing_from_pydantic: CLI options not in NvdParams
        - pydantic_missing_from_cli: NvdParams fields without CLI options
    """
    # CLI options that should be in NvdParams
    cli_for_pydantic = cli_options - CLI_ONLY_OPTIONS
    cli_missing = sorted(cli_for_pydantic - pydantic_fields)

    # NvdParams fields that should have CLI options
    pydantic_for_cli = pydantic_fields - PYDANTIC_NO_CLI_NEEDED
    pydantic_missing = sorted(pydantic_for_cli - cli_options)

    return cli_missing, pydantic_missing


def compare_params_to_schema(
    nf_params: dict[str, dict],
    schema: dict,
) -> tuple[list[str], list[str], list[str]]:
    """
    Compare Nextflow params against schema properties.

    Returns:
        (missing, extra, internal)
        - missing: Configurable params not in schema
        - extra: Schema properties not in nextflow.config (informational)
        - internal: Internal params correctly excluded from schema
    """
    schema_properties = set(schema.get("properties", {}).keys())

    configurable = {
        name for name, info in nf_params.items() if info["category"] == "configurable"
    }
    internal = {
        name for name, info in nf_params.items() if info["category"] == "internal"
    }

    # Missing from schema
    missing = sorted(configurable - schema_properties)

    # In schema but not in nextflow.config (might be future params or typos)
    all_nf_params = set(nf_params.keys())
    extra = sorted(schema_properties - all_nf_params)

    return missing, extra, sorted(internal)


def main() -> int:
    """Main entry point."""
    # Locate files relative to repo root
    repo_root = Path(__file__).parent.parent.parent
    config_path = repo_root / "nextflow.config"
    schema_path = repo_root / "schemas" / "nvd-params.latest.schema.json"

    print("=" * 70)
    print("NVD Schema Completeness Validation")
    print("=" * 70)
    print()

    # Check files exist
    if not config_path.exists():
        print(f"ERROR: nextflow.config not found at {config_path}")
        return 1

    if not schema_path.exists():
        print(f"ERROR: Schema not found at {schema_path}")
        return 1

    print(f"Config:  {config_path}")
    print(f"Schema:  {schema_path}")
    print()

    # Extract and compare
    print("Parsing nextflow.config...")
    nf_params = extract_nextflow_params(config_path)

    print("Loading schema...")
    schema = load_schema(schema_path)

    print()

    # Track overall pass/fail
    all_passed = True

    # =========================================================================
    # Check 1: nextflow.config → JSON schema coverage
    # =========================================================================
    print("-" * 70)
    print("CHECK 1: nextflow.config → JSON schema")
    print("-" * 70)
    print()

    missing_nf, extra_schema, internal = compare_params_to_schema(nf_params, schema)

    configurable_count = sum(
        1 for info in nf_params.values() if info["category"] == "configurable"
    )
    schema_count = len(schema.get("properties", {}))

    print(f"Configurable params in nextflow.config: {configurable_count}")
    print(f"Properties in schema:                   {schema_count}")
    print(f"Internal params (correctly excluded):   {len(internal)}")
    print()

    if missing_nf:
        print("MISSING FROM SCHEMA:")
        print()
        for param in missing_nf:
            info = nf_params[param]
            print(f"  - {param} (line {info['line']})")
            print(f"    Default: {info['value']}")
        print()
        print("Add these parameters to schemas/nvd-params.*.schema.json")
        print()
        all_passed = False
    else:
        print("✓ All configurable nextflow.config params are in schema")
        print()

    if extra_schema:
        print("EXTRA IN SCHEMA (informational):")
        print()
        for param in extra_schema:
            print(f"  - {param}")
        print()

    # =========================================================================
    # Check 2: NvdParams ↔ JSON schema sync
    # =========================================================================
    print("-" * 70)
    print("CHECK 2: NvdParams ↔ JSON schema")
    print("-" * 70)
    print()

    try:
        pydantic_fields = get_pydantic_fields()
        print(f"Fields in NvdParams:    {len(pydantic_fields)}")
        print(f"Properties in schema:   {schema_count}")
        print()

        missing_from_schema, missing_from_pydantic = compare_pydantic_to_schema(
            pydantic_fields, schema
        )

        if missing_from_schema:
            print("PYDANTIC FIELDS MISSING FROM SCHEMA:")
            print()
            for param in missing_from_schema:
                print(f"  - {param}")
            print()
            print("Add these to schemas/nvd-params.*.schema.json")
            print()
            all_passed = False
        else:
            print("✓ All NvdParams fields are in schema (or correctly excluded)")

        if missing_from_pydantic:
            print()
            print("SCHEMA PROPERTIES MISSING FROM NVDPARAMS:")
            print()
            for param in missing_from_pydantic:
                print(f"  - {param}")
            print()
            print("Add these to lib/py_nvd/models.py NvdParams class")
            print()
            all_passed = False
        else:
            print("✓ All schema properties are in NvdParams (or correctly excluded)")

        print()

        # =====================================================================
        # Check 3: CLI options ↔ NvdParams fields
        # =====================================================================
        print("-" * 70)
        print("CHECK 3: CLI options ↔ NvdParams fields")
        print("-" * 70)
        print()

        run_py_path = repo_root / "lib" / "py_nvd" / "cli" / "commands" / "run.py"
        if not run_py_path.exists():
            print(f"WARNING: run.py not found at {run_py_path}")
            print("         Skipping CLI sync check")
            print()
        else:
            cli_options = extract_cli_options(run_py_path)
            print(f"CLI options in run():   {len(cli_options)}")
            print(f"Fields in NvdParams:    {len(pydantic_fields)}")
            print()

            cli_missing, pydantic_missing = compare_cli_to_pydantic(
                cli_options, pydantic_fields
            )

            if cli_missing:
                print("CLI OPTIONS NOT IN NVDPARAMS:")
                print()
                for opt in cli_missing:
                    print(f"  - {opt}")
                print()
                print("Add these to lib/py_nvd/models.py NvdParams class,")
                print("or add to CLI_ONLY_OPTIONS if they're intentionally CLI-only")
                print()
                all_passed = False
            else:
                print(
                    "✓ All CLI options map to NvdParams fields (or correctly excluded)"
                )

            if pydantic_missing:
                print()
                print("NVDPARAMS FIELDS WITHOUT CLI OPTIONS:")
                print()
                for field in pydantic_missing:
                    print(f"  - {field}")
                print()
                print("Add CLI options to lib/py_nvd/cli/commands/run.py,")
                print("or add to PYDANTIC_NO_CLI_NEEDED if CLI option not needed")
                print()
                all_passed = False
            else:
                print("✓ All NvdParams fields have CLI options (or correctly excluded)")

            print()

    except ImportError as e:
        print(f"WARNING: Could not import NvdParams: {e}")
        print("         Skipping Pydantic and CLI sync checks")
        print()

    # =========================================================================
    # Final verdict
    # =========================================================================
    print("=" * 70)
    if all_passed:
        print("PASSED: All schema completeness checks passed")
        return 0
    else:
        print("FAILED: Schema completeness checks failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
