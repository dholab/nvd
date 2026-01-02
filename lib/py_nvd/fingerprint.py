"""
Pipeline fingerprint generation and verification.

This module provides functions to generate and verify blake3 hashes of
critical pipeline files (main.nf, nextflow.config) to ensure the CLI
is operating on the correct NVD pipeline, not some other Nextflow project.

The fingerprint is stored in lib/py_nvd/_fingerprint.json and should be
regenerated whenever main.nf or nextflow.config change (automated via
pre-commit hook or pixi task).
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import blake3

# Path to the fingerprint file (relative to this module)
FINGERPRINT_FILE = Path(__file__).parent / "_fingerprint.json"

# Files to include in the fingerprint
FINGERPRINTED_FILES = ["main.nf", "nextflow.config"]


def hash_file(path: Path) -> str:
    """
    Compute blake3 hash of a file.

    Args:
        path: Path to the file to hash.

    Returns:
        Hex-encoded blake3 hash string.
    """
    assert path.exists(), f"Cannot hash non-existent file: {path}"
    assert path.is_file(), f"Expected file, got directory: {path}"

    hasher = blake3.blake3()
    with open(path, "rb") as f:
        # Read in chunks for memory efficiency on large files
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)

    result = hasher.hexdigest()
    assert isinstance(result, str), "blake3 hexdigest should return str"
    assert len(result) == 64, f"blake3 hex digest should be 64 chars, got {len(result)}"
    return result


def generate_fingerprint(pipeline_root: Path) -> dict[str, str]:
    """
    Generate fingerprint hashes for pipeline files.

    Args:
        pipeline_root: Path to the NVD pipeline root directory.

    Returns:
        Dictionary mapping filenames to their blake3 hashes.

    Raises:
        FileNotFoundError: If any fingerprinted file is missing.
    """
    assert pipeline_root.is_dir(), f"Pipeline root must be a directory: {pipeline_root}"

    fingerprint = {}
    for filename in FINGERPRINTED_FILES:
        file_path = pipeline_root / filename
        if not file_path.exists():
            msg = (
                f"Cannot generate fingerprint: {filename} not found in {pipeline_root}"
            )
            raise FileNotFoundError(msg)
        fingerprint[filename] = hash_file(file_path)

    assert len(fingerprint) == len(FINGERPRINTED_FILES), "All files should be hashed"
    return fingerprint


def save_fingerprint(pipeline_root: Path, output_path: Path | None = None) -> Path:
    """
    Generate and save fingerprint to JSON file.

    Args:
        pipeline_root: Path to the NVD pipeline root directory.
        output_path: Where to save the fingerprint. Defaults to FINGERPRINT_FILE.

    Returns:
        Path to the saved fingerprint file.
    """
    output = output_path or FINGERPRINT_FILE
    fingerprint = generate_fingerprint(pipeline_root)

    assert output.parent.exists(), f"Output directory must exist: {output.parent}"

    with open(output, "w") as f:
        json.dump(fingerprint, f, indent=2)
        f.write("\n")  # Trailing newline

    assert output.exists(), f"Fingerprint file should exist after save: {output}"
    return output


def load_fingerprint(fingerprint_path: Path | None = None) -> dict[str, str] | None:
    """
    Load fingerprint from JSON file.

    Args:
        fingerprint_path: Path to fingerprint file. Defaults to FINGERPRINT_FILE.

    Returns:
        Dictionary of filename -> hash, or None if file doesn't exist.
    """
    path = fingerprint_path or FINGERPRINT_FILE
    if not path.exists():
        return None

    with open(path) as f:
        result = json.load(f)

    assert isinstance(result, dict), f"Fingerprint should be a dict, got {type(result)}"
    return result


def verify_pipeline(
    candidate: Path,
    fingerprint_path: Path | None = None,
    strict: bool = True,
) -> bool:
    """
    Verify a candidate directory is the authentic NVD pipeline.

    Args:
        candidate: Path to check.
        fingerprint_path: Path to fingerprint file. Defaults to FINGERPRINT_FILE.
        strict: If True, require fingerprint file to exist and all hashes to match.
                If False, return True when fingerprint file is missing (dev mode).

    Returns:
        True if the candidate is verified as the NVD pipeline.
    """
    assert candidate.is_dir(), f"Candidate must be a directory: {candidate}"

    expected = load_fingerprint(fingerprint_path)

    if expected is None:
        # No fingerprint file - in strict mode this is a failure,
        # in non-strict mode we assume dev environment and allow it
        return not strict

    return all(
        (candidate / filename).exists()
        and hash_file(candidate / filename) == expected_hash
        for filename, expected_hash in expected.items()
    )


def is_dev_mode() -> bool:
    """
    Check if running in development mode.

    Development mode is detected by:
    1. NVD_DEV_MODE environment variable is set
    2. The package is installed in editable mode (py_nvd source is in a git repo)

    Returns:
        True if in development mode.
    """
    if os.environ.get("NVD_DEV_MODE"):
        return True

    # Check if we're in an editable install (source dir has .git)
    source_dir = Path(__file__).parent.parent.parent  # lib/py_nvd -> lib -> repo root
    if (source_dir / ".git").exists():
        return True

    return False


def _find_pipeline_root_for_fingerprint() -> Path | None:
    """
    Locate pipeline root for fingerprint generation.

    Tries (in order):
    1. Walk up from this file's location
    2. Current working directory

    Returns:
        Path to pipeline root, or None if not found.
    """
    # Try walking up from this file
    current = Path(__file__).parent.parent.parent  # lib/py_nvd -> lib -> repo root
    if (current / "main.nf").exists():
        return current

    # Try cwd
    cwd = Path.cwd()
    if (cwd / "main.nf").exists():
        return cwd

    return None


def main() -> None:
    """
    CLI entry point for generating/updating the pipeline fingerprint.

    Usage:
        nvd-fingerprint          # Generate from repo root (auto-detected)
        nvd-fingerprint /path    # Generate from specified path
    """
    # Determine pipeline root - explicit path takes precedence
    if len(sys.argv) > 1:
        pipeline_root = Path(sys.argv[1]).resolve()
        if not pipeline_root.is_dir():
            print(f"Error: Not a directory: {pipeline_root}")
            sys.exit(1)
        if not (pipeline_root / "main.nf").exists():
            print(f"Error: No main.nf found in {pipeline_root}")
            print("Are you sure this is the NVD pipeline root?")
            sys.exit(1)
    else:
        pipeline_root = _find_pipeline_root_for_fingerprint()
        if pipeline_root is None:
            print("Error: Could not find pipeline root (main.nf not found)")
            print()
            print("Try one of:")
            print("  - Run from the NVD repository root directory")
            print("  - Pass the path explicitly: nvd-fingerprint /path/to/nvd")
            sys.exit(1)

    try:
        output = save_fingerprint(pipeline_root)
        fingerprint = load_fingerprint(output)
        print(f"Fingerprint saved to: {output}")
        print("Hashes:")
        for filename, file_hash in (fingerprint or {}).items():
            print(f"  {filename}: {file_hash[:16]}...")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
