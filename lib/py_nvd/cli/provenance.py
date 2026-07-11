"""Resolve source identity for Nextflow run metadata."""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

SOURCE_REVISION_ENV = "NVD_SOURCE_REVISION"
SOURCE_DIRTY_ENV = "NVD_SOURCE_DIRTY"
VCS_TIMEOUT_SECONDS = 5


@dataclass(frozen=True)
class SourceRevision:
    """A source snapshot identifier and its working-tree state."""

    revision: str | None
    dirty: bool | None


def _capture(command: list[str]) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(  # noqa: S603 - commands are fixed by this module
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=VCS_TIMEOUT_SECONDS,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None


def _resolve_jj(pipeline_root: Path) -> SourceRevision | None:
    if not (pipeline_root / ".jj").exists() or shutil.which("jj") is None:
        return None
    result = _capture(
        [
            "jj",
            "-R",
            str(pipeline_root),
            "log",
            "-r",
            "@",
            "--no-graph",
            "-T",
            "commit_id",
        ],
    )
    revision = result.stdout.strip() if result and result.returncode == 0 else ""
    return SourceRevision(revision=revision, dirty=None) if revision else None


def _resolve_git(pipeline_root: Path) -> SourceRevision | None:
    if not (pipeline_root / ".git").exists() or shutil.which("git") is None:
        return None
    revision_result = _capture(
        ["git", "-C", str(pipeline_root), "rev-parse", "HEAD"],
    )
    revision = (
        revision_result.stdout.strip()
        if revision_result and revision_result.returncode == 0
        else ""
    )
    if not revision:
        return None
    status_result = _capture(
        ["git", "-C", str(pipeline_root), "status", "--porcelain"],
    )
    dirty = (
        bool(status_result.stdout.strip())
        if status_result and status_result.returncode == 0
        else None
    )
    return SourceRevision(revision=revision, dirty=dirty)


def resolve_source_revision(pipeline_root: Path) -> SourceRevision:
    """Resolve the current source snapshot without requiring a particular VCS."""
    pipeline_root = Path(pipeline_root)
    return (
        _resolve_jj(pipeline_root)
        or _resolve_git(pipeline_root)
        or SourceRevision(revision=None, dirty=None)
    )


def nextflow_environment(
    pipeline_root: Path,
    *,
    environ: dict[str, str] | None = None,
) -> dict[str, str]:
    """Build a Nextflow environment containing source provenance."""
    environment = dict(os.environ if environ is None else environ)
    source = resolve_source_revision(pipeline_root)
    environment[SOURCE_REVISION_ENV] = source.revision or "unknown"
    environment[SOURCE_DIRTY_ENV] = (
        str(source.dirty).lower() if source.dirty is not None else "unknown"
    )
    return environment
