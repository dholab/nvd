from __future__ import annotations

import subprocess
from pathlib import Path

from py_nvd.cli import provenance


def test_resolve_source_revision_prefers_jj(
    monkeypatch: object,
    tmp_path: Path,
) -> None:
    tmp_path = Path(tmp_path)
    (tmp_path / ".jj").mkdir()
    monkeypatch.setattr(provenance.shutil, "which", lambda command: f"/bin/{command}")
    monkeypatch.setattr(
        provenance,
        "_capture",
        lambda _command: subprocess.CompletedProcess([], 0, "abc123\n", ""),
    )

    assert provenance.resolve_source_revision(tmp_path) == provenance.SourceRevision(
        revision="abc123",
        dirty=None,
    )


def test_resolve_source_revision_reports_git_worktree_state(
    monkeypatch: object,
    tmp_path: Path,
) -> None:
    tmp_path = Path(tmp_path)
    (tmp_path / ".git").mkdir()
    monkeypatch.setattr(
        provenance.shutil,
        "which",
        lambda command: "/bin/git" if command == "git" else None,
    )
    results = iter(
        (
            subprocess.CompletedProcess([], 0, "def456\n", ""),
            subprocess.CompletedProcess([], 0, " M modules/utils.nf\n", ""),
        ),
    )
    monkeypatch.setattr(provenance, "_capture", lambda _command: next(results))

    assert provenance.resolve_source_revision(tmp_path) == provenance.SourceRevision(
        revision="def456",
        dirty=True,
    )


def test_nextflow_environment_uses_unknown_without_checkout(tmp_path: Path) -> None:
    tmp_path = Path(tmp_path)
    environment = provenance.nextflow_environment(
        tmp_path,
        environ={"EXISTING": "value"},
    )

    assert environment == {
        "EXISTING": "value",
        "NVD_SOURCE_REVISION": "unknown",
        "NVD_SOURCE_DIRTY": "unknown",
    }
