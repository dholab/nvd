"""Refactor tripwires for Nextflow-to-repository references.

These tests are intentionally not a Nextflow parser. They catch a narrow class
of high-cost refactor mistakes: active workflow source invoking a repository
Python script that has been deleted or left syntactically broken.
"""

from __future__ import annotations

import py_compile
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

NEXTFLOW_SOURCE_GLOBS = (
    "main.nf",
    "modules/*.nf",
    "subworkflows/*.nf",
    "workflows/*.nf",
)

# External executables with .py names. Keep this list tiny; if it grows, prefer
# a real integration smoke test over turning this file into a static analyzer.
EXTERNAL_PY_COMMANDS = frozenset({"spades.py"})

BLOCK_COMMENT_RE = re.compile(r"/\*.*?\*/", re.DOTALL)
PYTHON_COMMAND_RE = re.compile(
    r"^\s*(?:python(?:3)?\s+)?(?P<script>[A-Za-z_][\w.-]*\.py)\b",
)


def active_nextflow_files() -> list[Path]:
    """Return active Nextflow source files checked by this tripwire."""
    paths: set[Path] = set()
    for pattern in NEXTFLOW_SOURCE_GLOBS:
        paths.update(ROOT.glob(pattern))
    return sorted(paths)


def referenced_python_commands(path: Path) -> list[tuple[int, str]]:
    """Return command-position Python scripts referenced by a Nextflow file."""
    text = path.read_text(encoding="utf-8")
    text = BLOCK_COMMENT_RE.sub("", text)

    references: list[tuple[int, str]] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.lstrip()
        if not stripped or stripped.startswith(("//", "#")):
            continue

        match = PYTHON_COMMAND_RE.match(line)
        if match is None:
            continue

        script = match.group("script")
        if script in EXTERNAL_PY_COMMANDS:
            continue

        references.append((line_number, script))
    return references


def test_nextflow_referenced_python_scripts_exist_and_compile(
    tmp_path: Path,
) -> None:
    """Active Nextflow Python script calls should point at valid bin scripts."""
    missing: list[str] = []
    referenced_scripts: set[Path] = set()

    for nextflow_file in active_nextflow_files():
        for line_number, script in referenced_python_commands(nextflow_file):
            script_path = ROOT / "bin" / script
            if script_path.is_file():
                referenced_scripts.add(script_path)
            else:
                missing.append(
                    f"{nextflow_file.relative_to(ROOT)}:{line_number}: {script}",
                )

    assert not missing, (
        "Nextflow references Python scripts that do not exist in bin/:\n"
        + "\n".join(missing)
    )

    for script_path in sorted(referenced_scripts):
        py_compile.compile(
            str(script_path),
            cfile=str(tmp_path / f"{script_path.stem}.pyc"),
            doraise=True,
        )
