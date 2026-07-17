"""Refactor tripwires for Nextflow-to-repository references.

These tests are intentionally not a Nextflow parser. They catch a narrow class
of high-cost refactor mistakes: active workflow source invoking a repository
Python script that has been deleted or left syntactically broken.
"""

from __future__ import annotations

import py_compile
import re
import subprocess
import sys
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


def test_resolve_read_inputs_starts_with_package_imports() -> None:
    """The read resolver script should start without source-tree path shims."""
    result = subprocess.run(  # noqa: S603
        [sys.executable, str(ROOT / "bin" / "resolve_read_inputs.py"), "--help"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0
    assert "Resolve NVD samplesheet read declarations" in result.stdout


def test_taxonomy_processes_forward_policy_arguments() -> None:
    """Taxonomy preflight and annotation processes should receive policy controls."""
    expected = {
        ROOT / "modules" / "utils.nf": ("ensure_taxonomy.py", "annotate_blast_lca.py"),
        ROOT / "modules" / "blast.nf": ("annotate_blast_results.py",),
    }

    missing: list[str] = []
    for path, commands in expected.items():
        text = path.read_text(encoding="utf-8")
        missing.extend(
            f"{path.relative_to(ROOT)}: missing {command}"
            for command in commands
            if command not in text
        )
        if "--taxonomy-mode" not in text:
            missing.append(f"{path.relative_to(ROOT)}: missing --taxonomy-mode")
        if "--taxonomy-max-age-days" not in text:
            missing.append(f"{path.relative_to(ROOT)}: missing --taxonomy-max-age-days")

    assert not missing, "\n".join(missing)


def test_skip_stage_params_gate_expensive_nextflow_processes() -> None:
    """Skip params should guard scheduling before SPAdes and MEGABLAST."""
    expected = {
        ROOT / "subworkflows" / "preprocess_contigs.nf": (
            "params.skip_assembly",
            "RUN_SPADES(",
        ),
        ROOT / "subworkflows" / "classify_with_megablast.nf": (
            "params.skip_blast",
            "MEGABLAST(",
        ),
    }

    missing: list[str] = []
    for path, fragments in expected.items():
        text = path.read_text(encoding="utf-8")
        missing.extend(
            f"{path.relative_to(ROOT)}: missing {fragment}"
            for fragment in fragments
            if fragment not in text
        )

    main_text = (ROOT / "workflows" / "nvd_main.nf").read_text(encoding="utf-8")
    missing.extend(
        f"workflows/nvd_main.nf: missing {fragment}"
        for fragment in (
            "requires_blast_db",
            "params.skip_assembly",
            "params.skip_blast",
        )
        if fragment not in main_text
    )

    assert not missing, "\n".join(missing)


def test_sourmash_reports_do_not_wait_for_runtime_taxonomy_normalization() -> None:
    """Per-sample sourmash reports should stream from each summary CSV."""
    reporting_module = (ROOT / "modules" / "reporting.nf").read_text(
        encoding="utf-8",
    )
    reporting_subworkflow = (ROOT / "subworkflows" / "reporting.nf").read_text(
        encoding="utf-8",
    )

    forbidden_fragments = (
        "BUILD_TAXONOMIC_PROFILE_NORMALIZATION_MAP",
        "NORMALIZE_TAXONOMIC_PROFILE_SUMMARY",
        "normalize_taxonomic_profile_summary.py",
        "taxonomic_profile.summary.normalized.csv",
    )
    for fragment in forbidden_fragments:
        assert fragment not in reporting_module
        assert fragment not in reporting_subworkflow

    assert "ch_sourmash_profile_summaries = ch_sourmash_tax_reports" in (
        reporting_subworkflow
    )
    assert "RENDER_TAXON_ABUNDANCE_SUNBURST(ch_sourmash_profile_summaries)" in (
        reporting_subworkflow
    )
    assert "RENDER_SOURMASH_SANKEY(ch_sourmash_profile_summaries)" in (
        reporting_subworkflow
    )
    assert "ch_merged_taxburst_input = ch_sourmash_profile_summaries" in (
        reporting_subworkflow
    )
