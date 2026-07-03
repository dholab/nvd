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
    """Skip params should guard scheduling before assembly and MEGABLAST."""
    expected = {
        ROOT / "workflows" / "nvd_main.nf": (
            "params.skip_assembly",
            "ch_reads_by_platform",
            "SHORT_READ_DENOVO_ASSEMBLY(",
            "LONG_READ_DENOVO_ENSEMBLY(",
        ),
        ROOT / "subworkflows" / "short_read_denovo_assembly.nf": (
            "RUN_SPADES(",
        ),
        ROOT / "subworkflows" / "long_read_denovo_ensembly.nf": (
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


def test_crumbs_coverage_process_streams_zero_inclusive_depth() -> None:
    """CRUMBS coverage summarization should consume a streamed depth -aa table."""
    text = (ROOT / "modules" / "samtools.nf").read_text(encoding="utf-8")
    expected_fragments = (
        "process SUMMARIZE_CONTIG_COVERAGE",
        "samtools depth -aa",
        "summarize_contig_coverage.py",
        "--depth-tsv -",
        "emit: coverage_summary",
    )

    missing = [fragment for fragment in expected_fragments if fragment not in text]

    assert not missing, "modules/samtools.nf missing: " + ", ".join(missing)


def test_crumbs_report_exporter_consumes_only_taxon_evidence() -> None:
    """CRUMBS terminal report exports should format the taxon evidence table."""
    crumbs_module = (ROOT / "modules" / "crumbs.nf").read_text(encoding="utf-8")
    crumbs_subworkflow = (ROOT / "subworkflows" / "crumbs_profiling.nf").read_text(
        encoding="utf-8",
    )

    expected_module_fragments = (
        "process EXPORT_CRUMBS_TAXONOMIC_REPORTS",
        "tuple val(sample_id), path(taxa_tsv)",
        "--taxa-tsv ${taxa_tsv}",
        "emit: krona",
        "emit: kreport",
    )
    missing = [
        fragment
        for fragment in expected_module_fragments
        if fragment not in crumbs_module
    ]

    assert not missing, "modules/crumbs.nf missing: " + ", ".join(missing)
    assert (
        "EXPORT_CRUMBS_TAXONOMIC_REPORTS(ESTIMATE_CRUMBS_PROFILE.out.taxa)"
        in crumbs_subworkflow
    )


def test_crumbs_taxburst_renders_from_taxon_evidence() -> None:
    """CRUMBS Taxburst should retain the complete taxon evidence hierarchy."""
    crumbs_module = (ROOT / "modules" / "crumbs.nf").read_text(encoding="utf-8")
    crumbs_subworkflow = (ROOT / "subworkflows" / "crumbs_profiling.nf").read_text(
        encoding="utf-8",
    )

    expected_module_fragments = (
        "process RENDER_CRUMBS_TAXBURST",
        "tuple val(sample_id), path(taxa_tsv)",
        "render_multisample_taxburst.py",
        "--input-format crumbs",
        '--summary "${sample_id}=${taxa_tsv}"',
        "--output-json",
        "emit: reports",
    )
    missing = [
        fragment
        for fragment in expected_module_fragments
        if fragment not in crumbs_module
    ]

    assert not missing, "modules/crumbs.nf missing: " + ", ".join(missing)
    assert (
        "RENDER_CRUMBS_TAXBURST(ESTIMATE_CRUMBS_PROFILE.out.taxa)"
        in crumbs_subworkflow
    )


def test_merged_crumbs_taxburst_groups_taxon_evidence() -> None:
    """Merged CRUMBS Taxburst should retain each sample's complete hierarchy."""
    crumbs_module = (ROOT / "modules" / "crumbs.nf").read_text(encoding="utf-8")
    crumbs_subworkflow = (ROOT / "subworkflows" / "crumbs_profiling.nf").read_text(
        encoding="utf-8",
    )

    expected_module_fragments = (
        "process RENDER_MERGED_CRUMBS_TAXBURST",
        "tuple val(sample_ids), path(taxa_tsvs)",
        "render_multisample_taxburst.py",
        "--input-format crumbs",
        "--output crumbs.taxburst.html",
        "emit: report",
    )
    missing = [
        fragment
        for fragment in expected_module_fragments
        if fragment not in crumbs_module
    ]

    assert not missing, "modules/crumbs.nf missing: " + ", ".join(missing)
    assert "ch_merged_taxburst_input = ESTIMATE_CRUMBS_PROFILE.out.taxa" in (
        crumbs_subworkflow
    )
    assert "RENDER_MERGED_CRUMBS_TAXBURST(ch_merged_taxburst_input)" in (
        crumbs_subworkflow
    )
