"""Focused Nextflow contract for BLAST risk-group annotation."""

from __future__ import annotations

import os
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def write_executable(path: Path, source: str) -> None:
    path.write_text(source, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_blast_risk_annotation_preserves_query_class(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_executable(
        bin_dir / "annotate_risk_groups.py",
        """#!/bin/sh
input=''
output=''
while [ "$#" -gt 0 ]; do
    case "$1" in
        --input) input="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        *) shift ;;
    esac
done
cp "$input" "$output"
""",
    )
    blast_tsv = tmp_path / "sample_A.short_assembly_contig.blast.merged_with_lca.tsv"
    blast_tsv.write_text("sample\tqseqid\ttaxid\nsample_A\tquery_1\t10244\n")
    risk_group_lookup = tmp_path / "risk_groups.tsv"
    risk_group_lookup.write_text("taxid\twho_risk_group\n10244\tRG2\n")
    taxonomy_dir = tmp_path / "taxonomy"
    taxonomy_dir.mkdir()

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ ANNOTATE_BLAST_RISK_GROUPS }} from '{ROOT / "modules" / "risk_groups"}'

workflow {{
    ANNOTATE_BLAST_RISK_GROUPS(
        Channel.of(tuple('sample_A', 'short_assembly_contig', file('{blast_tsv}'))),
        Channel.value(file('{risk_group_lookup}')),
        Channel.value('{taxonomy_dir}'),
    )
    ANNOTATE_BLAST_RISK_GROUPS.out.view {{ sample_id, query_class, annotated ->
        "RISK_RESULT: ${{sample_id}}|${{query_class}}|${{annotated.name}}"
    }}
}}
""",
        encoding="utf-8",
    )

    nextflow = shutil.which("nextflow")
    assert nextflow is not None
    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    environment["PATH"] = f"{bin_dir}{os.pathsep}{ROOT / 'bin'}{os.pathsep}{environment['PATH']}"
    completed = subprocess.run(  # noqa: S603
        [nextflow, "-C", "/dev/null", "run", str(workflow)],
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert (
        "RISK_RESULT: "
        "sample_A|short_assembly_contig|"
        "sample_A.short_assembly_contig.blast.with_risk_groups.tsv"
    ) in completed.stdout
