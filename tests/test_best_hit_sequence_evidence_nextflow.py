"""Nextflow contract tests for best-hit sequence evidence publishing."""

from __future__ import annotations

import os
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UTILS_MODULE = ROOT / "modules" / "utils"
RESULTS_CONFIG = ROOT / "conf" / "results.config"
QBT_HEADER = (
    "sample_id\tqseqid\tqlen\tbest_hit_reference_accession\t"
    "best_hit_reference_title\tbest_hit_query_start_1based\t"
    "best_hit_query_end_1based\tbest_hit_reference_length\t"
    "best_hit_reference_start_1based\tbest_hit_reference_end_1based\t"
    "best_hit_reference_strand"
)


def write_fake_blastdbcmd(bin_dir: Path, *, fail_if_called: bool = False) -> None:
    blastdbcmd = bin_dir / "blastdbcmd"
    if fail_if_called:
        blastdbcmd.write_text("#!/usr/bin/env bash\nexit 99\n", encoding="utf-8")
    else:
        blastdbcmd.write_text(
            """#!/usr/bin/env bash
set -euo pipefail
batch=""
target_only=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    -entry_batch) batch="$2"; shift 2 ;;
    -target_only) target_only=true; shift ;;
    *) shift ;;
  esac
done
if [[ "$target_only" != true ]]; then
  exit 42
fi
while IFS= read -r accession; do
  case "$accession" in
    REF_A) printf 'REF_A\tACGTACGT\n' ;;
    REF_B) printf 'REF_B\tTTTTGGGG\n' ;;
    *) exit 2 ;;
  esac
done < "$batch"
""",
            encoding="utf-8",
        )
    blastdbcmd.chmod(blastdbcmd.stat().st_mode | stat.S_IXUSR)


def run_nextflow(
    workflow: Path,
    tmp_path: Path,
    bin_dir: Path,
) -> subprocess.CompletedProcess[str]:
    config = tmp_path / "nextflow.config"
    config.write_text(
        f"""\
params.results = '{tmp_path / "results"}'
params.experimental = true
params.skip_blast = false
params.blast_db_prefix = 'mini'
params.no_enrichment = true
params.virus_index = null
params.virus_index_url = null
params.virus_reference_fasta = null
includeConfig '{RESULTS_CONFIG}'
""",
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment["PATH"] = (
        f"{bin_dir}{os.pathsep}{ROOT / 'bin'}{os.pathsep}{environment['PATH']}"
    )
    nextflow = shutil.which("nextflow")
    assert nextflow is not None
    return subprocess.run(  # noqa: S603
        [
            nextflow,
            "-C",
            str(config),
            "run",
            str(workflow),
            "-work-dir",
            str(tmp_path / "work"),
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )


def test_best_hit_sequence_evidence_process_stages_pairs_and_publishes(
    tmp_path: Path,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_blastdbcmd(bin_dir)

    qbt = tmp_path / "query_big_table.tsv"
    qbt.write_text(
        f"{QBT_HEADER}\n"
        "sample-A\tq1\t8\tREF_A\tAlpha reference\t1\t8\t8\t2\t5\tplus\n"
        "sample-B\tq1\t8\tREF_B\tBeta reference\t1\t8\t8\t7\t4\tminus\n",
        encoding="utf-8",
    )
    query_dir_a = tmp_path / "a"
    query_dir_b = tmp_path / "b"
    query_dir_a.mkdir()
    query_dir_b.mkdir()
    query_a = query_dir_a / "queries.fasta"
    query_b = query_dir_b / "queries.fasta"
    query_a.write_text(">q1\nAAAACCCC\n", encoding="utf-8")
    query_b.write_text(">q1\nGGGGTTTT\n", encoding="utf-8")
    blast_db = tmp_path / "blastdb"
    blast_db.mkdir()

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ EMIT_BEST_HIT_SEQUENCE_EVIDENCE }} from '{UTILS_MODULE}'

workflow {{
    EMIT_BEST_HIT_SEQUENCE_EVIDENCE(
        Channel.of(file('{qbt}')),
        Channel.value(['sample-A', 'sample-B']),
        Channel.of([file('{query_a}'), file('{query_b}')]),
        Channel.of(file('{blast_db}')),
    )
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, tmp_path, bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    evidence_dir = tmp_path / "results" / "nvd" / "11_best_hit_sequences"
    assert sorted(path.name for path in evidence_dir.iterdir()) == [
        "best_hit_placements.bed",
        "query_sequences.fasta",
        "selected_references.fasta",
    ]
    assert "sample-A|q1" in (evidence_dir / "query_sequences.fasta").read_text(
        encoding="utf-8",
    )
    assert "sample-B|q1" in (evidence_dir / "query_sequences.fasta").read_text(
        encoding="utf-8",
    )
    assert ">REF_A Alpha reference" in (
        evidence_dir / "selected_references.fasta"
    ).read_text(encoding="utf-8")
    assert ">REF_B Beta reference" in (
        evidence_dir / "selected_references.fasta"
    ).read_text(encoding="utf-8")


def test_best_hit_sequence_evidence_process_publishes_empty_qbt_outputs(
    tmp_path: Path,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_blastdbcmd(bin_dir, fail_if_called=True)

    qbt = tmp_path / "query_big_table.tsv"
    qbt.write_text(f"{QBT_HEADER}\n", encoding="utf-8")
    blast_db = tmp_path / "blastdb"
    blast_db.mkdir()
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ EMIT_BEST_HIT_SEQUENCE_EVIDENCE }} from '{UTILS_MODULE}'

workflow {{
    EMIT_BEST_HIT_SEQUENCE_EVIDENCE(
        Channel.of(file('{qbt}')),
        Channel.value([]),
        Channel.of([]),
        Channel.of(file('{blast_db}')),
    )
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, tmp_path, bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    evidence_dir = tmp_path / "results" / "nvd" / "11_best_hit_sequences"
    assert (evidence_dir / "query_sequences.fasta").read_text(encoding="utf-8") == ""
    assert (evidence_dir / "selected_references.fasta").read_text(
        encoding="utf-8",
    ) == ""
    assert (evidence_dir / "best_hit_placements.bed").read_text(encoding="utf-8") == ""
