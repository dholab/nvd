"""Tests for best-hit sequence evidence emission."""

from __future__ import annotations

import csv
import stat
from typing import TYPE_CHECKING

import pytest
from emit_best_hit_sequence_evidence import BestHitSequenceEvidenceError, main

if TYPE_CHECKING:
    from pathlib import Path


QBT_COLUMNS = [
    "sample_id",
    "qseqid",
    "qlen",
    "best_hit_reference_accession",
    "best_hit_reference_title",
    "best_hit_query_start_1based",
    "best_hit_query_end_1based",
    "best_hit_reference_length",
    "best_hit_reference_start_1based",
    "best_hit_reference_end_1based",
    "best_hit_reference_strand",
]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=QBT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_fake_blastdbcmd(path: Path, references: dict[str, str]) -> Path:
    script = path / "blastdbcmd"
    reference_cases = "\n".join(
        f"    {accession}) printf '%s\\t%s\\n' '{accession}' '{sequence}' ;;"
        for accession, sequence in references.items()
    )
    script.write_text(
        f"""#!/usr/bin/env bash
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
while IFS= read -r entry; do
  case "$entry" in
{reference_cases}
    *) exit 2 ;;
  esac
done < "$batch"
""",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def qbt_row(**overrides: str) -> dict[str, str]:
    row = {
        "sample_id": "sample-A",
        "qseqid": "q1",
        "qlen": "8",
        "best_hit_reference_accession": "REF_A",
        "best_hit_reference_title": "Alpha reference",
        "best_hit_query_start_1based": "1",
        "best_hit_query_end_1based": "8",
        "best_hit_reference_length": "16",
        "best_hit_reference_start_1based": "3",
        "best_hit_reference_end_1based": "10",
        "best_hit_reference_strand": "plus",
    }
    row.update(overrides)
    return row


def run_evidence(
    tmp_path: Path,
    *,
    qbt_rows: list[dict[str, str]],
    query_fastas: list[tuple[str, str]],
    references: dict[str, str],
) -> Path:
    qbt = tmp_path / "query_big_table.tsv"
    write_tsv(qbt, qbt_rows)
    output_dir = tmp_path / "out"
    blastdbcmd = write_fake_blastdbcmd(tmp_path, references)
    args = [
        "--query-big-table",
        str(qbt),
        "--blast-db",
        str(tmp_path / "blastdb"),
        "--blast-db-prefix",
        "mini",
        "--output-dir",
        str(output_dir),
        "--blastdbcmd-bin",
        str(blastdbcmd),
    ]
    for index, (sample_id, text) in enumerate(query_fastas, start=1):
        fasta = tmp_path / f"queries_{index}.fasta"
        fasta.write_text(text, encoding="utf-8")
        args.extend(["--query-fasta", sample_id, str(fasta)])

    main(args)
    return output_dir


def test_emits_literal_unwrapped_fastas_and_bed_with_deduplicated_references(
    tmp_path: Path,
) -> None:
    rows = [
        qbt_row(
            sample_id="sample-B",
            qseqid="shared_query",
            qlen="8",
            best_hit_reference_accession="REF_B",
            best_hit_reference_title="Beta reference",
            best_hit_reference_length="20",
            best_hit_reference_start_1based="20",
            best_hit_reference_end_1based="17",
            best_hit_reference_strand="minus",
        ),
        qbt_row(
            qseqid="shared_query",
            best_hit_reference_start_1based="3",
            best_hit_reference_end_1based="6",
        ),
        qbt_row(
            qseqid="unique_query",
            qlen="6",
            best_hit_query_end_1based="6",
            best_hit_reference_start_1based="10",
            best_hit_reference_end_1based="14",
        ),
    ]
    output_dir = run_evidence(
        tmp_path,
        qbt_rows=rows,
        query_fastas=[
            (
                "sample-A",
                ">attempted_only\nNNNN\n>attempted_only\nCCCC\n>shared_query attempted but represented\nAAAACCCC\n>unique_query\nGG\nGGTT\n",
            ),
            ("sample-B", ">shared_query\nTTTTGGGG\n"),
        ],
        references={"REF_A": "ACGTACGTACGTACGT", "REF_B": "TTTTCCCCAAAAGGGGNNNN"},
    )

    assert (output_dir / "query_sequences.fasta").read_text(encoding="utf-8") == (
        ">sample-A|shared_query\nAAAACCCC\n"
        ">sample-A|unique_query\nGGGGTT\n"
        ">sample-B|shared_query\nTTTTGGGG\n"
    )
    assert (output_dir / "selected_references.fasta").read_text(encoding="utf-8") == (
        ">REF_A Alpha reference\nACGTACGTACGTACGT\n"
        ">REF_B Beta reference\nTTTTCCCCAAAAGGGGNNNN\n"
    )
    assert (output_dir / "best_hit_placements.bed").read_text(encoding="utf-8") == (
        "REF_A\t2\t6\tsample-A|shared_query\t0\t+\n"
        "REF_A\t9\t14\tsample-A|unique_query\t0\t+\n"
        "REF_B\t16\t20\tsample-B|shared_query\t0\t-\n"
    )


def test_header_only_qbt_emits_three_empty_files_without_scanning_inputs(
    tmp_path: Path,
) -> None:
    qbt = tmp_path / "query_big_table.tsv"
    write_tsv(qbt, [])
    failing_blastdbcmd = tmp_path / "blastdbcmd"
    failing_blastdbcmd.write_text("#!/usr/bin/env bash\nexit 99\n", encoding="utf-8")
    failing_blastdbcmd.chmod(failing_blastdbcmd.stat().st_mode | stat.S_IXUSR)
    missing_query_fasta = tmp_path / "missing.fasta"
    output_dir = tmp_path / "out"

    main(
        [
            "--query-big-table",
            str(qbt),
            "--query-fasta",
            "sample-A",
            str(missing_query_fasta),
            "--blast-db",
            str(tmp_path),
            "--blast-db-prefix",
            "mini",
            "--output-dir",
            str(output_dir),
            "--blastdbcmd-bin",
            str(failing_blastdbcmd),
        ],
    )

    assert sorted(path.name for path in output_dir.iterdir()) == [
        "best_hit_placements.bed",
        "query_sequences.fasta",
        "selected_references.fasta",
    ]
    assert (output_dir / "query_sequences.fasta").read_text(encoding="utf-8") == ""
    assert (output_dir / "selected_references.fasta").read_text(encoding="utf-8") == ""
    assert (output_dir / "best_hit_placements.bed").read_text(encoding="utf-8") == ""


@pytest.mark.parametrize(
    ("rows", "query_text", "references", "match"),
    [
        ([qbt_row(qlen="9")], ">q1\nAAAACCCC\n", {"REF_A": "ACGTACGTACGTACGT"}, "qlen"),
        (
            [qbt_row(best_hit_reference_length="17")],
            ">q1\nAAAACCCC\n",
            {"REF_A": "ACGTACGTACGTACGT"},
            "best_hit_reference_length",
        ),
        (
            [qbt_row(best_hit_reference_end_1based="17")],
            ">q1\nAAAACCCC\n",
            {"REF_A": "ACGTACGTACGTACGT"},
            "exceeds best_hit_reference_length",
        ),
    ],
)
def test_cross_artifact_integrity_failures_are_loud(
    tmp_path: Path,
    rows: list[dict[str, str]],
    query_text: str,
    references: dict[str, str],
    match: str,
) -> None:
    with pytest.raises(BestHitSequenceEvidenceError, match=match):
        run_evidence(
            tmp_path,
            qbt_rows=rows,
            query_fastas=[("sample-A", query_text)],
            references=references,
        )


def test_missing_query_record_fails_loudly(tmp_path: Path) -> None:
    qbt = tmp_path / "query_big_table.tsv"
    write_tsv(qbt, [qbt_row(qseqid="missing")])
    query_fasta = tmp_path / "queries.fasta"
    query_fasta.write_text(">other\nAAAAAAAA\n", encoding="utf-8")

    with pytest.raises(BestHitSequenceEvidenceError, match="Missing query sequence"):
        main(
            [
                "--query-big-table",
                str(qbt),
                "--query-fasta",
                "sample-A",
                str(query_fasta),
                "--blast-db",
                str(tmp_path),
                "--blast-db-prefix",
                "mini",
                "--output-dir",
                str(tmp_path / "out"),
            ],
        )
