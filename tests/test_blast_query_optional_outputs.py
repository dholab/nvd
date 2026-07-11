"""Nextflow contract tests for empty BLAST query batches."""

from __future__ import annotations

import gzip
import os
import shutil
import sqlite3
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BLAST_MODULE = ROOT / "modules" / "blast"


def write_empty_query_lookup(path: Path) -> None:
    """Write the lookup schema consumed by SELECT_BLAST_QUERIES."""
    with sqlite3.connect(path) as connection:
        connection.execute(
            """
            create table query_sequences (
                qseqid text primary key,
                query_class text not null
            )
            """,
        )


def test_empty_query_producers_succeed_without_emitting_tuples(tmp_path: Path) -> None:
    """Empty contig and read inputs should close their output channels cleanly."""
    nextflow = shutil.which("nextflow")
    assert nextflow is not None

    empty_fasta = tmp_path / "water.target_filtered.fasta"
    empty_fasta.touch()
    lookup = tmp_path / "water.query_sequences.sqlite"
    write_empty_query_lookup(lookup)
    empty_fastq = tmp_path / "water.single_read.fastq.gz"
    with gzip.open(empty_fastq, "wb"):
        pass

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ SELECT_BLAST_QUERIES; NORMALIZE_READ_BLAST_QUERIES }} from '{BLAST_MODULE}'

workflow {{
    SELECT_BLAST_QUERIES(Channel.of(tuple(
        'water',
        'illumina',
        'single',
        file('{empty_fasta}'),
        file('{lookup}'),
    )))

    NORMALIZE_READ_BLAST_QUERIES(Channel.of(tuple(
        'water',
        'illumina',
        'single_read',
        file('{empty_fastq}'),
    )))

    SELECT_BLAST_QUERIES.out.query_fastas
        .mix(NORMALIZE_READ_BLAST_QUERIES.out.queries)
        .view {{ value -> throw new IllegalStateException("unexpected query tuple: ${{value}}") }}
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["PATH"] = f"{ROOT / 'bin'}{os.pathsep}{environment['PATH']}"
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
    assert "Missing output file(s)" not in diagnostics
