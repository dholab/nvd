"""Runtime contracts for BLAST samples without contig-derived products."""

from __future__ import annotations

import csv
import os
import shutil
import sqlite3
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UTILS_MODULE = ROOT / "modules" / "utils"
CRUMBS_MODULE = ROOT / "modules" / "crumbs"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_read_query_lookup(path: Path) -> Path:
    with sqlite3.connect(path) as connection:
        connection.execute(
            """
            create table query_sequences (
                qseqid text primary key,
                sample_id text not null,
                query_class text not null,
                producer text not null,
                source_id text not null,
                support_record_count integer not null,
                length integer not null,
                sha256 text not null
            )
            """,
        )
        connection.execute(
            """
            insert into query_sequences values (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "nvdReadQuery_water_000001",
                "water",
                "single_read",
                "source_read",
                "read-1",
                3,
                25,
                "sha",
            ),
        )
    return path


def test_read_only_sample_accepts_empty_count_and_coverage_lists(
    tmp_path: Path,
) -> None:
    nextflow = shutil.which("nextflow")
    assert nextflow is not None

    blast = write_tsv(
        tmp_path / "water.blast.tsv",
        [
            "sample",
            "qseqid",
            "qlen",
            "adjusted_taxid",
            "adjusted_taxid_name",
            "adjusted_taxid_rank",
            "adjustment_method",
        ],
        [
            {
                "sample": "water",
                "qseqid": "nvdReadQuery_water_000001",
                "qlen": 25,
                "adjusted_taxid": 10239,
                "adjusted_taxid_name": "Viruses",
                "adjusted_taxid_rank": "superkingdom",
                "adjustment_method": "dominant",
            },
        ],
    )
    lookup = write_read_query_lookup(tmp_path / "water.query_sequences.sqlite")
    taxonomy = write_tsv(
        tmp_path / "profile_taxonomy.tsv",
        ["taxon_id", "taxon_name", "rank", "taxpath", "taxpathsn", "rankpath"],
        [
            {
                "taxon_id": 10239,
                "taxon_name": "Viruses",
                "rank": "superkingdom",
                "taxpath": "1|10239",
                "taxpathsn": "root|Viruses",
                "rankpath": "no rank|superkingdom",
            },
        ],
    )

    lib = tmp_path / "lib"
    lib.mkdir()
    shutil.copy2(ROOT / "lib" / "NvdUtils.groovy", lib / "NvdUtils.groovy")

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ ADD_READ_COUNTS_TO_BLAST }} from '{UTILS_MODULE}'
include {{ ESTIMATE_CRUMBS_PROFILE }} from '{CRUMBS_MODULE}'

params.no_enrichment = true

workflow {{
    ADD_READ_COUNTS_TO_BLAST(
        Channel.of(tuple(
            'water',
            file('{blast}'),
            10,
            [],
            [file('{lookup}')],
        )),
        Channel.value('run-test'),
    )

    ESTIMATE_CRUMBS_PROFILE(
        ADD_READ_COUNTS_TO_BLAST.out.map {{ sample_id, final_blast ->
            tuple(sample_id, final_blast, [], file('{taxonomy}'))
        }}
    )

    ESTIMATE_CRUMBS_PROFILE.out.queries.view {{ query -> "READ_ONLY_QUERY: ${{query}}" }}
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
    assert "READ_ONLY_QUERY:" in completed.stdout
