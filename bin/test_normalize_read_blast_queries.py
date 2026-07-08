from __future__ import annotations

import hashlib
import sqlite3
from pathlib import Path

import pytest
from normalize_read_blast_queries import (
    ReadQueryNormalizationError,
    VsearchUnique,
    main,
    normalize_unique_reads,
    parse_vsearch_fasta,
    write_fasta,
    write_query_lookup,
)


def lookup_rows(path: Path) -> list[dict[str, object]]:
    with sqlite3.connect(path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(
            """
            select
                qseqid,
                sample_id,
                evidence_class,
                producer,
                source_id,
                support_record_count,
                length,
                sha256
            from query_sequences
            order by qseqid
            """,
        ).fetchall()
    return [dict(row) for row in rows]


@pytest.mark.parametrize(
    (
        "support_count_policy",
        "expected_support_counts",
        "expected_collapsed_support",
    ),
    [
        ("abundance", [1, 1, 1, 1, 3], 3),
        ("one", [1, 1, 1, 1, 1], 1),
    ],
)
def test_normalizes_distinct_casava_mates_without_losing_the_evidence_class(
    tmp_path: Path,
    support_count_policy: str,
    expected_support_counts: list[int],
    expected_collapsed_support: int,
) -> None:
    fastq = tmp_path / "water.fastq"
    fastq.write_text(
        """@LH00642:196:23L3KYLT3:2:2282:37185:27261 1:N:0:CAATCGGCTG+TTCCTACAGC
ACGT
+
IIII
@LH00642:196:23L3KYLT3:2:2282:37185:27261 2:N:0:CAATCGGCTG+TTCCTACAGC
TGCA
+
IIII
@cluster-b 1:N:0:INDEX
CCCC
+
IIII
@cluster-c 1:N:0:INDEX
CCCC
+
IIII
@cluster-d 2:N:0:INDEX
CCCC
+
IIII
@cluster-e 1:N:0:INDEX
GGGG
+
IIII
@cluster-f 2:N:0:INDEX
TTTT
+
IIII
""",
        encoding="utf-8",
    )

    main(
        [
            "--sample-id",
            "Water_20260707",
            "--evidence-class",
            "single_read",
            "--producer",
            "source_read",
            "--fastq",
            str(fastq),
            "--support-count-policy",
            support_count_policy,
            "--output-dir",
            str(tmp_path),
        ],
    )

    fasta = tmp_path / "Water_20260707.single_read.blast_queries.fasta"
    lookup = tmp_path / "Water_20260707.single_read.query_sequences.sqlite"
    rows = lookup_rows(lookup)
    fasta_headers = [
        line[1:].split()
        for line in fasta.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]
    fasta_qseqids = [fields[0] for fields in fasta_headers]
    fasta_source_ids = {
        fields[3].removeprefix("source_id=") for fields in fasta_headers
    }
    expected_query_count = 5

    assert len(rows) == expected_query_count
    assert {len(fields) for fields in fasta_headers} == {4}
    assert len(set(fasta_qseqids)) == expected_query_count
    assert fasta_qseqids == [row["qseqid"] for row in rows]
    sqlite_source_ids = {row["source_id"] for row in rows}
    assert fasta_source_ids == sqlite_source_ids
    assert sqlite_source_ids >= {
        "LH00642:196:23L3KYLT3:2:2282:37185:27261/1",
        "LH00642:196:23L3KYLT3:2:2282:37185:27261/2",
    }
    assert (
        sorted(row["support_record_count"] for row in rows) == expected_support_counts
    )
    collapsed_row = next(row for row in rows if row["source_id"] == "cluster-b")
    collapsed_sequence = b"CCCC"
    assert collapsed_row["support_record_count"] == expected_collapsed_support
    assert collapsed_row["length"] == len(collapsed_sequence)
    assert collapsed_row["sha256"] == hashlib.sha256(collapsed_sequence).hexdigest()


def test_preserves_existing_mate_suffixes_in_casava_labels(tmp_path: Path) -> None:
    vsearch_fasta = tmp_path / "uniques.fasta"
    vsearch_fasta.write_text(
        ">LH00642:196:23L3KYLT3:2:2282:37185:27261/1 1:N:0:INDEX;size=1\nACGT\n"
        ">LH00642:196:23L3KYLT3:2:2282:37185:27261/2 2:N:0:INDEX;size=1\nTGCA\n",
        encoding="utf-8",
    )

    with vsearch_fasta.open(encoding="utf-8") as handle:
        uniques = parse_vsearch_fasta(handle)

    assert [unique.source_id for unique in uniques] == [
        "LH00642:196:23L3KYLT3:2:2282:37185:27261/1",
        "LH00642:196:23L3KYLT3:2:2282:37185:27261/2",
    ]


def test_discards_non_casava_comments_from_source_ids(tmp_path: Path) -> None:
    vsearch_fasta = tmp_path / "uniques.fasta"
    vsearch_fasta.write_text(
        ">custom-read instrument metadata;size=1\nACGT\n",
        encoding="utf-8",
    )

    with vsearch_fasta.open(encoding="utf-8") as handle:
        uniques = parse_vsearch_fasta(handle)

    assert uniques[0].source_id == "custom-read"


def test_does_not_treat_a_casava_shaped_custom_comment_as_a_mate(
    tmp_path: Path,
) -> None:
    vsearch_fasta = tmp_path / "uniques.fasta"
    vsearch_fasta.write_text(
        ">custom-read 1:N:0:INDEX;size=1\nACGT\n",
        encoding="utf-8",
    )

    with vsearch_fasta.open(encoding="utf-8") as handle:
        uniques = parse_vsearch_fasta(handle)

    assert uniques[0].source_id == "custom-read"


def test_rejects_conflicting_suffix_and_casava_mate(tmp_path: Path) -> None:
    vsearch_fasta = tmp_path / "uniques.fasta"
    vsearch_fasta.write_text(
        ">LH00642:196:23L3KYLT3:2:2282:37185:27261/1 2:N:0:INDEX;size=1\nACGT\n",
        encoding="utf-8",
    )

    with (
        vsearch_fasta.open(encoding="utf-8") as handle,
        pytest.raises(
            ReadQueryNormalizationError,
            match="conflicting mate annotations",
        ),
    ):
        parse_vsearch_fasta(handle)


def test_rejects_reused_full_identifier_for_distinct_sequences(tmp_path: Path) -> None:
    fastq = tmp_path / "reused.fastq"
    fastq.write_text(
        """@LH00642:196:23L3KYLT3:2:2282:37185:27261 1:N:0:INDEX
ACGT
+
IIII
@LH00642:196:23L3KYLT3:2:2282:37185:27261 1:N:0:INDEX
TGCA
+
IIII
""",
        encoding="utf-8",
    )

    with pytest.raises(
        ReadQueryNormalizationError,
        match=(
            "duplicate vsearch representative source_id: "
            "'LH00642:196:23L3KYLT3:2:2282:37185:27261/1'"
        ),
    ):
        main(
            [
                "--sample-id",
                "sample-1",
                "--evidence-class",
                "single_read",
                "--producer",
                "source_read",
                "--fastq",
                str(fastq),
                "--support-count-policy",
                "abundance",
                "--output-dir",
                str(tmp_path),
            ],
        )


def test_rejects_unsafe_sample_id_before_creating_output_dir(tmp_path: Path) -> None:
    fastq = tmp_path / "empty.fastq"
    fastq.touch()
    output_dir = tmp_path / "output"

    with pytest.raises(
        ReadQueryNormalizationError,
        match="sample_id must contain only letters",
    ):
        main(
            [
                "--sample-id",
                "../escaped",
                "--evidence-class",
                "single_read",
                "--producer",
                "source_read",
                "--fastq",
                str(fastq),
                "--support-count-policy",
                "abundance",
                "--output-dir",
                str(output_dir),
                "--vsearch-bin",
                "definitely-not-vsearch",
            ],
        )

    assert not output_dir.exists()


def test_rejects_repeated_representative_source_ids() -> None:
    with pytest.raises(
        ReadQueryNormalizationError,
        match="duplicate vsearch representative source_id: 'readA/1'",
    ):
        normalize_unique_reads(
            [
                VsearchUnique(source_id="readA/1", abundance=1, sequence="ACGT"),
                VsearchUnique(source_id="readA/1", abundance=1, sequence="TGCA"),
            ],
            sample_id="sample-1",
            evidence_class="single_read",
            producer="source_read",
            support_count_policy="abundance",
        )


def test_parses_vsearch_fasta_size_annotations(tmp_path: Path) -> None:
    vsearch_fasta = tmp_path / "uniques.fasta"
    vsearch_fasta.write_text(
        ">readA/1;size=3\nACGTACGT\n>readC/1;size=1\nTTTTAAAA\n",
        encoding="utf-8",
    )

    with vsearch_fasta.open(encoding="utf-8") as handle:
        uniques = parse_vsearch_fasta(handle)

    assert uniques == [
        VsearchUnique(source_id="readA/1", abundance=3, sequence="ACGTACGT"),
        VsearchUnique(source_id="readC/1", abundance=1, sequence="TTTTAAAA"),
    ]


def test_writes_single_read_queries_with_abundance_support(tmp_path: Path) -> None:
    queries = normalize_unique_reads(
        [
            VsearchUnique(source_id="readA/1", abundance=3, sequence="ACGTACGT"),
            VsearchUnique(source_id="readC/1", abundance=1, sequence="TTTTAAAA"),
        ],
        sample_id="sample-1",
        evidence_class="single_read",
        producer="source_read",
        support_count_policy="abundance",
    )
    fasta = tmp_path / "sample-1.single_read.blast_queries.fasta"
    lookup = tmp_path / "sample-1.single_read.query_sequences.sqlite"

    write_fasta(fasta, queries)
    write_query_lookup(lookup, queries)

    assert fasta.read_text(encoding="utf-8") == (
        ">nvdReadQuery_sample-1_000001 evidence_class=single_read producer=source_read source_id=readA/1\n"
        "ACGTACGT\n"
        ">nvdReadQuery_sample-1_000002 evidence_class=single_read producer=source_read source_id=readC/1\n"
        "TTTTAAAA\n"
    )
    assert lookup_rows(lookup) == [
        {
            "qseqid": "nvdReadQuery_sample-1_000001",
            "sample_id": "sample-1",
            "evidence_class": "single_read",
            "producer": "source_read",
            "source_id": "readA/1",
            "support_record_count": 3,
            "length": 8,
            "sha256": hashlib.sha256(b"ACGTACGT").hexdigest(),
        },
        {
            "qseqid": "nvdReadQuery_sample-1_000002",
            "sample_id": "sample-1",
            "evidence_class": "single_read",
            "producer": "source_read",
            "source_id": "readC/1",
            "support_record_count": 1,
            "length": 8,
            "sha256": hashlib.sha256(b"TTTTAAAA").hexdigest(),
        },
    ]


def test_deduplicated_support_policy_uses_one_per_unique_query() -> None:
    queries = normalize_unique_reads(
        [VsearchUnique(source_id="readA/1", abundance=9, sequence="ACGT")],
        sample_id="sample-1",
        evidence_class="overlap_merged_pair",
        producer="bbmerge",
        support_count_policy="one",
    )

    assert queries[0].qseqid == "nvdMergeReadQuery_sample-1_000001"
    assert queries[0].support_record_count == 1


def test_rejects_vsearch_fasta_without_size_annotation(tmp_path: Path) -> None:
    vsearch_fasta = tmp_path / "bad.fasta"
    vsearch_fasta.write_text(">readA/1\nACGT\n", encoding="utf-8")

    with (
        vsearch_fasta.open(encoding="utf-8") as handle,
        pytest.raises(ReadQueryNormalizationError, match="lacks ';size='"),
    ):
        parse_vsearch_fasta(handle)
