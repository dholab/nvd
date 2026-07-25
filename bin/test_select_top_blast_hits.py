"""Behavioral tests for retaining top BLAST references."""

import csv
import sys
from pathlib import Path

import polars as pl
import pytest
from select_top_blast_hits import main


def blast_row(
    qseqid: str,
    sseqid: str,
    bitscore: int,
    staxids: str,
    stitle: str | None = None,
    *,
    saccver: str | None = None,
    qstart: int = 1,
    qend: int = 214,
    slen: int = 1_000,
    sstart: int = 101,
    send: int = 314,
    sstrand: str = "plus",
) -> str:
    """Build one headerless BLAST outfmt-6 row."""
    return "\t".join(
        [
            qseqid,
            "216",
            sseqid,
            stitle if stitle is not None else f"Reference {sseqid}",
            "214",
            "99.0",
            "1e-50",
            str(bitscore),
            f"Taxon {staxids}",
            staxids,
            saccver if saccver is not None else sseqid,
            str(qstart),
            str(qend),
            str(slen),
            str(sstart),
            str(send),
            sstrand,
        ],
    )


def run_selection(input_file: Path, output_file: Path, retention_count: int) -> None:
    """Run the retention CLI with a controlled argument list."""
    original_argv = sys.argv
    try:
        sys.argv = [
            "select_top_blast_hits.py",
            "--input-file",
            str(input_file),
            "--output-file",
            str(output_file),
            "--blast-retention-count",
            str(retention_count),
        ]
        main()
    finally:
        sys.argv = original_argv


def test_retention_preserves_boundary_ties_before_taxid_expansion(
    tmp_path: Path,
) -> None:
    """A multi-taxid reference and its boundary tie both survive top-five retention."""
    input_file = tmp_path / "blast.tsv"
    output_file = tmp_path / "retained.tsv"
    input_file.write_text(
        "\n".join(
            [
                blast_row("query-1", "ref-1", 100, "1"),
                blast_row("query-1", "ref-2", 99, "2"),
                blast_row("query-1", "ref-3", 98, "3"),
                blast_row("query-1", "ref-4", 97, "4"),
                blast_row("query-1", "ref-5a", 96, "5;6"),
                blast_row("query-1", "ref-5b", 96, "7"),
                blast_row("query-1", "ref-7", 95, "8"),
            ],
        )
        + "\n",
    )

    run_selection(input_file, output_file, 5)

    retained = pl.read_csv(output_file, separator="\t")
    assert set(retained["sseqid"]) == {
        "ref-1",
        "ref-2",
        "ref-3",
        "ref-4",
        "ref-5a",
        "ref-5b",
    }
    assert retained.height == 7
    assert set(retained.filter(pl.col("sseqid") == "ref-5a")["staxids"]) == {5, 6}


def test_retention_count_changes_the_reference_scope(tmp_path: Path) -> None:
    """The user-selected ordinal boundary changes retained references."""
    input_file = tmp_path / "blast.tsv"
    top_one = tmp_path / "top-one.tsv"
    top_two = tmp_path / "top-two.tsv"
    input_file.write_text(
        "\n".join(
            [
                blast_row("query-1", "ref-1", 100, "1"),
                blast_row("query-1", "ref-2", 99, "2"),
                blast_row("query-1", "ref-3", 98, "3"),
            ],
        )
        + "\n",
    )

    run_selection(input_file, top_one, 1)
    run_selection(input_file, top_two, 2)

    assert set(pl.read_csv(top_one, separator="\t")["sseqid"]) == {"ref-1"}
    assert set(pl.read_csv(top_two, separator="\t")["sseqid"]) == {
        "ref-1",
        "ref-2",
    }


def test_repeated_hsps_collapse_to_the_best_reference_row(tmp_path: Path) -> None:
    """Repeated HSPs consume one position and keep deterministic best evidence."""
    input_file = tmp_path / "blast.tsv"
    output_file = tmp_path / "retained.tsv"
    rows = [
        blast_row(
            "query-1",
            "ref-1",
            100,
            "1",
            saccver="NC_000001.1",
            qstart=3,
            qend=210,
            slen=2_000,
            sstart=900,
            send=693,
            sstrand="minus",
        ),
        blast_row(
            "query-1",
            "ref-1",
            90,
            "1",
            saccver="NC_000001.1",
            qstart=10,
            qend=100,
            slen=2_000,
            sstart=400,
            send=491,
        ),
        blast_row("query-1", "ref-2", 99, "2"),
        blast_row("query-1", "ref-3", 98, "3"),
    ]
    input_file.write_text("\n".join(rows) + "\n")

    run_selection(input_file, output_file, 2)

    retained = pl.read_csv(output_file, separator="\t")
    assert set(retained["sseqid"]) == {"ref-1", "ref-2"}
    assert retained.filter(pl.col("sseqid") == "ref-1")["bitscore"].to_list() == [
        100.0,
    ]
    best_reference = retained.filter(pl.col("sseqid") == "ref-1").row(
        0,
        named=True,
    )
    assert best_reference["saccver"] == "NC_000001.1"
    assert best_reference["qstart"] == 3
    assert best_reference["qend"] == 210
    assert best_reference["slen"] == 2_000
    assert best_reference["sstart"] == 900
    assert best_reference["send"] == 693
    assert best_reference["sstrand"] == "minus"


def test_fully_tied_hsps_use_coordinates_as_a_stable_tiebreaker(
    tmp_path: Path,
) -> None:
    input_file = tmp_path / "blast.tsv"
    output_file = tmp_path / "retained.tsv"
    input_file.write_text(
        "\n".join(
            [
                blast_row(
                    "query-1",
                    "ref-1",
                    100,
                    "1",
                    qstart=20,
                    qend=119,
                    sstart=500,
                    send=599,
                ),
                blast_row(
                    "query-1",
                    "ref-1",
                    100,
                    "1",
                    qstart=10,
                    qend=109,
                    sstart=700,
                    send=799,
                ),
            ],
        )
        + "\n",
    )

    run_selection(input_file, output_file, 1)

    [retained] = pl.read_csv(output_file, separator="\t").iter_rows(named=True)
    assert retained["qstart"] == 10
    assert retained["qend"] == 109
    assert retained["sstart"] == 700
    assert retained["send"] == 799


def test_literal_quotes_in_blast_titles_survive_selection(tmp_path: Path) -> None:
    input_file = tmp_path / "blast.tsv"
    output_file = tmp_path / "retained.tsv"
    titles = [
        'Leptotrichia buccalis strain RMA 2181" 16S ribosomal RNA gene',
        '"Literal quote at the start',
    ]
    input_file.write_text(
        "\n".join(
            blast_row(f"query-{index}", f"ref-{index}", 100, "9606", title)
            for index, title in enumerate(titles, start=1)
        )
        + "\n",
    )

    run_selection(input_file, output_file, 1)

    with output_file.open(newline="", encoding="utf-8") as handle:
        retained = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["stitle"] for row in retained] == titles
    assert all(None not in row and None not in row.values() for row in retained)


def test_malformed_numeric_blast_field_fails_explicitly(tmp_path: Path) -> None:
    input_file = tmp_path / "blast.tsv"
    output_file = tmp_path / "retained.tsv"
    malformed = blast_row("query-1", "ref-1", 100, "1").replace("\t100\t", "\tbad\t")
    input_file.write_text(malformed + "\n")

    with pytest.raises((pl.exceptions.ComputeError, pl.exceptions.SchemaError)):
        run_selection(input_file, output_file, 1)


def test_ragged_blast_row_fails_explicitly(tmp_path: Path) -> None:
    input_file = tmp_path / "blast.tsv"
    output_file = tmp_path / "retained.tsv"
    input_file.write_text(blast_row("query-1", "ref-1", 100, "1") + "\textra\n")

    with pytest.raises(ValueError, match="18 fields"):
        run_selection(input_file, output_file, 1)


def test_short_blast_row_fails_explicitly(tmp_path: Path) -> None:
    input_file = tmp_path / "blast.tsv"
    output_file = tmp_path / "retained.tsv"
    short_row = blast_row("query-1", "ref-1", 100, "1").rsplit("\t", 1)[0]
    input_file.write_text(short_row + "\n")

    with pytest.raises(ValueError, match="16 fields"):
        run_selection(input_file, output_file, 1)


def test_missing_taxid_lexeme_becomes_null_for_taxonomy_normalization(
    tmp_path: Path,
) -> None:
    input_file = tmp_path / "blast.tsv"
    output_file = tmp_path / "retained.tsv"
    input_file.write_text(blast_row("query-1", "ref-1", 100, "N/A") + "\n")

    run_selection(input_file, output_file, 1)

    retained = pl.read_csv(output_file, separator="\t")
    assert retained["staxids"].null_count() == 1
