import sys
from pathlib import Path

from concat_multiblast_hits import main

ANNOTATED_BLAST_COLUMNS = [
    "task",
    "sample",
    "qseqid",
    "qlen",
    "sseqid",
    "stitle",
    "length",
    "pident",
    "evalue",
    "bitscore",
    "sscinames",
    "staxids",
    "saccver",
    "qstart",
    "qend",
    "slen",
    "sstart",
    "send",
    "sstrand",
    "rank",
]


def test_header_only_inputs_preserve_lca_columns(
    tmp_path: Path,
    monkeypatch,
) -> None:
    header = "\t".join(
        [
            *ANNOTATED_BLAST_COLUMNS,
            "adjusted_taxid",
            "adjusted_taxid_name",
            "adjusted_taxid_rank",
            "adjustment_method",
        ],
    )
    batch_a = tmp_path / "sample.short_assembly_contig.blast.merged_with_lca.tsv"
    batch_b = tmp_path / "sample.long_assembly_contig.blast.merged_with_lca.tsv"
    output = tmp_path / "sample_blast.merged_with_lca.tsv"
    batch_a.write_text(f"{header}\n", encoding="utf-8")
    batch_b.write_text(f"{header}\n", encoding="utf-8")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "concat_multiblast_hits.py",
            "--blast-hits",
            str(batch_a),
            "--blast-hits",
            str(batch_b),
            "--output-file",
            str(output),
        ],
    )

    main()

    assert output.read_text(encoding="utf-8") == f"{header}\n"


def test_empty_inputs_use_the_expanded_annotated_blast_schema(
    tmp_path: Path,
    monkeypatch,
) -> None:
    batch_a = tmp_path / "sample.megablast.tsv"
    batch_b = tmp_path / "sample.blastn.tsv"
    output = tmp_path / "sample.combined.tsv"
    batch_a.touch()
    batch_b.touch()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "concat_multiblast_hits.py",
            "--blast-hits",
            str(batch_a),
            "--blast-hits",
            str(batch_b),
            "--output-file",
            str(output),
        ],
    )

    main()

    assert output.read_text(encoding="utf-8") == (
        "\t".join(ANNOTATED_BLAST_COLUMNS) + "\n"
    )
