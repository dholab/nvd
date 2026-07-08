import sys
from pathlib import Path

from concat_multiblast_hits import main


def test_header_only_inputs_preserve_lca_columns(
    tmp_path: Path,
    monkeypatch,
) -> None:
    header = "\t".join(
        [
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
            "rank",
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
