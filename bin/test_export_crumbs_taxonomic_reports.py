"""Tests for CRUMBS taxonomic report exports."""

from __future__ import annotations

import csv
import shutil
import subprocess
from typing import TYPE_CHECKING

import pytest
from export_crumbs_taxonomic_reports import CrumbsReportExportError, main

if TYPE_CHECKING:
    from pathlib import Path


TAXA_FIELDS = [
    "sample_id",
    "taxon_id",
    "taxon_name",
    "rank",
    "taxpath",
    "taxpathsn",
    "rankpath",
    "n_contigs",
    "total_contig_length",
    "total_covered_bases_1x",
    "total_raw_aligned_bases",
    "total_crumbs_score",
    "taxon_crumbs",
    "percentage_emitted",
    "n_zero_crumbs_contigs",
    "median_contig_breadth_1x",
    "min_contig_breadth_1x",
    "max_contig_breadth_1x",
    "fraction_crumbs_from_top_contig",
    "fraction_crumbs_from_low_breadth_contigs",
]


def write_tsv(path: Path, rows: list[dict[str, object]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TAXA_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        return next(reader)


def taxon_row(  # noqa: PLR0913 - explicit taxonomy fixtures are clearer here.
    *,
    taxon_id: int,
    taxon_name: str,
    rank: str,
    taxpath: str,
    taxpathsn: str,
    rankpath: str,
    percentage_emitted: float,
) -> dict[str, object]:
    return {
        "sample_id": "sample-1",
        "taxon_id": taxon_id,
        "taxon_name": taxon_name,
        "rank": rank,
        "taxpath": taxpath,
        "taxpathsn": taxpathsn,
        "rankpath": rankpath,
        "n_contigs": 1,
        "total_contig_length": 100,
        "total_covered_bases_1x": 100,
        "total_raw_aligned_bases": 100,
        "total_crumbs_score": 100,
        "taxon_crumbs": 1.0,
        "percentage_emitted": percentage_emitted,
        "n_zero_crumbs_contigs": 0,
        "median_contig_breadth_1x": 1.0,
        "min_contig_breadth_1x": 1.0,
        "max_contig_breadth_1x": 1.0,
        "fraction_crumbs_from_top_contig": 1.0,
        "fraction_crumbs_from_low_breadth_contigs": 0.0,
    }


def run_export(tmp_path: Path, taxa_tsv: Path, *, scale: int = 1_000_000_000) -> None:
    main(
        [
            "--sample-id",
            "sample-1",
            "--taxa-tsv",
            str(taxa_tsv),
            "--output-dir",
            str(tmp_path),
            "--profile-mass-scale",
            str(scale),
        ],
    )


def test_exports_krona_rows_with_taxburst_rank_columns(tmp_path: Path) -> None:
    taxa_tsv = write_tsv(
        tmp_path / "sample-1.crumbs.taxa.tsv",
        [
            taxon_row(
                taxon_id=9606,
                taxon_name="Homo sapiens",
                rank="species",
                taxpath="2759|7711|40674|9443|9604|9605|9606",
                taxpathsn="Eukaryota|Chordata|Mammalia|Primates|Hominidae|Homo|Homo sapiens",
                rankpath="superkingdom|phylum|class|order|family|genus|species",
                percentage_emitted=33.3333333333,
            ),
            taxon_row(
                taxon_id=10244,
                taxon_name="Monkeypox virus",
                rank="species",
                taxpath="10239|2732408|2732506|2732544|10240|10242|10244",
                taxpathsn="Viruses|Monodnaviria|Shotokuvirae|Cossaviricota|Quintoviricetes|Piccovirales|Monkeypox virus",
                rankpath="superkingdom|realm|kingdom|phylum|class|order|species",
                percentage_emitted=66.6666666667,
            ),
        ],
    )

    run_export(tmp_path, taxa_tsv)

    krona_tsv = tmp_path / "sample-1.crumbs.krona.tsv"
    rows = read_tsv(krona_tsv)
    assert read_tsv_header(krona_tsv) == [
        "fraction",
        "superkingdom",
        "phylum",
        "class",
        "order",
    ]
    assert rows[0] == {
        "fraction": "0.333333333333",
        "superkingdom": "Eukaryota",
        "phylum": "Chordata",
        "class": "Mammalia",
        "order": "Primates",
    }
    assert rows[1]["superkingdom"] == "Viruses"
    assert rows[1]["phylum"] == "Cossaviricota"
    assert rows[1]["class"] == "Quintoviricetes"
    assert rows[1]["order"] == "Piccovirales"


def test_krona_header_omits_trailing_blank_ranks_for_taxburst(
    tmp_path: Path,
) -> None:
    taxa_tsv = write_tsv(
        tmp_path / "sample-1.crumbs.taxa.tsv",
        [
            taxon_row(
                taxon_id=10258,
                taxon_name="Orf virus",
                rank="species",
                taxpath="10239|2732408|2732506|2732544|10240|10255|10258",
                taxpathsn="Viruses|Nucleocytoviricota|Pokkesviricetes|Chitovirales|Poxviridae|Parapoxvirus|Orf virus",
                rankpath="superkingdom|phylum|class|order|family|genus|species",
                percentage_emitted=100.0,
            ),
        ],
    )

    run_export(tmp_path, taxa_tsv)

    krona_tsv = tmp_path / "sample-1.crumbs.krona.tsv"
    assert read_tsv_header(krona_tsv) == [
        "fraction",
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    assert (
        krona_tsv.read_text(encoding="utf-8")
        .splitlines()[1]
        .endswith(
            "\tOrf virus",
        )
    )


def test_krona_assigns_rankless_taxa_to_unclassified_superkingdom(
    tmp_path: Path,
) -> None:
    taxa_tsv = write_tsv(
        tmp_path / "sample-1.crumbs.taxa.tsv",
        [
            taxon_row(
                taxon_id=1,
                taxon_name="root",
                rank="no rank",
                taxpath="1",
                taxpathsn="root",
                rankpath="no rank",
                percentage_emitted=100.0,
            ),
        ],
    )

    run_export(tmp_path, taxa_tsv)

    krona_tsv = tmp_path / "sample-1.crumbs.krona.tsv"
    assert read_tsv_header(krona_tsv) == ["fraction", "superkingdom"]
    assert read_tsv(krona_tsv) == [
        {"fraction": "1", "superkingdom": "unclassified"},
    ]

    taxburst = shutil.which("taxburst")
    assert taxburst is not None
    subprocess.run(  # noqa: S603 - exercise the installed TaxBurst CLI compatibility seam.
        [
            taxburst,
            "-F",
            "krona",
            str(krona_tsv),
            "-o",
            str(tmp_path / "sample-1.crumbs.taxburst.html"),
            "--save-json",
            str(tmp_path / "sample-1.crumbs.taxburst.json"),
        ],
        check=True,
    )


def test_exports_ancestor_expanded_kreport_with_scaled_profile_mass(
    tmp_path: Path,
) -> None:
    taxa_tsv = write_tsv(
        tmp_path / "sample-1.crumbs.taxa.tsv",
        [
            taxon_row(
                taxon_id=9605,
                taxon_name="Homo",
                rank="genus",
                taxpath="131567|2759|9605",
                taxpathsn="cellular organisms|Eukaryota|Homo",
                rankpath="no rank|superkingdom|genus",
                percentage_emitted=10.0,
            ),
            taxon_row(
                taxon_id=9606,
                taxon_name="Homo sapiens",
                rank="species",
                taxpath="131567|2759|9605|9606",
                taxpathsn="cellular organisms|Eukaryota|Homo|Homo sapiens",
                rankpath="no rank|superkingdom|genus|species",
                percentage_emitted=30.0,
            ),
        ],
    )

    run_export(tmp_path, taxa_tsv, scale=1_000)

    lines = (
        (tmp_path / "sample-1.crumbs.kreport.txt")
        .read_text(
            encoding="utf-8",
        )
        .splitlines()
    )
    assert lines == [
        "40.00\t400\t0\t-\t131567\tcellular organisms",
        "40.00\t400\t0\tD\t2759\t  Eukaryota",
        "40.00\t400\t100\tG\t9605\t    Homo",
        "30.00\t300\t300\tS\t9606\t      Homo sapiens",
    ]


def test_all_zero_profile_writes_empty_report_files(tmp_path: Path) -> None:
    taxa_tsv = write_tsv(
        tmp_path / "sample-1.crumbs.taxa.tsv",
        [
            taxon_row(
                taxon_id=9606,
                taxon_name="Homo sapiens",
                rank="species",
                taxpath="131567|2759|9605|9606",
                taxpathsn="cellular organisms|Eukaryota|Homo|Homo sapiens",
                rankpath="no rank|superkingdom|genus|species",
                percentage_emitted=0.0,
            ),
        ],
    )

    run_export(tmp_path, taxa_tsv)

    assert read_tsv(tmp_path / "sample-1.crumbs.krona.tsv") == []
    assert (tmp_path / "sample-1.crumbs.kreport.txt").read_text(encoding="utf-8") == ""


def test_misaligned_taxonomy_paths_fail_loudly(tmp_path: Path) -> None:
    taxa_tsv = write_tsv(
        tmp_path / "sample-1.crumbs.taxa.tsv",
        [
            taxon_row(
                taxon_id=9606,
                taxon_name="Homo sapiens",
                rank="species",
                taxpath="131567|2759|9605|9606",
                taxpathsn="cellular organisms|Eukaryota|Homo|Homo sapiens",
                rankpath="no rank|superkingdom|species",
                percentage_emitted=100.0,
            ),
        ],
    )

    with pytest.raises(
        CrumbsReportExportError,
        match="misaligned taxpath, taxpathsn, and rankpath",
    ):
        run_export(tmp_path, taxa_tsv)
