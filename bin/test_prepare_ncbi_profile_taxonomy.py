"""Tests for the NCBI CRUMBS profile taxonomy adapter."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

import pytest
from prepare_ncbi_profile_taxonomy import main
from py_nvd import taxonomy

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def minimal_taxdump(tmp_path: Path) -> Path:
    taxdump_dir = tmp_path / "taxdump"
    taxdump_dir.mkdir()

    (taxdump_dir / "nodes.dmp").write_text(
        """\
1\t|\t1\t|\tno rank\t|\t
131567\t|\t1\t|\tno rank\t|\t
2759\t|\t131567\t|\tsuperkingdom\t|\t
33208\t|\t2759\t|\tkingdom\t|\t
9604\t|\t33208\t|\tfamily\t|\t
207598\t|\t9604\t|\tsubfamily\t|\t
9605\t|\t207598\t|\tgenus\t|\t
9606\t|\t9605\t|\tspecies\t|\t
9596\t|\t207598\t|\tgenus\t|\t
9598\t|\t9596\t|\tspecies\t|\t
""",
        encoding="utf-8",
    )
    (taxdump_dir / "names.dmp").write_text(
        """\
1\t|\troot\t|\t\t|\tscientific name\t|
131567\t|\tcellular organisms\t|\t\t|\tscientific name\t|
2759\t|\tEukaryota\t|\t\t|\tscientific name\t|
33208\t|\tMetazoa\t|\t\t|\tscientific name\t|
9604\t|\tHominidae\t|\t\t|\tscientific name\t|
207598\t|\tHomininae\t|\t\t|\tscientific name\t|
9605\t|\tHomo\t|\t\t|\tscientific name\t|
9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|
9596\t|\tPan\t|\t\t|\tscientific name\t|
9598\t|\tPan troglodytes\t|\t\t|\tscientific name\t|
""",
        encoding="utf-8",
    )
    (taxdump_dir / "merged.dmp").write_text(
        """\
12345\t|\t9606\t|
""",
        encoding="utf-8",
    )

    taxonomy._build_sqlite_from_dmp(taxdump_dir)  # noqa: SLF001 - test fixture builder
    return taxdump_dir


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_observed_taxids_are_written_as_profile_taxonomy_rows(
    tmp_path: Path,
    minimal_taxdump: Path,
) -> None:
    taxids = tmp_path / "observed_taxids.txt"
    taxids.write_text("9606\n9598\n", encoding="utf-8")
    output = tmp_path / "profile_taxonomy.tsv"

    main(
        [
            "--taxids",
            str(taxids),
            "--taxonomy-dir",
            str(minimal_taxdump),
            "--output",
            str(output),
        ],
    )

    assert read_rows(output) == [
        {
            "taxon_id": "9606",
            "taxon_name": "Homo sapiens",
            "rank": "species",
            "taxpath": "131567|2759|33208|9604|207598|9605|9606",
            "taxpathsn": "cellular organisms|Eukaryota|Metazoa|Hominidae|Homininae|Homo|Homo sapiens",
        },
        {
            "taxon_id": "9598",
            "taxon_name": "Pan troglodytes",
            "rank": "species",
            "taxpath": "131567|2759|33208|9604|207598|9596|9598",
            "taxpathsn": "cellular organisms|Eukaryota|Metazoa|Hominidae|Homininae|Pan|Pan troglodytes",
        },
    ]


def test_merged_taxids_collapse_to_current_taxon_id(
    tmp_path: Path,
    minimal_taxdump: Path,
) -> None:
    taxids = tmp_path / "observed_taxids.txt"
    taxids.write_text("12345\n9606\n", encoding="utf-8")
    output = tmp_path / "profile_taxonomy.tsv"

    main(
        [
            "--taxids",
            str(taxids),
            "--taxonomy-dir",
            str(minimal_taxdump),
            "--output",
            str(output),
        ],
    )

    rows = read_rows(output)
    assert len(rows) == 1
    assert rows[0]["taxon_id"] == "9606"
    assert rows[0]["taxon_name"] == "Homo sapiens"


def test_unknown_taxid_fails_loudly(tmp_path: Path, minimal_taxdump: Path) -> None:
    taxids = tmp_path / "observed_taxids.txt"
    taxids.write_text("999999\n", encoding="utf-8")
    output = tmp_path / "profile_taxonomy.tsv"

    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "--taxids",
                str(taxids),
                "--taxonomy-dir",
                str(minimal_taxdump),
                "--output",
                str(output),
            ],
        )

    assert exc_info.value.code == 1
