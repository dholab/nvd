"""Behavior tests for WVDB sourmash reference preparation."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "prepare_wvdb_sourmash_inputs.py"
MANIFEST_HEADER = (
    "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,"
    "with_abundance,name,filename"
)
LINEAGE_HEADER = "ident,superkingdom,phylum,class,order,family,genus,species,strain"


def run_script(*args: str | Path) -> subprocess.CompletedProcess[str]:
    """Run the WVDB preparation script through its command-line interface."""
    return subprocess.run(  # noqa: S603
        [sys.executable, str(SCRIPT), *(str(arg) for arg in args)],
        check=False,
        cwd=ROOT,
        text=True,
        capture_output=True,
    )


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows as dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read diagnostic TSV rows as dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv_rows(path: Path, rows: list[list[str]]) -> None:
    """Write TSV rows to a fixture file."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def write_manifest(
    path: Path,
    rows: list[tuple[str, str]],
    *,
    scaled: str = "50",
) -> Path:
    """Write a minimal sourmash manifest from identifier/digest pairs."""
    lines = ["# SOURMASH-MANIFEST-VERSION: 1.0", MANIFEST_HEADER]
    lines.extend(
        f"signatures/{md5}.sig.gz,{md5},{md5[:8]},31,DNA,0,{scaled},10,False,{ident},ref.fa"
        for ident, md5 in rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_lineages(path: Path, rows: list[tuple[str, ...]]) -> Path:
    """Write minimal lineage rows for curation tests."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(LINEAGE_HEADER.split(","))
        writer.writerows(rows)
    return path


def run_curation(
    tmp_path: Path,
    *,
    manifests: tuple[Path, Path],
    lineages: tuple[Path, Path],
    policy: str,
) -> tuple[subprocess.CompletedProcess[str], Path, Path]:
    """Run lineage curation and return the process and output paths."""
    output = tmp_path / "combined.lineages.csv"
    diagnostics = tmp_path / "taxonomy-conflicts.tsv"
    result = run_script(
        "curate-lineages",
        "--reference-manifest-csv",
        manifests[0],
        "--wvdb-manifest-csv",
        manifests[1],
        "--reference-lineages-csv",
        lineages[0],
        "--wvdb-lineages-csv",
        lineages[1],
        "--taxonomy-conflict-policy",
        policy,
        "--lineages-csv",
        output,
        "--diagnostics-tsv",
        diagnostics,
    )
    return result, output, diagnostics


def write_wvdb_fixtures(tmp_path: Path) -> tuple[Path, Path]:
    """Write tiny WVDB-like FASTA and annotation fixtures."""
    fasta = tmp_path / "WVDB_v1.0.fasta"
    annotations = tmp_path / "WVDB_v1.0_annotations.tsv"

    fasta.write_text(
        """\
>contig_with_species original WVDB description
ACGTACGT
>contig_without_species
TGCATGCA
""",
        encoding="utf-8",
    )
    write_tsv_rows(
        annotations,
        [
            [
                "contig",
                "Phylum_gNd",
                "Class_gNd",
                "Order_gNd",
                "Family_gNd",
                "Phylum_RdRp",
                "Class_RdRp",
                "Order_RdRp",
                "Family_RdRp",
                "Genus_RdRp",
                "Species_RdRp",
                "Order_consensus",
                "Family_consensus",
                "Family_ICTV",
            ],
            [
                "contig_with_species",
                "Fallbackphylum",
                "Fallbackclass",
                "Fallbackorder",
                "Fallbackfamily",
                "RdRpphylum",
                "RdRpclass",
                "RdRporder",
                "RdRpfamily",
                "Knownvirus",
                "Knownvirus species",
                "Consensusorder",
                "Consensusfamily",
                "ICTVfamily",
            ],
            [
                "contig_without_species",
                "Fallbackonlyphylum",
                "Fallbackonlyclass",
                "Fallbackonlyorder",
                "Fallbackonlyfamily",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "ICTVfallbackfamily",
            ],
        ],
    )
    return fasta, annotations


def test_prepare_inputs_normalizes_fasta_and_maps_lineages(tmp_path: Path) -> None:
    """prepare-inputs writes normalized FASTA IDs and sourmash lineages."""
    fasta, annotations = write_wvdb_fixtures(tmp_path)
    normalized_fasta = tmp_path / "WVDB_v1.0.normalized.fasta"
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--normalized-fasta",
        normalized_fasta,
        "--lineages-csv",
        lineages_csv,
    )

    assert result.returncode == 0, result.stderr
    assert normalized_fasta.read_text(encoding="utf-8") == (
        ">WVDB|contig_with_species\nACGTACGT\n>WVDB|contig_without_species\nTGCATGCA\n"
    )

    rows = read_csv_rows(lineages_csv)
    assert rows == [
        {
            "ident": "WVDB|contig_with_species",
            "superkingdom": "Viruses",
            "phylum": "RdRpphylum",
            "class": "RdRpclass",
            "order": "Consensusorder",
            "family": "Consensusfamily",
            "genus": "Knownvirus",
            "species": "Knownvirus species",
        },
        {
            "ident": "WVDB|contig_without_species",
            "superkingdom": "Viruses",
            "phylum": "Fallbackonlyphylum",
            "class": "Fallbackonlyclass",
            "order": "Fallbackonlyorder",
            "family": "Fallbackonlyfamily",
            "genus": "",
            "species": "",
        },
    ]


def test_manifest_coverage_ignores_non_wvdb_references(tmp_path: Path) -> None:
    """check-manifest-coverage succeeds when WVDB rows are covered."""
    manifest = tmp_path / "manifest.csv"
    lineages = tmp_path / "lineages.csv"

    manifest.write_text(
        "# SOURMASH-MANIFEST-VERSION: 1.0\n"
        "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\n"
        "signatures/a.sig.gz,abc,abc,31,DNA,0,50,10,False,NC_005336.1,ref.fa\n"
        "signatures/b.sig.gz,def,def,31,DNA,0,50,10,False,WVDB|contig_1,ref.fa\n",
        encoding="utf-8",
    )
    lineages.write_text(
        "ident,superkingdom,phylum,class,order,family,genus,species\n"
        "WVDB|contig_1,Viruses,,,,,,\n",
        encoding="utf-8",
    )

    result = run_script(
        "check-manifest-coverage",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 0, result.stderr


def test_manifest_coverage_accepts_prepared_sourmash_taxonomy(
    tmp_path: Path,
) -> None:
    """check-manifest-coverage accepts sourmash tax prepare output."""
    manifest = tmp_path / "manifest.csv"
    lineages = tmp_path / "lineages.csv"

    manifest.write_text(
        "# SOURMASH-MANIFEST-VERSION: 1.0\n"
        "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\n"
        "signatures/b.sig.gz,def,def,31,DNA,0,50,10,False,WVDB|contig_1,ref.fa\n",
        encoding="utf-8",
    )
    lineages.write_text(
        "identifiers,superkingdom,phylum,class,order,family,genus,species,strain\n"
        "WVDB|contig_1,Viruses,,,,,,,\n",
        encoding="utf-8",
    )

    result = run_script(
        "check-manifest-coverage",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 0, result.stderr


def test_manifest_coverage_fails_on_missing_wvdb_lineage(tmp_path: Path) -> None:
    """check-manifest-coverage fails when a WVDB signature lacks a lineage."""
    manifest = tmp_path / "manifest.csv"
    lineages = tmp_path / "lineages.csv"

    manifest.write_text(
        "# SOURMASH-MANIFEST-VERSION: 1.0\n"
        "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\n"
        "signatures/b.sig.gz,def,def,31,DNA,0,50,10,False,WVDB|missing_contig,ref.fa\n",
        encoding="utf-8",
    )
    lineages.write_text(
        "ident,superkingdom,phylum,class,order,family,genus,species\n"
        "WVDB|other_contig,Viruses,,,,,,\n",
        encoding="utf-8",
    )

    result = run_script(
        "check-manifest-coverage",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "WVDB signatures are missing lineages" in result.stderr
    assert "WVDB|missing_contig" in result.stderr


def test_curate_lineages_most_specified_rewrites_conflicting_parents(
    tmp_path: Path,
) -> None:
    """The uniquely most-specified lineage should become canonical with a warning."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1 reference title", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            (
                "NCBI_1",
                "Viruses",
                "unclassified Viruses phylum",
                "unclassified Viruses class",
                "unclassified Viruses order",
                "unclassified Viruses family",
                "unclassified Viruses genus",
                "Shared virus",
                "ncbi-strain",
            ),
        ],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Picornaviridae",
                "",
                "Shared_virus",
                "wvdb-strain",
            ),
        ],
    )

    result, output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
        policy="most-specified",
    )

    assert result.returncode == 0, result.stderr
    output_rows = read_csv_rows(output)
    rows = {row["ident"]: row for row in output_rows}
    assert [row["ident"] for row in output_rows] == ["NCBI_1", "WVDB|1"]
    assert rows["NCBI_1"]["phylum"] == "Pisuviricota"
    assert rows["NCBI_1"]["species"] == "Shared_virus"
    assert rows["NCBI_1"]["strain"] == "ncbi-strain"
    assert rows["WVDB|1"]["strain"] == "wvdb-strain"
    assert "resolved 1 taxonomy conflict" in result.stderr
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert {row["normalized_name"] for row in diagnostic_rows} == {"shared virus"}
    assert {row["decision"] for row in diagnostic_rows} == {"winner", "rewritten"}


def test_curate_lineages_source_policy_overrides_specificity(tmp_path: Path) -> None:
    """An explicit NCBI policy should win even when WVDB is more specified."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            (
                "NCBI_1",
                "Viruses",
                "unclassified Viruses phylum",
                "",
                "",
                "",
                "",
                "Shared virus",
                "",
            ),
        ],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Picornaviridae",
                "",
                "Shared_virus",
                "",
            ),
        ],
    )

    for policy, expected_phylum, expected_species in (
        ("ncbi-wins", "unclassified Viruses phylum", "Shared virus"),
        ("wvdb-wins", "Pisuviricota", "Shared_virus"),
    ):
        result, output, _diagnostics = run_curation(
            tmp_path,
            manifests=(reference_manifest, wvdb_manifest),
            lineages=(reference_lineages, wvdb_lineages),
            policy=policy,
        )

        assert result.returncode == 0, result.stderr
        rows = {row["ident"]: row for row in read_csv_rows(output)}
        assert rows["WVDB|1"]["phylum"] == expected_phylum
        assert rows["WVDB|1"]["species"] == expected_species


def test_curate_lineages_fails_on_equally_specified_conflicts(
    tmp_path: Path,
) -> None:
    """Most-specified must not break incompatible ties by input order or frequency."""
    reference_manifest = write_manifest(tmp_path / "reference.manifest.csv", [])
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "a" * 32), ("WVDB|2", "b" * 32)],
    )
    reference_lineages = write_lineages(tmp_path / "reference.lineages.csv", [])
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Picornaviridae",
                "",
                "Shared_virus",
                "",
            ),
            (
                "WVDB|2",
                "Viruses",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Dicistroviridae",
                "",
                "Shared_virus",
                "",
            ),
        ],
    )

    result, _output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
        policy="most-specified",
    )

    assert result.returncode == 1
    assert "ambiguous taxonomy conflict" in result.stderr
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert {row["conflict_type"] for row in diagnostic_rows} == {"taxonomy_name"}
    assert {row["decision"] for row in diagnostic_rows} == {"ambiguous"}


def test_curate_lineages_fails_on_duplicate_signature_identity(
    tmp_path: Path,
) -> None:
    """Duplicate sketch identities remain a hard build error with evidence."""
    digest = "c" * 32
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", digest)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", digest)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [("NCBI_1", "Viruses", "Reference phylum", "", "", "", "", "Virus", "")],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [("WVDB|1", "Viruses", "WVDB phylum", "", "", "", "", "Virus", "")],
    )

    result, _output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
        policy="most-specified",
    )

    assert result.returncode == 1
    assert digest in result.stderr
    assert "NCBI_1" in result.stderr
    assert "WVDB|1" in result.stderr
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert [row["ident"] for row in diagnostic_rows] == ["NCBI_1", "WVDB|1"]
    assert {row["conflict_type"] for row in diagnostic_rows} == {"signature_identity"}


def test_curate_lineages_removes_rows_absent_from_manifests(tmp_path: Path) -> None:
    """Published taxonomy should contain only identifiers retained in the collection."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(tmp_path / "wvdb.manifest.csv", [])
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            ("NCBI_1", "Viruses", "Pisuviricota", "", "", "", "", "Virus 1", ""),
            ("NCBI_extra", "Viruses", "Pisuviricota", "", "", "", "", "Virus 2", ""),
        ],
    )
    wvdb_lineages = write_lineages(tmp_path / "wvdb.lineages.csv", [])

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
        policy="most-specified",
    )

    assert result.returncode == 0, result.stderr
    assert [row["ident"] for row in read_csv_rows(output)] == ["NCBI_1"]


def test_signature_identity_includes_sketch_parameters(tmp_path: Path) -> None:
    """Matching digests with different sketch parameters are distinct identities."""
    digest = "d" * 32
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", digest)],
        scaled="50",
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", digest)],
        scaled="100",
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [("NCBI_1", "Viruses", "Phylum A", "", "", "", "", "Virus A", "")],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [("WVDB|1", "Viruses", "Phylum B", "", "", "", "", "Virus B", "")],
    )

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
        policy="most-specified",
    )

    assert result.returncode == 0, result.stderr
    assert [row["ident"] for row in read_csv_rows(output)] == ["NCBI_1", "WVDB|1"]


def test_curate_lineages_rejects_duplicate_identifiers(tmp_path: Path) -> None:
    """Lineage identifiers must be unique across both taxonomy sources."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("SHARED", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [("SHARED", "Viruses", "Phylum A", "", "", "", "", "Virus A", "")],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            ("SHARED", "Viruses", "Phylum B", "", "", "", "", "Virus B", ""),
            ("WVDB|1", "Viruses", "Phylum B", "", "", "", "", "Virus C", ""),
        ],
    )

    result, _output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
        policy="most-specified",
    )

    assert result.returncode == 1
    assert "Duplicate lineage identifier: SHARED" in result.stderr


def test_curate_lineages_rejects_cross_rank_name_reuse(tmp_path: Path) -> None:
    """Taxburst requires node names to be globally unique across taxonomy ranks."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            (
                "NCBI_1",
                "Viruses",
                "Phylum A",
                "Class A",
                "Order A",
                "Family A",
                "Repeated taxon",
                "Species A",
                "",
            ),
        ],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Phylum B",
                "Class B",
                "Order B",
                "Family B",
                "Genus B",
                "Repeated_taxon",
                "",
            ),
        ],
    )

    result, _output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
        policy="most-specified",
    )

    assert result.returncode == 1
    assert "ambiguous taxonomy conflict" in result.stderr
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert {row["conflict_type"] for row in diagnostic_rows} == {"taxonomy_rank"}
