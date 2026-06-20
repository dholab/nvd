"""Generate the checked-in mini viral integration test fixtures.

This is a maintainer utility, not pipeline runtime code. NVD consumes the
generated FASTA, Deacon index, BLAST database, and samplesheet; it should not
learn how to build custom references during normal runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import subprocess
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = REPO_ROOT / "tests" / "data"
LOGGER = logging.getLogger(__name__)

REFERENCE_ACCESSIONS = (
    "NC_005336.1",  # Orf virus complete genome
    "NC_003310.1",  # Monkeypox virus complete genome, clade 1
    "NC_063383.1",  # Monkeypox virus complete genome, clade 2
)

REFERENCE_TAXIDS = {
    "NC_005336.1": 10258,
    "NC_003310.1": 10244,
    "NC_063383.1": 10244,
}

SRA_RUNS = (
    {
        "sample_id": "orf_virus_ov_pt001_2024",
        "srr": "ERR17363849",
        "platform": "illumina",
        "expected_organism": "Orf virus",
        "taxid": 10258,
        "size_mb": 33,
    },
    {
        "sample_id": "orf_virus_hs_pt001_2025",
        "srr": "ERR17363848",
        "platform": "illumina",
        "expected_organism": "Orf virus",
        "taxid": 10258,
        "size_mb": 89,
    },
    {
        "sample_id": "monkeypox_pt1020_2026",
        "srr": "ERR17356125",
        "platform": "illumina",
        "expected_organism": "Monkeypox virus",
        "taxid": 10244,
        "size_mb": 57,
    },
    {
        "sample_id": "nigerian_orf_srr38321624",
        "srr": "SRR38321624",
        "platform": "illumina",
        "expected_organism": "Orf virus",
        "taxid": 10258,
        "size_mb": 210,
    },
)


@dataclass(frozen=True)
class FixturePaths:
    data_dir: Path
    reference_fasta: Path
    blast_taxid_map: Path
    samplesheet: Path
    manifest: Path
    deacon_index: Path
    blast_prefix: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate mini viral references for NVD integration tests.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Directory where generated fixture files are written.",
    )
    return parser.parse_args()


def fixture_paths(data_dir: Path) -> FixturePaths:
    return FixturePaths(
        data_dir=data_dir,
        reference_fasta=data_dir / "mini_virus_reference.fasta",
        blast_taxid_map=data_dir / "mini_virus_blast_taxid_map.tsv",
        samplesheet=data_dir / "integration_sra_samplesheet.csv",
        manifest=data_dir / "reference.manifest.json",
        deacon_index=data_dir / "mini_virus_deacon.k31w1.idx",
        blast_prefix=data_dir / "mini_virus_blast",
    )


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    LOGGER.info("+ %s", " ".join(command))
    return subprocess.run(  # noqa: S603
        command,
        check=True,
        text=True,
        capture_output=True,
    )


def fetch_text(url: str) -> str:
    request = urllib.request.Request(  # noqa: S310
        url,
        headers={"User-Agent": "nvd-integration-fixture-generator"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
        return response.read().decode("utf-8")


def fetch_reference_fasta(paths: FixturePaths) -> str:
    accession_csv = ",".join(REFERENCE_ACCESSIONS)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession_csv}&rettype=fasta&retmode=text"
    )
    fasta = fetch_text(url)
    missing = [acc for acc in REFERENCE_ACCESSIONS if acc not in fasta]
    if missing:
        message = f"Fetched FASTA did not contain expected accessions: {missing}"
        raise RuntimeError(message)
    paths.reference_fasta.write_text(fasta, encoding="utf-8")
    return url


def write_samplesheet(paths: FixturePaths) -> None:
    lines = ["sample_id,srr,platform,fastq1,fastq2"]
    lines.extend(
        f"{run_info['sample_id']},{run_info['srr']},{run_info['platform']},,"
        for run_info in SRA_RUNS
    )
    paths.samplesheet.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_blast_taxid_map(paths: FixturePaths) -> None:
    lines = [f"{accession}\t{taxid}" for accession, taxid in REFERENCE_TAXIDS.items()]
    paths.blast_taxid_map.write_text("\n".join(lines) + "\n", encoding="utf-8")


def remove_previous_outputs(paths: FixturePaths) -> None:
    for path in (
        paths.reference_fasta,
        paths.samplesheet,
        paths.manifest,
        paths.deacon_index,
    ):
        if path.exists():
            path.unlink()
    for path in paths.data_dir.glob(f"{paths.blast_prefix.name}.*"):
        path.unlink()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def manifest_command(command: list[str]) -> list[str]:
    return [
        repo_relative(Path(part)) if part.startswith(str(REPO_ROOT)) else part
        for part in command
    ]


def collect_files(paths: FixturePaths) -> list[Path]:
    generated = [
        paths.reference_fasta,
        paths.samplesheet,
        paths.deacon_index,
        *sorted(paths.data_dir.glob(f"{paths.blast_prefix.name}.*")),
    ]
    return [path for path in generated if path.exists()]


def command_version(command: list[str]) -> str:
    try:
        completed = run(command)
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        return f"unavailable: {error}"
    return (completed.stdout + completed.stderr).strip()


def build_manifest(
    paths: FixturePaths,
    reference_url: str,
    commands: list[list[str]],
) -> dict[str, Any]:
    files = collect_files(paths)
    return {
        "description": "Mini viral references and SRA samples for NVD integration tests.",
        "maintainer_note": (
            "Regenerate with tests/scripts/generate_integration_fixtures.py; "
            "NVD runtime code should consume these files, not build them."
        ),
        "reference_accessions": list(REFERENCE_ACCESSIONS),
        "reference_fasta_url": reference_url,
        "sra_runs": list(SRA_RUNS),
        "commands": [manifest_command(command) for command in commands],
        "tool_versions": {
            "makeblastdb": command_version(["makeblastdb", "-version"]),
            "blastdbcmd": command_version(["blastdbcmd", "-version"]),
            "deacon": command_version(["deacon", "--version"]),
        },
        "files": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in files
        },
    }


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    paths = fixture_paths(args.data_dir.resolve())
    paths.data_dir.mkdir(parents=True, exist_ok=True)
    remove_previous_outputs(paths)

    reference_url = fetch_reference_fasta(paths)
    write_samplesheet(paths)
    write_blast_taxid_map(paths)

    commands = [
        [
            "makeblastdb",
            "-in",
            str(paths.reference_fasta),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-taxid_map",
            str(paths.blast_taxid_map),
            "-out",
            str(paths.blast_prefix),
        ],
        [
            "blastdbcmd",
            "-db",
            str(paths.blast_prefix),
            "-entry",
            "all",
        ],
        [
            "deacon",
            "index",
            "build",
            "-k",
            "31",
            "-w",
            "1",
            "-o",
            str(paths.deacon_index),
            str(paths.reference_fasta),
        ],
    ]

    run(commands[0])
    run(commands[1])
    run(commands[2])

    manifest = build_manifest(paths, reference_url, commands)
    paths.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    LOGGER.info("Wrote %s", paths.manifest)


if __name__ == "__main__":
    main()
