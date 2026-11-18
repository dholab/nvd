"""Build the synthetic wastewater metagenome fixtures for the NVD e2e bundle.

This is a maintainer utility, not pipeline runtime code. It simulates Illumina
paired-end reads with InSilicoSeq from published RefSeq genomes, so the test
dataset is defined by a short recipe in this file rather than by opaque FASTQ
files committed to the repository.

Why synthetic reads at all
--------------------------
NVD's purpose is finding a low-abundance virus inside a large environmental
metagenome. Real public data cannot give us that in a small file:

* Hybrid-capture libraries (every SRA run the e2e used before this) are 97-100%
  target virus, so they test the pipeline's plumbing and nothing else.
* A low-abundance sample is inherently large. Assembly needs roughly 10-30x
  depth, so the reads a target needs scale with its genome length, while the
  total reads scale as "target reads / target fraction". Orf virus is 139,962 bp
  and monkeypox is about 197,000 bp, so a 1%-abundance Orf sample needs on the
  order of a million read pairs - over 100 MB.

Simulating solves both problems at once: we choose the abundance, we know the
ground truth exactly, and the recipe is a few hundred lines instead of gigabytes.

Where the community composition comes from
------------------------------------------
The two communities below are modelled on a real 18-sample wastewater run
(``rapid-ww-metagenomics`` run ``rww_20261029_32720``, Illumina NovaSeq X). The
properties reproduced here are:

* Plant viruses dominate the virome. Tobamoviruses (pepper mild mottle virus and
  tobacco mosaic virus) were the most abundant and most frequently detected
  species in that run.
* Human enteric viruses sit below them: astrovirus, cosavirus, norovirus and
  mastadenovirus were all called repeatedly.
* crAssphage-type gut phage is a reliable sewage component.
* Large DNA viruses appear only as traces. That run detected Poxviridae
  (molluscum contagiosum virus) in 6 of 18 samples on a total of 21 reads,
  covering at most 0.4% of the genome. The high-complexity community below puts
  Orf virus and monkeypox at the same kind of level deliberately.
* Bacteria make up the bulk of what remains after host and rRNA depletion.

Determinism
-----------
InSilicoSeq produces identical reads for a given seed, but Python's gzip writes
a modification time into the header, so the compressed bytes differ between
runs. We therefore let InSilicoSeq write plain FASTQ and compress it ourselves
with the timestamp zeroed, which makes the published ``.fastq.gz`` files
byte-reproducible and lets one sha256 serve both download integrity and rebuild
verification.

Regenerate with::

    pixi run -e dev python tests/scripts/build_synthetic_metagenomes.py \\
        --out-dir .scratch/e2e-bundle/synthetic
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]
LOGGER = logging.getLogger(__name__)

#: Seed for every random draw, so a rebuild reproduces the published reads.
DEFAULT_SEED = 42

#: InSilicoSeq error model. The real run this is modelled on was NovaSeq X.
DEFAULT_ERROR_MODEL = "NovaSeq"

#: Read length the NovaSeq error model emits, used for the coverage arithmetic.
NOVASEQ_READ_LENGTH = 151

#: Gzip level for the published files. Only the level affects size; the
#: timestamp is always zeroed, see ``compress_deterministically``.
GZIP_LEVEL = 9

NCBI_EFETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={accession}&rettype=fasta&retmode=text"
)

#: NCBI allows roughly three requests a second without an API key, and answers
#: HTTP 429 above that. Pause between fetches and back off when throttled.
NCBI_REQUEST_DELAY_SECONDS = 0.5
NCBI_MAX_ATTEMPTS = 5

#: HTTP status NCBI answers with when the request rate is too high.
HTTP_TOO_MANY_REQUESTS = 429


@dataclass(frozen=True)
class ReferenceGenome:
    """A published genome that reads are simulated from.

    Attributes:
        accession: RefSeq nucleotide accession, also the FASTA record id and the
            key InSilicoSeq matches read counts against.
        taxon_name: Scientific name, for the documentation table.
        role: What this genome represents in the community. One of ``target``,
            ``plant_virus``, ``enteric_virus``, ``phage``, ``host`` or
            ``bacterium``.
        genome_length: Length in bases, verified against NCBI. Used to turn read
            counts into expected coverage without downloading anything.
        taxid: NCBI taxonomy id. Only the targets need one, because only they
            are expected to produce BLAST hits against the mini database.
    """

    accession: str
    taxon_name: str
    role: str
    genome_length: int
    taxid: int | None = None


#: Every genome either community draws from. Lengths were verified by fetching
#: each record from NCBI; see ``verify_declared_genome_lengths``.
REFERENCE_GENOMES: Mapping[str, ReferenceGenome] = {
    # --- targets: the only viruses in the e2e's mini BLAST database ---------
    "NC_005336.1": ReferenceGenome(
        accession="NC_005336.1",
        taxon_name="Orf virus",
        role="target",
        genome_length=139_962,
        taxid=10_258,
    ),
    "NC_063383.1": ReferenceGenome(
        accession="NC_063383.1",
        taxon_name="Monkeypox virus",
        role="target",
        genome_length=197_209,
        taxid=10_244,
    ),
    # --- plant viruses: the dominant signal in real wastewater --------------
    "NC_003630.1": ReferenceGenome(
        accession="NC_003630.1",
        taxon_name="Pepper mild mottle virus",
        role="plant_virus",
        genome_length=6_357,
    ),
    "NC_001367.1": ReferenceGenome(
        accession="NC_001367.1",
        taxon_name="Tobacco mosaic virus",
        role="plant_virus",
        genome_length=6_395,
    ),
    # --- human enteric viruses ---------------------------------------------
    "NC_039477.1": ReferenceGenome(
        accession="NC_039477.1",
        taxon_name="Norovirus GII",
        role="enteric_virus",
        genome_length=7_567,
    ),
    "NC_001943.1": ReferenceGenome(
        accession="NC_001943.1",
        taxon_name="Human astrovirus",
        role="enteric_virus",
        genome_length=6_813,
    ),
    "NC_012800.1": ReferenceGenome(
        accession="NC_012800.1",
        taxon_name="Cosavirus A",
        role="enteric_virus",
        genome_length=7_632,
    ),
    "NC_001405.1": ReferenceGenome(
        accession="NC_001405.1",
        taxon_name="Human adenovirus C",
        role="enteric_virus",
        genome_length=35_937,
    ),
    # --- gut phage ----------------------------------------------------------
    "NC_024711.1": ReferenceGenome(
        accession="NC_024711.1",
        taxon_name="Carjivirus communis (crAssphage)",
        role="phage",
        genome_length=97_065,
    ),
    # --- host ---------------------------------------------------------------
    "NC_012920.1": ReferenceGenome(
        accession="NC_012920.1",
        taxon_name="Homo sapiens mitochondrion",
        role="host",
        genome_length=16_569,
    ),
    # --- bacterial bulk -----------------------------------------------------
    "NC_000913.3": ReferenceGenome(
        accession="NC_000913.3",
        taxon_name="Escherichia coli K-12 MG1655",
        role="bacterium",
        genome_length=4_641_652,
    ),
    "NC_002516.2": ReferenceGenome(
        accession="NC_002516.2",
        taxon_name="Pseudomonas aeruginosa PAO1",
        role="bacterium",
        genome_length=6_264_404,
    ),
}


@dataclass(frozen=True)
class SyntheticCommunity:
    """One simulated sample: which genomes, and how many reads from each.

    Attributes:
        sample_id: Samplesheet id, also the output file prefix.
        complexity: ``low`` or ``high``. Low complexity means few genomes and a
            target deep enough to assemble; high complexity means many genomes
            and targets at trace level.
        description: One sentence for the README provenance table.
        read_counts: Accession to **total** reads from that genome. InSilicoSeq
            counts both mates, so 1,000 here means 500 pairs.
        target_accessions: The genomes the e2e asserts BLAST hits for.
    """

    sample_id: str
    complexity: str
    description: str
    read_counts: Mapping[str, int]
    target_accessions: tuple[str, ...] = field(default=())


#: The low-complexity community: five genomes, Orf deep enough that SPAdes
#: produces contigs. This covers the assembly route.
LOW_COMPLEXITY_COMMUNITY = SyntheticCommunity(
    sample_id="synthetic_ww_lowcx",
    complexity="low",
    description=(
        "Low-complexity synthetic wastewater: five genomes with Orf virus at "
        "about 6x depth, so assembly and contig BLAST are exercised."
    ),
    read_counts={
        "NC_005336.1": 6_000,  # Orf virus, ~6.5x over 140 kb
        "NC_000913.3": 60_000,  # E. coli bulk
        "NC_003630.1": 16_000,  # pepper mild mottle virus
        "NC_024711.1": 12_000,  # crAssphage
        "NC_012920.1": 6_000,  # human mitochondrion
    },
    target_accessions=("NC_005336.1",),
)

#: The high-complexity community: twelve genomes with the composition of a real
#: sewage virome, and both targets at trace level. This covers detection from
#: unassembled reads, which is how a large DNA virus actually shows up in
#: wastewater.
HIGH_COMPLEXITY_COMMUNITY = SyntheticCommunity(
    sample_id="synthetic_ww_highcx",
    complexity="high",
    description=(
        "High-complexity synthetic wastewater: twelve genomes modelled on a "
        "real NovaSeq X sewage run, with Orf virus and monkeypox at under 0.5% "
        "of reads, the level real poxviruses appear at."
    ),
    read_counts={
        "NC_005336.1": 500,  # Orf virus, ~0.5x: detectable, not assemblable
        "NC_063383.1": 400,  # monkeypox clade IIb, ~0.3x
        "NC_003630.1": 12_000,  # pepper mild mottle virus, dominant
        "NC_001367.1": 6_000,  # tobacco mosaic virus
        "NC_039477.1": 3_000,  # norovirus GII
        "NC_001943.1": 2_000,  # human astrovirus
        "NC_012800.1": 2_000,  # cosavirus A
        "NC_001405.1": 2_500,  # human adenovirus C
        "NC_024711.1": 8_000,  # crAssphage
        "NC_012920.1": 3_000,  # human mitochondrion
        "NC_000913.3": 40_000,  # E. coli bulk
        "NC_002516.2": 30_000,  # Pseudomonas aeruginosa bulk
    },
    target_accessions=("NC_005336.1", "NC_063383.1"),
)

COMMUNITIES: tuple[SyntheticCommunity, ...] = (
    LOW_COMPLEXITY_COMMUNITY,
    HIGH_COMPLEXITY_COMMUNITY,
)


# ---------------------------------------------------------------------------
# Planning: pure functions over the declarations above
# ---------------------------------------------------------------------------


def community_by_sample_id(sample_id: str) -> SyntheticCommunity:
    """Return the declared community with this sample id.

    Raises:
        KeyError: if no community declares that sample id.
    """
    for community in COMMUNITIES:
        if community.sample_id == sample_id:
            return community
    message = f"{sample_id} is not a declared synthetic community"
    raise KeyError(message)


def total_simulated_reads(community: SyntheticCommunity) -> int:
    """Total reads across every genome, counting both mates of a pair."""
    return sum(community.read_counts.values())


def target_read_fraction(community: SyntheticCommunity, accession: str) -> float:
    """Share of the sample's reads that come from one genome."""
    return community.read_counts[accession] / total_simulated_reads(community)


def expected_target_coverage(
    community: SyntheticCommunity,
    accession: str,
    read_length: int = NOVASEQ_READ_LENGTH,
) -> float:
    """Expected depth over a genome, before deduplication.

    Simple bases-over-length arithmetic. It is what decides whether a target can
    assemble: below roughly 5x, SPAdes will not build a usable contig, which is
    why the low- and high-complexity communities are deliberately on opposite
    sides of that line.
    """
    genome = REFERENCE_GENOMES[accession]
    simulated_bases = community.read_counts[accession] * read_length
    return simulated_bases / genome.genome_length


def write_readcount_file(
    community: SyntheticCommunity,
    destination: Path,
) -> Path:
    """Write the InSilicoSeq ``--readcount_file``: ``accession<TAB>reads``."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        f"{accession}\t{reads}" for accession, reads in community.read_counts.items()
    ]
    destination.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return destination


# ---------------------------------------------------------------------------
# Reference genomes
# ---------------------------------------------------------------------------


def fetch_reference_fasta(accession: str, cache_dir: Path) -> Path:
    """Download one genome from NCBI, reusing an already-downloaded copy."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = cache_dir / f"{accession}.fasta"
    if destination.is_file() and destination.stat().st_size > 0:
        LOGGER.debug("Reusing cached %s", destination)
        return destination

    url = NCBI_EFETCH_URL.format(accession=accession)
    LOGGER.info("Fetching %s from NCBI", accession)
    partial = destination.with_suffix(".fasta.part")

    payload = b""
    for attempt in range(1, NCBI_MAX_ATTEMPTS + 1):
        time.sleep(NCBI_REQUEST_DELAY_SECONDS)
        try:
            with urllib.request.urlopen(url) as response:  # noqa: S310 - fixed URL
                payload = response.read()
        except urllib.error.HTTPError as error:
            # 429 means we asked too fast; anything else is not worth retrying.
            if error.code != HTTP_TOO_MANY_REQUESTS or attempt == NCBI_MAX_ATTEMPTS:
                raise
            backoff = NCBI_REQUEST_DELAY_SECONDS * 2**attempt
            LOGGER.warning(
                "NCBI throttled %s (attempt %d); retrying in %.1fs",
                accession,
                attempt,
                backoff,
            )
            time.sleep(backoff)
        else:
            break

    partial.write_bytes(payload)
    if not partial.read_text(encoding="utf-8").startswith(">"):
        partial.unlink(missing_ok=True)
        message = f"NCBI did not return FASTA for {accession}"
        raise RuntimeError(message)
    partial.rename(destination)
    return destination


def measure_fasta_length(fasta: Path) -> int:
    """Total sequence length in a FASTA file, ignoring headers."""
    return sum(
        len(line.strip())
        for line in fasta.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )


def verify_declared_genome_lengths(accessions: Iterable[str], cache_dir: Path) -> None:
    """Fail loudly if a downloaded genome is not the length we declared.

    The coverage arithmetic, and therefore the whole design of the two
    communities, rests on these lengths. If NCBI revises a record we want a
    clear error rather than a silently wrong abundance.
    """
    for accession in accessions:
        fasta = fetch_reference_fasta(accession, cache_dir)
        observed = measure_fasta_length(fasta)
        declared = REFERENCE_GENOMES[accession].genome_length
        if observed != declared:
            message = (
                f"{accession} is {observed} bp but REFERENCE_GENOMES declares "
                f"{declared} bp. Update the declaration and re-check the "
                f"community read budgets."
            )
            raise RuntimeError(message)


def concatenate_reference_fastas(
    accessions: Sequence[str],
    cache_dir: Path,
    destination: Path,
) -> Path:
    """Write one multi-FASTA holding every genome a community draws from."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        for accession in accessions:
            fasta = fetch_reference_fasta(accession, cache_dir)
            handle.write(fasta.read_text(encoding="utf-8").rstrip("\n") + "\n")
    return destination


# ---------------------------------------------------------------------------
# Simulation and deterministic packaging
# ---------------------------------------------------------------------------


# Six parameters, but each is an independent InSilicoSeq input; grouping them
# into a config object would add a layer without removing a decision.
def build_iss_command(  # noqa: PLR0913
    genomes_fasta: Path,
    readcount_file: Path,
    output_prefix: Path,
    *,
    seed: int = DEFAULT_SEED,
    error_model: str = DEFAULT_ERROR_MODEL,
    cpus: int = 4,
) -> list[str]:
    """Return the InSilicoSeq command line, with every path made absolute.

    The command runs with its working directory set to the output directory, so
    that InSilicoSeq's scratch files land there rather than wherever the builder
    was invoked from. That means a relative path from the caller would no longer
    resolve once the process starts, so each one is resolved here.
    """
    return [
        "iss",
        "generate",
        "--genomes",
        str(genomes_fasta.resolve()),
        "--readcount_file",
        str(readcount_file.resolve()),
        "--model",
        error_model,
        "--seed",
        str(seed),
        "--cpus",
        str(cpus),
        "--quiet",
        "--output",
        str(output_prefix.resolve()),
    ]


def simulate_community_reads(  # noqa: PLR0913
    genomes_fasta: Path,
    readcount_file: Path,
    output_prefix: Path,
    *,
    seed: int = DEFAULT_SEED,
    error_model: str = DEFAULT_ERROR_MODEL,
    cpus: int = 4,
) -> tuple[Path, Path]:
    """Run InSilicoSeq and return the plain R1/R2 FASTQ paths.

    Deliberately does not pass ``--compress``: we compress afterwards with the
    timestamp zeroed so the published files are byte-reproducible.
    """
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    resolved_prefix = output_prefix.resolve()
    command = build_iss_command(
        genomes_fasta,
        readcount_file,
        output_prefix,
        seed=seed,
        error_model=error_model,
        cpus=cpus,
    )
    LOGGER.info("Simulating reads: %s", " ".join(command))
    subprocess.run(command, check=True, cwd=resolved_prefix.parent)  # noqa: S603

    output_prefix = resolved_prefix
    read_one = output_prefix.with_name(f"{output_prefix.name}_R1.fastq")
    read_two = output_prefix.with_name(f"{output_prefix.name}_R2.fastq")
    for path in (read_one, read_two):
        if not path.is_file():
            message = f"InSilicoSeq did not write {path}"
            raise RuntimeError(message)

    # InSilicoSeq leaves per-cpu temporary VCFs behind even without
    # --store_mutations; drop them so the output directory is just the reads.
    for leftover in output_prefix.parent.glob(f"{output_prefix.name}.iss.tmp.*"):
        leftover.unlink()

    return read_one, read_two


def compress_deterministically(
    source: Path,
    destination: Path,
    level: int = GZIP_LEVEL,
) -> Path:
    """Gzip a file with the header timestamp zeroed.

    A gzip header carries two fields that vary between otherwise identical
    runs: the modification time and the original filename. Both change the
    compressed bytes, and therefore the sha256, even when the reads are the
    same. ``mtime=0`` with an empty ``filename`` is the library equivalent of
    ``gzip -n``.

    Note that the filename has to be suppressed explicitly: passing a path to
    ``GzipFile(filename=...)`` writes that name into the header, so compressing
    through differently-named temporary files would otherwise produce different
    bytes for the same reads.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    with (
        source.open("rb") as raw,
        partial.open("wb") as sink,
        gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=sink,
            compresslevel=level,
            mtime=0,
        ) as compressed,
    ):
        shutil.copyfileobj(raw, compressed)
    partial.rename(destination)
    return destination


def sha256_of_file(path: Path) -> str:
    """Hex sha256 of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_fastq_reads(path: Path) -> int:
    """Number of FASTQ records, accepting plain or gzipped input."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        return sum(1 for index, _ in enumerate(handle) if index % 4 == 0)


def tool_versions() -> dict[str, str]:
    """Record the versions that produced the reads, for the manifest."""
    versions: dict[str, str] = {}
    for name, command in (("insilicoseq", ["iss", "--version"]),):
        try:
            completed = subprocess.run(  # noqa: S603
                command,
                check=True,
                capture_output=True,
                text=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            versions[name] = "unknown"
        else:
            versions[name] = (completed.stdout or completed.stderr).strip()
    return versions


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------


# Six parameters: the community plus the three directories and two knobs the
# CLI exposes. All are keyword-only below, so call sites stay readable.
def build_community(  # noqa: PLR0913
    community: SyntheticCommunity,
    *,
    cache_dir: Path,
    out_dir: Path,
    work_dir: Path,
    seed: int = DEFAULT_SEED,
    cpus: int = 4,
) -> dict[str, Any]:
    """Simulate one community and return its manifest entry."""
    accessions = sorted(community.read_counts)
    verify_declared_genome_lengths(accessions, cache_dir)

    genomes_fasta = concatenate_reference_fastas(
        accessions,
        cache_dir,
        work_dir / f"{community.sample_id}.genomes.fasta",
    )
    readcount_file = write_readcount_file(
        community,
        work_dir / f"{community.sample_id}.readcounts.tsv",
    )
    read_one, read_two = simulate_community_reads(
        genomes_fasta,
        readcount_file,
        work_dir / community.sample_id,
        seed=seed,
        cpus=cpus,
    )

    published: list[dict[str, Any]] = []
    for mate, plain in (("R1", read_one), ("R2", read_two)):
        destination = out_dir / f"{community.sample_id}_{mate}.fastq.gz"
        compress_deterministically(plain, destination)
        published.append(
            {
                "name": destination.name,
                "mate": mate,
                "bytes": destination.stat().st_size,
                "sha256": sha256_of_file(destination),
                "reads": count_fastq_reads(destination),
                "url": None,  # filled in after the Zenodo upload
            },
        )
        plain.unlink()

    composition = [
        {
            "accession": accession,
            "taxon_name": REFERENCE_GENOMES[accession].taxon_name,
            "role": REFERENCE_GENOMES[accession].role,
            "genome_length": REFERENCE_GENOMES[accession].genome_length,
            "reads": reads,
            "read_fraction": round(target_read_fraction(community, accession), 6),
            "expected_coverage": round(
                expected_target_coverage(community, accession),
                3,
            ),
        }
        for accession, reads in sorted(
            community.read_counts.items(),
            key=lambda item: item[1],
            reverse=True,
        )
    ]

    return {
        "sample_id": community.sample_id,
        "platform": "illumina",
        "complexity": community.complexity,
        "description": community.description,
        "method": "insilicoseq",
        "seed": seed,
        "error_model": DEFAULT_ERROR_MODEL,
        "read_length": NOVASEQ_READ_LENGTH,
        "total_reads": total_simulated_reads(community),
        "targets": [
            {
                "accession": accession,
                "taxon_name": REFERENCE_GENOMES[accession].taxon_name,
                "taxid": REFERENCE_GENOMES[accession].taxid,
                "reads": community.read_counts[accession],
                "read_fraction": round(
                    target_read_fraction(community, accession),
                    6,
                ),
                "expected_coverage": round(
                    expected_target_coverage(community, accession),
                    3,
                ),
            }
            for accession in community.target_accessions
        ],
        "composition": composition,
        "files": published,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build the synthetic wastewater metagenome fixtures for the NVD "
            "end-to-end test bundle."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / ".scratch" / "e2e-bundle" / "synthetic",
        help="Where the published .fastq.gz files are written.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=REPO_ROOT / ".scratch" / "e2e-bundle" / "genomes",
        help="Where downloaded RefSeq genomes are cached.",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=REPO_ROOT / ".scratch" / "e2e-bundle" / "work",
        help="Scratch space for multi-FASTAs and uncompressed reads.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument(
        "--only",
        action="append",
        default=None,
        metavar="SAMPLE_ID",
        help="Build just this community; repeatable.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    selected = (
        COMMUNITIES
        if not args.only
        else tuple(community_by_sample_id(sample_id) for sample_id in args.only)
    )

    for directory in (args.out_dir, args.cache_dir, args.work_dir):
        directory.mkdir(parents=True, exist_ok=True)

    entries = [
        build_community(
            community,
            cache_dir=args.cache_dir,
            out_dir=args.out_dir,
            work_dir=args.work_dir,
            seed=args.seed,
            cpus=args.cpus,
        )
        for community in selected
    ]

    manifest = {
        "description": (
            "Synthetic wastewater metagenome fixtures for the NVD end-to-end "
            "test, simulated with InSilicoSeq from published RefSeq genomes."
        ),
        "seed": args.seed,
        "tool_versions": tool_versions(),
        "samples": entries,
    }
    manifest_path = args.out_dir / "synthetic_metagenomes.manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    LOGGER.info("Wrote %s", manifest_path)

    for entry in entries:
        total_bytes = sum(item["bytes"] for item in entry["files"])
        LOGGER.info(
            "%s: %d reads, %.1f MB, targets %s",
            entry["sample_id"],
            entry["total_reads"],
            total_bytes / 1_000_000,
            ", ".join(
                f"{target['taxon_name']} {target['read_fraction'] * 100:.2f}%"
                for target in entry["targets"]
            ),
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
