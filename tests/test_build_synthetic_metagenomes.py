"""Tests for the synthetic wastewater metagenome fixture builder.

These cover the pure planning functions (read budgets, coverage arithmetic,
readcount files, deterministic compression). The InSilicoSeq invocation itself
is exercised by the ``slow`` test at the bottom, because simulating reads takes
tens of seconds and needs the ``dev`` Pixi environment.
"""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import pytest

from tests.scripts.build_synthetic_metagenomes import (
    COMMUNITIES,
    REFERENCE_GENOMES,
    SyntheticCommunity,
    build_iss_command,
    community_by_sample_id,
    compress_deterministically,
    count_fastq_reads,
    expected_target_coverage,
    sha256_of_file,
    target_read_fraction,
    total_simulated_reads,
    write_readcount_file,
)

# A tiny stand-in community so the arithmetic tests do not depend on the real
# read budgets, which we expect to be tuned over time.
EXAMPLE_COMMUNITY = SyntheticCommunity(
    sample_id="example_sample",
    complexity="low",
    description="Two-genome example used by the unit tests.",
    read_counts={"NC_005336.1": 1_000, "NC_003630.1": 9_000},
    target_accessions=("NC_005336.1",),
)


def test_total_simulated_reads_sums_every_genome() -> None:
    assert total_simulated_reads(EXAMPLE_COMMUNITY) == 10_000


def test_target_read_fraction_is_share_of_all_reads() -> None:
    assert target_read_fraction(EXAMPLE_COMMUNITY, "NC_005336.1") == pytest.approx(0.1)


def test_expected_target_coverage_uses_genome_length() -> None:
    # 1,000 reads x 151 bp over the 139,962 bp Orf genome.
    coverage = expected_target_coverage(
        EXAMPLE_COMMUNITY,
        "NC_005336.1",
        read_length=151,
    )
    assert coverage == pytest.approx(1_000 * 151 / 139_962, rel=1e-6)


def test_write_readcount_file_emits_one_tab_separated_row_per_genome(
    tmp_path: Path,
) -> None:
    destination = write_readcount_file(
        EXAMPLE_COMMUNITY,
        tmp_path / "counts.tsv",
    )
    rows = destination.read_text(encoding="utf-8").strip().splitlines()
    assert sorted(rows) == sorted(
        ["NC_005336.1\t1000", "NC_003630.1\t9000"],
    )


def test_compress_deterministically_is_byte_identical_across_runs(
    tmp_path: Path,
) -> None:
    source = tmp_path / "reads.fastq"
    source.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")

    first = compress_deterministically(source, tmp_path / "first.fastq.gz")
    second = compress_deterministically(source, tmp_path / "second.fastq.gz")

    assert first.read_bytes() == second.read_bytes()
    assert sha256_of_file(first) == sha256_of_file(second)
    # And it must still be readable gzip holding the original content.
    with gzip.open(first, "rt", encoding="utf-8") as handle:
        assert handle.read() == "@r1\nACGT\n+\nIIII\n"


def test_compress_deterministically_omits_the_gzip_timestamp(
    tmp_path: Path,
) -> None:
    """A recorded sha256 must not change just because the clock moved."""
    source = tmp_path / "reads.fastq"
    source.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    compressed = compress_deterministically(source, tmp_path / "reads.fastq.gz")

    # Bytes 4-8 of a gzip member are the modification time, little-endian.
    assert compressed.read_bytes()[4:8] == b"\x00\x00\x00\x00"


def test_build_iss_command_makes_every_path_absolute(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """InSilicoSeq runs with cwd set elsewhere, so relative paths would break."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "work").mkdir()

    command = build_iss_command(
        Path("work/genomes.fasta"),
        Path("work/counts.tsv"),
        Path("work/sample"),
        seed=42,
    )

    path_arguments = [
        command[command.index(flag) + 1]
        for flag in ("--genomes", "--readcount_file", "--output")
    ]
    assert all(Path(argument).is_absolute() for argument in path_arguments)
    assert command[command.index("--seed") + 1] == "42"


def test_count_fastq_reads_counts_four_line_records(tmp_path: Path) -> None:
    plain = tmp_path / "reads.fastq"
    plain.write_text("@a\nACGT\n+\nIIII\n@b\nTGCA\n+\nIIII\n", encoding="utf-8")
    assert count_fastq_reads(plain) == 2

    compressed = compress_deterministically(plain, tmp_path / "reads.fastq.gz")
    assert count_fastq_reads(compressed) == 2


class TestDeclaredCommunities:
    """Guard rails on the real community definitions."""

    def test_two_communities_are_declared_one_of_each_complexity(self) -> None:
        assert len(COMMUNITIES) == 2
        assert {community.complexity for community in COMMUNITIES} == {"low", "high"}

    def test_every_accession_has_a_declared_reference_genome(self) -> None:
        for community in COMMUNITIES:
            unknown = set(community.read_counts) - set(REFERENCE_GENOMES)
            assert not unknown, f"{community.sample_id} references {unknown}"

    def test_every_target_accession_is_part_of_its_community(self) -> None:
        for community in COMMUNITIES:
            missing = set(community.target_accessions) - set(community.read_counts)
            assert not missing, f"{community.sample_id} lacks reads for {missing}"

    def test_community_sample_ids_are_unique(self) -> None:
        sample_ids = [community.sample_id for community in COMMUNITIES]
        assert len(sample_ids) == len(set(sample_ids))

    def test_low_complexity_community_targets_assembly_depth(self) -> None:
        """The low-complexity sample must assemble, so it needs real depth."""
        community = community_by_sample_id("synthetic_ww_lowcx")
        for accession in community.target_accessions:
            assert expected_target_coverage(community, accession) >= 5.0

    def test_high_complexity_community_keeps_targets_at_trace_level(self) -> None:
        """The high-complexity sample mirrors real sewage: under 1% of reads."""
        community = community_by_sample_id("synthetic_ww_highcx")
        for accession in community.target_accessions:
            assert target_read_fraction(community, accession) < 0.01

    def test_high_complexity_community_is_more_diverse(self) -> None:
        low = community_by_sample_id("synthetic_ww_lowcx")
        high = community_by_sample_id("synthetic_ww_highcx")
        assert len(high.read_counts) > len(low.read_counts)

    def test_community_by_sample_id_rejects_an_unknown_id(self) -> None:
        with pytest.raises(KeyError, match="no_such_sample"):
            community_by_sample_id("no_such_sample")


@pytest.mark.slow
def test_insilicoseq_is_available_and_deterministic(tmp_path: Path) -> None:
    """Pin the behaviour the builder relies on: a seed fixes the reads."""
    genome = tmp_path / "genome.fasta"
    genome.write_text(">NC_TEST.1 test\n" + "ACGTACGTGA" * 200 + "\n", encoding="utf-8")
    counts = tmp_path / "counts.tsv"
    counts.write_text("NC_TEST.1\t100\n", encoding="utf-8")

    outputs = []
    for attempt in ("first", "second"):
        prefix = tmp_path / attempt
        subprocess.run(  # noqa: S603
            [
                "iss",
                "generate",
                "--genomes",
                str(genome),
                "--readcount_file",
                str(counts),
                "--model",
                "NovaSeq",
                "--seed",
                "42",
                "--cpus",
                "1",
                "--quiet",
                "--output",
                str(prefix),
            ],
            check=True,
            cwd=tmp_path,
        )
        outputs.append((prefix.parent / f"{attempt}_R1.fastq").read_bytes())

    assert outputs[0] == outputs[1]
