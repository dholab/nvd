"""Tests for FASTX profiling sidecars."""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from pathlib import Path

import pytest
from profile_fastx import RawProfile, main, parse_bbtools_histogram, run_reformat

FASTQ_RECORD_COUNT = 3
FASTQ_TOTAL_BASES = 33
FASTA_RECORD_COUNT = 3
FASTA_TOTAL_BASES = 210


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_fastq_gz(path: Path, records: list[tuple[str, str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for name, sequence, quality in records:
            handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in records),
        encoding="utf-8",
    )


def test_parse_bbtools_histogram_skips_headers_and_zero_counts(tmp_path: Path) -> None:
    histogram = tmp_path / "histogram.tsv"
    histogram.write_text(
        "#Quality\tcount1\tfraction1\n0\t0\t0.00000\n39\t2\t1.00000\n",
        encoding="utf-8",
    )

    assert parse_bbtools_histogram(histogram, expected_axis="Quality") == Counter(
        {39: 2},
    )


def test_parse_bbtools_histogram_combines_paired_quality_counts(tmp_path: Path) -> None:
    histogram = tmp_path / "histogram.tsv"
    histogram.write_text(
        "#Quality\tcount1\tfraction1\tcount2\tfraction2\n"
        "20\t1\t0.33333\t3\t0.42857\n"
        "39\t2\t0.66667\t4\t0.57143\n",
        encoding="utf-8",
    )

    assert parse_bbtools_histogram(histogram, expected_axis="Quality") == Counter(
        {39: 6, 20: 4},
    )


def test_parse_bbtools_histogram_requires_the_expected_header(tmp_path: Path) -> None:
    missing_header = tmp_path / "missing-header.tsv"
    missing_header.write_text("39\t2\t1.00000\n", encoding="utf-8")
    wrong_axis = tmp_path / "wrong-axis.tsv"
    wrong_axis.write_text("#Length\tCount\n39\t2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="precedes its header"):
        parse_bbtools_histogram(missing_header, expected_axis="Quality")
    with pytest.raises(ValueError, match="expected Quality"):
        parse_bbtools_histogram(wrong_axis, expected_axis="Quality")


def test_run_reformat_profiles_fastq_records_as_singletons(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_file = tmp_path / "reads.fastq.gz"
    input_file.touch()
    observed_command: list[str] = []

    monkeypatch.setattr(
        "profile_fastx.shutil.which",
        lambda _name: "/tools/reformat.sh",
    )

    def fake_run(command: list[str], *, check: bool) -> None:
        observed_command.extend(command)
        assert check is True
        for argument, contents in (
            ("lhist=", "#Length\tCount\n100\t2\n"),
            ("aqhist=", "#Quality\tcount1\tfraction1\n40\t2\t1.00000\n"),
        ):
            output_path = next(
                Path(value.removeprefix(argument))
                for value in command
                if value.startswith(argument)
            )
            output_path.write_text(contents, encoding="utf-8")

    monkeypatch.setattr("profile_fastx.subprocess.run", fake_run)

    profile = run_reformat(input_file, tmpdir=tmp_path)

    assert "out=/dev/null" in observed_command
    assert "int=f" in observed_command
    assert profile == RawProfile(
        length_counts=Counter({100: 2}),
        quality_counts=Counter({40: 2}),
    )


def test_profile_fastx_counts_gzipped_fastq_and_writes_histograms(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastq = tmp_path / "reads.fastq.gz"
    output_prefix = tmp_path / "sample-1.preprocessed_single_read"
    write_fastq_gz(
        fastq,
        [
            ("read-1", "ACGTACGT", "IIIIIIII"),
            ("read-2", "ACGTACGTACGT", "555555555555"),
            ("read-3", "ACGTACGTACGTA", "?????????????"),
        ],
    )
    monkeypatch.setattr(
        "profile_fastx.run_reformat",
        lambda _input_file, **_kwargs: RawProfile(
            length_counts=Counter({8: 1, 12: 1, 13: 1}),
            quality_counts=Counter({40: 1, 20: 1, 30: 1}),
        ),
    )

    main(
        [
            "--sample-id",
            "sample-1",
            "--stage",
            "preprocessed_single_read",
            "--input",
            str(fastq),
            "--output-prefix",
            str(output_prefix),
            "--length-bin-size",
            "5",
            "--quality-bin-size",
            "10",
            "--threshold",
            "min_read_length:length:50",
            "--threshold",
            "min_read_quality:quality:20",
        ],
    )

    profile = read_json(Path(f"{output_prefix}.fastx_profile.json"))
    assert profile["sample_id"] == "sample-1"
    assert profile["stage"] == "preprocessed_single_read"
    assert profile["format"] == "fastq"
    assert profile["sequence_count"] == FASTQ_RECORD_COUNT
    assert profile["total_bases"] == FASTQ_TOTAL_BASES
    assert profile["length"] == {"min": 8, "max": 13, "mean": 11.0, "bin_size": 5}
    assert profile["quality"] == {
        "mean": 30.0,
        "bin_size": 10,
        "encoding": "phred33_integer_average",
    }
    assert profile["thresholds"] == [
        {
            "name": "min_read_length",
            "axis": "length",
            "value": 50.0,
            "sequence_count_at_or_above": 0,
            "bases_at_or_above": 0,
        },
        {
            "name": "min_read_quality",
            "axis": "quality",
            "value": 20.0,
            "sequence_count_at_or_above": 3,
            "bases_at_or_above": None,
        },
    ]
    assert (
        Path(f"{output_prefix}.sequence_count.txt").read_text(encoding="utf-8") == "3\n"
    )
    assert read_tsv(Path(f"{output_prefix}.length_histogram.tsv")) == [
        {
            "sample_id": "sample-1",
            "stage": "preprocessed_single_read",
            "bin_start": "5",
            "bin_end": "9",
            "sequence_count": "1",
        },
        {
            "sample_id": "sample-1",
            "stage": "preprocessed_single_read",
            "bin_start": "10",
            "bin_end": "14",
            "sequence_count": "2",
        },
    ]
    assert read_tsv(Path(f"{output_prefix}.quality_histogram.tsv")) == [
        {
            "sample_id": "sample-1",
            "stage": "preprocessed_single_read",
            "bin_start": "20",
            "bin_end": "29",
            "sequence_count": "1",
        },
        {
            "sample_id": "sample-1",
            "stage": "preprocessed_single_read",
            "bin_start": "30",
            "bin_end": "39",
            "sequence_count": "1",
        },
        {
            "sample_id": "sample-1",
            "stage": "preprocessed_single_read",
            "bin_start": "40",
            "bin_end": "49",
            "sequence_count": "1",
        },
    ]


def test_profile_fastx_counts_sequences_and_bases_above_length_thresholds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastq = tmp_path / "reads.fastq.gz"
    output_prefix = tmp_path / "sample-1.preprocessed_single_read"
    write_fastq_gz(fastq, [("read-1", "AC", "II")])
    monkeypatch.setattr(
        "profile_fastx.run_reformat",
        lambda _input_file, **_kwargs: RawProfile(
            length_counts=Counter({199: 1, 200: 2, 999: 1, 1000: 3, 1200: 1}),
            quality_counts=Counter({40: 8}),
        ),
    )

    main(
        [
            "--sample-id",
            "sample-1",
            "--stage",
            "preprocessed_single_read",
            "--input",
            str(fastq),
            "--output-prefix",
            str(output_prefix),
            "--threshold",
            "metamdbg_min_read_overlap:length:200",
            "--threshold",
            "myloasm_min_read_length:length:1000",
            "--threshold",
            "metaflye_min_read_length:length:1001",
        ],
    )

    profile = read_json(Path(f"{output_prefix}.fastx_profile.json"))
    assert profile["thresholds"] == [
        {
            "name": "metamdbg_min_read_overlap",
            "axis": "length",
            "value": 200.0,
            "sequence_count_at_or_above": 7,
            "bases_at_or_above": 5599,
        },
        {
            "name": "myloasm_min_read_length",
            "axis": "length",
            "value": 1000.0,
            "sequence_count_at_or_above": 4,
            "bases_at_or_above": 4200,
        },
        {
            "name": "metaflye_min_read_length",
            "axis": "length",
            "value": 1001.0,
            "sequence_count_at_or_above": 1,
            "bases_at_or_above": 1200,
        },
    ]


def test_profile_fastx_sums_multiple_fasta_inputs_without_quality(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    output_prefix = tmp_path / "sample-1.filtered_contigs"
    write_fasta(first, [("contig-1", "A" * 80)])
    write_fasta(second, [("contig-2", "C" * 10), ("contig-3", "G" * 120)])
    raw_profiles = {
        first: RawProfile(length_counts=Counter({80: 1}), quality_counts=Counter()),
        second: RawProfile(
            length_counts=Counter({10: 1, 120: 1}),
            quality_counts=Counter(),
        ),
    }
    monkeypatch.setattr(
        "profile_fastx.run_reformat",
        lambda input_file, **_kwargs: raw_profiles[input_file],
    )

    main(
        [
            "--sample-id",
            "sample-1",
            "--stage",
            "filtered_contigs",
            "--input",
            str(first),
            str(second),
            "--output-prefix",
            str(output_prefix),
            "--length-bin-size",
            "50",
        ],
    )

    profile = read_json(Path(f"{output_prefix}.fastx_profile.json"))
    assert profile["format"] == "fasta_or_empty"
    assert profile["sequence_count"] == FASTA_RECORD_COUNT
    assert profile["total_bases"] == FASTA_TOTAL_BASES
    assert profile["quality"] == {"mean": None, "bin_size": 1, "encoding": None}
    assert not Path(f"{output_prefix}.quality_histogram.tsv").exists()
    assert read_tsv(Path(f"{output_prefix}.length_histogram.tsv")) == [
        {
            "sample_id": "sample-1",
            "stage": "filtered_contigs",
            "bin_start": "0",
            "bin_end": "49",
            "sequence_count": "1",
        },
        {
            "sample_id": "sample-1",
            "stage": "filtered_contigs",
            "bin_start": "50",
            "bin_end": "99",
            "sequence_count": "1",
        },
        {
            "sample_id": "sample-1",
            "stage": "filtered_contigs",
            "bin_start": "100",
            "bin_end": "149",
            "sequence_count": "1",
        },
    ]


def test_profile_fastx_uses_bbtools_integer_average_quality_bins(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastq = tmp_path / "reads.fastq.gz"
    output_prefix = tmp_path / "sample-1.preprocessed_single_read"
    write_fastq_gz(fastq, [("read-1", "AC", "IH")])
    monkeypatch.setattr(
        "profile_fastx.run_reformat",
        lambda _input_file, **_kwargs: RawProfile(
            length_counts=Counter({2: 1}),
            quality_counts=Counter({39: 1}),
        ),
    )

    main(
        [
            "--sample-id",
            "sample-1",
            "--stage",
            "preprocessed_single_read",
            "--input",
            str(fastq),
            "--output-prefix",
            str(output_prefix),
        ],
    )

    profile = read_json(Path(f"{output_prefix}.fastx_profile.json"))
    assert profile["quality"] == {
        "mean": 39.0,
        "bin_size": 1,
        "encoding": "phred33_integer_average",
    }
    assert read_tsv(Path(f"{output_prefix}.quality_histogram.tsv")) == [
        {
            "sample_id": "sample-1",
            "stage": "preprocessed_single_read",
            "bin_start": "39",
            "bin_end": "39",
            "sequence_count": "1",
        },
    ]


def test_profile_fastx_defaults_to_input_basename_for_single_input(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastq = tmp_path / "sample-1.single_read.filtered.fastq.gz"
    write_fastq_gz(fastq, [("read-1", "ACGT", "IIII")])
    monkeypatch.setattr(
        "profile_fastx.run_reformat",
        lambda _input_file, **_kwargs: RawProfile(
            length_counts=Counter({4: 1}),
            quality_counts=Counter({40: 1}),
        ),
    )

    main(
        [
            "--sample-id",
            "sample-1",
            "--stage",
            "preprocessed_single_read",
            "--input",
            str(fastq),
        ],
    )

    assert (tmp_path / "sample-1.single_read.filtered.fastx_profile.json").exists()
    assert (tmp_path / "sample-1.single_read.filtered.length_histogram.tsv").exists()
    assert (tmp_path / "sample-1.single_read.filtered.quality_histogram.tsv").exists()
    assert (tmp_path / "sample-1.single_read.filtered.sequence_count.txt").exists()
