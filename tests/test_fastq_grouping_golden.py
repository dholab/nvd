"""Golden seam tests for generated grouped FASTQ inputs."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from py_nvd.cli.commands.samplesheet import (
    FastqGrouping,
    scan_fastq_directory,
    write_samplesheet,
)
from py_nvd.read_inputs import ResolvedSample, read_rows, resolve_rows
from stream_fastqs_to_deacon import DeaconStreamConfig, run_deacon_stream

Runner = Callable[[list[str]], subprocess.CompletedProcess[str]]


def write_fastq_gz(path: Path, read_names: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for name in read_names:
            handle.write(f"@{name}\nACGT\n+\n!!!!\n")
    return path


def command_io(command: list[str]) -> tuple[Path, Path, tuple[Path, ...]]:
    output = Path(command[command.index("--output") + 1])
    summary = Path(command[command.index("--summary") + 1])
    index = command.index("filter") + 1
    while command[index].startswith("--"):
        index += 2
    return output, summary, tuple(Path(item) for item in command[index + 1 :])


def open_maybe_gzip(path: Path, mode: str) -> IO[Any]:
    if path.name.endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def read_fastq_names(path: Path) -> list[str]:
    names: list[str] = []
    with open_maybe_gzip(path, "rt") as handle:
        while True:
            header = handle.readline().strip()
            if not header:
                return names
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                message = f"truncated FASTQ record in {path}"
                raise RuntimeError(message)
            names.append(header.removeprefix("@"))


def paired_read_name(r1: str, r2: str) -> str | None:
    if not r1.endswith("/1") or not r2.endswith("/2"):
        return None
    left = r1.removesuffix("/1")
    right = r2.removesuffix("/2")
    return left if left == right else None


def pair_checking_deacon_runner(report: Path) -> Runner:
    def run(command: list[str]) -> subprocess.CompletedProcess[str]:
        output, summary, inputs = command_io(command)
        if len(inputs) != 2:  # noqa: PLR2004
            return subprocess.CompletedProcess(command, 1)

        r1_records = read_fastq_names(inputs[0])
        r2_records = read_fastq_names(inputs[1])
        if len(r1_records) != len(r2_records):
            return subprocess.CompletedProcess(command, 1)

        pairs: list[list[str]] = []
        for r1, r2 in zip(r1_records, r2_records, strict=True):
            name = paired_read_name(r1, r2)
            if name is None:
                return subprocess.CompletedProcess(command, 1)
            pairs.append([name, name])

        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("", encoding="utf-8")
        summary.write_text(json.dumps({"pairs": pairs}), encoding="utf-8")
        report.write_text(json.dumps(pairs), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    return run


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def stream_config(
    tmp_path: Path,
    *,
    r1: tuple[Path, ...],
    r2: tuple[Path, ...],
) -> DeaconStreamConfig:
    index = tmp_path / "fake-deacon.idx"
    index.write_text("fake index", encoding="utf-8")
    return DeaconStreamConfig(
        sample_id="Sample",
        index=index,
        output=tmp_path / "out" / "Sample.fastq.gz",
        summary=tmp_path / "out" / "Sample.deacon.json",
        reads=(),
        r1=r1,
        r2=r2,
        threads=1,
        abs_threshold=1,
        rel_threshold=0.0,
        deacon_bin="deacon",
    )


@dataclass(frozen=True)
class FastqGroupingObservation:
    warnings: list[str]
    resolution_warnings: tuple[str, ...]
    samples: tuple[ResolvedSample, ...]
    expected_r1: tuple[Path, ...]
    expected_r2: tuple[Path, ...]
    streamed_pairs: list[list[str]]
    expected_pairs: list[list[str]]


def test_generated_illumina_lane_grouping_streams_paired_lanes_in_order(
    tmp_path: Path,
) -> None:
    fastq_dir = tmp_path / "fastqs"
    r1_l002 = write_fastq_gz(
        fastq_dir / "Sample_S1_L002_R1_001.fastq.gz",
        ["L002_A/1"],
    )
    r1_l001 = write_fastq_gz(
        fastq_dir / "Sample_S1_L001_R1_001.fastq.gz",
        ["L001_A/1", "L001_B/1"],
    )
    r2_l002 = write_fastq_gz(
        fastq_dir / "Sample_S1_L002_R2_001.fastq.gz",
        ["L002_A/2"],
    )
    r2_l001 = write_fastq_gz(
        fastq_dir / "Sample_S1_L001_R2_001.fastq.gz",
        ["L001_A/2", "L001_B/2"],
    )

    samples, warnings = scan_fastq_directory(
        fastq_dir,
        sanitize=True,
        group_by=FastqGrouping.ILLUMINA_LANES,
    )
    samplesheet = tmp_path / "samplesheet.csv"
    write_samplesheet(samples, samplesheet, platform="illumina")

    assert warnings == []
    assert read_csv_rows(samplesheet) == [
        {
            "sample_id": "Sample",
            "srr": "",
            "platform": "illumina",
            "fastq1": "",
            "fastq2": "",
            "fastq1_glob": str(
                fastq_dir / "Sample_S1_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz",
            ),
            "fastq2_glob": str(
                fastq_dir / "Sample_S1_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9].fastq.gz",
            ),
        },
    ]

    resolution = resolve_rows(read_rows(samplesheet))
    resolved = resolution.samples[0]

    assert resolved.sample_id == "Sample"
    assert resolved.r1 == (r1_l001, r1_l002)
    assert resolved.r2 == (r2_l001, r2_l002)

    report = tmp_path / "pairs.json"
    run_deacon_stream(
        stream_config(tmp_path, r1=resolved.r1, r2=resolved.r2),
        runner=pair_checking_deacon_runner(report),
    )

    assert json.loads(report.read_text(encoding="utf-8")) == [
        ["L001_A", "L001_A"],
        ["L001_B", "L001_B"],
        ["L002_A", "L002_A"],
    ]


lane_chunk_bundles = st.lists(
    st.tuples(
        st.integers(min_value=1, max_value=999),
        st.integers(min_value=1, max_value=999),
        st.integers(min_value=1, max_value=5),
    ),
    min_size=1,
    max_size=16,
    unique_by=lambda item: (item[0], item[1]),
)

stress_lane_chunk_bundles = st.lists(
    st.tuples(
        st.integers(min_value=1, max_value=999),
        st.integers(min_value=1, max_value=999),
        st.integers(min_value=1, max_value=10),
    ),
    min_size=1,
    max_size=32,
    unique_by=lambda item: (item[0], item[1]),
)


def run_grouped_illumina_lane_case(
    bundles: list[tuple[int, int, int]],
) -> FastqGroupingObservation:
    with tempfile.TemporaryDirectory(prefix="nvd-lane-concat-property.") as tmp:
        case_dir = Path(tmp)
        fastq_dir = case_dir / "fastqs"

        for lane, chunk, record_count in bundles:
            lane_label = f"L{lane:03d}"
            chunk_label = f"{chunk:03d}"
            record_names = [
                f"{lane_label}_{chunk_label}_read{record_index}"
                for record_index in range(record_count)
            ]
            write_fastq_gz(
                fastq_dir / f"Sample_S1_{lane_label}_R1_{chunk_label}.fastq.gz",
                [f"{name}/1" for name in record_names],
            )
            write_fastq_gz(
                fastq_dir / f"Sample_S1_{lane_label}_R2_{chunk_label}.fastq.gz",
                [f"{name}/2" for name in record_names],
            )

        samples, warnings = scan_fastq_directory(
            fastq_dir,
            sanitize=True,
            group_by=FastqGrouping.ILLUMINA_LANES,
        )
        samplesheet = case_dir / "samplesheet.csv"
        write_samplesheet(samples, samplesheet, platform="illumina")

        resolution = resolve_rows(read_rows(samplesheet))
        resolved = resolution.samples[0]

        ordered_bundles = sorted(bundles)
        expected_r1 = tuple(
            fastq_dir / f"Sample_S1_L{lane:03d}_R1_{chunk:03d}.fastq.gz"
            for lane, chunk, _record_count in ordered_bundles
        )
        expected_r2 = tuple(
            fastq_dir / f"Sample_S1_L{lane:03d}_R2_{chunk:03d}.fastq.gz"
            for lane, chunk, _record_count in ordered_bundles
        )
        expected_pairs = [
            [f"L{lane:03d}_{chunk:03d}_read{record_index}"] * 2
            for lane, chunk, record_count in ordered_bundles
            for record_index in range(record_count)
        ]

        assert warnings == []
        assert resolved.sample_id == "Sample"
        assert resolved.platform == "illumina"
        assert resolved.source == "paired_globs"
        assert resolved.r1 == expected_r1
        assert resolved.r2 == expected_r2

        report = case_dir / "pairs.json"
        run_deacon_stream(
            stream_config(case_dir, r1=resolved.r1, r2=resolved.r2),
            runner=pair_checking_deacon_runner(report),
        )

        return FastqGroupingObservation(
            warnings=warnings,
            resolution_warnings=resolution.warnings,
            samples=resolution.samples,
            expected_r1=expected_r1,
            expected_r2=expected_r2,
            streamed_pairs=json.loads(report.read_text(encoding="utf-8")),
            expected_pairs=expected_pairs,
        )


@given(bundles=lane_chunk_bundles)
@settings(
    max_examples=150,
    deadline=None,
)
def test_illumina_lane_grouping_resolution_and_streaming_preserve_pair_order(
    bundles: list[tuple[int, int, int]],
) -> None:
    case = run_grouped_illumina_lane_case(bundles)
    assert case.warnings == []
    assert case.resolution_warnings == ()
    assert len(case.samples) == 1
    resolved = case.samples[0]
    assert resolved.sample_id == "Sample"
    assert resolved.platform == "illumina"
    assert resolved.source == "paired_globs"
    assert resolved.r1 == case.expected_r1
    assert resolved.r2 == case.expected_r2
    assert case.streamed_pairs == case.expected_pairs


@pytest.mark.slow
@given(bundles=stress_lane_chunk_bundles)
@settings(
    max_examples=500,
    deadline=None,
)
def test_illumina_lane_grouping_resolution_and_streaming_preserve_pair_order_stress(
    bundles: list[tuple[int, int, int]],
) -> None:
    case = run_grouped_illumina_lane_case(bundles)
    assert case.warnings == []
    assert case.resolution_warnings == ()
    assert len(case.samples) == 1
    resolved = case.samples[0]
    assert resolved.sample_id == "Sample"
    assert resolved.platform == "illumina"
    assert resolved.source == "paired_globs"
    assert resolved.r1 == case.expected_r1
    assert resolved.r2 == case.expected_r2
    assert case.streamed_pairs == case.expected_pairs
