"""Tests for the Deacon FIFO streaming helper."""

from __future__ import annotations

import gzip
import json
import lzma
import shutil
import stat
import subprocess
import time
from collections.abc import Callable
from pathlib import Path
from typing import IO, Any

import pytest
from stream_fastqs_to_deacon import (
    DeaconStreamConfig,
    StreamError,
    config_from_args,
    parse_args,
    run_deacon_stream,
)

HANG_GUARD_SECONDS = 5
PAIRED_INPUT_COUNT = 2
ROOT = Path(__file__).resolve().parents[1]
INTEGRATION_DATA = ROOT / "tests" / "data"
Runner = Callable[[list[str]], subprocess.CompletedProcess[str]]

def write_fastq(path: Path, read_names: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(f"@{name}\nACGT\n+\n!!!!\n" for name in read_names)
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    elif path.name.endswith(".xz"):
        with lzma.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def write_path_list(path: Path, files: tuple[Path, ...]) -> Path:
    path.write_text("\n".join(str(file) for file in files) + "\n", encoding="utf-8")
    return path


def command_io(command: list[str]) -> tuple[Path, Path, Path, tuple[Path, ...]]:
    output = Path(command[command.index("--output") + 1])
    summary = Path(command[command.index("--summary") + 1])
    index = command.index("filter") + 1
    while command[index].startswith("--"):
        index += 2
    return output, summary, Path(command[index]), tuple(Path(item) for item in command[index + 1 :])


def open_maybe_gzip(path: Path, mode: str) -> IO[Any]:
    if path.name.endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def read_fastq_names(path: Path) -> list[str]:
    names = []
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


def canonical_pair_name(name: str) -> str:
    return name.removesuffix("/1").removesuffix("/2")


def copying_deacon_runner() -> Runner:
    def run(command: list[str]) -> subprocess.CompletedProcess[str]:
        try:
            output, summary, index, inputs = command_io(command)
            output.parent.mkdir(parents=True, exist_ok=True)
            with output.open("wb") as out:
                for input_path in inputs:
                    with open_maybe_gzip(input_path, "rb") as inp:
                        out.write(inp.read())
            summary.write_text(
                json.dumps({"index": str(index), "inputs": [str(item) for item in inputs]}),
                encoding="utf-8",
            )
            return subprocess.CompletedProcess(command, 0)
        except Exception as error:  # noqa: BLE001
            return subprocess.CompletedProcess(command, 1, stderr=str(error))

    return run


def stat_recording_deacon_runner(report: Path) -> Runner:
    def run(command: list[str]) -> subprocess.CompletedProcess[str]:
        output, summary, _index, inputs = command_io(command)
        records = []
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("wb") as out:
            for input_path in inputs:
                mode = input_path.stat().st_mode
                records.append(
                    {
                        "path": str(input_path),
                        "fifo": stat.S_ISFIFO(mode),
                        "regular": stat.S_ISREG(mode),
                    },
                )
                with input_path.open("rb") as inp:
                    out.write(inp.read())
        report.write_text(json.dumps(records), encoding="utf-8")
        summary.write_text(json.dumps({"inputs": [str(item) for item in inputs]}), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    return run


def failing_deacon_runner() -> Runner:
    return lambda command: subprocess.CompletedProcess(command, 7)


def first_fifo_only_success_deacon_runner() -> Runner:
    def run(command: list[str]) -> subprocess.CompletedProcess[str]:
        _output, _summary, _index, inputs = command_io(command)
        with inputs[0].open("rb") as handle:
            handle.read()
        return subprocess.CompletedProcess(command, 0)

    return run


def partial_read_failing_deacon_runner() -> Runner:
    def run(command: list[str]) -> subprocess.CompletedProcess[str]:
        _output, _summary, _index, inputs = command_io(command)
        for input_path in inputs:
            with input_path.open("rb") as handle:
                handle.read(16)
        return subprocess.CompletedProcess(command, 7)

    return run


def pair_checking_deacon_runner(report: Path) -> Runner:
    def run(command: list[str]) -> subprocess.CompletedProcess[str]:
        output, summary, _index, inputs = command_io(command)
        if len(inputs) != PAIRED_INPUT_COUNT:
            return subprocess.CompletedProcess(command, 1)
        r1_records = read_fastq_names(inputs[0])
        r2_records = read_fastq_names(inputs[1])
        if len(r1_records) != len(r2_records):
            return subprocess.CompletedProcess(command, 1)
        pairs = []
        for r1, r2 in zip(r1_records, r2_records, strict=True):
            left = canonical_pair_name(r1)
            right = canonical_pair_name(r2)
            if left != right:
                return subprocess.CompletedProcess(command, 1)
            pairs.append([left, right])
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("", encoding="utf-8")
        summary.write_text(json.dumps({"pairs": pairs}), encoding="utf-8")
        report.write_text(json.dumps(pairs), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    return run


def config(
    tmp_path: Path,
    deacon_bin: str | Path = "deacon",
    *,
    reads: tuple[Path, ...] = (),
    r1: tuple[Path, ...] = (),
    r2: tuple[Path, ...] = (),
) -> DeaconStreamConfig:
    index = tmp_path / "index.dcn"
    index.write_text("fake index", encoding="utf-8")
    return DeaconStreamConfig(
        sample_id="S1",
        index=index,
        output=tmp_path / "out" / "S1.fastq",
        summary=tmp_path / "out" / "S1.json",
        reads=reads,
        r1=r1,
        r2=r2,
        threads=1,
        abs_threshold=1,
        rel_threshold=0.0,
        deacon_bin=str(deacon_bin),
    )


def test_single_end_gzip_bundle_streams_in_order(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A", "L1_B"])
    lane2 = write_fastq(tmp_path / "lane2.fastq.gz", ["L2_C"])

    run_deacon_stream(config(tmp_path, reads=(lane1, lane2)), runner=copying_deacon_runner())

    output = (tmp_path / "out" / "S1.fastq").read_text(encoding="utf-8")
    assert output.index("@L1_A") < output.index("@L1_B") < output.index("@L2_C")


def test_plain_fastq_bundle_streams_in_order(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq", ["L1_A"])
    lane2 = write_fastq(tmp_path / "lane2.fastq", ["L2_B"])

    run_deacon_stream(config(tmp_path, reads=(lane1, lane2)), runner=copying_deacon_runner())

    output = (tmp_path / "out" / "S1.fastq").read_text(encoding="utf-8")
    assert output.index("@L1_A") < output.index("@L2_B")


def test_single_end_gzip_inputs_are_fastq_gz_fifos_not_regular_files(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"])
    lane2 = write_fastq(tmp_path / "lane2.fastq.gz", ["L2_B"])
    report = tmp_path / "input_stats.json"

    run_deacon_stream(
        config(tmp_path, reads=(lane1, lane2)),
        runner=stat_recording_deacon_runner(report),
    )

    records = json.loads(report.read_text(encoding="utf-8"))
    assert records == [{"path": records[0]["path"], "fifo": True, "regular": False}]
    assert Path(records[0]["path"]).name == "reads.fastq.gz"
    assert not Path(records[0]["path"]).exists()


def test_paired_gzip_inputs_are_fastq_gz_fifos_not_regular_files(tmp_path: Path) -> None:
    r1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r2 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.gz", ["L1_A/2"])
    report = tmp_path / "input_stats.json"

    run_deacon_stream(config(tmp_path, r1=(r1,), r2=(r2,)), runner=stat_recording_deacon_runner(report))

    records = json.loads(report.read_text(encoding="utf-8"))
    assert [(Path(record["path"]).name, record["fifo"], record["regular"]) for record in records] == [
        ("r1.fastq.gz", True, False),
        ("r2.fastq.gz", True, False),
    ]
    assert not Path(records[0]["path"]).exists()
    assert not Path(records[1]["path"]).exists()


def test_paired_bundle_streams_both_mates_in_order(tmp_path: Path) -> None:
    r1_lane1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r1_lane2 = write_fastq(tmp_path / "Sample_L002_R1_001.fastq.gz", ["L2_B/1"])
    r2_lane1 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.gz", ["L1_A/2"])
    r2_lane2 = write_fastq(tmp_path / "Sample_L002_R2_001.fastq.gz", ["L2_B/2"])

    run_deacon_stream(
        config(
            tmp_path,
            r1=(r1_lane1, r1_lane2),
            r2=(r2_lane1, r2_lane2),
        ),
        runner=copying_deacon_runner(),
    )

    output = (tmp_path / "out" / "S1.fastq").read_text(encoding="utf-8")
    assert output.index("@L1_A/1") < output.index("@L2_B/1")
    assert output.index("@L1_A/2") < output.index("@L2_B/2")


def test_paired_gzip_fifo_preserves_record_ordinality_across_lanes(tmp_path: Path) -> None:
    r1_lane1 = write_fastq(
        tmp_path / "Sample_L001_R1_001.fastq.gz",
        ["L001_A/1", "L001_B/1"],
    )
    r2_lane1 = write_fastq(
        tmp_path / "Sample_L001_R2_001.fastq.gz",
        ["L001_A/2", "L001_B/2"],
    )
    r1_lane2 = write_fastq(tmp_path / "Sample_L002_R1_001.fastq.gz", ["L002_C/1"])
    r2_lane2 = write_fastq(tmp_path / "Sample_L002_R2_001.fastq.gz", ["L002_C/2"])
    report = tmp_path / "pairs.json"

    run_deacon_stream(
        config(
            tmp_path,
            r1=(r1_lane1, r1_lane2),
            r2=(r2_lane1, r2_lane2),
        ),
        runner=pair_checking_deacon_runner(report),
    )

    assert json.loads(report.read_text(encoding="utf-8")) == [
        ["L001_A", "L001_A"],
        ["L001_B", "L001_B"],
        ["L002_C", "L002_C"],
    ]


def test_paired_gzip_fifo_detects_swapped_r2_lane_order(tmp_path: Path) -> None:
    r1_lane1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L001_A/1"])
    r2_lane1 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.gz", ["L001_A/2"])
    r1_lane2 = write_fastq(tmp_path / "Sample_L002_R1_001.fastq.gz", ["L002_C/1"])
    r2_lane2 = write_fastq(tmp_path / "Sample_L002_R2_001.fastq.gz", ["L002_C/2"])

    with pytest.raises(StreamError, match="deacon filter failed"):
        run_deacon_stream(
            config(
                tmp_path,
                r1=(r1_lane1, r1_lane2),
                r2=(r2_lane2, r2_lane1),
            ),
            runner=pair_checking_deacon_runner(tmp_path / "pairs.json"),
        )


def test_missing_producer_file_fails_even_if_consumer_would_succeed(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"])
    missing = tmp_path / "missing.fastq.gz"

    with pytest.raises(StreamError, match="missing files"):
        run_deacon_stream(config(tmp_path, reads=(lane1, missing)), runner=copying_deacon_runner())


def test_corrupt_gzip_consumer_failure_is_reported(tmp_path: Path) -> None:
    corrupt = tmp_path / "corrupt.fastq.gz"
    corrupt.write_bytes(b"not actually gzip")

    with pytest.raises(StreamError, match="deacon filter failed"):
        run_deacon_stream(config(tmp_path, reads=(corrupt,)), runner=copying_deacon_runner())


def test_truncated_gzip_consumer_failure_after_partial_stream_is_reported(tmp_path: Path) -> None:
    valid = write_fastq(tmp_path / "valid.fastq.gz", ["L1_A", "L1_B"])
    truncated = tmp_path / "truncated.fastq.gz"
    original = valid.read_bytes()
    truncated.write_bytes(original[:-8])

    with pytest.raises(StreamError, match="deacon filter failed"):
        run_deacon_stream(config(tmp_path, reads=(truncated,)), runner=copying_deacon_runner())


def test_deacon_failure_before_fifo_open_does_not_hang(tmp_path: Path) -> None:
    reads = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"])

    with pytest.raises(StreamError, match="deacon filter failed with exit code 7"):
        run_deacon_stream(config(tmp_path, reads=(reads,)), runner=failing_deacon_runner())


def test_deacon_failure_after_reading_from_fifo_does_not_hang(tmp_path: Path) -> None:
    reads = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A", "L1_B"])

    started = time.monotonic()
    with pytest.raises(StreamError, match="deacon filter failed with exit code 7"):
        run_deacon_stream(config(tmp_path, reads=(reads,)), runner=partial_read_failing_deacon_runner())

    assert time.monotonic() - started < HANG_GUARD_SECONDS


def test_deacon_success_before_opening_all_plain_paired_fifos_does_not_hang(tmp_path: Path) -> None:
    r1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq", ["L1_A/1"])
    r2 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq", ["L1_A/2"])

    started = time.monotonic()
    with pytest.raises(StreamError, match="no reader opened FIFO"):
        run_deacon_stream(
            config(tmp_path, r1=(r1,), r2=(r2,)),
            runner=first_fifo_only_success_deacon_runner(),
        )

    assert time.monotonic() - started < HANG_GUARD_SECONDS


def test_missing_deacon_executable_fails_before_streaming(tmp_path: Path) -> None:
    reads = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"])
    missing_deacon = tmp_path / "definitely-not-deacon"

    with pytest.raises(StreamError, match="deacon executable was not found"):
        run_deacon_stream(config(tmp_path, missing_deacon, reads=(reads,)))


def test_non_executable_deacon_path_fails_before_streaming(tmp_path: Path) -> None:
    reads = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"])
    not_executable = tmp_path / "not_deacon"
    not_executable.write_text("not executable", encoding="utf-8")

    with pytest.raises(StreamError, match="deacon executable was not found"):
        run_deacon_stream(config(tmp_path, not_executable, reads=(reads,)))


def test_mixed_compression_is_rejected(tmp_path: Path) -> None:
    plain = write_fastq(tmp_path / "lane1.fastq", ["L1_A"])
    gz = write_fastq(tmp_path / "lane2.fastq.gz", ["L2_B"])

    with pytest.raises(StreamError, match="one compression type"):
        run_deacon_stream(config(tmp_path, reads=(plain, gz)), runner=copying_deacon_runner())


def test_xz_bundle_streams_in_order(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq.xz", ["L1_A"])
    lane2 = write_fastq(tmp_path / "lane2.fastq.xz", ["L2_B"])

    run_deacon_stream(config(tmp_path, reads=(lane1, lane2)), runner=copying_deacon_runner())

    output = (tmp_path / "out" / "S1.fastq").read_text(encoding="utf-8")
    assert output.index("@L1_A") < output.index("@L2_B")


def test_paired_list_length_mismatch_is_rejected(tmp_path: Path) -> None:
    r1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r2 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.gz", ["L1_A/2"])
    extra_r2 = write_fastq(tmp_path / "Sample_L002_R2_001.fastq.gz", ["L2_B/2"])

    with pytest.raises(StreamError, match="R1/R2 file list length mismatch"):
        run_deacon_stream(config(tmp_path, r1=(r1,), r2=(r2, extra_r2)), runner=copying_deacon_runner())


def test_paired_mixed_compression_across_sample_is_rejected(tmp_path: Path) -> None:
    r1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r2 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.xz", ["L1_A/2"])

    with pytest.raises(StreamError, match="one compression type"):
        run_deacon_stream(config(tmp_path, r1=(r1,), r2=(r2,)), runner=copying_deacon_runner())


def test_zst_bundle_streams_raw_bytes_through_fastq_zst_fifo(tmp_path: Path) -> None:
    zst = tmp_path / "lane1.fastq.zst"
    zst.write_bytes(b"placeholder")
    report = tmp_path / "input_stats.json"

    run_deacon_stream(config(tmp_path, reads=(zst,)), runner=stat_recording_deacon_runner(report))

    records = json.loads(report.read_text(encoding="utf-8"))
    assert [(Path(record["path"]).name, record["fifo"], record["regular"]) for record in records] == [
        ("reads.fastq.zst", True, False),
    ]
    assert (tmp_path / "out" / "S1.fastq").read_bytes() == b"placeholder"


def test_paired_zst_inputs_are_fastq_zst_fifos(tmp_path: Path) -> None:
    r1 = tmp_path / "Sample_L001_R1_001.fastq.zst"
    r2 = tmp_path / "Sample_L001_R2_001.fastq.zst"
    r1.write_bytes(b"r1-zstd-bytes")
    r2.write_bytes(b"r2-zstd-bytes")
    report = tmp_path / "input_stats.json"

    run_deacon_stream(config(tmp_path, r1=(r1,), r2=(r2,)), runner=stat_recording_deacon_runner(report))

    records = json.loads(report.read_text(encoding="utf-8"))
    assert [(Path(record["path"]).name, record["fifo"], record["regular"]) for record in records] == [
        ("r1.fastq.zst", True, False),
        ("r2.fastq.zst", True, False),
    ]


def test_cli_list_files_are_supported(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    reads = (write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"]),)
    reads_list = write_path_list(tmp_path / "reads.txt", reads)
    argv = [
        "--sample-id",
        "S1",
        "--index",
        str(tmp_path / "index.dcn"),
        "--output",
        str(tmp_path / "out.fastq.gz"),
        "--summary",
        str(tmp_path / "summary.json"),
        "--reads-list",
        str(reads_list),
    ]
    (tmp_path / "index.dcn").write_text("fake", encoding="utf-8")

    parsed = parse_args(argv)
    loaded = config_from_args(parsed)

    assert loaded.reads == reads
    assert capsys.readouterr().err == ""


@pytest.fixture(scope="module")
def deacon_bin() -> str:
    deacon = shutil.which("deacon")
    if deacon is None:
        pytest.skip("deacon executable is not available")
    return deacon


def build_deacon_index(tmp_path: Path, deacon: str) -> Path:
    reference = tmp_path / "reference.fasta"
    reference.write_text(
        ">toy\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
        encoding="utf-8",
    )
    index = tmp_path / "toy.k3w1.idx"
    subprocess.run(  # noqa: S603
        [deacon, "index", "build", "-k", "3", "-w", "1", "-q", "-o", str(index), str(reference)],
        check=True,
    )
    return index


def concatenate_gzip_to_plain(output: Path, files: tuple[Path, ...]) -> Path:
    with output.open("w", encoding="utf-8") as out:
        for path in files:
            with gzip.open(path, "rt", encoding="utf-8") as inp:
                out.write(inp.read())
    return output


def run_direct_deacon(
    tmp_path: Path,
    deacon: str,
    index: Path,
    inputs: tuple[Path, ...],
) -> Path:
    output = tmp_path / "control" / "out.fastq"
    summary = tmp_path / "control" / "summary.json"
    output.parent.mkdir(parents=True)
    subprocess.run(  # noqa: S603
        [
            deacon,
            "filter",
            "--threads",
            "1",
            "--abs-threshold",
            "1",
            "--rel-threshold",
            "0.0",
            "--summary",
            str(summary),
            "--output",
            str(output),
            str(index),
            *(str(path) for path in inputs),
        ],
        check=True,
    )
    return output


def real_deacon_config(  # noqa: PLR0913
    tmp_path: Path,
    deacon: str,
    index: Path,
    *,
    reads: tuple[Path, ...] = (),
    r1: tuple[Path, ...] = (),
    r2: tuple[Path, ...] = (),
) -> DeaconStreamConfig:
    return DeaconStreamConfig(
        sample_id="real",
        index=index,
        output=tmp_path / "real" / "out.fastq",
        summary=tmp_path / "real" / "summary.json",
        reads=reads,
        r1=r1,
        r2=r2,
        threads=1,
        abs_threshold=1,
        rel_threshold=0.0,
        deacon_bin=deacon,
    )


def test_real_deacon_reads_single_end_fifo_stream(tmp_path: Path, deacon_bin: str) -> None:
    index = build_deacon_index(tmp_path, deacon_bin)
    lane1 = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A", "L1_B"])
    lane2 = write_fastq(tmp_path / "lane2.fastq.gz", ["L2_C"])
    control_input = concatenate_gzip_to_plain(tmp_path / "control.fastq", (lane1, lane2))
    control_output = run_direct_deacon(tmp_path, deacon_bin, index, (control_input,))

    run_deacon_stream(real_deacon_config(tmp_path, deacon_bin, index, reads=(lane1, lane2)))

    output_path = tmp_path / "real" / "out.fastq"
    assert output_path.read_bytes() == control_output.read_bytes()
    output = output_path.read_text(encoding="utf-8")
    assert output.index("@L1_A") < output.index("@L1_B") < output.index("@L2_C")
    assert (tmp_path / "real" / "summary.json").is_file()


def test_real_deacon_reads_paired_fifo_stream(tmp_path: Path, deacon_bin: str) -> None:
    index = build_deacon_index(tmp_path, deacon_bin)
    r1_lane1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r1_lane2 = write_fastq(tmp_path / "Sample_L002_R1_001.fastq.gz", ["L2_B/1"])
    r2_lane1 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.gz", ["L1_A/2"])
    r2_lane2 = write_fastq(tmp_path / "Sample_L002_R2_001.fastq.gz", ["L2_B/2"])
    control_r1 = concatenate_gzip_to_plain(tmp_path / "control_R1.fastq", (r1_lane1, r1_lane2))
    control_r2 = concatenate_gzip_to_plain(tmp_path / "control_R2.fastq", (r2_lane1, r2_lane2))
    control_output = run_direct_deacon(tmp_path, deacon_bin, index, (control_r1, control_r2))

    run_deacon_stream(
        real_deacon_config(
            tmp_path,
            deacon_bin,
            index,
            r1=(r1_lane1, r1_lane2),
            r2=(r2_lane1, r2_lane2),
        ),
    )

    output_path = tmp_path / "real" / "out.fastq"
    assert output_path.read_bytes() == control_output.read_bytes()
    output = output_path.read_text(encoding="utf-8")
    assert output.index("@L1_A/1") < output.index("@L1_A/2")
    assert output.index("@L2_B/1") < output.index("@L2_B/2")
    assert output.index("@L1_A/2") < output.index("@L2_B/1")


def test_real_deacon_streams_larger_paired_bundle_in_pair_order(
    tmp_path: Path,
    deacon_bin: str,
) -> None:
    index = INTEGRATION_DATA / "mini_virus_deacon.k31w1.idx"
    r1_lane1 = INTEGRATION_DATA / "hits_only_R1.fastq.gz"
    r2_lane1 = INTEGRATION_DATA / "hits_only_R2.fastq.gz"
    r1_lane2 = INTEGRATION_DATA / "water_plus_hits_R1.fastq.gz"
    r2_lane2 = INTEGRATION_DATA / "water_plus_hits_R2.fastq.gz"
    fixtures = (index, r1_lane1, r2_lane1, r1_lane2, r2_lane2)
    missing = [path for path in fixtures if not path.is_file()]
    if missing:
        pytest.skip("integration FASTQ fixtures are unavailable")

    control_r1 = concatenate_gzip_to_plain(tmp_path / "control_R1.fastq", (r1_lane1, r1_lane2))
    control_r2 = concatenate_gzip_to_plain(tmp_path / "control_R2.fastq", (r2_lane1, r2_lane2))
    control_output = run_direct_deacon(tmp_path, deacon_bin, index, (control_r1, control_r2))

    run_deacon_stream(
        real_deacon_config(
            tmp_path,
            deacon_bin,
            index,
            r1=(r1_lane1, r1_lane2),
            r2=(r2_lane1, r2_lane2),
        ),
    )

    output_path = tmp_path / "real" / "out.fastq"
    assert output_path.read_bytes() == control_output.read_bytes()
    assert (tmp_path / "real" / "summary.json").is_file()
