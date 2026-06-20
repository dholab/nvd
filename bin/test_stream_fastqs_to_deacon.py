"""Tests for the Deacon FIFO streaming helper."""

from __future__ import annotations

import gzip
import json
import lzma
import shutil
import subprocess
import time
from pathlib import Path

import pytest
from stream_fastqs_to_deacon import (
    DeaconStreamConfig,
    StreamError,
    config_from_args,
    parse_args,
    run_deacon_stream,
)

HANG_GUARD_SECONDS = 5

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


def fake_deacon(path: Path) -> Path:
    path.write_text(
        """#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

args = sys.argv[1:]
output = Path(args[args.index('--output') + 1])
summary = Path(args[args.index('--summary') + 1])
index = args.index('filter') + 1
while args[index].startswith('--'):
    flag = args[index]
    index += 2
index_path = Path(args[index])
inputs = [Path(item) for item in args[index + 1:]]
output.parent.mkdir(parents=True, exist_ok=True)
with output.open('wb') as out:
    for input_path in inputs:
        with input_path.open('rb') as inp:
            out.write(inp.read())
summary.write_text(json.dumps({'index': str(index_path), 'inputs': [str(item) for item in inputs]}), encoding='utf-8')
""",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def stat_recording_deacon(path: Path, report: Path) -> Path:
    path.write_text(
        f"""#!/usr/bin/env python3
from __future__ import annotations

import json
import stat
import sys
from pathlib import Path

args = sys.argv[1:]
output = Path(args[args.index('--output') + 1])
summary = Path(args[args.index('--summary') + 1])
index = args.index('filter') + 1
while args[index].startswith('--'):
    index += 2
inputs = [Path(item) for item in args[index + 1:]]
records = []
output.parent.mkdir(parents=True, exist_ok=True)
with output.open('wb') as out:
    for input_path in inputs:
        mode = input_path.stat().st_mode
        records.append({{'path': str(input_path), 'fifo': stat.S_ISFIFO(mode), 'regular': stat.S_ISREG(mode)}})
        with input_path.open('rb') as inp:
            out.write(inp.read())
Path({str(report)!r}).write_text(json.dumps(records), encoding='utf-8')
summary.write_text(json.dumps({{'inputs': [str(item) for item in inputs]}}), encoding='utf-8')
""",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def failing_deacon(path: Path) -> Path:
    path.write_text(
        """#!/usr/bin/env python3
from __future__ import annotations

import sys

sys.exit(7)
""",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def marker_deacon(path: Path, marker: Path) -> Path:
    path.write_text(
        f"""#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

Path({str(marker)!r}).write_text('started', encoding='utf-8')
""",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def first_fifo_only_success_deacon(path: Path) -> Path:
    path.write_text(
        """#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

args = sys.argv[1:]
index = args.index('filter') + 1
while args[index].startswith('--'):
    index += 2
inputs = [Path(item) for item in args[index + 1:]]
with inputs[0].open('rb') as handle:
    handle.read()
sys.exit(0)
""",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def partial_read_failing_deacon(path: Path) -> Path:
    path.write_text(
        """#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

args = sys.argv[1:]
index = args.index('filter') + 1
while args[index].startswith('--'):
    index += 2
inputs = [Path(item) for item in args[index + 1:]]
for input_path in inputs:
    with input_path.open('rb') as handle:
        handle.read(16)
sys.exit(7)
""",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def config(
    tmp_path: Path,
    deacon_bin: Path,
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
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    run_deacon_stream(config(tmp_path, deacon, reads=(lane1, lane2)))

    output = (tmp_path / "out" / "S1.fastq").read_text(encoding="utf-8")
    assert output.index("@L1_A") < output.index("@L1_B") < output.index("@L2_C")


def test_plain_fastq_bundle_streams_in_order(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq", ["L1_A"])
    lane2 = write_fastq(tmp_path / "lane2.fastq", ["L2_B"])
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    run_deacon_stream(config(tmp_path, deacon, reads=(lane1, lane2)))

    output = (tmp_path / "out" / "S1.fastq").read_text(encoding="utf-8")
    assert output.index("@L1_A") < output.index("@L2_B")


def test_deacon_inputs_are_fifos_not_regular_files(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"])
    lane2 = write_fastq(tmp_path / "lane2.fastq.gz", ["L2_B"])
    report = tmp_path / "input_stats.json"
    deacon = stat_recording_deacon(tmp_path / "stat_deacon.py", report)

    run_deacon_stream(config(tmp_path, deacon, reads=(lane1, lane2)))

    records = json.loads(report.read_text(encoding="utf-8"))
    assert records == [{"path": records[0]["path"], "fifo": True, "regular": False}]
    assert not Path(records[0]["path"]).exists()


def test_paired_bundle_streams_both_mates_in_order(tmp_path: Path) -> None:
    r1_lane1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r1_lane2 = write_fastq(tmp_path / "Sample_L002_R1_001.fastq.gz", ["L2_B/1"])
    r2_lane1 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.gz", ["L1_A/2"])
    r2_lane2 = write_fastq(tmp_path / "Sample_L002_R2_001.fastq.gz", ["L2_B/2"])
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    run_deacon_stream(
        config(
            tmp_path,
            deacon,
            r1=(r1_lane1, r1_lane2),
            r2=(r2_lane1, r2_lane2),
        ),
    )

    output = (tmp_path / "out" / "S1.fastq").read_text(encoding="utf-8")
    assert output.index("@L1_A/1") < output.index("@L2_B/1")
    assert output.index("@L1_A/2") < output.index("@L2_B/2")


def test_missing_producer_file_fails_even_if_consumer_would_succeed(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"])
    missing = tmp_path / "missing.fastq.gz"
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    with pytest.raises(StreamError, match="missing files"):
        run_deacon_stream(config(tmp_path, deacon, reads=(lane1, missing)))


def test_corrupt_gzip_producer_failure_is_reported(tmp_path: Path) -> None:
    corrupt = tmp_path / "corrupt.fastq.gz"
    corrupt.write_bytes(b"not actually gzip")
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    with pytest.raises(StreamError, match="FASTQ stream producer failed"):
        run_deacon_stream(config(tmp_path, deacon, reads=(corrupt,)))


def test_truncated_gzip_failure_after_partial_stream_is_reported(tmp_path: Path) -> None:
    valid = write_fastq(tmp_path / "valid.fastq.gz", ["L1_A", "L1_B"])
    truncated = tmp_path / "truncated.fastq.gz"
    original = valid.read_bytes()
    truncated.write_bytes(original[:-8])
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    with pytest.raises(StreamError, match="FASTQ stream producer failed"):
        run_deacon_stream(config(tmp_path, deacon, reads=(truncated,)))


def test_deacon_failure_before_fifo_open_does_not_hang(tmp_path: Path) -> None:
    reads = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A"])
    deacon = failing_deacon(tmp_path / "failing_deacon.py")

    with pytest.raises(StreamError, match="deacon filter failed with exit code 7"):
        run_deacon_stream(config(tmp_path, deacon, reads=(reads,)))


def test_deacon_failure_after_reading_from_fifo_does_not_hang(tmp_path: Path) -> None:
    reads = write_fastq(tmp_path / "lane1.fastq.gz", ["L1_A", "L1_B"])
    deacon = partial_read_failing_deacon(tmp_path / "partial_read_failing.py")

    started = time.monotonic()
    with pytest.raises(StreamError, match="deacon filter failed with exit code 7"):
        run_deacon_stream(config(tmp_path, deacon, reads=(reads,)))

    assert time.monotonic() - started < HANG_GUARD_SECONDS


def test_deacon_success_before_opening_all_paired_fifos_does_not_hang(tmp_path: Path) -> None:
    r1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r2 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.gz", ["L1_A/2"])
    deacon = first_fifo_only_success_deacon(tmp_path / "first_fifo_only.py")

    started = time.monotonic()
    with pytest.raises(StreamError, match="no reader opened FIFO"):
        run_deacon_stream(config(tmp_path, deacon, r1=(r1,), r2=(r2,)))

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
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    with pytest.raises(StreamError, match="one compression type"):
        run_deacon_stream(config(tmp_path, deacon, reads=(plain, gz)))


def test_xz_bundle_streams_in_order(tmp_path: Path) -> None:
    lane1 = write_fastq(tmp_path / "lane1.fastq.xz", ["L1_A"])
    lane2 = write_fastq(tmp_path / "lane2.fastq.xz", ["L2_B"])
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    run_deacon_stream(config(tmp_path, deacon, reads=(lane1, lane2)))

    output = (tmp_path / "out" / "S1.fastq").read_text(encoding="utf-8")
    assert output.index("@L1_A") < output.index("@L2_B")


def test_paired_list_length_mismatch_is_rejected(tmp_path: Path) -> None:
    r1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r2 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.gz", ["L1_A/2"])
    extra_r2 = write_fastq(tmp_path / "Sample_L002_R2_001.fastq.gz", ["L2_B/2"])
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    with pytest.raises(StreamError, match="R1/R2 file list length mismatch"):
        run_deacon_stream(config(tmp_path, deacon, r1=(r1,), r2=(r2, extra_r2)))


def test_paired_mixed_compression_across_sample_is_rejected(tmp_path: Path) -> None:
    r1 = write_fastq(tmp_path / "Sample_L001_R1_001.fastq.gz", ["L1_A/1"])
    r2 = write_fastq(tmp_path / "Sample_L001_R2_001.fastq.xz", ["L1_A/2"])
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    with pytest.raises(StreamError, match="one compression type"):
        run_deacon_stream(config(tmp_path, deacon, r1=(r1,), r2=(r2,)))


def test_zst_input_is_rejected_until_supported(tmp_path: Path) -> None:
    zst = tmp_path / "lane1.fastq.zst"
    zst.write_bytes(b"placeholder")
    deacon = fake_deacon(tmp_path / "fake_deacon.py")

    with pytest.raises(StreamError, match="zst input streaming is not implemented"):
        run_deacon_stream(config(tmp_path, deacon, reads=(zst,)))


def test_zst_input_is_rejected_before_deacon_starts(tmp_path: Path) -> None:
    zst = tmp_path / "lane1.fastq.zst"
    zst.write_bytes(b"placeholder")
    marker = tmp_path / "started.txt"
    deacon = marker_deacon(tmp_path / "marker_deacon.py", marker)

    with pytest.raises(StreamError, match="zst input streaming is not implemented"):
        run_deacon_stream(config(tmp_path, deacon, reads=(zst,)))

    assert not marker.exists()


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
