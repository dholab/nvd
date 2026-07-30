"""Behavioral tests for streamed SRA source selection and Deacon execution."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import stat
import threading
from collections.abc import Iterator  # noqa: TC003 - keep imports at module frontmatter
from contextlib import contextmanager
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path  # noqa: TC003 - keep imports at module frontmatter

import pytest
from stream_sra_to_deacon import (
    BackendProcessError,
    IntegrityError,
    StreamConfig,
    UnsupportedLayoutError,
    configure_logging,
    execute,
    run,
)

HTTP_OK = 200
REQUESTS_WITH_ONE_RETRY = 2


@dataclass
class EnaTestServer:
    file_report_url: str
    object_host: str
    routes: dict[str, bytes]
    request_counts: dict[str, int]

    def set_report(self, payload: object) -> None:
        self.routes["/filereport"] = json.dumps(payload).encode()


@contextmanager
def ena_server(
    payload: object,
    objects: dict[str, bytes] | None = None,
    report_statuses: list[int] | None = None,
    omit_content_length: set[str] | None = None,
) -> Iterator[EnaTestServer]:
    routes = {"/filereport": json.dumps(payload).encode()}
    routes.update(objects or {})
    statuses = list(report_statuses or [200])
    paths_without_length = omit_content_length or set()
    request_counts: dict[str, int] = {}

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            path, _separator, _query = self.path.partition("?")
            request_counts[path] = request_counts.get(path, 0) + 1
            body = routes.get(path)
            if body is None:
                self.send_error(404)
                return
            status = (
                statuses.pop(0)
                if path == "/filereport" and len(statuses) > 1
                else statuses[0]
            )
            if status != HTTP_OK:
                self.send_error(status)
                return
            self.send_response(HTTP_OK)
            if path not in paths_without_length:
                self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            try:
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                return

        def log_message(self, _format: str, *_args: object) -> None:
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever)
    thread.start()
    try:
        host, port = server.server_address
        object_host = f"{host}:{port}"
        yield EnaTestServer(
            file_report_url=f"http://{object_host}/filereport",
            object_host=object_host,
            routes=routes,
            request_counts=request_counts,
        )
    finally:
        server.shutdown()
        thread.join()
        server.server_close()


def write_executable(path: Path, source: str) -> Path:
    path.write_text(source, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def raw_fastq_paths(tmp_path: Path) -> list[Path]:
    return [
        path
        for pattern in ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
        for path in tmp_path.rglob(pattern)
        if not path.name.endswith(".target_enriched.fastq.gz")
    ]


def paired_report(
    server: EnaTestServer,
    r1: bytes,
    r2: bytes,
    *,
    r1_md5: str | None = None,
    r1_bytes: int | None = None,
) -> list[dict[str, str]]:
    return [
        {
            "run_accession": "SRR_TEST",
            "instrument_platform": "ILLUMINA",
            "library_layout": "PAIRED",
            "fastq_ftp": (
                f"{server.object_host}/SRR_TEST_1.fastq.gz;"
                f"{server.object_host}/SRR_TEST_2.fastq.gz"
            ),
            "fastq_md5": (
                f"{r1_md5 or hashlib.md5(r1).hexdigest()};"  # noqa: S324 - transfer checksum
                f"{hashlib.md5(r2).hexdigest()}"  # noqa: S324 - transfer checksum
            ),
            "fastq_bytes": f"{r1_bytes or len(r1)};{len(r2)}",
        },
    ]


def write_sra_tool_fakes(bin_dir: Path) -> None:
    write_executable(
        bin_dir / "prefetch",
        """#!/usr/bin/env python3
from pathlib import Path
import sys

accession = sys.argv[-1]
archive = Path(accession)
archive.mkdir(parents=True, exist_ok=True)
(archive / f"{accession}.sra").write_bytes(b"archive")
""",
    )
    write_executable(
        bin_dir / "vdb-validate",
        """#!/usr/bin/env python3
raise SystemExit(0)
""",
    )
    write_executable(
        bin_dir / "fastq-dump",
        """#!/usr/bin/env python3
import sys

sys.stdout.write(
    "@SRR_TEST.1.1 source length=4\\nACG1\\n+\\nIIII\\n"
    "@SRR_TEST.1.2 source length=4\\nTGC2\\n+\\nIIII\\n"
)
""",
    )
    write_executable(
        bin_dir / "deacon",
        """#!/usr/bin/env python3
import gzip
import json
import os
import stat
import sys
from pathlib import Path

args = sys.argv[1:]
if os.environ.get("FAKE_DEACON_EXIT_EARLY"):
    raise SystemExit(0)
output = Path(args[args.index("--output") + 1])
summary = Path(args[args.index("--summary") + 1])
inputs = args[-2:]
report_path = os.environ.get("FAKE_DEACON_INPUT_REPORT")
if report_path and inputs != ["-", "-"]:
    records = [
        {"path": path, "fifo": stat.S_ISFIFO(Path(path).stat().st_mode)}
        for path in inputs
    ]
    Path(report_path).write_text(json.dumps(records), encoding="utf-8")
if inputs == ["-", "-"]:
    payload = sys.stdin.buffer.read()
else:
    payload = b"".join(gzip.decompress(Path(path).read_bytes()) for path in inputs)
output.write_bytes(payload)
summary.write_text(json.dumps({"seqs_in": 2, "seqs_out": 2}), encoding="utf-8")
""",
    )


def stream_config(
    tmp_path: Path,
    file_report_url: str,
    *,
    ena_object_host: str = "ftp.sra.ebi.ac.uk",
    ena_object_scheme: str = "https",
) -> StreamConfig:
    index = tmp_path / "target.idx"
    index.write_bytes(b"index")
    return StreamConfig(
        accession="SRR_TEST",
        sra_directory=tmp_path / "sra",
        index=index,
        output=tmp_path / "sample.target_enriched.fastq.gz",
        summary=tmp_path / "sample.deacon_filter.json",
        threads=2,
        abs_threshold=1,
        rel_threshold=0.0,
        deplete=False,
        ena_file_report_url=file_report_url,
        ena_object_host=ena_object_host,
        ena_object_scheme=ena_object_scheme,
        retry_delay_seconds=0.0,
    )


def test_run_selects_sra_when_ena_has_no_compatible_fastqs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_sra_tool_fakes(bin_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    report = [
        {
            "run_accession": "SRR_TEST",
            "instrument_platform": "ILLUMINA",
            "library_layout": "PAIRED",
            "fastq_ftp": "",
            "fastq_md5": "",
            "fastq_bytes": "",
        },
    ]

    configure_logging()
    with ena_server(report) as server:
        result = execute(stream_config(tmp_path, server.file_report_url))

    event = capsys.readouterr().err
    assert result.backend == "fastq-dump"
    assert result.selection_reason == "ena_fastq_absent"
    assert event.count("nvd.sra_stream ") == 1
    assert "backend=fastq-dump" in event
    assert "selection_reason=ena_fastq_absent" in event
    assert "outcome=success" in event
    assert (
        (tmp_path / "sample.target_enriched.fastq.gz")
        .read_text(
            encoding="utf-8",
        )
        .startswith("@SRR_TEST.1.1")
    )
    assert json.loads(
        (tmp_path / "sample.deacon_filter.json").read_text(encoding="utf-8"),
    ) == {"seqs_in": 2, "seqs_out": 2}
    assert raw_fastq_paths(tmp_path) == []


def test_run_streams_compatible_ena_pair_through_compressed_fifos(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_sra_tool_fakes(bin_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    input_report = tmp_path / "deacon-inputs.json"
    monkeypatch.setenv("FAKE_DEACON_INPUT_REPORT", str(input_report))
    r1 = gzip.compress(b"@SRR_TEST.1 source/1\nACGT\n+\nIIII\n")
    r2 = gzip.compress(b"@SRR_TEST.1 source/2\nTGCA\n+\nIIII\n")

    object_routes = {
        "/SRR_TEST_1.fastq.gz": r1,
        "/SRR_TEST_2.fastq.gz": r2,
    }
    with ena_server({}, object_routes) as server:
        server.set_report(paired_report(server, r1, r2))
        configure_logging()
        result = execute(
            stream_config(
                tmp_path,
                server.file_report_url,
                ena_object_host=server.object_host,
                ena_object_scheme="http",
            ),
        )

    output = (tmp_path / "sample.target_enriched.fastq.gz").read_text(
        encoding="utf-8",
    )
    event = capsys.readouterr().err
    deacon_inputs = json.loads(input_report.read_text(encoding="utf-8"))
    assert result.backend == "ena"
    assert result.selection_reason == "compatible_pair"
    assert [record["fifo"] for record in deacon_inputs] == [True, True]
    assert event.count("nvd.sra_stream ") == 1
    assert "backend=ena" in event
    assert "outcome=success" in event
    assert "source/1" in output
    assert "source/2" in output
    assert not (tmp_path / "sra").exists()
    assert raw_fastq_paths(tmp_path) == []


def test_committed_ena_integrity_failure_does_not_fall_back(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_sra_tool_fakes(bin_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    r1 = gzip.compress(b"@SRR_TEST.1 source/1\nACGT\n+\nIIII\n")
    r2 = gzip.compress(b"@SRR_TEST.1 source/2\nTGCA\n+\nIIII\n")
    object_routes = {
        "/SRR_TEST_1.fastq.gz": r1,
        "/SRR_TEST_2.fastq.gz": r2,
    }

    configure_logging()
    with ena_server({}, object_routes) as server:
        server.set_report(paired_report(server, r1, r2, r1_md5="0" * 32))
        with pytest.raises(IntegrityError, match="r1 compressed MD5 mismatch"):
            execute(
                stream_config(
                    tmp_path,
                    server.file_report_url,
                    ena_object_host=server.object_host,
                    ena_object_scheme="http",
                ),
            )

    event = capsys.readouterr().err
    assert not (tmp_path / "sra").exists()
    assert event.count("nvd.sra_stream ") == 1
    assert "backend=ena" in event
    assert "outcome=error" in event
    assert "error_type=IntegrityError" in event


def test_committed_ena_stops_when_deacon_exits_before_opening_fifos(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_sra_tool_fakes(bin_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_DEACON_EXIT_EARLY", "1")
    r1 = gzip.compress(b"@SRR_TEST.1 source/1\nACGT\n+\nIIII\n")
    r2 = gzip.compress(b"@SRR_TEST.1 source/2\nTGCA\n+\nIIII\n")
    objects = {
        "/SRR_TEST_1.fastq.gz": r1,
        "/SRR_TEST_2.fastq.gz": r2,
    }

    with ena_server({}, objects) as server:
        server.set_report(paired_report(server, r1, r2))
        config = stream_config(
            tmp_path,
            server.file_report_url,
            ena_object_host=server.object_host,
            ena_object_scheme="http",
        )
        with pytest.raises(BackendProcessError, match="before both ENA objects"):
            run(config)

    assert not (tmp_path / "sra").exists()


def test_run_selects_sra_when_ena_lookup_remains_unavailable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_sra_tool_fakes(bin_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")

    with ena_server({}, report_statuses=[503, 503]) as server:
        result = run(stream_config(tmp_path, server.file_report_url))

    assert result.backend == "fastq-dump"
    assert result.selection_reason == "ena_lookup_unavailable"
    assert server.request_counts["/filereport"] == REQUESTS_WITH_ONE_RETRY


def test_run_rejects_single_end_layout_from_successful_ena_report(
    tmp_path: Path,
) -> None:
    report = [
        {
            "run_accession": "SRR_TEST",
            "instrument_platform": "ILLUMINA",
            "library_layout": "SINGLE",
            "fastq_ftp": "",
            "fastq_md5": "",
            "fastq_bytes": "",
        },
    ]

    with (
        ena_server(report) as server,
        pytest.raises(UnsupportedLayoutError, match="requires a paired run"),
    ):
        run(stream_config(tmp_path, server.file_report_url))

    assert not (tmp_path / "sra").exists()


def test_run_selects_sra_when_ena_object_fails_preflight(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_sra_tool_fakes(bin_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    r1 = gzip.compress(b"@SRR_TEST.1 source/1\nACGT\n+\nIIII\n")
    r2 = gzip.compress(b"@SRR_TEST.1 source/2\nTGCA\n+\nIIII\n")
    object_routes = {
        "/SRR_TEST_1.fastq.gz": r1,
        "/SRR_TEST_2.fastq.gz": r2,
    }

    with ena_server({}, object_routes) as server:
        server.set_report(paired_report(server, r1, r2, r1_bytes=len(r1) + 1))
        result = run(
            stream_config(
                tmp_path,
                server.file_report_url,
                ena_object_host=server.object_host,
                ena_object_scheme="http",
            ),
        )

    assert result.backend == "fastq-dump"
    assert result.selection_reason == "ena_objects_unavailable"
    assert server.request_counts["/SRR_TEST_1.fastq.gz"] == REQUESTS_WITH_ONE_RETRY
    assert server.request_counts["/SRR_TEST_2.fastq.gz"] == REQUESTS_WITH_ONE_RETRY
    assert (
        (tmp_path / "sample.target_enriched.fastq.gz")
        .read_text(encoding="utf-8")
        .startswith("@SRR_TEST.1.1")
    )


def test_run_selects_sra_when_ena_object_has_no_content_length(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_sra_tool_fakes(bin_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    r1 = gzip.compress(b"@SRR_TEST.1 source/1\nACGT\n+\nIIII\n")
    r2 = gzip.compress(b"@SRR_TEST.1 source/2\nTGCA\n+\nIIII\n")
    objects = {
        "/SRR_TEST_1.fastq.gz": r1,
        "/SRR_TEST_2.fastq.gz": r2,
    }

    with ena_server(
        {},
        objects,
        omit_content_length={"/SRR_TEST_1.fastq.gz"},
    ) as server:
        server.set_report(paired_report(server, r1, r2))
        result = run(
            stream_config(
                tmp_path,
                server.file_report_url,
                ena_object_host=server.object_host,
                ena_object_scheme="http",
            ),
        )

    assert result.backend == "fastq-dump"
    assert result.selection_reason == "ena_objects_unavailable"
    assert server.request_counts["/SRR_TEST_1.fastq.gz"] == REQUESTS_WITH_ONE_RETRY
