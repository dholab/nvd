#!/usr/bin/env python3
"""Prefer paired ENA FASTQs while streaming one SRA run into Deacon."""

from __future__ import annotations

import argparse
import errno
import hashlib
import http.client
import json
import os
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path
from threading import Event
from typing import BinaryIO, Literal, cast

from loguru import logger

ENA_FILE_REPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FILE_REPORT_FIELDS = (
    "run_accession",
    "instrument_platform",
    "library_layout",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
)
HTTP_TOO_MANY_REQUESTS = 429
HTTP_SERVER_ERROR_MIN = 500
COPY_BUFFER_BYTES = 1024 * 1024
ENA_PAIR_FILE_COUNT = 2
HTTP_OK = 200
MD5_DIGEST_BYTES = 16
FIFO_OPEN_TIMEOUT_SECONDS = 30.0
PROCESS_TERMINATION_TIMEOUT_SECONDS = 5.0
SUPERVISOR_POLL_SECONDS = 0.05


class SraStreamError(RuntimeError):
    """Expected failure while streaming one run into Deacon."""

    def __init__(
        self,
        message: str,
        *,
        backend: str,
        selection_reason: str,
        component: str,
    ) -> None:
        super().__init__(message)
        self.backend = backend
        self.selection_reason = selection_reason
        self.component = component


class UnsupportedLayoutError(SraStreamError):
    """Raised when a run is outside the paired-Illumina contract."""

    def __init__(self, message: str) -> None:
        super().__init__(
            message,
            backend="none",
            selection_reason="unsupported_layout",
            component="ena_discovery",
        )


class BackendProcessError(SraStreamError):
    """Raised when an external backend process fails."""


class EnaTransferError(SraStreamError):
    """Raised when an ENA object cannot be streamed completely."""

    def __init__(self, message: str, component: str) -> None:
        super().__init__(
            message,
            backend="ena",
            selection_reason="compatible_pair",
            component=component,
        )


class IntegrityError(SraStreamError):
    """Raised when transferred bytes do not match ENA metadata."""

    def __init__(self, message: str, component: str) -> None:
        super().__init__(
            message,
            backend="ena",
            selection_reason="compatible_pair",
            component=component,
        )


@dataclass(frozen=True)
class StreamConfig:
    accession: str
    sra_directory: Path
    index: Path
    output: Path
    summary: Path
    threads: int
    abs_threshold: int
    rel_threshold: float
    deplete: bool
    ena_file_report_url: str = ENA_FILE_REPORT_URL
    ena_object_host: str = "ftp.sra.ebi.ac.uk"
    ena_object_scheme: str = "https"
    http_timeout_seconds: float = 30.0
    retry_delay_seconds: float = 1.0


@dataclass(frozen=True)
class EnaFastq:
    url: str
    expected_bytes: int
    expected_md5: bytes


@dataclass(frozen=True)
class EnaPair:
    fastqs: tuple[EnaFastq, EnaFastq]


@dataclass(frozen=True)
class Sra:
    reason: str


@dataclass(frozen=True)
class StreamResult:
    backend: Literal["ena", "fastq-dump"]
    selection_reason: str
    transferred_bytes: tuple[int, int] | None = None


def normalize_ena_url(raw_url: str, config: StreamConfig) -> str | None:
    value = raw_url.strip()
    if value.startswith("ftp://"):
        host_path = value.removeprefix("ftp://")
    elif value.startswith("https://"):
        parsed = urllib.parse.urlparse(value)
        if parsed.netloc != config.ena_object_host:
            return None
        return value if config.ena_object_scheme == "https" else None
    elif "://" in value:
        return None
    else:
        host_path = value
    host, separator, _path = host_path.partition("/")
    if not separator or host != config.ena_object_host:
        return None
    return f"{config.ena_object_scheme}://{host_path}"


def ena_pair_from_row(  # noqa: PLR0911 - guard clauses mirror ENA eligibility rules
    row: dict[str, object],
    config: StreamConfig,
) -> EnaPair | Sra:
    raw_urls = str(row.get("fastq_ftp") or "")
    if not raw_urls:
        return Sra("ena_fastq_absent")
    values = (
        raw_urls.split(";"),
        str(row.get("fastq_md5") or "").split(";"),
        str(row.get("fastq_bytes") or "").split(";"),
    )
    if any(len(items) != ENA_PAIR_FILE_COUNT for items in values):
        return Sra("ena_layout_incompatible")

    fastqs: dict[str, EnaFastq] = {}
    for raw_url, raw_md5, raw_bytes in zip(*values, strict=True):
        url = normalize_ena_url(raw_url, config)
        try:
            expected_md5 = bytes.fromhex(raw_md5)
            expected_bytes = int(raw_bytes)
        except ValueError:
            return Sra("ena_metadata_incomplete")
        if url is None or len(expected_md5) != MD5_DIGEST_BYTES or expected_bytes < 1:
            return Sra("ena_metadata_incomplete")
        basename = Path(urllib.parse.urlparse(url).path).name
        if basename.endswith("_1.fastq.gz"):
            mate = "r1"
        elif basename.endswith("_2.fastq.gz"):
            mate = "r2"
        else:
            return Sra("ena_layout_incompatible")
        if mate in fastqs:
            return Sra("ena_layout_incompatible")
        fastqs[mate] = EnaFastq(url, expected_bytes, expected_md5)
    if set(fastqs) != {"r1", "r2"}:
        return Sra("ena_layout_incompatible")
    return EnaPair((fastqs["r1"], fastqs["r2"]))


def classify_ena_report(payload: object, config: StreamConfig) -> EnaPair | Sra:
    if not isinstance(payload, list) or len(payload) != 1:
        return Sra("ena_metadata_incomplete")
    if not isinstance(payload[0], dict):
        return Sra("ena_metadata_incomplete")
    row = cast("dict[str, object]", payload[0])
    if row.get("run_accession") != config.accession:
        return Sra("ena_metadata_incomplete")
    platform = row.get("instrument_platform")
    layout = row.get("library_layout")
    if not platform or not layout:
        return Sra("ena_metadata_incomplete")
    if platform != "ILLUMINA" or layout != "PAIRED":
        requirement = "an Illumina run" if platform != "ILLUMINA" else "a paired run"
        observed = platform if platform != "ILLUMINA" else layout
        message = (
            f"stream_sra requires {requirement}; ENA reports {observed!r} "
            f"for {config.accession}"
        )
        raise UnsupportedLayoutError(message)
    return ena_pair_from_row(row, config)


def query_ena(config: StreamConfig) -> EnaPair | Sra:
    query = urllib.parse.urlencode(
        {
            "accession": config.accession,
            "result": "read_run",
            "format": "json",
            "fields": ",".join(ENA_FILE_REPORT_FIELDS),
        },
    )
    url = f"{config.ena_file_report_url}?{query}"
    for attempt in range(2):
        try:
            with urllib.request.urlopen(  # noqa: S310 - fixed or test-injected endpoint
                url,
                timeout=config.http_timeout_seconds,
            ) as response:
                return classify_ena_report(json.load(response), config)
        except urllib.error.HTTPError as error:
            retryable = (
                error.code == HTTP_TOO_MANY_REQUESTS
                or error.code >= HTTP_SERVER_ERROR_MIN
            )
            if attempt == 0 and retryable:
                time.sleep(config.retry_delay_seconds)
                continue
            return Sra("ena_lookup_unavailable")
        except (json.JSONDecodeError, UnicodeDecodeError):
            return Sra("ena_metadata_incomplete")
        except (
            http.client.HTTPException,
            urllib.error.URLError,
            TimeoutError,
            OSError,
        ):
            if attempt == 0:
                time.sleep(config.retry_delay_seconds)
    return Sra("ena_lookup_unavailable")


def close_responses(responses: tuple[http.client.HTTPResponse, ...]) -> None:
    for response in responses:
        with suppress(OSError):
            response.close()


def preflight_ena(
    config: StreamConfig,
    pair: EnaPair,
) -> tuple[http.client.HTTPResponse, http.client.HTTPResponse] | None:
    for attempt in range(2):
        opened: list[http.client.HTTPResponse] = []
        keep_open = False
        try:
            for fastq in pair.fastqs:
                response = urllib.request.urlopen(  # noqa: S310 - URL host validated above
                    fastq.url,
                    timeout=config.http_timeout_seconds,
                )
                opened.append(response)
            if len(opened) == ENA_PAIR_FILE_COUNT and all(
                response.status == HTTP_OK
                and response.headers.get("Content-Encoding") in (None, "identity")
                and response.headers.get("Content-Length") is not None
                and int(response.headers.get("Content-Length", ""))
                == fastq.expected_bytes
                and urllib.parse.urlparse(response.url)[:2]
                == urllib.parse.urlparse(fastq.url)[:2]
                for response, fastq in zip(opened, pair.fastqs, strict=True)
            ):
                keep_open = True
                return opened[0], opened[1]
        except (
            TypeError,
            ValueError,
            http.client.HTTPException,
            urllib.error.HTTPError,
            urllib.error.URLError,
            TimeoutError,
            OSError,
        ):
            pass
        finally:
            if not keep_open:
                close_responses(tuple(opened))
        if attempt == 0:
            time.sleep(config.retry_delay_seconds)
    return None


def deacon_command(config: StreamConfig, inputs: tuple[str, str]) -> list[str]:
    command = [
        "deacon",
        "filter",
        "--threads",
        str(config.threads),
        "--abs-threshold",
        str(config.abs_threshold),
        "--rel-threshold",
        str(config.rel_threshold),
        "--summary",
        str(config.summary),
        "--output",
        str(config.output),
        str(config.index),
        *inputs,
    ]
    if config.deplete:
        command.insert(2, "--deplete")
    return command


def stop_process(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    process.terminate()
    try:
        process.wait(timeout=PROCESS_TERMINATION_TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()


def run_checked(command: list[str], selection: Sra, cwd: Path | None = None) -> None:
    try:
        subprocess.run(command, cwd=cwd, check=True)  # noqa: S603
    except (OSError, subprocess.CalledProcessError) as error:
        message = f"backend process failed: {' '.join(command)}"
        raise BackendProcessError(
            message,
            backend="fastq-dump",
            selection_reason=selection.reason,
            component=command[0],
        ) from error


def run_sra(config: StreamConfig, selection: Sra) -> StreamResult:
    config.sra_directory.mkdir(parents=True, exist_ok=True)
    config.output.parent.mkdir(parents=True, exist_ok=True)
    config.summary.parent.mkdir(parents=True, exist_ok=True)
    run_checked(
        ["prefetch", "--max-size", "u", "--output-directory", ".", config.accession],
        selection,
        config.sra_directory,
    )
    archive = config.sra_directory / config.accession
    run_checked(["vdb-validate", str(archive)], selection)
    fastq_command = [
        "fastq-dump",
        "--split-spot",
        "--skip-technical",
        "--readids",
        "--stdout",
        str(archive),
    ]
    try:
        fastq = subprocess.Popen(  # noqa: S603
            fastq_command,
            stdout=subprocess.PIPE,
        )
    except OSError as error:
        message = "could not launch fastq-dump"
        raise BackendProcessError(
            message,
            backend="fastq-dump",
            selection_reason=selection.reason,
            component="fastq-dump",
        ) from error
    if fastq.stdout is None:
        stop_process(fastq)
        message = "fastq-dump stdout pipe was not created"
        raise BackendProcessError(
            message,
            backend="fastq-dump",
            selection_reason=selection.reason,
            component="fastq-dump",
        )
    try:
        deacon = subprocess.Popen(  # noqa: S603
            deacon_command(config, ("-", "-")),
            stdin=fastq.stdout,
        )
    except OSError as error:
        stop_process(fastq)
        message = "could not launch Deacon for SRA Toolkit input"
        raise BackendProcessError(
            message,
            backend="fastq-dump",
            selection_reason=selection.reason,
            component="deacon",
        ) from error
    finally:
        fastq.stdout.close()

    while fastq.poll() is None and deacon.poll() is None:
        time.sleep(SUPERVISOR_POLL_SECONDS)
    if fastq.poll() not in (None, 0):
        stop_process(deacon)
    elif deacon.poll() is not None and fastq.poll() is None:
        stop_process(fastq)
    fastq_returncode = fastq.wait()
    deacon_returncode = deacon.wait()
    if fastq_returncode != 0 or deacon_returncode != 0:
        message = (
            "streaming backend failed: "
            f"fastq-dump={fastq_returncode}, deacon={deacon_returncode}"
        )
        raise BackendProcessError(
            message,
            backend="fastq-dump",
            selection_reason=selection.reason,
            component="pipeline",
        )
    return StreamResult("fastq-dump", selection.reason)


def open_fifo_for_write(fifo: Path, stop: Event) -> BinaryIO:
    deadline = time.monotonic() + FIFO_OPEN_TIMEOUT_SECONDS
    while not stop.is_set() and time.monotonic() < deadline:
        try:
            descriptor = os.open(fifo, os.O_WRONLY | os.O_NONBLOCK)
            os.set_blocking(descriptor, True)
            return os.fdopen(descriptor, "wb", buffering=0)
        except OSError as error:
            if error.errno != errno.ENXIO:
                raise
            stop.wait(0.01)
    message = f"Deacon did not open ENA FIFO: {fifo}"
    raise EnaTransferError(message, fifo.stem)


def transfer_to_fifo(
    name: str,
    fastq: EnaFastq,
    response: http.client.HTTPResponse,
    fifo: Path,
    stop: Event,
) -> int:
    transferred = 0
    digest = hashlib.md5(usedforsecurity=False)
    buffer = bytearray(COPY_BUFFER_BYTES)
    view = memoryview(buffer)
    try:
        with open_fifo_for_write(fifo, stop) as output:
            while not stop.is_set():
                count = response.readinto(buffer)
                if not count:
                    break
                digest.update(view[:count])
                remaining = view[:count]
                while remaining:
                    written = output.write(remaining)
                    if written is None or written < 1:
                        message = "short write while forwarding ENA FASTQ"
                        raise EnaTransferError(message, name)
                    remaining = remaining[written:]
                transferred += count
    except (OSError, ValueError, http.client.HTTPException) as error:
        message = f"ENA {name} transfer stopped after {transferred} compressed bytes"
        raise EnaTransferError(message, name) from error
    if transferred != fastq.expected_bytes:
        message = (
            f"ENA {name} byte count mismatch: expected {fastq.expected_bytes}, "
            f"observed {transferred}"
        )
        raise IntegrityError(message, name)
    if digest.digest() != fastq.expected_md5:
        message = f"ENA {name} compressed MD5 mismatch"
        raise IntegrityError(message, name)
    return transferred


def first_failed(futures: tuple[Future[int], ...]) -> Future[int] | None:
    return next(
        (
            future
            for future in futures
            if future.done() and future.exception() is not None
        ),
        None,
    )


def run_ena(
    config: StreamConfig,
    pair: EnaPair,
    responses: tuple[http.client.HTTPResponse, http.client.HTTPResponse],
) -> StreamResult:
    config.output.parent.mkdir(parents=True, exist_ok=True)
    config.summary.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{config.accession}.ena-input.") as tmp:
        fifos = (Path(tmp) / "r1.fastq.gz", Path(tmp) / "r2.fastq.gz")
        for fifo in fifos:
            os.mkfifo(fifo)
        try:
            deacon = subprocess.Popen(  # noqa: S603
                deacon_command(config, (str(fifos[0]), str(fifos[1]))),
            )
        except OSError as error:
            message = "could not launch Deacon for ENA input"
            raise BackendProcessError(
                message,
                backend="ena",
                selection_reason="compatible_pair",
                component="deacon",
            ) from error

        stop = Event()
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = tuple(
                executor.submit(
                    transfer_to_fifo,
                    name,
                    fastq,
                    response,
                    fifo,
                    stop,
                )
                for name, fastq, response, fifo in zip(
                    ("r1", "r2"),
                    pair.fastqs,
                    responses,
                    fifos,
                    strict=True,
                )
            )
            while True:
                all_done = all(future.done() for future in futures)
                failed = first_failed(futures)
                if failed is not None:
                    stop.set()
                    close_responses(responses)
                    stop_process(deacon)
                    failed.result()
                if all_done:
                    break
                if deacon.poll() is not None:
                    stop.set()
                    close_responses(responses)
                    stop_process(deacon)
                    message = (
                        f"deacon filter exited {deacon.returncode} before both ENA "
                        "objects completed"
                    )
                    raise BackendProcessError(
                        message,
                        backend="ena",
                        selection_reason="compatible_pair",
                        component="deacon",
                    )
                stop.wait(SUPERVISOR_POLL_SECONDS)
            transferred = futures[0].result(), futures[1].result()

        deacon_returncode = deacon.wait()
        if deacon_returncode != 0:
            message = f"deacon filter exited {deacon_returncode} while reading ENA"
            raise BackendProcessError(
                message,
                backend="ena",
                selection_reason="compatible_pair",
                component="deacon",
            )
    return StreamResult("ena", "compatible_pair", transferred)


def run(config: StreamConfig) -> StreamResult:
    selection = query_ena(config)
    if isinstance(selection, Sra):
        return run_sra(config, selection)
    responses = preflight_ena(config, selection)
    if responses is None:
        return run_sra(config, Sra("ena_objects_unavailable"))
    try:
        return run_ena(config, selection, responses)
    finally:
        close_responses(responses)


def render_stream_event(fields: dict[str, object]) -> str:
    rendered = " ".join(f"{name}={value}" for name, value in fields.items())
    return f"nvd.sra_stream {rendered}"


def execute(config: StreamConfig) -> StreamResult:
    started = time.monotonic()
    try:
        result = run(config)
    except SraStreamError as error:
        logger.error(
            render_stream_event(
                {
                    "accession": config.accession,
                    "backend": error.backend,
                    "selection_reason": error.selection_reason,
                    "outcome": "error",
                    "duration_ms": round((time.monotonic() - started) * 1000),
                    "error_type": type(error).__name__,
                    "component": error.component,
                },
            ),
        )
        raise
    fields: dict[str, object] = {
        "accession": config.accession,
        "backend": result.backend,
        "selection_reason": result.selection_reason,
        "outcome": "success",
        "duration_ms": round((time.monotonic() - started) * 1000),
    }
    if result.transferred_bytes is not None:
        fields["transferred_r1_bytes"], fields["transferred_r2_bytes"] = (
            result.transferred_bytes
        )
    logger.info(render_stream_event(fields))
    return result


def configure_logging() -> None:
    logger.remove()
    logger.add(sys.stderr, format="{message}", level="INFO")


def parse_args(argv: list[str] | None = None) -> StreamConfig:
    parser = argparse.ArgumentParser(
        description="Prefer paired ENA FASTQs while streaming one SRA run into Deacon.",
    )
    parser.add_argument("--accession", required=True)
    parser.add_argument("--sra-directory", required=True, type=Path)
    parser.add_argument("--index", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--abs-threshold", required=True, type=int)
    parser.add_argument("--rel-threshold", required=True, type=float)
    parser.add_argument("--deplete", action="store_true")
    args = parser.parse_args(argv)
    return StreamConfig(
        accession=args.accession,
        sra_directory=args.sra_directory,
        index=args.index,
        output=args.output,
        summary=args.summary,
        threads=args.threads,
        abs_threshold=args.abs_threshold,
        rel_threshold=args.rel_threshold,
        deplete=args.deplete,
    )


def main() -> None:
    configure_logging()
    execute(parse_args())


if __name__ == "__main__":
    main()
