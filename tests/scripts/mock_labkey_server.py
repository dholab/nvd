"""In-memory mock of the LabKey query API and WebDAV surface.

The pipeline's upload/validate scripts talk to a LabKey server through the
``labkey`` Python client (``APIWrapper``) and to WebDAV through a small urllib
client (``bin/webdav_CLIent.py``). This module stands up just enough of both to
run the LIMS-enabled end-to-end test without a real server.

Request shapes were learned by reading the installed ``labkey`` client:

* ``APIWrapper(domain, container, api_key=...)`` forces ``use_ssl=True`` and
  builds URLs as ``https://<domain>/<container>/<controller>-<action>``. The
  domain must therefore be a bare ``host:port`` and the endpoint must speak
  TLS. We serve HTTPS with a self-signed cert for ``127.0.0.1`` and expect the
  caller to point ``SSL_CERT_FILE`` / ``REQUESTS_CA_BUNDLE`` at that cert so the
  client (requests) and the urllib WebDAV client both trust it.
* Before any request the client GETs ``login-whoami.api`` once per session and
  reads a ``CSRF`` field from the JSON body.
* ``select_rows`` POSTs form-encoded data to ``query-getQuery.api`` with
  ``schemaName``, ``query.queryName`` and filter params named
  ``query.<column>~<op>`` (e.g. ``query.experiment~eq``).
* ``insert_rows`` POSTs JSON to ``query-insertRows.api`` with
  ``{"schemaName", "queryName", "rows": [...]}``.
* ``delete_rows`` POSTs JSON to ``query-deleteRows.api`` with the same shape;
  ``validate_labkey.py`` inserts a probe row then deletes it by ``Key``.
* The WebDAV ``upload`` command issues a plain ``PUT`` of the file bytes.

The mock records every insert operation in :attr:`MockLabKeyServer.inserts`
(a list of ``{"schema", "table", "rows"}``) and serves ``select_rows`` from the
net stored state, so a second run sees the first run's inserted combos as
already present and skips them.
"""

from __future__ import annotations

import json
import socket
import ssl
import subprocess
import threading
import urllib.parse
from collections import defaultdict
from contextlib import closing, contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterator


def _free_port(host: str) -> int:
    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        sock.bind((host, 0))
        return sock.getsockname()[1]


class MockLabKeyServer:
    """Shared in-memory state plus URLs the pipeline params should point at."""

    def __init__(self, host: str, port: int, cert_file: Path) -> None:
        self.host = host
        self.port = port
        self.cert_file = cert_file
        self._lock = threading.Lock()
        self._key = 0
        # Net stored rows per table (list-style LabKey query), each carrying a
        # server-assigned integer "Key".
        self.stored: dict[str, list[dict[str, Any]]] = defaultdict(list)
        # Append-only log of insert operations: {"schema", "table", "rows"}.
        self.inserts: list[dict[str, Any]] = []
        # Append-only log of WebDAV PUT target paths.
        self.webdav_puts: list[str] = []

    @property
    def base_url(self) -> str:
        """Bare ``host:port`` for ``--labkey_server`` (client prepends https)."""
        return f"{self.host}:{self.port}"

    @property
    def server_url(self) -> str:
        return f"https://{self.host}:{self.port}"

    @property
    def webdav_url(self) -> str:
        """Full HTTPS URL for ``--labkey_webdav`` (urllib client keeps scheme)."""
        return f"https://{self.host}:{self.port}/webdav/"

    # --- operations invoked by the request handler (all lock-guarded) ---------

    def record_insert(
        self,
        schema: str,
        table: str,
        rows: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
        with self._lock:
            self.inserts.append({"schema": schema, "table": table, "rows": rows})
            stored_rows = []
            for row in rows:
                self._key += 1
                stored = {**row, "Key": self._key}
                self.stored[table].append(stored)
                stored_rows.append(stored)
            return stored_rows

    def record_delete(
        self,
        table: str,
        rows: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
        with self._lock:
            keys_to_drop = {row.get("Key") for row in rows if row.get("Key") is not None}
            if keys_to_drop:
                self.stored[table] = [
                    stored
                    for stored in self.stored[table]
                    if stored.get("Key") not in keys_to_drop
                ]
            return rows

    def select(self, table: str, filters: dict[str, str]) -> list[dict[str, Any]]:
        with self._lock:
            return [
                dict(row)
                for row in self.stored[table]
                if all(str(row.get(col)) == str(val) for col, val in filters.items())
            ]

    def record_put(self, path: str) -> None:
        with self._lock:
            self.webdav_puts.append(path)


def _make_handler(server: MockLabKeyServer) -> type[BaseHTTPRequestHandler]:  # noqa: C901
    class Handler(BaseHTTPRequestHandler):
        # Silence per-request stderr logging.
        def log_message(self, *_args: object) -> None:
            return

        def _send_json(self, payload: dict[str, Any], status: int = 200) -> None:
            body = json.dumps(payload).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _read_body(self) -> bytes:
            length = int(self.headers.get("Content-Length", 0) or 0)
            return self.rfile.read(length) if length else b""

        # -- GET: CSRF whoami handshake (and lenient WebDAV existence checks) ---
        def do_GET(self) -> None:
            if self.path.endswith("login-whoami.api"):
                self._send_json({"CSRF": "mock-csrf-token"})
                return
            # Any other GET (e.g. a WebDAV existence probe) reports "not found"
            # so the client falls through to creating/uploading.
            self.send_response(404)
            self.end_headers()

        def do_POST(self) -> None:
            body = self._read_body()
            if self.path.endswith("query-getQuery.api"):
                self._handle_select(body)
            elif self.path.endswith("query-insertRows.api"):
                self._handle_insert(body)
            elif self.path.endswith("query-deleteRows.api"):
                self._handle_delete(body)
            else:
                self._send_json({"rows": []})

        # -- WebDAV verbs -----------------------------------------------------
        def do_PUT(self) -> None:
            self._read_body()  # drain uploaded bytes
            server.record_put(self.path)
            self.send_response(201)
            self.end_headers()

        def do_MKCOL(self) -> None:
            self.send_response(201)
            self.end_headers()

        def do_PROPFIND(self) -> None:
            self.send_response(207)
            self.end_headers()

        # -- query handlers ---------------------------------------------------
        def _handle_select(self, body: bytes) -> None:
            params = urllib.parse.parse_qs(body.decode())
            table = params.get("query.queryName", [""])[0]
            schema = params.get("schemaName", [""])[0]
            filters: dict[str, str] = {}
            for key, values in params.items():
                if "~" not in key:
                    continue
                column = key.split("~", 1)[0].removeprefix("query.")
                filters[column] = values[0]
            rows = server.select(table, filters)
            self._send_json(
                {
                    "schemaName": schema,
                    "queryName": table,
                    "rowCount": len(rows),
                    "rows": rows,
                },
            )

        def _handle_insert(self, body: bytes) -> None:
            payload = json.loads(body.decode() or "{}")
            schema = payload.get("schemaName", "")
            table = payload.get("queryName", "")
            rows = payload.get("rows", [])
            stored_rows = server.record_insert(schema, table, rows)
            self._send_json(
                {
                    "schemaName": schema,
                    "queryName": table,
                    "command": "insert",
                    "rowsAffected": len(stored_rows),
                    "rows": stored_rows,
                },
            )

        def _handle_delete(self, body: bytes) -> None:
            payload = json.loads(body.decode() or "{}")
            schema = payload.get("schemaName", "")
            table = payload.get("queryName", "")
            rows = payload.get("rows", [])
            deleted = server.record_delete(table, rows)
            self._send_json(
                {
                    "schemaName": schema,
                    "queryName": table,
                    "command": "delete",
                    "rowsAffected": len(deleted),
                    "rows": deleted,
                },
            )

    return Handler


def _generate_self_signed_cert(directory: Path) -> tuple[Path, Path]:
    """Write a self-signed cert/key for 127.0.0.1 using the openssl CLI."""
    cert_file = directory / "mock_labkey_cert.pem"
    key_file = directory / "mock_labkey_key.pem"
    subprocess.run(  # noqa: S603
        [  # noqa: S607
            "openssl",
            "req",
            "-x509",
            "-newkey",
            "rsa:2048",
            "-nodes",
            "-keyout",
            str(key_file),
            "-out",
            str(cert_file),
            "-days",
            "2",
            "-subj",
            "/CN=127.0.0.1",
            "-addext",
            "subjectAltName=IP:127.0.0.1,DNS:localhost",
        ],
        check=True,
        capture_output=True,
    )
    return cert_file, key_file


@contextmanager
def mock_labkey_server(host: str = "127.0.0.1") -> Iterator[MockLabKeyServer]:
    """Run the mock LabKey/WebDAV endpoint for the duration of the context.

    Yields a :class:`MockLabKeyServer` whose ``base_url`` feeds
    ``--labkey_server`` (bare host:port), ``webdav_url`` feeds
    ``--labkey_webdav``, and ``cert_file`` should be exported via
    ``SSL_CERT_FILE`` / ``REQUESTS_CA_BUNDLE`` so the upload clients trust the
    self-signed TLS certificate.
    """
    with TemporaryDirectory(prefix="mock_labkey_") as tmp:
        tmp_path = Path(tmp)
        cert_file, key_file = _generate_self_signed_cert(tmp_path)

        port = _free_port(host)
        state = MockLabKeyServer(host, port, cert_file)
        httpd = ThreadingHTTPServer((host, port), _make_handler(state))

        ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
        ssl_context.load_cert_chain(certfile=str(cert_file), keyfile=str(key_file))
        httpd.socket = ssl_context.wrap_socket(httpd.socket, server_side=True)

        thread = threading.Thread(target=httpd.serve_forever, daemon=True)
        thread.start()
        try:
            yield state
        finally:
            httpd.shutdown()
            httpd.server_close()
            thread.join(timeout=5)
