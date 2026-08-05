"""Unit tests for the mock LabKey/WebDAV server.

These exercise the mock through the very clients the pipeline uses (the
``labkey`` APIWrapper and the urllib WebDAV client), so a green run here proves
the request-shape matching is correct independently of the full pipeline.
"""

from __future__ import annotations

import os

from tests.scripts.mock_labkey_server import mock_labkey_server


def test_insert_then_select_sees_recorded_row() -> None:
    from labkey.api_wrapper import APIWrapper

    with mock_labkey_server() as mock:
        os.environ["SSL_CERT_FILE"] = str(mock.cert_file)
        os.environ["REQUESTS_CA_BUNDLE"] = str(mock.cert_file)
        try:
            api = APIWrapper(mock.base_url, "test", api_key="dummy")

            # Nothing recorded yet: the combo is absent.
            empty = api.query.select_rows(
                schema_name="lists",
                query_name="hits",
                filter_array=[],
            )
            assert empty["rows"] == []

            api.query.insert_rows(
                schema_name="lists",
                query_name="hits",
                rows=[
                    {"experiment": 1, "sample_id": "s1", "query_class": "single_read"},
                ],
            )

            # The insert op is logged with its rows carrying query_class.
            assert [op["table"] for op in mock.inserts] == ["hits"]
            assert mock.inserts[0]["rows"][0]["query_class"] == "single_read"

            # A filtered select now sees the recorded combo (dedup ledger).
            from labkey.query import QueryFilter

            present = api.query.select_rows(
                schema_name="lists",
                query_name="hits",
                filter_array=[
                    QueryFilter("experiment", 1, "eq"),
                    QueryFilter("sample_id", "s1", "eq"),
                    QueryFilter("query_class", "single_read", "eq"),
                ],
                max_rows=1,
            )
            assert present["rows"], "recorded combo should be visible to select_rows"

            # A different combo is still absent.
            other = api.query.select_rows(
                schema_name="lists",
                query_name="hits",
                filter_array=[QueryFilter("sample_id", "other", "eq")],
                max_rows=1,
            )
            assert other["rows"] == []
        finally:
            os.environ.pop("SSL_CERT_FILE", None)
            os.environ.pop("REQUESTS_CA_BUNDLE", None)


def test_validate_insert_then_delete_leaves_no_stored_row() -> None:
    """validate_labkey inserts a probe row then deletes it by Key."""
    from labkey.api_wrapper import APIWrapper

    with mock_labkey_server() as mock:
        os.environ["SSL_CERT_FILE"] = str(mock.cert_file)
        os.environ["REQUESTS_CA_BUNDLE"] = str(mock.cert_file)
        try:
            api = APIWrapper(mock.base_url, "test", api_key="dummy")
            result = api.query.insert_rows("lists", "hits", [{"Experiment": 123456789}])
            inserted = result["rows"][0]
            assert "Key" in inserted
            api.query.delete_rows("lists", "hits", [inserted])
            assert mock.stored["hits"] == []
        finally:
            os.environ.pop("SSL_CERT_FILE", None)
            os.environ.pop("REQUESTS_CA_BUNDLE", None)


def test_webdav_put_is_recorded() -> None:
    import base64
    import ssl
    import urllib.request

    with mock_labkey_server() as mock:
        token = base64.b64encode(b"apikey:dummy").decode()
        request = urllib.request.Request(  # noqa: S310
            mock.webdav_url + "1/sample/nvd/file.gz",
            data=b"payload",
            method="PUT",
            headers={"Authorization": f"Basic {token}"},
        )
        # Trust the mock's self-signed cert via an explicit context rather than
        # the SSL_CERT_FILE env var: env-based trust is only honored if set
        # before the process's first SSL context is created, which is unreliable
        # in a shared pytest process (passes on macOS, fails on Linux CI).
        context = ssl.create_default_context(cafile=str(mock.cert_file))
        with urllib.request.urlopen(request, context=context) as response:  # noqa: S310
            assert response.getcode() == 201
        assert any(path.endswith("file.gz") for path in mock.webdav_puts)
