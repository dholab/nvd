"""
Tests for py_nvd.db module.

Tests for database connection, schema management, and graceful degradation helpers.
"""

import sqlite3
from pathlib import Path

import pytest

from py_nvd.db import (
    SchemaMismatchError,
    StateUnavailableError,
    format_state_warning,
)


class TestStateUnavailableError:
    """Tests for StateUnavailableError exception."""

    def test_error_attributes(self, tmp_path):
        """Verify exception has correct attributes."""
        original = ValueError("original error")
        error = StateUnavailableError(
            db_path=tmp_path / "state.sqlite",
            operation="test operation",
            reason="test reason",
            original_error=original,
        )

        assert error.db_path == tmp_path / "state.sqlite"
        assert error.operation == "test operation"
        assert error.reason == "test reason"
        assert error.original_error is original

    def test_error_attributes_with_string_path(self, tmp_path):
        """Verify exception handles string path correctly."""
        error = StateUnavailableError(
            db_path=str(tmp_path / "state.sqlite"),
            operation="test operation",
            reason="test reason",
        )

        assert error.db_path == tmp_path / "state.sqlite"
        assert isinstance(error.db_path, Path)

    def test_error_message_format(self, tmp_path):
        """Verify error message is helpful."""
        error = StateUnavailableError(
            db_path=tmp_path / "state.sqlite",
            operation="Registering run",
            reason="Permission denied",
        )

        message = str(error)
        assert "Registering run" in message
        assert "Permission denied" in message
        assert "--sync" in message
        assert str(tmp_path / "state.sqlite") in message

    def test_error_without_original_error(self, tmp_path):
        """Verify exception works without original_error."""
        error = StateUnavailableError(
            db_path=tmp_path / "state.sqlite",
            operation="test",
            reason="test reason",
        )

        assert error.original_error is None


class TestFormatStateWarning:
    """Tests for format_state_warning() function."""

    def test_format_warning_permission_error(self, tmp_path):
        """Verify explanation for PermissionError."""
        error = PermissionError("Permission denied")
        warning = format_state_warning(
            operation="Registering run",
            context="Run 'test_run' with 5 samples",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["Duplicate detection DISABLED"],
        )

        assert "STATE DATABASE UNAVAILABLE" in warning
        assert "PermissionError" in warning
        assert "not readable/writable" in warning
        assert str(tmp_path / "state.sqlite") in warning

    def test_format_warning_database_error(self, tmp_path):
        """Verify explanation for sqlite3.DatabaseError."""
        error = sqlite3.DatabaseError("file is not a database")
        warning = format_state_warning(
            operation="Opening database",
            context="Connecting to state database",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["State tracking DISABLED"],
        )

        assert "DatabaseError" in warning
        assert "corrupted" in warning

    def test_format_warning_operational_error_readonly(self, tmp_path):
        """Verify explanation for readonly OperationalError."""
        error = sqlite3.OperationalError("attempt to write a readonly database")
        warning = format_state_warning(
            operation="Writing to database",
            context="Inserting run record",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["Run registration DISABLED"],
        )

        assert "OperationalError" in warning
        assert "read-only filesystem" in warning

    def test_format_warning_operational_error_locked(self, tmp_path):
        """Verify explanation for locked OperationalError."""
        error = sqlite3.OperationalError("database is locked")
        warning = format_state_warning(
            operation="Acquiring lock",
            context="Locking sample for processing",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["Sample locking DISABLED"],
        )

        assert "OperationalError" in warning
        assert "locked by another process" in warning

    def test_format_warning_operational_error_generic(self, tmp_path):
        """Verify explanation for generic OperationalError."""
        error = sqlite3.OperationalError("some other error")
        warning = format_state_warning(
            operation="Database operation",
            context="Performing query",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["Operation DISABLED"],
        )

        assert "OperationalError" in warning
        assert "failed unexpectedly" in warning

    def test_format_warning_file_not_found(self, tmp_path):
        """Verify explanation for FileNotFoundError."""
        error = FileNotFoundError("No such file or directory")
        warning = format_state_warning(
            operation="Opening database",
            context="Connecting to state database",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["State tracking DISABLED"],
        )

        assert "FileNotFoundError" in warning
        assert "does not exist" in warning

    def test_format_warning_schema_mismatch(self, tmp_path):
        """Verify explanation for SchemaMismatchError."""
        error = SchemaMismatchError(
            db_path=tmp_path / "state.sqlite",
            current_version=1,
            expected_version=2,
        )
        warning = format_state_warning(
            operation="Opening database",
            context="Checking schema version",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["Database needs migration"],
        )

        assert "SchemaMismatchError" in warning
        assert "version (1)" in warning
        assert "version (2)" in warning

    def test_format_warning_unknown_error(self, tmp_path):
        """Verify fallback explanation for unknown errors."""
        error = RuntimeError("Something unexpected happened")
        warning = format_state_warning(
            operation="Unknown operation",
            context="Doing something",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["Unknown consequence"],
        )

        assert "RuntimeError" in warning
        assert "unexpected error" in warning

    def test_format_warning_includes_consequences(self, tmp_path):
        """Verify consequences are listed in warning."""
        error = PermissionError("Permission denied")
        warning = format_state_warning(
            operation="Registering run",
            context="Run 'test_run'",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=[
                "Duplicate detection DISABLED",
                "Sample locking DISABLED",
                "Run provenance tracking DISABLED",
            ],
        )

        assert "Duplicate detection DISABLED" in warning
        assert "Sample locking DISABLED" in warning
        assert "Run provenance tracking DISABLED" in warning

    def test_format_warning_includes_resolution(self, tmp_path):
        """Verify resolution instructions are included."""
        error = PermissionError("Permission denied")
        warning = format_state_warning(
            operation="Test",
            context="Test context",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["Test consequence"],
        )

        assert "Resolution:" in warning
        assert "--sync" in warning
        assert "file permissions" in warning

    def test_format_warning_structure(self, tmp_path):
        """Verify overall warning structure."""
        error = ValueError("test error")
        warning = format_state_warning(
            operation="Test operation",
            context="Test context details",
            error=error,
            db_path=tmp_path / "state.sqlite",
            consequences=["Consequence 1", "Consequence 2"],
        )

        # Check section headers
        assert "=" * 10 in warning  # Separator lines
        assert "WARNING: STATE DATABASE UNAVAILABLE" in warning
        assert "Context:" in warning
        assert "Operation:" in warning
        assert "Error:" in warning
        assert "Consequence:" in warning
        assert "Resolution:" in warning

    def test_format_warning_with_string_path(self, tmp_path):
        """Verify warning works with string path."""
        error = ValueError("test")
        warning = format_state_warning(
            operation="Test",
            context="Test",
            error=error,
            db_path=str(tmp_path / "state.sqlite"),
            consequences=["Test"],
        )

        assert str(tmp_path / "state.sqlite") in warning
