"""
Database connection and initialization.

Provides context managers for state database connections with automatic
schema initialization and version checking.

NVD Home Directory:
    All NVD data is stored under ~/.nvd/ by default:
        ~/.nvd/
        ├── user.config         # Nextflow configuration
        ├── state.sqlite        # Run/sample/upload tracking
        └── taxdump/            # NCBI taxonomy data

Path Resolution (hierarchical fallback):
    Config:   NVD_CONFIG env var > ~/.nvd/user.config
    State:    NVD_STATE_DIR env var > ~/.nvd/
    Taxonomy: NVD_TAXONOMY_DB env var > ~/.nvd/taxdump/

Schema Migration Safety:
    When the database schema version doesn't match the expected version,
    the database must be recreated. This is a DESTRUCTIVE operation that
    loses all existing state (runs, samples, uploads). To prevent accidental
    data loss:

    1. Interactive mode (TTY): User is prompted for confirmation with a
       5-minute timeout. If no response, the operation is aborted.
    2. Non-interactive mode: A SchemaMismatchError is raised.
    3. Explicit opt-in: Pass allow_destructive_update=True to skip prompts.

    In all cases where destruction proceeds, a timestamped backup is created.
"""

from __future__ import annotations

import os
import select
import shutil
import sqlite3
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Generator

# The NVD home directory - all NVD data lives here by default
DEFAULT_NVD_HOME = Path.home() / ".nvd"

# Environment variables for overriding default paths
ENV_VAR_CONFIG = "NVD_CONFIG"
ENV_VAR = "NVD_STATE_DIR"  # Keep short name for backward compatibility
ENV_VAR_TAXONOMY = "NVD_TAXONOMY_DB"

# Default paths (all relative to NVD home)
DEFAULT_CONFIG_PATH = DEFAULT_NVD_HOME / "user.config"
DEFAULT_STATE_DIR = DEFAULT_NVD_HOME

EXPECTED_VERSION = 2

# Timeout for interactive confirmation prompt (seconds)
CONFIRMATION_TIMEOUT_SECONDS = 300  # 5 minutes


def utc_now_iso() -> str:
    """
    Return current UTC time as ISO8601 string with Z suffix.

    All timestamps in the NVD state database should use this function
    to ensure consistent UTC storage. This enables correct cross-timezone
    comparisons and avoids daylight saving time discontinuities.

    Returns:
        ISO8601 timestamp string like "2026-01-26T15:30:00Z"
    """
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


class SchemaMismatchError(Exception):
    """
    Raised when database schema version doesn't match and destructive update not allowed.

    This error indicates that the existing state database has an incompatible schema
    and cannot be used without being recreated. Recreating the database is destructive
    and requires explicit user consent.

    Attributes:
        db_path: Path to the database file
        current_version: Schema version in the existing database
        expected_version: Schema version required by the code
    """

    def __init__(
        self,
        db_path: Path,
        current_version: int,
        expected_version: int,
    ):
        self.db_path = db_path
        self.current_version = current_version
        self.expected_version = expected_version
        super().__init__(
            f"Schema mismatch: database at {db_path} has version {current_version}, "
            f"expected {expected_version}. Recreating the database would destroy all "
            f"existing state (runs, samples, uploads). To proceed, either:\n"
            f"  1. Run interactively and confirm when prompted\n"
            f"  2. Pass --update-db-destructive to explicitly allow destruction\n"
            f"  3. Manually backup and delete: {db_path}"
        )


class ResourceUnavailableError(Exception):
    """
    Base class for errors when a required resource is unavailable.

    This is the parent class for StateUnavailableError and TaxonomyUnavailableError.
    It provides a common interface for errors that occur when --sync is requested
    but the underlying resource (database, taxonomy, etc.) cannot be accessed.

    Attributes:
        resource_path: Path to the resource (database file, directory, etc.)
        operation: What operation was being attempted
        reason: Human-readable explanation of why the resource is unavailable
        original_error: The underlying exception that caused the failure
    """

    # Subclasses should override these
    resource_type: str = "Resource"
    recovery_hint: str = "remove the --sync flag"

    def __init__(
        self,
        resource_path: Path | str,
        operation: str,
        reason: str,
        original_error: Exception | None = None,
    ):
        self.resource_path = (
            Path(resource_path) if isinstance(resource_path, str) else resource_path
        )
        self.operation = operation
        self.reason = reason
        self.original_error = original_error
        super().__init__(
            f"{self.resource_type} unavailable during '{operation}': {reason}\n"
            f"Path: {resource_path}\n"
            f"The --sync flag requires this resource to be available.\n"
            f"To continue without it, {self.recovery_hint}."
        )


class StateUnavailableError(ResourceUnavailableError):
    """
    Raised when the state database is unavailable and --sync was requested.

    This error indicates that the state database could not be accessed (due to
    permission errors, corruption, missing files, etc.) and the user explicitly
    requested synchronization with the coordination system via --sync.

    Without --sync, this condition results in a warning and graceful degradation
    instead of an exception.

    Attributes:
        db_path: Path to the database file (or intended path)
        operation: What operation was being attempted
        reason: Human-readable explanation of why the database is unavailable
        original_error: The underlying exception that caused the failure
    """

    resource_type = "State database"
    recovery_hint = "remove the --sync flag"

    def __init__(
        self,
        db_path: Path | str,
        operation: str,
        reason: str,
        original_error: Exception | None = None,
    ):
        super().__init__(
            resource_path=db_path,
            operation=operation,
            reason=reason,
            original_error=original_error,
        )

    @property
    def db_path(self) -> Path:
        """Alias for resource_path for backward compatibility."""
        return self.resource_path


def format_state_warning(
    *,
    operation: str,
    context: str,
    error: Exception,
    db_path: Path | str,
    consequences: list[str],
) -> str:
    """
    Format a detailed warning message for state database unavailability.

    This produces a "wide event" style warning with full context about what
    was being attempted, what went wrong, and what the consequences are.

    Args:
        operation: Short description of the operation (e.g., "Registering run")
        context: Detailed context about what was being done
        error: The exception that occurred
        db_path: Path to the state database
        consequences: List of consequences of continuing without state

    Returns:
        Formatted multi-line warning string suitable for logging

    Example:
        >>> warning = format_state_warning(
        ...     operation="Registering run",
        ...     context="Run 'friendly_turing' with 12 samples (sample_set_id: abc123)",
        ...     error=PermissionError("Permission denied"),
        ...     db_path=Path("/home/user/.nvd/state.sqlite"),
        ...     consequences=[
        ...         "Duplicate detection DISABLED",
        ...         "Sample locking DISABLED",
        ...     ],
        ... )
    """
    # Validate inputs
    assert operation, "operation must not be empty"
    assert context, "context must not be empty"
    assert consequences, "consequences must not be empty"

    # Determine error type and provide helpful explanation
    error_type = type(error).__name__
    error_msg = str(error)

    # Provide more helpful explanations for common errors
    # Note: Order matters! OperationalError is a subclass of DatabaseError,
    # so we must check for OperationalError first.
    if isinstance(error, PermissionError):
        explanation = (
            "The database file exists but is not readable/writable by the current user."
        )
    elif isinstance(error, sqlite3.OperationalError):
        if "readonly" in error_msg.lower():
            explanation = "The database or its directory is on a read-only filesystem."
        elif "locked" in error_msg.lower():
            explanation = "The database is locked by another process."
        else:
            explanation = "A database operation failed unexpectedly."
    elif isinstance(error, sqlite3.DatabaseError):
        explanation = "The database file appears to be corrupted or is not a valid SQLite database."
    elif isinstance(error, FileNotFoundError):
        explanation = "The database file or its parent directory does not exist."
    elif isinstance(error, SchemaMismatchError):
        explanation = (
            f"The database schema version ({error.current_version}) doesn't match "
            f"the expected version ({error.expected_version})."
        )
    else:
        explanation = "An unexpected error occurred while accessing the database."

    # Build consequence lines
    consequence_lines = "\n".join(f"  - {c}" for c in consequences)

    return f"""
================================================================================
WARNING: STATE DATABASE UNAVAILABLE
================================================================================

Context:
  Operation: {operation}
  Details: {context}

Error:
  Type: {error_type}
  Message: {error_msg}
  Path: {db_path}
  Explanation: {explanation}

Consequence:
  Continuing WITHOUT coordination. This means:
{consequence_lines}

Resolution:
  To fix the underlying issue, check file permissions and database integrity.
  To require coordination (fail instead of continuing), add --sync to the command.

================================================================================
""".strip()


def get_config_path(explicit_path: Path | str | None = None) -> Path:
    """
    Resolve the NVD config file path using hierarchical fallback.

    Priority:
        1. Explicit path argument (from CLI)
        2. NVD_CONFIG environment variable
        3. Default: ~/.nvd/user.config

    Unlike state/taxonomy directories, this does NOT create the file
    if it doesn't exist - config files should be created explicitly.

    Args:
        explicit_path: Optional explicit path from CLI

    Returns:
        Resolved Path to config file
    """
    if explicit_path is not None:
        return Path(explicit_path)
    if ENV_VAR_CONFIG in os.environ:
        return Path(os.environ[ENV_VAR_CONFIG])
    return DEFAULT_CONFIG_PATH


def get_state_dir(explicit_path: Path | str | None = None) -> Path:
    """
    Resolve the NVD state directory using hierarchical fallback.

    Priority:
        1. Explicit path argument (from CLI or Nextflow param)
        2. NVD_STATE_DIR environment variable
        3. Default: ~/.nvd/

    The directory is created if it doesn't exist.

    Args:
        explicit_path: Optional explicit path from CLI/Nextflow

    Returns:
        Resolved Path to state directory
    """
    if explicit_path is not None:
        state_dir = Path(explicit_path)
    elif ENV_VAR in os.environ:
        state_dir = Path(os.environ[ENV_VAR])
    else:
        state_dir = DEFAULT_STATE_DIR

    state_dir.mkdir(parents=True, exist_ok=True)
    return state_dir


def get_state_db_path(state_dir: Path | str | None = None) -> Path:
    """Get path to the state database."""
    return get_state_dir(state_dir) / "state.sqlite"


def get_taxdump_dir(state_dir: Path | str | None = None) -> Path:
    """
    Get path to the taxdump directory containing .dmp files and taxonomy.sqlite.

    Priority:
        1. NVD_TAXONOMY_DB environment variable (for shared cluster installs)
        2. {state_dir}/taxdump (default)

    The NVD_TAXONOMY_DB override is useful for cluster environments where
    taxonomy data is pre-downloaded to a shared location.
    """
    if ENV_VAR_TAXONOMY in os.environ:
        return Path(os.environ[ENV_VAR_TAXONOMY])
    return get_state_dir(state_dir) / "taxdump"


def get_taxonomy_db_path(state_dir: Path | str | None = None) -> Path:
    """Get path to the taxonomy SQLite database (inside taxdump directory)."""
    return get_taxdump_dir(state_dir) / "taxonomy.sqlite"


def get_hits_dir(state_dir: Path | str | None = None) -> Path:
    """
    Get path to the hits directory containing parquet files.

    The hits directory uses Hive-style partitioning by month:
        Uncompacted: {state_dir}/hits/month=NULL/{sample_set_id}/{sample_id}/data.parquet
        Compacted:   {state_dir}/hits/month=YYYY-MM/data.parquet

    The directory is created if it doesn't exist.

    Args:
        state_dir: Optional explicit state directory. If None, resolves
                   via NVD_STATE_DIR env var, then ~/.nvd/

    Returns:
        Path to the hits directory
    """
    hits_dir = get_state_dir(state_dir) / "hits"
    hits_dir.mkdir(parents=True, exist_ok=True)
    return hits_dir


def _init_schema(conn: sqlite3.Connection) -> None:
    """Initialize database schema from schema.sql."""
    schema = files("py_nvd").joinpath("schema.sql").read_text()
    conn.executescript(schema)


def _get_version(conn: sqlite3.Connection) -> int:
    """Get the schema version from the database."""
    return conn.execute("PRAGMA user_version").fetchone()[0]


def _check_version(conn: sqlite3.Connection) -> bool:
    """Check if database schema version matches expected."""
    return _get_version(conn) == EXPECTED_VERSION


def _configure_connection(conn: sqlite3.Connection) -> None:
    """Apply standard connection configuration.

    Uses DELETE journal mode instead of WAL for network filesystem compatibility.
    WAL mode requires shared memory (mmap) which doesn't work reliably on NFS,
    GPFS, Lustre, and similar network filesystems. DELETE mode uses a traditional
    rollback journal that is slower but works correctly over the network.
    """
    conn.execute("PRAGMA journal_mode=DELETE")
    conn.execute("PRAGMA busy_timeout=30000")
    conn.row_factory = sqlite3.Row


def _get_backup_path(db_path: Path) -> Path:
    """Generate a timestamped backup path for the database."""
    timestamp = datetime.now(tz=timezone.utc).strftime("%Y%m%d_%H%M%S")
    return db_path.with_suffix(f".backup_{timestamp}.sqlite")


def _create_backup(db_path: Path) -> Path:
    """
    Create a timestamped backup of the database.

    Returns the path to the backup file.
    """
    backup_path = _get_backup_path(db_path)
    shutil.copy2(db_path, backup_path)
    return backup_path


def create_backup(state_dir: Path | str | None = None) -> Path:
    """
    Create a timestamped backup of the state database.

    This is the public API for creating backups, intended for use after
    successful pipeline runs to preserve a known-good state.

    Args:
        state_dir: Optional explicit state directory. If None, resolves
                   via NVD_STATE_DIR env var, then ~/.nvd/

    Returns:
        Path to the backup file.

    Raises:
        FileNotFoundError: If the database doesn't exist.
    """
    db_path = get_state_db_path(state_dir)
    if not db_path.exists():
        msg = f"Database not found: {db_path}"
        raise FileNotFoundError(msg)
    return _create_backup(db_path)


def _get_db_stats(db_path: Path) -> dict:
    """
    Get statistics about the database for display to user.

    Returns dict with file_size, version, and row counts per table.
    """
    stats = {
        "file_size": db_path.stat().st_size,
        "version": None,
        "tables": {},
    }

    conn = sqlite3.connect(db_path, timeout=5.0)
    try:
        stats["version"] = _get_version(conn)

        # Get list of tables
        tables = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
        ).fetchall()

        for (table_name,) in tables:
            count = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]  # noqa: S608
            stats["tables"][table_name] = count
    finally:
        conn.close()

    return stats


def _format_file_size(size_bytes: int) -> str:
    """Format bytes as human-readable size."""
    size = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def _prompt_for_confirmation(db_path: Path, timeout_seconds: int) -> bool:
    """
    Prompt user for confirmation before destructive schema update.

    Displays database statistics and waits for user input.
    Returns True if user confirms, False if user declines or timeout.

    Only works when stdin is a TTY. Returns False immediately if not.
    """
    if not sys.stdin.isatty():
        return False

    # Gather stats to show user what they're about to lose
    try:
        stats = _get_db_stats(db_path)
    except Exception:
        stats = {
            "file_size": db_path.stat().st_size,
            "version": "unknown",
            "tables": {},
        }

    # Build warning message
    print("\n" + "=" * 70, file=sys.stderr)
    print(
        "WARNING: DATABASE SCHEMA MISMATCH - DESTRUCTIVE UPDATE REQUIRED",
        file=sys.stderr,
    )
    print("=" * 70, file=sys.stderr)
    print(f"\nDatabase: {db_path}", file=sys.stderr)
    print(f"File size: {_format_file_size(stats['file_size'])}", file=sys.stderr)
    print(f"Current schema version: {stats['version']}", file=sys.stderr)
    print(f"Required schema version: {EXPECTED_VERSION}", file=sys.stderr)

    if stats["tables"]:
        print("\nData that will be PERMANENTLY LOST:", file=sys.stderr)
        for table, count in stats["tables"].items():
            print(f"  - {table}: {count} rows", file=sys.stderr)

    print("\nA backup will be created before deletion.", file=sys.stderr)
    print(f"You have {timeout_seconds // 60} minutes to respond.", file=sys.stderr)
    print(
        "\nType 'yes' to proceed, or anything else to abort: ", file=sys.stderr, end=""
    )
    sys.stderr.flush()

    # Wait for input with timeout (Unix-only via select)
    # On Windows, this will fall back to blocking input
    try:
        if sys.platform != "win32":
            ready, _, _ = select.select([sys.stdin], [], [], timeout_seconds)
            if not ready:
                print("\n\nTimeout reached. Aborting.", file=sys.stderr)
                return False
        response = input().strip().lower()
    except EOFError:
        return False

    return response == "yes"


def _handle_schema_mismatch(
    db_path: Path,
    current_version: int,
    allow_destructive_update: bool,
) -> Path | None:
    """
    Handle a schema version mismatch.

    Args:
        db_path: Path to the existing database
        current_version: Schema version in the existing database
        allow_destructive_update: If True, proceed without prompting

    Returns:
        Path to backup file if destruction was approved, None otherwise

    Raises:
        SchemaMismatchError: If destruction not allowed and not interactive
    """
    if allow_destructive_update:
        # Explicit opt-in: backup and proceed
        backup_path = _create_backup(db_path)
        print(
            f"Schema mismatch (v{current_version} -> v{EXPECTED_VERSION}). "
            f"Backup created: {backup_path}",
            file=sys.stderr,
        )
        db_path.unlink()
        return backup_path

    # Try interactive confirmation
    if sys.stdin.isatty():
        if _prompt_for_confirmation(db_path, CONFIRMATION_TIMEOUT_SECONDS):
            backup_path = _create_backup(db_path)
            print(f"\nBackup created: {backup_path}", file=sys.stderr)
            print("Proceeding with schema update...\n", file=sys.stderr)
            db_path.unlink()
            return backup_path
        else:
            print("\nAborted. Database unchanged.", file=sys.stderr)
            raise SchemaMismatchError(db_path, current_version, EXPECTED_VERSION)

    # Non-interactive and no explicit opt-in: refuse
    raise SchemaMismatchError(db_path, current_version, EXPECTED_VERSION)


@contextmanager
def connect(
    state_dir: Path | str | None = None,
    allow_destructive_update: bool = False,
) -> Generator[sqlite3.Connection, None, None]:
    """
    Context manager for state database connections.

    Creates the database and schema if missing. If the schema version doesn't
    match, behavior depends on allow_destructive_update and whether stdin is
    a TTY:

    - allow_destructive_update=True: Backup and recreate automatically
    - TTY available: Prompt user for confirmation (5-minute timeout)
    - Non-interactive: Raise SchemaMismatchError

    In all cases where the database is recreated, a timestamped backup is
    created first.

    Args:
        state_dir: Optional explicit state directory. If None, resolves
                   via NVD_STATE_DIR env var, then ~/.nvd/
        allow_destructive_update: If True, allow schema mismatch to trigger
                                  automatic backup and recreation without
                                  prompting. Use this when the CLI flag
                                  --update-db-destructive is passed.

    Raises:
        SchemaMismatchError: If schema mismatch detected and user doesn't
                             consent to destructive update.

    Usage:
        with connect() as conn:
            conn.execute("SELECT ...")

        # Or with explicit path
        with connect(state_dir="/scratch/nvd") as conn:
            conn.execute("SELECT ...")

        # Or allowing destructive updates (from CLI with --update-db-destructive)
        with connect(allow_destructive_update=True) as conn:
            conn.execute("SELECT ...")
    """
    db_path = get_state_db_path(state_dir)

    # Handle schema version mismatch
    if db_path.exists():
        conn = sqlite3.connect(db_path, timeout=30.0)
        if not _check_version(conn):
            current_version = _get_version(conn)
            conn.close()
            # This may raise SchemaMismatchError if user doesn't consent
            _handle_schema_mismatch(db_path, current_version, allow_destructive_update)
            # If we get here, the old database was backed up and deleted
            conn = sqlite3.connect(db_path, timeout=30.0)
            _configure_connection(conn)
            _init_schema(conn)
        else:
            _configure_connection(conn)
    else:
        conn = sqlite3.connect(db_path, timeout=30.0)
        _configure_connection(conn)
        _init_schema(conn)

    try:
        yield conn
    finally:
        conn.close()
