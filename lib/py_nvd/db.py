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
    """Apply standard connection configuration."""
    conn.execute("PRAGMA journal_mode=WAL")
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
