"""Tests for state database import/export functionality."""

import json
import shutil
import tarfile
from pathlib import Path

import pytest
from typer.testing import CliRunner

from py_nvd import __version__
from py_nvd.cli.app import app
from py_nvd.cli.commands.state import (
    _build_export_manifest,
    _do_merge_import,
    _do_replace_import,
    _find_conflicts,
    _get_table_counts,
    _print_conflicts,
    _validate_archive,
)
from py_nvd.db import EXPECTED_VERSION, connect
from py_nvd.state import (
    compute_sample_set_id,
    register_processed_sample,
    register_run,
)

runner = CliRunner()


def _make_archive(tmp_path: Path, state_dir: Path) -> Path:
    """Create a valid .nvd-state archive from an existing state directory.

    Builds the manifest from the live database, bundles it with the
    state.sqlite file into a tar.gz, and returns the archive path.
    """
    db_path = state_dir / "state.sqlite"
    archive_path = tmp_path / "test.nvd-state"

    with connect(state_dir) as conn:
        manifest = _build_export_manifest(conn)

    staging = tmp_path / "archive_staging"
    staging.mkdir()
    shutil.copy2(db_path, staging / "state.sqlite")
    (staging / "manifest.json").write_text(json.dumps(manifest, default=str))

    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(staging / "state.sqlite", arcname="state.sqlite")
        tar.add(staging / "manifest.json", arcname="manifest.json")

    shutil.rmtree(staging)
    return archive_path


def _seed_lab_data(
    state_dir: Path,
    run_id: str,
    sample_ids: list[str],
) -> str:
    """Seed a state database with a run and processed samples.

    Returns the sample_set_id for the registered run.
    """
    sample_set_id = compute_sample_set_id(sample_ids)
    register_run(run_id, sample_set_id, state_dir=state_dir)
    for sid in sample_ids:
        register_processed_sample(sid, sample_set_id, run_id, state_dir=state_dir)
    return sample_set_id


class TestGetTableCounts:
    """Tests for _get_table_counts helper."""

    def test_empty_database(self, tmp_path: Path) -> None:
        """All counts are zero in a freshly initialized database."""
        with connect(tmp_path) as conn:
            counts = _get_table_counts(conn)

        assert all(v == 0 for v in counts.values())
        assert "runs" in counts
        assert "processed_samples" in counts
        assert "uploads" in counts

    def test_reflects_seeded_data(self, tmp_path: Path) -> None:
        """Counts reflect data inserted into the database."""
        _seed_lab_data(tmp_path, "run_001", ["s1", "s2", "s3"])

        with connect(tmp_path) as conn:
            counts = _get_table_counts(conn)

        assert counts["runs"] == 1
        assert counts["processed_samples"] == 3


class TestBuildExportManifest:
    """Tests for _build_export_manifest helper."""

    def test_contains_required_fields(self, tmp_path: Path) -> None:
        """Manifest contains all fields required by _validate_archive."""
        with connect(tmp_path) as conn:
            manifest = _build_export_manifest(conn)

        for field in ("schema_version", "exported_at", "nvd_version", "counts"):
            assert field in manifest

    def test_schema_version_matches_expected(self, tmp_path: Path) -> None:
        """Manifest schema_version matches the current EXPECTED_VERSION."""
        with connect(tmp_path) as conn:
            manifest = _build_export_manifest(conn)

        assert manifest["schema_version"] == EXPECTED_VERSION

    def test_nvd_version_matches_package(self, tmp_path: Path) -> None:
        """Manifest nvd_version matches the installed package version."""
        with connect(tmp_path) as conn:
            manifest = _build_export_manifest(conn)

        assert manifest["nvd_version"] == __version__

    def test_counts_reflect_data(self, tmp_path: Path) -> None:
        """Manifest counts reflect the actual database contents."""
        _seed_lab_data(tmp_path, "run_001", ["s1", "s2"])

        with connect(tmp_path) as conn:
            manifest = _build_export_manifest(conn)

        assert manifest["counts"]["runs"] == 1
        assert manifest["counts"]["processed_samples"] == 2


class TestValidateArchive:
    """Tests for _validate_archive — archive extraction and validation."""

    def test_extracts_and_returns_manifest(self, tmp_path: Path) -> None:
        """Valid archive is extracted and manifest is returned."""
        _seed_lab_data(tmp_path / "source", "run_001", ["s1"])
        archive_path = _make_archive(tmp_path, tmp_path / "source")

        manifest, tmpdir_path = _validate_archive(archive_path)

        try:
            assert "schema_version" in manifest
            assert "counts" in manifest
            assert (tmpdir_path / "state.sqlite").exists()
            assert (tmpdir_path / "manifest.json").exists()
        finally:
            shutil.rmtree(tmpdir_path)

    def test_fails_when_archive_missing(self, tmp_path: Path) -> None:
        """Raises SystemExit when archive file does not exist."""
        with pytest.raises(SystemExit):
            _validate_archive(tmp_path / "nonexistent.nvd-state")

    def test_fails_when_state_sqlite_missing(self, tmp_path: Path) -> None:
        """Raises SystemExit when archive lacks state.sqlite."""
        archive_path = tmp_path / "bad.nvd-state"
        manifest = {
            "schema_version": EXPECTED_VERSION,
            "exported_at": "now",
            "nvd_version": __version__,
            "counts": {},
        }

        staging = tmp_path / "staging"
        staging.mkdir()
        (staging / "manifest.json").write_text(json.dumps(manifest))

        with tarfile.open(archive_path, "w:gz") as tar:
            tar.add(staging / "manifest.json", arcname="manifest.json")

        shutil.rmtree(staging)

        with pytest.raises(SystemExit):
            _validate_archive(archive_path)

    def test_fails_when_manifest_missing(self, tmp_path: Path) -> None:
        """Raises SystemExit when archive lacks manifest.json."""
        # Create a DB to bundle
        with connect(tmp_path / "source"):
            pass
        archive_path = tmp_path / "bad.nvd-state"

        staging = tmp_path / "staging"
        staging.mkdir()
        shutil.copy2(tmp_path / "source" / "state.sqlite", staging / "state.sqlite")

        with tarfile.open(archive_path, "w:gz") as tar:
            tar.add(staging / "state.sqlite", arcname="state.sqlite")

        shutil.rmtree(staging)

        with pytest.raises(SystemExit):
            _validate_archive(archive_path)

    def test_fails_on_malformed_manifest(self, tmp_path: Path) -> None:
        """Raises SystemExit when manifest.json is not valid JSON."""
        with connect(tmp_path / "source"):
            pass
        archive_path = tmp_path / "bad.nvd-state"

        staging = tmp_path / "staging"
        staging.mkdir()
        shutil.copy2(tmp_path / "source" / "state.sqlite", staging / "state.sqlite")
        (staging / "manifest.json").write_text("not json {{{")

        with tarfile.open(archive_path, "w:gz") as tar:
            tar.add(staging / "state.sqlite", arcname="state.sqlite")
            tar.add(staging / "manifest.json", arcname="manifest.json")

        shutil.rmtree(staging)

        with pytest.raises(SystemExit):
            _validate_archive(archive_path)

    def test_fails_on_missing_manifest_fields(self, tmp_path: Path) -> None:
        """Raises SystemExit when manifest is missing required fields."""
        with connect(tmp_path / "source"):
            pass
        archive_path = tmp_path / "bad.nvd-state"

        staging = tmp_path / "staging"
        staging.mkdir()
        shutil.copy2(tmp_path / "source" / "state.sqlite", staging / "state.sqlite")
        (staging / "manifest.json").write_text(json.dumps({"partial": True}))

        with tarfile.open(archive_path, "w:gz") as tar:
            tar.add(staging / "state.sqlite", arcname="state.sqlite")
            tar.add(staging / "manifest.json", arcname="manifest.json")

        shutil.rmtree(staging)

        with pytest.raises(SystemExit):
            _validate_archive(archive_path)


class TestFindConflicts:
    """Tests for _find_conflicts — PK collision detection between databases."""

    def test_no_conflicts_with_disjoint_data(self, tmp_path: Path) -> None:
        """Two databases with completely different data produce no conflicts."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        _seed_lab_data(local_dir, "run_local", ["s1", "s2"])
        _seed_lab_data(backup_dir, "run_backup", ["s3", "s4"])

        with connect(local_dir) as conn:
            conflicts = _find_conflicts(conn, backup_dir / "state.sqlite")

        assert conflicts == {}

    def test_no_conflicts_when_pk_matches_with_identical_data(
        self,
        tmp_path: Path,
    ) -> None:
        """Matching PKs with identical display columns are not conflicts."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        _seed_lab_data(local_dir, "run_001", ["s1"])
        _seed_lab_data(backup_dir, "run_001", ["s1"])

        with connect(local_dir) as conn:
            conflicts = _find_conflicts(conn, local_dir / "state.sqlite")

        assert conflicts == {}

    def test_detects_run_conflict(self, tmp_path: Path) -> None:
        """Detects conflict when same run_id has different status."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        sample_set_id = compute_sample_set_id(["s1"])
        register_run("shared_run", sample_set_id, state_dir=local_dir)

        # In the backup, create the same run_id but with a different sample_set_id
        # so the started_at or status will differ
        other_set_id = compute_sample_set_id(["s2"])
        register_run("shared_run", other_set_id, state_dir=backup_dir)

        # Manually change the status in the backup to force a display-column diff
        with connect(backup_dir) as conn:
            conn.execute(
                "UPDATE runs SET status = 'completed' WHERE run_id = 'shared_run'",
            )
            conn.commit()

        with connect(local_dir) as conn:
            conflicts = _find_conflicts(conn, backup_dir / "state.sqlite")

        assert "runs" in conflicts
        assert len(conflicts["runs"]) == 1
        assert conflicts["runs"][0]["pk"] == {"run_id": "shared_run"}

    def test_detects_composite_pk_conflict(self, tmp_path: Path) -> None:
        """Detects conflict on a table with composite primary key."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        sample_set_id = compute_sample_set_id(["s1"])
        register_run("run_001", sample_set_id, state_dir=local_dir)
        register_run("run_001", sample_set_id, state_dir=backup_dir)

        register_processed_sample("s1", sample_set_id, "run_001", state_dir=local_dir)
        register_processed_sample("s1", sample_set_id, "run_001", state_dir=backup_dir)

        # Change status in backup to create a conflict on (sample_id, sample_set_id)
        with connect(backup_dir) as conn:
            conn.execute(
                "UPDATE processed_samples SET status = 'failed' WHERE sample_id = 's1'",
            )
            conn.commit()

        with connect(local_dir) as conn:
            conflicts = _find_conflicts(conn, backup_dir / "state.sqlite")

        assert "processed_samples" in conflicts
        assert len(conflicts["processed_samples"]) == 1
        conflict = conflicts["processed_samples"][0]
        assert conflict["pk"]["sample_id"] == "s1"
        assert conflict["local"]["status"] == "completed"
        assert conflict["incoming"]["status"] == "failed"

    def test_empty_databases_no_conflicts(self, tmp_path: Path) -> None:
        """Two empty databases produce no conflicts."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        # Just initialize both databases (connect creates schema)
        with connect(local_dir):
            pass
        with connect(backup_dir):
            pass

        with connect(local_dir) as conn:
            conflicts = _find_conflicts(conn, backup_dir / "state.sqlite")

        assert conflicts == {}


class TestPrintConflicts:
    """Tests for _print_conflicts — conflict report formatting."""

    def test_json_output_structure(self, capsys: pytest.CaptureFixture[str]) -> None:
        """JSON output contains error, conflicts, and conflict_count."""
        conflicts = {
            "runs": [
                {
                    "pk": {"run_id": "run_001"},
                    "local": {"started_at": "2024-01-01", "status": "running"},
                    "incoming": {"started_at": "2024-06-01", "status": "completed"},
                },
            ],
        }

        _print_conflicts(conflicts, json_output=True)

        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert data["error"] == "Merge aborted due to conflicts"
        assert data["conflict_count"] == 1
        assert "runs" in data["conflicts"]

    def test_human_output_contains_details(
        self,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Human-readable output includes table name, PK, and local/incoming values."""
        conflicts = {
            "runs": [
                {
                    "pk": {"run_id": "run_001"},
                    "local": {"started_at": "2024-01-01", "status": "running"},
                    "incoming": {"started_at": "2024-06-01", "status": "completed"},
                },
            ],
        }

        _print_conflicts(conflicts, json_output=False)

        captured = capsys.readouterr()
        output = captured.out
        assert "run_001" in output
        assert "running" in output
        assert "completed" in output


class TestDoMergeImport:
    """Tests for _do_merge_import — INSERT WHERE NOT EXISTS semantics."""

    def test_imports_new_rows(self, tmp_path: Path) -> None:
        """New rows from backup are inserted into the local database."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        _seed_lab_data(local_dir, "run_local", ["s1"])
        _seed_lab_data(backup_dir, "run_backup", ["s2"])

        with connect(local_dir) as conn:
            imported = _do_merge_import(conn, backup_dir / "state.sqlite")

        assert imported["runs"] == 1
        assert imported["processed_samples"] == 1

    def test_skips_existing_rows(self, tmp_path: Path) -> None:
        """Rows that already exist in local are not duplicated."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        _seed_lab_data(local_dir, "run_001", ["s1"])
        _seed_lab_data(backup_dir, "run_001", ["s1"])

        with connect(local_dir) as conn:
            imported = _do_merge_import(conn, backup_dir / "state.sqlite")

        assert imported["runs"] == 0
        assert imported["processed_samples"] == 0

    def test_mixed_new_and_existing(self, tmp_path: Path) -> None:
        """Only new rows are imported when backup has both new and existing data."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        shared_set = compute_sample_set_id(["s1"])
        register_run("shared_run", shared_set, state_dir=local_dir)
        register_processed_sample("s1", shared_set, "shared_run", state_dir=local_dir)

        # Backup has the same run plus a new one
        register_run("shared_run", shared_set, state_dir=backup_dir)
        register_processed_sample("s1", shared_set, "shared_run", state_dir=backup_dir)
        new_set = compute_sample_set_id(["s2"])
        register_run("new_run", new_set, state_dir=backup_dir)
        register_processed_sample("s2", new_set, "new_run", state_dir=backup_dir)

        with connect(local_dir) as conn:
            imported = _do_merge_import(conn, backup_dir / "state.sqlite")

        assert imported["runs"] == 1  # new_run
        assert imported["processed_samples"] == 1  # s2

    def test_returns_per_table_counts(self, tmp_path: Path) -> None:
        """Return dict has an entry for every table in _IMPORT_TABLES."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        with connect(local_dir):
            pass
        with connect(backup_dir):
            pass

        with connect(local_dir) as conn:
            imported = _do_merge_import(conn, backup_dir / "state.sqlite")

        assert "runs" in imported
        assert "processed_samples" in imported
        assert "uploads" in imported
        assert "databases" in imported
        assert "taxonomy_versions" in imported
        assert "sra_cache" in imported
        assert "sample_locks" in imported

    def test_local_data_unchanged_after_merge(self, tmp_path: Path) -> None:
        """Existing local rows are not modified by the merge."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        _seed_lab_data(local_dir, "run_local", ["s1"])
        _seed_lab_data(backup_dir, "run_backup", ["s2"])

        with connect(local_dir) as conn:
            # Capture local state before merge
            local_run_before = conn.execute(
                "SELECT * FROM runs WHERE run_id = 'run_local'",
            ).fetchone()

            _do_merge_import(conn, backup_dir / "state.sqlite")

            # Verify local run is unchanged
            local_run_after = conn.execute(
                "SELECT * FROM runs WHERE run_id = 'run_local'",
            ).fetchone()

        assert dict(local_run_before) == dict(local_run_after)


class TestDoReplaceImport:
    """Tests for _do_replace_import — database file swap."""

    def test_replaces_database(self, tmp_path: Path) -> None:
        """Backup database replaces the current one."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        _seed_lab_data(local_dir, "run_local", ["s1"])
        _seed_lab_data(backup_dir, "run_backup", ["s2"])

        _do_replace_import(backup_dir / "state.sqlite", local_dir / "state.sqlite")

        # Local DB should now contain backup's data
        with connect(local_dir) as conn:
            runs = conn.execute("SELECT run_id FROM runs").fetchall()

        run_ids = {r["run_id"] for r in runs}
        assert "run_backup" in run_ids
        assert "run_local" not in run_ids

    def test_bak_cleaned_up_on_success(self, tmp_path: Path) -> None:
        """The .bak file is removed after successful replace."""
        local_dir = tmp_path / "local"
        backup_dir = tmp_path / "backup"

        with connect(local_dir):
            pass
        with connect(backup_dir):
            pass

        _do_replace_import(backup_dir / "state.sqlite", local_dir / "state.sqlite")

        assert not (local_dir / "state.sqlite.bak").exists()

    def test_creates_db_when_current_missing(self, tmp_path: Path) -> None:
        """Replace works even when no current database exists."""
        local_dir = tmp_path / "local"
        local_dir.mkdir(parents=True)
        backup_dir = tmp_path / "backup"

        _seed_lab_data(backup_dir, "run_backup", ["s1"])

        _do_replace_import(backup_dir / "state.sqlite", local_dir / "state.sqlite")

        assert (local_dir / "state.sqlite").exists()
        with connect(local_dir) as conn:
            runs = conn.execute("SELECT run_id FROM runs").fetchall()
        assert len(runs) == 1


class TestExportImportRoundTrip:
    """End-to-end round-trip tests: export from one lab, import into another."""

    def test_merge_round_trip(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Data exported from lab A can be merge-imported into lab B."""
        lab_a = tmp_path / "lab_a"
        lab_b = tmp_path / "lab_b"

        # Seed lab A
        _seed_lab_data(lab_a, "run_a", ["s1", "s2"])

        # Export from lab A
        monkeypatch.setenv("NVD_STATE_DIR", str(lab_a))
        archive = tmp_path / "export.nvd-state"
        result = runner.invoke(app, ["state", "export", str(archive)])
        assert result.exit_code == 0
        assert archive.exists()

        # Seed lab B with different data
        _seed_lab_data(lab_b, "run_b", ["s3", "s4"])

        # Import into lab B
        monkeypatch.setenv("NVD_STATE_DIR", str(lab_b))
        result = runner.invoke(app, ["state", "import", str(archive), "--yes"])
        assert result.exit_code == 0

        # Verify lab B now has both labs' data
        with connect(lab_b) as conn:
            runs = conn.execute("SELECT run_id FROM runs").fetchall()
            samples = conn.execute("SELECT sample_id FROM processed_samples").fetchall()

        run_ids = {r["run_id"] for r in runs}
        sample_ids = {s["sample_id"] for s in samples}
        assert run_ids == {"run_a", "run_b"}
        assert sample_ids == {"s1", "s2", "s3", "s4"}

    def test_replace_round_trip(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Replace import discards local data and adopts the archive."""
        lab_a = tmp_path / "lab_a"
        lab_b = tmp_path / "lab_b"

        _seed_lab_data(lab_a, "run_a", ["s1"])

        # Export from lab A
        monkeypatch.setenv("NVD_STATE_DIR", str(lab_a))
        archive = tmp_path / "export.nvd-state"
        result = runner.invoke(app, ["state", "export", str(archive)])
        assert result.exit_code == 0

        # Seed lab B with different data
        _seed_lab_data(lab_b, "run_b", ["s2"])

        # Replace import into lab B
        monkeypatch.setenv("NVD_STATE_DIR", str(lab_b))
        result = runner.invoke(
            app,
            ["state", "import", str(archive), "--replace", "--yes"],
        )
        assert result.exit_code == 0

        # Lab B should only have lab A's data now
        with connect(lab_b) as conn:
            runs = conn.execute("SELECT run_id FROM runs").fetchall()

        run_ids = {r["run_id"] for r in runs}
        assert run_ids == {"run_a"}
        assert "run_b" not in run_ids

    def test_merge_aborts_on_conflict(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Merge import aborts cleanly when conflicts are detected."""
        lab_a = tmp_path / "lab_a"
        lab_b = tmp_path / "lab_b"

        # Both labs have the same run_id but different data
        sample_set = compute_sample_set_id(["s1"])
        register_run("shared_run", sample_set, state_dir=lab_a)
        register_run("shared_run", sample_set, state_dir=lab_b)

        # Change status in lab_a to create a conflict
        with connect(lab_a) as conn:
            conn.execute(
                "UPDATE runs SET status = 'completed' WHERE run_id = 'shared_run'",
            )
            conn.commit()

        # Export from lab A
        monkeypatch.setenv("NVD_STATE_DIR", str(lab_a))
        archive = tmp_path / "export.nvd-state"
        result = runner.invoke(app, ["state", "export", str(archive)])
        assert result.exit_code == 0

        # Try to merge into lab B — should fail
        monkeypatch.setenv("NVD_STATE_DIR", str(lab_b))
        result = runner.invoke(app, ["state", "import", str(archive), "--yes"])
        assert result.exit_code != 0

        # Lab B data should be unchanged
        with connect(lab_b) as conn:
            run = conn.execute(
                "SELECT status FROM runs WHERE run_id = 'shared_run'",
            ).fetchone()
        assert run["status"] == "running"


class TestExportCommand:
    """Tests for the nvd state export CLI command."""

    def test_json_outputs_manifest(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """--json outputs manifest without creating an archive."""
        _seed_lab_data(tmp_path, "run_001", ["s1"])
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))

        result = runner.invoke(app, ["state", "export", "--json"])
        assert result.exit_code == 0

        data = json.loads(result.stdout)
        assert "schema_version" in data
        assert "counts" in data
        assert data["counts"]["runs"] == 1

    def test_rejects_json_with_path(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Cannot specify both --json and a path."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        with connect(tmp_path):
            pass

        result = runner.invoke(
            app,
            ["state", "export", "--json", str(tmp_path / "out.nvd-state")],
        )
        assert result.exit_code != 0

    def test_requires_path_without_json(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Path is required when --json is not specified."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        with connect(tmp_path):
            pass

        result = runner.invoke(app, ["state", "export"])
        assert result.exit_code != 0

    def test_appends_extension(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Auto-appends .nvd-state extension if missing."""
        _seed_lab_data(tmp_path, "run_001", ["s1"])
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))

        result = runner.invoke(app, ["state", "export", str(tmp_path / "backup")])
        assert result.exit_code == 0
        assert (tmp_path / "backup.nvd-state").exists()

    def test_refuses_overwrite_without_force(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Refuses to overwrite an existing archive without --force."""
        _seed_lab_data(tmp_path, "run_001", ["s1"])
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))

        archive = tmp_path / "backup.nvd-state"
        archive.touch()

        result = runner.invoke(app, ["state", "export", str(archive)])
        assert result.exit_code != 0

    def test_force_overwrites(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """--force allows overwriting an existing archive."""
        _seed_lab_data(tmp_path, "run_001", ["s1"])
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))

        archive = tmp_path / "backup.nvd-state"
        archive.touch()

        result = runner.invoke(app, ["state", "export", str(archive), "--force"])
        assert result.exit_code == 0
        assert archive.stat().st_size > 0

    def test_archive_contains_expected_files(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Archive contains state.sqlite and manifest.json."""
        _seed_lab_data(tmp_path, "run_001", ["s1"])
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))

        archive = tmp_path / "backup.nvd-state"
        result = runner.invoke(app, ["state", "export", str(archive)])
        assert result.exit_code == 0

        with tarfile.open(archive, "r:gz") as tar:
            names = tar.getnames()
        assert "state.sqlite" in names
        assert "manifest.json" in names


class TestImportCommand:
    """Tests for the nvd state import CLI command."""

    def test_schema_mismatch_aborts(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Import aborts when archive schema version doesn't match."""
        # Create an archive with a wrong schema version
        with connect(tmp_path / "source"):
            pass

        staging = tmp_path / "staging"
        staging.mkdir()
        shutil.copy2(tmp_path / "source" / "state.sqlite", staging / "state.sqlite")
        manifest = {
            "schema_version": 9999,
            "exported_at": "2024-01-01",
            "nvd_version": __version__,
            "counts": {},
        }
        (staging / "manifest.json").write_text(json.dumps(manifest))

        archive = tmp_path / "bad_version.nvd-state"
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(staging / "state.sqlite", arcname="state.sqlite")
            tar.add(staging / "manifest.json", arcname="manifest.json")
        shutil.rmtree(staging)

        # Try to import
        target = tmp_path / "target"
        with connect(target):
            pass
        monkeypatch.setenv("NVD_STATE_DIR", str(target))

        result = runner.invoke(app, ["state", "import", str(archive), "--yes"])
        assert result.exit_code != 0

    def test_merge_requires_existing_db(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Merge import fails when local database does not exist."""
        source = tmp_path / "source"
        _seed_lab_data(source, "run_001", ["s1"])
        archive = _make_archive(tmp_path, source)

        # Point at a directory with no database
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        monkeypatch.setenv("NVD_STATE_DIR", str(empty_dir))

        result = runner.invoke(app, ["state", "import", str(archive), "--yes"])
        assert result.exit_code != 0

    def test_merge_json_output(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Merge import with --json produces structured output."""
        lab_a = tmp_path / "lab_a"
        lab_b = tmp_path / "lab_b"

        _seed_lab_data(lab_a, "run_a", ["s1"])
        archive = _make_archive(tmp_path, lab_a)

        _seed_lab_data(lab_b, "run_b", ["s2"])
        monkeypatch.setenv("NVD_STATE_DIR", str(lab_b))

        result = runner.invoke(
            app,
            ["state", "import", str(archive), "--yes", "--json"],
        )
        assert result.exit_code == 0

        data = json.loads(result.stdout)
        assert data["mode"] == "merge"
        assert data["success"] is True
        assert "imported" in data

    def test_replace_json_output(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Replace import with --json produces structured output."""
        lab_a = tmp_path / "lab_a"
        lab_b = tmp_path / "lab_b"

        _seed_lab_data(lab_a, "run_a", ["s1"])
        archive = _make_archive(tmp_path, lab_a)

        _seed_lab_data(lab_b, "run_b", ["s2"])
        monkeypatch.setenv("NVD_STATE_DIR", str(lab_b))

        result = runner.invoke(
            app,
            ["state", "import", str(archive), "--replace", "--yes", "--json"],
        )
        assert result.exit_code == 0

        data = json.loads(result.stdout)
        assert data["mode"] == "replace"
        assert data["success"] is True

    def test_import_nonexistent_archive(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Import fails cleanly when archive file does not exist."""
        monkeypatch.setenv("NVD_STATE_DIR", str(tmp_path))
        with connect(tmp_path):
            pass

        result = runner.invoke(
            app,
            ["state", "import", str(tmp_path / "ghost.nvd-state"), "--yes"],
        )
        assert result.exit_code != 0
