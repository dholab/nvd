"""Tests for py_nvd.state module."""

import tempfile
from pathlib import Path

import pytest

from py_nvd.state import (
    check_taxonomy_drift,
    check_upload,
    complete_sample,
    compute_sample_set_id,
    get_processed_sample,
    get_run,
    get_run_by_sample_set,
    get_sample_history,
    get_samples_for_run,
    get_taxonomy_version,
    get_upload,
    get_uploads_for_sample,
    hash_upload_content,
    list_runs,
    list_samples,
    list_uploads,
    record_taxonomy_version,
    record_upload,
    register_processed_sample,
    register_run,
    was_sample_ever_processed,
    was_sample_ever_uploaded,
)


class TestComputeSampleSetId:
    """Tests for compute_sample_set_id()."""

    def test_deterministic(self):
        """Same input produces same output."""
        samples = ["sample_a", "sample_b", "sample_c"]
        id1 = compute_sample_set_id(samples)
        id2 = compute_sample_set_id(samples)
        assert id1 == id2

    def test_order_independent(self):
        """Order of samples doesn't affect the ID."""
        id1 = compute_sample_set_id(["a", "b", "c"])
        id2 = compute_sample_set_id(["c", "a", "b"])
        id3 = compute_sample_set_id(["b", "c", "a"])
        assert id1 == id2 == id3

    def test_duplicates_ignored(self):
        """Duplicate sample IDs are deduplicated."""
        id1 = compute_sample_set_id(["a", "b", "c"])
        id2 = compute_sample_set_id(["a", "a", "b", "b", "c", "c"])
        assert id1 == id2

    def test_different_samples_different_id(self):
        """Different samples produce different IDs."""
        id1 = compute_sample_set_id(["a", "b", "c"])
        id2 = compute_sample_set_id(["a", "b", "d"])
        assert id1 != id2

    def test_returns_16_char_hex(self):
        """Returns a 16-character hex string."""
        result = compute_sample_set_id(["sample1", "sample2"])
        assert len(result) == 16
        assert all(c in "0123456789abcdef" for c in result)

    def test_empty_list(self):
        """Empty list produces a valid ID."""
        result = compute_sample_set_id([])
        assert len(result) == 16


class TestRegisterRun:
    """Tests for register_run() and related functions."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_register_new_run(self, temp_state_dir):
        """Registering a new run succeeds."""
        sample_set_id = compute_sample_set_id(["s1", "s2"])
        run = register_run("run_001", sample_set_id, state_dir=temp_state_dir)

        assert run is not None
        assert run.run_id == "run_001"
        assert run.sample_set_id == sample_set_id
        assert run.status == "running"
        assert run.experiment_id is None

    def test_register_with_experiment_id(self, temp_state_dir):
        """Registering with optional experiment_id stores it."""
        sample_set_id = compute_sample_set_id(["s1", "s2"])
        run = register_run(
            "run_001",
            sample_set_id,
            experiment_id=12345,
            state_dir=temp_state_dir,
        )

        assert run is not None
        assert run.experiment_id == 12345

    def test_duplicate_sample_set_returns_none(self, temp_state_dir):
        """Attempting to register same sample_set_id returns None."""
        sample_set_id = compute_sample_set_id(["s1", "s2"])

        run1 = register_run("run_001", sample_set_id, state_dir=temp_state_dir)
        run2 = register_run("run_002", sample_set_id, state_dir=temp_state_dir)

        assert run1 is not None
        assert run2 is None  # Duplicate sample set rejected

    def test_different_sample_sets_allowed(self, temp_state_dir):
        """Different sample sets can be registered."""
        set1 = compute_sample_set_id(["s1", "s2"])
        set2 = compute_sample_set_id(["s3", "s4"])

        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)

        assert run1 is not None
        assert run2 is not None

    def test_get_run(self, temp_state_dir):
        """get_run retrieves a registered run."""
        sample_set_id = compute_sample_set_id(["s1"])
        register_run("run_001", sample_set_id, state_dir=temp_state_dir)

        run = get_run("run_001", state_dir=temp_state_dir)
        assert run is not None
        assert run.run_id == "run_001"

    def test_get_run_not_found(self, temp_state_dir):
        """get_run returns None for unknown run_id."""
        run = get_run("nonexistent", state_dir=temp_state_dir)
        assert run is None

    def test_get_run_by_sample_set(self, temp_state_dir):
        """get_run_by_sample_set retrieves by sample_set_id."""
        sample_set_id = compute_sample_set_id(["s1", "s2"])
        register_run("run_001", sample_set_id, state_dir=temp_state_dir)

        run = get_run_by_sample_set(sample_set_id, state_dir=temp_state_dir)
        assert run is not None
        assert run.run_id == "run_001"

    def test_list_runs(self, temp_state_dir):
        """list_runs returns all registered runs."""
        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])

        register_run("run_001", set1, state_dir=temp_state_dir)
        register_run("run_002", set2, state_dir=temp_state_dir)

        runs = list_runs(state_dir=temp_state_dir)
        assert len(runs) == 2
        run_ids = {r.run_id for r in runs}
        assert run_ids == {"run_001", "run_002"}

    def test_list_runs_filter_by_status(self, temp_state_dir):
        """list_runs filters by status."""
        from py_nvd.state import complete_run

        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])
        set3 = compute_sample_set_id(["s3"])

        register_run("run_001", set1, state_dir=temp_state_dir)
        register_run("run_002", set2, state_dir=temp_state_dir)
        register_run("run_003", set3, state_dir=temp_state_dir)

        # Complete one, fail another, leave third running
        complete_run("run_001", status="completed", state_dir=temp_state_dir)
        complete_run("run_002", status="failed", state_dir=temp_state_dir)

        # Filter by status
        completed = list_runs(status="completed", state_dir=temp_state_dir)
        assert len(completed) == 1
        assert completed[0].run_id == "run_001"

        failed = list_runs(status="failed", state_dir=temp_state_dir)
        assert len(failed) == 1
        assert failed[0].run_id == "run_002"

        running = list_runs(status="running", state_dir=temp_state_dir)
        assert len(running) == 1
        assert running[0].run_id == "run_003"

    def test_list_runs_filter_by_limit(self, temp_state_dir):
        """list_runs respects limit parameter."""
        for i in range(5):
            set_id = compute_sample_set_id([f"s{i}"])
            register_run(f"run_{i:03d}", set_id, state_dir=temp_state_dir)

        runs = list_runs(limit=3, state_dir=temp_state_dir)
        assert len(runs) == 3

        runs = list_runs(limit=1, state_dir=temp_state_dir)
        assert len(runs) == 1

    def test_list_runs_filter_by_since_date(self, temp_state_dir):
        """list_runs filters by since with a date."""
        from datetime import date

        from py_nvd.db import connect

        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])

        # Insert runs with explicit timestamps for testing
        with connect(temp_state_dir) as conn:
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, started_at, status)
                   VALUES (?, ?, ?, ?)""",
                ("run_old", set1, "2024-01-15 10:00:00", "completed"),
            )
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, started_at, status)
                   VALUES (?, ?, ?, ?)""",
                ("run_new", set2, "2024-12-20 10:00:00", "completed"),
            )
            conn.commit()

        # Filter since 2024-06-01
        runs = list_runs(since=date(2024, 6, 1), state_dir=temp_state_dir)
        assert len(runs) == 1
        assert runs[0].run_id == "run_new"

        # Filter since 2024-01-01 (should get both)
        runs = list_runs(since=date(2024, 1, 1), state_dir=temp_state_dir)
        assert len(runs) == 2

    def test_list_runs_filter_by_since_datetime(self, temp_state_dir):
        """list_runs filters by since with a datetime."""
        from datetime import datetime

        from py_nvd.db import connect

        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])

        # Insert runs with explicit timestamps
        with connect(temp_state_dir) as conn:
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, started_at, status)
                   VALUES (?, ?, ?, ?)""",
                ("run_morning", set1, "2024-12-20 09:00:00", "completed"),
            )
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, started_at, status)
                   VALUES (?, ?, ?, ?)""",
                ("run_afternoon", set2, "2024-12-20 14:00:00", "completed"),
            )
            conn.commit()

        # Filter since noon
        runs = list_runs(
            since=datetime(2024, 12, 20, 12, 0, 0), state_dir=temp_state_dir
        )
        assert len(runs) == 1
        assert runs[0].run_id == "run_afternoon"

    def test_list_runs_combined_filters(self, temp_state_dir):
        """list_runs combines multiple filters."""
        from datetime import date

        from py_nvd.db import connect
        from py_nvd.state import complete_run

        # Create runs with different dates and statuses
        with connect(temp_state_dir) as conn:
            for i, (run_id, set_id, started, status) in enumerate(
                [
                    ("run_a", compute_sample_set_id(["a"]), "2024-01-10", "completed"),
                    ("run_b", compute_sample_set_id(["b"]), "2024-06-15", "completed"),
                    ("run_c", compute_sample_set_id(["c"]), "2024-06-20", "failed"),
                    ("run_d", compute_sample_set_id(["d"]), "2024-12-01", "completed"),
                    ("run_e", compute_sample_set_id(["e"]), "2024-12-15", "completed"),
                ]
            ):
                conn.execute(
                    """INSERT INTO runs (run_id, sample_set_id, started_at, status)
                       VALUES (?, ?, ?, ?)""",
                    (run_id, set_id, started, status),
                )
            conn.commit()

        # Completed runs since June, limit 2
        runs = list_runs(
            status="completed",
            since=date(2024, 6, 1),
            limit=2,
            state_dir=temp_state_dir,
        )
        assert len(runs) == 2
        # Should be most recent first: run_e, run_d
        assert runs[0].run_id == "run_e"
        assert runs[1].run_id == "run_d"


class TestProcessedSamples:
    """Tests for processed sample tracking functions."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def registered_run(self, temp_state_dir):
        """Register a run and return (run, sample_set_id, temp_state_dir)."""
        sample_set_id = compute_sample_set_id(["s1", "s2", "s3"])
        run = register_run("run_001", sample_set_id, state_dir=temp_state_dir)
        return run, sample_set_id, temp_state_dir

    def test_register_processed_sample(self, registered_run):
        """Registering a processed sample succeeds."""
        run, sample_set_id, state_dir = registered_run

        sample = register_processed_sample(
            sample_id="s1",
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            blast_db_version="core-nt_2025-01-01",
            stat_db_version="v1.0",
            state_dir=state_dir,
        )

        assert sample is not None
        assert sample.sample_id == "s1"
        assert sample.sample_set_id == sample_set_id
        assert sample.run_id == "run_001"
        assert sample.status == "processing"
        assert sample.blast_db_version == "core-nt_2025-01-01"
        assert sample.stat_db_version == "v1.0"

    def test_get_processed_sample(self, registered_run):
        """get_processed_sample retrieves a registered sample."""
        run, sample_set_id, state_dir = registered_run

        register_processed_sample(
            sample_id="s1",
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            state_dir=state_dir,
        )

        sample = get_processed_sample("s1", sample_set_id, state_dir=state_dir)
        assert sample is not None
        assert sample.sample_id == "s1"

    def test_get_processed_sample_not_found(self, temp_state_dir):
        """get_processed_sample returns None for unknown sample."""
        sample = get_processed_sample(
            "nonexistent",
            "fake_set",
            state_dir=temp_state_dir,
        )
        assert sample is None

    def test_get_samples_for_run(self, registered_run):
        """get_samples_for_run returns all samples for a run."""
        run, sample_set_id, state_dir = registered_run

        for sid in ["s1", "s2", "s3"]:
            register_processed_sample(
                sample_id=sid,
                sample_set_id=sample_set_id,
                run_id=run.run_id,
                state_dir=state_dir,
            )

        samples = get_samples_for_run(run.run_id, state_dir=state_dir)
        assert len(samples) == 3
        sample_ids = {s.sample_id for s in samples}
        assert sample_ids == {"s1", "s2", "s3"}

    def test_was_sample_ever_processed(self, registered_run):
        """was_sample_ever_processed returns True for processed samples."""
        run, sample_set_id, state_dir = registered_run

        # Not processed yet
        assert was_sample_ever_processed("s1", state_dir=state_dir) is False

        # Process it
        register_processed_sample(
            sample_id="s1",
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            state_dir=state_dir,
        )

        # Now it's processed
        assert was_sample_ever_processed("s1", state_dir=state_dir) is True
        # Other samples still not processed
        assert was_sample_ever_processed("s2", state_dir=state_dir) is False

    def test_complete_sample(self, registered_run):
        """complete_sample updates the sample status."""
        run, sample_set_id, state_dir = registered_run

        register_processed_sample(
            sample_id="s1",
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            state_dir=state_dir,
        )

        # Initially processing
        sample = get_processed_sample("s1", sample_set_id, state_dir=state_dir)
        assert sample is not None
        assert sample.status == "processing"

        # Mark as completed
        complete_sample("s1", sample_set_id, status="completed", state_dir=state_dir)

        sample = get_processed_sample("s1", sample_set_id, state_dir=state_dir)
        assert sample is not None
        assert sample.status == "completed"

    def test_complete_sample_failed(self, registered_run):
        """complete_sample can mark a sample as failed."""
        run, sample_set_id, state_dir = registered_run

        register_processed_sample(
            sample_id="s1",
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            state_dir=state_dir,
        )

        complete_sample("s1", sample_set_id, status="failed", state_dir=state_dir)

        sample = get_processed_sample("s1", sample_set_id, state_dir=state_dir)
        assert sample is not None
        assert sample.status == "failed"

    def test_list_samples(self, temp_state_dir):
        """list_samples returns all processed samples."""
        set1 = compute_sample_set_id(["s1", "s2"])
        set2 = compute_sample_set_id(["s3"])

        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        register_processed_sample("s1", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s3", set2, run2.run_id, state_dir=temp_state_dir)

        samples = list_samples(state_dir=temp_state_dir)
        assert len(samples) == 3
        sample_ids = {s.sample_id for s in samples}
        assert sample_ids == {"s1", "s2", "s3"}

    def test_list_samples_filter_by_run(self, temp_state_dir):
        """list_samples filters by run_id."""
        set1 = compute_sample_set_id(["s1", "s2"])
        set2 = compute_sample_set_id(["s3"])

        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        register_processed_sample("s1", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s3", set2, run2.run_id, state_dir=temp_state_dir)

        samples = list_samples(run_id="run_001", state_dir=temp_state_dir)
        assert len(samples) == 2
        sample_ids = {s.sample_id for s in samples}
        assert sample_ids == {"s1", "s2"}

    def test_list_samples_filter_by_status(self, temp_state_dir):
        """list_samples filters by status."""
        set1 = compute_sample_set_id(["s1", "s2", "s3"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run.run_id, state_dir=temp_state_dir)
        register_processed_sample("s3", set1, run.run_id, state_dir=temp_state_dir)

        complete_sample("s1", set1, status="completed", state_dir=temp_state_dir)
        complete_sample("s2", set1, status="failed", state_dir=temp_state_dir)
        # s3 remains "processing"

        completed = list_samples(status="completed", state_dir=temp_state_dir)
        assert len(completed) == 1
        assert completed[0].sample_id == "s1"

        failed = list_samples(status="failed", state_dir=temp_state_dir)
        assert len(failed) == 1
        assert failed[0].sample_id == "s2"

        processing = list_samples(status="processing", state_dir=temp_state_dir)
        assert len(processing) == 1
        assert processing[0].sample_id == "s3"

    def test_list_samples_filter_by_limit(self, temp_state_dir):
        """list_samples respects limit parameter."""
        set1 = compute_sample_set_id(["s1", "s2", "s3", "s4", "s5"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        for sid in ["s1", "s2", "s3", "s4", "s5"]:
            register_processed_sample(sid, set1, run.run_id, state_dir=temp_state_dir)

        samples = list_samples(limit=3, state_dir=temp_state_dir)
        assert len(samples) == 3

    def test_list_samples_combined_filters(self, temp_state_dir):
        """list_samples combines multiple filters."""
        set1 = compute_sample_set_id(["s1", "s2"])
        set2 = compute_sample_set_id(["s3", "s4"])

        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        register_processed_sample("s1", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s3", set2, run2.run_id, state_dir=temp_state_dir)
        register_processed_sample("s4", set2, run2.run_id, state_dir=temp_state_dir)

        complete_sample("s1", set1, status="completed", state_dir=temp_state_dir)
        complete_sample("s3", set2, status="completed", state_dir=temp_state_dir)

        # Completed samples from run_001
        samples = list_samples(
            run_id="run_001", status="completed", state_dir=temp_state_dir
        )
        assert len(samples) == 1
        assert samples[0].sample_id == "s1"

    def test_get_sample_history(self, temp_state_dir):
        """get_sample_history returns all processing records for a sample across runs."""
        # Sample s1 appears in two different runs
        set1 = compute_sample_set_id(["s1", "s2"])
        set2 = compute_sample_set_id(["s1", "s3"])  # s1 is in both

        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        register_processed_sample("s1", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s1", set2, run2.run_id, state_dir=temp_state_dir)
        register_processed_sample("s3", set2, run2.run_id, state_dir=temp_state_dir)

        # s1 should have 2 records
        history = get_sample_history("s1", state_dir=temp_state_dir)
        assert len(history) == 2
        run_ids = {h.run_id for h in history}
        assert run_ids == {"run_001", "run_002"}

        # s2 should have 1 record
        history = get_sample_history("s2", state_dir=temp_state_dir)
        assert len(history) == 1
        assert history[0].run_id == "run_001"

    def test_get_sample_history_not_found(self, temp_state_dir):
        """get_sample_history returns empty list for unknown sample."""
        history = get_sample_history("nonexistent", state_dir=temp_state_dir)
        assert history == []


class TestUploads:
    """Tests for upload tracking functions."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def processed_sample(self, temp_state_dir):
        """Register a run and processed sample, return context."""
        sample_set_id = compute_sample_set_id(["s1"])
        run = register_run("run_001", sample_set_id, state_dir=temp_state_dir)
        assert run is not None
        register_processed_sample(
            sample_id="s1",
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            state_dir=temp_state_dir,
        )
        return {
            "sample_id": "s1",
            "sample_set_id": sample_set_id,
            "run_id": run.run_id,
            "state_dir": temp_state_dir,
        }

    def test_record_upload(self, processed_sample):
        """record_upload creates an upload record."""
        ctx = processed_sample

        upload = record_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            target_metadata={"experiment_id": 42, "labkey_row_id": 100},
            state_dir=ctx["state_dir"],
        )

        assert upload is not None
        assert upload.sample_id == "s1"
        assert upload.upload_type == "blast"
        assert upload.upload_target == "labkey"
        assert upload.content_hash == "abc123"
        assert upload.target_metadata == {"experiment_id": 42, "labkey_row_id": 100}

    def test_get_upload(self, processed_sample):
        """get_upload retrieves a specific upload."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            state_dir=ctx["state_dir"],
        )

        upload = get_upload(
            ctx["sample_id"],
            ctx["sample_set_id"],
            "blast",
            "labkey",
            state_dir=ctx["state_dir"],
        )
        assert upload is not None
        assert upload.content_hash == "abc123"

    def test_get_upload_not_found(self, temp_state_dir):
        """get_upload returns None for unknown upload."""
        upload = get_upload(
            "s1",
            "fake_set",
            "blast",
            "labkey",
            state_dir=temp_state_dir,
        )
        assert upload is None

    def test_check_upload_not_uploaded(self, processed_sample):
        """check_upload returns 'not_uploaded' for new uploads."""
        ctx = processed_sample

        result = check_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            state_dir=ctx["state_dir"],
        )
        assert result == "not_uploaded"

    def test_check_upload_already_uploaded(self, processed_sample):
        """check_upload returns 'already_uploaded' for same content."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            state_dir=ctx["state_dir"],
        )

        result = check_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",  # Same hash
            state_dir=ctx["state_dir"],
        )
        assert result == "already_uploaded"

    def test_check_upload_content_changed(self, processed_sample):
        """check_upload returns 'content_changed' for different content."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            state_dir=ctx["state_dir"],
        )

        result = check_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="different_hash",  # Different hash
            state_dir=ctx["state_dir"],
        )
        assert result == "content_changed"

    def test_get_uploads_for_sample(self, processed_sample):
        """get_uploads_for_sample returns all uploads for a sample."""
        ctx = processed_sample

        # Upload to multiple targets
        record_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="hash1",
            state_dir=ctx["state_dir"],
        )
        record_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast_fasta",
            upload_target="labkey",
            content_hash="hash2",
            state_dir=ctx["state_dir"],
        )

        uploads = get_uploads_for_sample(ctx["sample_id"], state_dir=ctx["state_dir"])
        assert len(uploads) == 2
        upload_types = {u.upload_type for u in uploads}
        assert upload_types == {"blast", "blast_fasta"}

    def test_was_sample_ever_uploaded_false(self, processed_sample):
        """was_sample_ever_uploaded returns False for un-uploaded samples."""
        ctx = processed_sample
        assert (
            was_sample_ever_uploaded(ctx["sample_id"], state_dir=ctx["state_dir"])
            is False
        )

    def test_was_sample_ever_uploaded_true(self, processed_sample):
        """was_sample_ever_uploaded returns True for uploaded samples."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="hash1",
            state_dir=ctx["state_dir"],
        )

        assert (
            was_sample_ever_uploaded(ctx["sample_id"], state_dir=ctx["state_dir"])
            is True
        )

    def test_was_sample_ever_uploaded_with_filters(self, processed_sample):
        """was_sample_ever_uploaded respects upload_type and upload_target filters."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx["sample_id"],
            sample_set_id=ctx["sample_set_id"],
            upload_type="blast",
            upload_target="labkey",
            content_hash="hash1",
            state_dir=ctx["state_dir"],
        )

        # Uploaded to labkey
        assert (
            was_sample_ever_uploaded(
                ctx["sample_id"],
                upload_target="labkey",
                state_dir=ctx["state_dir"],
            )
            is True
        )

        # Not uploaded to globus
        assert (
            was_sample_ever_uploaded(
                ctx["sample_id"],
                upload_target="globus",
                state_dir=ctx["state_dir"],
            )
            is False
        )

        # Uploaded blast type
        assert (
            was_sample_ever_uploaded(
                ctx["sample_id"],
                upload_type="blast",
                state_dir=ctx["state_dir"],
            )
            is True
        )

        # Not uploaded gottcha2 type
        assert (
            was_sample_ever_uploaded(
                ctx["sample_id"],
                upload_type="gottcha2",
                state_dir=ctx["state_dir"],
            )
            is False
        )

    def test_cross_run_deduplication(self, temp_state_dir):
        """was_sample_ever_uploaded detects uploads across different runs."""
        # First run
        set1 = compute_sample_set_id(["s1", "s2"])
        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run1 is not None
        register_processed_sample("s1", set1, run1.run_id, state_dir=temp_state_dir)
        record_upload(
            sample_id="s1",
            sample_set_id=set1,
            upload_type="blast",
            upload_target="labkey",
            content_hash="hash1",
            state_dir=temp_state_dir,
        )

        # Second run with overlapping sample
        set2 = compute_sample_set_id(["s1", "s3"])  # s1 is in both runs
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run2 is not None
        register_processed_sample("s1", set2, run2.run_id, state_dir=temp_state_dir)

        # s1 was uploaded in run1, should be detected in run2 context
        assert (
            was_sample_ever_uploaded(
                "s1",
                upload_type="blast",
                upload_target="labkey",
                state_dir=temp_state_dir,
            )
            is True
        )

        # s3 was never uploaded
        assert (
            was_sample_ever_uploaded(
                "s3",
                upload_target="labkey",
                state_dir=temp_state_dir,
            )
            is False
        )

    def test_list_uploads(self, temp_state_dir):
        """list_uploads returns all uploads."""
        # Set up two samples with uploads
        set1 = compute_sample_set_id(["s1", "s2"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run.run_id, state_dir=temp_state_dir)

        record_upload("s1", set1, "blast", "labkey", "hash1", state_dir=temp_state_dir)
        record_upload(
            "s1", set1, "blast_fasta", "labkey", "hash2", state_dir=temp_state_dir
        )
        record_upload("s2", set1, "blast", "labkey", "hash3", state_dir=temp_state_dir)

        uploads = list_uploads(state_dir=temp_state_dir)
        assert len(uploads) == 3

    def test_list_uploads_filter_by_sample(self, temp_state_dir):
        """list_uploads filters by sample_id."""
        set1 = compute_sample_set_id(["s1", "s2"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run.run_id, state_dir=temp_state_dir)

        record_upload("s1", set1, "blast", "labkey", "hash1", state_dir=temp_state_dir)
        record_upload("s2", set1, "blast", "labkey", "hash2", state_dir=temp_state_dir)

        uploads = list_uploads(sample_id="s1", state_dir=temp_state_dir)
        assert len(uploads) == 1
        assert uploads[0].sample_id == "s1"

    def test_list_uploads_filter_by_type(self, temp_state_dir):
        """list_uploads filters by upload_type."""
        set1 = compute_sample_set_id(["s1"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)

        record_upload("s1", set1, "blast", "labkey", "hash1", state_dir=temp_state_dir)
        record_upload(
            "s1", set1, "blast_fasta", "labkey", "hash2", state_dir=temp_state_dir
        )
        record_upload(
            "s1", set1, "gottcha2", "labkey", "hash3", state_dir=temp_state_dir
        )

        uploads = list_uploads(upload_type="blast", state_dir=temp_state_dir)
        assert len(uploads) == 1
        assert uploads[0].upload_type == "blast"

    def test_list_uploads_filter_by_target(self, temp_state_dir):
        """list_uploads filters by upload_target."""
        set1 = compute_sample_set_id(["s1"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)

        record_upload("s1", set1, "blast", "labkey", "hash1", state_dir=temp_state_dir)
        record_upload(
            "s1", set1, "blast_fasta", "local", "hash2", state_dir=temp_state_dir
        )

        uploads = list_uploads(upload_target="labkey", state_dir=temp_state_dir)
        assert len(uploads) == 1
        assert uploads[0].upload_target == "labkey"

    def test_list_uploads_with_limit(self, temp_state_dir):
        """list_uploads respects limit parameter."""
        set1 = compute_sample_set_id(["s1"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)

        record_upload("s1", set1, "blast", "labkey", "hash1", state_dir=temp_state_dir)
        record_upload(
            "s1", set1, "blast_fasta", "labkey", "hash2", state_dir=temp_state_dir
        )
        record_upload(
            "s1", set1, "gottcha2", "labkey", "hash3", state_dir=temp_state_dir
        )

        uploads = list_uploads(limit=2, state_dir=temp_state_dir)
        assert len(uploads) == 2


class TestHashUploadContent:
    """Tests for hash_upload_content()."""

    def test_deterministic(self):
        """Same input produces same hash."""
        data = {"key": "value", "number": 42}
        hash1 = hash_upload_content(data)
        hash2 = hash_upload_content(data)
        assert hash1 == hash2

    def test_order_independent_for_dicts(self):
        """Dict key order doesn't affect hash (sorted keys)."""
        data1 = {"b": 2, "a": 1}
        data2 = {"a": 1, "b": 2}
        assert hash_upload_content(data1) == hash_upload_content(data2)

    def test_different_data_different_hash(self):
        """Different data produces different hashes."""
        hash1 = hash_upload_content({"key": "value1"})
        hash2 = hash_upload_content({"key": "value2"})
        assert hash1 != hash2

    def test_handles_lists(self):
        """Can hash lists of dicts."""
        data = [{"id": 1}, {"id": 2}]
        result = hash_upload_content(data)
        assert len(result) == 64  # SHA256 hex length


class TestCrossRunDeduplicationE2E:
    """
    End-to-end tests for cross-run upload deduplication.

    These tests simulate the workflow used by labkey_upload_blast_results.py:
    1. First run processes and uploads samples A, B, C
    2. Second run has overlapping samples B, C, D
    3. Only sample D should be uploaded (B, C already uploaded)
    """

    @pytest.fixture
    def temp_state_dir(self, tmp_path):
        """Create a temporary state directory."""
        state_dir = tmp_path / "state"
        state_dir.mkdir()
        return state_dir

    def test_upload_workflow_first_run(self, temp_state_dir):
        """First run: all samples are uploaded."""
        samples = ["sample_A", "sample_B", "sample_C"]
        sample_set_id = compute_sample_set_id(samples)
        run = register_run("run_001", sample_set_id, state_dir=temp_state_dir)
        assert run is not None

        # Register all samples as processed
        for sample_id in samples:
            register_processed_sample(
                sample_id=sample_id,
                sample_set_id=sample_set_id,
                run_id=run.run_id,
                state_dir=temp_state_dir,
            )

        # Simulate upload script's deduplication check
        samples_to_upload = []
        samples_to_skip = []
        for sample_id in samples:
            if was_sample_ever_uploaded(
                sample_id,
                upload_type="blast",
                upload_target="labkey",
                state_dir=temp_state_dir,
            ):
                samples_to_skip.append(sample_id)
            else:
                samples_to_upload.append(sample_id)

        # First run: nothing should be skipped
        assert samples_to_skip == []
        assert set(samples_to_upload) == set(samples)

        # Record uploads
        for sample_id in samples_to_upload:
            record_upload(
                sample_id=sample_id,
                sample_set_id=sample_set_id,
                upload_type="blast",
                upload_target="labkey",
                content_hash=f"hash_{sample_id}",
                state_dir=temp_state_dir,
            )

    def test_upload_workflow_second_run_with_overlap(self, temp_state_dir):
        """Second run: overlapping samples are skipped."""
        # First run: upload A, B, C
        samples_run1 = ["sample_A", "sample_B", "sample_C"]
        set1 = compute_sample_set_id(samples_run1)
        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run1 is not None

        for sample_id in samples_run1:
            register_processed_sample(
                sample_id=sample_id,
                sample_set_id=set1,
                run_id=run1.run_id,
                state_dir=temp_state_dir,
            )
            record_upload(
                sample_id=sample_id,
                sample_set_id=set1,
                upload_type="blast",
                upload_target="labkey",
                content_hash=f"hash_{sample_id}",
                state_dir=temp_state_dir,
            )

        # Second run: B, C, D (B and C overlap with run1)
        samples_run2 = ["sample_B", "sample_C", "sample_D"]
        set2 = compute_sample_set_id(samples_run2)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run2 is not None

        for sample_id in samples_run2:
            register_processed_sample(
                sample_id=sample_id,
                sample_set_id=set2,
                run_id=run2.run_id,
                state_dir=temp_state_dir,
            )

        # Simulate upload script's deduplication check
        samples_to_upload = []
        samples_to_skip = []
        for sample_id in samples_run2:
            if was_sample_ever_uploaded(
                sample_id,
                upload_type="blast",
                upload_target="labkey",
                state_dir=temp_state_dir,
            ):
                samples_to_skip.append(sample_id)
            else:
                samples_to_upload.append(sample_id)

        # B and C should be skipped (already uploaded in run1)
        assert set(samples_to_skip) == {"sample_B", "sample_C"}
        # Only D should be uploaded
        assert samples_to_upload == ["sample_D"]

    def test_dedup_respects_upload_type(self, temp_state_dir):
        """Deduplication is per upload_type - blast and gottcha2 are independent."""
        sample_id = "sample_X"
        sample_set_id = compute_sample_set_id([sample_id])
        run = register_run("run_001", sample_set_id, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample(
            sample_id=sample_id,
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            state_dir=temp_state_dir,
        )

        # Upload only BLAST results
        record_upload(
            sample_id=sample_id,
            sample_set_id=sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="blast_hash",
            state_dir=temp_state_dir,
        )

        # BLAST should be detected as uploaded
        assert was_sample_ever_uploaded(
            sample_id,
            upload_type="blast",
            upload_target="labkey",
            state_dir=temp_state_dir,
        )

        # GOTTCHA2 should NOT be detected (different upload_type)
        assert not was_sample_ever_uploaded(
            sample_id,
            upload_type="gottcha2",
            upload_target="labkey",
            state_dir=temp_state_dir,
        )

    def test_dedup_respects_upload_target(self, temp_state_dir):
        """Deduplication is per upload_target - labkey and local are independent."""
        sample_id = "sample_Y"
        sample_set_id = compute_sample_set_id([sample_id])
        run = register_run("run_001", sample_set_id, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample(
            sample_id=sample_id,
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            state_dir=temp_state_dir,
        )

        # Upload to LabKey only
        record_upload(
            sample_id=sample_id,
            sample_set_id=sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="labkey_hash",
            state_dir=temp_state_dir,
        )

        # LabKey should be detected
        assert was_sample_ever_uploaded(
            sample_id,
            upload_type="blast",
            upload_target="labkey",
            state_dir=temp_state_dir,
        )

        # Local should NOT be detected (different target)
        assert not was_sample_ever_uploaded(
            sample_id,
            upload_type="blast",
            upload_target="local",
            state_dir=temp_state_dir,
        )


class TestTaxonomyDrift:
    """Tests for taxonomy version tracking and drift detection."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_record_and_get_taxonomy_version(self, temp_state_dir):
        """Can record and retrieve taxonomy version for a run."""
        # Register a run first (required for foreign key)
        sample_set_id = compute_sample_set_id(["s1"])
        register_run("run_001", sample_set_id, state_dir=temp_state_dir)

        # Record taxonomy version
        record_taxonomy_version(
            run_id="run_001",
            file_hash="abc123def456",
            ncbi_rebuild_timestamp="2024-01-15",
            file_size=500_000_000,
            downloaded_at="2024-01-20T10:30:00",
            state_dir=temp_state_dir,
        )

        # Retrieve it
        version = get_taxonomy_version("run_001", state_dir=temp_state_dir)
        assert version is not None
        assert version.run_id == "run_001"
        assert version.file_hash == "abc123def456"
        assert version.ncbi_rebuild_timestamp == "2024-01-15"
        assert version.file_size == 500_000_000
        assert version.downloaded_at == "2024-01-20T10:30:00"

    def test_get_taxonomy_version_not_found(self, temp_state_dir):
        """Returns None for non-existent run."""
        version = get_taxonomy_version("nonexistent_run", state_dir=temp_state_dir)
        assert version is None

    def test_same_version_no_drift(self, temp_state_dir):
        """Same taxonomy hash means no drift."""
        # Register two runs
        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])
        register_run("run_001", set1, state_dir=temp_state_dir)
        register_run("run_002", set2, state_dir=temp_state_dir)

        # Record same taxonomy version for both
        record_taxonomy_version(
            run_id="run_001",
            file_hash="same_hash_abc123",
            ncbi_rebuild_timestamp="2024-01-15",
            file_size=500_000_000,
            downloaded_at="2024-01-20T10:30:00",
            state_dir=temp_state_dir,
        )
        record_taxonomy_version(
            run_id="run_002",
            file_hash="same_hash_abc123",  # Same hash
            ncbi_rebuild_timestamp="2024-01-15",
            file_size=500_000_000,
            downloaded_at="2024-01-25T14:00:00",  # Different download time is OK
            state_dir=temp_state_dir,
        )

        # No drift detected
        assert check_taxonomy_drift("run_001", "run_002", temp_state_dir) is False

    def test_different_version_drift_detected(self, temp_state_dir):
        """Different taxonomy hash means drift detected."""
        # Register two runs
        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])
        register_run("run_001", set1, state_dir=temp_state_dir)
        register_run("run_002", set2, state_dir=temp_state_dir)

        # Record different taxonomy versions
        record_taxonomy_version(
            run_id="run_001",
            file_hash="old_taxonomy_hash",
            ncbi_rebuild_timestamp="2024-01-15",
            file_size=500_000_000,
            downloaded_at="2024-01-20T10:30:00",
            state_dir=temp_state_dir,
        )
        record_taxonomy_version(
            run_id="run_002",
            file_hash="new_taxonomy_hash",  # Different hash
            ncbi_rebuild_timestamp="2024-02-15",
            file_size=510_000_000,
            downloaded_at="2024-02-20T10:30:00",
            state_dir=temp_state_dir,
        )

        # Drift detected
        assert check_taxonomy_drift("run_001", "run_002", temp_state_dir) is True

    def test_missing_version_assumes_drift(self, temp_state_dir):
        """Missing version info conservatively assumes drift."""
        # Register only one run with taxonomy
        set1 = compute_sample_set_id(["s1"])
        register_run("run_001", set1, state_dir=temp_state_dir)

        record_taxonomy_version(
            run_id="run_001",
            file_hash="abc123",
            ncbi_rebuild_timestamp="2024-01-15",
            file_size=500_000_000,
            downloaded_at="2024-01-20T10:30:00",
            state_dir=temp_state_dir,
        )

        # run_002 doesn't exist - should assume drift (conservative)
        assert check_taxonomy_drift("run_001", "run_002", temp_state_dir) is True

    def test_both_missing_assumes_drift(self, temp_state_dir):
        """Both versions missing conservatively assumes drift."""
        # Neither run has taxonomy recorded
        assert check_taxonomy_drift("run_001", "run_002", temp_state_dir) is True

    def test_drift_check_is_symmetric(self, temp_state_dir):
        """Drift check gives same result regardless of argument order."""
        # Register two runs with different taxonomy
        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])
        register_run("run_001", set1, state_dir=temp_state_dir)
        register_run("run_002", set2, state_dir=temp_state_dir)

        record_taxonomy_version(
            run_id="run_001",
            file_hash="hash_a",
            ncbi_rebuild_timestamp=None,
            file_size=500_000_000,
            downloaded_at="2024-01-20",
            state_dir=temp_state_dir,
        )
        record_taxonomy_version(
            run_id="run_002",
            file_hash="hash_b",
            ncbi_rebuild_timestamp=None,
            file_size=500_000_000,
            downloaded_at="2024-02-20",
            state_dir=temp_state_dir,
        )

        # Order shouldn't matter
        assert check_taxonomy_drift(
            "run_001", "run_002", temp_state_dir
        ) == check_taxonomy_drift("run_002", "run_001", temp_state_dir)


class TestSchemaMigrationSafety:
    """Tests for schema migration safety in db.py."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_new_database_created_without_prompts(self, temp_state_dir):
        """Creating a new database doesn't require any prompts."""
        from py_nvd.db import connect

        db_path = temp_state_dir / "state.sqlite"
        assert not db_path.exists()

        with connect(temp_state_dir) as conn:
            # Should succeed without any prompts
            version = conn.execute("PRAGMA user_version").fetchone()[0]
            assert version == 2  # EXPECTED_VERSION

        assert db_path.exists()

    def test_matching_schema_version_proceeds_normally(self, temp_state_dir):
        """Database with matching schema version works without prompts."""
        from py_nvd.db import connect

        # Create database with correct schema
        with connect(temp_state_dir) as conn:
            conn.execute(
                "INSERT INTO runs (run_id, sample_set_id, started_at, status) "
                "VALUES ('test_run', 'abc123', datetime('now'), 'running')"
            )
            conn.commit()

        # Second connection should work fine
        with connect(temp_state_dir) as conn:
            row = conn.execute(
                "SELECT * FROM runs WHERE run_id = 'test_run'"
            ).fetchone()
            assert row is not None

    def test_schema_mismatch_raises_error_in_non_interactive(self, temp_state_dir):
        """Schema mismatch raises SchemaMismatchError when not interactive."""
        import sqlite3

        from py_nvd.db import EXPECTED_VERSION, SchemaMismatchError, connect

        db_path = temp_state_dir / "state.sqlite"

        # Create a database with wrong schema version
        conn = sqlite3.connect(db_path)
        conn.execute(f"PRAGMA user_version = {EXPECTED_VERSION + 1}")
        conn.execute("CREATE TABLE dummy (id INTEGER)")
        conn.commit()
        conn.close()

        # In non-interactive mode (tests), should raise error
        with pytest.raises(SchemaMismatchError) as exc_info:
            with connect(temp_state_dir):
                pass

        assert exc_info.value.current_version == EXPECTED_VERSION + 1
        assert exc_info.value.expected_version == EXPECTED_VERSION
        assert exc_info.value.db_path == db_path

    def test_allow_destructive_update_creates_backup(self, temp_state_dir):
        """allow_destructive_update=True creates backup before deletion."""
        import sqlite3

        from py_nvd.db import EXPECTED_VERSION, connect

        db_path = temp_state_dir / "state.sqlite"

        # Create a database with wrong schema version and some data
        conn = sqlite3.connect(db_path)
        conn.execute(f"PRAGMA user_version = {EXPECTED_VERSION - 1}")
        conn.execute("CREATE TABLE old_data (value TEXT)")
        conn.execute("INSERT INTO old_data VALUES ('important')")
        conn.commit()
        conn.close()

        original_size = db_path.stat().st_size

        # Connect with allow_destructive_update=True
        with connect(temp_state_dir, allow_destructive_update=True) as conn:
            # Should have new schema
            version = conn.execute("PRAGMA user_version").fetchone()[0]
            assert version == EXPECTED_VERSION

        # Backup should exist
        backups = list(temp_state_dir.glob("state.backup_*.sqlite"))
        assert len(backups) == 1

        # Backup should contain the old data
        backup_conn = sqlite3.connect(backups[0])
        old_version = backup_conn.execute("PRAGMA user_version").fetchone()[0]
        assert old_version == EXPECTED_VERSION - 1
        row = backup_conn.execute("SELECT value FROM old_data").fetchone()
        assert row[0] == "important"
        backup_conn.close()

    def test_allow_destructive_update_preserves_data_in_backup(self, temp_state_dir):
        """Backup contains all original data before schema migration."""
        import sqlite3

        from py_nvd.db import EXPECTED_VERSION, connect

        db_path = temp_state_dir / "state.sqlite"

        # Create database with old schema and multiple tables of data
        conn = sqlite3.connect(db_path)
        conn.execute(f"PRAGMA user_version = {EXPECTED_VERSION - 1}")
        conn.execute("CREATE TABLE runs (run_id TEXT PRIMARY KEY, data TEXT)")
        conn.execute("CREATE TABLE samples (sample_id TEXT, run_id TEXT)")
        for i in range(10):
            conn.execute(f"INSERT INTO runs VALUES ('run_{i}', 'data_{i}')")
            conn.execute(f"INSERT INTO samples VALUES ('sample_{i}', 'run_{i}')")
        conn.commit()
        conn.close()

        # Migrate
        with connect(temp_state_dir, allow_destructive_update=True):
            pass

        # Verify backup has all data
        backups = list(temp_state_dir.glob("state.backup_*.sqlite"))
        backup_conn = sqlite3.connect(backups[0])
        run_count = backup_conn.execute("SELECT COUNT(*) FROM runs").fetchone()[0]
        sample_count = backup_conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        assert run_count == 10
        assert sample_count == 10
        backup_conn.close()

    def test_schema_mismatch_error_message_is_helpful(self, temp_state_dir):
        """SchemaMismatchError message tells user what to do."""
        import sqlite3

        from py_nvd.db import EXPECTED_VERSION, SchemaMismatchError, connect

        db_path = temp_state_dir / "state.sqlite"

        # Create database with wrong version
        conn = sqlite3.connect(db_path)
        conn.execute(f"PRAGMA user_version = 999")
        conn.close()

        with pytest.raises(SchemaMismatchError) as exc_info:
            with connect(temp_state_dir):
                pass

        error_msg = str(exc_info.value)
        assert "999" in error_msg  # Current version
        assert str(EXPECTED_VERSION) in error_msg  # Expected version
        assert "--update-db-destructive" in error_msg  # CLI flag
        assert str(db_path) in error_msg  # Path to database


class TestDatabasePathLookup:
    """Tests for get_database_by_path() and get_databases_by_path()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def temp_db_path(self, temp_state_dir):
        """Create a temporary 'database' directory to register."""
        db_dir = temp_state_dir / "fake_blast_db"
        db_dir.mkdir()
        return db_dir

    def test_get_database_by_path_not_found(self, temp_state_dir, temp_db_path):
        """Looking up unregistered path returns (None, None)."""
        from py_nvd.state import get_database_by_path

        db, warning = get_database_by_path("blast", temp_db_path, temp_state_dir)

        assert db is None
        assert warning is None

    def test_get_database_by_path_single_match(self, temp_state_dir, temp_db_path):
        """Looking up registered path returns (database, None)."""
        from py_nvd.state import get_database_by_path, register_database

        register_database("blast", "v1.0", str(temp_db_path), state_dir=temp_state_dir)

        db, warning = get_database_by_path("blast", temp_db_path, temp_state_dir)

        assert db is not None
        assert db.version == "v1.0"
        assert db.db_type == "blast"
        assert warning is None

    def test_get_database_by_path_multiple_versions_returns_newest(
        self, temp_state_dir, temp_db_path
    ):
        """When multiple versions share a path, returns newest with warning."""
        import time

        from py_nvd.state import get_database_by_path, register_database

        # Register first version
        register_database("blast", "v1.0", str(temp_db_path), state_dir=temp_state_dir)
        # Small delay to ensure different timestamps
        time.sleep(0.1)
        # Register second version at same path
        register_database("blast", "v2.0", str(temp_db_path), state_dir=temp_state_dir)

        db, warning = get_database_by_path("blast", temp_db_path, temp_state_dir)

        assert db is not None
        assert db.version == "v2.0"  # Newest
        assert warning is not None
        assert "v1.0" in warning
        assert "v2.0" in warning
        assert "Multiple versions" in warning

    def test_get_database_by_path_canonicalizes_path(
        self, temp_state_dir, temp_db_path
    ):
        """Path is canonicalized before lookup."""
        from py_nvd.state import get_database_by_path, register_database

        # Register with absolute path
        register_database(
            "blast", "v1.0", str(temp_db_path.resolve()), state_dir=temp_state_dir
        )

        # Look up with relative path (from temp_state_dir)
        relative_path = temp_db_path.relative_to(temp_state_dir.parent)
        # We need to be in the right directory for relative paths to work
        # Instead, test with a symlink
        symlink_path = temp_state_dir / "symlink_to_db"
        symlink_path.symlink_to(temp_db_path)

        db, warning = get_database_by_path("blast", symlink_path, temp_state_dir)

        assert db is not None
        assert db.version == "v1.0"

    def test_get_database_by_path_different_types_independent(
        self, temp_state_dir, temp_db_path
    ):
        """Different db_types at same path are independent."""
        from py_nvd.state import get_database_by_path, register_database

        register_database(
            "blast", "blast-v1", str(temp_db_path), state_dir=temp_state_dir
        )
        register_database(
            "gottcha2", "gottcha-v1", str(temp_db_path), state_dir=temp_state_dir
        )

        blast_db, _ = get_database_by_path("blast", temp_db_path, temp_state_dir)
        gottcha_db, _ = get_database_by_path("gottcha2", temp_db_path, temp_state_dir)

        assert blast_db is not None
        assert blast_db.version == "blast-v1"
        assert gottcha_db is not None
        assert gottcha_db.version == "gottcha-v1"

    def test_get_databases_by_path_returns_all_versions(
        self, temp_state_dir, temp_db_path
    ):
        """get_databases_by_path returns all versions at a path."""
        import time

        from py_nvd.state import get_databases_by_path, register_database

        register_database("blast", "v1.0", str(temp_db_path), state_dir=temp_state_dir)
        time.sleep(0.1)
        register_database("blast", "v2.0", str(temp_db_path), state_dir=temp_state_dir)
        time.sleep(0.1)
        register_database("blast", "v3.0", str(temp_db_path), state_dir=temp_state_dir)

        databases = get_databases_by_path("blast", temp_db_path, temp_state_dir)

        assert len(databases) == 3
        # Should be ordered newest first
        assert databases[0].version == "v3.0"
        assert databases[1].version == "v2.0"
        assert databases[2].version == "v1.0"

    def test_get_databases_by_path_empty_when_not_found(
        self, temp_state_dir, temp_db_path
    ):
        """get_databases_by_path returns empty list when path not registered."""
        from py_nvd.state import get_databases_by_path

        databases = get_databases_by_path("blast", temp_db_path, temp_state_dir)

        assert databases == []

    def test_get_database_by_version_still_works(self, temp_state_dir, temp_db_path):
        """Renamed get_database_by_version still works correctly."""
        from py_nvd.state import get_database_by_version, register_database

        register_database("blast", "v1.0", str(temp_db_path), state_dir=temp_state_dir)

        # Lookup by specific version
        db = get_database_by_version("blast", "v1.0", temp_state_dir)
        assert db is not None
        assert db.version == "v1.0"

        # Lookup latest (version=None)
        db = get_database_by_version("blast", None, temp_state_dir)
        assert db is not None
        assert db.version == "v1.0"

        # Lookup non-existent version
        db = get_database_by_version("blast", "v999", temp_state_dir)
        assert db is None


class TestResolveDatabaseVersions:
    """Tests for resolve_database_versions() and _resolve_single_database()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def temp_db_paths(self, temp_state_dir):
        """Create temporary 'database' directories for each type."""
        paths = {}
        for db_type in ["blast", "gottcha2", "stat"]:
            db_dir = temp_state_dir / f"fake_{db_type}_db"
            db_dir.mkdir()
            paths[db_type] = db_dir
        return paths

    # =========================================================================
    # Case 1: No paths provided - nothing to resolve
    # =========================================================================

    def test_no_paths_returns_empty_resolution(self, temp_state_dir):
        """When no paths provided, returns empty resolution."""
        from py_nvd.state import resolve_database_versions

        resolution = resolve_database_versions(state_dir=temp_state_dir)

        assert resolution.blast_db_version is None
        assert resolution.gottcha2_db_version is None
        assert resolution.stat_db_version is None
        assert resolution.warnings == []
        assert resolution.auto_registered == []

    def test_version_without_path_passes_through(self, temp_state_dir):
        """Version provided without path is returned as-is (unusual but valid)."""
        from py_nvd.state import resolve_database_versions

        resolution = resolve_database_versions(
            blast_db_version="v1.0",
            state_dir=temp_state_dir,
        )

        # Version passes through even without path
        assert resolution.blast_db_version == "v1.0"
        assert resolution.warnings == []
        assert resolution.auto_registered == []

    # =========================================================================
    # Case 2: Path and version both provided - auto-register
    # =========================================================================

    def test_path_and_version_auto_registers_new(self, temp_state_dir, temp_db_paths):
        """Path + version auto-registers when not in registry."""
        from py_nvd.state import get_database_by_version, resolve_database_versions

        resolution = resolve_database_versions(
            blast_db=temp_db_paths["blast"],
            blast_db_version="core-nt_2025-01",
            state_dir=temp_state_dir,
        )

        # Version is resolved
        assert resolution.blast_db_version == "core-nt_2025-01"
        assert resolution.warnings == []

        # Was auto-registered
        assert len(resolution.auto_registered) == 1
        assert resolution.auto_registered[0][0] == "blast"
        assert resolution.auto_registered[0][1] == "core-nt_2025-01"

        # Verify it's in the database
        db = get_database_by_version("blast", "core-nt_2025-01", temp_state_dir)
        assert db is not None
        assert db.version == "core-nt_2025-01"

    def test_path_and_version_same_as_registered_is_noop(
        self, temp_state_dir, temp_db_paths
    ):
        """Path + version matching registry is a no-op."""
        from py_nvd.state import register_database, resolve_database_versions

        # Pre-register
        register_database(
            "blast",
            "core-nt_2025-01",
            str(temp_db_paths["blast"]),
            state_dir=temp_state_dir,
        )

        resolution = resolve_database_versions(
            blast_db=temp_db_paths["blast"],
            blast_db_version="core-nt_2025-01",
            state_dir=temp_state_dir,
        )

        # Version is resolved
        assert resolution.blast_db_version == "core-nt_2025-01"
        # No warnings, no new registrations
        assert resolution.warnings == []
        assert resolution.auto_registered == []

    def test_path_and_version_different_from_registered_warns_and_updates(
        self, temp_state_dir, temp_db_paths
    ):
        """Path + version different from registry warns and updates."""
        from py_nvd.state import (
            get_database_by_version,
            register_database,
            resolve_database_versions,
        )

        # Pre-register with old version
        register_database(
            "blast",
            "core-nt_2024-01",
            str(temp_db_paths["blast"]),
            state_dir=temp_state_dir,
        )

        resolution = resolve_database_versions(
            blast_db=temp_db_paths["blast"],
            blast_db_version="core-nt_2025-01",  # Different version
            state_dir=temp_state_dir,
        )

        # Uses provided version
        assert resolution.blast_db_version == "core-nt_2025-01"

        # Warns about new registration (old one preserved)
        assert len(resolution.warnings) == 1
        assert "registering new version" in resolution.warnings[0].lower()
        assert "core-nt_2024-01" in resolution.warnings[0]
        assert "core-nt_2025-01" in resolution.warnings[0]

        # Was registered
        assert len(resolution.auto_registered) == 1

        # Verify new version is in database
        db = get_database_by_version("blast", "core-nt_2025-01", temp_state_dir)
        assert db is not None

        # Verify old version is ALSO still in database (history preserved)
        old_db = get_database_by_version("blast", "core-nt_2024-01", temp_state_dir)
        assert old_db is not None

    # =========================================================================
    # Case 3: Path provided without version - auto-resolve from registry
    # =========================================================================

    def test_path_without_version_resolves_from_registry(
        self, temp_state_dir, temp_db_paths
    ):
        """Path without version resolves from registry."""
        from py_nvd.state import register_database, resolve_database_versions

        # Pre-register
        register_database(
            "blast",
            "core-nt_2025-01",
            str(temp_db_paths["blast"]),
            state_dir=temp_state_dir,
        )

        resolution = resolve_database_versions(
            blast_db=temp_db_paths["blast"],
            # No blast_db_version provided
            state_dir=temp_state_dir,
        )

        # Resolved from registry
        assert resolution.blast_db_version == "core-nt_2025-01"
        assert resolution.warnings == []
        assert resolution.auto_registered == []

    def test_path_without_version_not_registered_warns(
        self, temp_state_dir, temp_db_paths
    ):
        """Path without version, not in registry, warns."""
        from py_nvd.state import resolve_database_versions

        resolution = resolve_database_versions(
            blast_db=temp_db_paths["blast"],
            # No blast_db_version provided, not registered
            state_dir=temp_state_dir,
        )

        # Cannot resolve
        assert resolution.blast_db_version is None

        # Warns about unregistered path
        assert len(resolution.warnings) == 1
        assert "not registered" in resolution.warnings[0]
        assert "blast" in resolution.warnings[0]
        assert "nvd state database register" in resolution.warnings[0]

        # No registrations
        assert resolution.auto_registered == []

    # =========================================================================
    # Multiple database types
    # =========================================================================

    def test_resolves_all_database_types(self, temp_state_dir, temp_db_paths):
        """Resolves all three database types independently."""
        from py_nvd.state import resolve_database_versions

        resolution = resolve_database_versions(
            blast_db=temp_db_paths["blast"],
            blast_db_version="blast-v1",
            gottcha2_db=temp_db_paths["gottcha2"],
            gottcha2_db_version="gottcha-v1",
            stat_index=temp_db_paths["stat"],
            stat_db_version="stat-v1",
            state_dir=temp_state_dir,
        )

        assert resolution.blast_db_version == "blast-v1"
        assert resolution.gottcha2_db_version == "gottcha-v1"
        assert resolution.stat_db_version == "stat-v1"
        assert resolution.warnings == []
        assert len(resolution.auto_registered) == 3

    def test_mixed_resolution_scenarios(self, temp_state_dir, temp_db_paths):
        """Different resolution scenarios for different db types."""
        from py_nvd.state import register_database, resolve_database_versions

        # Pre-register BLAST
        register_database(
            "blast",
            "blast-v1",
            str(temp_db_paths["blast"]),
            state_dir=temp_state_dir,
        )

        resolution = resolve_database_versions(
            # BLAST: path only, should resolve from registry
            blast_db=temp_db_paths["blast"],
            # GOTTCHA2: path + version, should auto-register
            gottcha2_db=temp_db_paths["gottcha2"],
            gottcha2_db_version="gottcha-v1",
            # STAT: path only, not registered, should warn
            stat_index=temp_db_paths["stat"],
            state_dir=temp_state_dir,
        )

        # BLAST resolved from registry
        assert resolution.blast_db_version == "blast-v1"

        # GOTTCHA2 auto-registered
        assert resolution.gottcha2_db_version == "gottcha-v1"
        assert any(r[0] == "gottcha2" for r in resolution.auto_registered)

        # STAT unresolved with warning
        assert resolution.stat_db_version is None
        assert any("stat" in w and "not registered" in w for w in resolution.warnings)

    # =========================================================================
    # Edge cases
    # =========================================================================

    def test_path_canonicalization_in_resolution(self, temp_state_dir, temp_db_paths):
        """Paths are canonicalized during resolution."""
        from py_nvd.state import register_database, resolve_database_versions

        # Register with absolute path
        register_database(
            "blast",
            "v1.0",
            str(temp_db_paths["blast"].resolve()),
            state_dir=temp_state_dir,
        )

        # Create symlink
        symlink = temp_state_dir / "blast_symlink"
        symlink.symlink_to(temp_db_paths["blast"])

        # Resolve using symlink
        resolution = resolve_database_versions(
            blast_db=symlink,
            state_dir=temp_state_dir,
        )

        # Should resolve via symlink
        assert resolution.blast_db_version == "v1.0"
        assert resolution.warnings == []

    def test_multiple_versions_at_path_uses_newest_with_warning(
        self, temp_state_dir, temp_db_paths
    ):
        """When multiple versions at path, uses newest and warns."""
        from py_nvd.state import register_database, resolve_database_versions

        # Register multiple versions at same path
        register_database(
            "blast", "v1.0", str(temp_db_paths["blast"]), state_dir=temp_state_dir
        )
        register_database(
            "blast", "v2.0", str(temp_db_paths["blast"]), state_dir=temp_state_dir
        )

        resolution = resolve_database_versions(
            blast_db=temp_db_paths["blast"],
            state_dir=temp_state_dir,
        )

        # Uses newest
        assert resolution.blast_db_version == "v2.0"

        # Warns about multiple versions
        assert len(resolution.warnings) == 1
        assert "Multiple versions" in resolution.warnings[0]
