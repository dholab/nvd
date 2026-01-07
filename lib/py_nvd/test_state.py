"""Tests for py_nvd.state module."""

import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

from py_nvd.db import EXPECTED_VERSION, SchemaMismatchError, connect
from py_nvd.models import Run
from py_nvd.state import (
    check_taxonomy_drift,
    check_upload,
    complete_run,
    complete_sample,
    compute_sample_set_id,
    get_database_by_path,
    get_database_by_version,
    get_databases_by_path,
    get_processed_sample,
    get_run,
    get_run_by_sample_set,
    get_sample_history,
    get_samples_for_run,
    get_samples_needing_upload,
    get_taxonomy_version,
    get_upload,
    get_uploaded_sample_ids,
    get_uploads_for_sample,
    hash_upload_content,
    list_runs,
    list_samples,
    list_uploads,
    mark_sample_uploaded,
    record_taxonomy_version,
    record_upload,
    register_database,
    register_processed_sample,
    register_run,
    resolve_database_versions,
    was_sample_ever_processed,
    was_sample_ever_uploaded,
)


@dataclass
class RunContext:
    """Context for a registered run fixture."""

    run: Run
    sample_set_id: str
    state_dir: Path


@dataclass
class ProcessedSampleContext:
    """Context for a processed sample fixture."""

    sample_id: str
    sample_set_id: str
    run_id: str
    state_dir: Path


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
    def registered_run(self, temp_state_dir) -> RunContext:
        """Register a run and return RunContext."""
        sample_set_id = compute_sample_set_id(["s1", "s2", "s3"])
        run = register_run("run_001", sample_set_id, state_dir=temp_state_dir)
        assert run is not None
        return RunContext(
            run=run, sample_set_id=sample_set_id, state_dir=temp_state_dir
        )

    def test_register_processed_sample(self, registered_run: RunContext):
        """Registering a processed sample succeeds with 'completed' status.

        Note: register_processed_sample is an alias for mark_sample_completed.
        Samples are now inserted as 'completed' (not 'processing') because
        this is called after REGISTER_HITS succeeds.
        """
        ctx = registered_run

        sample = register_processed_sample(
            sample_id="s1",
            sample_set_id=ctx.sample_set_id,
            run_id=ctx.run.run_id,
            blast_db_version="core-nt_2025-01-01",
            stat_db_version="v1.0",
            state_dir=ctx.state_dir,
        )

        assert sample is not None
        assert sample.sample_id == "s1"
        assert sample.sample_set_id == ctx.sample_set_id
        assert sample.run_id == "run_001"
        assert sample.status == "completed"  # Now starts as completed
        assert sample.blast_db_version == "core-nt_2025-01-01"
        assert sample.stat_db_version == "v1.0"

    def test_get_processed_sample(self, registered_run: RunContext):
        """get_processed_sample retrieves a registered sample."""
        ctx = registered_run

        register_processed_sample(
            sample_id="s1",
            sample_set_id=ctx.sample_set_id,
            run_id=ctx.run.run_id,
            state_dir=ctx.state_dir,
        )

        sample = get_processed_sample("s1", ctx.sample_set_id, state_dir=ctx.state_dir)
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

    def test_get_samples_for_run(self, registered_run: RunContext):
        """get_samples_for_run returns all samples for a run."""
        ctx = registered_run

        for sid in ["s1", "s2", "s3"]:
            register_processed_sample(
                sample_id=sid,
                sample_set_id=ctx.sample_set_id,
                run_id=ctx.run.run_id,
                state_dir=ctx.state_dir,
            )

        samples = get_samples_for_run(ctx.run.run_id, state_dir=ctx.state_dir)
        assert len(samples) == 3
        sample_ids = {s.sample_id for s in samples}
        assert sample_ids == {"s1", "s2", "s3"}

    def test_was_sample_ever_processed(self, registered_run: RunContext):
        """was_sample_ever_processed returns True for processed samples."""
        ctx = registered_run

        # Not processed yet
        assert was_sample_ever_processed("s1", state_dir=ctx.state_dir) is False

        # Process it
        register_processed_sample(
            sample_id="s1",
            sample_set_id=ctx.sample_set_id,
            run_id=ctx.run.run_id,
            state_dir=ctx.state_dir,
        )

        # Now it's processed
        assert was_sample_ever_processed("s1", state_dir=ctx.state_dir) is True
        # Other samples still not processed
        assert was_sample_ever_processed("s2", state_dir=ctx.state_dir) is False

    def test_complete_sample(self, registered_run: RunContext):
        """complete_sample can transition sample to 'uploaded' status.

        Note: Samples now start as 'completed' (after REGISTER_HITS).
        The complete_sample function is used to transition to 'uploaded'
        after LabKey upload succeeds.
        """
        ctx = registered_run

        register_processed_sample(
            sample_id="s1",
            sample_set_id=ctx.sample_set_id,
            run_id=ctx.run.run_id,
            state_dir=ctx.state_dir,
        )

        # Initially completed (after REGISTER_HITS)
        sample = get_processed_sample("s1", ctx.sample_set_id, state_dir=ctx.state_dir)
        assert sample is not None
        assert sample.status == "completed"

        # Mark as uploaded (after LabKey upload)
        complete_sample(
            "s1", ctx.sample_set_id, status="uploaded", state_dir=ctx.state_dir
        )

        sample = get_processed_sample("s1", ctx.sample_set_id, state_dir=ctx.state_dir)
        assert sample is not None
        assert sample.status == "uploaded"

    def test_complete_sample_failed(self, registered_run: RunContext):
        """complete_sample can mark a sample as failed."""
        ctx = registered_run

        register_processed_sample(
            sample_id="s1",
            sample_set_id=ctx.sample_set_id,
            run_id=ctx.run.run_id,
            state_dir=ctx.state_dir,
        )

        complete_sample(
            "s1", ctx.sample_set_id, status="failed", state_dir=ctx.state_dir
        )

        sample = get_processed_sample("s1", ctx.sample_set_id, state_dir=ctx.state_dir)
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
        """list_samples filters by status.

        Note: Samples now start as 'completed' (not 'processing').
        Valid statuses are: completed, uploaded, failed.
        """
        set1 = compute_sample_set_id(["s1", "s2", "s3"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run.run_id, state_dir=temp_state_dir)
        register_processed_sample("s3", set1, run.run_id, state_dir=temp_state_dir)

        # All start as 'completed', transition some to other states
        complete_sample("s1", set1, status="uploaded", state_dir=temp_state_dir)
        complete_sample("s2", set1, status="failed", state_dir=temp_state_dir)
        # s3 remains "completed"

        uploaded = list_samples(status="uploaded", state_dir=temp_state_dir)
        assert len(uploaded) == 1
        assert uploaded[0].sample_id == "s1"

        failed = list_samples(status="failed", state_dir=temp_state_dir)
        assert len(failed) == 1
        assert failed[0].sample_id == "s2"

        completed = list_samples(status="completed", state_dir=temp_state_dir)
        assert len(completed) == 1
        assert completed[0].sample_id == "s3"

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
        """list_samples combines multiple filters.

        Note: Samples now start as 'completed'. We transition some to
        'uploaded' to test filtering.
        """
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

        # Transition some to 'uploaded'
        complete_sample("s1", set1, status="uploaded", state_dir=temp_state_dir)
        complete_sample("s3", set2, status="uploaded", state_dir=temp_state_dir)

        # Uploaded samples from run_001
        samples = list_samples(
            run_id="run_001", status="uploaded", state_dir=temp_state_dir
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


class TestPerSampleStateTracking:
    """Tests for per-sample state tracking functions (mark_sample_uploaded, get_samples_needing_upload, get_uploaded_sample_ids).

    These functions support the per-sample resume feature, enabling:
    - Marking samples as 'uploaded' after LabKey upload succeeds
    - Querying which samples need upload (status='completed')
    - Checking which samples have already been uploaded (for run-start gating)
    """

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_mark_sample_uploaded_success(self, temp_state_dir):
        """mark_sample_uploaded transitions sample from completed to uploaded."""
        from py_nvd.state import mark_sample_uploaded

        set1 = compute_sample_set_id(["s1"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        # Register sample (starts as 'completed')
        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)
        sample = get_processed_sample("s1", set1, state_dir=temp_state_dir)
        assert sample is not None
        assert sample.status == "completed"

        # Mark as uploaded
        result = mark_sample_uploaded("s1", set1, state_dir=temp_state_dir)
        assert result is True

        # Verify status changed
        sample = get_processed_sample("s1", set1, state_dir=temp_state_dir)
        assert sample is not None
        assert sample.status == "uploaded"

    def test_mark_sample_uploaded_not_found(self, temp_state_dir):
        """mark_sample_uploaded returns False for non-existent sample."""
        from py_nvd.state import mark_sample_uploaded

        result = mark_sample_uploaded(
            "nonexistent", "fake_set", state_dir=temp_state_dir
        )
        assert result is False

    def test_get_samples_needing_upload(self, temp_state_dir):
        """get_samples_needing_upload returns only samples with status='completed'."""
        from py_nvd.state import get_samples_needing_upload, mark_sample_uploaded

        set1 = compute_sample_set_id(["s1", "s2", "s3"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        # Register three samples (all start as 'completed')
        register_processed_sample("s1", set1, run.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run.run_id, state_dir=temp_state_dir)
        register_processed_sample("s3", set1, run.run_id, state_dir=temp_state_dir)

        # All three should need upload initially
        needing_upload = get_samples_needing_upload(set1, state_dir=temp_state_dir)
        assert len(needing_upload) == 3
        assert {s.sample_id for s in needing_upload} == {"s1", "s2", "s3"}

        # Mark s1 as uploaded
        mark_sample_uploaded("s1", set1, state_dir=temp_state_dir)

        # Now only s2 and s3 should need upload
        needing_upload = get_samples_needing_upload(set1, state_dir=temp_state_dir)
        assert len(needing_upload) == 2
        assert {s.sample_id for s in needing_upload} == {"s2", "s3"}

        # Mark s2 as failed
        complete_sample("s2", set1, status="failed", state_dir=temp_state_dir)

        # Now only s3 should need upload (failed samples don't need upload)
        needing_upload = get_samples_needing_upload(set1, state_dir=temp_state_dir)
        assert len(needing_upload) == 1
        assert needing_upload[0].sample_id == "s3"

    def test_get_samples_needing_upload_empty_set(self, temp_state_dir):
        """get_samples_needing_upload returns empty list for unknown sample set."""
        from py_nvd.state import get_samples_needing_upload

        needing_upload = get_samples_needing_upload(
            "nonexistent_set", state_dir=temp_state_dir
        )
        assert needing_upload == []

    def test_get_uploaded_sample_ids(self, temp_state_dir):
        """get_uploaded_sample_ids returns IDs of samples with status='uploaded'."""
        from py_nvd.state import get_uploaded_sample_ids, mark_sample_uploaded

        set1 = compute_sample_set_id(["s1", "s2"])
        set2 = compute_sample_set_id(["s3", "s4"])

        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        # Register samples in two different sample sets
        register_processed_sample("s1", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s2", set1, run1.run_id, state_dir=temp_state_dir)
        register_processed_sample("s3", set2, run2.run_id, state_dir=temp_state_dir)
        register_processed_sample("s4", set2, run2.run_id, state_dir=temp_state_dir)

        # Initially none are uploaded
        uploaded = get_uploaded_sample_ids(
            ["s1", "s2", "s3", "s4"], state_dir=temp_state_dir
        )
        assert uploaded == set()

        # Mark s1 and s3 as uploaded
        mark_sample_uploaded("s1", set1, state_dir=temp_state_dir)
        mark_sample_uploaded("s3", set2, state_dir=temp_state_dir)

        # Check uploaded samples
        uploaded = get_uploaded_sample_ids(
            ["s1", "s2", "s3", "s4"], state_dir=temp_state_dir
        )
        assert uploaded == {"s1", "s3"}

        # Check subset of samples
        uploaded = get_uploaded_sample_ids(["s1", "s2"], state_dir=temp_state_dir)
        assert uploaded == {"s1"}

        # Check samples not in query
        uploaded = get_uploaded_sample_ids(["s5", "s6"], state_dir=temp_state_dir)
        assert uploaded == set()

    def test_get_uploaded_sample_ids_empty_list(self, temp_state_dir):
        """get_uploaded_sample_ids returns empty set for empty input."""
        from py_nvd.state import get_uploaded_sample_ids

        uploaded = get_uploaded_sample_ids([], state_dir=temp_state_dir)
        assert uploaded == set()

    def test_get_uploaded_sample_ids_cross_sample_set(self, temp_state_dir):
        """get_uploaded_sample_ids finds uploads across different sample sets.

        This is important for the run-start check: if a sample was uploaded
        in a previous run (different sample_set_id), we should still detect it.
        """
        from py_nvd.state import get_uploaded_sample_ids, mark_sample_uploaded

        # Same sample in two different sample sets (reprocessed)
        set1 = compute_sample_set_id(["s1", "s2"])
        set2 = compute_sample_set_id(["s1", "s3"])  # s1 appears in both

        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        # Register s1 in first run and mark uploaded
        register_processed_sample("s1", set1, run1.run_id, state_dir=temp_state_dir)
        mark_sample_uploaded("s1", set1, state_dir=temp_state_dir)

        # Register s1 in second run (completed, not uploaded yet)
        register_processed_sample("s1", set2, run2.run_id, state_dir=temp_state_dir)

        # Query should find s1 as uploaded (from first run)
        uploaded = get_uploaded_sample_ids(["s1"], state_dir=temp_state_dir)
        assert uploaded == {"s1"}


class TestUploads:
    """Tests for upload tracking functions."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def processed_sample(self, temp_state_dir) -> ProcessedSampleContext:
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
        return ProcessedSampleContext(
            sample_id="s1",
            sample_set_id=sample_set_id,
            run_id=run.run_id,
            state_dir=temp_state_dir,
        )

    def test_record_upload(self, processed_sample: ProcessedSampleContext):
        """record_upload creates an upload record."""
        ctx = processed_sample

        upload = record_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            target_metadata={"experiment_id": 42, "labkey_row_id": 100},
            state_dir=ctx.state_dir,
        )

        assert upload is not None
        assert upload.sample_id == "s1"
        assert upload.upload_type == "blast"
        assert upload.upload_target == "labkey"
        assert upload.content_hash == "abc123"
        assert upload.target_metadata == {"experiment_id": 42, "labkey_row_id": 100}

    def test_get_upload(self, processed_sample: ProcessedSampleContext):
        """get_upload retrieves a specific upload."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            state_dir=ctx.state_dir,
        )

        upload = get_upload(
            ctx.sample_id,
            ctx.sample_set_id,
            "blast",
            "labkey",
            state_dir=ctx.state_dir,
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

    def test_check_upload_not_uploaded(self, processed_sample: ProcessedSampleContext):
        """check_upload returns 'not_uploaded' for new uploads."""
        ctx = processed_sample

        result = check_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            state_dir=ctx.state_dir,
        )
        assert result == "not_uploaded"

    def test_check_upload_already_uploaded(
        self, processed_sample: ProcessedSampleContext
    ):
        """check_upload returns 'already_uploaded' for same content."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            state_dir=ctx.state_dir,
        )

        result = check_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",  # Same hash
            state_dir=ctx.state_dir,
        )
        assert result == "already_uploaded"

    def test_check_upload_content_changed(
        self, processed_sample: ProcessedSampleContext
    ):
        """check_upload returns 'content_changed' for different content."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="abc123",
            state_dir=ctx.state_dir,
        )

        result = check_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="different_hash",  # Different hash
            state_dir=ctx.state_dir,
        )
        assert result == "content_changed"

    def test_get_uploads_for_sample(self, processed_sample: ProcessedSampleContext):
        """get_uploads_for_sample returns all uploads for a sample."""
        ctx = processed_sample

        # Upload to multiple targets
        record_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="hash1",
            state_dir=ctx.state_dir,
        )
        record_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast_fasta",
            upload_target="labkey",
            content_hash="hash2",
            state_dir=ctx.state_dir,
        )

        uploads = get_uploads_for_sample(ctx.sample_id, state_dir=ctx.state_dir)
        assert len(uploads) == 2
        upload_types = {u.upload_type for u in uploads}
        assert upload_types == {"blast", "blast_fasta"}

    def test_was_sample_ever_uploaded_false(
        self, processed_sample: ProcessedSampleContext
    ):
        """was_sample_ever_uploaded returns False for un-uploaded samples."""
        ctx = processed_sample
        assert was_sample_ever_uploaded(ctx.sample_id, state_dir=ctx.state_dir) is False

    def test_was_sample_ever_uploaded_true(
        self, processed_sample: ProcessedSampleContext
    ):
        """was_sample_ever_uploaded returns True for uploaded samples."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="hash1",
            state_dir=ctx.state_dir,
        )

        assert was_sample_ever_uploaded(ctx.sample_id, state_dir=ctx.state_dir) is True

    def test_was_sample_ever_uploaded_with_filters(
        self, processed_sample: ProcessedSampleContext
    ):
        """was_sample_ever_uploaded respects upload_type and upload_target filters."""
        ctx = processed_sample

        record_upload(
            sample_id=ctx.sample_id,
            sample_set_id=ctx.sample_set_id,
            upload_type="blast",
            upload_target="labkey",
            content_hash="hash1",
            state_dir=ctx.state_dir,
        )

        # Uploaded to labkey
        assert (
            was_sample_ever_uploaded(
                ctx.sample_id,
                upload_target="labkey",
                state_dir=ctx.state_dir,
            )
            is True
        )

        # Not uploaded to globus
        assert (
            was_sample_ever_uploaded(
                ctx.sample_id,
                upload_target="globus",
                state_dir=ctx.state_dir,
            )
            is False
        )

        # Uploaded blast type
        assert (
            was_sample_ever_uploaded(
                ctx.sample_id,
                upload_type="blast",
                state_dir=ctx.state_dir,
            )
            is True
        )

        # Not uploaded gottcha2 type
        assert (
            was_sample_ever_uploaded(
                ctx.sample_id,
                upload_type="gottcha2",
                state_dir=ctx.state_dir,
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
        db_path = temp_state_dir / "state.sqlite"
        assert not db_path.exists()

        with connect(temp_state_dir) as conn:
            # Should succeed without any prompts
            version = conn.execute("PRAGMA user_version").fetchone()[0]
            assert version == 2  # EXPECTED_VERSION

        assert db_path.exists()

    def test_matching_schema_version_proceeds_normally(self, temp_state_dir):
        """Database with matching schema version works without prompts."""
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
        db, warning = get_database_by_path("blast", temp_db_path, temp_state_dir)

        assert db is None
        assert warning is None

    def test_get_database_by_path_single_match(self, temp_state_dir, temp_db_path):
        """Looking up registered path returns (database, None)."""
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
        databases = get_databases_by_path("blast", temp_db_path, temp_state_dir)

        assert databases == []

    def test_get_database_by_version_still_works(self, temp_state_dir, temp_db_path):
        """Renamed get_database_by_version still works correctly."""
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

    def test_no_paths_returns_empty_resolution(self, temp_state_dir):
        """When no paths provided, returns empty resolution."""
        resolution = resolve_database_versions(state_dir=temp_state_dir)

        assert resolution.blast_db_version is None
        assert resolution.gottcha2_db_version is None
        assert resolution.stat_db_version is None
        assert resolution.warnings == []
        assert resolution.auto_registered == []

    def test_version_without_path_passes_through(self, temp_state_dir):
        """Version provided without path is returned as-is (unusual but valid)."""
        resolution = resolve_database_versions(
            blast_db_version="v1.0",
            state_dir=temp_state_dir,
        )

        # Version passes through even without path
        assert resolution.blast_db_version == "v1.0"
        assert resolution.warnings == []
        assert resolution.auto_registered == []

    def test_path_and_version_auto_registers_new(self, temp_state_dir, temp_db_paths):
        """Path + version auto-registers when not in registry."""
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

    def test_path_without_version_resolves_from_registry(
        self, temp_state_dir, temp_db_paths
    ):
        """Path without version resolves from registry."""
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

    def test_resolves_all_database_types(self, temp_state_dir, temp_db_paths):
        """Resolves all three database types independently."""
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

    def test_path_canonicalization_in_resolution(self, temp_state_dir, temp_db_paths):
        """Paths are canonicalized during resolution."""
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


class TestSampleLocks:
    """Tests for sample lock functions.

    Sample locks prevent duplicate processing across concurrent runs.
    They use TTL-based expiration and machine fingerprinting.
    """

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def registered_run(self, temp_state_dir) -> RunContext:
        """Register a run and return RunContext."""
        sample_set_id = compute_sample_set_id(["s1", "s2", "s3"])
        run = register_run("run_001", sample_set_id, state_dir=temp_state_dir)
        assert run is not None
        return RunContext(
            run=run, sample_set_id=sample_set_id, state_dir=temp_state_dir
        )

    def test_acquire_sample_lock_basic(self, registered_run: RunContext):
        """Basic lock acquisition succeeds and returns None."""
        from py_nvd.state import acquire_sample_lock, get_lock

        ctx = registered_run

        # Acquire lock
        result = acquire_sample_lock("s1", ctx.run.run_id, state_dir=ctx.state_dir)

        # None means success (no blocking lock)
        assert result is None

        # Lock should exist
        lock = get_lock("s1", state_dir=ctx.state_dir)
        assert lock is not None
        assert lock.sample_id == "s1"
        assert lock.run_id == ctx.run.run_id
        assert not lock.is_expired

    def test_acquire_sample_lock_blocked_by_another_run(self, temp_state_dir):
        """Lock acquisition fails when blocked by another run's active lock."""
        from py_nvd.state import acquire_sample_lock

        # Register two runs
        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])
        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        # Run 1 acquires lock
        result1 = acquire_sample_lock("s1", "run_001", state_dir=temp_state_dir)
        assert result1 is None  # Success

        # Run 2 tries to acquire same sample - should be blocked
        result2 = acquire_sample_lock("s1", "run_002", state_dir=temp_state_dir)
        assert result2 is not None  # Returns blocking lock
        assert result2.run_id == "run_001"
        assert result2.sample_id == "s1"

    def test_acquire_sample_lock_overwrites_expired(self, temp_state_dir):
        """Lock acquisition succeeds when existing lock is expired."""
        from datetime import datetime, timedelta, timezone

        from py_nvd.state import acquire_sample_lock, get_lock

        # Register a run
        set1 = compute_sample_set_id(["s1"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        # Insert an expired lock directly
        now = datetime.now(timezone.utc)
        expired = now - timedelta(hours=1)
        with connect(temp_state_dir) as conn:
            conn.execute(
                """INSERT INTO sample_locks 
                   (sample_id, run_id, hostname, username, locked_at, expires_at)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    "s1",
                    "old_run",
                    "old_host",
                    "old_user",
                    now.isoformat(),
                    expired.isoformat(),
                ),
            )
            conn.commit()

        # New run should be able to acquire (steal expired lock)
        result = acquire_sample_lock("s1", "run_001", state_dir=temp_state_dir)
        assert result is None  # Success

        # Lock should now be owned by run_001
        lock = get_lock("s1", state_dir=temp_state_dir)
        assert lock is not None
        assert lock.run_id == "run_001"

    def test_acquire_sample_lock_same_run_same_fingerprint_refreshes(
        self, registered_run: RunContext
    ):
        """Same run + same fingerprint refreshes TTL (legitimate resume)."""
        from py_nvd.state import acquire_sample_lock, get_lock

        ctx = registered_run

        # First acquisition
        result1 = acquire_sample_lock("s1", ctx.run.run_id, state_dir=ctx.state_dir)
        assert result1 is None

        lock1 = get_lock("s1", state_dir=ctx.state_dir)
        assert lock1 is not None
        expires1 = lock1.expires_at

        # Small delay to ensure different timestamp
        import time

        time.sleep(0.1)

        # Second acquisition (same run, same machine) - should refresh TTL
        result2 = acquire_sample_lock("s1", ctx.run.run_id, state_dir=ctx.state_dir)
        assert result2 is None  # Success (refresh)

        lock2 = get_lock("s1", state_dir=ctx.state_dir)
        assert lock2 is not None
        # TTL should be refreshed (expires_at should be later)
        assert lock2.expires_at >= expires1

    def test_acquire_sample_lock_same_run_different_fingerprint_active_conflict(
        self, temp_state_dir
    ):
        """Same run + different fingerprint + active lock = conflict."""
        from datetime import datetime, timedelta, timezone

        from py_nvd.state import acquire_sample_lock

        # Register a run
        set1 = compute_sample_set_id(["s1"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        # Insert a lock with different fingerprint (simulating another machine)
        now = datetime.now(timezone.utc)
        expires = now + timedelta(hours=72)
        with connect(temp_state_dir) as conn:
            conn.execute(
                """INSERT INTO sample_locks 
                   (sample_id, run_id, hostname, username, locked_at, expires_at)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    "s1",
                    "run_001",
                    "other_host",
                    "other_user",
                    now.isoformat(),
                    expires.isoformat(),
                ),
            )
            conn.commit()

        # Same run_id but different machine - should conflict
        result = acquire_sample_lock("s1", "run_001", state_dir=temp_state_dir)
        assert result is not None  # Conflict
        assert result.run_id == "run_001"
        assert result.hostname == "other_host"
        assert result.username == "other_user"

    def test_acquire_sample_lock_same_run_different_fingerprint_expired_acquires(
        self, temp_state_dir
    ):
        """Same run + different fingerprint + expired = crash recovery."""
        from datetime import datetime, timedelta, timezone

        from py_nvd.state import acquire_sample_lock, get_lock

        # Register a run
        set1 = compute_sample_set_id(["s1"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        # Insert an expired lock with different fingerprint
        now = datetime.now(timezone.utc)
        expired = now - timedelta(hours=1)
        with connect(temp_state_dir) as conn:
            conn.execute(
                """INSERT INTO sample_locks 
                   (sample_id, run_id, hostname, username, locked_at, expires_at)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    "s1",
                    "run_001",
                    "crashed_host",
                    "crashed_user",
                    now.isoformat(),
                    expired.isoformat(),
                ),
            )
            conn.commit()

        # Same run_id, different machine, but expired - should acquire (crash recovery)
        result = acquire_sample_lock("s1", "run_001", state_dir=temp_state_dir)
        assert result is None  # Success

        # Lock should now have current machine's fingerprint
        lock = get_lock("s1", state_dir=temp_state_dir)
        assert lock is not None
        assert lock.run_id == "run_001"
        # Fingerprint should be updated to current machine
        assert lock.hostname != "crashed_host"

    def test_acquire_sample_locks_batch_mixed_results(self, temp_state_dir):
        """Batch acquisition returns acquired and conflicts separately."""
        from datetime import datetime, timedelta, timezone

        from py_nvd.state import acquire_sample_locks

        # Register two runs
        set1 = compute_sample_set_id(["s1", "s2", "s3"])
        set2 = compute_sample_set_id(["s4"])
        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        # Run 1 locks s1 and s2
        now = datetime.now(timezone.utc)
        expires = now + timedelta(hours=72)
        with connect(temp_state_dir) as conn:
            for sample_id in ["s1", "s2"]:
                conn.execute(
                    """INSERT INTO sample_locks 
                       (sample_id, run_id, hostname, username, locked_at, expires_at)
                       VALUES (?, ?, ?, ?, ?, ?)""",
                    (
                        sample_id,
                        "run_001",
                        "host1",
                        "user1",
                        now.isoformat(),
                        expires.isoformat(),
                    ),
                )
            conn.commit()

        # Run 2 tries to acquire s1, s2, s3 (s1, s2 blocked, s3 available)
        acquired, conflicts = acquire_sample_locks(
            ["s1", "s2", "s3"], "run_002", state_dir=temp_state_dir
        )

        assert acquired == ["s3"]
        assert len(conflicts) == 2
        conflict_samples = {c.sample_id for c in conflicts}
        assert conflict_samples == {"s1", "s2"}

    def test_release_sample_lock_releases_own(self, registered_run: RunContext):
        """release_sample_lock releases lock held by this run."""
        from py_nvd.state import acquire_sample_lock, get_lock, release_sample_lock

        ctx = registered_run

        # Acquire lock
        acquire_sample_lock("s1", ctx.run.run_id, state_dir=ctx.state_dir)
        assert get_lock("s1", state_dir=ctx.state_dir) is not None

        # Release lock
        released = release_sample_lock("s1", ctx.run.run_id, state_dir=ctx.state_dir)
        assert released is True

        # Lock should be gone
        assert get_lock("s1", state_dir=ctx.state_dir) is None

    def test_release_sample_lock_wont_release_others(self, temp_state_dir):
        """release_sample_lock won't release another run's lock."""
        from py_nvd.state import acquire_sample_lock, get_lock, release_sample_lock

        # Register two runs
        set1 = compute_sample_set_id(["s1"])
        set2 = compute_sample_set_id(["s2"])
        run1 = register_run("run_001", set1, state_dir=temp_state_dir)
        run2 = register_run("run_002", set2, state_dir=temp_state_dir)
        assert run1 is not None
        assert run2 is not None

        # Run 1 acquires lock
        acquire_sample_lock("s1", "run_001", state_dir=temp_state_dir)

        # Run 2 tries to release - should fail
        released = release_sample_lock("s1", "run_002", state_dir=temp_state_dir)
        assert released is False

        # Lock should still exist
        lock = get_lock("s1", state_dir=temp_state_dir)
        assert lock is not None
        assert lock.run_id == "run_001"

    def test_release_all_locks_for_run(self, registered_run: RunContext):
        """release_all_locks_for_run releases all locks for a run."""
        from py_nvd.state import (
            acquire_sample_lock,
            get_locks_for_run,
            release_all_locks_for_run,
        )

        ctx = registered_run

        # Acquire multiple locks
        for sample_id in ["s1", "s2", "s3"]:
            acquire_sample_lock(sample_id, ctx.run.run_id, state_dir=ctx.state_dir)

        assert len(get_locks_for_run(ctx.run.run_id, state_dir=ctx.state_dir)) == 3

        # Release all
        count = release_all_locks_for_run(ctx.run.run_id, state_dir=ctx.state_dir)
        assert count == 3

        # All locks should be gone
        assert len(get_locks_for_run(ctx.run.run_id, state_dir=ctx.state_dir)) == 0

    def test_get_lock_returns_none_for_expired(self, temp_state_dir):
        """get_lock returns None for expired locks by default."""
        from datetime import datetime, timedelta, timezone

        from py_nvd.state import get_lock

        # Register a run
        set1 = compute_sample_set_id(["s1"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        # Insert an expired lock
        now = datetime.now(timezone.utc)
        expired = now - timedelta(hours=1)
        with connect(temp_state_dir) as conn:
            conn.execute(
                """INSERT INTO sample_locks 
                   (sample_id, run_id, hostname, username, locked_at, expires_at)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                ("s1", "run_001", "host", "user", now.isoformat(), expired.isoformat()),
            )
            conn.commit()

        # Default: expired locks are not returned
        assert get_lock("s1", state_dir=temp_state_dir) is None

        # With include_expired: expired locks are returned
        lock = get_lock("s1", state_dir=temp_state_dir, include_expired=True)
        assert lock is not None
        assert lock.is_expired

    def test_cleanup_expired_locks_removes_only_expired(self, temp_state_dir):
        """cleanup_expired_locks removes only expired locks."""
        from datetime import datetime, timedelta

        from py_nvd.state import cleanup_expired_locks, list_locks

        # Register a run
        set1 = compute_sample_set_id(["s1", "s2"])
        run = register_run("run_001", set1, state_dir=temp_state_dir)
        assert run is not None

        # Use local timezone to match cleanup_expired_locks implementation
        now = datetime.now().astimezone()
        expired = now - timedelta(hours=1)
        active = now + timedelta(hours=72)

        # Insert one expired and one active lock
        with connect(temp_state_dir) as conn:
            conn.execute(
                """INSERT INTO sample_locks 
                   (sample_id, run_id, hostname, username, locked_at, expires_at)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                ("s1", "run_001", "host", "user", now.isoformat(), expired.isoformat()),
            )
            conn.execute(
                """INSERT INTO sample_locks 
                   (sample_id, run_id, hostname, username, locked_at, expires_at)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                ("s2", "run_001", "host", "user", now.isoformat(), active.isoformat()),
            )
            conn.commit()

        # Both locks exist (including expired)
        all_locks = list_locks(state_dir=temp_state_dir, include_expired=True)
        assert len(all_locks) == 2

        # Cleanup
        removed = cleanup_expired_locks(state_dir=temp_state_dir)
        assert removed == 1

        # Only active lock remains
        remaining = list_locks(state_dir=temp_state_dir, include_expired=True)
        assert len(remaining) == 1
        assert remaining[0].sample_id == "s2"


class TestEffectiveRunStatus:
    """Tests for run status with TTL-based failure detection."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_get_effective_run_status_returns_failed_for_stale(self, temp_state_dir):
        """get_effective_run_status returns 'failed' for stale 'running' runs."""
        from py_nvd.state import get_effective_run_status

        # Insert a run that started long ago
        with connect(temp_state_dir) as conn:
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, started_at, status)
                   VALUES (?, ?, ?, ?)""",
                ("old_run", "set123", "2020-01-01T00:00:00", "running"),
            )
            conn.commit()

        run = get_run("old_run", state_dir=temp_state_dir)
        assert run is not None
        assert run.status == "running"

        # Effective status should be 'failed' (way past TTL)
        effective = get_effective_run_status(run, ttl_hours=72)
        assert effective == "failed"

    def test_get_effective_run_status_returns_running_for_fresh(self, temp_state_dir):
        """get_effective_run_status returns 'running' for fresh runs."""
        from py_nvd.state import get_effective_run_status

        set1 = compute_sample_set_id(["s1"])
        run = register_run("fresh_run", set1, state_dir=temp_state_dir)
        assert run is not None
        assert run.status == "running"

        # Effective status should still be 'running' (just started)
        effective = get_effective_run_status(run, ttl_hours=72)
        assert effective == "running"

    def test_get_effective_run_status_returns_actual_for_completed(
        self, temp_state_dir
    ):
        """get_effective_run_status returns actual status for non-running runs."""
        from py_nvd.state import get_effective_run_status

        set1 = compute_sample_set_id(["s1"])
        run = register_run("completed_run", set1, state_dir=temp_state_dir)
        assert run is not None

        # Mark as completed
        complete_run("completed_run", status="completed", state_dir=temp_state_dir)

        run = get_run("completed_run", state_dir=temp_state_dir)
        assert run is not None
        assert run.status == "completed"

        # Effective status should be 'completed' (not affected by TTL)
        effective = get_effective_run_status(run, ttl_hours=72)
        assert effective == "completed"

    def test_mark_stale_runs_failed_updates_only_stale(self, temp_state_dir):
        """mark_stale_runs_failed updates only stale 'running' runs."""
        from datetime import datetime

        from py_nvd.state import mark_stale_runs_failed

        # Insert runs with different ages
        with connect(temp_state_dir) as conn:
            # Very old run (should be marked failed)
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, started_at, status)
                   VALUES (?, ?, ?, ?)""",
                ("old_run", "set1", "2020-01-01T00:00:00", "running"),
            )
            # Fresh run (should stay running)
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, started_at, status)
                   VALUES (?, ?, ?, ?)""",
                ("fresh_run", "set2", datetime.now().isoformat(), "running"),
            )
            # Already completed (should stay completed)
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, started_at, completed_at, status)
                   VALUES (?, ?, ?, ?, ?)""",
                (
                    "done_run",
                    "set3",
                    "2020-01-01T00:00:00",
                    "2020-01-02T00:00:00",
                    "completed",
                ),
            )
            conn.commit()

        # Mark stale runs as failed
        count = mark_stale_runs_failed(ttl_hours=72, state_dir=temp_state_dir)
        assert count == 1  # Only old_run

        # Verify statuses
        old = get_run("old_run", state_dir=temp_state_dir)
        assert old is not None
        assert old.status == "failed"

        fresh = get_run("fresh_run", state_dir=temp_state_dir)
        assert fresh is not None
        assert fresh.status == "running"

        done = get_run("done_run", state_dir=temp_state_dir)
        assert done is not None
        assert done.status == "completed"
