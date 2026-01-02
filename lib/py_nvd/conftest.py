"""Shared pytest fixtures for py_nvd tests."""

import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def temp_state_dir():
    """Create a temporary directory for state database."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)
