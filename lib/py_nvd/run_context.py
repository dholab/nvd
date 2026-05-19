"""Stateless run-context helpers."""

from __future__ import annotations

import hashlib


def compute_sample_set_id(sample_ids: list[str]) -> str:
    """Compute a deterministic 16-character identifier for a sample set."""
    canonical = sorted(set(sample_ids))
    content = "\n".join(canonical)
    return hashlib.sha256(content.encode()).hexdigest()[:16]
