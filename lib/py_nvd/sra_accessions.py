"""Small local checks for SRA accession-shaped strings."""

from __future__ import annotations

import re

SRA_ACCESSION = re.compile(r"^[SED]RR\d+$", re.IGNORECASE)


def looks_like_sra_accession(value: str) -> bool:
    """Return true when a string has the common SRR/ERR/DRR accession shape."""
    return SRA_ACCESSION.fullmatch(value) is not None
