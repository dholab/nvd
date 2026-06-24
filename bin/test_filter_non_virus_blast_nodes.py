"""Tests for filtering annotated BLAST rows to viral query groups."""

from __future__ import annotations

import pandas as pd
from filter_non_virus_blast_nodes import contains_non_phage_viruses


def group_with(*, rank: str, stitle: str = "viral reference") -> pd.DataFrame:
    """Build a minimal qseqid group for filter predicate tests."""
    return pd.DataFrame(
        {
            "qseqid": ["contig1"],
            "rank": [rank],
            "stitle": [stitle],
        },
    )


def test_superkingdom_viruses_lineage_is_kept() -> None:
    """Legacy superkingdom lineage strings are still viral."""
    group = group_with(rank="root:cellular organisms; superkingdom:Viruses")

    assert contains_non_phage_viruses(group)


def test_acellular_root_viruses_lineage_is_kept() -> None:
    """Current taxonomy strings can identify Viruses at acellular root."""
    group = group_with(rank="acellular root:Viruses; realm:Riboviria")

    assert contains_non_phage_viruses(group)


def test_non_viral_lineage_is_rejected() -> None:
    """Non-viral taxonomy strings should not keep the query group."""
    group = group_with(rank="domain:Bacteria; phylum:Pseudomonadota")

    assert not contains_non_phage_viruses(group)


def test_viral_phage_only_group_is_rejected() -> None:
    """Phage-only viral query groups keep the existing exclusion behavior."""
    group = group_with(
        rank="acellular root:Viruses; realm:Duplodnaviria",
        stitle="Escherichia phage lambda",
    )

    assert not contains_non_phage_viruses(group)


def test_substring_lookalike_is_not_viral() -> None:
    """The viral predicate should match exact token names, not substrings."""
    group = group_with(rank="note:NotVirusesMaybe; domain:Bacteria")

    assert not contains_non_phage_viruses(group)
