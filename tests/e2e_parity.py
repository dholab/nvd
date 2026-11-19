"""Compare what a reduced e2e run detects against a recorded full-size run.

The e2e samplesheet was cut from five SRA runs to three, one per distinct code
path. That is only safe if the three remaining runs still detect the same
organisms the five did. These functions make that a mechanical check instead of
a judgement call: ``record_e2e_parity.py`` freezes a summary of a full run into
``tests/data/e2e_parity.json``, and the e2e compares against it.

The comparison is deliberately asymmetric. Detecting *fewer* things than the
full run is expected -- that is the whole point of dropping two samples, and the
e2e separately asserts that each sample still finds the organism it must. What
would signal a real behaviour change is the reduced run detecting something the
full run never did, so extra taxids are errors.

Reference accessions are only warnings. Which genome inside a taxon wins a BLAST
hit depends on the BLAST database build, so an accession swap within an
unchanged taxon is routine drift rather than a regression.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

DATA = Path(__file__).resolve().parent / "data"
PARITY_PATH = DATA / "e2e_parity.json"

# Columns of 12_experiment_summary/experiment_blast_results.tsv that the
# comparison reads. saccver (not sseqid) is the versioned accession the rest of
# the pipeline keys on.
SAMPLE_COLUMN = "sample"
TAXID_COLUMN = "adjusted_taxid"
REFERENCE_COLUMN = "saccver"


def summarize_detections(
    rows: Iterable[Mapping[str, str]],
) -> dict[str, dict[str, list[str]]]:
    """Reduce blast result rows to the detections each sample made.

    Thousands of hit rows collapse to a handful of taxids per sample, which is
    the level the parity check cares about. Blank values are dropped because a
    sample with no hit still contributes a row with empty taxonomy fields.
    """
    taxids_by_sample: dict[str, set[str]] = {}
    references_by_sample: dict[str, set[str]] = {}

    for row in rows:
        sample_id = (row.get(SAMPLE_COLUMN) or "").strip()
        if not sample_id:
            continue
        taxids_by_sample.setdefault(sample_id, set())
        references_by_sample.setdefault(sample_id, set())

        taxid = (row.get(TAXID_COLUMN) or "").strip()
        if taxid:
            taxids_by_sample[sample_id].add(taxid)

        reference = (row.get(REFERENCE_COLUMN) or "").strip()
        if reference:
            references_by_sample[sample_id].add(reference)

    return {
        sample_id: {
            "taxids": sorted(taxids_by_sample[sample_id]),
            "references": sorted(references_by_sample[sample_id]),
        }
        for sample_id in sorted(taxids_by_sample)
    }


def find_parity_violations(
    small: Mapping[str, Mapping[str, list[str]]],
    full: Mapping[str, Mapping[str, list[str]]],
    sample_map: Mapping[str, str],
) -> tuple[list[str], list[str]]:
    """Check a reduced run against a recorded full run.

    ``sample_map`` maps a sample id in the reduced run to its counterpart in the
    recording, and defines the comparison's scope: samples outside it (the
    synthetic metagenomes, which have no full-size counterpart and carry their
    own assertions) must not be passed in.

    Returns ``(errors, warnings)`` as human-readable messages.
    """
    errors: list[str] = []
    warnings: list[str] = []

    errors.extend(
        f"Sample {small_sample_id!r} has no entry in the parity sample map, so "
        f"there is nothing to compare it against. Add a mapping or keep it out "
        f"of the parity check."
        for small_sample_id in sorted(small)
        if small_sample_id not in sample_map
    )

    for small_sample_id, full_sample_id in sorted(sample_map.items()):
        if full_sample_id not in full:
            errors.append(
                f"Sample {small_sample_id!r} maps to {full_sample_id!r}, which the "
                f"recorded full run does not contain.",
            )
            continue
        if small_sample_id not in small:
            errors.append(
                f"Sample {small_sample_id!r} produced no rows, so it cannot be "
                f"compared against {full_sample_id!r}.",
            )
            continue

        extra_taxids = set(small[small_sample_id]["taxids"]) - set(
            full[full_sample_id]["taxids"],
        )
        if extra_taxids:
            errors.append(
                f"Sample {small_sample_id!r} detected taxid(s) "
                f"{sorted(extra_taxids)} that the full run ({full_sample_id!r}) "
                f"never detected.",
            )

        extra_references = set(small[small_sample_id]["references"]) - set(
            full[full_sample_id]["references"],
        )
        if extra_references:
            warnings.append(
                f"Sample {small_sample_id!r} hit reference(s) "
                f"{sorted(extra_references)} absent from the full run "
                f"({full_sample_id!r}). Expected after a BLAST database rebuild as "
                f"long as the taxids match.",
            )

    return errors, warnings


def load_parity(path: Path = PARITY_PATH) -> dict[str, Any]:
    """Read the recorded full-run summary."""
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)
