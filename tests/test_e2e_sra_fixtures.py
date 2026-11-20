"""The downsampled SRA fixtures must stay committed, intact and documented.

These fixtures are what let the default end-to-end test run offline. Nothing
else notices if one goes missing or drifts: the e2e would simply assert against
different reads, or a contributor would find an undocumented FASTQ in
``tests/data``. These tests make both failures loud.

The full SRA runs are still downloaded and exercised by the ``e2e_full`` tier,
which is what keeps the ``FETCH_FASTQ``/sracha path under test. That division is
deliberate, and the manifest's description records it.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"
MANIFEST_PATH = DATA / "e2e_sra_fixtures.manifest.json"

# Below this, a fixture risks losing its target virus entirely and the e2e would
# assert nothing. The measured assembly floor is around 1,000 viral pairs.
MINIMUM_VIRAL_READS = 1000

# The seed every fixture was sampled with. Changing it changes every checksum.
EXPECTED_DOWNSAMPLE_SEED = 42


def load_manifest() -> dict:
    """Read the SRA fixture provenance manifest."""
    with MANIFEST_PATH.open(encoding="utf-8") as handle:
        return json.load(handle)


def sha256_of_file(path: Path) -> str:
    """Return the hex sha256 of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_every_manifest_fixture_is_committed_and_matches_its_checksum() -> None:
    """A missing or drifted fixture must fail here, not silently change the e2e."""
    for fixture in load_manifest()["fixtures"]:
        for file_entry in fixture["files"]:
            path = DATA / file_entry["name"]
            assert path.is_file(), f"{file_entry['name']} is not committed to tests/data"
            assert path.stat().st_size == file_entry["bytes"], (
                f"{file_entry['name']} is {path.stat().st_size} bytes, "
                f"manifest says {file_entry['bytes']}"
            )
            assert sha256_of_file(path) == file_entry["sha256"], (
                f"{file_entry['name']} does not match its pinned sha256"
            )


def test_every_committed_sra_fixture_is_documented() -> None:
    """An undocumented FASTQ in tests/data is how fixture provenance gets lost.

    The ``water*`` fixtures already have no recorded source accession, which is
    exactly the situation this guards against repeating.
    """
    documented = {
        file_entry["name"]
        for fixture in load_manifest()["fixtures"]
        for file_entry in fixture["files"]
    }
    committed = {path.name for path in DATA.glob("*_s42*.fastq.gz")}
    assert committed == documented, (
        "Committed downsampled fixtures and manifest entries disagree. "
        f"Only on disk: {sorted(committed - documented)}. "
        f"Only in the manifest: {sorted(documented - committed)}."
    )


def test_each_fixture_keeps_enough_viral_reads_to_assert_on() -> None:
    """A fraction lowered too far would leave the e2e asserting against nothing."""
    for fixture in load_manifest()["fixtures"]:
        assert fixture["viral_reads_deacon"] >= MINIMUM_VIRAL_READS, (
            f"{fixture['sample_id']} keeps only {fixture['viral_reads_deacon']} "
            f"viral reads, below the {MINIMUM_VIRAL_READS} needed to assert on"
        )
        assert fixture["viral_reads_deacon"] <= fixture["reads"], (
            f"{fixture['sample_id']} reports more viral reads than total reads"
        )


def test_the_manifest_pins_the_reproducibility_parameters() -> None:
    """Seed, fraction and seqkit version together determine the bytes."""
    manifest = load_manifest()
    assert manifest["seed"] == EXPECTED_DOWNSAMPLE_SEED
    assert manifest["tool_versions"]["seqkit"] == "seqkit v2.13.0"
    for fixture in manifest["fixtures"]:
        assert 0 < fixture["fraction"] <= 1, (
            f"{fixture['sample_id']} has an out-of-range fraction"
        )
        assert fixture["source_accession"], (
            f"{fixture['sample_id']} has no source accession to re-create it from"
        )


def test_the_long_read_code_path_is_covered_offline() -> None:
    """Losing the only ONT fixture would silently drop long-read assembly coverage.

    The synthetic metagenomes are paired Illumina, so this fixture is the only
    thing exercising the long-read route in the offline tier.
    """
    platforms = {fixture["platform"] for fixture in load_manifest()["fixtures"]}
    assert "ont" in platforms, (
        "No ONT fixture is committed, so the offline e2e covers no long-read path."
    )
