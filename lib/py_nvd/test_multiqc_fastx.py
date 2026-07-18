"""Tests for typed FASTX report profiles."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import yaml

from py_nvd.multiqc_domains import LineGraphSection, build_domain_sections
from py_nvd.multiqc_fastx import (
    FASTX_MAX_POINTS,
    InvalidFastxProfile,
    ParsedFastxProfile,
    collect_fastx_profiles,
    histogram_series,
    rebin_series,
    series_key,
)
from py_nvd.multiqc_packages import Domain, ProfileReceipt, ReportPackage, write_package
from py_nvd.multiqc_report import write_domain_sections

if TYPE_CHECKING:
    from pathlib import Path


def test_profile_parses_once_and_histogram_is_rebinned_for_presentation(
    tmp_path: Path,
) -> None:
    profile = tmp_path / "profile.json"
    profile.write_text(
        json.dumps(
            {
                "sample_id": "sample_A",
                "stage": "post_depletion",
                "format": "fastq",
                "sequence_count": 1001,
                "total_bases": 10010,
                "length": {
                    "min": 10,
                    "max": 10,
                    "mean": 10,
                    "n50": 10,
                    "l50": 501,
                },
                "quality": {"mean": 30},
                "thresholds": [],
            },
        ),
        encoding="utf-8",
    )
    histogram = tmp_path / "length.tsv"
    histogram.write_text(
        "sample_id\tstage\tbin_start\tbin_end\tsequence_count\tproducer_note\n"
        + "\n".join(
            f"sample_A\tpost_depletion\t{index}\t{index}\t1\ttrusted"
            for index in range(1001)
        )
        + "\n",
        encoding="utf-8",
    )
    receipt = ProfileReceipt(
        domain=Domain.FASTX,
        sample_id="sample_A",
        stage="post_depletion",
    )
    package_dir = write_package(receipt, (profile, histogram), output_root=tmp_path)
    [parsed_profile] = collect_fastx_profiles(
        (ReportPackage(package_dir, receipt),),
        domain="fastx",
    )

    assert isinstance(parsed_profile, ParsedFastxProfile)
    [series] = histogram_series((parsed_profile,), quality=False).values()
    assert len(series) <= 1000
    assert sum(series.values()) == 1001


def test_malformed_profile_becomes_a_local_invalid_profile(tmp_path: Path) -> None:
    profile = tmp_path / "profile.json"
    profile.write_text("{}", encoding="utf-8")
    histogram = tmp_path / "length.tsv"
    histogram.write_text(
        "sample_id\tstage\tbin_start\tbin_end\tsequence_count\n",
        encoding="utf-8",
    )
    receipt = ProfileReceipt(
        domain=Domain.FASTX,
        sample_id="sample_A",
        stage="post_depletion",
    )
    package_dir = write_package(receipt, (profile, histogram), output_root=tmp_path)

    [parsed_profile] = collect_fastx_profiles(
        (ReportPackage(package_dir, receipt),),
        domain="fastx",
    )

    assert isinstance(parsed_profile, InvalidFastxProfile)


def test_histogram_from_another_sample_becomes_a_local_invalid_profile(
    tmp_path: Path,
) -> None:
    profile = tmp_path / "profile.json"
    profile.write_text(
        json.dumps(
            {
                "sample_id": "sample_A",
                "stage": "post_depletion",
                "format": "fastq",
                "sequence_count": 1,
                "total_bases": 10,
                "length": {"min": 10, "max": 10, "mean": 10, "n50": 10, "l50": 1},
                "quality": {"mean": 30},
                "thresholds": [],
            },
        ),
        encoding="utf-8",
    )
    histogram = tmp_path / "length.tsv"
    histogram.write_text(
        "sample_id\tstage\tbin_start\tbin_end\tsequence_count\n"
        "sample_B\tpost_depletion\t10\t10\t1\n",
        encoding="utf-8",
    )
    receipt = ProfileReceipt(
        domain=Domain.FASTX,
        sample_id="sample_A",
        stage="post_depletion",
    )
    package_dir = write_package(receipt, (profile, histogram), output_root=tmp_path)

    [parsed_profile] = collect_fastx_profiles(
        (ReportPackage(package_dir, receipt),),
        domain="fastx",
    )

    assert isinstance(parsed_profile, InvalidFastxProfile)
    assert (
        str(parsed_profile.error)
        == "FASTX histogram identity does not match its receipt"
    )


def test_rebin_preserves_sparse_counts() -> None:
    points = dict.fromkeys(range(0, 4000, 2), 1)
    rebinned = rebin_series(points)
    assert len(rebinned) <= 1000
    assert sum(rebinned.values()) == sum(points.values())


def test_series_keys_remain_distinct_after_multiqc_sample_name_cleaning() -> None:
    first = ProfileReceipt(
        domain=Domain.FASTX,
        sample_id="sample_A",
        stage="single_read",
        query_class="single_read",
        producer="single_read",
    )
    second = first.model_copy(update={"sample_id": "sample_B"})

    assert series_key(first) == "sample_A | single_read"
    assert series_key(second) == "sample_B | single_read"


def test_domain_sections_split_length_histograms_by_fastx_stage(
    tmp_path: Path,
) -> None:
    packages = []
    for index, stage in enumerate(
        ("single_read", "overlap_merged_pair", "filtered_contigs"),
        1,
    ):
        sample_id = f"sample_{index}"
        profile = tmp_path / f"{stage}.json"
        profile.write_text(
            json.dumps(
                {
                    "sample_id": sample_id,
                    "stage": stage,
                    "sequence_count": 1,
                    "total_bases": 100,
                    "length": {
                        "min": 100,
                        "max": 100,
                        "mean": 100,
                        "n50": 100,
                        "l50": 1,
                    },
                    "quality": {"mean": None},
                    "thresholds": [],
                },
            ),
            encoding="utf-8",
        )
        histogram = tmp_path / f"{stage}.tsv"
        histogram.write_text(
            "sample_id\tstage\tbin_start\tbin_end\tsequence_count\n"
            f"{sample_id}\t{stage}\t100\t149\t1\n",
            encoding="utf-8",
        )
        receipt = ProfileReceipt(
            domain=Domain.FASTX,
            sample_id=sample_id,
            stage=stage,
        )
        package_dir = write_package(
            receipt,
            (profile, histogram),
            output_root=tmp_path,
        )
        packages.append(ReportPackage(package_dir, receipt))

    sections = build_domain_sections(
        packages=tuple(packages),
        sample_ids=(),
        target_enrichment_enabled=False,
        depletion_enabled=False,
        assembly_enabled=False,
    )

    expected = {
        "nvd_fastx_single_read_length_distribution": "sample_1 | single_read",
        "nvd_fastx_overlap_merged_pair_length_distribution": (
            "sample_2 | overlap_merged_pair"
        ),
        "nvd_fastx_filtered_contigs_length_distribution": (
            "sample_3 | filtered_contigs"
        ),
    }
    assert set(sections) == {"nvd_fastx_profiles", *expected}
    for section_id, series_name in expected.items():
        section = sections[section_id]
        assert isinstance(section, LineGraphSection)
        assert section.data == {series_name: {100: 1}}

    rendered = tmp_path / "rendered"
    rendered.mkdir()
    write_domain_sections(rendered, sections)
    assert not (rendered / "nvd_fastx_length_distribution_mqc.yaml").exists()
    for section_id in expected:
        filename = f"{section_id}_mqc.yaml"
        path = rendered / filename
        assert path.is_file()
        content = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert content["pconfig"] == {
            "xlab": "Length bin start",
            "ylab": "Sequence count",
            "style": "lines",
            "smooth_points": FASTX_MAX_POINTS,
            "xlog": False,
            "ylog": False,
            "logswitch": False,
        }
