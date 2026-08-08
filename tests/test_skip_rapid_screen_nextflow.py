from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_main_workflow_gates_screening_and_its_dependents_together() -> None:
    source = (ROOT / "workflows" / "nvd_main.nf").read_text(encoding="utf-8")

    assert (
        "def rapid_screening_enabled = "
        "params.experimental == true && !params.skip_rapid_screen"
    ) in source
    assert (
        "if (rapid_screening_enabled) {\n    rapid_screening = RAPID_SCREENING("
    ) in source
    assert "SAMPLE_SIMILARITY_QC(rapid_screening.query_sketches)" in source
    assert "if (rapid_screening_enabled) {\n    RAPID_SCREENING_EVAL(" in source


def test_reporting_separates_sourmash_outputs_from_other_experimental_outputs() -> None:
    source = (ROOT / "subworkflows" / "reporting.nf").read_text(encoding="utf-8")

    assert (
        "def rapid_screening_enabled = "
        "params.experimental == true && !params.skip_rapid_screen"
    ) in source
    assert "if (params.experimental == true) {\n        BUILD_SEQUENCE_FLOW" in source
    assert (
        "if (rapid_screening_enabled) {\n        ch_sourmash_profile_summaries"
    ) in source
    assert (
        "if (rapid_screening_enabled) {\n        ch_reporting_terminal_outputs"
    ) in source
    for output_name in (
        "taxon_abundance_sunbursts",
        "merged_taxon_abundance_sunburst",
        "sourmash_sankey_reports",
    ):
        assert f"{output_name} = rapid_screening_enabled ?" in source
