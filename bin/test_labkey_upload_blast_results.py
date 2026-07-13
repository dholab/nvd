"""Tests for preserving taxonomy absence at the LabKey upload boundary."""

import polars as pl
from labkey_upload_blast_results import dataframe_to_records, validate_dataframe


def test_unavailable_staxid_remains_absent_in_labkey_record() -> None:
    logs: list[str] = []
    frame = pl.DataFrame(
        {
            "experiment": [1],
            "sample_id": ["sample-1"],
            "qseqid": ["query-1"],
            "staxids": [None],
            "adjusted_taxid": [9606],
        },
        schema_overrides={"staxids": pl.Int64},
    )

    [record] = dataframe_to_records(validate_dataframe(frame, "blast.csv", logs))

    assert record["staxids"] == ""
    assert record["staxids"] != 0
