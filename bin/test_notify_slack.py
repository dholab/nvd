"""Behavior tests for the state-free NVD run-completion notification."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

from notify_slack import (
    CompletionRequest,
    NotificationLineage,
    SlackMessage,
    build_labkey_experiment_url,
    build_report,
    post_lineage,
    render_lineage,
)

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


def write_tsv(path: Path, rows: list[dict[str, str]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def rendered_text(blocks: tuple[dict[str, object], ...]) -> str:
    return json.dumps(blocks, sort_keys=True, ensure_ascii=False)


def final_row(**values: str) -> dict[str, str]:
    return {
        "task": "megablast",
        "sample": values["sample"],
        "qseqid": values["qseqid"],
        "qlen": values["qlen"],
        "sseqid": values["sseqid"],
        "adjusted_taxid": values["assigned_taxid"],
        "adjusted_taxid_name": values["assigned_name"],
        "adjusted_taxid_rank": values["assigned_rank"],
        "adjustment_method": values.get("adjustment_method", "dominant"),
        "who_risk_group": values["risk_group"],
        "who_risk_group_source_taxid": values["source_taxid"],
        "who_risk_group_source_name": values["source_name"],
        "query_class": values["query_class"],
    }


def test_completion_report_collapses_hit_rows_and_groups_final_risk_assignments(
    tmp_path: Path,
) -> None:
    results = write_tsv(
        tmp_path / "experiment_blast_results.tsv",
        [
            final_row(
                sample="sample-1",
                qseqid="query-1",
                sseqid="reference-a",
                qlen="600",
                assigned_taxid="123",
                assigned_name="Example virus",
                assigned_rank="species",
                risk_group="RG3",
                source_taxid="122",
                source_name="Example family",
                query_class="short_assembly_contig",
            ),
            # A second retained reference row for the same query is not a second
            # query or taxonomic assignment.
            final_row(
                sample="sample-1",
                qseqid="query-1",
                sseqid="reference-b",
                qlen="600",
                assigned_taxid="123",
                assigned_name="Example virus",
                assigned_rank="species",
                risk_group="RG3",
                source_taxid="122",
                source_name="Example family",
                query_class="short_assembly_contig",
            ),
            final_row(
                sample="sample-1",
                qseqid="query-2",
                sseqid="reference-c",
                qlen="400",
                assigned_taxid="123",
                assigned_name="Example virus",
                assigned_rank="species",
                risk_group="RG3",
                source_taxid="122",
                source_name="Example family",
                query_class="single_read",
                adjustment_method="lca",
            ),
            final_row(
                sample="sample-2",
                qseqid="query-3",
                sseqid="reference-d",
                qlen="50",
                assigned_taxid="999",
                assigned_name="Exceptional virus",
                assigned_rank="species",
                risk_group="RG4",
                source_taxid="999",
                source_name="Exceptional virus",
                query_class="single_read",
            ),
            final_row(
                sample="sample-3",
                qseqid="query-4",
                sseqid="reference-e",
                qlen="500",
                assigned_taxid="222",
                assigned_name="RG2 virus",
                assigned_rank="species",
                risk_group="RG2",
                source_taxid="222",
                source_name="RG2 virus",
                query_class="short_assembly_contig",
            ),
        ],
    )
    request = CompletionRequest(
        run_id="blue_koala",
        experiment_id="12345",
        sample_set_id="abc123",
        results_path=results,
        labkey_url="https://example.org/experiment/12345",
    )

    report = build_report(request)

    assert [assignment.risk_group for assignment in report.assignments] == [
        "RG4",
        "RG3",
    ]
    rg3_assignment = report.assignments[1]
    assert rg3_assignment.sample_id == "sample-1"
    assert rg3_assignment.taxon_name == "Example virus"
    expected_supporting_queries = 2
    expected_query_span = 1000
    assert rg3_assignment.supporting_query_count == expected_supporting_queries
    assert rg3_assignment.total_query_span == expected_query_span
    assert dict(rg3_assignment.query_class_counts) == {
        "short_assembly_contig": 1,
        "single_read": 1,
    }
    assert rg3_assignment.assignment_methods == ("dominant", "lca")

    lineage = render_lineage(report)
    parent_text = rendered_text(lineage.parent.blocks)
    assert "NVD run complete" in parent_text
    assert "✅" not in parent_text
    assert "risk context" not in parent_text.lower()
    assert "RG4: *1 sample · 1 taxonomic assignment*" in parent_text
    assert "RG3: *1 sample · 1 taxonomic assignment*" in parent_text
    assert "View this experiment in LabKey" in parent_text
    assert "RG4: 1 sample · 1 taxonomic assignment" in lineage.parent.text
    reply_text = " ".join(rendered_text(reply.blocks) for reply in lineage.replies)
    assert "Example family" in reply_text
    assert "Example virus" in reply_text
    assert "2 supporting queries" in reply_text
    assert "short assembly contig x1" in reply_text
    assert "single read x1" in reply_text
    assert "1,000 bp" in reply_text
    assert "dominant, lca" in reply_text
    assert "Example family" in " ".join(reply.text for reply in lineage.replies)
    assert "2 supporting queries" in " ".join(reply.text for reply in lineage.replies)


def test_completion_without_priority_assignments_still_sends_parent(
    tmp_path: Path,
) -> None:
    results = write_tsv(
        tmp_path / "experiment_blast_results.tsv",
        [
            final_row(
                sample="sample-1",
                qseqid="query-1",
                sseqid="reference-a",
                qlen="500",
                assigned_taxid="222",
                assigned_name="RG2 virus",
                assigned_rank="species",
                risk_group="RG2",
                source_taxid="222",
                source_name="RG2 virus",
                query_class="short_assembly_contig",
            ),
        ],
    )
    report = build_report(
        CompletionRequest(
            run_id="blue_koala",
            experiment_id=None,
            sample_set_id="abc123",
            results_path=results,
            labkey_url=None,
        ),
    )

    lineage = render_lineage(report)

    assert not report.assignments
    assert not lineage.replies
    parent_text = rendered_text(lineage.parent.blocks)
    assert "no RG3/RG4-associated assignments" in parent_text
    assert "LabKey" not in parent_text
    assert "Experiment" not in parent_text
    assert "None" not in lineage.parent.text


def test_final_assignments_are_paginated_without_omissions(tmp_path: Path) -> None:
    results = write_tsv(
        tmp_path / "experiment_blast_results.tsv",
        [
            final_row(
                sample="sample-1",
                qseqid=f"query-{index:02d}",
                sseqid=f"reference-{index:02d}",
                qlen="50",
                assigned_taxid=str(10_000 + index),
                assigned_name=f"Signal {index:02d}",
                assigned_rank="species",
                risk_group="RG3",
                source_taxid="122",
                source_name="Example family",
                query_class="single_read",
            )
            for index in range(49)
        ],
    )
    report = build_report(
        CompletionRequest(
            run_id="blue_koala",
            experiment_id="12345",
            sample_set_id="abc123",
            results_path=results,
            labkey_url=None,
        ),
    )

    lineage = render_lineage(report)

    expected_reply_count = 2
    slack_block_limit = 50
    assert len(lineage.replies) == expected_reply_count
    assert all(len(reply.blocks) <= slack_block_limit for reply in lineage.replies)
    reply_text = " ".join(rendered_text(reply.blocks) for reply in lineage.replies)
    for index in range(49):
        assert reply_text.count(f"Signal {index:02d}") == 1


def test_final_assignments_paginate_accessible_fallback_text(tmp_path: Path) -> None:
    results = write_tsv(
        tmp_path / "experiment_blast_results.tsv",
        [
            final_row(
                sample="sample-1",
                qseqid=f"query-{index:02d}",
                sseqid=f"reference-{index:02d}",
                qlen="50",
                assigned_taxid=str(20_000 + index),
                assigned_name=f"Signal {index:02d} " + "x" * 1_400,
                assigned_rank="species",
                risk_group="RG3",
                source_taxid="122",
                source_name="Example family",
                query_class="single_read",
            )
            for index in range(30)
        ],
    )
    report = build_report(
        CompletionRequest(
            run_id="blue_koala",
            experiment_id=None,
            sample_set_id="abc123",
            results_path=results,
            labkey_url=None,
        ),
    )

    lineage = render_lineage(report)

    expected_reply_count = 2
    slack_message_text_limit = 40_000
    assert len(lineage.replies) == expected_reply_count
    assert all(len(reply.text) < slack_message_text_limit for reply in lineage.replies)


def test_labkey_experiment_url_uses_the_uploaded_experiment_field() -> None:
    url = build_labkey_experiment_url(
        server="labkey.example.org",
        project="dho/projects/Infinite Path",
        results_list="BLAST meta hits",
        experiment_id="12345",
    )

    assert url == (
        "https://labkey.example.org/dho/projects/Infinite%20Path/list-grid.view?"
        "name=BLAST+meta+hits&query.experiment~eq=12345"
    )


def test_post_lineage_posts_parent_then_thread_replies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("SLACK_BOT_TOKEN", "xoxb-test")
    client = MagicMock()
    client.chat_postMessage.return_value = {
        "ok": True,
        "channel": "C123",
        "ts": "123.456",
    }
    lineage = NotificationLineage(
        parent=SlackMessage(text="parent fallback", blocks=()),
        replies=(SlackMessage(text="reply fallback", blocks=()),),
    )

    with patch("notify_slack.WebClient", return_value=client):
        assert post_lineage("C123", lineage)

    expected_post_count = 2
    assert client.chat_postMessage.call_count == expected_post_count
    assert client.chat_postMessage.call_args_list[0].kwargs == {
        "channel": "C123",
        "text": "parent fallback",
        "blocks": [],
        "unfurl_links": False,
        "unfurl_media": False,
    }
    assert client.chat_postMessage.call_args_list[1].kwargs == {
        "channel": "C123",
        "thread_ts": "123.456",
        "text": "reply fallback",
        "blocks": [],
        "unfurl_links": False,
        "unfurl_media": False,
    }


def test_post_lineage_skips_without_token(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SLACK_BOT_TOKEN", raising=False)
    lineage = NotificationLineage(
        parent=SlackMessage(text="parent fallback", blocks=()),
        replies=(),
    )

    assert not post_lineage("C123", lineage)
