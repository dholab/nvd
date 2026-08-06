"""Behavior tests for rapid-screening Slack notifications."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

from notify_rapid_screening_slack import (
    AlertRequest,
    NotificationLineage,
    SlackMessage,
    build_alert,
    parse_args,
    post_lineage,
    render_lineage,
)

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


def write_csv(path: Path, rows: list[dict[str, str]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def rendered_text(blocks: tuple[dict[str, object], ...]) -> str:
    return json.dumps(blocks, sort_keys=True, ensure_ascii=False)


def test_build_alert_groups_all_configured_risk_assignments_and_orders_evidence(
    tmp_path: Path,
) -> None:
    risk_summary = write_csv(
        tmp_path / "sample-1.with-risk-groups.csv",
        [
            {
                "rank": "species",
                "lineage": "Viruses;RG3 source;Lower share",
                "fraction": "0.01",
                "f_weighted_at_rank": "0.1",
                "bp_match_at_rank": "50",
                "who_risk_group": "RG3",
                "who_risk_group_source_taxid": "123",
                "who_risk_group_source_name": "RG3 source",
            },
            {
                "rank": "species",
                "lineage": "Viruses;RG4 source;Single hash signal",
                "fraction": "0.0002",
                "f_weighted_at_rank": "0.001",
                "bp_match_at_rank": "50",
                "who_risk_group": "RG4",
                "who_risk_group_source_taxid": "999",
                "who_risk_group_source_name": "RG4 source",
            },
            {
                "rank": "species",
                "lineage": "Viruses;RG3 source;Higher share",
                "fraction": "0.02",
                "f_weighted_at_rank": "0.2",
                "bp_match_at_rank": "100",
                "who_risk_group": "RG3",
                "who_risk_group_source_taxid": "123",
                "who_risk_group_source_name": "RG3 source",
            },
            {
                "rank": "species",
                "lineage": "Viruses;RG2 source;Excluded by current policy",
                "fraction": "0.4",
                "f_weighted_at_rank": "0.5",
                "bp_match_at_rank": "2000",
                "who_risk_group": "RG2",
                "who_risk_group_source_taxid": "222",
                "who_risk_group_source_name": "RG2 source",
            },
            {
                "rank": "species",
                "lineage": "Viruses;RG1 source;Valid non-alerting context",
                "fraction": "0.1",
                "f_weighted_at_rank": "0.1",
                "bp_match_at_rank": "500",
                "who_risk_group": "RG1",
                "who_risk_group_source_taxid": "111",
                "who_risk_group_source_name": "RG1 source",
            },
        ],
    )
    enrichment = tmp_path / "sample-1.deacon_filter.json"
    enrichment.write_text(
        json.dumps(
            {
                "seqs_in": 100_000,
                "seqs_out": 125,
                "seqs_out_proportion": 0.00125,
                "bp_in": 15_000_000,
                "bp_out": 18_750,
                "bp_out_proportion": 0.00125,
            },
        ),
        encoding="utf-8",
    )

    alert = build_alert(
        AlertRequest(
            sample_id="sample-1",
            run_id="blue_koala",
            risk_summary_path=risk_summary,
            enrichment_path=enrichment,
            target_enrichment_enabled=True,
            minimum_risk_group="RG3",
        ),
    )

    assert alert is not None
    assert alert.enrichment is not None
    expected_retained_sequences = 125
    assert alert.enrichment.sequences_retained == expected_retained_sequences
    assert [context.risk_group for context in alert.risk_contexts] == ["RG4", "RG3"]
    assert [context.source_name for context in alert.risk_contexts] == [
        "RG4 source",
        "RG3 source",
    ]
    assert [
        assignment.taxon_name for assignment in alert.risk_contexts[1].assignments
    ] == ["Higher share", "Lower share"]
    expected_single_hash_bp = 50
    assert (
        alert.risk_contexts[0].assignments[0].estimated_distinct_bp
        == expected_single_hash_bp
    )
    expected_distinct_bp = 100
    assert (
        alert.risk_contexts[1].assignments[0].estimated_distinct_bp
        == expected_distinct_bp
    )

    lineage = render_lineage(alert)
    parent_text = rendered_text(lineage.parent.blocks)
    assert "🔴 RG4 · Preliminary rapid-screening alert" in parent_text
    assert "sample-1" in parent_text
    assert "blue_koala" in parent_text
    assert "1 risk context · 1 taxonomic assignment" in parent_text
    assert "1 risk context · 2 taxonomic assignments" in parent_text
    assert "125 / 100,000" in parent_text
    assert "https://github.com/dholab/nvd/issues/NNN" in parent_text
    assert [reply.text.split(maxsplit=1)[0] for reply in lineage.replies] == [
        "RG4-associated",
        "RG3-associated",
    ]
    reply_text = " ".join(rendered_text(reply.blocks) for reply in lineage.replies)
    assert "Single hash signal" in reply_text
    assert "Estimated distinct support: *50 bp*" in reply_text
    assert "Higher share" in reply_text
    assert "Lower share" in reply_text


def test_render_lineage_paginates_without_omitting_assignments(
    tmp_path: Path,
) -> None:
    risk_summary = write_csv(
        tmp_path / "sample-2.with-risk-groups.csv",
        [
            {
                "rank": "species",
                "lineage": f"Viruses;RG3 source;Signal {index:02d}",
                "fraction": "0.00002",
                "f_weighted_at_rank": "0.00002",
                "bp_match_at_rank": "50",
                "who_risk_group": "RG3",
                "who_risk_group_source_taxid": "123",
                "who_risk_group_source_name": "RG3 source",
            }
            for index in range(49)
        ],
    )
    request = AlertRequest(
        sample_id="sample-2",
        run_id="blue_koala",
        risk_summary_path=risk_summary,
        enrichment_path=tmp_path / "unused-enrichment.json",
        target_enrichment_enabled=False,
        minimum_risk_group="RG3",
    )

    alert = build_alert(request)

    assert alert is not None
    assert alert.enrichment is None
    lineage = render_lineage(alert)
    assert "Target enrichment:* not enabled" in rendered_text(lineage.parent.blocks)
    expected_reply_count = 2
    slack_block_limit = 50
    assert len(lineage.replies) == expected_reply_count
    assert all(len(reply.blocks) <= slack_block_limit for reply in lineage.replies)
    reply_text = " ".join(rendered_text(reply.blocks) for reply in lineage.replies)
    for index in range(49):
        assert reply_text.count(f"*Signal {index:02d}*") == 1

    no_alert = build_alert(
        AlertRequest(
            sample_id=request.sample_id,
            run_id=request.run_id,
            risk_summary_path=request.risk_summary_path,
            enrichment_path=request.enrichment_path,
            target_enrichment_enabled=request.target_enrichment_enabled,
            minimum_risk_group="RG4",
        ),
    )
    assert no_alert is None


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

    with patch("notify_rapid_screening_slack.WebClient", return_value=client):
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


def test_cli_accepts_rg2_as_the_minimum_risk_group(tmp_path: Path) -> None:
    placeholder = tmp_path / "placeholder"

    args = parse_args(
        [
            "--sample-id",
            "sample-1",
            "--run-id",
            "blue_koala",
            "--channel",
            "C123",
            "--risk-summary",
            str(placeholder),
            "--enrichment",
            str(placeholder),
            "--minimum-risk-group",
            "RG2",
        ],
    )

    assert args.minimum_risk_group == "RG2"
