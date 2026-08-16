#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["loguru", "slack_sdk"]
# ///

"""Send one Slack notification lineage for a rapid-screening sample."""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import sys
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import TYPE_CHECKING, Literal, TypeAlias, cast

from loguru import logger
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError, SlackClientError
from slack_sdk.http_retry.builtin_handlers import RateLimitErrorRetryHandler

if TYPE_CHECKING:
    from collections.abc import Sequence

    from slack_sdk.web import SlackResponse

RiskGroup: TypeAlias = Literal["RG1", "RG2", "RG3", "RG4"]
MinimumRiskGroup: TypeAlias = Literal["RG2", "RG3", "RG4"]
RISK_GROUP_PRIORITY = {"RG1": 1, "RG2": 2, "RG3": 3, "RG4": 4}
NOTIFICATION_MINIMUMS = ("RG2", "RG3", "RG4")
RISK_GROUP_MARKER = {"RG1": "⚪", "RG2": "🟡", "RG3": "🟠", "RG4": "🔴"}
EVIDENCE_BROWSER_ISSUE_URL = "https://github.com/dholab/nvd/issues/NNN"
THREAD_ASSIGNMENTS_PER_PAGE = 48
SlackBlock: TypeAlias = dict[str, object]


@dataclass(frozen=True)
class EnrichmentContext:
    sequences_input: int
    sequences_retained: int
    sequences_retained_proportion: Decimal
    bases_input: int
    bases_retained: int
    bases_retained_proportion: Decimal


@dataclass(frozen=True)
class TaxonomicAssignment:
    taxon_name: str
    lineage: str
    rank: str
    weighted_fraction: Decimal
    estimated_distinct_bp: int


@dataclass(frozen=True)
class RiskContext:
    risk_group: RiskGroup
    source_taxid: int
    source_name: str
    assignments: tuple[TaxonomicAssignment, ...]


@dataclass(frozen=True)
class RapidAlert:
    sample_id: str
    run_id: str
    target_enrichment_enabled: bool
    enrichment: EnrichmentContext | None
    risk_contexts: tuple[RiskContext, ...]


@dataclass(frozen=True)
class AlertRequest:
    sample_id: str
    run_id: str
    risk_summary_path: Path
    enrichment_path: Path
    target_enrichment_enabled: bool
    minimum_risk_group: MinimumRiskGroup


@dataclass(frozen=True)
class SlackMessage:
    text: str
    blocks: tuple[SlackBlock, ...]


@dataclass(frozen=True)
class NotificationLineage:
    parent: SlackMessage
    replies: tuple[SlackMessage, ...]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Send one rapid-screening Slack notification lineage.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--channel", required=True)
    parser.add_argument("--risk-summary", type=Path, required=True)
    parser.add_argument("--enrichment", type=Path, required=True)
    parser.add_argument("--target-enrichment-enabled", action="store_true")
    parser.add_argument(
        "--minimum-risk-group",
        choices=NOTIFICATION_MINIMUMS,
        required=True,
    )
    parser.add_argument("-v", "--verbose", action="count", default=0)
    return parser.parse_args(argv)


def configure_logging(verbosity: int) -> None:
    logger.remove()
    level = "WARNING" if verbosity == 0 else "INFO" if verbosity == 1 else "DEBUG"
    logger.add(sys.stderr, level=level, format="<level>{level: <8}</level> | {message}")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_enrichment(path: Path) -> EnrichmentContext:
    summary = json.loads(path.read_text(encoding="utf-8"))
    return EnrichmentContext(
        sequences_input=int(summary["seqs_in"]),
        sequences_retained=int(summary["seqs_out"]),
        sequences_retained_proportion=Decimal(str(summary["seqs_out_proportion"])),
        bases_input=int(summary["bp_in"]),
        bases_retained=int(summary["bp_out"]),
        bases_retained_proportion=Decimal(str(summary["bp_out_proportion"])),
    )


def parse_risk_group(value: str) -> RiskGroup:
    if value not in RISK_GROUP_PRIORITY:
        message = f"Unsupported rapid-screening risk group: {value!r}"
        raise ValueError(message)
    return cast("RiskGroup", value)


def parse_minimum_risk_group(value: str) -> MinimumRiskGroup:
    if value not in NOTIFICATION_MINIMUMS:
        message = f"Unsupported minimum notification risk group: {value!r}"
        raise ValueError(message)
    return cast("MinimumRiskGroup", value)


def build_alert(request: AlertRequest) -> RapidAlert | None:
    """Build a complete alert model, or None when no configured rows qualify."""
    minimum_priority = RISK_GROUP_PRIORITY[request.minimum_risk_group]
    grouped: dict[
        tuple[RiskGroup, int, str],
        list[TaxonomicAssignment],
    ] = {}
    for row in read_csv_rows(request.risk_summary_path):
        risk_group_text = (row.get("who_risk_group") or "").strip()
        if not risk_group_text:
            continue
        risk_group = parse_risk_group(risk_group_text)
        if RISK_GROUP_PRIORITY[risk_group] < minimum_priority:
            continue
        source_taxid = int(row["who_risk_group_source_taxid"])
        source_name = row["who_risk_group_source_name"].strip()
        lineage = row["lineage"].strip()
        rank = row["rank"].strip()
        weighted_fraction = Decimal(row["f_weighted_at_rank"])
        estimated_distinct_bp = int(row["bp_match_at_rank"])
        taxon_name = lineage.rsplit(";", maxsplit=1)[-1].strip()
        grouped.setdefault((risk_group, source_taxid, source_name), []).append(
            TaxonomicAssignment(
                taxon_name=taxon_name,
                lineage=lineage,
                rank=rank,
                weighted_fraction=weighted_fraction,
                estimated_distinct_bp=estimated_distinct_bp,
            ),
        )

    if not grouped:
        return None

    contexts: list[RiskContext] = []
    for (risk_group, source_taxid, source_name), assignments in grouped.items():
        assignments.sort(
            key=lambda assignment: (
                -assignment.weighted_fraction,
                -assignment.estimated_distinct_bp,
                assignment.taxon_name,
            ),
        )
        contexts.append(
            RiskContext(
                risk_group=risk_group,
                source_taxid=source_taxid,
                source_name=source_name,
                assignments=tuple(assignments),
            ),
        )
    contexts.sort(
        key=lambda context: (
            -RISK_GROUP_PRIORITY[context.risk_group],
            -max(assignment.weighted_fraction for assignment in context.assignments),
            context.source_name,
            context.source_taxid,
        ),
    )
    enrichment = (
        read_enrichment(request.enrichment_path)
        if request.target_enrichment_enabled
        else None
    )
    return RapidAlert(
        sample_id=request.sample_id,
        run_id=request.run_id,
        target_enrichment_enabled=request.target_enrichment_enabled,
        enrichment=enrichment,
        risk_contexts=tuple(contexts),
    )


def escape_slack_text(value: str) -> str:
    return html.escape(value, quote=False)


def format_percentage(fraction: Decimal) -> str:
    return f"{format((fraction * 100).normalize(), 'f')}%"


def plural(count: int, singular: str) -> str:
    return f"{count} {singular if count == 1 else f'{singular}s'}"


def section(text: str) -> SlackBlock:
    return {
        "type": "section",
        "text": {"type": "mrkdwn", "text": text, "verbatim": True},
    }


def render_parent(alert: RapidAlert) -> SlackMessage:
    highest_risk_group = alert.risk_contexts[0].risk_group
    context_counts = Counter(context.risk_group for context in alert.risk_contexts)
    assignment_counts = Counter(
        context.risk_group
        for context in alert.risk_contexts
        for _assignment in context.assignments
    )
    count_lines = [
        (
            f"*{risk_group}:* "
            f"{plural(context_counts[risk_group], 'risk context')} · "
            f"{plural(assignment_counts[risk_group], 'taxonomic assignment')}"
        )
        for risk_group in ("RG4", "RG3", "RG2")
        if context_counts[risk_group]
    ]

    if alert.target_enrichment_enabled and alert.enrichment is not None:
        enrichment = alert.enrichment
        enrichment_text = "\n".join(
            [
                "*Target enrichment*",
                "Reads retained: "
                f"*{enrichment.sequences_retained:,} / "
                f"{enrichment.sequences_input:,} "
                f"({format_percentage(enrichment.sequences_retained_proportion)})*",
                "Bases retained: "
                f"*{enrichment.bases_retained:,} bp / "
                f"{enrichment.bases_input:,} bp "
                f"({format_percentage(enrichment.bases_retained_proportion)})*",
            ],
        )
    elif alert.target_enrichment_enabled:
        enrichment_text = "*Target enrichment:* statistics unavailable"
    else:
        enrichment_text = "*Target enrichment:* not enabled"

    sample_id = escape_slack_text(alert.sample_id)
    run_id = escape_slack_text(alert.run_id)
    blocks: tuple[SlackBlock, ...] = (
        {
            "type": "header",
            "text": {
                "type": "plain_text",
                "text": (
                    f"{RISK_GROUP_MARKER[highest_risk_group]} {highest_risk_group} · "
                    "Preliminary rapid-screening alert"
                ),
            },
        },
        {
            "type": "section",
            "fields": [
                {"type": "mrkdwn", "text": f"*Sample*\n`{sample_id}`"},
                {"type": "mrkdwn", "text": f"*Run*\n`{run_id}`"},
            ],
        },
        section("\n".join(count_lines)),
        section(enrichment_text),
        {
            "type": "context",
            "elements": [
                {
                    "type": "mrkdwn",
                    "text": (
                        "Maximum-sensitivity k-mer screen. Later BLAST/LCA analysis "
                        "remains in progress. Complete screening evidence follows in "
                        "this thread."
                    ),
                },
                {
                    "type": "mrkdwn",
                    "text": (
                        f"<{EVIDENCE_BROWSER_ISSUE_URL}|Direct evidence links are "
                        "planned; follow dholab/nvd#NNN for updates.>"
                    ),
                },
            ],
        },
    )
    return SlackMessage(
        text=(
            f"{highest_risk_group} preliminary rapid-screening alert for sample "
            f"{alert.sample_id} in run {alert.run_id}; complete evidence is in the thread."
        ),
        blocks=blocks,
    )


def render_assignment(
    context: RiskContext,
    assignment: TaxonomicAssignment,
) -> SlackBlock:
    taxon_name = escape_slack_text(assignment.taxon_name)
    rank = escape_slack_text(assignment.rank)
    source_name = escape_slack_text(context.source_name)
    lineage = escape_slack_text(assignment.lineage)
    return section(
        "\n".join(
            [
                f"*{taxon_name}*",
                f"Assignment rank: *{rank}*",
                (
                    f"Risk-group source: *{context.risk_group} · {source_name}* "
                    f"(TaxID `{context.source_taxid}`)"
                ),
                (
                    "Assigned sequence-data share: "
                    f"*{format_percentage(assignment.weighted_fraction)}*"
                ),
                (
                    "Estimated distinct support: "
                    f"*{assignment.estimated_distinct_bp:,} bp*"
                ),
                f"Lineage: {lineage}",
            ],
        ),
    )


def render_replies(alert: RapidAlert) -> tuple[SlackMessage, ...]:
    replies: list[SlackMessage] = []
    for risk_group in ("RG4", "RG3", "RG2"):
        assignments = [
            (context, assignment)
            for context in alert.risk_contexts
            if context.risk_group == risk_group
            for assignment in context.assignments
        ]
        if not assignments:
            continue
        pages = [
            assignments[index : index + THREAD_ASSIGNMENTS_PER_PAGE]
            for index in range(0, len(assignments), THREAD_ASSIGNMENTS_PER_PAGE)
        ]
        for page_index, page in enumerate(pages, 1):
            page_suffix = (
                f" · page {page_index} of {len(pages)}" if len(pages) > 1 else ""
            )
            blocks: tuple[SlackBlock, ...] = (
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": f"{risk_group}-associated screening evidence{page_suffix}",
                    },
                },
                *(
                    render_assignment(context, assignment)
                    for context, assignment in page
                ),
            )
            replies.append(
                SlackMessage(
                    text=(
                        f"{risk_group}-associated screening evidence for "
                        f"{alert.sample_id}{page_suffix}"
                    ),
                    blocks=blocks,
                ),
            )
    return tuple(replies)


def render_lineage(alert: RapidAlert) -> NotificationLineage:
    return NotificationLineage(
        parent=render_parent(alert),
        replies=render_replies(alert),
    )


def post_message(
    client: WebClient,
    *,
    channel: str,
    message: SlackMessage,
    thread_ts: str | None = None,
) -> SlackResponse:
    if thread_ts is None:
        return client.chat_postMessage(
            channel=channel,
            text=message.text,
            blocks=list(message.blocks),
            unfurl_links=False,
            unfurl_media=False,
        )
    return client.chat_postMessage(
        channel=channel,
        thread_ts=thread_ts,
        text=message.text,
        blocks=list(message.blocks),
        unfurl_links=False,
        unfurl_media=False,
    )


def post_lineage(channel: str, lineage: NotificationLineage) -> bool:
    """Post one best-effort parent-and-replies lineage to Slack."""
    token = os.environ.get("SLACK_BOT_TOKEN")
    if not token:
        logger.warning("SLACK_BOT_TOKEN not set; skipping rapid-screening notification")
        return False

    client = WebClient(token=token)
    client.retry_handlers.append(RateLimitErrorRetryHandler(max_retry_count=3))
    try:
        parent = post_message(client, channel=channel, message=lineage.parent)
        if not parent.get("ok", False):
            logger.warning("Slack API returned a non-ok parent response: {}", parent)
            return False
        parent_channel = str(parent["channel"])
        parent_timestamp = str(parent["ts"])
        for reply in lineage.replies:
            response = post_message(
                client,
                channel=parent_channel,
                message=reply,
                thread_ts=parent_timestamp,
            )
            if not response.get("ok", False):
                logger.warning(
                    "Slack API returned a non-ok reply response: {}",
                    response,
                )
                return False
    except (SlackApiError, SlackClientError, OSError) as error:
        logger.warning("Rapid-screening Slack notification failed: {}", error)
        return False
    return True


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    configure_logging(args.verbose)
    alert = build_alert(
        AlertRequest(
            sample_id=args.sample_id,
            run_id=args.run_id,
            risk_summary_path=args.risk_summary,
            enrichment_path=args.enrichment,
            target_enrichment_enabled=args.target_enrichment_enabled,
            minimum_risk_group=parse_minimum_risk_group(args.minimum_risk_group),
        ),
    )
    if alert is None:
        logger.info(
            "No rapid-screening assignments met the {} notification policy for {}",
            args.minimum_risk_group,
            args.sample_id,
        )
        return
    post_lineage(args.channel, render_lineage(alert))


if __name__ == "__main__":
    main()
