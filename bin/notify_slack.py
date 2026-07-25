#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["loguru", "slack_sdk"]
# ///

"""Send a state-free Slack notification when an NVD run completes."""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import os
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal, TypeAlias, cast
from urllib.parse import quote, urlencode

from loguru import logger
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError, SlackClientError
from slack_sdk.http_retry.builtin_handlers import RateLimitErrorRetryHandler

if TYPE_CHECKING:
    from collections.abc import Sequence

    from slack_sdk.web import SlackResponse

RiskGroup: TypeAlias = Literal["RG1", "RG2", "RG3", "RG4"]
RISK_GROUP_PRIORITY = {"RG1": 1, "RG2": 2, "RG3": 3, "RG4": 4}
NOTIFIED_RISK_GROUPS: tuple[RiskGroup, ...] = ("RG4", "RG3")
THREAD_ASSIGNMENTS_PER_PAGE = 48
SLACK_SECTION_TEXT_LIMIT = 3_000
SLACK_MESSAGE_TEXT_BUDGET = 38_000
SlackBlock: TypeAlias = dict[str, object]


@dataclass(frozen=True)
class CompletionRequest:
    run_id: str
    experiment_id: str | None
    sample_set_id: str
    results_path: Path
    labkey_url: str | None


@dataclass(frozen=True)
class QueryAssignment:
    sample_id: str
    query_id: str
    assigned_taxid: int
    taxon_name: str
    taxon_rank: str
    risk_group: RiskGroup
    risk_source_taxid: int
    risk_source_name: str
    query_class: str
    query_span: int
    assignment_method: str


@dataclass(frozen=True)
class TaxonomicAssignment:
    risk_group: RiskGroup
    risk_source_taxid: int
    risk_source_name: str
    sample_id: str
    assigned_taxid: int
    taxon_name: str
    taxon_rank: str
    supporting_query_count: int
    query_class_counts: tuple[tuple[str, int], ...]
    total_query_span: int
    assignment_methods: tuple[str, ...]


@dataclass(frozen=True)
class CompletionReport:
    run_id: str
    experiment_id: str | None
    sample_set_id: str
    labkey_url: str | None
    assignments: tuple[TaxonomicAssignment, ...]


@dataclass(frozen=True)
class SlackMessage:
    text: str
    blocks: tuple[SlackBlock, ...]


@dataclass(frozen=True)
class NotificationLineage:
    parent: SlackMessage
    replies: tuple[SlackMessage, ...]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Send an NVD completion notification")
    parser.add_argument("--run-id", required=True)
    experiment = parser.add_mutually_exclusive_group()
    experiment.add_argument("--experiment-id")
    experiment.add_argument("--experiment-id-base64", help=argparse.SUPPRESS)
    parser.add_argument("--channel", required=True)
    parser.add_argument("--sample-set-id", required=True)
    parser.add_argument("--experiment-results", type=Path, required=True)
    parser.add_argument("--labkey-server")
    parser.add_argument("--labkey-project")
    parser.add_argument("--labkey-results-list")
    parser.add_argument("--labkey-config-base64", help=argparse.SUPPRESS)
    parser.add_argument("-v", "--verbose", action="count", default=0)
    return parser.parse_args(argv)


def configure_logging(verbosity: int) -> None:
    logger.remove()
    level = "WARNING" if verbosity == 0 else "INFO" if verbosity == 1 else "DEBUG"
    logger.add(sys.stderr, level=level, format="<level>{level: <8}</level> | {message}")


def parse_risk_group(value: str) -> RiskGroup:
    if value not in RISK_GROUP_PRIORITY:
        message = f"Unsupported final-analysis risk group: {value!r}"
        raise ValueError(message)
    return cast("RiskGroup", value)


def build_labkey_experiment_url(
    *,
    server: str,
    project: str,
    results_list: str,
    experiment_id: str,
) -> str:
    """Build the provisional LabKey list-grid URL for one experiment."""
    base_url = server.rstrip("/")
    if not base_url.startswith(("http://", "https://")):
        base_url = f"https://{base_url}"
    project_path = quote(project.strip("/"), safe="/")
    query = urlencode(
        {
            "name": results_list,
            "query.experiment~eq": experiment_id,
        },
    )
    return f"{base_url}/{project_path}/list-grid.view?{query}"


def decode_base64(value: str) -> str:
    return base64.b64decode(value, validate=True).decode("utf-8")


def experiment_id_from_args(args: argparse.Namespace) -> str | None:
    if args.experiment_id_base64:
        return decode_base64(args.experiment_id_base64)
    return args.experiment_id


def optional_labkey_url(
    args: argparse.Namespace,
    *,
    experiment_id: str | None,
) -> str | None:
    if args.labkey_config_base64:
        config = json.loads(decode_base64(args.labkey_config_base64))
        values = (
            config.get("server"),
            config.get("project"),
            config.get("results_list"),
        )
    else:
        values = (args.labkey_server, args.labkey_project, args.labkey_results_list)
    if not any(values):
        return None
    if not all(values):
        message = (
            "LabKey completion links require --labkey-server, --labkey-project, "
            "and --labkey-results-list together"
        )
        raise ValueError(message)
    if experiment_id is None:
        message = "LabKey completion links require an experiment ID"
        raise ValueError(message)
    server, project, results_list = values
    return build_labkey_experiment_url(
        server=server,
        project=project,
        results_list=results_list,
        experiment_id=experiment_id,
    )


def read_result_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def query_assignment(row: dict[str, str]) -> QueryAssignment:
    return QueryAssignment(
        sample_id=row["sample"].strip(),
        query_id=row["qseqid"].strip(),
        assigned_taxid=int(row["adjusted_taxid"]),
        taxon_name=row["adjusted_taxid_name"].strip(),
        taxon_rank=row["adjusted_taxid_rank"].strip(),
        risk_group=parse_risk_group(row["who_risk_group"].strip()),
        risk_source_taxid=int(row["who_risk_group_source_taxid"]),
        risk_source_name=row["who_risk_group_source_name"].strip(),
        query_class=row["query_class"].strip(),
        query_span=int(row["qlen"]),
        assignment_method=row["adjustment_method"].strip(),
    )


def deduplicate_queries(
    rows: Sequence[dict[str, str]],
    *,
    minimum_priority: int,
) -> tuple[QueryAssignment, ...]:
    by_query: dict[tuple[str, str], QueryAssignment] = {}
    for row in rows:
        risk_group_text = (row.get("who_risk_group") or "").strip()
        if not risk_group_text:
            continue
        risk_group = parse_risk_group(risk_group_text)
        if RISK_GROUP_PRIORITY[risk_group] < minimum_priority:
            continue
        assignment = query_assignment(row)
        key = (assignment.sample_id, assignment.query_id)
        previous = by_query.setdefault(key, assignment)
        if previous != assignment:
            message = (
                "Conflicting final taxonomic assignments for sample/query "
                f"{assignment.sample_id!r}/{assignment.query_id!r}"
            )
            raise ValueError(message)
    return tuple(by_query.values())


def build_report(request: CompletionRequest) -> CompletionReport:
    queries = deduplicate_queries(
        read_result_rows(request.results_path),
        minimum_priority=RISK_GROUP_PRIORITY["RG3"],
    )
    grouped: dict[
        tuple[RiskGroup, int, str, str, int, str, str],
        list[QueryAssignment],
    ] = {}
    for query in queries:
        key = (
            query.risk_group,
            query.risk_source_taxid,
            query.risk_source_name,
            query.sample_id,
            query.assigned_taxid,
            query.taxon_name,
            query.taxon_rank,
        )
        grouped.setdefault(key, []).append(query)

    assignments: list[TaxonomicAssignment] = []
    for key, supporting_queries in grouped.items():
        (
            risk_group,
            source_taxid,
            source_name,
            sample_id,
            assigned_taxid,
            taxon_name,
            taxon_rank,
        ) = key
        query_class_counts = Counter(query.query_class for query in supporting_queries)
        assignments.append(
            TaxonomicAssignment(
                risk_group=risk_group,
                risk_source_taxid=source_taxid,
                risk_source_name=source_name,
                sample_id=sample_id,
                assigned_taxid=assigned_taxid,
                taxon_name=taxon_name,
                taxon_rank=taxon_rank,
                supporting_query_count=len(supporting_queries),
                query_class_counts=tuple(sorted(query_class_counts.items())),
                total_query_span=sum(query.query_span for query in supporting_queries),
                assignment_methods=tuple(
                    sorted({query.assignment_method for query in supporting_queries}),
                ),
            ),
        )
    assignments.sort(
        key=lambda assignment: (
            -RISK_GROUP_PRIORITY[assignment.risk_group],
            assignment.risk_source_name,
            assignment.risk_source_taxid,
            assignment.sample_id,
            -assignment.supporting_query_count,
            -assignment.total_query_span,
            assignment.taxon_name,
            assignment.assigned_taxid,
        ),
    )
    return CompletionReport(
        run_id=request.run_id,
        experiment_id=request.experiment_id,
        sample_set_id=request.sample_set_id,
        labkey_url=request.labkey_url,
        assignments=tuple(assignments),
    )


def escape_slack_text(value: str) -> str:
    return html.escape(value, quote=False)


def plural(count: int, singular: str, plural_form: str | None = None) -> str:
    return f"{count} {singular if count == 1 else plural_form or f'{singular}s'}"


def section(text: str) -> SlackBlock:
    if len(text) > SLACK_SECTION_TEXT_LIMIT:
        message = (
            "Slack section text exceeds the 3,000-character limit: "
            f"{len(text):,} characters"
        )
        raise ValueError(message)
    return {
        "type": "section",
        "text": {"type": "mrkdwn", "text": text, "verbatim": True},
    }


def render_parent(report: CompletionReport) -> SlackMessage:
    samples_by_risk: dict[RiskGroup, set[str]] = {
        risk_group: set() for risk_group in NOTIFIED_RISK_GROUPS
    }
    assignments_by_risk: Counter[RiskGroup] = Counter()
    for assignment in report.assignments:
        assignments_by_risk[assignment.risk_group] += 1
        samples_by_risk[assignment.risk_group].add(assignment.sample_id)
    risk_lines = [
        (
            f"{risk_group}: *{plural(len(samples_by_risk[risk_group]), 'sample')} · "
            f"{plural(assignments_by_risk[risk_group], 'taxonomic assignment')}*"
        )
        for risk_group in NOTIFIED_RISK_GROUPS
    ]
    plain_risk_lines = [
        (
            f"{risk_group}: {plural(len(samples_by_risk[risk_group]), 'sample')} · "
            f"{plural(assignments_by_risk[risk_group], 'taxonomic assignment')}"
        )
        for risk_group in NOTIFIED_RISK_GROUPS
    ]
    if not report.assignments:
        no_assignments = (
            "Final BLAST/LCA analysis produced no RG3/RG4-associated assignments."
        )
    else:
        no_assignments = ""

    identity_fields: list[dict[str, str]] = [
        {
            "type": "mrkdwn",
            "text": f"*Run*\n`{escape_slack_text(report.run_id)}`",
        },
    ]
    if report.experiment_id is not None:
        identity_fields.append(
            {
                "type": "mrkdwn",
                "text": (f"*Experiment*\n`{escape_slack_text(report.experiment_id)}`"),
            },
        )
    identity_fields.append(
        {
            "type": "mrkdwn",
            "text": f"*Sample set*\n`{escape_slack_text(report.sample_set_id)}`",
        },
    )
    blocks: list[SlackBlock] = [
        {
            "type": "header",
            "text": {"type": "plain_text", "text": "NVD run complete"},
        },
        {
            "type": "section",
            "fields": identity_fields,
        },
        section("*Final risk-group assignments*\n" + "\n".join(risk_lines)),
    ]
    if no_assignments:
        blocks.append(section(no_assignments))
    if report.labkey_url:
        blocks.append(
            section(
                f"<{report.labkey_url}|View this experiment in LabKey>",
            ),
        )
    blocks.append(
        {
            "type": "context",
            "elements": [
                {
                    "type": "mrkdwn",
                    "text": (
                        "Final BLAST/LCA analysis is available. Taxonomic assignments "
                        "are analytical results and should not be interpreted as "
                        "diagnostic confirmation."
                    ),
                },
            ],
        },
    )
    fallback_lines = [
        "NVD run complete",
        f"Run: {report.run_id}",
        f"Sample set: {report.sample_set_id}",
        "Final risk-group assignments",
        *plain_risk_lines,
    ]
    if report.experiment_id is not None:
        fallback_lines.insert(2, f"Experiment: {report.experiment_id}")
    if no_assignments:
        fallback_lines.append(no_assignments)
    if report.labkey_url:
        fallback_lines.append(f"View this experiment in LabKey: {report.labkey_url}")
    fallback_lines.append(
        "Taxonomic assignments are analytical results and should not be "
        "interpreted as diagnostic confirmation.",
    )
    return SlackMessage(text="\n".join(fallback_lines), blocks=tuple(blocks))


def render_assignment(
    assignment: TaxonomicAssignment,
) -> SlackBlock:
    query_classes = " · ".join(
        f"{escape_slack_text(query_class.replace('_', ' '))} x{count}"
        for query_class, count in assignment.query_class_counts
    )
    methods = ", ".join(
        escape_slack_text(method) for method in assignment.assignment_methods
    )
    return section(
        "\n".join(
            [
                (
                    f"*{escape_slack_text(assignment.sample_id)} · "
                    f"{escape_slack_text(assignment.taxon_name)}*"
                ),
                (
                    f"Assignment: *{escape_slack_text(assignment.taxon_rank)}* · "
                    f"TaxID `{assignment.assigned_taxid}`"
                ),
                (
                    f"Risk-group source: *{assignment.risk_group} · "
                    f"{escape_slack_text(assignment.risk_source_name)}* "
                    f"(TaxID `{assignment.risk_source_taxid}`)"
                ),
                (
                    f"Support: *{plural(assignment.supporting_query_count, 'supporting query', 'supporting queries')}*"
                    f" · {query_classes}"
                ),
                f"Total assigned query span: *{assignment.total_query_span:,} bp*",
                f"Assignment methods: *{methods}*",
            ],
        ),
    )


def assignment_fallback(assignment: TaxonomicAssignment) -> str:
    query_classes = ", ".join(
        f"{query_class.replace('_', ' ')} x{count}"
        for query_class, count in assignment.query_class_counts
    )
    return "\n".join(
        [
            f"{assignment.sample_id} · {assignment.taxon_name}",
            f"Assignment: {assignment.taxon_rank} · TaxID {assignment.assigned_taxid}",
            (
                f"Risk-group source: {assignment.risk_group} · "
                f"{assignment.risk_source_name} "
                f"(TaxID {assignment.risk_source_taxid})"
            ),
            (
                f"Support: {plural(assignment.supporting_query_count, 'supporting query', 'supporting queries')}"
                f" · {query_classes}"
            ),
            f"Total assigned query span: {assignment.total_query_span:,} bp",
            f"Assignment methods: {', '.join(assignment.assignment_methods)}",
        ],
    )


def paginate_assignments(
    assignments: Sequence[TaxonomicAssignment],
) -> tuple[tuple[TaxonomicAssignment, ...], ...]:
    pages: list[tuple[TaxonomicAssignment, ...]] = []
    page: list[TaxonomicAssignment] = []
    fallback_length = 0
    for assignment in assignments:
        assignment_length = len(assignment_fallback(assignment))
        separator_length = 2 if page else 0
        if page and (
            len(page) >= THREAD_ASSIGNMENTS_PER_PAGE
            or fallback_length + separator_length + assignment_length
            > SLACK_MESSAGE_TEXT_BUDGET
        ):
            pages.append(tuple(page))
            page = []
            fallback_length = 0
            separator_length = 0
        page.append(assignment)
        fallback_length += separator_length + assignment_length
    if page:
        pages.append(tuple(page))
    return tuple(pages)


def render_replies(report: CompletionReport) -> tuple[SlackMessage, ...]:
    replies: list[SlackMessage] = []
    for risk_group in NOTIFIED_RISK_GROUPS:
        assignments = [
            assignment
            for assignment in report.assignments
            if assignment.risk_group == risk_group
        ]
        pages = paginate_assignments(assignments)
        for page_index, page in enumerate(pages, 1):
            page_suffix = (
                f" · page {page_index} of {len(pages)}" if len(pages) > 1 else ""
            )
            blocks: tuple[SlackBlock, ...] = (
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": (
                            f"{risk_group}-associated final analysis assignments"
                            f"{page_suffix}"
                        ),
                    },
                },
                *(render_assignment(assignment) for assignment in page),
            )
            title = f"{risk_group}-associated final analysis assignments{page_suffix}"
            replies.append(
                SlackMessage(
                    text=title
                    + "\n\n"
                    + "\n\n".join(
                        assignment_fallback(assignment) for assignment in page
                    ),
                    blocks=blocks,
                ),
            )
    return tuple(replies)


def render_lineage(report: CompletionReport) -> NotificationLineage:
    return NotificationLineage(
        parent=render_parent(report),
        replies=render_replies(report),
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
        logger.warning("SLACK_BOT_TOKEN not set; skipping completion notification")
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
        logger.warning("Completion Slack notification failed: {}", error)
        return False
    return True


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    configure_logging(args.verbose)
    experiment_id = experiment_id_from_args(args)
    report = build_report(
        CompletionRequest(
            run_id=args.run_id,
            experiment_id=experiment_id,
            sample_set_id=args.sample_set_id,
            results_path=args.experiment_results,
            labkey_url=optional_labkey_url(args, experiment_id=experiment_id),
        ),
    )
    post_lineage(args.channel, render_lineage(report))


if __name__ == "__main__":
    main()
