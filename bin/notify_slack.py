#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
#     "slack_sdk",
# ]
# ///

"""
Slack notification for NVD run completion.

Sends a formatted message to a Slack channel with run statistics.
Designed to fail gracefully (warn but don't exit non-zero).

Usage:
    notify_slack.py \
        --run-id "run_2026-01-10_abc123" \
        --experiment-id "exp_2026_001" \
        --channel "C0123456789" \
        --state-dir /path/to/state \
        --sample-set-id "abc123def456" \
        --labkey-url "https://dholk.primate.wisc.edu/..."

Environment:
    SLACK_BOT_TOKEN: Bot OAuth token (xoxb-...)

Exit codes:
    0: Always (notification failure is non-fatal to avoid failing the pipeline)
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

from loguru import logger
from py_nvd import state
from py_nvd.cli.utils import format_duration
from py_nvd.hits import (
    get_highlights_string,
    get_novel_taxa,
    get_run_report,
    get_top_movers,
)
from py_nvd.models import RunReport
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError, SlackClientError


def configure_logging(verbosity: int) -> None:
    """Configure loguru logging based on verbosity level."""
    logger.remove()

    if verbosity == 0:
        level = "WARNING"
    elif verbosity == 1:
        level = "INFO"
    else:
        level = "DEBUG"

    logger.add(
        sys.stderr,
        level=level,
        format="<level>{level: <8}</level> | {message}",
    )


def get_report_safe(state_dir: Path, sample_set_id: str) -> RunReport | None:
    """
    Safely get run report with graceful error handling.

    Returns RunReport on success, None on any error.
    All errors are logged as warnings but never raised.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    try:
        report = get_run_report(sample_set_id, state_dir)
        if report is None:
            logger.warning(
                "No hits found for sample set | sample_set_id={id}",
                id=sample_set_id,
            )
        return report

    except FileNotFoundError as e:
        logger.warning(
            "State directory or parquet files not found | path={path} | error={err}",
            path=state_dir,
            err=str(e),
        )
        return None

    except PermissionError as e:
        logger.warning(
            "Permission denied reading state | path={path} | error={err}",
            path=state_dir,
            err=str(e),
        )
        return None

    except Exception:
        logger.opt(exception=True).warning(
            "Unexpected error getting run report | sample_set_id={id}",
            id=sample_set_id,
        )
        return None


def get_highlights_safe(
    state_dir: Path,
    sample_set_id: str,
    limit: int = 5,
) -> str:
    """
    Safely get highlights string with graceful error handling.

    Returns formatted string on success, empty string on any error.
    """
    assert sample_set_id, "sample_set_id cannot be empty"
    assert limit > 0, "limit must be positive"

    try:
        result = get_highlights_string(sample_set_id, limit=limit, state_dir=state_dir)
        assert isinstance(result, str), "get_highlights_string must return str"
        return result
    except Exception:
        logger.opt(exception=True).warning(
            "Failed to get highlights | sample_set_id={id}",
            id=sample_set_id,
        )
        return ""


# Minimum number of historical runs required before showing "top movers"
MIN_BASELINE_RUNS = 3

# Minimum absolute prevalence change (percentage points) to be considered "notable"
MIN_CHANGE_THRESHOLD = 25.0

# Maximum number of movers to show in the message
MAX_MOVERS_SHOWN = 3


def get_cross_run_context(
    state_dir: Path,
    sample_set_id: str,
) -> tuple[int | None, list[tuple[str, float]]]:
    """
    Get cross-run analytics context for the Slack message.

    Returns:
        Tuple of (novel_taxa_count, notable_movers)
        - novel_taxa_count: Number of first-time taxa, or None on error
        - notable_movers: List of (taxon_name, change_pct) for taxa with
          |change| >= MIN_CHANGE_THRESHOLD, limited to MAX_MOVERS_SHOWN

    All errors are logged as warnings but never raised.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    novel_count: int | None = None
    notable_movers: list[tuple[str, float]] = []

    # Get novel taxa count
    try:
        novel_taxa = get_novel_taxa(sample_set_id, state_dir)
        novel_count = len(novel_taxa)
    except Exception:
        logger.opt(exception=True).warning(
            "Failed to get novel taxa | sample_set_id={id}",
            id=sample_set_id,
        )

    # Get top movers (only if enough baseline data)
    try:
        movers = get_top_movers(sample_set_id, limit=10, state_dir=state_dir)
        for m in movers:
            if m.baseline_run_count < MIN_BASELINE_RUNS:
                break
            if abs(m.change_pct) >= MIN_CHANGE_THRESHOLD:
                notable_movers.append((m.name, m.change_pct))
            if len(notable_movers) >= MAX_MOVERS_SHOWN:
                break
    except Exception:
        logger.opt(exception=True).warning(
            "Failed to get top movers | sample_set_id={id}",
            id=sample_set_id,
        )

    assert novel_count is None or novel_count >= 0, "novel_count must be non-negative"
    assert len(notable_movers) <= MAX_MOVERS_SHOWN, "notable_movers exceeds limit"

    return novel_count, notable_movers


def format_cross_run_section(
    novel_count: int | None,
    notable_movers: list[tuple[str, float]],
) -> str:
    """
    Format the cross-run context section for the Slack message.

    Returns empty string if there's nothing interesting to report.
    """
    lines: list[str] = []

    if novel_count is not None and novel_count > 0:
        taxa_word = "taxon" if novel_count == 1 else "taxa"
        lines.append(f"- {novel_count} novel {taxa_word} (first time seen)")

    if notable_movers:
        formatted = ", ".join(
            f"{name} ({change:+.0f}%)" for name, change in notable_movers
        )
        lines.append(f"- Notable: {formatted}")

    if not lines:
        return ""

    result = "*Cross-run context:*\n" + "\n".join(lines)
    assert result.startswith("*Cross-run context:*"), "format invariant violated"
    return result


def build_message(
    experiment_id: str,
    run_id: str,
    report: RunReport | None,
    highlights: str,
    cross_run_section: str,
    labkey_url: str,
    duration_str: str | None = None,
) -> str:
    """
    Build the Slack mrkdwn message.

    Args:
        experiment_id: Experiment identifier
        run_id: Nextflow run ID
        report: RunReport from get_run_report(), or None if unavailable
        highlights: Top findings string from get_highlights_string()
        cross_run_section: Formatted cross-run context (may be empty)
        labkey_url: URL to LabKey results
        duration_str: Formatted duration string, or None

    Returns:
        Formatted Slack message in mrkdwn format.
    """
    assert experiment_id, "experiment_id cannot be empty"
    assert run_id, "run_id cannot be empty"
    assert labkey_url, "labkey_url cannot be empty"

    # Build duration line if available
    duration_line = f"*Duration:* {duration_str}\n" if duration_str else ""

    # If we have a report, use the rich format
    if report is not None:
        samples_with_hits_pct = (
            f" ({report.samples_with_viral_hits} with hits)"
            if report.samples_with_viral_hits > 0
            else ""
        )

        results_section = f"""*Results:*
- {report.samples_analyzed:,} samples analyzed{samples_with_hits_pct}
- {report.unique_hits:,} unique contigs identified
- {report.viral_taxa_found:,} taxa detected"""

        # Add highlights if available
        highlights_line = f"\n\n*Top findings:* {highlights}" if highlights else ""

        # Add cross-run section if available
        cross_run_block = f"\n\n{cross_run_section}" if cross_run_section else ""

        message = f"""*NVD Run Complete* :white_check_mark:

*Experiment:* `{experiment_id}`
*Run:* `{run_id}`
{duration_line}
{results_section}{highlights_line}{cross_run_block}

<{labkey_url}|View on LabKey>"""

    else:
        # Fallback: minimal message when report unavailable
        message = f"""*NVD Run Complete* :white_check_mark:

*Experiment:* `{experiment_id}`
*Run:* `{run_id}`
{duration_line}
_Run statistics unavailable_

<{labkey_url}|View on LabKey>"""

    assert message.startswith("*NVD Run Complete*"), "Message format invariant violated"

    return message


# Slack API errors that are transient and should be retried
RETRYABLE_SLACK_ERRORS = frozenset(
    {
        "rate_limited",
        "ratelimited",
        "request_timeout",
        "service_unavailable",
        "fatal_error",
    },
)

# Maximum number of retry attempts for transient errors
MAX_RETRIES = 2


def send_notification(
    channel: str,
    message: str,
    token: str | None = None,
) -> bool:
    """
    Send notification to Slack with retry logic for transient errors.

    Retries up to MAX_RETRIES times for rate limiting, timeouts, and
    service unavailability. Uses exponential backoff between retries.

    Returns True on success, False on failure.
    Never raises exceptions - logs warnings instead.
    """
    assert channel, "channel cannot be empty"
    assert message, "message cannot be empty"

    token = token or os.environ.get("SLACK_BOT_TOKEN")

    if not token:
        logger.warning("SLACK_BOT_TOKEN not set, skipping notification")
        return False

    client = WebClient(token=token)

    for attempt in range(MAX_RETRIES + 1):
        try:
            response = client.chat_postMessage(
                channel=channel,
                text=message,
                mrkdwn=True,
            )
            logger.info(
                "Slack notification sent | channel={ch} | ts={ts}",
                ch=channel,
                ts=response["ts"],
            )
            return True

        except SlackApiError as e:
            error_code = e.response.get("error", "unknown")

            # Rate limited - use Retry-After header
            if error_code in ("rate_limited", "ratelimited"):
                retry_after = int(e.response.headers.get("Retry-After", 30))
                if attempt < MAX_RETRIES:
                    logger.warning(
                        "Slack rate limited, waiting {sec}s | attempt={att}/{max}",
                        sec=retry_after,
                        att=attempt + 1,
                        max=MAX_RETRIES,
                    )
                    time.sleep(retry_after)
                    continue

            # Other retryable errors - exponential backoff
            if error_code in RETRYABLE_SLACK_ERRORS and attempt < MAX_RETRIES:
                wait_time = 2**attempt
                logger.warning(
                    "Slack API error '{code}', retrying in {wait}s | attempt={att}/{max}",
                    code=error_code,
                    wait=wait_time,
                    att=attempt + 1,
                    max=MAX_RETRIES,
                )
                time.sleep(wait_time)
                continue

            # Permanent error or retries exhausted - log details and fail
            logger.warning(
                "Slack API error | channel={ch} | error={code}",
                ch=channel,
                code=error_code,
            )

            if error_code == "channel_not_found":
                logger.warning(
                    "Channel '{ch}' not found. Check channel ID.", ch=channel,
                )
            elif error_code == "not_in_channel":
                logger.warning(
                    "Bot not in channel '{ch}'. Invite the bot first.", ch=channel,
                )
            elif error_code == "invalid_auth":
                logger.warning("Invalid Slack token. Check SLACK_BOT_TOKEN.")
            elif error_code == "token_revoked":
                logger.warning("Slack token has been revoked.")
            elif error_code == "missing_scope":
                logger.warning("Slack token missing required scope (chat:write).")
            elif error_code == "msg_too_long":
                logger.warning("Message exceeds Slack's 40,000 character limit.")

            return False

        except SlackClientError as e:
            logger.warning(
                "Slack client error | channel={ch} | error={err}",
                ch=channel,
                err=str(e),
            )
            return False

        except (TimeoutError, OSError, ConnectionError) as e:
            if attempt < MAX_RETRIES:
                wait_time = 2**attempt
                logger.warning(
                    "Network error, retrying in {wait}s | attempt={att}/{max} | error={err}",
                    wait=wait_time,
                    att=attempt + 1,
                    max=MAX_RETRIES,
                    err=str(e),
                )
                time.sleep(wait_time)
                continue

            logger.warning(
                "Network error sending Slack notification | channel={ch} | error={err}",
                ch=channel,
                err=str(e),
            )
            return False

        except Exception:
            logger.opt(exception=True).warning(
                "Unexpected error sending Slack notification | channel={ch}",
                ch=channel,
            )
            return False

    # Unreachable: loop always returns on success or exhausted retries
    return False


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Send Slack notification for NVD run completion",
    )
    parser.add_argument(
        "--run-id",
        type=str,
        required=True,
        help="Nextflow run ID (workflow.runName)",
    )
    parser.add_argument(
        "--experiment-id",
        type=str,
        required=True,
        help="Experiment identifier",
    )
    parser.add_argument(
        "--channel",
        type=str,
        required=True,
        help="Slack channel ID (e.g., C0123456789)",
    )
    parser.add_argument(
        "--state-dir",
        type=Path,
        required=True,
        help="Path to state directory containing state.sqlite",
    )
    parser.add_argument(
        "--sample-set-id",
        type=str,
        required=True,
        help="Sample set ID for run-specific stats",
    )
    parser.add_argument(
        "--labkey-url",
        type=str,
        required=True,
        help="LabKey results URL",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (-v for INFO, -vv for DEBUG)",
    )

    args = parser.parse_args()

    configure_logging(args.verbose)

    logger.info(
        "Sending Slack notification | run={run} | sample_set={ss}",
        run=args.run_id,
        ss=args.sample_set_id,
    )

    # Get run report (comprehensive stats)
    report = get_report_safe(args.state_dir, args.sample_set_id)
    if report:
        logger.debug(
            "Run report | samples={s} | hits={h} | taxa={t}",
            s=report.samples_analyzed,
            h=report.unique_hits,
            t=report.viral_taxa_found,
        )
    else:
        logger.debug("Run report unavailable")

    # Get highlights string (top taxa one-liner)
    highlights = get_highlights_safe(args.state_dir, args.sample_set_id)
    if highlights:
        logger.debug("Highlights: {h}", h=highlights)

    # Get cross-run context (novel taxa, top movers)
    novel_count, notable_movers = get_cross_run_context(
        args.state_dir, args.sample_set_id,
    )
    cross_run_section = format_cross_run_section(novel_count, notable_movers)
    if cross_run_section:
        logger.debug("Cross-run context available")

    # Get run duration from state database
    duration_str: str | None = None
    try:
        run = state.get_run(args.run_id, args.state_dir)
        if run and run.duration_seconds is not None:
            duration_str = format_duration(run.duration_seconds)
            logger.debug("Run duration: {d}", d=duration_str)
    except Exception as e:
        logger.warning("Failed to get run duration | error={err}", err=str(e))

    # Build and send message
    message = build_message(
        experiment_id=args.experiment_id,
        run_id=args.run_id,
        report=report,
        highlights=highlights,
        cross_run_section=cross_run_section,
        labkey_url=args.labkey_url,
        duration_str=duration_str,
    )

    logger.debug("Message:\n{msg}", msg=message)

    success = send_notification(args.channel, message)

    # Always exit 0 - we don't want to fail the pipeline
    if success:
        print(f"Slack notification sent for run {args.run_id}")
    else:
        logger.warning("Slack notification failed, but continuing pipeline")
        print("Warning: Slack notification failed (non-fatal)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
