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
from pathlib import Path

from loguru import logger
from py_nvd.hits import get_hit_stats, get_stats_for_sample_set
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


def get_run_stats(state_dir: Path, sample_set_id: str) -> dict[str, int]:
    """
    Query run-specific stats from the state database.

    Returns dict with keys: samples_analyzed, hits_identified
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    stats = get_stats_for_sample_set(sample_set_id, state_dir)

    result = {
        "samples_analyzed": stats.unique_samples if stats else 0,
        "hits_identified": stats.total_observations if stats else 0,
    }

    assert result["samples_analyzed"] >= 0, "samples_analyzed must be non-negative"
    assert result["hits_identified"] >= 0, "hits_identified must be non-negative"

    return result


def format_byte_size(size_bytes: int) -> str:
    """Format byte count as human-readable string (KB/MB/GB)."""
    assert size_bytes >= 0, "size_bytes must be non-negative"

    if size_bytes >= 1024**3:
        return f"{size_bytes / 1024**3:.1f} GB"
    elif size_bytes >= 1024**2:
        return f"{size_bytes / 1024**2:.1f} MB"
    else:
        return f"{size_bytes / 1024:.1f} KB"


def get_cumulative_stats(state_dir: Path) -> dict[str, int | str]:
    """
    Query cumulative stats from the state database.

    Returns dict with keys: total_samples (int), total_hits (int), state_size (str)
    """
    hits_stats = get_hit_stats(state_dir)

    # Calculate state directory size
    size_bytes = 0
    state_path = Path(state_dir)
    if state_path.exists():
        for f in state_path.rglob("*"):
            if f.is_file():
                size_bytes += f.stat().st_size

    result: dict[str, int | str] = {
        "total_samples": hits_stats.unique_samples,
        "total_hits": hits_stats.total_hits,
        "state_size": format_byte_size(size_bytes),
    }

    assert isinstance(result["total_samples"], int), "total_samples must be int"
    assert isinstance(result["total_hits"], int), "total_hits must be int"
    assert result["total_samples"] >= 0, "total_samples must be non-negative"
    assert result["total_hits"] >= 0, "total_hits must be non-negative"

    return result


def build_message(
    experiment_id: str,
    run_id: str,
    samples_analyzed: int,
    hits_identified: int,
    labkey_url: str,
    cumulative: dict[str, int | str],
) -> str:
    """Build the Slack mrkdwn message."""
    assert experiment_id, "experiment_id cannot be empty"
    assert run_id, "run_id cannot be empty"
    assert labkey_url, "labkey_url cannot be empty"
    assert samples_analyzed >= 0, "samples_analyzed must be non-negative"
    assert hits_identified >= 0, "hits_identified must be non-negative"

    # Format numeric values with commas for readability
    total_samples = cumulative["total_samples"]
    total_hits = cumulative["total_hits"]
    total_samples_str = (
        f"{total_samples:,}" if isinstance(total_samples, int) else str(total_samples)
    )
    total_hits_str = (
        f"{total_hits:,}" if isinstance(total_hits, int) else str(total_hits)
    )

    message = f"""*NVD Run Complete* :white_check_mark:

*Experiment:* `{experiment_id}`
*Run:* `{run_id}`
*Samples analyzed:* {samples_analyzed}
*Hits identified (this run):* {hits_identified}

*Results:* <{labkey_url}|View Hits on LabKey>

---
_Cumulative stats:_
- Total samples in database: {total_samples_str}
- Total unique hits: {total_hits_str}
- State directory: {cumulative["state_size"]}"""

    assert message.startswith("*NVD Run Complete*"), "Message format invariant violated"

    return message


def send_notification(
    channel: str,
    message: str,
    token: str | None = None,
) -> bool:
    """
    Send notification to Slack.

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

    try:
        response = client.chat_postMessage(
            channel=channel,
            text=message,
            mrkdwn=True,
        )
        logger.info(f"Slack notification sent: ts={response['ts']}")
        return True

    except SlackApiError as e:
        error_code = e.response.get("error", "unknown")
        logger.warning(f"Slack API error: {error_code}")

        if error_code == "channel_not_found":
            logger.warning(f"Channel '{channel}' not found. Check channel ID.")
        elif error_code == "not_in_channel":
            logger.warning(f"Bot not in channel '{channel}'. Invite the bot first.")
        elif error_code == "invalid_auth":
            logger.warning("Invalid Slack token. Check SLACK_BOT_TOKEN.")

        return False

    except SlackClientError as e:
        logger.warning(f"Slack client error: {e}")
        return False

    except Exception as e:
        logger.warning(
            f"Unexpected error sending Slack notification: {e}", exc_info=True
        )
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

    logger.info(f"Sending Slack notification for run {args.run_id}")

    # Get run-specific stats from state database
    try:
        run_stats = get_run_stats(args.state_dir, args.sample_set_id)
        logger.debug(f"Run stats: {run_stats}")
    except Exception as e:
        logger.warning(f"Failed to get run stats: {e}")
        run_stats = {"samples_analyzed": 0, "hits_identified": 0}

    # Get cumulative stats from state database
    cumulative: dict[str, int | str]
    try:
        cumulative = get_cumulative_stats(args.state_dir)
        logger.debug(f"Cumulative stats: {cumulative}")
    except Exception as e:
        logger.warning(f"Failed to get cumulative stats: {e}")
        cumulative = {
            "total_samples": "N/A",
            "total_hits": "N/A",
            "state_size": "N/A",
        }

    # Build and send message
    message = build_message(
        experiment_id=args.experiment_id,
        run_id=args.run_id,
        samples_analyzed=run_stats["samples_analyzed"],
        hits_identified=run_stats["hits_identified"],
        labkey_url=args.labkey_url,
        cumulative=cumulative,
    )

    logger.debug(f"Message:\n{message}")

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
