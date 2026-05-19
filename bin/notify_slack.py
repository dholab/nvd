#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
#     "slack_sdk",
# ]
# ///

"""Send a state-free Slack notification for NVD run completion."""

from __future__ import annotations

import argparse
import os
import sys

from loguru import logger
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError, SlackClientError


def configure_logging(verbosity: int) -> None:
    logger.remove()
    level = "WARNING" if verbosity == 0 else "INFO" if verbosity == 1 else "DEBUG"
    logger.add(sys.stderr, level=level, format="<level>{level: <8}</level> | {message}")


def build_message(
    *,
    run_id: str,
    experiment_id: str,
    sample_set_id: str,
    labkey_url: str,
) -> str:
    """Build the completion message body."""
    return "\n".join(
        [
            ":white_check_mark: *NVD run completed*",
            f"*Run:* `{run_id}`",
            f"*Experiment:* `{experiment_id}`",
            f"*Sample set:* `{sample_set_id}`",
            f"*LabKey results:* {labkey_url}",
        ],
    )


def send_message(channel: str, message: str) -> bool:
    """Send a Slack message, returning False rather than raising on failures."""
    token = os.environ.get("SLACK_BOT_TOKEN")
    if not token:
        logger.warning("SLACK_BOT_TOKEN not set; skipping Slack notification")
        return False

    try:
        client = WebClient(token=token)
        response = client.chat_postMessage(channel=channel, text=message, mrkdwn=True)
    except (SlackApiError, SlackClientError, OSError) as exc:
        logger.warning(f"Slack notification failed: {exc}")
        return False

    ok = bool(response.get("ok", False))
    if not ok:
        logger.warning(f"Slack API returned non-ok response: {response}")
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(description="Send NVD completion notification")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--channel", required=True)
    parser.add_argument("--sample-set-id", required=True)
    parser.add_argument("--labkey-url", required=True)
    parser.add_argument("-v", "--verbose", action="count", default=0)
    args = parser.parse_args()

    configure_logging(args.verbose)
    message = build_message(
        run_id=args.run_id,
        experiment_id=args.experiment_id,
        sample_set_id=args.sample_set_id,
        labkey_url=args.labkey_url,
    )
    send_message(args.channel, message)


if __name__ == "__main__":
    main()
