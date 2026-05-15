"""
Tests for notify_slack.py Slack notification script.

The graceful-failure design means error paths are the feature—
they need thorough coverage to ensure silent failures stay silent.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Final
from unittest.mock import MagicMock, patch

import pytest
from notify_slack import (
    build_message,
    format_cross_run_section,
    get_highlights_safe,
    get_report_safe,
    send_notification,
)
from py_nvd.models import RunReport
from slack_sdk.errors import SlackApiError, SlackClientError

if TYPE_CHECKING:
    from pathlib import Path

FAKE_SLACK_BOT_TOKEN: Final = "xoxb-test"  # noqa: S105 - fake Slack-shaped token for exercising WebClient wiring.
FAKE_ENV_SLACK_BOT_TOKEN: Final = "xoxb-from-env"  # noqa: S105 - fake Slack-shaped token for env fallback coverage.


def make_report(**overrides: object) -> RunReport:
    """Build a minimal RunReport for Slack message tests."""
    values = {
        "sample_set_id": "sample_set_123",
        "run_date": "2026-01-10",
        "samples_analyzed": 42,
        "samples_with_viral_hits": 7,
        "samples_negative": 35,
        "unique_hits": 11,
        "viral_taxa_found": 3,
        "median_contig_length": 750,
        "contigs_over_500bp": 8,
        "top_findings": (),
        "sample_summaries": (),
    }
    values.update(overrides)
    return RunReport(**values)


class TestBuildMessage:
    """Tests for Slack message formatting."""

    def test_formats_all_fields_correctly(self) -> None:
        """Report, highlights, and cross-run context appear in the message."""
        message = build_message(
            experiment_id="exp_001",
            run_id="run_2026-01-10",
            report=make_report(),
            highlights="Adenovirus in 3 samples",
            cross_run_section="*Cross-run context:*\n- 2 novel taxa (first time seen)",
            labkey_url="https://example.com/results",
            duration_str="12m 34s",
        )

        assert "exp_001" in message
        assert "run_2026-01-10" in message
        assert "42" in message
        assert "11" in message
        assert "3" in message
        assert "Adenovirus in 3 samples" in message
        assert "2 novel taxa" in message
        assert "12m 34s" in message
        assert "https://example.com/results" in message

    def test_handles_zero_report_counts(self) -> None:
        """Zero values in a report display correctly."""
        message = build_message(
            experiment_id="exp_001",
            run_id="run_001",
            report=make_report(
                samples_analyzed=0,
                samples_with_viral_hits=0,
                samples_negative=0,
                unique_hits=0,
                viral_taxa_found=0,
            ),
            highlights="",
            cross_run_section="",
            labkey_url="https://example.com",
        )

        assert "0 samples analyzed" in message
        assert "0 unique contigs identified" in message
        assert "0 taxa detected" in message

    def test_labkey_url_in_link_format(self) -> None:
        """URL is wrapped in Slack link syntax."""
        message = build_message(
            experiment_id="exp",
            run_id="run",
            report=make_report(samples_analyzed=1, unique_hits=1, viral_taxa_found=1),
            highlights="",
            cross_run_section="",
            labkey_url="https://dholk.example.com/results",
        )

        assert "<https://dholk.example.com/results|View on LabKey>" in message

    def test_starts_with_header(self) -> None:
        """Message starts with the expected header."""
        message = build_message(
            experiment_id="exp",
            run_id="run",
            report=make_report(samples_analyzed=1, unique_hits=1, viral_taxa_found=1),
            highlights="",
            cross_run_section="",
            labkey_url="https://example.com",
        )

        assert message.startswith("*NVD Run Complete*")

    def test_handles_missing_report(self) -> None:
        """Fallback message is used when run statistics are unavailable."""
        message = build_message(
            experiment_id="exp",
            run_id="run",
            report=None,
            highlights="",
            cross_run_section="",
            labkey_url="https://example.com",
        )

        assert "Run statistics unavailable" in message

    def test_assertion_on_empty_experiment_id(self) -> None:
        """Empty experiment_id raises assertion."""
        with pytest.raises(AssertionError, match="experiment_id cannot be empty"):
            build_message(
                experiment_id="",
                run_id="run",
                report=make_report(),
                highlights="",
                cross_run_section="",
                labkey_url="https://example.com",
            )

    def test_assertion_on_empty_run_id(self) -> None:
        """Empty run_id raises assertion."""
        with pytest.raises(AssertionError, match="run_id cannot be empty"):
            build_message(
                experiment_id="exp",
                run_id="",
                report=make_report(),
                highlights="",
                cross_run_section="",
                labkey_url="https://example.com",
            )

    def test_assertion_on_empty_labkey_url(self) -> None:
        """Empty labkey_url raises assertion."""
        with pytest.raises(AssertionError, match="labkey_url cannot be empty"):
            build_message(
                experiment_id="exp",
                run_id="run",
                report=make_report(),
                highlights="",
                cross_run_section="",
                labkey_url="",
            )


class TestReportHelpers:
    """Tests for report helper fallbacks."""

    def test_get_report_safe_returns_report(self, tmp_path: Path) -> None:
        """get_report_safe returns the report from get_run_report."""
        report = make_report()
        with patch("notify_slack.get_run_report", return_value=report):
            assert get_report_safe(tmp_path, "sample_set_123") == report

    def test_get_report_safe_returns_none_on_exception(self, tmp_path: Path) -> None:
        """get_report_safe swallows read errors and returns None."""
        with patch("notify_slack.get_run_report", side_effect=RuntimeError("boom")):
            assert get_report_safe(tmp_path, "sample_set_123") is None

    def test_get_highlights_safe_returns_empty_string_on_exception(
        self,
        tmp_path: Path,
    ) -> None:
        """get_highlights_safe swallows errors and returns an empty string."""
        with patch(
            "notify_slack.get_highlights_string",
            side_effect=RuntimeError("boom"),
        ):
            assert get_highlights_safe(tmp_path, "sample_set_123") == ""


class TestFormatCrossRunSection:
    """Tests for cross-run section formatting."""

    def test_empty_context_returns_empty_string(self) -> None:
        """No novel taxa or movers produce no section."""
        assert format_cross_run_section(None, []) == ""

    def test_formats_novel_taxa(self) -> None:
        """Novel taxa counts are included."""
        assert "2 novel taxa" in format_cross_run_section(2, [])

    def test_formats_notable_movers(self) -> None:
        """Notable movers are included with signed percentages."""
        section = format_cross_run_section(None, [("Adenoviridae", 42.0)])
        assert "Adenoviridae (+42%)" in section


class TestSendNotification:
    """Tests for Slack API interaction and error handling."""

    def test_returns_false_when_no_token(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Returns False and logs warning when SLACK_BOT_TOKEN not set."""
        monkeypatch.delenv("SLACK_BOT_TOKEN", raising=False)

        result = send_notification("C123", "test message", token=None)

        assert result is False

    def test_returns_true_on_success(self) -> None:
        """Returns True when message posts successfully."""
        mock_client = MagicMock()
        mock_client.chat_postMessage.return_value = {"ts": "1234567890.123456"}

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification(
                "C123", "test message", token=FAKE_SLACK_BOT_TOKEN
            )

        assert result is True
        mock_client.chat_postMessage.assert_called_once_with(
            channel="C123",
            text="test message",
            mrkdwn=True,
        )

    def test_handles_channel_not_found(self) -> None:
        """Returns False with helpful message for channel_not_found."""
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_response.get.return_value = "channel_not_found"
        mock_client.chat_postMessage.side_effect = SlackApiError(
            message="channel_not_found",
            response=mock_response,
        )

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification(
                "C123", "test message", token=FAKE_SLACK_BOT_TOKEN
            )

        assert result is False

    def test_handles_not_in_channel(self) -> None:
        """Returns False with helpful message for not_in_channel."""
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_response.get.return_value = "not_in_channel"
        mock_client.chat_postMessage.side_effect = SlackApiError(
            message="not_in_channel",
            response=mock_response,
        )

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification(
                "C123", "test message", token=FAKE_SLACK_BOT_TOKEN
            )

        assert result is False

    def test_handles_invalid_auth(self) -> None:
        """Returns False with helpful message for invalid_auth."""
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_response.get.return_value = "invalid_auth"
        mock_client.chat_postMessage.side_effect = SlackApiError(
            message="invalid_auth",
            response=mock_response,
        )

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification(
                "C123", "test message", token=FAKE_SLACK_BOT_TOKEN
            )

        assert result is False

    def test_handles_slack_client_error(self) -> None:
        """Returns False for SlackClientError (network issues)."""
        mock_client = MagicMock()
        mock_client.chat_postMessage.side_effect = SlackClientError("Network error")

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification(
                "C123", "test message", token=FAKE_SLACK_BOT_TOKEN
            )

        assert result is False

    def test_handles_unexpected_exception(self) -> None:
        """Returns False and doesn't raise for unexpected errors."""
        mock_client = MagicMock()
        mock_client.chat_postMessage.side_effect = RuntimeError("Unexpected!")

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification(
                "C123", "test message", token=FAKE_SLACK_BOT_TOKEN
            )

        assert result is False

    def test_assertion_on_empty_channel(self) -> None:
        """Empty channel raises assertion."""
        with pytest.raises(AssertionError, match="channel cannot be empty"):
            send_notification("", "test message", token=FAKE_SLACK_BOT_TOKEN)

    def test_assertion_on_empty_message(self) -> None:
        """Empty message raises assertion."""
        with pytest.raises(AssertionError, match="message cannot be empty"):
            send_notification("C123", "", token=FAKE_SLACK_BOT_TOKEN)

    def test_uses_env_token_when_not_provided(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Uses SLACK_BOT_TOKEN from environment when token param is None."""
        monkeypatch.setenv("SLACK_BOT_TOKEN", FAKE_ENV_SLACK_BOT_TOKEN)

        mock_client = MagicMock()
        mock_client.chat_postMessage.return_value = {"ts": "123"}

        with patch(
            "notify_slack.WebClient",
            return_value=mock_client,
        ) as mock_webclient:
            result = send_notification("C123", "test", token=None)

        assert result is True
        mock_webclient.assert_called_once_with(token=FAKE_ENV_SLACK_BOT_TOKEN)
