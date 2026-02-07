"""
Tests for notify_slack.py Slack notification script.

The graceful-failure design means error paths are the feature—
they need thorough coverage to ensure silent failures stay silent.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
from notify_slack import (
    build_message,
    format_byte_size,
    get_cumulative_stats,
    get_run_stats,
    send_notification,
)
from slack_sdk.errors import SlackApiError, SlackClientError


class TestFormatByteSize:
    """Tests for byte size formatting."""

    def test_formats_zero_as_kb(self):
        """Zero bytes displays as 0.0 KB."""
        assert format_byte_size(0) == "0.0 KB"

    def test_formats_small_size_as_kb(self):
        """Sizes under 1MB display as KB."""
        assert format_byte_size(500 * 1024) == "500.0 KB"

    def test_formats_size_as_mb(self):
        """Sizes 1MB-1GB display as MB."""
        assert format_byte_size(50 * 1024**2) == "50.0 MB"

    def test_formats_size_as_gb(self):
        """Sizes over 1GB display as GB."""
        assert format_byte_size(2 * 1024**3) == "2.0 GB"

    def test_boundary_exactly_1mb(self):
        """Exactly 1MB (1024**2 bytes) shows as MB not KB."""
        assert format_byte_size(1024**2) == "1.0 MB"

    def test_boundary_exactly_1gb(self):
        """Exactly 1GB (1024**3 bytes) shows as GB not MB."""
        assert format_byte_size(1024**3) == "1.0 GB"

    def test_boundary_just_under_1mb(self):
        """Just under 1MB shows as KB."""
        assert format_byte_size(1024**2 - 1) == "1024.0 KB"

    def test_fractional_gb(self):
        """Fractional GB values display with one decimal."""
        assert format_byte_size(int(2.5 * 1024**3)) == "2.5 GB"

    def test_assertion_on_negative(self):
        """Negative byte count raises assertion."""
        with pytest.raises(AssertionError, match="non-negative"):
            format_byte_size(-1)


class TestBuildMessage:
    """Tests for Slack message formatting."""

    def test_formats_all_fields_correctly(self):
        """All parameters appear in the output message."""
        message = build_message(
            experiment_id="exp_001",
            run_id="run_2026-01-10",
            samples_analyzed=42,
            hits_identified=7,
            labkey_url="https://example.com/results",
            cumulative={
                "total_samples": 1000,
                "total_hits": 500,
                "state_size": "2.3 GB",
            },
        )

        assert "exp_001" in message
        assert "run_2026-01-10" in message
        assert "42" in message
        assert "7" in message
        assert "https://example.com/results" in message
        assert "1,000" in message  # formatted with comma
        assert "500" in message
        assert "2.3 GB" in message

    def test_handles_zero_counts(self):
        """Zero samples/hits display correctly."""
        message = build_message(
            experiment_id="exp_001",
            run_id="run_001",
            samples_analyzed=0,
            hits_identified=0,
            labkey_url="https://example.com",
            cumulative={
                "total_samples": 0,
                "total_hits": 0,
                "state_size": "0.0 KB",
            },
        )

        assert "*Samples analyzed:* 0" in message
        assert "*Hits identified (this run):* 0" in message
        assert "Total samples in database: 0" in message

    def test_labkey_url_in_link_format(self):
        """URL is wrapped in Slack link syntax."""
        message = build_message(
            experiment_id="exp",
            run_id="run",
            samples_analyzed=1,
            hits_identified=1,
            labkey_url="https://dholk.example.com/results",
            cumulative={"total_samples": 1, "total_hits": 1, "state_size": "1 KB"},
        )

        assert "<https://dholk.example.com/results|View Hits on LabKey>" in message

    def test_starts_with_header(self):
        """Message starts with the expected header."""
        message = build_message(
            experiment_id="exp",
            run_id="run",
            samples_analyzed=1,
            hits_identified=1,
            labkey_url="https://example.com",
            cumulative={"total_samples": 1, "total_hits": 1, "state_size": "1 KB"},
        )

        assert message.startswith("*NVD Run Complete*")

    def test_handles_string_cumulative_values(self):
        """Handles N/A strings in cumulative stats (fallback case)."""
        message = build_message(
            experiment_id="exp",
            run_id="run",
            samples_analyzed=1,
            hits_identified=1,
            labkey_url="https://example.com",
            cumulative={
                "total_samples": "N/A",
                "total_hits": "N/A",
                "state_size": "N/A",
            },
        )

        assert "Total samples in database: N/A" in message
        assert "Total unique hits: N/A" in message
        assert "State directory: N/A" in message

    def test_assertion_on_empty_experiment_id(self):
        """Empty experiment_id raises assertion."""
        with pytest.raises(AssertionError, match="experiment_id cannot be empty"):
            build_message(
                experiment_id="",
                run_id="run",
                samples_analyzed=1,
                hits_identified=1,
                labkey_url="https://example.com",
                cumulative={"total_samples": 1, "total_hits": 1, "state_size": "1 KB"},
            )

    def test_assertion_on_empty_run_id(self):
        """Empty run_id raises assertion."""
        with pytest.raises(AssertionError, match="run_id cannot be empty"):
            build_message(
                experiment_id="exp",
                run_id="",
                samples_analyzed=1,
                hits_identified=1,
                labkey_url="https://example.com",
                cumulative={"total_samples": 1, "total_hits": 1, "state_size": "1 KB"},
            )

    def test_assertion_on_negative_samples(self):
        """Negative samples_analyzed raises assertion."""
        with pytest.raises(AssertionError, match="non-negative"):
            build_message(
                experiment_id="exp",
                run_id="run",
                samples_analyzed=-1,
                hits_identified=1,
                labkey_url="https://example.com",
                cumulative={"total_samples": 1, "total_hits": 1, "state_size": "1 KB"},
            )


class TestGetRunStats:
    """Tests for run-specific stats extraction."""

    def test_returns_zeros_when_stats_is_none(self, tmp_path):
        """When get_stats_for_sample_set returns None, returns zeros."""
        with patch(
            "notify_slack.get_stats_for_sample_set",
            return_value=None,
        ) as mock_get:
            result = get_run_stats(tmp_path, "nonexistent_sample_set")

        assert result == {"samples_analyzed": 0, "hits_identified": 0}
        mock_get.assert_called_once_with("nonexistent_sample_set", tmp_path)

    def test_extracts_correct_fields_from_stats(self, tmp_path):
        """Maps unique_samples → samples_analyzed, total_observations → hits_identified."""
        mock_stats = MagicMock()
        mock_stats.unique_samples = 42
        mock_stats.total_observations = 7

        with patch("notify_slack.get_stats_for_sample_set", return_value=mock_stats):
            result = get_run_stats(tmp_path, "sample_set_123")

        assert result == {"samples_analyzed": 42, "hits_identified": 7}

    def test_assertion_on_empty_sample_set_id(self, tmp_path):
        """Empty sample_set_id raises assertion."""
        with pytest.raises(AssertionError, match="sample_set_id cannot be empty"):
            get_run_stats(tmp_path, "")


class TestGetCumulativeStats:
    """Tests for cumulative stats and size formatting."""

    def test_returns_stats_from_hit_stats(self, tmp_path):
        """Returns unique_samples and total_hits from get_hit_stats."""
        mock_stats = MagicMock()
        mock_stats.unique_samples = 1000
        mock_stats.total_hits = 500

        with patch("notify_slack.get_hit_stats", return_value=mock_stats):
            result = get_cumulative_stats(tmp_path)

        assert result["total_samples"] == 1000
        assert result["total_hits"] == 500

    def test_calculates_directory_size(self, tmp_path):
        """Calculates state directory size from files."""
        # Create some test files
        (tmp_path / "file1.txt").write_text("x" * 1000)
        (tmp_path / "file2.txt").write_text("y" * 2000)
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        (subdir / "file3.txt").write_text("z" * 500)

        mock_stats = MagicMock()
        mock_stats.unique_samples = 0
        mock_stats.total_hits = 0

        with patch("notify_slack.get_hit_stats", return_value=mock_stats):
            result = get_cumulative_stats(tmp_path)

        # Total should be 3500 bytes = 3.4 KB
        assert result["state_size"] == "3.4 KB"

    def test_handles_empty_directory(self, tmp_path):
        """Empty state directory shows 0.0 KB."""
        mock_stats = MagicMock()
        mock_stats.unique_samples = 0
        mock_stats.total_hits = 0

        with patch("notify_slack.get_hit_stats", return_value=mock_stats):
            result = get_cumulative_stats(tmp_path)

        assert result["state_size"] == "0.0 KB"

    def test_handles_nonexistent_directory(self, tmp_path):
        """Nonexistent state directory shows 0.0 KB."""
        nonexistent = tmp_path / "does_not_exist"

        mock_stats = MagicMock()
        mock_stats.unique_samples = 0
        mock_stats.total_hits = 0

        with patch("notify_slack.get_hit_stats", return_value=mock_stats):
            result = get_cumulative_stats(nonexistent)

        assert result["state_size"] == "0.0 KB"


class TestSendNotification:
    """Tests for Slack API interaction and error handling."""

    def test_returns_false_when_no_token(self, monkeypatch):
        """Returns False and logs warning when SLACK_BOT_TOKEN not set."""
        monkeypatch.delenv("SLACK_BOT_TOKEN", raising=False)

        result = send_notification("C123", "test message", token=None)

        assert result is False

    def test_returns_true_on_success(self, monkeypatch):
        """Returns True when message posts successfully."""
        mock_client = MagicMock()
        mock_client.chat_postMessage.return_value = {"ts": "1234567890.123456"}

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification("C123", "test message", token="xoxb-test")

        assert result is True
        mock_client.chat_postMessage.assert_called_once_with(
            channel="C123",
            text="test message",
            mrkdwn=True,
        )

    def test_handles_channel_not_found(self):
        """Returns False with helpful message for channel_not_found."""
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_response.get.return_value = "channel_not_found"
        mock_client.chat_postMessage.side_effect = SlackApiError(
            message="channel_not_found",
            response=mock_response,
        )

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification("C123", "test message", token="xoxb-test")

        assert result is False

    def test_handles_not_in_channel(self):
        """Returns False with helpful message for not_in_channel."""
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_response.get.return_value = "not_in_channel"
        mock_client.chat_postMessage.side_effect = SlackApiError(
            message="not_in_channel",
            response=mock_response,
        )

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification("C123", "test message", token="xoxb-test")

        assert result is False

    def test_handles_invalid_auth(self):
        """Returns False with helpful message for invalid_auth."""
        mock_client = MagicMock()
        mock_response = MagicMock()
        mock_response.get.return_value = "invalid_auth"
        mock_client.chat_postMessage.side_effect = SlackApiError(
            message="invalid_auth",
            response=mock_response,
        )

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification("C123", "test message", token="xoxb-test")

        assert result is False

    def test_handles_slack_client_error(self):
        """Returns False for SlackClientError (network issues)."""
        mock_client = MagicMock()
        mock_client.chat_postMessage.side_effect = SlackClientError("Network error")

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification("C123", "test message", token="xoxb-test")

        assert result is False

    def test_handles_unexpected_exception(self):
        """Returns False and doesn't raise for unexpected errors."""
        mock_client = MagicMock()
        mock_client.chat_postMessage.side_effect = RuntimeError("Unexpected!")

        with patch("notify_slack.WebClient", return_value=mock_client):
            result = send_notification("C123", "test message", token="xoxb-test")

        assert result is False

    def test_assertion_on_empty_channel(self):
        """Empty channel raises assertion."""
        with pytest.raises(AssertionError, match="channel cannot be empty"):
            send_notification("", "test message", token="xoxb-test")

    def test_assertion_on_empty_message(self):
        """Empty message raises assertion."""
        with pytest.raises(AssertionError, match="message cannot be empty"):
            send_notification("C123", "", token="xoxb-test")

    def test_uses_env_token_when_not_provided(self, monkeypatch):
        """Uses SLACK_BOT_TOKEN from environment when token param is None."""
        monkeypatch.setenv("SLACK_BOT_TOKEN", "xoxb-from-env")

        mock_client = MagicMock()
        mock_client.chat_postMessage.return_value = {"ts": "123"}

        with patch(
            "notify_slack.WebClient",
            return_value=mock_client,
        ) as mock_webclient:
            result = send_notification("C123", "test", token=None)

        assert result is True
        mock_webclient.assert_called_once_with(token="xoxb-from-env")
