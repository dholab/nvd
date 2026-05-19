"""Tests for the stateless Slack notification helper."""

from __future__ import annotations

from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

from notify_slack import build_message, send_message

if TYPE_CHECKING:
    import pytest


def test_build_message_includes_run_context() -> None:
    message = build_message(
        run_id="run_001",
        experiment_id="exp_001",
        sample_set_id="abc123",
        labkey_url="https://example.com/results",
    )

    assert "run_001" in message
    assert "exp_001" in message
    assert "abc123" in message
    assert "https://example.com/results" in message


def test_send_message_skips_without_token(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SLACK_BOT_TOKEN", raising=False)
    assert not send_message("C123", "hello")


def test_send_message_uses_slack_client(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SLACK_BOT_TOKEN", "xoxb-test")
    fake_client = MagicMock()
    fake_client.chat_postMessage.return_value = {"ok": True}

    with patch("notify_slack.WebClient", return_value=fake_client):
        assert send_message("C123", "hello")

    fake_client.chat_postMessage.assert_called_once_with(
        channel="C123",
        text="hello",
        mrkdwn=True,
    )
