"""Tests for py_nvd.cli.utils module."""

import pytest

from py_nvd.cli.utils import format_duration


class TestFormatDuration:
    """Tests for format_duration function."""

    def test_zero_seconds(self) -> None:
        """Zero seconds formats correctly."""
        assert format_duration(0) == "0s"

    def test_seconds_only(self) -> None:
        """Values under 60 seconds show seconds only."""
        assert format_duration(1) == "1s"
        assert format_duration(30) == "30s"
        assert format_duration(59) == "59s"

    def test_minutes_and_seconds(self) -> None:
        """Values between 1-60 minutes show minutes and seconds."""
        assert format_duration(60) == "1m 0s"
        assert format_duration(61) == "1m 1s"
        assert format_duration(125) == "2m 5s"
        assert format_duration(3599) == "59m 59s"

    def test_hours_and_minutes(self) -> None:
        """Values >= 1 hour show hours and minutes (no seconds)."""
        assert format_duration(3600) == "1h 0m"
        assert format_duration(3660) == "1h 1m"
        assert format_duration(3725) == "1h 2m"
        assert format_duration(7200) == "2h 0m"

    def test_boundary_59_to_60_seconds(self) -> None:
        """Boundary between seconds-only and minutes+seconds."""
        assert format_duration(59) == "59s"
        assert format_duration(60) == "1m 0s"

    def test_boundary_3599_to_3600_seconds(self) -> None:
        """Boundary between minutes+seconds and hours+minutes."""
        assert format_duration(3599) == "59m 59s"
        assert format_duration(3600) == "1h 0m"

    def test_fractional_seconds_truncated(self) -> None:
        """Fractional seconds are truncated, not rounded."""
        assert format_duration(0.9) == "0s"
        assert format_duration(1.9) == "1s"
        assert format_duration(59.9) == "59s"

    def test_fractional_minutes_truncated(self) -> None:
        """Fractional values in minute range are truncated."""
        assert format_duration(60.9) == "1m 0s"
        assert format_duration(119.9) == "1m 59s"

    def test_fractional_hours_truncated(self) -> None:
        """Fractional values in hour range are truncated."""
        assert format_duration(3600.9) == "1h 0m"
        assert format_duration(3659.9) == "1h 0m"  # 59.9s truncates to 0m

    def test_large_values_24_hours(self) -> None:
        """24 hours (86400 seconds) formats correctly."""
        assert format_duration(86400) == "24h 0m"

    def test_large_values_100_hours(self) -> None:
        """100 hours formats correctly."""
        assert format_duration(360000) == "100h 0m"

    def test_large_values_with_minutes(self) -> None:
        """Large values with non-zero minutes."""
        assert format_duration(90000) == "25h 0m"
        assert format_duration(90060) == "25h 1m"

    def test_negative_raises_assertion_error(self) -> None:
        """Negative duration raises AssertionError."""
        with pytest.raises(AssertionError, match="must be non-negative"):
            format_duration(-1)

    def test_negative_fractional_raises_assertion_error(self) -> None:
        """Negative fractional duration raises AssertionError."""
        with pytest.raises(AssertionError, match="must be non-negative"):
            format_duration(-0.5)

    def test_nan_raises_assertion_error(self) -> None:
        """NaN raises AssertionError."""
        with pytest.raises(AssertionError, match="cannot be NaN"):
            format_duration(float("nan"))

    def test_positive_infinity_raises_assertion_error(self) -> None:
        """Positive infinity raises AssertionError."""
        with pytest.raises(AssertionError, match="cannot be infinite"):
            format_duration(float("inf"))

    def test_negative_infinity_raises_assertion_error(self) -> None:
        """Negative infinity raises AssertionError."""
        with pytest.raises(AssertionError, match="cannot be infinite"):
            format_duration(float("-inf"))

    def test_integer_input(self) -> None:
        """Integer input works correctly."""
        assert format_duration(60) == "1m 0s"
        assert format_duration(3600) == "1h 0m"

    def test_float_input(self) -> None:
        """Float input works correctly."""
        assert format_duration(60.0) == "1m 0s"
        assert format_duration(3600.0) == "1h 0m"
