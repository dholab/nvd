"""
Secrets management for the NVD CLI.

Wraps `nextflow secrets` to provide a friendlier interface for managing
pipeline secrets, particularly the LabKey API key used for uploading results.

The pipeline requires a secret named 'LABKEY_API_KEY' for LabKey uploads.

Note: Nextflow's secrets implementation stores secrets as plain JSON in
~/.nextflow/secrets/store.json with 600 permissions. The secrets are NOT
encrypted. Additionally, secrets are copied into .command.env files for
each task that uses them - a significant security concern for shared
filesystems.
"""

from __future__ import annotations

import getpass
import subprocess
import sys

import typer

from py_nvd.cli.prompts import (
    PROMPT_TIMEOUT_SECONDS,
    PromptTimeoutError,
    confirm_with_timeout,
    is_interactive,
)
from py_nvd.cli.utils import console, error, info, success, warning

# The secret name used by the pipeline for LabKey API key
LABKEY_SECRET_NAME = "LABKEY_API_KEY"  # noqa: S105

# The secret name used by the pipeline for Slack notifications
SLACK_SECRET_NAME = "SLACK_BOT_TOKEN"  # noqa: S105

# Number of trailing characters to show when masking a secret value
_MASK_VISIBLE_CHARS = 4

# All secret names we recognize
KNOWN_SECRETS = {
    LABKEY_SECRET_NAME: "LabKey API key for uploading results",
    SLACK_SECRET_NAME: "Slack bot token for run notifications (xoxb-...)",
}


def _run_nextflow_secrets(
    args: list[str],
    *,
    capture: bool = True,
) -> subprocess.CompletedProcess[str]:
    """
    Run a nextflow secrets command.

    Args:
        args: Arguments to pass to `nextflow secrets`
        capture: Whether to capture stdout/stderr

    Returns:
        CompletedProcess with result

    Raises:
        SystemExit: If nextflow is not found (via error())
    """
    cmd = ["nextflow", "secrets", *args]
    try:
        return subprocess.run(  # noqa: S603
            cmd,
            capture_output=capture,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        # error() calls sys.exit() and never returns
        error("Nextflow not found. Please ensure nextflow is installed and in PATH.")


def _list_secrets() -> list[str]:
    """Get list of all secret names from Nextflow."""
    result = _run_nextflow_secrets(["list"])
    if result.returncode != 0:
        return []

    output = result.stdout.strip()

    # Handle "no secrets available" message
    if not output or "no secrets available" in output.lower():
        return []

    # Parse output - nextflow secrets list outputs one name per line
    lines = output.split("\n")
    # Filter out any header/empty lines
    return [line.strip() for line in lines if line.strip() and not line.startswith("-")]


def _secret_exists(name: str) -> bool:
    """Check if a secret exists."""
    return name in _list_secrets()


def _get_secret(name: str) -> str | None:
    """Get a secret value by name."""
    result = _run_nextflow_secrets(["get", name])
    if result.returncode != 0:
        return None
    return result.stdout.strip()


def _set_secret(name: str, value: str) -> bool:
    """Set a secret value. Returns True on success."""
    result = _run_nextflow_secrets(["set", name, value])
    return result.returncode == 0


def _delete_secret(name: str) -> bool:
    """Delete a secret. Returns True on success."""
    result = _run_nextflow_secrets(["delete", name])
    return result.returncode == 0


secrets_app = typer.Typer(
    name="secrets",
    help="Manage pipeline secrets (LabKey API keys, etc.)",
    no_args_is_help=True,
)


@secrets_app.command("list")
@secrets_app.command("ls", hidden=True)  # Alias
def list_secrets(
    show_values: bool = typer.Option(
        False,
        "--show-values",
        "-v",
        help="Show secret values (use with caution)",
    ),
) -> None:
    """
    List all configured secrets.

    Shows which secrets are configured and their status. Use --show-values
    to display the actual secret values (not recommended in shared terminals).
    """
    console.print()

    secrets = _list_secrets()

    if not secrets:
        info("No secrets configured.")
        console.print()
        info("Run [cyan]nvd secrets set[/cyan] to configure the LabKey API key.")
        return

    console.print("[bold]Configured Secrets:[/bold]")
    console.print()

    for name in secrets:
        description = KNOWN_SECRETS.get(name, "(unknown)")

        if show_values:
            value = _get_secret(name)
            # Mask all but last few characters for display
            if value and len(value) > _MASK_VISIBLE_CHARS:
                masked = (
                    "*" * (len(value) - _MASK_VISIBLE_CHARS)
                    + value[-_MASK_VISIBLE_CHARS:]
                )
            else:
                masked = "****"
            console.print(f"  [cyan]{name}[/cyan] = {masked}")
            console.print(f"    [dim]{description}[/dim]")
        else:
            console.print(f"  [cyan]{name}[/cyan] - {description}")

    console.print()

    # Check if the required secret is configured
    has_labkey_secret = LABKEY_SECRET_NAME in secrets

    if not has_labkey_secret:
        warning(f"'{LABKEY_SECRET_NAME}' is not configured. LabKey uploads will fail.")
        console.print()
        info("To configure, run:")
        console.print(f"  [cyan]nvd secrets set {LABKEY_SECRET_NAME}[/cyan]")
    else:
        success(
            f"'{LABKEY_SECRET_NAME}' is configured - pipeline ready for LabKey uploads",
        )


@secrets_app.command("set")
def set_secret(
    name: str = typer.Argument(
        LABKEY_SECRET_NAME,
        help=f"Secret name (default: {LABKEY_SECRET_NAME})",
    ),
    value: str | None = typer.Argument(
        None,
        help="Secret value (will prompt securely if not provided)",
    ),
) -> None:
    """
    Set a secret value.

    If the value is not provided as an argument, you will be prompted to
    enter it securely (input will be hidden). This is the recommended
    approach to avoid the secret appearing in shell history.

    The pipeline requires a secret named 'LABKEY_API_KEY' for LabKey uploads.

    Examples:

        # Set LabKey API key (will prompt for value securely)
        nvd secrets set

        # Set with value directly (appears in shell history - not recommended)
        nvd secrets set LABKEY_API_KEY "your-api-key"
    """
    console.print()

    # Prompt for value if not provided
    if value is None:
        if not is_interactive():
            error(
                "Cannot prompt for secret value in non-interactive mode. Provide value as argument.",
            )

        console.print(f"[bold]Setting secret:[/bold] {name}")
        console.print()

        try:
            # Use getpass for secure input (no echo)
            value = getpass.getpass(f"Enter value for '{name}': ")
            if not value:
                error("Secret value cannot be empty.")

            # Confirm by entering again
            value_confirm = getpass.getpass(f"Confirm value for '{name}': ")
            if value != value_confirm:
                error("Values do not match.")

        except KeyboardInterrupt:
            console.print("\n[yellow]Cancelled[/yellow]")
            raise typer.Abort from None

    # Check if overwriting
    if _secret_exists(name):
        warning(f"Secret '{name}' already exists and will be overwritten.")
        if is_interactive():
            try:
                if not confirm_with_timeout(
                    console,
                    "Continue?",
                    default=True,
                    timeout_seconds=PROMPT_TIMEOUT_SECONDS,
                ):
                    info("Aborted")
                    raise typer.Abort from None
            except PromptTimeoutError:
                error("Confirmation timed out.")
            except KeyboardInterrupt:
                console.print("\n[yellow]Cancelled[/yellow]")
                raise typer.Abort from None

    # Set the secret
    if _set_secret(name, value):
        success(f"Secret '{name}' has been set.")
    else:
        error(f"Failed to set secret '{name}'.")


@secrets_app.command("get")
def get_secret(
    name: str = typer.Argument(..., help="Secret name to retrieve"),
) -> None:
    """
    Get a secret value (DEPRECATED).

    This command is deprecated and will be removed in a future release.
    Use 'nvd secrets list --show-values' for a masked view instead.

    Displays the secret value in plain text to stdout. Requires interactive
    confirmation with no way to bypass - this is intentional to discourage
    use in scripts.
    """
    console.print()

    # Deprecation warning
    warning(
        "[bold]This command is deprecated[/bold] and will be removed in a future release.",
    )
    info("Use [cyan]nvd secrets list --show-values[/cyan] for a masked view instead.")
    console.print()

    # Must be interactive - no bypass allowed
    if not is_interactive():
        error(
            "This command requires interactive confirmation and cannot be used "
            "in scripts or non-interactive environments.",
        )

    if not _secret_exists(name):
        error(f"Secret '{name}' not found.")

    # Require explicit confirmation
    warning("This will print the secret value in plain text to stdout.")
    warning("The value may be visible in terminal history, screen recordings, or logs.")
    console.print()

    try:
        if not confirm_with_timeout(
            console,
            "Are you sure you want to display the secret?",
            default=False,
            timeout_seconds=PROMPT_TIMEOUT_SECONDS,
        ):
            info("Aborted.")
            raise typer.Abort from None
    except PromptTimeoutError:
        console.print()
        error(
            f"Confirmation timed out after {PROMPT_TIMEOUT_SECONDS // 60} minutes. "
            "No secret was displayed.",
        )
    except KeyboardInterrupt:
        console.print("\n[yellow]Cancelled[/yellow]")
        raise typer.Abort from None

    console.print()

    value = _get_secret(name)
    if value is None:
        error(f"Failed to retrieve secret '{name}'.")

    console.print(f"[bold]{name}[/bold] = {value}")


@secrets_app.command("delete")
@secrets_app.command("rm", hidden=True)  # Alias
def delete_secret(
    name: str = typer.Argument(..., help="Secret name to delete"),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Delete without confirmation",
    ),
) -> None:
    """
    Delete a secret.

    Permanently removes a secret from the Nextflow secrets store.
    """
    console.print()

    if not _secret_exists(name):
        error(f"Secret '{name}' not found.")

    if not force:
        warning(f"This will permanently delete secret '{name}'.")
        if is_interactive():
            try:
                if not confirm_with_timeout(
                    console,
                    "Are you sure?",
                    default=False,
                    timeout_seconds=PROMPT_TIMEOUT_SECONDS,
                ):
                    info("Aborted")
                    raise typer.Abort from None
            except PromptTimeoutError:
                error("Confirmation timed out.")
            except KeyboardInterrupt:
                console.print("\n[yellow]Cancelled[/yellow]")
                raise typer.Abort from None
        else:
            error("Use --force to delete in non-interactive mode.")

    if _delete_secret(name):
        success(f"Secret '{name}' has been deleted.")
    else:
        error(f"Failed to delete secret '{name}'.")


@secrets_app.command("check")
def check_secrets() -> None:
    """
    Check if required secrets are configured.

    Verifies that the LABKEY_API_KEY secret is configured for pipeline
    uploads to LabKey.
    """
    console.print()
    console.print("[bold]Checking Pipeline Secrets[/bold]")
    console.print()

    secrets = _list_secrets()

    has_labkey = LABKEY_SECRET_NAME in secrets

    if has_labkey:
        success(f"'{LABKEY_SECRET_NAME}' is configured")
        console.print()
        success("Secret configuration is valid - pipeline ready for LabKey uploads")
    else:
        console.print(
            f"[red]✗[/red] '{LABKEY_SECRET_NAME}' is [bold]not configured[/bold]",
        )
        console.print()
        info("LabKey uploads will fail without this secret.")
        console.print()
        info("To configure, run:")
        console.print(f"  [cyan]nvd secrets set {LABKEY_SECRET_NAME}[/cyan]")
        sys.exit(1)
