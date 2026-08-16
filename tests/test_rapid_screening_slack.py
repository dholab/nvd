"""Nextflow contract tests for rapid-screening Slack notifications."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NOTIFICATIONS_MODULE = ROOT / "modules" / "notifications.nf"


def write_csv(path: Path, rows: list[dict[str, str]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_notifier_uses_hardcoded_rg3_policy_without_contacting_slack(
    tmp_path: Path,
) -> None:
    nextflow = shutil.which("nextflow")
    assert nextflow is not None
    risk_summary = write_csv(
        tmp_path / "sample-1.with-risk-groups.csv",
        [
            {
                "rank": "species",
                "lineage": "Viruses;RG2 source;RG2 signal",
                "fraction": "0.1",
                "f_weighted_at_rank": "0.1",
                "bp_match_at_rank": "500",
                "who_risk_group": "RG2",
                "who_risk_group_source_taxid": "222",
                "who_risk_group_source_name": "RG2 source",
            },
        ],
    )
    enrichment = tmp_path / "sample-1.deacon_filter.json"
    enrichment.write_text("{}\n", encoding="utf-8")
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

params.experimental = true
params.slack_enabled = true
params.slack_channel = 'C123'

include {{ NOTIFY_RAPID_SCREENING_SLACK }} from '{NOTIFICATIONS_MODULE}'

workflow {{
    NOTIFY_RAPID_SCREENING_SLACK(
        Channel.of(tuple(
            'sample-1',
            file('{risk_summary}'),
            file('{enrichment}'),
        )),
        false,
        'blue_koala',
    )
    NOTIFY_RAPID_SCREENING_SLACK.out.done.view()
}}
""",
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment.pop("SLACK_BOT_TOKEN", None)
    environment["PATH"] = f"{ROOT / 'bin'}{os.pathsep}{environment['PATH']}"
    environment["NXF_HOME"] = str(tmp_path / "nxf-home")
    secret_setup = subprocess.run(  # noqa: S603
        [nextflow, "secrets", "set", "SLACK_BOT_TOKEN", "xoxb-test"],
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert secret_setup.returncode == 0, secret_setup.stderr

    completed = subprocess.run(  # noqa: S603
        [nextflow, "-C", "/dev/null", "run", str(workflow)],
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert "[sample-1, true]" in diagnostics
    command_files = list(tmp_path.glob("work/**/.command.sh"))
    assert len(command_files) == 1
    command = command_files[0].read_text(encoding="utf-8")
    assert "--minimum-risk-group RG3" in command
