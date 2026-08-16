"""Nextflow contract test for the run-completion Slack notification."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NOTIFICATIONS_MODULE = ROOT / "modules" / "notifications.nf"
REPORTING_SUBWORKFLOW = ROOT / "subworkflows" / "reporting.nf"
LIMS_INTEGRATION_SUBWORKFLOW = ROOT / "subworkflows" / "lims_integration.nf"


def test_reporting_labkey_completion_gate_uses_emitted_lims_output() -> None:
    reporting_source = REPORTING_SUBWORKFLOW.read_text(encoding="utf-8")
    lims_source = LIMS_INTEGRATION_SUBWORKFLOW.read_text(encoding="utf-8")

    lims_emit_block = lims_source.split("    emit:\n", maxsplit=1)[1]
    lims_outputs = set(
        re.findall(r"^    ([A-Za-z_][A-Za-z0-9_]*)\s=", lims_emit_block, re.MULTILINE)
    )
    reporting_lims_refs = set(
        re.findall(r"LIMS_INTEGRATION\.out\.([A-Za-z_][A-Za-z0-9_]*)", reporting_source)
    )

    assert "LIMS_INTEGRATION.out.uploads_done" in reporting_source
    assert reporting_lims_refs <= lims_outputs


@pytest.mark.parametrize("labkey_enabled", [False, True])
def test_completion_notifier_waits_for_final_gate_with_optional_labkey(
    tmp_path: Path,
    *,
    labkey_enabled: bool,
) -> None:
    nextflow = shutil.which("nextflow")
    assert nextflow is not None
    results = tmp_path / "experiment_blast_results.tsv"
    results.write_text(
        "sample\tqseqid\twho_risk_group\n",
        encoding="utf-8",
    )
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_notifier = fake_bin / "notify_slack.py"
    fake_notifier.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    fake_notifier.chmod(0o755)
    labkey_config = (
        """\
params.labkey = true
params.labkey_server = 'labkey.example.org'
params.labkey_project_name = 'dho/projects/example'
params.labkey_blast_meta_hits_list = "O'Connor BLAST hits"
"""
        if labkey_enabled
        else "params.labkey = false\n"
    )
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

params.slack_enabled = true
params.slack_channel = 'C123'
params.experiment_id = {"'12345'" if labkey_enabled else "null"}
{labkey_config}

include {{ NOTIFY_RUN_COMPLETION_SLACK }} from '{NOTIFICATIONS_MODULE}'

process FINAL_GATE {{
    output:
    val true

    script:
    '''
    sleep 1
    '''
}}

process FINAL_RESULTS {{
    input:
    val ready
    path source_results, stageAs: 'source.tsv'

    output:
    path 'experiment_blast_results.tsv'

    script:
    \"\"\"
    cp "${{source_results}}" experiment_blast_results.tsv
    \"\"\"
}}

workflow {{
    FINAL_GATE()
    FINAL_RESULTS(FINAL_GATE.out, Channel.value(file('{results}')))
    NOTIFY_RUN_COMPLETION_SLACK(
        FINAL_RESULTS.out,
        Channel.value('abc123'),
        'blue_koala',
    )
    NOTIFY_RUN_COMPLETION_SLACK.out.done.view()
}}
""",
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment.pop("SLACK_BOT_TOKEN", None)
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
    assert "true" in diagnostics
    command_files = list(tmp_path.glob("work/**/.command.sh"))
    notification_commands = [
        path.read_text(encoding="utf-8")
        for path in command_files
        if "notify_slack.py" in path.read_text(encoding="utf-8")
    ]
    assert len(notification_commands) == 1
    command = notification_commands[0]
    assert "--experiment-results experiment_blast_results.tsv" in command
    assert "--sample-set-id 'abc123'" in command
    if labkey_enabled:
        assert "--experiment-id-base64" in command
        assert "--labkey-config-base64" in command
    else:
        assert "--experiment-id-base64" not in command
        assert "--labkey-config-base64" not in command
