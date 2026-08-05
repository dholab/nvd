import csv
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent
SCRIPT = BIN / "prepare_blast_labkey.py"


def test_query_class_column_survives_into_labkey_csv(tmp_path: Path) -> None:
    tsv = tmp_path / "in.tsv"
    tsv.write_text("sample\tqseqid\tquery_class\ttask\trank\n" "s1\tq1\tsingle_read\tmegablast\tspecies:x\n")
    out = tmp_path / "out.csv"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--blast-csv", str(tsv), "--output", str(out),
         "--meta", "s1", "--experiment-id", "7"],
        check=True,
    )
    [row] = list(csv.DictReader(out.open()))
    assert row["experiment"] == "7"
    assert row["query_class"] == "single_read"      # passthrough, not renamed
    assert row["sample_id"] == "s1"                  # 'sample' -> 'sample_id' rename
    assert row["blast_task"] == "megablast"          # 'task' -> 'blast_task' rename
