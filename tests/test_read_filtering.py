"""Focused behavioral tests for read quality, length, and entropy filtering."""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BBMAP_MODULE = ROOT / "modules" / "bbmap"
PREPROCESS_READS = ROOT / "subworkflows" / "preprocess_reads"

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None or shutil.which("bbduk.sh") is None,
    reason="FILTER_READS behavior requires Nextflow and locked BBTools",
)

HIGH_COMPLEXITY = "ACGTTGCAGATCGTAGCTAGGCTAACGATCGTACGGCATTCGATGCTAGCATCGTACGTA"
LOW_COMPLEXITY = "T" * len(HIGH_COMPLEXITY)
DINUCLEOTIDE_REPEAT = "AC" * (len(HIGH_COMPLEXITY) // 2)
ILLUMINA_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
SECOND_MATE = 2


def write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    """Write named sequence and quality triples as gzip-compressed FASTQ."""
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for name, sequence, quality in records:
            handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


def run_filter(
    tmp_path: Path,
    records: list[tuple[str, str, str]],
    *,
    filter_reads: bool,
    filter_low_complexity_reads: bool,
    read_structure: str = "single",
    min_read_length: int = 50,
    min_read_entropy: float = 0.5,
) -> str:
    """Run the public FILTER_READS Nextflow process and return its FASTQ text."""
    reads = tmp_path / "reads.fastq.gz"
    write_fastq(reads, records)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

params.filter_reads = {str(filter_reads).lower()}
params.filter_low_complexity_reads = {str(filter_low_complexity_reads).lower()}
params.min_read_entropy = {min_read_entropy}
params.min_read_length = {min_read_length}
params.max_read_length = null
params.min_read_quality_illumina = 20
params.min_read_quality_nanopore = 20

include {{ FILTER_READS }} from '{BBMAP_MODULE}'

workflow {{
    meta = [id: 'sample_A', platform: 'illumina', read_structure: '{read_structure}', query_class: 'single_read']
    FILTER_READS(Channel.of(tuple(meta, file('{reads}'))))
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    completed = subprocess.run(  # noqa: S603
        ["nextflow", "-C", "/dev/null", "run", str(workflow)],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=45,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    outputs = list(
        (tmp_path / "work").glob("**/sample_A.single_read.filtered.fastq.gz")
    )
    assert len(outputs) == 1, outputs
    with gzip.open(outputs[0], "rt", encoding="utf-8") as handle:
        return handle.read()


def run_dedup(
    tmp_path: Path,
    records: list[tuple[str, str, str]],
) -> str:
    """Run the public DEDUP_WITH_CLUMPIFY process and return its FASTQ text."""
    reads = tmp_path / "reads.fastq.gz"
    write_fastq(reads, records)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ DEDUP_WITH_CLUMPIFY }} from '{BBMAP_MODULE}'

workflow {{
    meta = [id: 'sample_A', platform: 'illumina', read_structure: 'single', query_class: 'single_read']
    DEDUP_WITH_CLUMPIFY(Channel.of(tuple(meta, file('{reads}'))))
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    completed = subprocess.run(  # noqa: S603
        ["nextflow", "-C", "/dev/null", "run", str(workflow)],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=45,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    outputs = list(
        (tmp_path / "work").glob("**/sample_A.single_read.dedup.fastq.gz")
    )
    assert len(outputs) == 1, outputs
    with gzip.open(outputs[0], "rt", encoding="utf-8") as handle:
        return handle.read()


@pytest.mark.skipif(
    shutil.which("clumpify.sh") is None,
    reason="DEDUP_WITH_CLUMPIFY behavior requires locked BBTools",
)
def test_clumpify_deduplicates_only_exact_sequences_independent_of_name(
    tmp_path: Path,
) -> None:
    """Different names collapse only when their read sequences are exact duplicates."""
    near_match = HIGH_COMPLEXITY[:-1] + ("A" if HIGH_COMPLEXITY[-1] != "A" else "C")
    output = run_dedup(
        tmp_path,
        [
            ("first-name", HIGH_COMPLEXITY, "I" * len(HIGH_COMPLEXITY)),
            ("second-name", HIGH_COMPLEXITY, "I" * len(HIGH_COMPLEXITY)),
            ("near-match", near_match, "I" * len(near_match)),
        ],
    )
    sequences = output.splitlines()[1::4]

    assert sorted(sequences) == sorted([HIGH_COMPLEXITY, near_match])


def run_adapter_trim(
    tmp_path: Path,
    records: list[tuple[str, str, str]],
    *,
    read_structure: str,
) -> str:
    """Run the public TRIM_ADAPTERS process and return its FASTQ text."""
    reads = tmp_path / "reads.fastq.gz"
    write_fastq(reads, records)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ TRIM_ADAPTERS }} from '{BBMAP_MODULE}'

workflow {{
    meta = [id: 'sample_A', platform: 'illumina', read_structure: '{read_structure}', query_class: 'single_read']
    TRIM_ADAPTERS(Channel.of(tuple(meta, file('{reads}'))))
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    completed = subprocess.run(  # noqa: S603
        ["nextflow", "-C", "/dev/null", "run", str(workflow)],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=45,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    outputs = list(
        (tmp_path / "work").glob("**/sample_A.single_read.trimmed.fastq.gz")
    )
    assert len(outputs) == 1, outputs
    with gzip.open(outputs[0], "rt", encoding="utf-8") as handle:
        return handle.read()


@pytest.mark.parametrize("adapter_mate", [1, SECOND_MATE])
def test_adapter_trimming_removes_complete_interleaved_pair_when_either_mate_matches(
    tmp_path: Path,
    adapter_mate: int,
) -> None:
    """An adapter-only mate removes its complete SRA-style interleaved pair."""
    records = [
        (
            "spot.1",
            ILLUMINA_ADAPTER if adapter_mate == 1 else HIGH_COMPLEXITY,
            "I"
            * (len(ILLUMINA_ADAPTER) if adapter_mate == 1 else len(HIGH_COMPLEXITY)),
        ),
        (
            "spot.2",
            ILLUMINA_ADAPTER if adapter_mate == SECOND_MATE else HIGH_COMPLEXITY,
            "I"
            * (
                len(ILLUMINA_ADAPTER)
                if adapter_mate == SECOND_MATE
                else len(HIGH_COMPLEXITY)
            ),
        ),
    ]

    output = run_adapter_trim(tmp_path, records, read_structure="interleaved")

    assert output == ""


def test_adapter_trimming_keeps_single_read_batches_unpaired(tmp_path: Path) -> None:
    """SRA-style suffixes do not make independent single reads into a pair."""
    output = run_adapter_trim(
        tmp_path,
        [
            ("spot.1", ILLUMINA_ADAPTER, "I" * len(ILLUMINA_ADAPTER)),
            ("spot.2", HIGH_COMPLEXITY, "I" * len(HIGH_COMPLEXITY)),
        ],
        read_structure="single",
    )

    assert "@spot.1\n" not in output
    assert "@spot.2\n" in output


@pytest.mark.skipif(
    any(shutil.which(tool) is None for tool in ("clumpify.sh", "deacon")),
    reason="preprocessing order requires locked Clumpify and Deacon",
)
def test_preprocessing_deduplicates_after_adapter_trimming(tmp_path: Path) -> None:
    """Reads made identical by adapter removal collapse before preprocessing ends."""
    lib_dir = tmp_path / "lib"
    lib_dir.mkdir()
    shutil.copy2(ROOT / "lib" / "NvdUtils.groovy", lib_dir)
    shutil.copy2(ROOT / "lib" / "NvdReporting.groovy", lib_dir)
    assets_dir = tmp_path / "assets"
    assets_dir.mkdir()
    shutil.copy2(ROOT / "assets" / "empty_deacon.k31w1.idx", assets_dir)

    reads = tmp_path / "reads.fastq.gz"
    write_fastq(
        reads,
        [
            (
                "adapter-bearing",
                HIGH_COMPLEXITY + ILLUMINA_ADAPTER,
                "I" * (len(HIGH_COMPLEXITY) + len(ILLUMINA_ADAPTER)),
            ),
            ("already-trimmed", HIGH_COMPLEXITY, "I" * len(HIGH_COMPLEXITY)),
        ],
    )

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

params.no_enrichment = true
params.virus_index = null
params.virus_index_url = null
params.virus_reference_fasta = null
params.virus_abs_threshold = 1
params.virus_rel_threshold = 0.0
params.merge_pairs = false
params.dedup = false
params.dedup_seq = true
params.trim_adapters = true
params.host_index = null
params.host_index_url = null
params.host_contaminants_fasta = null
params.filter_reads = false
params.filter_low_complexity_reads = false
params.min_read_quality_illumina = 20
params.min_read_quality_nanopore = 12
params.min_read_length = 50
params.min_consecutive_bases = 200

include {{ PREPROCESS_READS }} from '{PREPROCESS_READS}'

workflow {{
    PREPROCESS_READS(
        Channel.of(tuple(
            [
                id: 'sample_A',
                platform: 'illumina',
                read_mode: 'single',
                r1_count: 1,
                deacon_read_structure: 'single',
            ],
            [file('{reads}')],
        )),
        Channel.empty(),
    )

    PREPROCESS_READS.out.read_batches.view {{ _id, _platform, _structure, _query_class, output ->
        "FINAL_READS:${{output.name}}"
    }}
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["PATH"] = (
        f"{ROOT / 'bin'}{os.pathsep}{environment['PATH']}"
    )
    environment["NXF_ANSI_LOG"] = "false"
    completed = subprocess.run(  # noqa: S603
        ["nextflow", "-C", "/dev/null", "run", str(workflow)],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=90,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert "FINAL_READS:sample_A.single_read.dedup.fastq.gz" in completed.stdout
    [output] = [
        path
        for path in (tmp_path / "work").glob(
            "**/sample_A.single_read.dedup.fastq.gz",
        )
        if not path.is_symlink()
    ]
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        records = handle.read().splitlines()
    assert len(records) == 4
    assert records[1] == HIGH_COMPLEXITY


def test_entropy_filter_does_not_enable_quality_filtering(tmp_path: Path) -> None:
    """Entropy-only filtering keeps a complex read despite low base quality."""
    output = run_filter(
        tmp_path,
        [
            ("complex-low-quality", HIGH_COMPLEXITY, "!" * len(HIGH_COMPLEXITY)),
            ("homopolymer-high-quality", LOW_COMPLEXITY, "I" * len(LOW_COMPLEXITY)),
            (
                "dinucleotide-high-quality",
                DINUCLEOTIDE_REPEAT,
                "I" * len(DINUCLEOTIDE_REPEAT),
            ),
        ],
        filter_reads=False,
        filter_low_complexity_reads=True,
        min_read_length=len(HIGH_COMPLEXITY) + 1,
    )

    assert "@complex-low-quality\n" in output
    assert "@homopolymer-high-quality\n" not in output
    assert "@dinucleotide-high-quality\n" not in output


def test_quality_filter_does_not_enable_entropy_filtering(tmp_path: Path) -> None:
    """Quality-only filtering keeps a high-quality low-complexity read."""
    output = run_filter(
        tmp_path,
        [
            ("complex-low-quality", HIGH_COMPLEXITY, "!" * len(HIGH_COMPLEXITY)),
            ("homopolymer-high-quality", LOW_COMPLEXITY, "I" * len(LOW_COMPLEXITY)),
        ],
        filter_reads=True,
        filter_low_complexity_reads=False,
    )

    assert "@complex-low-quality\n" not in output
    assert "@homopolymer-high-quality\n" in output


def test_combined_filters_can_emit_valid_empty_fastq(tmp_path: Path) -> None:
    """An all-rejected batch remains a valid empty gzip output."""
    output = run_filter(
        tmp_path,
        [
            ("complex-low-quality", HIGH_COMPLEXITY, "!" * len(HIGH_COMPLEXITY)),
            ("homopolymer-high-quality", LOW_COMPLEXITY, "I" * len(LOW_COMPLEXITY)),
        ],
        filter_reads=True,
        filter_low_complexity_reads=True,
    )

    assert output == ""


def test_entropy_filter_removes_complete_interleaved_pair(tmp_path: Path) -> None:
    """An entropy failure removes both mates from an intact interleaved pair."""
    output = run_filter(
        tmp_path,
        [
            ("pair/1", HIGH_COMPLEXITY, "I" * len(HIGH_COMPLEXITY)),
            ("pair/2", LOW_COMPLEXITY, "I" * len(LOW_COMPLEXITY)),
        ],
        filter_reads=False,
        filter_low_complexity_reads=True,
        read_structure="interleaved",
    )

    assert output == ""


@pytest.mark.parametrize(
    ("min_read_entropy", "expected_ids"),
    [
        (
            0.5,
            {"000113", "000139", "000141", "000144", "000146", "control"},
        ),
        (0.9, {"000144", "control"}),
    ],
    ids=["moderate-default", "strict-compatibility"],
)
def test_entropy_thresholds_filter_observed_repeat_rich_reads(
    tmp_path: Path,
    min_read_entropy: float,
    expected_ids: set[str],
) -> None:
    """The moderate default removes strong repeats while 0.9 remains available."""
    observed = [
        (
            "000110",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        ),
        (
            "000113",
            "TCCTCTCTCTCTTCCTCTCTCTCTTCCTCTTCTACTTTTGTTGACGACGAAGAAATCGAGAGACAGATACGAGAGGAGCTGAGGTGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA",
        ),
        (
            "000138",
            "TAATCCGCCCACACACACACATACACATACACACACCCACACCCACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
        ),
        (
            "000139",
            "TGGCTGATGGCGTCATGCGTTCGTTTGTTATTTTGTTGTTGTTTTTCGTTGTTTTCTTTTGTTGGTTCAAGGATACGATGGTGGTCTTCCTTGTGTTTTTCTCTTTTTTTCTTTTTTTCTTTTTGTTTTCATCTTTTTTTTTTTTTT",
        ),
        (
            "000141",
            "TTGAAAGCAAGACGAAGCAAAACAGAATGACGGGAAAAGAGGCTCTTCCAAGAAGAAACCAGGAAAAACGCACGTGGTTTGTTGCTTTTTTCGTTTTGGCGCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        ),
        (
            "000144",
            "TGTGTGATGTGTGCGTACGTCTGTGCGTCTGACATCCCTAGGGTCGGCATCAGCCTGTGCTTGCCCGCCGCGTGCGCCCAGGACATCGAGGAGGGCTACGAGTTCCTCCTCCCTTACCTCGCGTACGCACACACACACACACACACA",
        ),
        (
            "000145",
            "TCCCGTCCATTATCACTCATTTTTTGCTCATTGTTGTTGTTGTTGTTGCTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG",
        ),
        (
            "000146",
            "TCGGAGAAATTTATTGGCTAAATGCCAGTATAAAATGGGCGAAGTGTCAAACCTTAAAAAGATAGAAACGAATGAAAAAAGGAAAGAACTGACAGTTGAGTGAAAGAATGAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        ),
        (
            "control",
            "ACGTTGCAGATCGTAGCTAGGCTAACGATCGTACGGCATTCGATGCTAGCATCGTACGTATCGATGGCATACGTAGCTGATCGTACGATCGTAGCTAGCATGCTACGATCGTACGTAGCATCGATGCTAGCATCGTAGCTAGCATGCA",
        ),
    ]
    output = run_filter(
        tmp_path,
        [(name, sequence, "I" * len(sequence)) for name, sequence in observed],
        filter_reads=False,
        filter_low_complexity_reads=True,
        min_read_entropy=min_read_entropy,
    )
    retained_ids = {
        line.removeprefix("@")
        for index, line in enumerate(output.splitlines())
        if index % 4 == 0
    }

    assert retained_ids == expected_ids


def test_either_filter_family_schedules_filter_reads() -> None:
    """Either read-filter family should schedule the shared FILTER_READS process."""
    preprocessing = (ROOT / "subworkflows" / "preprocess_reads.nf").read_text(
        encoding="utf-8",
    )

    assert "params.filter_reads || params.filter_low_complexity_reads" in preprocessing
