#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["loguru>=0.7.3", "polars>=1.27.1"]
# ///
"""Prepare WVDB inputs for a combined sourmash reference database."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

import polars as pl
from loguru import logger

WVDB_ID_PREFIX = "WVDB|"
TAXONOMY_RANKS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
REFERENCE_TAXONOMY_RANKS = (*TAXONOMY_RANKS, "strain")
CURATION_TAXONOMY_RANKS = REFERENCE_TAXONOMY_RANKS
TAXPATH_COLUMN = "taxpath"
LINEAGE_COLUMNS = ("ident", *TAXONOMY_RANKS, TAXPATH_COLUMN)
COMBINED_LINEAGE_COLUMNS = ("ident", *REFERENCE_TAXONOMY_RANKS, TAXPATH_COLUMN)
CURATED_LINEAGE_COLUMNS = COMBINED_LINEAGE_COLUMNS
LINEAGE_IDENTIFIER_COLUMNS = ("ident", "identifiers")
SIGNATURE_IDENTITY_COLUMNS = (
    "md5",
    "ksize",
    "moltype",
    "num",
    "scaled",
    "with_abundance",
)
TAXONOMY_CONFLICT_POLICIES = (
    "error",
    "most-specified",
    "ncbi-wins",
    "wvdb-wins",
)
UNKNOWN_TAXA = frozenset({"", "unknown", "na", "n/a", "none"})
DIAGNOSTIC_COLUMNS = (
    "conflict_type",
    "rank",
    "normalized_name",
    "policy",
    "decision",
    "source",
    "ident",
    "signature_identity",
    "specificity",
    "lineage",
)


@dataclass(frozen=True)
class CurationPaths:
    """Input and output paths for combined-reference taxonomy curation."""

    reference_manifest_csv: Path
    wvdb_manifest_csv: Path
    reference_lineages_csv: Path
    wvdb_lineages_csv: Path
    lineages_csv: Path
    diagnostics_tsv: Path


PLACEHOLDER_SENTINELS = frozenset({"unclassified", "unknown", "na", "n/a", "none"})
VIRUSES_TAXID = "10239"
WVDB_LOCAL_TAXID_BASE = 900_000_000_000
WVDB_LOCAL_TAXID_MODULUS = 99_999_999_999


@dataclass(frozen=True)
class ReferenceLineageIndex:
    """Indexes used to map WVDB lineage names onto reference taxpaths."""

    rows: list[dict[str, str]]
    rows_by_rank_name: dict[tuple[str, str], list[dict[str, str]]]
    prefix_taxids: dict[tuple[tuple[str, str], ...], str]


def configure_logging() -> None:
    """Configure Loguru for build progress on stderr."""
    logger.remove()
    logger.add(
        sys.stderr,
        format="<level>{level: <8}</level> | {message}",
        level="INFO",
    )


def ensure_parent(path: Path) -> None:
    """Create the parent directory for an output path when needed."""
    path.parent.mkdir(parents=True, exist_ok=True)


def cleaned_string(name: str) -> pl.Expr:
    """Return a stripped string column expression with blanks converted to null."""
    value = pl.col(name).cast(pl.Utf8).str.strip_chars()
    return pl.when(value.str.len_chars() > 0).then(value).otherwise(None)


def informative_string(name: str) -> pl.Expr:
    """Return a stripped source value, with unknown sentinels treated as null."""
    value = cleaned_string(name)
    return (
        pl.when(value.str.to_lowercase().is_in(PLACEHOLDER_SENTINELS))
        .then(None)
        .otherwise(value)
    )


def sentinel_string(name: str) -> pl.Expr:
    """Return a stripped source value only when it is an unknown sentinel."""
    value = cleaned_string(name)
    return (
        pl.when(value.str.to_lowercase().is_in(PLACEHOLDER_SENTINELS))
        .then(value)
        .otherwise(None)
    )


def first_nonempty(available_columns: set[str], *columns: str) -> pl.Expr:
    """Return the first informative value among columns present in a frame.

    Unknown sentinels such as ``Unclassified`` fall through to later informative
    source columns. If no informative value exists, the first sentinel is kept so
    downstream lineage construction can preserve the explicit unknown rank.
    """
    present_columns = [name for name in columns if name in available_columns]
    values = [informative_string(name) for name in present_columns]
    sentinel_values = [sentinel_string(name) for name in present_columns]
    if not values:
        return pl.lit("")
    return pl.coalesce([*values, *sentinel_values]).fill_null("")


def normalize_fasta_ids(source: Path, destination: Path) -> int:
    """Write a FASTA with WVDB-prefixed, whitespace-free sequence identifiers."""
    ensure_parent(destination)
    seen: set[str] = set()

    with (
        source.open(encoding="utf-8") as src,
        destination.open("w", encoding="utf-8") as dst,
    ):
        for line_number, line in enumerate(src, start=1):
            if not line.startswith(">"):
                dst.write(line)
                continue

            raw_id = line[1:].strip().split(maxsplit=1)[0]
            if not raw_id:
                msg = f"Empty FASTA identifier at line {line_number} in {source}"
                raise ValueError(msg)

            normalized_id = f"{WVDB_ID_PREFIX}{raw_id}"
            if normalized_id in seen:
                msg = f"Duplicate FASTA identifier after normalization: {normalized_id}"
                raise ValueError(msg)
            seen.add(normalized_id)

            dst.write(f">{normalized_id}\n")

    return len(seen)


def scan_csv(path: Path, *, separator: str = ",", skip_rows: int = 0) -> pl.LazyFrame:
    """Read a CSV-like file lazily with all columns as strings."""
    return pl.scan_csv(
        path,
        separator=separator,
        skip_rows=skip_rows,
        infer_schema_length=0,
        null_values=[""],
    )


def require_lazy_columns(
    frame: pl.LazyFrame,
    *,
    path: Path,
    columns: set[str],
) -> pl.LazyFrame:
    """Fail if a lazy CSV frame lacks required columns."""
    missing = columns.difference(frame.collect_schema().names())
    if missing:
        msg = f"{path} is missing required column(s): {', '.join(sorted(missing))}"
        raise ValueError(msg)
    return frame


def normalize_lineage_identifier_column(
    frame: pl.LazyFrame,
    *,
    path: Path,
) -> pl.LazyFrame:
    """Expose raw and prepared sourmash lineage identifiers as ident."""
    available_columns = set(frame.collect_schema().names())
    for column in LINEAGE_IDENTIFIER_COLUMNS:
        if column not in available_columns:
            continue
        if column == "ident":
            return frame
        return frame.rename({column: "ident"})

    expected = ", ".join(LINEAGE_IDENTIFIER_COLUMNS)
    msg = f"{path} is missing required lineage identifier column (expected one of: {expected})"
    raise ValueError(msg)


def scan_annotations(annotations_tsv: Path) -> pl.LazyFrame:
    """Read the WVDB annotation TSV lazily as strings."""
    return scan_csv(
        annotations_tsv,
        separator="\t",
    ).pipe(require_lazy_columns, path=annotations_tsv, columns={"contig"})


def scan_sourmash_manifest(manifest_csv: Path) -> pl.LazyFrame:
    """Read a sourmash manifest CSV, skipping its metadata preamble."""
    return scan_csv(manifest_csv, skip_rows=1).pipe(
        require_lazy_columns,
        path=manifest_csv,
        columns={"name"},
    )


def require_unique_identifiers(lineages: pl.DataFrame) -> pl.DataFrame:
    """Fail if lineage identifiers are duplicated."""
    duplicates = (
        lineages.group_by("ident")
        .len()
        .filter(pl.col("len") > 1)
        .select("ident")
        .get_column("ident")
        .to_list()
    )
    if duplicates:
        examples = ", ".join(duplicates[:10])
        msg = f"Duplicate WVDB lineage identifiers: {examples}"
        raise ValueError(msg)
    return lineages


def split_taxpath(taxpath: str) -> list[str]:
    """Split a sourmash taxpath while preserving internal blank ranks."""
    taxids = clean_cell(taxpath).split("|") if clean_cell(taxpath) else []
    while taxids and not taxids[-1]:
        taxids.pop()
    return taxids


def require_reference_lineage_columns(
    rows: list[dict[str, str]],
    *,
    path: Path,
) -> list[dict[str, str]]:
    """Fail if reference lineages lack columns needed for taxpath mapping."""
    if not rows:
        msg = f"{path} contains no reference lineage rows"
        raise ValueError(msg)

    columns = set(rows[0])
    missing = {TAXPATH_COLUMN, *TAXONOMY_RANKS}.difference(columns)
    if missing:
        msg = f"{path} is missing required column(s): {', '.join(sorted(missing))}"
        raise ValueError(msg)
    return rows


def load_reference_lineages(path: Path | None) -> list[dict[str, str]]:
    """Load reference lineages used to reuse known taxpaths."""
    if path is None:
        return []
    return require_reference_lineage_columns(read_csv_dicts(path), path=path)


def reference_row_taxids(row: dict[str, str]) -> dict[str, str]:
    """Return reference taxids keyed by taxonomy rank."""
    taxids = split_taxpath(row.get(TAXPATH_COLUMN, ""))
    return {
        rank: taxids[index]
        for index, rank in enumerate(REFERENCE_TAXONOMY_RANKS)
        if index < len(taxids) and taxids[index]
    }


def reference_prefix(row: dict[str, str], *, through_rank: str) -> tuple[str, ...]:
    """Return a reference lineage's names through a rank."""
    rank_index = TAXONOMY_RANKS.index(through_rank)
    return tuple(clean_cell(row.get(rank)) for rank in TAXONOMY_RANKS[: rank_index + 1])


def non_placeholder_value(row: dict[str, str], rank: str) -> str:
    """Return a rank value only when it is an informative taxon name."""
    value = clean_cell(row.get(rank))
    if is_placeholder_sentinel(value, rank=rank):
        return ""
    return value


def reference_rows_by_rank_name(
    reference_rows: list[dict[str, str]],
) -> dict[tuple[str, str], list[dict[str, str]]]:
    """Index reference rows by exact rank/name values."""
    rows_by_rank_name: defaultdict[tuple[str, str], list[dict[str, str]]] = defaultdict(
        list,
    )
    for row in reference_rows:
        for rank in TAXONOMY_RANKS:
            name = clean_cell(row.get(rank))
            if name:
                rows_by_rank_name[(rank, name)].append(row)
    return dict(rows_by_rank_name)


def build_reference_lineage_index(
    reference_rows: list[dict[str, str]],
) -> ReferenceLineageIndex:
    """Build reusable lookup indexes for reference lineage matching."""
    return ReferenceLineageIndex(
        rows=reference_rows,
        rows_by_rank_name=reference_rows_by_rank_name(reference_rows),
        prefix_taxids=reference_prefix_taxids(reference_rows),
    )


def candidate_matches_provided_values(
    candidate: dict[str, str],
    values: dict[str, str],
) -> bool:
    """Return whether a reference row agrees with all informative WVDB ranks."""
    return all(
        not value or clean_cell(candidate.get(rank)) == value
        for rank, value in values.items()
    )


def unique_reference_prefix_candidate(
    row: dict[str, str],
    reference_index: ReferenceLineageIndex,
) -> tuple[str, dict[str, str]] | None:
    """Return a unique reference lineage prefix compatible with a WVDB row."""
    if not reference_index.rows:
        return None

    informative_values = {
        rank: non_placeholder_value(row, rank)
        for rank in TAXONOMY_RANKS
        if non_placeholder_value(row, rank)
    }
    if not informative_values:
        return None

    for deepest_rank in reversed(TAXONOMY_RANKS):
        deepest_name = informative_values.get(deepest_rank)
        if not deepest_name:
            continue

        candidate_values = {
            rank: value
            for rank, value in informative_values.items()
            if TAXONOMY_RANKS.index(rank) <= TAXONOMY_RANKS.index(deepest_rank)
        }
        candidates = [
            candidate
            for candidate in reference_index.rows_by_rank_name.get(
                (deepest_rank, deepest_name),
                [],
            )
            if candidate_matches_provided_values(candidate, candidate_values)
        ]
        if not candidates:
            continue

        unique_prefixes: dict[tuple[str, ...], dict[str, str]] = {}
        for candidate in candidates:
            unique_prefixes.setdefault(
                reference_prefix(candidate, through_rank=deepest_rank),
                candidate,
            )

        if len(unique_prefixes) == 1:
            return deepest_rank, next(iter(unique_prefixes.values()))

    return None


def backfill_wvdb_lineage_from_reference(
    row: dict[str, str],
    reference_index: ReferenceLineageIndex,
) -> dict[str, str]:
    """Fill missing or placeholder WVDB ancestors from a unique reference path."""
    match = unique_reference_prefix_candidate(row, reference_index)
    if match is None:
        return row

    deepest_rank, candidate = match
    output = dict(row)
    for rank in TAXONOMY_RANKS[: TAXONOMY_RANKS.index(deepest_rank) + 1]:
        value = clean_cell(output.get(rank))
        if value and not is_placeholder_sentinel(value, rank=rank):
            continue
        reference_value = clean_cell(candidate.get(rank))
        if reference_value:
            output[rank] = reference_value
    return output


def backfill_wvdb_lineages_from_reference(
    lineages: pl.DataFrame,
    reference_index: ReferenceLineageIndex,
) -> pl.DataFrame:
    """Backfill WVDB missing ancestors when the reference path is unique."""
    if not reference_index.rows:
        logger.info("reference_backfill.skipped reason=no_reference_lineages")
        return lineages

    rows = []
    backfilled_rows = 0
    for row in lineages.iter_rows(named=True):
        input_row = {str(key): clean_cell(value) for key, value in row.items()}
        output_row = backfill_wvdb_lineage_from_reference(input_row, reference_index)
        if any(input_row.get(rank) != output_row.get(rank) for rank in TAXONOMY_RANKS):
            backfilled_rows += 1
        rows.append(output_row)

    logger.info(
        "reference_backfill.done rows={} backfilled_rows={}",
        len(rows),
        backfilled_rows,
    )
    return pl.DataFrame(rows).select(lineages.columns)


def reference_prefix_taxids(
    reference_rows: list[dict[str, str]],
) -> dict[tuple[tuple[str, str], ...], str]:
    """Return exact unambiguous reference prefix-to-taxid mappings."""
    taxids_by_prefix: defaultdict[tuple[tuple[str, str], ...], set[str]] = defaultdict(
        set,
    )
    for row in reference_rows:
        row_taxids = reference_row_taxids(row)
        prefix: list[tuple[str, str]] = []
        for rank in TAXONOMY_RANKS:
            name = clean_cell(row.get(rank))
            taxid = row_taxids.get(rank)
            if not name:
                continue
            prefix.append((rank, name))
            if taxid:
                taxids_by_prefix[tuple(prefix)].add(taxid)

    mappings = {
        prefix: next(iter(taxids))
        for prefix, taxids in taxids_by_prefix.items()
        if len(taxids) == 1
    }
    mappings[(("superkingdom", "Viruses"),)] = VIRUSES_TAXID
    return mappings


def wvdb_local_taxid(prefix: tuple[tuple[str, str], ...]) -> str:
    """Return a stable numeric WVDB-local taxid for a lineage prefix."""
    canonical = "|".join(f"{rank}={name}" for rank, name in prefix)
    digest = hashlib.sha1(canonical.encode("utf-8"), usedforsecurity=False)
    offset = int(digest.hexdigest()[:12], 16) % WVDB_LOCAL_TAXID_MODULUS
    return str(WVDB_LOCAL_TAXID_BASE + offset)


def taxpath_for_row(
    row: dict[str, str],
    reference_taxids: dict[tuple[tuple[str, str], ...], str],
    *,
    reference_nodes: set[tuple[tuple[str, str], ...]] | None = None,
    local_nodes: set[tuple[tuple[str, str], ...]] | None = None,
) -> str:
    """Return the sourmash taxpath for a finalized WVDB lineage row."""
    prefix: list[tuple[str, str]] = []
    taxids: list[str] = []
    for rank in TAXONOMY_RANKS:
        name = clean_cell(row.get(rank))
        if not name:
            continue
        prefix.append((rank, name))
        prefix_key = tuple(prefix)
        taxid = reference_taxids.get(prefix_key)
        if taxid is not None:
            if reference_nodes is not None:
                reference_nodes.add(prefix_key)
        else:
            taxid = wvdb_local_taxid(prefix_key)
            if local_nodes is not None:
                local_nodes.add(prefix_key)
        taxids.append(taxid)
    return "|".join(taxids)


def add_wvdb_taxpaths(
    lineages: pl.DataFrame,
    reference_index: ReferenceLineageIndex,
) -> pl.DataFrame:
    """Add BioBoxes-compatible taxpaths to WVDB lineages."""
    rows = []
    reference_nodes: set[tuple[tuple[str, str], ...]] = set()
    local_nodes: set[tuple[tuple[str, str], ...]] = set()
    for row in lineages.iter_rows(named=True):
        output_row = {str(key): clean_cell(value) for key, value in row.items()}
        output_row[TAXPATH_COLUMN] = taxpath_for_row(
            output_row,
            reference_index.prefix_taxids,
            reference_nodes=reference_nodes,
            local_nodes=local_nodes,
        )
        rows.append(output_row)
    logger.info(
        "taxpath_generation.done rows={} reference_nodes={} local_nodes={}",
        len(rows),
        len(reference_nodes),
        len(local_nodes),
    )
    return pl.DataFrame(rows).select(LINEAGE_COLUMNS)


def clean_cell(value: object) -> str:
    """Return a stripped string value, treating nulls as empty strings."""
    if value is None:
        return ""
    return str(value).strip()


def bare_placeholder(rank: str) -> str:
    """Return the canonical pre-context placeholder name for a taxonomy rank."""
    return f"unclassified {rank}"


def is_bare_placeholder(name: str) -> bool:
    """Return whether a name is one of the canonical pre-context placeholders."""
    return name in {bare_placeholder(rank) for rank in TAXONOMY_RANKS}


def is_placeholder_sentinel(value: object, *, rank: str) -> bool:
    """Return whether a source value means this rank is explicitly unknown."""
    text = clean_cell(value)
    if not text:
        return False
    lowered = text.lower()
    return lowered in PLACEHOLDER_SENTINELS or lowered == bare_placeholder(rank)


def has_lower_rank_value(values: dict[str, str], rank_index: int) -> bool:
    """Return whether any lower rank has a source value."""
    return any(values[rank] for rank in TAXONOMY_RANKS[rank_index + 1 :])


def nearest_non_placeholder_context(parent: tuple[str, ...]) -> str:
    """Return the nearest classified ancestor name, falling back to root."""
    for name in reversed(parent):
        if not is_bare_placeholder(name):
            return name
    return "root"


def placeholder_contexts(parent: tuple[str, ...]) -> list[str]:
    """Return increasingly specific parent contexts for a placeholder node."""
    contexts: list[str] = []

    def add(context: str) -> None:
        if context and context not in contexts:
            contexts.append(context)

    add(nearest_non_placeholder_context(parent))
    for depth in range(1, len(parent) + 1):
        add(" > ".join(parent[-depth:]))
    if not contexts:
        add("root")
    return contexts


def short_hash(prefix: tuple[str, ...]) -> str:
    """Return a short deterministic suffix for pathological placeholder clashes."""
    digest = hashlib.sha1(";".join(prefix).encode("utf-8"), usedforsecurity=False)
    return digest.hexdigest()[:8]


def choose_placeholder_names(
    placeholder_nodes: dict[tuple[str, ...], str],
    *,
    reserved_names: set[str],
) -> dict[tuple[str, ...], str]:
    """Choose stable, contextual display names for WVDB placeholder nodes."""
    nodes_by_placeholder: defaultdict[str, list[tuple[str, ...]]] = defaultdict(list)
    for prefix in placeholder_nodes:
        nodes_by_placeholder[prefix[-1]].append(prefix)

    replacements: dict[tuple[str, ...], str] = {}
    reserved = set(reserved_names)
    for placeholder, prefixes in nodes_by_placeholder.items():
        contexts = {prefix: placeholder_contexts(prefix[:-1]) for prefix in prefixes}
        context_indexes = dict.fromkeys(prefixes, 0)

        while True:
            names = {
                prefix: f"{placeholder} [under {contexts[prefix][context_indexes[prefix]]}]"
                for prefix in prefixes
            }
            counts = Counter(names.values())
            colliding = {
                prefix
                for prefix, name in names.items()
                if counts[name] > 1 or name in reserved
            }
            if not colliding:
                replacements.update(names)
                reserved.update(names.values())
                break

            advanced = False
            for prefix in colliding:
                if context_indexes[prefix] < len(contexts[prefix]) - 1:
                    context_indexes[prefix] += 1
                    advanced = True

            if not advanced:
                resolved_names = {
                    prefix: f"{name} #{short_hash(prefix)}"
                    if prefix in colliding
                    else name
                    for prefix, name in names.items()
                }
                replacements.update(resolved_names)
                reserved.update(resolved_names.values())
                break

    return replacements


def contextualize_wvdb_placeholders(lineages: pl.DataFrame) -> pl.DataFrame:
    """Convert WVDB unknown-rank markers into stable contextual lineage names."""
    output_rows: list[dict[str, str]] = []
    prefixes_by_row: list[dict[str, tuple[str, ...]]] = []
    placeholder_nodes: dict[tuple[str, ...], str] = {}
    real_names: set[str] = set()

    for row in lineages.iter_rows(named=True):
        values = {rank: clean_cell(row.get(rank)) for rank in TAXONOMY_RANKS}
        output_row = {str(key): clean_cell(value) for key, value in row.items()}
        prefixes: dict[str, tuple[str, ...]] = {}
        lineage_prefix: list[str] = []

        for rank_index, rank in enumerate(TAXONOMY_RANKS):
            value = values[rank]
            is_placeholder = is_placeholder_sentinel(value, rank=rank) or (
                not value and has_lower_rank_value(values, rank_index)
            )

            if is_placeholder:
                label = bare_placeholder(rank)
            elif value:
                label = value
            else:
                label = ""

            output_row[rank] = label
            if not label:
                continue

            lineage_prefix.append(label)
            prefix = tuple(lineage_prefix)
            prefixes[rank] = prefix
            if is_placeholder:
                placeholder_nodes.setdefault(prefix, rank)
            else:
                real_names.add(label)

        output_rows.append(output_row)
        prefixes_by_row.append(prefixes)

    replacements = choose_placeholder_names(
        placeholder_nodes,
        reserved_names=real_names,
    )
    for output_row, prefixes in zip(output_rows, prefixes_by_row, strict=True):
        for rank, prefix in prefixes.items():
            replacement = replacements.get(prefix)
            if replacement is not None:
                output_row[rank] = replacement

    logger.info(
        "placeholder_context.done rows={} placeholder_nodes={} replacements={}",
        len(output_rows),
        len(placeholder_nodes),
        len(replacements),
    )
    return pl.DataFrame(output_rows).select(lineages.columns)


def select_wvdb_lineages(annotations: pl.LazyFrame) -> pl.LazyFrame:
    """Project WVDB annotation columns into sourmash lineage columns."""
    available_columns = set(annotations.collect_schema().names())
    return annotations.select(
        (
            pl.lit(WVDB_ID_PREFIX) + pl.col("contig").cast(pl.Utf8).str.strip_chars()
        ).alias("ident"),
        pl.lit("Viruses").alias("superkingdom"),
        first_nonempty(available_columns, "Phylum_RdRp", "Phylum_gNd").alias("phylum"),
        first_nonempty(available_columns, "Class_RdRp", "Class_gNd").alias("class"),
        first_nonempty(
            available_columns,
            "Order_consensus",
            "Order_RdRp",
            "Order_gNd",
        ).alias("order"),
        first_nonempty(
            available_columns,
            "Family_consensus",
            "Family_RdRp",
            "Family_gNd",
            "Family_ICTV",
        ).alias("family"),
        first_nonempty(available_columns, "Genus_RdRp").alias("genus"),
        first_nonempty(available_columns, "Species_RdRp").alias("species"),
    )


def build_wvdb_lineages(
    annotations: pl.LazyFrame,
    *,
    reference_lineages_csv: Path | None = None,
) -> pl.DataFrame:
    """Build WVDB sourmash lineages as a validated DataFrame."""
    started = perf_counter()
    reference_rows = load_reference_lineages(reference_lineages_csv)
    logger.info(
        "reference_lineages.loaded path={} rows={}",
        reference_lineages_csv or "none",
        len(reference_rows),
    )
    reference_index = build_reference_lineage_index(reference_rows)
    logger.info(
        "reference_lineages.indexed rank_name_keys={} prefix_taxids={}",
        len(reference_index.rows_by_rank_name),
        len(reference_index.prefix_taxids),
    )

    lineages = annotations.pipe(select_wvdb_lineages).collect()
    logger.info(
        "wvdb_lineages.projected rows={} columns={}",
        lineages.height,
        len(lineages.columns),
    )
    lineages = lineages.pipe(backfill_wvdb_lineages_from_reference, reference_index)
    lineages = lineages.pipe(contextualize_wvdb_placeholders)
    lineages = lineages.pipe(add_wvdb_taxpaths, reference_index)
    lineages = lineages.pipe(require_unique_identifiers)
    logger.info(
        "wvdb_lineages.built rows={} elapsed_seconds={:.2f}",
        lineages.height,
        perf_counter() - started,
    )
    return lineages


def write_wvdb_lineages(
    annotations_tsv: Path,
    lineages_csv: Path,
    *,
    reference_lineages_csv: Path | None = None,
) -> None:
    """Write sourmash-compatible WVDB lineages from the WVDB annotation TSV."""
    ensure_parent(lineages_csv)
    build_wvdb_lineages(
        scan_annotations(annotations_tsv),
        reference_lineages_csv=reference_lineages_csv,
    ).write_csv(lineages_csv)


def prepare_inputs(
    *,
    fasta: Path,
    annotations_tsv: Path,
    normalized_fasta: Path,
    lineages_csv: Path,
    reference_lineages_csv: Path | None = None,
) -> None:
    """Write normalized FASTA and sourmash lineage CSV outputs."""
    started = perf_counter()
    logger.info(
        "prepare_inputs.start fasta={} annotations_tsv={} reference_lineages_csv={} normalized_fasta={} lineages_csv={}",
        fasta,
        annotations_tsv,
        reference_lineages_csv or "none",
        normalized_fasta,
        lineages_csv,
    )
    sequence_count = normalize_fasta_ids(fasta, normalized_fasta)
    logger.info(
        "normalize_fasta_ids.done source={} destination={} sequences={}",
        fasta,
        normalized_fasta,
        sequence_count,
    )
    write_wvdb_lineages(
        annotations_tsv,
        lineages_csv,
        reference_lineages_csv=reference_lineages_csv,
    )
    logger.info(
        "prepare_inputs.done lineages_csv={} elapsed_seconds={:.2f}",
        lineages_csv,
        perf_counter() - started,
    )


def normalized_lineage_row(row: dict[str, str], *, path: Path) -> dict[str, str]:
    """Project a lineage row into the combined publishable schema."""
    ident = clean_cell(row.get("ident") or row.get("identifiers"))
    if not ident:
        msg = f"{path} contains a lineage row without ident/identifiers"
        raise ValueError(msg)

    taxpath = clean_cell(row.get(TAXPATH_COLUMN))
    if not taxpath:
        msg = f"{path} contains lineage {ident} without a taxpath"
        raise ValueError(msg)

    output = dict.fromkeys(COMBINED_LINEAGE_COLUMNS, "")
    output["ident"] = ident
    output[TAXPATH_COLUMN] = taxpath
    for rank in REFERENCE_TAXONOMY_RANKS:
        output[rank] = clean_cell(row.get(rank))
    return output


def combine_lineages(
    *,
    reference_lineages_csv: Path,
    wvdb_lineages_csv: Path,
    lineages_csv: Path,
) -> None:
    """Combine reference and WVDB lineages without stripping taxpaths."""
    started = perf_counter()
    logger.info(
        "combine_lineages.start reference_lineages_csv={} wvdb_lineages_csv={} lineages_csv={}",
        reference_lineages_csv,
        wvdb_lineages_csv,
        lineages_csv,
    )
    ensure_parent(lineages_csv)
    reference_rows = read_csv_dicts(reference_lineages_csv)
    wvdb_rows = read_csv_dicts(wvdb_lineages_csv)
    rows = [
        normalized_lineage_row(row, path=reference_lineages_csv)
        for row in reference_rows
    ]
    rows.extend(
        normalized_lineage_row(row, path=wvdb_lineages_csv) for row in wvdb_rows
    )

    with lineages_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COMBINED_LINEAGE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    logger.info(
        "combine_lineages.done reference_rows={} wvdb_rows={} output_rows={} elapsed_seconds={:.2f}",
        len(reference_rows),
        len(wvdb_rows),
        len(rows),
        perf_counter() - started,
    )


def find_missing_wvdb_lineages(
    manifest: pl.LazyFrame,
    lineages: pl.LazyFrame,
) -> pl.DataFrame:
    """Return WVDB manifest names that do not have matching lineages."""
    return (
        manifest.filter(pl.col("name").str.starts_with(WVDB_ID_PREFIX))
        .select("name")
        .join(lineages.select("ident"), left_on="name", right_on="ident", how="anti")
        .collect()
    )


def fail_on_missing_wvdb_lineages(missing: pl.DataFrame) -> pl.DataFrame:
    """Fail when the missing-lineage frame is not empty."""

    if missing.height:
        examples = ", ".join(missing.get_column("name").head(10).to_list())
        msg = f"{missing.height} WVDB signatures are missing lineages. Examples: {examples}"
        raise ValueError(msg)
    return missing


def check_manifest_coverage(manifest_csv: Path, lineages_csv: Path) -> None:
    """Verify that every WVDB signature in a sourmash manifest has a lineage."""
    manifest = scan_sourmash_manifest(manifest_csv)
    lineages = scan_csv(lineages_csv).pipe(
        normalize_lineage_identifier_column,
        path=lineages_csv,
    )

    find_missing_wvdb_lineages(manifest, lineages).pipe(fail_on_missing_wvdb_lineages)


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    """Read a CSV file into string-valued dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return [
            {key: clean_cell(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def read_manifest_dicts(path: Path, *, source: str) -> list[dict[str, str]]:
    """Read sourmash manifest rows and attach their curation source."""
    with path.open(newline="", encoding="utf-8") as handle:
        first_line = handle.readline()
        if not first_line.startswith("# SOURMASH-MANIFEST-VERSION:"):
            handle.seek(0)
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or ())
        required = {"name", *SIGNATURE_IDENTITY_COLUMNS}
        missing = required.difference(columns)
        if missing:
            msg = f"{path} is missing required column(s): {', '.join(sorted(missing))}"
            raise ValueError(msg)
        records = []
        for row in reader:
            record = {key: clean_cell(value) for key, value in row.items()}
            record["ident"] = record["name"].split(maxsplit=1)[0]
            record["source"] = source
            records.append(record)
        return records


def read_lineage_dicts(path: Path, *, source: str) -> list[dict[str, str]]:
    """Read lineage rows into one curation schema with source metadata."""
    rows = []
    for raw in read_csv_dicts(path):
        ident = clean_cell(raw.get("ident") or raw.get("identifiers"))
        if not ident:
            msg = f"{path} contains a lineage row without ident/identifiers"
            raise ValueError(msg)
        row = {rank: clean_cell(raw.get(rank)) for rank in CURATION_TAXONOMY_RANKS}
        row["ident"] = ident
        row[TAXPATH_COLUMN] = clean_cell(raw.get(TAXPATH_COLUMN))
        row["source"] = source
        rows.append(row)
    return rows


def lineage_by_identifier(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Index unique lineage records by identifier."""
    by_ident: dict[str, dict[str, str]] = {}
    for row in rows:
        ident = row["ident"]
        if ident in by_ident:
            msg = f"Duplicate lineage identifier: {ident}"
            raise ValueError(msg)
        by_ident[ident] = row
    return by_ident


def signature_identity(record: dict[str, str]) -> tuple[str, ...]:
    """Return the complete manifest-level identity of a sourmash sketch."""
    return tuple(record[column] for column in SIGNATURE_IDENTITY_COLUMNS)


def signature_identity_text(record: dict[str, str]) -> str:
    """Format a sketch identity for diagnostics."""
    return ";".join(
        f"{column}={record[column]}" for column in SIGNATURE_IDENTITY_COLUMNS
    )


def normalized_taxon_name(value: object) -> str:
    """Normalize WVDB/NCBI spelling differences for taxonomy comparison."""
    return " ".join(clean_cell(value).replace("_", " ").casefold().split())


def is_informative_taxon(value: object) -> bool:
    """Return whether a rank value provides more than unknown classification."""
    normalized = normalized_taxon_name(value)
    return normalized not in UNKNOWN_TAXA and not normalized.startswith("unclassified")


def lineage_prefix(row: dict[str, str], rank_index: int) -> tuple[str, ...]:
    """Return a normalized lineage prefix through the selected rank."""
    return tuple(
        normalized_taxon_name(row[rank])
        for rank in CURATION_TAXONOMY_RANKS[: rank_index + 1]
    )


def specificity(row: dict[str, str], rank_index: int) -> tuple[int, int]:
    """Score informative ranks through the taxon whose placement is contested."""
    informative = [
        index
        for index, rank in enumerate(CURATION_TAXONOMY_RANKS[: rank_index + 1])
        if is_informative_taxon(row[rank])
    ]
    return len(informative), max(informative, default=-1)


def lineage_text(row: dict[str, str]) -> str:
    """Format a lineage for actionable diagnostics."""
    return " > ".join(row[rank] or "<blank>" for rank in CURATION_TAXONOMY_RANKS)


def write_diagnostics(path: Path, rows: list[dict[str, str]]) -> None:
    """Write all automatic and unresolved curation decisions."""
    ensure_parent(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=DIAGNOSTIC_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def fail_on_duplicate_signature_identities(
    manifest_rows: list[dict[str, str]],
    lineages: dict[str, dict[str, str]],
    diagnostics: list[dict[str, str]],
) -> None:
    """Reject duplicate sketch identities instead of relying on input ordering."""
    grouped: defaultdict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for record in manifest_rows:
        grouped[signature_identity(record)].append(record)

    duplicate_groups = [records for records in grouped.values() if len(records) > 1]
    for records in duplicate_groups:
        for record in records:
            lineage = lineages.get(record["ident"], {})
            diagnostics.append(
                {
                    "conflict_type": "signature_identity",
                    "rank": "",
                    "normalized_name": "",
                    "policy": "error",
                    "decision": "ambiguous",
                    "source": record["source"],
                    "ident": record["ident"],
                    "signature_identity": signature_identity_text(record),
                    "specificity": "",
                    "lineage": lineage_text(lineage) if lineage else "<missing>",
                },
            )

    if duplicate_groups:
        records = duplicate_groups[0]
        msg = (
            "Duplicate sourmash signature identity: "
            f"{signature_identity_text(records[0])}; identifiers="
            f"{', '.join(record['ident'] for record in records)}; lineages="
            f"{' || '.join(lineage_text(lineages[record['ident']]) for record in records if record['ident'] in lineages)}"
        )
        raise ValueError(msg)


def preferred_source(policy: str) -> str | None:
    """Return the source selected by an explicit source-precedence policy."""
    if policy == "ncbi-wins":
        return "ncbi"
    if policy == "wvdb-wins":
        return "wvdb"
    return None


def select_canonical_lineage(
    rows: list[dict[str, str]],
    *,
    rank_index: int,
    policy: str,
) -> dict[str, str] | None:
    """Select one canonical lineage or return None for an unresolved conflict."""
    if policy == "error":
        return None

    candidates = rows
    source = preferred_source(policy)
    if source and any(row["source"] == source for row in rows):
        candidates = [row for row in rows if row["source"] == source]

    best_score = max(specificity(row, rank_index) for row in candidates)
    best = [row for row in candidates if specificity(row, rank_index) == best_score]
    prefixes = {lineage_prefix(row, rank_index) for row in best}
    if len(prefixes) != 1:
        return None
    return min(best, key=lambda row: row["ident"])


def taxonomy_conflicts_at_rank(
    rows: list[dict[str, str]],
    rank_index: int,
) -> list[tuple[str, list[dict[str, str]]]]:
    """Return taxon-name groups assigned to incompatible parents at one rank."""
    rank = CURATION_TAXONOMY_RANKS[rank_index]
    grouped: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        name = normalized_taxon_name(row[rank])
        if name:
            grouped[name].append(row)
    return [
        (name, candidates)
        for name, candidates in grouped.items()
        if len({lineage_prefix(row, rank_index) for row in candidates}) > 1
    ]


def append_taxonomy_diagnostics(
    diagnostics: list[dict[str, str]],
    candidates: list[dict[str, str]],
    *,
    conflict: tuple[int, str, str],
    winner: dict[str, str] | None,
) -> None:
    """Record every candidate participating in one taxonomy-name conflict."""
    rank_index, name, policy = conflict
    rank = CURATION_TAXONOMY_RANKS[rank_index]
    winner_prefix = lineage_prefix(winner, rank_index) if winner else None
    for row in candidates:
        decision = (
            "ambiguous"
            if winner is None
            else "winner"
            if lineage_prefix(row, rank_index) == winner_prefix
            else "rewritten"
        )
        diagnostics.append(
            {
                "conflict_type": "taxonomy_name",
                "rank": rank,
                "normalized_name": name,
                "policy": policy,
                "decision": decision,
                "source": row["source"],
                "ident": row["ident"],
                "signature_identity": "",
                "specificity": ":".join(
                    str(value) for value in specificity(row, rank_index)
                ),
                "lineage": lineage_text(row),
            },
        )


def rewrite_taxpath_prefix(
    row: dict[str, str],
    winner: dict[str, str],
    rank_index: int,
) -> None:
    """Keep taxids aligned when curation rewrites a lineage prefix."""
    taxids = split_taxpath(row.get(TAXPATH_COLUMN, ""))
    winner_taxids = split_taxpath(winner.get(TAXPATH_COLUMN, ""))
    if not taxids or not winner_taxids:
        return

    width = len(CURATION_TAXONOMY_RANKS)
    taxids.extend([""] * (width - len(taxids)))
    winner_taxids.extend([""] * (width - len(winner_taxids)))
    for index in range(rank_index + 1):
        taxids[index] = winner_taxids[index]
    while taxids and not taxids[-1]:
        taxids.pop()
    row[TAXPATH_COLUMN] = "|".join(taxids)


def curate_taxonomy_names(
    rows: list[dict[str, str]],
    *,
    policy: str,
    diagnostics: list[dict[str, str]],
) -> tuple[int, list[str]]:
    """Canonicalize conflicting taxon placements from broad to specific ranks."""
    resolved = 0
    unresolved: list[str] = []
    for rank_index, rank in enumerate(CURATION_TAXONOMY_RANKS[1:], start=1):
        for name, candidates in taxonomy_conflicts_at_rank(rows, rank_index):
            winner = select_canonical_lineage(
                candidates,
                rank_index=rank_index,
                policy=policy,
            )
            if winner is None:
                unresolved.append(f"{rank}={name}")

            append_taxonomy_diagnostics(
                diagnostics,
                candidates,
                conflict=(rank_index, name, policy),
                winner=winner,
            )

            if winner is None:
                continue

            resolved += 1
            for row in candidates:
                for ancestor_rank in CURATION_TAXONOMY_RANKS[: rank_index + 1]:
                    row[ancestor_rank] = winner[ancestor_rank]
                rewrite_taxpath_prefix(row, winner, rank_index)

    return resolved, unresolved


def find_cross_rank_taxon_names(
    rows: list[dict[str, str]],
    *,
    policy: str,
    diagnostics: list[dict[str, str]],
) -> list[str]:
    """Record names reused at multiple ranks, which Taxburst cannot represent."""
    occurrences: defaultdict[str, list[tuple[str, dict[str, str]]]] = defaultdict(list)
    for row in rows:
        for rank in CURATION_TAXONOMY_RANKS:
            name = normalized_taxon_name(row[rank])
            if name:
                occurrences[name].append((rank, row))

    conflicts = []
    for name, candidates in occurrences.items():
        ranks = {rank for rank, _row in candidates}
        if len(ranks) <= 1:
            continue
        conflicts.append(f"{name} ({', '.join(sorted(ranks))})")
        for rank, row in candidates:
            diagnostics.append(
                {
                    "conflict_type": "taxonomy_rank",
                    "rank": rank,
                    "normalized_name": name,
                    "policy": policy,
                    "decision": "ambiguous",
                    "source": row["source"],
                    "ident": row["ident"],
                    "signature_identity": "",
                    "specificity": "",
                    "lineage": lineage_text(row),
                },
            )
    return conflicts


def write_curated_lineages(path: Path, rows: list[dict[str, str]]) -> None:
    """Write only the lineage columns consumed by sourmash taxonomy commands."""
    ensure_parent(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CURATED_LINEAGE_COLUMNS,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {column: row[column] for column in CURATED_LINEAGE_COLUMNS} for row in rows
        )


def curate_lineages(
    paths: CurationPaths,
    *,
    policy: str,
) -> None:
    """Validate sketch identities and canonicalize conflicting taxonomy names."""
    manifest_rows = [
        *read_manifest_dicts(paths.reference_manifest_csv, source="ncbi"),
        *read_manifest_dicts(paths.wvdb_manifest_csv, source="wvdb"),
    ]
    lineage_rows = [
        *read_lineage_dicts(paths.reference_lineages_csv, source="ncbi"),
        *read_lineage_dicts(paths.wvdb_lineages_csv, source="wvdb"),
    ]
    lineages = lineage_by_identifier(lineage_rows)
    missing = [
        record["ident"] for record in manifest_rows if record["ident"] not in lineages
    ]
    if missing:
        msg = f"{len(missing)} signatures are missing lineages. Examples: {', '.join(missing[:10])}"
        raise ValueError(msg)

    diagnostics: list[dict[str, str]] = []
    try:
        fail_on_duplicate_signature_identities(manifest_rows, lineages, diagnostics)
    except ValueError:
        write_diagnostics(paths.diagnostics_tsv, diagnostics)
        raise

    retained_rows = [lineages[record["ident"]] for record in manifest_rows]
    resolved, unresolved = curate_taxonomy_names(
        retained_rows,
        policy=policy,
        diagnostics=diagnostics,
    )
    cross_rank_conflicts = find_cross_rank_taxon_names(
        retained_rows,
        policy=policy,
        diagnostics=diagnostics,
    )
    write_diagnostics(paths.diagnostics_tsv, diagnostics)
    if resolved:
        warnings.warn(
            f"resolved {resolved} taxonomy conflict(s) using policy {policy}; "
            f"evidence={paths.diagnostics_tsv}",
            RuntimeWarning,
            stacklevel=2,
        )
    if unresolved or cross_rank_conflicts:
        examples = [*unresolved, *cross_rank_conflicts]
        msg = (
            f"{len(examples)} ambiguous taxonomy conflict(s) under policy {policy}. "
            f"Examples: {', '.join(examples[:10])}. Evidence: {paths.diagnostics_tsv}"
        )
        raise ValueError(msg)

    write_curated_lineages(paths.lineages_csv, retained_rows)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Prepare WVDB inputs for a combined NCBI Virus + WVDB sourmash reference.",
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    prepare = subcommands.add_parser(
        "prepare-inputs",
        help="Normalize WVDB FASTA IDs and write WVDB sourmash lineages.",
    )
    prepare.add_argument("--fasta", type=Path, required=True)
    prepare.add_argument("--annotations-tsv", type=Path, required=True)
    prepare.add_argument(
        "--reference-lineages-csv",
        type=Path,
        help="Reference sourmash lineages CSV used to reuse exact taxpaths.",
    )
    prepare.add_argument("--normalized-fasta", type=Path, required=True)
    prepare.add_argument("--lineages-csv", type=Path, required=True)

    combine = subcommands.add_parser(
        "combine-lineages",
        help="Combine reference and WVDB lineages while preserving taxpaths.",
    )
    combine.add_argument("--reference-lineages-csv", type=Path, required=True)
    combine.add_argument("--wvdb-lineages-csv", type=Path, required=True)
    combine.add_argument("--lineages-csv", type=Path, required=True)

    check = subcommands.add_parser(
        "check-manifest-coverage",
        help="Fail if any WVDB signatures in a sourmash manifest lack lineages.",
    )
    check.add_argument("--manifest-csv", type=Path, required=True)
    check.add_argument("--lineages-csv", type=Path, required=True)

    curate = subcommands.add_parser(
        "curate-lineages",
        help="Validate signature identity and resolve conflicting taxonomy names.",
    )
    curate.add_argument("--reference-manifest-csv", type=Path, required=True)
    curate.add_argument("--wvdb-manifest-csv", type=Path, required=True)
    curate.add_argument("--reference-lineages-csv", type=Path, required=True)
    curate.add_argument("--wvdb-lineages-csv", type=Path, required=True)
    curate.add_argument(
        "--taxonomy-conflict-policy",
        choices=TAXONOMY_CONFLICT_POLICIES,
        required=True,
    )
    curate.add_argument("--lineages-csv", type=Path, required=True)
    curate.add_argument("--diagnostics-tsv", type=Path, required=True)

    return parser


def main() -> None:
    """Run the command-line interface."""
    configure_logging()
    parser = build_parser()
    args = parser.parse_args()

    try:
        if args.command == "prepare-inputs":
            prepare_inputs(
                fasta=args.fasta,
                annotations_tsv=args.annotations_tsv,
                normalized_fasta=args.normalized_fasta,
                lineages_csv=args.lineages_csv,
                reference_lineages_csv=args.reference_lineages_csv,
            )
        elif args.command == "combine-lineages":
            combine_lineages(
                reference_lineages_csv=args.reference_lineages_csv,
                wvdb_lineages_csv=args.wvdb_lineages_csv,
                lineages_csv=args.lineages_csv,
            )
        elif args.command == "check-manifest-coverage":
            check_manifest_coverage(args.manifest_csv, args.lineages_csv)
        elif args.command == "curate-lineages":
            curate_lineages(
                CurationPaths(
                    reference_manifest_csv=args.reference_manifest_csv,
                    wvdb_manifest_csv=args.wvdb_manifest_csv,
                    reference_lineages_csv=args.reference_lineages_csv,
                    wvdb_lineages_csv=args.wvdb_lineages_csv,
                    lineages_csv=args.lineages_csv,
                    diagnostics_tsv=args.diagnostics_tsv,
                ),
                policy=args.taxonomy_conflict_policy,
            )
    except ValueError as exc:
        parser.exit(1, f"error: {exc}\n")


if __name__ == "__main__":
    main()
