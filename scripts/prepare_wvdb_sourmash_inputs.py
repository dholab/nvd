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
UNKNOWN_TAXA = frozenset({"", "unknown", "na", "n/a", "none"})
DIAGNOSTIC_COLUMNS = (
    "conflict_type",
    "rank",
    "normalized_name",
    "resolution",
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


def expression_is_informative(value: pl.Expr) -> pl.Expr:
    """Return whether a string expression contains a classified value."""
    normalized = value.cast(pl.Utf8).str.strip_chars().str.to_lowercase()
    return (
        (normalized.str.len_chars() > 0)
        & ~normalized.is_in(PLACEHOLDER_SENTINELS)
        & ~normalized.str.starts_with("unclassified")
    )


def prefer_informative(primary: pl.Expr, fallback: pl.Expr) -> pl.Expr:
    """Prefer classified primary evidence without discarding explicit sentinels."""
    return (
        pl.when(expression_is_informative(primary))
        .then(primary)
        .when(expression_is_informative(fallback))
        .then(fallback)
        .when(primary.cast(pl.Utf8).str.strip_chars().str.len_chars() > 0)
        .then(primary)
        .otherwise(fallback)
    )


def lineage_source_score(values: dict[str, pl.Expr]) -> tuple[pl.Expr, pl.Expr]:
    """Return deepest classified rank and classified-rank count for one source."""
    informative = {
        rank: expression_is_informative(value) for rank, value in values.items()
    }
    depth = pl.max_horizontal(
        pl.when(informative[rank]).then(pl.lit(index)).otherwise(pl.lit(0))
        for index, rank in enumerate(("phylum", "class", "order", "family"), start=1)
    )
    count = pl.sum_horizontal(
        informative[rank].cast(pl.UInt8)
        for rank in ("phylum", "class", "order", "family")
    )
    return depth, count


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


def wvdb_local_taxid(
    prefix: tuple[tuple[str, str], ...],
    *,
    collision_attempt: int = 0,
) -> str:
    """Return a stable numeric WVDB-local taxid for a lineage prefix."""
    if any("|" in name or "=" in name for _rank, name in prefix):
        encoded = "".join(
            f"{len(rank)}:{rank}{len(name)}:{name}" for rank, name in prefix
        )
        canonical = f"v2:{encoded}"
    else:
        canonical = "|".join(f"{rank}={name}" for rank, name in prefix)
    if collision_attempt:
        canonical = f"{canonical}\x00collision_attempt={collision_attempt}"
    digest = hashlib.sha1(canonical.encode("utf-8"), usedforsecurity=False)
    offset = int(digest.hexdigest()[:12], 16) % WVDB_LOCAL_TAXID_MODULUS
    return str(WVDB_LOCAL_TAXID_BASE + offset)


def is_wvdb_local_taxid(taxid: str) -> bool:
    """Return whether a numeric taxid belongs to NVD's synthetic range."""
    if not taxid.isdecimal():
        return False
    value = int(taxid)
    return (
        WVDB_LOCAL_TAXID_BASE
        <= value
        < (WVDB_LOCAL_TAXID_BASE + WVDB_LOCAL_TAXID_MODULUS)
    )


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
    rdrp = {
        rank: first_nonempty(available_columns, f"{rank.title()}_RdRp")
        for rank in ("phylum", "class", "order", "family")
    }
    gnd = {
        rank: first_nonempty(available_columns, f"{rank.title()}_gNd")
        for rank in ("phylum", "class", "order", "family")
    }
    rdrp_depth, rdrp_count = lineage_source_score(rdrp)
    gnd_depth, gnd_count = lineage_source_score(gnd)
    use_rdrp = (rdrp_depth > gnd_depth) | (
        (rdrp_depth == gnd_depth) & (rdrp_count >= gnd_count)
    )
    selected = {
        rank: pl.when(use_rdrp).then(rdrp[rank]).otherwise(gnd[rank])
        for rank in ("phylum", "class", "order", "family")
    }
    order = prefer_informative(
        first_nonempty(available_columns, "Order_consensus"),
        selected["order"],
    )
    family = prefer_informative(
        first_nonempty(available_columns, "Family_consensus"),
        prefer_informative(
            selected["family"],
            first_nonempty(available_columns, "Family_ICTV"),
        ),
    )
    return annotations.select(
        (
            pl.lit(WVDB_ID_PREFIX) + pl.col("contig").cast(pl.Utf8).str.strip_chars()
        ).alias("ident"),
        pl.lit("Viruses").alias("superkingdom"),
        selected["phylum"].alias("phylum"),
        selected["class"].alias("class"),
        order.alias("order"),
        family.alias("family"),
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


@dataclass(frozen=True)
class TaxpathValidation:
    """Collected positional and global taxpath invariant violations."""

    short: list[str]
    long: list[str]
    missing_taxids: list[str]
    orphan_taxids: list[str]
    non_numeric: list[str]
    taxids_by_prefix: dict[tuple[tuple[str, str], ...], set[str]]
    prefixes_by_taxid: dict[str, set[tuple[tuple[str, str], ...]]]


def require_publishable_lineage_schema(lineages_csv: Path) -> None:
    """Require every identifier, rank, and taxpath column in the artifact."""
    with lineages_csv.open(newline="", encoding="utf-8") as handle:
        columns = set(csv.DictReader(handle).fieldnames or ())
    missing_columns = {"ident", TAXPATH_COLUMN, *REFERENCE_TAXONOMY_RANKS}.difference(
        columns,
    )
    if missing_columns:
        msg = f"missing required lineage columns: {', '.join(sorted(missing_columns))}"
        raise ValueError(msg)


def identifier_validation_errors(
    manifest_rows: list[dict[str, str]],
    rows: list[dict[str, str]],
) -> list[str]:
    """Return manifest and lineage identifier cardinality errors."""
    manifest_identifier_counts = Counter(row["ident"] for row in manifest_rows)
    manifest_identifiers = set(manifest_identifier_counts)
    duplicate_manifest_identifiers = sorted(
        ident for ident, count in manifest_identifier_counts.items() if count > 1
    )
    identifier_counts = Counter(
        clean_cell(row.get("ident") or row.get("identifiers")) for row in rows
    )
    blank_identifiers = identifier_counts.get("", 0)
    lineage_identifiers = set(identifier_counts).difference({""})
    manifest_identifiers_without_lineages = sorted(
        manifest_identifiers - lineage_identifiers,
    )
    lineage_identifiers_absent_from_manifest = sorted(
        lineage_identifiers - manifest_identifiers,
    )
    duplicate_identifiers = sorted(
        ident for ident, count in identifier_counts.items() if count > 1
    )
    errors: list[str] = []
    if duplicate_identifiers:
        errors.append(
            f"duplicate identifiers: {len(duplicate_identifiers)}; examples: "
            f"{', '.join(duplicate_identifiers[:10])}",
        )
    if blank_identifiers:
        errors.append(f"blank identifiers: {blank_identifiers}")
    if duplicate_manifest_identifiers:
        errors.append(
            "duplicate manifest identifiers: "
            f"{len(duplicate_manifest_identifiers)}; examples: "
            f"{', '.join(duplicate_manifest_identifiers[:10])}",
        )
    if manifest_identifiers_without_lineages:
        errors.append(
            "manifest identifiers without lineages: "
            f"{len(manifest_identifiers_without_lineages)}; examples: "
            f"{', '.join(manifest_identifiers_without_lineages[:10])}",
        )
    if lineage_identifiers_absent_from_manifest:
        errors.append(
            "lineage identifiers absent from manifest: "
            f"{len(lineage_identifiers_absent_from_manifest)}; examples: "
            f"{', '.join(lineage_identifiers_absent_from_manifest[:10])}",
        )
    return errors


def inspect_taxpaths(rows: list[dict[str, str]]) -> TaxpathValidation:
    """Collect positional and global taxonomy-node mapping evidence."""
    width = len(REFERENCE_TAXONOMY_RANKS)
    short_taxpaths: list[str] = []
    long_taxpaths: list[str] = []
    populated_ranks_without_taxids: list[str] = []
    taxids_without_populated_ranks: list[str] = []
    non_numeric_taxids: list[str] = []
    taxids_by_prefix: defaultdict[tuple[tuple[str, str], ...], set[str]] = defaultdict(
        set,
    )
    prefixes_by_taxid: defaultdict[str, set[tuple[tuple[str, str], ...]]] = defaultdict(
        set,
    )
    for row in rows:
        ident = clean_cell(row.get("ident") or row.get("identifiers"))
        taxpath = clean_cell(row.get(TAXPATH_COLUMN))
        taxids = taxpath.split("|") if taxpath else []
        if len(taxids) < width:
            short_taxpaths.append(ident)
        elif len(taxids) > width:
            long_taxpaths.append(ident)
        positional_taxids = [*taxids[:width], *([""] * max(0, width - len(taxids)))]
        prefix: list[tuple[str, str]] = []
        for index, rank in enumerate(REFERENCE_TAXONOMY_RANKS):
            name = clean_cell(row.get(rank))
            taxid = positional_taxids[index]
            if name and not positional_taxids[index]:
                populated_ranks_without_taxids.append(f"{ident}:{rank}")
            elif not name and positional_taxids[index]:
                taxids_without_populated_ranks.append(f"{ident}:{rank}")
            if taxid and not taxid.isdecimal():
                non_numeric_taxids.append(
                    f"{ident}:{rank}={taxid}",
                )
            if name:
                prefix.append((rank, name))
                if taxid:
                    prefix_key = tuple(prefix)
                    taxids_by_prefix[prefix_key].add(taxid)
                    prefixes_by_taxid[taxid].add(prefix_key)
    return TaxpathValidation(
        short=short_taxpaths,
        long=long_taxpaths,
        missing_taxids=populated_ranks_without_taxids,
        orphan_taxids=taxids_without_populated_ranks,
        non_numeric=non_numeric_taxids,
        taxids_by_prefix=dict(taxids_by_prefix),
        prefixes_by_taxid=dict(prefixes_by_taxid),
    )


def taxpath_validation_errors(validation: TaxpathValidation) -> list[str]:
    """Format every failed positional or global taxpath invariant."""
    prefixes_with_conflicting_taxids = {
        prefix: taxids
        for prefix, taxids in validation.taxids_by_prefix.items()
        if len(taxids) > 1
    }
    taxids_with_conflicting_prefixes = {
        taxid: prefixes
        for taxid, prefixes in validation.prefixes_by_taxid.items()
        if len(prefixes) > 1
    }
    errors: list[str] = []
    if validation.short:
        errors.append(
            f"short taxpaths: {len(validation.short)}; examples: "
            f"{', '.join(validation.short[:10])}",
        )
    if validation.long:
        errors.append(
            f"long taxpaths: {len(validation.long)}; examples: "
            f"{', '.join(validation.long[:10])}",
        )
    if validation.missing_taxids:
        errors.append(
            "populated ranks with missing taxids: "
            f"{len(validation.missing_taxids)}; examples: "
            f"{', '.join(validation.missing_taxids[:10])}",
        )
    if validation.orphan_taxids:
        errors.append(
            "taxids without populated ranks: "
            f"{len(validation.orphan_taxids)}; examples: "
            f"{', '.join(validation.orphan_taxids[:10])}",
        )
    if validation.non_numeric:
        errors.append(
            f"non-numeric taxids: {len(validation.non_numeric)}; examples: "
            f"{', '.join(validation.non_numeric[:10])}",
        )
    if prefixes_with_conflicting_taxids:
        errors.append(
            "lineage prefixes with conflicting taxids: "
            f"{len(prefixes_with_conflicting_taxids)}",
        )
    if taxids_with_conflicting_prefixes:
        errors.append(
            "taxids with conflicting lineage prefixes: "
            f"{len(taxids_with_conflicting_prefixes)}",
        )
    return errors


def semantic_taxonomy_errors(rows: list[dict[str, str]]) -> list[str]:
    """Return incompatible same-rank placements and cross-rank name reuse."""
    same_rank_conflicts = 0
    ranks_by_name: defaultdict[str, set[str]] = defaultdict(set)
    for rank_index, rank in enumerate(REFERENCE_TAXONOMY_RANKS):
        prefixes_by_name: defaultdict[str, set[tuple[str, ...]]] = defaultdict(set)
        for row in rows:
            name = normalized_taxon_name(row[rank])
            if not name:
                continue
            ranks_by_name[name].add(rank)
            prefixes_by_name[name].add(comparison_lineage_prefix(row, rank_index))
        same_rank_conflicts += sum(
            len(prefixes) > 1 for prefixes in prefixes_by_name.values()
        )
    cross_rank_conflicts = sum(len(ranks) > 1 for ranks in ranks_by_name.values())
    errors = []
    if same_rank_conflicts:
        errors.append(f"same-rank taxonomy conflicts: {same_rank_conflicts}")
    if cross_rank_conflicts:
        errors.append(f"cross-rank taxonomy conflicts: {cross_rank_conflicts}")
    return errors


def validate_lineages(manifest_csv: Path, lineages_csv: Path) -> None:
    """Validate positional BioBoxes invariants for a production lineage artifact."""
    manifest_rows = read_manifest_dicts(manifest_csv, source="combined")
    require_publishable_lineage_schema(lineages_csv)
    rows = read_csv_dicts(lineages_csv)
    errors = [
        *identifier_validation_errors(manifest_rows, rows),
        *taxpath_validation_errors(inspect_taxpaths(rows)),
        *semantic_taxonomy_errors(rows),
    ]
    if errors:
        raise ValueError("\n".join(["Lineage validation failed:", *errors]))


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
        row = {rank: clean_cell(raw.get(rank)) for rank in REFERENCE_TAXONOMY_RANKS}
        for rank in REFERENCE_TAXONOMY_RANKS:
            row[f"source_{rank}"] = row[rank]
        row["ident"] = ident
        row["leaf_taxid"] = clean_cell(raw.get("taxid"))
        row[TAXPATH_COLUMN] = clean_cell(raw.get(TAXPATH_COLUMN))
        row["source_taxpath"] = row[TAXPATH_COLUMN]
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
        for rank in REFERENCE_TAXONOMY_RANKS[: rank_index + 1]
    )


def display_lineage_prefix(row: dict[str, str], rank_index: int) -> tuple[str, ...]:
    """Return the exact displayed lineage prefix through the selected rank."""
    return tuple(
        clean_cell(row[rank]) for rank in REFERENCE_TAXONOMY_RANKS[: rank_index + 1]
    )


def comparison_taxon_name(value: object) -> str:
    """Return a semantic comparison key for classified and placeholder taxa."""
    if not is_informative_taxon(value):
        return "<uninformative>"
    return normalized_taxon_name(value)


def comparison_lineage_prefix(
    row: dict[str, str],
    rank_index: int,
) -> tuple[str, ...]:
    """Return a lineage prefix with equivalent placeholders canonicalized."""
    return tuple(
        comparison_taxon_name(row[rank])
        for rank in REFERENCE_TAXONOMY_RANKS[: rank_index + 1]
    )


def first_divergent_rank(prefixes: set[tuple[str, ...]]) -> int:
    """Return the first rank at which otherwise comparable prefixes differ."""
    width = len(next(iter(prefixes)))
    for rank_index in range(width):
        if len({prefix[rank_index] for prefix in prefixes}) > 1:
            return rank_index
    msg = "Cannot find a divergent rank in equivalent lineage prefixes"
    raise ValueError(msg)


def specificity(row: dict[str, str], rank_index: int) -> tuple[int, int]:
    """Score informative ranks through the taxon whose placement is contested."""
    informative = [
        index
        for index, rank in enumerate(REFERENCE_TAXONOMY_RANKS[: rank_index + 1])
        if is_informative_taxon(row[rank])
    ]
    return len(informative), max(informative, default=-1)


def lineage_text(row: dict[str, str]) -> str:
    """Format a lineage for actionable diagnostics."""
    return " > ".join(row[rank] or "<blank>" for rank in REFERENCE_TAXONOMY_RANKS)


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


def diagnose_duplicate_signature_identities(
    manifest_rows: list[dict[str, str]],
    lineages: dict[str, dict[str, str]],
    diagnostics: list[dict[str, str]],
) -> str | None:
    """Record duplicate sketch identities and return an error message when found."""
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
                    "resolution": "rejected-duplicate-signature",
                    "decision": "ambiguous",
                    "source": record["source"],
                    "ident": record["ident"],
                    "signature_identity": signature_identity_text(record),
                    "specificity": "",
                    "lineage": lineage_text(lineage) if lineage else "<missing>",
                },
            )

    if not duplicate_groups:
        return None

    records = duplicate_groups[0]
    return (
        "Duplicate sourmash signature identity: "
        f"{signature_identity_text(records[0])}; identifiers="
        f"{', '.join(record['ident'] for record in records)}; lineages="
        f"{' || '.join(lineage_text(lineages[record['ident']]) for record in records if record['ident'] in lineages)}"
    )


def select_canonical_lineage(
    rows: list[dict[str, str]],
    *,
    rank_index: int,
) -> tuple[dict[str, str] | None, str]:
    """Select a defensible canonical lineage or request contextualization."""
    best_score = max(specificity(row, rank_index) for row in rows)
    best = [row for row in rows if specificity(row, rank_index) == best_score]
    prefixes = {comparison_lineage_prefix(row, rank_index) for row in best}
    if len(prefixes) == 1:
        resolution = (
            "selected-most-specific"
            if len(best) < len(rows)
            else "canonicalized-placeholder"
        )
        return (
            min(best, key=lambda row: (row["source"] != "ncbi", row["ident"])),
            resolution,
        )

    divergent_rank = first_divergent_rank(prefixes)
    informative_ncbi_prefixes = {
        comparison_lineage_prefix(row, rank_index)
        for row in best
        if row["source"] == "ncbi"
        and is_informative_taxon(row[REFERENCE_TAXONOMY_RANKS[divergent_rank]])
    }
    if len(informative_ncbi_prefixes) != 1:
        return None, "contextualized-ambiguity"
    ncbi_prefix = next(iter(informative_ncbi_prefixes))
    return (
        min(
            (
                row
                for row in best
                if comparison_lineage_prefix(row, rank_index) == ncbi_prefix
            ),
            key=lambda row: (row["source"] != "ncbi", row["ident"]),
        ),
        "selected-informative-ncbi",
    )


def taxonomy_conflicts_at_rank(
    rows: list[dict[str, str]],
    rank_index: int,
) -> list[tuple[str, list[dict[str, str]]]]:
    """Return taxon-name groups assigned to incompatible parents at one rank."""
    rank = REFERENCE_TAXONOMY_RANKS[rank_index]
    grouped: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        name = normalized_taxon_name(row[rank])
        if name:
            grouped[name].append(row)
    return [
        (name, candidates)
        for name, candidates in grouped.items()
        if len({display_lineage_prefix(row, rank_index) for row in candidates}) > 1
    ]


def append_taxonomy_diagnostics(
    diagnostics: list[dict[str, str]],
    candidates: list[dict[str, str]],
    *,
    conflict: tuple[int, str, str],
    winner: dict[str, str],
) -> None:
    """Record every candidate participating in one taxonomy-name conflict."""
    rank_index, name, resolution = conflict
    rank = REFERENCE_TAXONOMY_RANKS[rank_index]
    winner_prefix = display_lineage_prefix(winner, rank_index)
    for row in candidates:
        decision = (
            "winner"
            if display_lineage_prefix(row, rank_index) == winner_prefix
            else "rewritten"
        )
        diagnostics.append(
            {
                "conflict_type": "taxonomy_name",
                "rank": rank,
                "normalized_name": name,
                "resolution": resolution,
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


def contextual_taxon_base(candidates: list[dict[str, str]], rank: str) -> str:
    """Choose a stable display name for one ambiguous taxon-name group."""
    ncbi_names = [row[rank] for row in candidates if row["source"] == "ncbi"]
    names = ncbi_names or [row[rank] for row in candidates]
    return min(names, key=lambda name: (name.count("_"), name.casefold(), name))


def informative_parent_contexts(
    row: dict[str, str],
    rank_index: int,
) -> list[str]:
    """Return increasingly specific informative parent contexts."""
    parents = [
        row[rank]
        for rank in REFERENCE_TAXONOMY_RANKS[:rank_index]
        if is_informative_taxon(row[rank])
    ]
    return [" > ".join(parents[-depth:]) for depth in range(1, len(parents) + 1)]


def contextualize_taxonomy_conflict(
    candidates: list[dict[str, str]],
    *,
    rank_index: int,
    reserved_names: set[str],
) -> None:
    """Give each incompatible parent placement a stable contextual taxon name."""
    rank = REFERENCE_TAXONOMY_RANKS[rank_index]
    base = contextual_taxon_base(candidates, rank)
    rows_by_parent: defaultdict[tuple[str, ...], list[dict[str, str]]] = defaultdict(
        list,
    )
    for row in candidates:
        rows_by_parent[lineage_prefix(row, rank_index - 1)].append(row)

    contexts = {
        parent: informative_parent_contexts(rows[0], rank_index)
        for parent, rows in rows_by_parent.items()
    }
    context_indexes = dict.fromkeys(rows_by_parent, 0)
    while True:
        labels = {
            parent: f"{base} [under {contexts[parent][context_indexes[parent]]}]"
            for parent in rows_by_parent
        }
        normalized_counts = Counter(
            normalized_taxon_name(label) for label in labels.values()
        )
        colliding = {
            parent
            for parent, label in labels.items()
            if normalized_counts[normalized_taxon_name(label)] > 1
            or normalized_taxon_name(label) in reserved_names
        }
        if not colliding:
            break

        advanced = False
        for parent in colliding:
            if context_indexes[parent] < len(contexts[parent]) - 1:
                context_indexes[parent] += 1
                advanced = True
        if not advanced:
            labels = {
                parent: f"{label} #{short_hash(parent)}"
                if parent in colliding
                else label
                for parent, label in labels.items()
            }
            break

    for parent, rows in rows_by_parent.items():
        label = labels[parent]
        reserved_names.add(normalized_taxon_name(label))
        for row in rows:
            row[rank] = label


def append_contextualization_diagnostics(
    diagnostics: list[dict[str, str]],
    candidates: list[dict[str, str]],
    *,
    rank_index: int,
    name: str,
) -> None:
    """Record loss-preserving renames for an incompatible taxon placement."""
    rank = REFERENCE_TAXONOMY_RANKS[rank_index]
    diagnostics.extend(
        [
            {
                "conflict_type": "taxonomy_name",
                "rank": rank,
                "normalized_name": name,
                "resolution": "contextualized-ambiguity",
                "decision": "contextualized",
                "source": row["source"],
                "ident": row["ident"],
                "signature_identity": "",
                "specificity": ":".join(
                    str(value) for value in specificity(row, rank_index)
                ),
                "lineage": lineage_text(row),
            }
            for row in candidates
        ],
    )


def source_taxid_candidates(
    rows: list[dict[str, str]],
) -> tuple[
    set[tuple[tuple[str, str], ...]],
    dict[tuple[tuple[str, str], ...], list[tuple[str, str]]],
]:
    """Collect finalized lineage prefixes and their source taxid evidence."""
    all_prefixes: set[tuple[tuple[str, str], ...]] = set()
    candidates: defaultdict[tuple[tuple[str, str], ...], list[tuple[str, str]]] = (
        defaultdict(list)
    )
    width = len(REFERENCE_TAXONOMY_RANKS)
    for row in rows:
        source_taxids = split_taxpath(row.get("source_taxpath", ""))
        if len(source_taxids) > width:
            msg = (
                f"{row['ident']}: source taxpath has {len(source_taxids)} positions; "
                f"expected at most {width}"
            )
            raise ValueError(msg)
        source_taxids.extend([""] * (width - len(source_taxids)))
        leaf_taxid = clean_cell(row.get("leaf_taxid"))
        if leaf_taxid and not leaf_taxid.isdecimal():
            msg = f"non-numeric taxid: {leaf_taxid}"
            raise ValueError(msg)
        deepest_named_index = max(
            (
                index
                for index, rank in enumerate(REFERENCE_TAXONOMY_RANKS)
                if clean_cell(row[rank])
            ),
            default=-1,
        )
        prefix: list[tuple[str, str]] = []
        for index, rank in enumerate(REFERENCE_TAXONOMY_RANKS):
            name = clean_cell(row[rank])
            if not name:
                continue
            prefix.append((rank, name))
            all_prefixes.add(tuple(prefix))
            source_name = clean_cell(row[f"source_{rank}"])
            same_taxon = bool(source_name) and (
                comparison_taxon_name(source_name) == comparison_taxon_name(name)
            )
            synthetic_wvdb_evidence = row["source"] == "wvdb" and is_wvdb_local_taxid(
                source_taxids[index],
            )
            if source_taxids[index] and same_taxon and not synthetic_wvdb_evidence:
                if not source_taxids[index].isdecimal():
                    msg = f"non-numeric taxid: {source_taxids[index]}"
                    raise ValueError(msg)
                candidates[tuple(prefix)].append((row["source"], source_taxids[index]))
            if (
                row["source"] == "ncbi"
                and index == deepest_named_index
                and leaf_taxid
                and leaf_taxid not in source_taxids
                and same_taxon
            ):
                candidates[tuple(prefix)].append((row["source"], leaf_taxid))
    return all_prefixes, dict(candidates)


def select_source_taxid(
    prefix: tuple[tuple[str, str], ...],
    values: list[tuple[str, str]],
    diagnostics: list[dict[str, str]],
) -> str | None:
    """Select unambiguous native evidence for one finalized lineage prefix."""
    all_taxids = {taxid for _source, taxid in values}
    ncbi_taxids = {taxid for source, taxid in values if source == "ncbi"}
    lineage = " > ".join(f"{rank}={name}" for rank, name in prefix)
    if len(ncbi_taxids) > 1:
        diagnostics.append(
            {
                "conflict_type": "taxid",
                "rank": prefix[-1][0],
                "normalized_name": normalized_taxon_name(prefix[-1][1]),
                "resolution": "synthesized-conflicting-native-taxids",
                "decision": "synthesized",
                "source": "ncbi",
                "ident": "",
                "signature_identity": "",
                "specificity": "",
                "lineage": f"{lineage}; taxids={','.join(sorted(ncbi_taxids))}",
            },
        )
        return None
    if ncbi_taxids:
        return next(iter(ncbi_taxids))
    if len(all_taxids) > 1:
        msg = (
            f"conflicting taxids for lineage prefix {lineage}: "
            f"{', '.join(sorted(all_taxids))}"
        )
        raise ValueError(msg)
    return next(iter(all_taxids), None)


def fail_on_conflicting_taxid_prefixes(
    taxids_by_prefix: dict[tuple[tuple[str, str], ...], str],
) -> defaultdict[str, set[tuple[tuple[str, str], ...]]]:
    """Return the inverse mapping after rejecting reused native taxids."""
    prefixes_by_taxid: defaultdict[str, set[tuple[tuple[str, str], ...]]] = defaultdict(
        set,
    )
    for prefix, taxid in taxids_by_prefix.items():
        prefixes_by_taxid[taxid].add(prefix)
    conflicts = {
        taxid: prefixes
        for taxid, prefixes in prefixes_by_taxid.items()
        if len(prefixes) > 1
    }
    if not conflicts:
        return prefixes_by_taxid
    taxid = min(conflicts)
    prefixes = sorted(conflicts[taxid])
    examples = " || ".join(
        " > ".join(f"{rank}={name}" for rank, name in prefix) for prefix in prefixes[:2]
    )
    msg = f"taxid {taxid} represents conflicting lineage prefixes: {examples}"
    raise ValueError(msg)


def allocate_synthetic_taxids(
    prefixes: set[tuple[tuple[str, str], ...]],
    taxids_by_prefix: dict[tuple[tuple[str, str], ...], str],
    prefixes_by_taxid: defaultdict[str, set[tuple[tuple[str, str], ...]]],
) -> None:
    """Assign deterministic unused local IDs to prefixes without native evidence."""
    for prefix in sorted(prefixes.difference(taxids_by_prefix)):
        collision_attempt = 0
        while True:
            taxid = wvdb_local_taxid(prefix, collision_attempt=collision_attempt)
            if taxid not in prefixes_by_taxid:
                taxids_by_prefix[prefix] = taxid
                prefixes_by_taxid[taxid].add(prefix)
                break
            collision_attempt += 1


def assign_materialized_taxpaths(
    rows: list[dict[str, str]],
    taxids_by_prefix: dict[tuple[tuple[str, str], ...], str],
) -> None:
    """Serialize finalized rank names into fixed-width positional taxpaths."""
    for row in rows:
        prefix: list[tuple[str, str]] = []
        output_taxids: list[str] = []
        for rank in REFERENCE_TAXONOMY_RANKS:
            name = clean_cell(row[rank])
            if not name:
                output_taxids.append("")
                continue
            prefix.append((rank, name))
            output_taxids.append(taxids_by_prefix[tuple(prefix)])
        row[TAXPATH_COLUMN] = "|".join(output_taxids)


def materialize_curated_taxpaths(
    rows: list[dict[str, str]],
    diagnostics: list[dict[str, str]],
) -> None:
    """Align every final rank name with one stable positional taxid."""
    all_prefixes, candidates = source_taxid_candidates(rows)

    taxids_by_prefix: dict[tuple[tuple[str, str], ...], str] = {}
    for prefix in all_prefixes:
        taxid = select_source_taxid(prefix, candidates.get(prefix, []), diagnostics)
        if taxid is not None:
            taxids_by_prefix[prefix] = taxid
    prefixes_by_taxid = fail_on_conflicting_taxid_prefixes(taxids_by_prefix)
    allocate_synthetic_taxids(all_prefixes, taxids_by_prefix, prefixes_by_taxid)
    assign_materialized_taxpaths(rows, taxids_by_prefix)


def curate_taxonomy_names(
    rows: list[dict[str, str]],
    *,
    diagnostics: list[dict[str, str]],
) -> int:
    """Canonicalize conflicting taxon placements from broad to specific ranks."""
    resolved = 0
    reserved_names = {
        normalized_taxon_name(row[rank])
        for row in rows
        for rank in REFERENCE_TAXONOMY_RANKS
        if row[rank]
    }
    for rank_index, _rank in enumerate(REFERENCE_TAXONOMY_RANKS[1:], start=1):
        for name, candidates in taxonomy_conflicts_at_rank(rows, rank_index):
            winner, resolution = select_canonical_lineage(
                candidates,
                rank_index=rank_index,
            )
            if winner is None:
                contextualize_taxonomy_conflict(
                    candidates,
                    rank_index=rank_index,
                    reserved_names=reserved_names,
                )
                append_contextualization_diagnostics(
                    diagnostics,
                    candidates,
                    rank_index=rank_index,
                    name=name,
                )
                resolved += 1
                continue

            append_taxonomy_diagnostics(
                diagnostics,
                candidates,
                conflict=(rank_index, name, resolution),
                winner=winner,
            )
            resolved += 1
            for row in candidates:
                for ancestor_rank in REFERENCE_TAXONOMY_RANKS[: rank_index + 1]:
                    row[ancestor_rank] = winner[ancestor_rank]

    return resolved


def find_cross_rank_taxon_names(
    rows: list[dict[str, str]],
    *,
    diagnostics: list[dict[str, str]],
) -> list[str]:
    """Record names reused at multiple ranks, which Taxburst cannot represent."""
    occurrences: defaultdict[str, list[tuple[str, dict[str, str]]]] = defaultdict(list)
    for row in rows:
        for rank in REFERENCE_TAXONOMY_RANKS:
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
                    "resolution": "rejected-cross-rank-reuse",
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
) -> None:
    """Validate sketch identities and canonicalize conflicting taxonomy names."""
    manifest_rows = [
        *read_manifest_dicts(paths.reference_manifest_csv, source="ncbi"),
        *read_manifest_dicts(paths.wvdb_manifest_csv, source="wvdb"),
    ]
    manifest_identifier_counts = Counter(row["ident"] for row in manifest_rows)
    duplicate_manifest_identifiers = sorted(
        ident for ident, count in manifest_identifier_counts.items() if count > 1
    )
    if duplicate_manifest_identifiers:
        msg = f"Duplicate manifest identifier: {duplicate_manifest_identifiers[0]}"
        raise ValueError(msg)
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
    duplicate_error = diagnose_duplicate_signature_identities(
        manifest_rows,
        lineages,
        diagnostics,
    )
    if duplicate_error is not None:
        write_diagnostics(paths.diagnostics_tsv, diagnostics)
        raise ValueError(duplicate_error)

    retained_rows = [lineages[record["ident"]] for record in manifest_rows]
    resolved = curate_taxonomy_names(
        retained_rows,
        diagnostics=diagnostics,
    )
    cross_rank_conflicts = find_cross_rank_taxon_names(
        retained_rows,
        diagnostics=diagnostics,
    )
    if resolved:
        logger.warning(
            "resolved {} taxonomy conflict(s); evidence={}",
            resolved,
            paths.diagnostics_tsv,
        )
    if cross_rank_conflicts:
        write_diagnostics(paths.diagnostics_tsv, diagnostics)
        msg = (
            f"{len(cross_rank_conflicts)} cross-rank taxonomy name conflict(s). "
            f"Examples: {', '.join(cross_rank_conflicts[:10])}. "
            f"Evidence: {paths.diagnostics_tsv}"
        )
        raise ValueError(msg)

    materialize_curated_taxpaths(retained_rows, diagnostics)
    write_diagnostics(paths.diagnostics_tsv, diagnostics)
    write_curated_lineages(paths.lineages_csv, retained_rows)


def write_checksums(artifacts: list[Path], output: Path) -> None:
    """Write a SHA-256 ready marker for one validated artifact set."""
    lines = []
    for artifact in artifacts:
        with artifact.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        lines.append(f"{digest}  {artifact.name}\n")
    ensure_parent(output)
    output.write_text("".join(lines), encoding="utf-8")


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

    validate = subcommands.add_parser(
        "validate-lineages",
        help="Validate a combined lineage artifact for BioBoxes publication.",
    )
    validate.add_argument("--manifest-csv", type=Path, required=True)
    validate.add_argument("--lineages-csv", type=Path, required=True)

    checksums = subcommands.add_parser(
        "write-checksums",
        help="Write SHA-256 checksums after an artifact set passes validation.",
    )
    checksums.add_argument("--output", type=Path, required=True)
    checksums.add_argument("artifacts", type=Path, nargs="+")

    curate = subcommands.add_parser(
        "curate-lineages",
        help="Validate signature identity and resolve conflicting taxonomy names.",
    )
    curate.add_argument("--reference-manifest-csv", type=Path, required=True)
    curate.add_argument("--wvdb-manifest-csv", type=Path, required=True)
    curate.add_argument("--reference-lineages-csv", type=Path, required=True)
    curate.add_argument("--wvdb-lineages-csv", type=Path, required=True)
    curate.add_argument("--lineages-csv", type=Path, required=True)
    curate.add_argument("--diagnostics-tsv", type=Path, required=True)

    return parser


def main() -> None:
    """Run the command-line interface."""
    configure_logging()
    parser = build_parser()
    args = parser.parse_args()

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
    elif args.command == "validate-lineages":
        validate_lineages(args.manifest_csv, args.lineages_csv)
    elif args.command == "write-checksums":
        write_checksums(args.artifacts, args.output)
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
        )


if __name__ == "__main__":
    main()
