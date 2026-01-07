-- NVD State Database Schema
-- 
-- This is a local cache/state database, not a production database.
-- If the schema version changes, the database is deleted and recreated.
--
-- Phase 5 (2024-12-28): Redesigned sample and upload tracking.
-- - Replaced `samples` with `processed_samples` (tracks processing with provenance)
-- - Replaced `labkey_uploads` with `uploads` (target-agnostic, enables cross-run dedup)
-- - Key insight: indexing uploads by sample_id enables "was this sample ever uploaded?"
--   queries across all runs, solving the overlapping sample set problem.

PRAGMA user_version = 2;

-- Run tracking
CREATE TABLE IF NOT EXISTS runs (
    run_id TEXT PRIMARY KEY,
    sample_set_id TEXT NOT NULL UNIQUE,  -- SHA256[:16] of sorted sample IDs (idempotency key)
    experiment_id INTEGER,                -- Optional LabKey linkage
    started_at TEXT NOT NULL,
    completed_at TEXT,
    status TEXT NOT NULL CHECK (status IN ('running', 'completed', 'failed'))
);

-- Processed samples with provenance
-- Tracks which samples were processed, when, and with what database versions.
-- Keyed by (sample_id, sample_set_id) to allow intentional reprocessing with different samples.
CREATE TABLE IF NOT EXISTS processed_samples (
    sample_id TEXT NOT NULL,
    sample_set_id TEXT NOT NULL,
    run_id TEXT NOT NULL REFERENCES runs(run_id),
    processed_at TEXT NOT NULL,
    blast_db_version TEXT,
    stat_db_version TEXT,
    taxonomy_hash TEXT,
    status TEXT NOT NULL CHECK (status IN ('completed', 'uploaded', 'failed')),
    PRIMARY KEY (sample_id, sample_set_id)
);

-- Index for "was this sample ever processed?" queries across all runs
CREATE INDEX IF NOT EXISTS idx_processed_samples_sample_id 
    ON processed_samples(sample_id);

-- Reference database registry
CREATE TABLE IF NOT EXISTS databases (
    db_type TEXT NOT NULL CHECK (db_type IN ('blast', 'stat', 'gottcha2', 'hostile')),
    version TEXT NOT NULL,
    path TEXT NOT NULL,
    checksum TEXT,
    registered_at TEXT NOT NULL,
    PRIMARY KEY (db_type, version)
);

-- Index for looking up databases by path (for version auto-resolution)
CREATE INDEX IF NOT EXISTS idx_databases_path ON databases(db_type, path);

-- Upload tracking (target-agnostic)
-- Tracks where results were uploaded with content hashing for idempotency.
-- Keyed by (sample_id, sample_set_id, upload_type, upload_target) to track each upload.
-- Target-specific metadata (e.g., experiment_id, labkey_row_id) stored as JSON.
CREATE TABLE IF NOT EXISTS uploads (
    sample_id TEXT NOT NULL,
    sample_set_id TEXT NOT NULL,
    upload_type TEXT NOT NULL CHECK (upload_type IN ('blast', 'blast_fasta', 'gottcha2', 'gottcha2_fasta')),
    upload_target TEXT NOT NULL,          -- 'labkey', 'local', 'globus', etc.
    content_hash TEXT NOT NULL,           -- SHA256 of uploaded data for idempotency
    uploaded_at TEXT NOT NULL,
    target_metadata TEXT,                 -- JSON for target-specific data (experiment_id, row_id, etc.)
    PRIMARY KEY (sample_id, sample_set_id, upload_type, upload_target),
    FOREIGN KEY (sample_id, sample_set_id) REFERENCES processed_samples(sample_id, sample_set_id)
);

-- Index for "was this sample ever uploaded?" queries across all runs
-- This is the key index that enables cross-run sample deduplication
CREATE INDEX IF NOT EXISTS idx_uploads_sample_id 
    ON uploads(sample_id);

-- Index for checking uploads by target (e.g., all LabKey uploads)
CREATE INDEX IF NOT EXISTS idx_uploads_target 
    ON uploads(upload_target);

-- Taxonomy version tracking (for reproducibility)
-- Records which taxonomy database version was used for each run.
-- Enables detection of taxonomy drift between runs.
CREATE TABLE IF NOT EXISTS taxonomy_versions (
    run_id TEXT PRIMARY KEY REFERENCES runs(run_id),
    file_hash TEXT NOT NULL,           -- SHA256 of gettax.sqlite file
    ncbi_rebuild_timestamp TEXT,       -- From rebuilds2 table in gettax.sqlite
    file_size INTEGER,                 -- File size in bytes
    downloaded_at TEXT NOT NULL,       -- When we downloaded this version
    recorded_at TEXT NOT NULL          -- When this record was created
);

-- SRA download cache
CREATE TABLE IF NOT EXISTS sra_cache (
    srr_accession TEXT PRIMARY KEY,
    download_path TEXT NOT NULL,
    downloaded_at TEXT NOT NULL,
    file_count INTEGER,
    total_bytes INTEGER
);

-- Run parameter presets
-- Stores named configurations of Nextflow parameters for reuse
CREATE TABLE IF NOT EXISTS presets (
    name TEXT PRIMARY KEY,
    description TEXT,
    params TEXT NOT NULL,  -- JSON object of param_name -> value
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);

-- Index for listing presets by creation date
CREATE INDEX IF NOT EXISTS idx_presets_created_at ON presets(created_at);

-- Unique biological sequences (deduplicated by canonical sequence hash)
-- A "hit" is a contig that had BLAST results. The hit_key is computed from
-- the canonical sequence (lexicographically smaller of seq vs reverse-complement),
-- making it strand-agnostic and deterministic across runs.
CREATE TABLE IF NOT EXISTS hits (
    hit_key TEXT PRIMARY KEY,              -- BLAKE3 hash of canonical sequence (32 hex chars)
    sequence_length INTEGER NOT NULL,      -- original length in bp
    sequence_compressed BLOB NOT NULL,     -- 2-bit + zlib compressed (with N positions)
    gc_content REAL NOT NULL,              -- GC fraction (0.0-1.0)
    first_seen_date TEXT NOT NULL          -- ISO8601
);

-- Records of when/where each hit was observed (many-to-one with hits)
-- Enables queries like "which samples have seen this hit?" and 
-- "what hits appeared in this run?"
CREATE TABLE IF NOT EXISTS hit_observations (
    id INTEGER PRIMARY KEY,
    hit_key TEXT NOT NULL REFERENCES hits(hit_key),
    sample_set_id TEXT NOT NULL,           -- links to runs table
    sample_id TEXT NOT NULL,               -- source sample
    run_date TEXT NOT NULL,                -- ISO8601
    contig_id TEXT,                        -- original SPAdes contig ID (for traceability)
    UNIQUE(hit_key, sample_set_id, sample_id)
);

CREATE INDEX IF NOT EXISTS idx_hit_observations_sample_set 
    ON hit_observations(sample_set_id);
CREATE INDEX IF NOT EXISTS idx_hit_observations_sample 
    ON hit_observations(sample_id);
CREATE INDEX IF NOT EXISTS idx_hit_observations_hit_key 
    ON hit_observations(hit_key);

-- Sample-level processing locks to prevent duplicate work across concurrent runs.
-- Locks are acquired at run start and released on completion/upload.
-- TTL-based expiration handles crashed runs (default 72 hours).
-- Machine fingerprint (hostname + username) enables:
--   - Legitimate resume: same run_id + same fingerprint → refresh TTL
--   - Conflict detection: same run_id + different fingerprint + active → error
--   - Crash recovery: same run_id + different fingerprint + expired → allow
CREATE TABLE IF NOT EXISTS sample_locks (
    sample_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES runs(run_id),
    hostname TEXT NOT NULL,
    username TEXT NOT NULL,
    locked_at TEXT NOT NULL,
    expires_at TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_sample_locks_run_id 
    ON sample_locks(run_id);
CREATE INDEX IF NOT EXISTS idx_sample_locks_expires_at 
    ON sample_locks(expires_at);
