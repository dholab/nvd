/*
 * Check run state and register the run atomically.
 *
 * This process prevents duplicate processing of the same sample set.
 * It reads sample IDs from the samplesheet, computes a deterministic
 * sample_set_id, and attempts to register the run. If the sample set
 * was already processed, the pipeline fails fast with a clear error.
 *
 * Also registers each sample in processed_samples table with provenance
 * (database versions) for tracking.
 *
 * Input: 
 *   - samplesheet: path to samplesheet CSV
 *   - state_dir: path to state directory (for SQLite database)
 * Output: 
 *   - ready: val true (gates downstream processes, backward compatible)
 *   - sample_set_id: the computed sample set ID (for upload tracking)
 */
process CHECK_RUN_STATE {

    label "low"
    cache false  // Always run this check

    input:
    path samplesheet
    path state_dir

    output:
    val true, emit: ready
    stdout emit: sample_set_id

    script:
    def exp_id = params.experiment_id ?: "None"
    def blast_db_ver = params.blast_db_version ? "\"${params.blast_db_version}\"" : "None"
    def stat_db_ver = params.stat_db_version ? "\"${params.stat_db_version}\"" : "None"
    """
    #!/usr/bin/env python3
    import csv
    import sys
    from py_nvd.state import (
        compute_sample_set_id, 
        register_run, 
        get_run_by_sample_set,
        register_processed_sample,
    )

    # State directory is passed from Nextflow and staged/mounted automatically
    state_dir = "${state_dir}"

    # Read sample IDs from samplesheet
    sample_ids = []
    with open("${samplesheet}") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row.get("sample_id", "").strip()
            if sample_id and not sample_id.startswith("#"):
                sample_ids.append(sample_id)

    if not sample_ids:
        print("ERROR: No sample IDs found in samplesheet", file=sys.stderr)
        sys.exit(1)

    # Deduplicate and compute sample set ID
    sample_ids = sorted(set(sample_ids))
    sample_set_id = compute_sample_set_id(sample_ids)
    experiment_id = ${exp_id}

    # Check if this sample set was already processed
    existing = get_run_by_sample_set(sample_set_id, state_dir=state_dir)
    if existing:
        print(f"ERROR: Sample set already processed", file=sys.stderr)
        print(f"  Previous run: {existing.run_id}", file=sys.stderr)
        print(f"  Sample set ID: {sample_set_id}", file=sys.stderr)
        print(f"  Status: {existing.status}", file=sys.stderr)
        if existing.experiment_id:
            print(f"  Experiment ID: {existing.experiment_id}", file=sys.stderr)
        sys.exit(1)

    # Register this run
    run = register_run("${workflow.runName}", sample_set_id, experiment_id, state_dir=state_dir)
    if not run:
        print(f"ERROR: Failed to register run", file=sys.stderr)
        print(f"  Sample set {sample_set_id} may have been claimed by concurrent run", file=sys.stderr)
        sys.exit(1)

    # Register each sample with provenance
    blast_db_version = ${blast_db_ver}
    stat_db_version = ${stat_db_ver}
    
    for sid in sample_ids:
        register_processed_sample(
            sample_id=sid,
            sample_set_id=sample_set_id,
            run_id="${workflow.runName}",
            blast_db_version=blast_db_version,
            stat_db_version=stat_db_version,
            state_dir=state_dir,
        )

    # Log to stderr (doesn't affect stdout output)
    print(f"Registered run: ${workflow.runName}", file=sys.stderr)
    print(f"  State dir: {state_dir}", file=sys.stderr)
    print(f"  Sample set ID: {sample_set_id}", file=sys.stderr)
    print(f"  Samples: {len(sample_ids)}", file=sys.stderr)
    print(f"  BLAST DB: {blast_db_version}", file=sys.stderr)
    print(f"  STAT DB: {stat_db_version}", file=sys.stderr)
    if experiment_id:
        print(f"  Experiment ID: {experiment_id}", file=sys.stderr)

    # Output sample_set_id to stdout for Nextflow to capture
    # IMPORTANT: end="" prevents trailing newline - stdout emit captures raw text
    # Do not add any other print() calls to stdout in this script
    print(sample_set_id, end="")
    """
}

process ANNOTATE_LEAST_COMMON_ANCESTORS {

    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2
    
    input:
    tuple val(sample_id), path(all_blast_hits)

    output:
    tuple val(sample_id), path("${sample_id}_blast.merged_with_lca.tsv")

    script:
    """
    annotate_blast_lca.py -i ${all_blast_hits} -o ${sample_id}_blast.merged_with_lca.tsv
    """

}
