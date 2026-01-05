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
 *   - tuple of (samplesheet, state_dir) - combined to ensure they travel together
 * Output: 
 *   - ready: val true (gates downstream processes, for combine() with other validation gates)
 *   - run_context: tuple of (sample_set_id, state_dir) - bundled for downstream upload processes
 *
 * Design note: We emit run_context as a tuple so that sample_set_id and state_dir
 * travel together to downstream processes. This prevents the Nextflow footgun where
 * separate queue channels get consumed by the first process, leaving nothing for others.
 */
process CHECK_RUN_STATE {

    label "low"
    cache false  // Always run this check

    input:
    tuple path(samplesheet), path(state_dir)

    output:
    val true, emit: ready
    tuple stdout, path(state_dir), emit: run_context

    script:
    def exp_arg = params.experiment_id ? "--experiment-id ${params.experiment_id}" : ""
    def blast_db_arg = params.blast_db_version ? "--blast-db-version '${params.blast_db_version}'" : ""
    def stat_db_arg = params.stat_db_version ? "--stat-db-version '${params.stat_db_version}'" : ""
    """
    check_run_state.py \\
        --verbose \\
        --samplesheet ${samplesheet} \\
        --run-id '${workflow.runName}' \\
        --state-dir ${state_dir} \\
        ${exp_arg} \\
        ${blast_db_arg} \\
        ${stat_db_arg}
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

/*
 * Register BLAST hits with idempotent keys in the state database.
 *
 * This process computes deterministic hit keys from contig sequences,
 * enabling cross-run and cross-sample deduplication. The hit key is
 * derived from the canonical sequence (lexicographically smaller of
 * the sequence and its reverse complement).
 *
 * Input:
 *   - tuple of (sample_id, contigs, blast_results, sample_set_id, state_dir)
 * Output:
 *   - tuple of (sample_id, log_file)
 *
 * Note: This is currently a "dead end" - output is not consumed by
 * downstream processes. The hits are registered in the state database
 * for later querying via `nvd hits export`.
 */
process REGISTER_HITS {

    tag "${sample_id}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(contigs), path(blast_results), val(sample_set_id), path(state_dir)

    output:
    tuple val(sample_id), path("${sample_id}_hits_registered.log")

    script:
    def blast_db_arg = params.blast_db_version ? "--blast-db-version '${params.blast_db_version}'" : ""
    def stat_db_arg = params.stat_db_version ? "--stat-db-version '${params.stat_db_version}'" : ""
    """
    register_hits.py \\
        --contigs ${contigs} \\
        --blast-results ${blast_results} \\
        --state-dir ${state_dir} \\
        --sample-set-id '${sample_set_id}' \\
        --sample-id '${sample_id}' \\
        --run-id '${workflow.runName}' \\
        ${blast_db_arg} \\
        ${stat_db_arg} \\
        -v \\
        > ${sample_id}_hits_registered.log 2>&1
    """
}

/*
 * Mark a pipeline run as completed or failed.
 *
 * This process should be called at the end of the workflow, gated on
 * completion of all sample processing. It updates the run status in
 * the state database.
 *
 * Input:
 *   - ready: Gate signal (true when all processing is complete)
 *   - state_dir: Path to state directory
 *   - status: "completed" or "failed"
 *
 * Output:
 *   - done: val true (can be used to gate downstream processes)
 *
 * Note: This process is intentionally not cached (cache false) because
 * the run status should always be updated at the end of a workflow run.
 */
process COMPLETE_RUN {

    tag "${workflow.runName}"
    label "low"
    cache false  // Always run this at workflow end

    input:
    val ready  // Gate: all processing complete
    path state_dir
    val status  // "completed" or "failed"

    output:
    val true, emit: done

    script:
    """
    complete_run.py \\
        --run-id '${workflow.runName}' \\
        --state-dir ${state_dir} \\
        --status '${status}' \\
        -v
    """
}
