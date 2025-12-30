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
