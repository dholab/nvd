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
 *   - tuple of (samplesheet, state_dir, taxonomy_dir) - combined to ensure they travel together
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
    tuple path(samplesheet), val(state_dir), val(taxonomy_dir)

    output:
    val true, emit: ready
    tuple stdout, val(state_dir), emit: run_context

    script:
    def exp_arg = params.experiment_id ? "--experiment-id ${params.experiment_id}" : ""
    def blast_db_arg = params.blast_db_version ? "--blast-db-version '${params.blast_db_version}'" : ""
    def stat_db_arg = params.stat_db_version ? "--stat-db-version '${params.stat_db_version}'" : ""
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
    """
    check_run_state.py \\
        --verbose \\
        --samplesheet ${samplesheet} \\
        --run-id '${workflow.runName}' \\
        ${state_dir_arg} \\
        ${taxonomy_dir_arg} \\
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
    val state_dir
    val taxonomy_dir

    output:
    tuple val(sample_id), path("${sample_id}_blast.merged_with_lca.tsv")

    script:
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
    """
    annotate_blast_lca.py -i ${all_blast_hits} -o ${sample_id}_blast.merged_with_lca.tsv ${state_dir_arg} ${taxonomy_dir_arg}
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
 *   - tuple of (sample_id, contigs, blast_results, sample_set_id, state_dir, hits_dir)
 *   - state_dir may be null in stateless mode
 *   - hits_dir is the base directory for parquet output (state_dir/hits or results/hits)
 * Output:
 *   - tuple of (sample_id, log_file, parquet_file)
 *
 * Note: This is currently a "dead end" - output is not consumed by
 * downstream processes. The hits are registered in the state database
 * for later querying via `nvd hits export`.
 */
process REGISTER_HITS {

    tag "${sample_id}"
    label "low"
    maxForks 1

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    // Publish parquet file - runs on every resume even if task is cached
    // Uses Hive-partitioned structure: {hits_dir}/month=NULL/{sample_set_id}/{sample_id}/data.parquet
    // hits_dir is either state_dir/hits (stateful) or results/hits (stateless)
    publishDir "${hits_dir}/month=NULL/${sample_set_id}/${sample_id}", mode: 'copy', pattern: "data.parquet"

    input:
    tuple val(sample_id), path(contigs), path(blast_results), val(sample_set_id), val(state_dir), val(hits_dir)

    output:
    tuple val(sample_id), path("${sample_id}_hits_registered.log"), path("data.parquet")

    script:
    assert hits_dir : "hits_dir cannot be null - this indicates a workflow configuration error"
    def blast_db_arg = params.blast_db_version ? "--blast-db-version '${params.blast_db_version}'" : ""
    def stat_db_arg = params.stat_db_version ? "--stat-db-version '${params.stat_db_version}'" : ""
    def labkey_arg = params.labkey ? "--labkey" : ""
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    """
    set -o pipefail
    register_hits.py \\
        --contigs ${contigs} \\
        --blast-results ${blast_results} \\
        ${state_dir_arg} \\
        --sample-set-id '${sample_set_id}' \\
        --sample-id '${sample_id}' \\
        --run-id '${workflow.runName}' \\
        --output data.parquet \\
        ${blast_db_arg} \\
        ${stat_db_arg} \\
        ${labkey_arg} \\
        -v 2>&1 | tee ${sample_id}_hits_registered.log
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
 *   - state_dir: Path to state directory (may be null in stateless mode)
 *   - status: "completed" or "failed"
 *
 * Output:
 *   - done: val true (can be used to gate downstream processes)
 *
 * Note: This process is intentionally not cached (cache false) because
 * the run status should always be updated at the end of a workflow run.
 * In stateless mode (state_dir is null), this is a no-op.
 */
process COMPLETE_RUN {

    tag "${workflow.runName}"
    label "low"
    cache false  // Always run this at workflow end

    input:
    val ready  // Gate: all processing complete
    val state_dir
    val status  // "completed" or "failed"

    output:
    val true, emit: done

    script:
    def backup_arg = status == "completed" ? "--backup" : ""
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    """
    complete_run.py \\
        --run-id '${workflow.runName}' \\
        ${state_dir_arg} \\
        --status '${status}' \\
        ${backup_arg} \\
        -v
    """
}

/*
 * Send Slack notification for run completion.
 *
 * This process runs after COMPLETE_RUN and sends a summary to Slack.
 * It is designed to never fail the pipeline - errors are logged as warnings.
 *
 * IMPORTANT: This process runs on the local executor (not submitted to the
 * cluster) because it needs access to the state database and py_nvd modules
 * in the local pixi environment. This is configured in:
 *   - conf/chtc-template.config (for CHTC users via `nvd setup`)
 *   - Any custom user.config for other HPC environments
 *
 * Requires:
 *   - params.slack_enabled = true
 *   - params.slack_channel set to a valid channel ID
 *   - params.experiment_id set
 *   - SLACK_BOT_TOKEN secret configured via `nvd secrets`
 *
 * Input:
 *   - ready: Gate signal from COMPLETE_RUN
 *   - run_context: tuple of (sample_set_id, state_dir)
 *   - labkey_url: URL to LabKey results
 */
process NOTIFY_SLACK {

    tag "${workflow.runName}"
    label "low"
    cache false  // Always run at workflow end

    // Never fail the pipeline due to Slack issues
    errorStrategy 'ignore'

    // Make SLACK_BOT_TOKEN available as environment variable
    secret 'SLACK_BOT_TOKEN'

    input:
    val ready  // Gate: COMPLETE_RUN finished
    tuple val(sample_set_id), val(state_dir)
    val labkey_url

    output:
    val true, emit: done

    when:
    params.slack_enabled && params.slack_channel

    script:
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    """
    notify_slack.py \\
        --run-id '${workflow.runName}' \\
        --experiment-id '${params.experiment_id}' \\
        --channel '${params.slack_channel}' \\
        ${state_dir_arg} \\
        --sample-set-id '${sample_set_id}' \\
        --labkey-url '${labkey_url}' \\
        -v || echo "Slack notification failed (non-fatal)"
    """
}
