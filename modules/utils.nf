/*
 * Compute the deterministic run context without writing workflow state.
 */
process COMPUTE_RUN_CONTEXT {

    label "low"
    cache false

    input:
    tuple path(samplesheet), val(state_dir)

    output:
    val true, emit: ready
    tuple stdout, val(state_dir), emit: run_context

    script:
    """
    compute_run_context.py \
        --verbose \
        --samplesheet ${samplesheet}
    """
}

/*
 * Ensure the NCBI taxonomy cache is ready before distributed annotation workers run.
 */
process ENSURE_TAXONOMY {

    label "low"
    cache false

    input:
    val state_dir
    val taxonomy_dir

    output:
    val taxonomy_dir, emit: taxonomy_dir

    script:
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
    """
    ensure_taxonomy.py \
        --verbose \
        ${state_dir_arg} \
        ${taxonomy_dir_arg}
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

process ADD_READ_COUNTS_TO_BLAST {
    /*
     * Append total_reads column to merged BLAST results.
     */

    tag "${sample_id}"
    label "low"

    input:
    tuple val(sample_id), path(blast_tsv), val(total_reads)

    output:
    tuple val(sample_id), path("${sample_id}_blast.final.tsv")

    script:
    """
    awk -v reads="${total_reads}" 'BEGIN{OFS="\t"} NR==1{print \$0, "total_reads"} NR>1{print \$0, reads}' \
        ${blast_tsv} > ${sample_id}_blast.final.tsv
    """
}

process CONCATENATE_EXPERIMENT_BLAST {
    /*
     * Concatenate all per-sample _blast.final.tsv files into a single
     * experiment-level TSV. Runs unconditionally (does not require LabKey).
     * Header is taken from the first file; subsequent files contribute
     * data rows only.
     */

    label "low"

    input:
    path "sample_results/*"

    output:
    path "experiment_blast_results.tsv", emit: concatenated_tsv

    script:
    """
    #!/usr/bin/env python3
    from pathlib import Path

    files = sorted(Path("sample_results").glob("*.tsv"))

    if not files:
        Path("experiment_blast_results.tsv").write_text("")
    else:
        with open("experiment_blast_results.tsv", "w") as out:
            # Write header from first file
            with open(files[0]) as f:
                header = f.readline()
                out.write(header)
                for line in f:
                    out.write(line)
            # Append data rows from remaining files
            for tsv in files[1:]:
                with open(tsv) as f:
                    f.readline()  # skip header
                    for line in f:
                        out.write(line)
    """
}

/*
 * Send a minimal stateless Slack notification for run completion.
 */
process NOTIFY_SLACK {

    tag "${workflow.runName}"
    label "low"
    cache false

    errorStrategy 'ignore'

    secret 'SLACK_BOT_TOKEN'

    input:
    val ready
    tuple val(sample_set_id), val(state_dir)
    val labkey_url

    output:
    val true, emit: done

    when:
    params.slack_enabled && params.slack_channel

    script:
    """
    notify_slack.py \
        --run-id '${workflow.runName}' \
        --experiment-id '${params.experiment_id}' \
        --channel '${params.slack_channel}' \
        --sample-set-id '${sample_set_id}' \
        --labkey-url '${labkey_url}' \
        -v || echo "Slack notification failed (non-fatal)"
    """
}
