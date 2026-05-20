include { ADD_READ_COUNTS_TO_BLAST; CONCATENATE_EXPERIMENT_BLAST } from "../modules/utils"
include { NOTIFY_SLACK } from "../modules/utils"
include { LIMS_INTEGRATION } from "./lims_integration"

workflow REPORTING {
    take:
    ch_blast_results
    ch_read_counts
    ch_contig_sequences
    ch_contig_read_counts
    ch_run_ready
    ch_run_context
    run_id

    main:
    if (params.labkey) {
        NvdUtils.validateLabkeyBlast(params)
    }

    // Add total read counts to the final merged BLAST results so they are
    // always present in the published TSV, regardless of LabKey.
    ch_split_read_counts = ch_read_counts
        .multiMap { sample_id, total_reads ->
            for_blast: tuple(sample_id, total_reads)
            for_labkey: tuple(sample_id, total_reads)
        }

    ch_blast_with_counts = ch_blast_results
        .join(ch_split_read_counts.for_blast, by: 0)

    ADD_READ_COUNTS_TO_BLAST(ch_blast_with_counts)

    ch_split_blast_results = ADD_READ_COUNTS_TO_BLAST.out
        .multiMap { sample_id, blast_tsv ->
            for_summary: tuple(sample_id, blast_tsv)
            for_labkey_trigger: tuple(sample_id, blast_tsv)
            for_labkey_upload: tuple(sample_id, blast_tsv)
            for_emit: tuple(sample_id, blast_tsv)
        }

    // Concatenate all per-sample final BLAST results into a single experiment-level TSV.
    // Runs unconditionally so every run produces an experiment summary, not just LabKey runs.
    CONCATENATE_EXPERIMENT_BLAST(
        ch_split_blast_results.for_summary.map { _sample_id, tsv -> tsv }.collect()
    )

    LIMS_INTEGRATION(
        ch_split_blast_results.for_labkey_upload,
        ch_contig_sequences,
        ch_split_read_counts.for_labkey,
        params.experiment_id,
        run_id,
        ch_contig_read_counts,
        ch_run_ready,
    )

    ch_labkey_url = channel.value(
        "https://${params.labkey_server}/${params.labkey_project_name}/list-grid.view?name=${params.labkey_blast_meta_hits_list}"
    )

    ch_slack_trigger = params.slack_enabled && params.slack_channel && params.labkey
        ? LIMS_INTEGRATION.out.registered
        : channel.empty()

    NOTIFY_SLACK(
        ch_slack_trigger,
        ch_run_context,
        ch_labkey_url,
    )

    emit:
    blast_results = ch_split_blast_results.for_emit
    experiment_blast = CONCATENATE_EXPERIMENT_BLAST.out.concatenated_tsv
    labkey_log = LIMS_INTEGRATION.out.upload_log
    final_labkey_log = LIMS_INTEGRATION.out.final_labkey_log
    labkey_registered = LIMS_INTEGRATION.out.registered
}
