include { ADD_READ_COUNTS_TO_BLAST; CONCATENATE_EXPERIMENT_BLAST; VIRUS_ENRICHMENT_REPORT } from "../modules/utils"
include { NOTIFY_SLACK } from "../modules/utils"
include { BUILD_TAXONOMIC_PROFILE_NORMALIZATION_MAP; NORMALIZE_TAXONOMIC_PROFILE_SUMMARY; RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST; RENDER_TAXON_ABUNDANCE_SUNBURST; RENDER_SOURMASH_SANKEY } from "../modules/reporting"
include { LIMS_INTEGRATION } from "./lims_integration"

workflow REPORTING {
    take:
    ch_blast_results
    ch_read_counts
    ch_contig_sequences
    ch_contig_read_counts
    ch_virus_enrichment_stats
    ch_run_ready
    ch_run_context
    ch_sourmash_tax_reports
    run_id

    main:
    if (params.labkey) {
        NvdUtils.validateLabkeyBlast(params)
    }

    // Enrich BLAST results with all pipeline metadata (mapped_reads, total_reads,
    // blast_db_version, nextflow_run_id) so the published TSV is complete
    // regardless of whether LabKey is enabled.
    ch_blast_finalize = ch_blast_results
        .join(ch_read_counts, by: 0)
        .join(ch_contig_read_counts, by: 0)

    ADD_READ_COUNTS_TO_BLAST(ch_blast_finalize, run_id)

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

    VIRUS_ENRICHMENT_REPORT(
        ch_virus_enrichment_stats.map { _sample_id, json -> json }.collect()
    )

    ch_sourmash_tax_report_files = ch_sourmash_tax_reports
        .map { _sample_id, _platform, _read_structure, tax_reports -> tax_reports }
        .collect()
        .filter { tax_reports -> tax_reports }

    BUILD_TAXONOMIC_PROFILE_NORMALIZATION_MAP(ch_sourmash_tax_report_files)

    ch_sourmash_tax_reports_with_normalization_map = ch_sourmash_tax_reports
        .combine(BUILD_TAXONOMIC_PROFILE_NORMALIZATION_MAP.out.map)

    NORMALIZE_TAXONOMIC_PROFILE_SUMMARY(ch_sourmash_tax_reports_with_normalization_map)
    RENDER_TAXON_ABUNDANCE_SUNBURST(NORMALIZE_TAXONOMIC_PROFILE_SUMMARY.out.summary)
    RENDER_SOURMASH_SANKEY(NORMALIZE_TAXONOMIC_PROFILE_SUMMARY.out.summary)

    ch_merged_taxburst_input = NORMALIZE_TAXONOMIC_PROFILE_SUMMARY.out.summary
        .map { sample_id, _platform, _read_structure, profile_summary -> tuple("all", sample_id, profile_summary) }
        .groupTuple()
        .map { _key, sample_ids, profile_summaries -> tuple(sample_ids, profile_summaries) }

    RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST(ch_merged_taxburst_input)

    LIMS_INTEGRATION(
        ch_split_blast_results.for_labkey_upload,
        ch_contig_sequences,
        params.experiment_id,
        run_id,
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
    virus_enrichment_report = VIRUS_ENRICHMENT_REPORT.out.summary_tsv
    taxon_abundance_sunbursts = RENDER_TAXON_ABUNDANCE_SUNBURST.out.reports
    merged_taxon_abundance_sunburst = RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST.out.report
    sourmash_sankey_reports = RENDER_SOURMASH_SANKEY.out.report
    labkey_log = LIMS_INTEGRATION.out.upload_log
    final_labkey_log = LIMS_INTEGRATION.out.final_labkey_log
    labkey_registered = LIMS_INTEGRATION.out.registered
}
