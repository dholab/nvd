include { ADD_READ_COUNTS_TO_BLAST; CONCATENATE_EXPERIMENT_BLAST; VIRUS_ENRICHMENT_REPORT } from "../modules/utils"
include { NOTIFY_SLACK } from "../modules/utils"
include { RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST; RENDER_TAXON_ABUNDANCE_SUNBURST; RENDER_SOURMASH_SANKEY } from "../modules/reporting"
include { CRUMBS_PROFILING } from "./crumbs_profiling"
include { LIMS_INTEGRATION } from "./lims_integration"

workflow REPORTING {
    take:
    ch_blast_results
    ch_read_counts
    ch_contig_sequences
    ch_contig_read_counts
    ch_filtered_bam
    ch_virus_enrichment_stats
    ch_taxonomy_dir
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

    ch_sourmash_profile_summaries = ch_sourmash_tax_reports
        .map { sample_id, platform, read_structure, tax_reports ->
            def report_files = tax_reports instanceof List ? tax_reports : [tax_reports]
            def summarized_csv = report_files.find { report -> report.name.endsWith(".summarized.csv") }
            assert summarized_csv : "Missing taxonomic profile summary CSV for ${sample_id}."
            tuple(sample_id, platform, read_structure, summarized_csv)
        }

    RENDER_TAXON_ABUNDANCE_SUNBURST(ch_sourmash_profile_summaries)
    RENDER_SOURMASH_SANKEY(ch_sourmash_profile_summaries)

    ch_merged_taxburst_input = ch_sourmash_profile_summaries
        .map { sample_id, _platform, _read_structure, profile_summary -> tuple("all", sample_id, profile_summary) }
        .groupTuple()
        .map { _key, sample_ids, profile_summaries -> tuple(sample_ids, profile_summaries) }

    RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST(ch_merged_taxburst_input)

    CRUMBS_PROFILING(
        ch_split_blast_results.for_emit,
        ch_filtered_bam,
        ch_taxonomy_dir,
    )

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
    crumbs_contigs = CRUMBS_PROFILING.out.contigs
    crumbs_taxa = CRUMBS_PROFILING.out.taxa
    crumbs_bioboxes_profile = CRUMBS_PROFILING.out.bioboxes_profile
    crumbs_qc = CRUMBS_PROFILING.out.qc
    crumbs_profile_taxonomy = CRUMBS_PROFILING.out.profile_taxonomy
}
