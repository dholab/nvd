include { ADD_READ_COUNTS_TO_BLAST; CONCATENATE_SAMPLE_BLAST_RESULTS; BUILD_QUERY_BIG_TABLE; BUILD_TAXON_BIG_TABLE; CONCATENATE_QUERY_BIG_TABLE; CONCATENATE_TAXON_BIG_TABLE; CONCATENATE_EXPERIMENT_BLAST_RESULTS; TARGET_ENRICHMENT_REPORT } from "../modules/utils"
include { BUILD_SEQUENCE_FLOW; RENDER_CONTIG_ALIGNMENT_PLOTS; RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST; RENDER_TAXON_ABUNDANCE_SUNBURST; RENDER_SOURMASH_SANKEY } from "../modules/reporting"
include { CRUMBS_PROFILING } from "./crumbs_profiling"
include { LIMS_INTEGRATION } from "./lims_integration"
include { RENDER_CONTIG_COVERAGE_HISTOGRAM } from "../modules/samtools"
include { ANNOTATE_BLAST_RISK_GROUPS } from "../modules/risk_groups"

workflow REPORTING {
    take:
    ch_blast_results
    ch_read_counts
    ch_contig_sequences  // tuple(sample_id, platform, read_structure, fasta, lookup)
    ch_query_fastas      // tuple(sample_id, platform, query_class, batch_fasta, lookup) from PREPARE_BLAST_QUERIES.out.queries
    ch_query_lookups
    ch_contig_read_counts
    ch_filtered_bam
    ch_no_contigs
    ch_target_enrichment_stats
    ch_taxonomy_dir
    ch_run_ready
    ch_run_context
    ch_sourmash_tax_reports
    ch_risk_group_lookup
    ch_sequence_flow_evidence
    run_id

    main:
    if (params.labkey) {
        NvdUtils.validateLabkeyBlast(params)
    }

    ch_contig_sequence_parts = ch_contig_sequences.multiMap { sample_id, _platform, _read_structure, fasta, _lookup ->
        for_lims: tuple(sample_id, fasta)
        for_alignment_plots: tuple(sample_id, fasta)
    }

    // Per-read-type query FASTAs (contig, merged, single) — the sequences
    // that were actually BLASTed — for eager per-read-type LIMS upload.
    ch_query_fastas_for_lims = ch_query_fastas.map { sample_id, _platform, query_class, batch_fasta, _lookup ->
        tuple(sample_id, query_class, batch_fasta)
    }

    ch_contig_alignment_plot_inputs = ch_contig_sequence_parts.for_alignment_plots
        .join(ch_filtered_bam, by: 0)
        .map { sample_id, fasta, bam, bai -> tuple(sample_id, fasta, bam, bai) }

    // Enrich BLAST results with all pipeline metadata (mapped_reads, total_reads,
    // blast_db_version, nextflow_run_id) so the published TSV is complete
    // regardless of whether LabKey is enabled. This now runs per (sample_id,
    // query_class) batch — before batches are stacked — so each read type's
    // hits can be uploaded eagerly downstream.
    ANNOTATE_BLAST_RISK_GROUPS(
        ch_blast_results,
        ch_risk_group_lookup,
        ch_taxonomy_dir,
    )

    ch_blast_finalize = ANNOTATE_BLAST_RISK_GROUPS.out
        .combine(ch_read_counts, by: 0)
        .combine(ch_contig_read_counts, by: 0)
        .combine(ch_query_lookups, by: 0)
        .map { sample_id, query_class, batch_tsv, total_reads, contig_counts, lookups ->
            tuple(sample_id, query_class, batch_tsv, total_reads, contig_counts, lookups)
        }

    def target_enrichment_enabled = NvdUtils.targetEnrichmentEnabled(params)
    ch_blast_db_version = Channel.value(params.blast_db_version)
    ch_virus_index_version = Channel.value(
        target_enrichment_enabled ? params.virus_index_version : "not_used"
    )
    ADD_READ_COUNTS_TO_BLAST(
        ch_blast_finalize,
        run_id,
        ch_blast_db_version,
        ch_virus_index_version,
    )

    // Collapse query-class partitions into the canonical per-sample BLAST result.
    CONCATENATE_SAMPLE_BLAST_RESULTS(
        ADD_READ_COUNTS_TO_BLAST.out
            .map { sample_id, _query_class, final_tsv -> tuple(sample_id, final_tsv) }
            .groupTuple()
    )

    RENDER_CONTIG_COVERAGE_HISTOGRAM(ch_filtered_bam)
    RENDER_CONTIG_ALIGNMENT_PLOTS(ch_contig_alignment_plot_inputs)

    ch_sample_blast_results = CONCATENATE_SAMPLE_BLAST_RESULTS.out
        .multiMap { sample_id, blast_tsv ->
            for_summary: tuple(sample_id, blast_tsv)
            for_big_table: tuple(sample_id, blast_tsv)
            for_emit: tuple(sample_id, blast_tsv)
            for_lims: tuple(sample_id, blast_tsv)
        }

    // Concatenate all per-sample final BLAST results into a single experiment-level TSV.
    // Runs unconditionally so every run produces an experiment summary, not just LabKey runs.
    CONCATENATE_EXPERIMENT_BLAST_RESULTS(
        ch_sample_blast_results.for_summary.map { _sample_id, tsv -> tsv }.collect()
    )

    if (params.experimental == true) {
        BUILD_SEQUENCE_FLOW(ch_sequence_flow_evidence.collect())

        CRUMBS_PROFILING(
            ch_sample_blast_results.for_emit,
            ch_filtered_bam,
            ch_no_contigs,
            ch_taxonomy_dir,
        )

        ch_query_big_table_inputs = ch_sample_blast_results.for_big_table
            .join(CRUMBS_PROFILING.out.queries, by: 0)
            .map { sample_id, blast_tsv, crumbs_tsv -> tuple(sample_id, blast_tsv, crumbs_tsv) }

        BUILD_QUERY_BIG_TABLE(ch_query_big_table_inputs)

        CONCATENATE_QUERY_BIG_TABLE(
            BUILD_QUERY_BIG_TABLE.out.map { _sample_id, tsv -> tsv }.collect()
        )

        ch_taxon_big_table_inputs = BUILD_QUERY_BIG_TABLE.out
            .join(CRUMBS_PROFILING.out.taxa, by: 0)
            .map { sample_id, query_big_table, crumbs_taxa_tsv -> tuple(sample_id, query_big_table, crumbs_taxa_tsv) }

        BUILD_TAXON_BIG_TABLE(ch_taxon_big_table_inputs)

        CONCATENATE_TAXON_BIG_TABLE(
            BUILD_TAXON_BIG_TABLE.out.map { _sample_id, tsv -> tsv }.collect()
        )
    }

    TARGET_ENRICHMENT_REPORT(
        ch_target_enrichment_stats.map { _sample_id, json -> json }.collect()
    )

    if (params.experimental == true) {
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

    }

    LIMS_INTEGRATION(
        ADD_READ_COUNTS_TO_BLAST.out,
        ch_sample_blast_results.for_lims,
        ch_contig_sequence_parts.for_lims,
        ch_query_fastas_for_lims,
        params.experiment_id,
        run_id,
        ch_run_ready,
    )

    // Completion means every reporting leaf has settled. Keep this terminal
    // inventory local to REPORTING so adding a report does not spread another
    // synchronization channel through the composition root.
    ch_reporting_terminal_outputs = CONCATENATE_EXPERIMENT_BLAST_RESULTS.out.concatenated_tsv
        .mix(RENDER_CONTIG_COVERAGE_HISTOGRAM.out.histogram)
        .mix(RENDER_CONTIG_ALIGNMENT_PLOTS.out.plots)
        .mix(TARGET_ENRICHMENT_REPORT.out.summary_tsv)

    if (params.experimental == true) {
        ch_reporting_terminal_outputs = ch_reporting_terminal_outputs
            .mix(BUILD_SEQUENCE_FLOW.out.sequence_flow)
            .mix(CONCATENATE_QUERY_BIG_TABLE.out.concatenated_tsv)
            .mix(CONCATENATE_TAXON_BIG_TABLE.out.concatenated_tsv)
            .mix(CRUMBS_PROFILING.out.krona)
            .mix(CRUMBS_PROFILING.out.taxburst)
            .mix(CRUMBS_PROFILING.out.merged_taxburst)
            .mix(RENDER_TAXON_ABUNDANCE_SUNBURST.out.reports)
            .mix(RENDER_SOURMASH_SANKEY.out.report)
            .mix(RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST.out.report)
    }

    if (params.labkey) {
        ch_reporting_terminal_outputs = ch_reporting_terminal_outputs
            .mix(LIMS_INTEGRATION.out.final_labkey_log)
    }

    ch_reporting_ready = ch_reporting_terminal_outputs
        .collect()
        .map { _outputs -> true }
    ch_labkey_ready = params.labkey
        ? LIMS_INTEGRATION.out.uploads_done
        : channel.value(true)
    ch_completed_results = CONCATENATE_EXPERIMENT_BLAST_RESULTS.out.concatenated_tsv
        .combine(ch_reporting_ready)
        .combine(ch_labkey_ready)
        .map { experiment_results, _reporting_ready, _labkey_ready -> experiment_results }

    emit:
    blast_results = ch_sample_blast_results.for_emit
    query_big_tables = params.experimental ? BUILD_QUERY_BIG_TABLE.out : channel.empty()
    query_big_table = params.experimental ? CONCATENATE_QUERY_BIG_TABLE.out.concatenated_tsv : channel.empty()
    taxon_big_tables = params.experimental ? BUILD_TAXON_BIG_TABLE.out : channel.empty()
    taxon_big_table = params.experimental ? CONCATENATE_TAXON_BIG_TABLE.out.concatenated_tsv : channel.empty()
    experiment_blast = CONCATENATE_EXPERIMENT_BLAST_RESULTS.out.concatenated_tsv
    completed_results = ch_completed_results
    sequence_flow = params.experimental ? BUILD_SEQUENCE_FLOW.out.sequence_flow : channel.empty()
    target_enrichment_report = TARGET_ENRICHMENT_REPORT.out.summary_tsv
    taxon_abundance_sunbursts = params.experimental ? RENDER_TAXON_ABUNDANCE_SUNBURST.out.reports : channel.empty()
    merged_taxon_abundance_sunburst = params.experimental ? RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST.out.report : channel.empty()
    sourmash_sankey_reports = params.experimental ? RENDER_SOURMASH_SANKEY.out.report : channel.empty()
    labkey_log = LIMS_INTEGRATION.out.upload_log
    final_labkey_log = LIMS_INTEGRATION.out.final_labkey_log
    labkey_uploads_done = LIMS_INTEGRATION.out.uploads_done
    crumbs_queries = params.experimental ? CRUMBS_PROFILING.out.queries : channel.empty()
    crumbs_taxa = params.experimental ? CRUMBS_PROFILING.out.taxa : channel.empty()
    crumbs_bioboxes_profile = params.experimental ? CRUMBS_PROFILING.out.bioboxes_profile : channel.empty()
    crumbs_qc = params.experimental ? CRUMBS_PROFILING.out.qc : channel.empty()
    crumbs_krona = params.experimental ? CRUMBS_PROFILING.out.krona : channel.empty()
    crumbs_kreport = params.experimental ? CRUMBS_PROFILING.out.kreport : channel.empty()
    crumbs_taxburst = params.experimental ? CRUMBS_PROFILING.out.taxburst : channel.empty()
    merged_crumbs_taxburst = params.experimental ? CRUMBS_PROFILING.out.merged_taxburst : channel.empty()
    crumbs_profile_taxonomy = params.experimental ? CRUMBS_PROFILING.out.profile_taxonomy : channel.empty()
}
