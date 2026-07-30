include { BUILD_RAPID_SCREENING_EVAL_DB ; RENDER_RAPID_SCREENING_EVAL } from "../modules/rapid_screening_eval"

workflow RAPID_SCREENING_EVAL {
    take:
    ch_sample_read_counts
    ch_sourmash_gather_csv
    ch_sourmash_tax_reports
    ch_sourmash_lineages
    ch_blast_final_tsv
    ch_crumbs_taxa
    ch_crumbs_contigs
    run_id

    main:
    ch_sample_read_counts_tsv = ch_sample_read_counts
        .map { sample_id, total_reads -> "${sample_id}\t${total_reads}" }
        .collectFile(name: "sample_read_counts.tsv", newLine: true, sort: true)

    ch_sourmash_gather_paths = ch_sourmash_gather_csv
        .map { _sample_id, _platform, _read_structure, gather_csv -> gather_csv }
        .collect()

    ch_sourmash_tax_summaries = ch_sourmash_tax_reports
        .map { sample_id, _platform, _read_structure, tax_reports ->
            def report_files = tax_reports instanceof List ? tax_reports : [tax_reports]
            def summarized_csv = report_files.find { report -> report.name.endsWith(".summarized.csv") }
            assert summarized_csv : "Missing sourmash tax summary CSV for ${sample_id}."
            summarized_csv
        }
        .collect()

    ch_blast_paths = ch_blast_final_tsv
        .map { _sample_id, blast_tsv -> blast_tsv }
        .collect()

    ch_crumbs_taxa_paths = ch_crumbs_taxa
        .map { _sample_id, taxa_tsv -> taxa_tsv }
        .collect()

    ch_crumbs_contig_paths = ch_crumbs_contigs
        .map { _sample_id, contigs_tsv -> contigs_tsv }
        .collect()

    BUILD_RAPID_SCREENING_EVAL_DB(
        ch_sample_read_counts_tsv,
        ch_sourmash_gather_paths,
        ch_sourmash_tax_summaries,
        ch_sourmash_lineages,
        ch_blast_paths,
        ch_crumbs_taxa_paths,
        ch_crumbs_contig_paths,
        run_id,
    )

    RENDER_RAPID_SCREENING_EVAL(
        BUILD_RAPID_SCREENING_EVAL_DB.out.database,
        BUILD_RAPID_SCREENING_EVAL_DB.out.followup_tsv,
        BUILD_RAPID_SCREENING_EVAL_DB.out.without_followup_tsv,
    )

    emit:
    database = BUILD_RAPID_SCREENING_EVAL_DB.out.database
    followup_tsv = BUILD_RAPID_SCREENING_EVAL_DB.out.followup_tsv
    without_followup_tsv = BUILD_RAPID_SCREENING_EVAL_DB.out.without_followup_tsv
    report = RENDER_RAPID_SCREENING_EVAL.out.report
}
