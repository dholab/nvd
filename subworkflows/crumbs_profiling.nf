include { SUMMARIZE_CONTIG_COVERAGE } from "../modules/samtools"
include {
    COLLECT_CRUMBS_TAXIDS ;
    PREPARE_NCBI_PROFILE_TAXONOMY ;
    ESTIMATE_CRUMBS_PROFILE ;
    EXPORT_CRUMBS_TAXONOMIC_REPORTS ;
    RENDER_CRUMBS_TAXBURST ;
    RENDER_MERGED_CRUMBS_TAXBURST
} from "../modules/crumbs"

workflow CRUMBS_PROFILING {
    take:
    ch_blast_results
    ch_filtered_bam
    ch_taxonomy_dir

    main:
    ch_taxid_inputs = ch_blast_results
        .map { _sample_id, blast_tsv -> tuple("crumbs", blast_tsv) }
        .groupTuple()

    COLLECT_CRUMBS_TAXIDS(ch_taxid_inputs)

    PREPARE_NCBI_PROFILE_TAXONOMY(
        COLLECT_CRUMBS_TAXIDS.out.taxids,
        ch_taxonomy_dir,
    )

    SUMMARIZE_CONTIG_COVERAGE(ch_filtered_bam)

    ch_profile_taxonomy = PREPARE_NCBI_PROFILE_TAXONOMY.out.profile_taxonomy
        .map { _profile_id, profile_taxonomy_tsv -> profile_taxonomy_tsv }

    // The estimator accepts zero or one coverage file. Existing mapback
    // samples contribute the singleton case; explicit no-contig samples add
    // the empty case when no-contig routing is activated.
    ch_coverage_inputs = SUMMARIZE_CONTIG_COVERAGE.out.coverage_summary.map { sample_id, coverage_tsv ->
        tuple(sample_id, [coverage_tsv])
    }

    ch_profile_inputs = ch_blast_results
        .join(ch_coverage_inputs, by: 0)
        .combine(ch_profile_taxonomy)
        .map { sample_id, blast_tsv, coverage_files, profile_taxonomy_tsv ->
            tuple(sample_id, blast_tsv, coverage_files, profile_taxonomy_tsv)
        }

    ESTIMATE_CRUMBS_PROFILE(ch_profile_inputs)
    EXPORT_CRUMBS_TAXONOMIC_REPORTS(ESTIMATE_CRUMBS_PROFILE.out.taxa)
    RENDER_CRUMBS_TAXBURST(ESTIMATE_CRUMBS_PROFILE.out.taxa)

    ch_merged_taxburst_input = ESTIMATE_CRUMBS_PROFILE.out.taxa
        .map { sample_id, taxa_tsv -> tuple("all", sample_id, taxa_tsv) }
        .groupTuple()
        .map { _key, sample_ids, taxa_tsvs -> tuple(sample_ids, taxa_tsvs) }

    RENDER_MERGED_CRUMBS_TAXBURST(ch_merged_taxburst_input)

    emit:
    queries = ESTIMATE_CRUMBS_PROFILE.out.queries
    taxa = ESTIMATE_CRUMBS_PROFILE.out.taxa
    bioboxes_profile = ESTIMATE_CRUMBS_PROFILE.out.bioboxes_profile
    qc = ESTIMATE_CRUMBS_PROFILE.out.qc
    krona = EXPORT_CRUMBS_TAXONOMIC_REPORTS.out.krona
    kreport = EXPORT_CRUMBS_TAXONOMIC_REPORTS.out.kreport
    taxburst = RENDER_CRUMBS_TAXBURST.out.reports
    merged_taxburst = RENDER_MERGED_CRUMBS_TAXBURST.out.report
    profile_taxonomy = PREPARE_NCBI_PROFILE_TAXONOMY.out.profile_taxonomy
}
