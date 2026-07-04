include { SUMMARIZE_CONTIG_COVERAGE } from "../modules/samtools"
include {
    COLLECT_CRUMBS_TAXIDS ;
    PREPARE_NCBI_PROFILE_TAXONOMY ;
    ESTIMATE_CRUMBS_PROFILE ;
    EXPORT_CRUMBS_TAXONOMIC_REPORTS ;
    RENDER_CRUMBS_TAXBURST
} from "../modules/crumbs"

workflow CRUMBS_PROFILING {
    take:
    ch_blast_results
    ch_filtered_bam
    ch_taxonomy_dir

    main:
    ch_requested_blast = ch_blast_results.filter { _sample_id, _blast_tsv ->
        params.experimental == true
    }

    ch_requested_bam = ch_filtered_bam.filter { _sample_id, _bam, _bai ->
        params.experimental == true
    }

    ch_taxid_inputs = ch_requested_blast
        .map { _sample_id, blast_tsv -> tuple("crumbs", blast_tsv) }
        .groupTuple()

    COLLECT_CRUMBS_TAXIDS(ch_taxid_inputs)

    PREPARE_NCBI_PROFILE_TAXONOMY(
        COLLECT_CRUMBS_TAXIDS.out.taxids,
        ch_taxonomy_dir,
    )

    SUMMARIZE_CONTIG_COVERAGE(ch_requested_bam)

    ch_profile_taxonomy = PREPARE_NCBI_PROFILE_TAXONOMY.out.profile_taxonomy
        .map { _profile_id, profile_taxonomy_tsv -> profile_taxonomy_tsv }

    ch_profile_inputs = ch_requested_blast
        .join(SUMMARIZE_CONTIG_COVERAGE.out.coverage_summary, by: 0)
        .combine(ch_profile_taxonomy)
        .map { sample_id, blast_tsv, coverage_tsv, profile_taxonomy_tsv ->
            tuple(sample_id, blast_tsv, coverage_tsv, profile_taxonomy_tsv)
        }

    ESTIMATE_CRUMBS_PROFILE(ch_profile_inputs)
    EXPORT_CRUMBS_TAXONOMIC_REPORTS(ESTIMATE_CRUMBS_PROFILE.out.taxa)
    RENDER_CRUMBS_TAXBURST(ESTIMATE_CRUMBS_PROFILE.out.taxa)

    emit:
    contigs = ESTIMATE_CRUMBS_PROFILE.out.contigs
    taxa = ESTIMATE_CRUMBS_PROFILE.out.taxa
    bioboxes_profile = ESTIMATE_CRUMBS_PROFILE.out.bioboxes_profile
    qc = ESTIMATE_CRUMBS_PROFILE.out.qc
    krona = EXPORT_CRUMBS_TAXONOMIC_REPORTS.out.krona
    kreport = EXPORT_CRUMBS_TAXONOMIC_REPORTS.out.kreport
    taxburst = RENDER_CRUMBS_TAXBURST.out.reports
    profile_taxonomy = PREPARE_NCBI_PROFILE_TAXONOMY.out.profile_taxonomy
}
