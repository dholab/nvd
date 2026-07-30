process BUILD_RAPID_SCREENING_EVAL_DB {

    label "low"

    input:
    path sample_read_counts_tsv
    path sourmash_gather_csvs
    path sourmash_tax_summaries
    path lineages_csv
    path blast_tsvs
    path crumbs_taxa_tsvs
    path crumbs_contigs_tsvs
    val run_id

    output:
    path "rapid_screening_eval.duckdb", emit: database
    path "exports/screening_signal_followup_by_sample_rank.tsv", emit: followup_tsv
    path "exports/screening_signals_without_same_rank_followup.tsv", emit: without_followup_tsv

    script:
    def asList = { value -> value instanceof List ? value : [value] }
    def pathArgs = { flag, paths ->
        asList(paths).findAll { it }.collect { path -> "${flag} '${path}'" }.join(" ")
    }
    def sourmashGatherArgs = pathArgs("--sourmash-gather-csv", sourmash_gather_csvs)
    def sourmashTaxArgs = pathArgs("--sourmash-tax-summary", sourmash_tax_summaries)
    def blastArgs = pathArgs("--blast-tsv", blast_tsvs)
    def crumbsTaxaArgs = pathArgs("--crumbs-taxa-tsv", crumbs_taxa_tsvs)
    def crumbsContigsArgs = pathArgs("--crumbs-contigs-tsv", crumbs_contigs_tsvs)
    """
    evaluate_rapid_screening.py \
        --sample-read-counts ${sample_read_counts_tsv} \
        --lineages-csv ${lineages_csv} \
        ${sourmashGatherArgs} \
        ${sourmashTaxArgs} \
        ${blastArgs} \
        ${crumbsTaxaArgs} \
        ${crumbsContigsArgs} \
        --output-database rapid_screening_eval.duckdb \
        --output-dir .
    """
}

process RENDER_RAPID_SCREENING_EVAL {

    label "low"

    input:
    path database
    path followup_tsv
    path without_followup_tsv

    output:
    path "rapid_screening_eval.html", emit: report

    script:
    """
    render_rapid_screening_eval.py \
        --database ${database} \
        --followup-tsv ${followup_tsv} \
        --without-followup-tsv ${without_followup_tsv} \
        --output rapid_screening_eval.html
    """
}
