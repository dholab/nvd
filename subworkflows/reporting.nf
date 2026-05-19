include { ADD_READ_COUNTS_TO_BLAST; CONCATENATE_EXPERIMENT_BLAST } from "../modules/utils"

workflow REPORTING {
    take:
    ch_blast_results
    ch_read_counts

    main:
    // Add total read counts to the final merged BLAST results so they are
    // always present in the published TSV, regardless of LabKey.
    ch_blast_with_counts = ch_blast_results
        .join(ch_read_counts, by: 0)

    ADD_READ_COUNTS_TO_BLAST(ch_blast_with_counts)

    // Concatenate all per-sample final BLAST results into a single experiment-level TSV.
    // Runs unconditionally so every run produces an experiment summary, not just LabKey runs.
    CONCATENATE_EXPERIMENT_BLAST(
        ADD_READ_COUNTS_TO_BLAST.out.map { _sample_id, tsv -> tsv }.collect()
    )

    emit:
    blast_results = ADD_READ_COUNTS_TO_BLAST.out
    experiment_blast = CONCATENATE_EXPERIMENT_BLAST.out.concatenated_tsv
}
