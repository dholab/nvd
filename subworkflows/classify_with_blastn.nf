include {
    MEGABLAST ;
    ANNOTATE_MEGABLAST_RESULTS ;
    FILTER_NON_VIRUS_MEGABLAST_NODES ;
    REMOVE_MEGABLAST_MAPPED_CONTIGS ;
    BLASTN_CLASSIFY ;
    ANNOTATE_BLASTN_RESULTS ;
    FILTER_NON_VIRUS_BLASTN_NODES ;
    EXTRACT_UNCLASSIFIED_CONTIGS ;
    MERGE_FILTERED_BLAST_RESULTS
} from "../modules/blast"

workflow CLASSIFY_WITH_BLASTN {
    take:
    ch_filtered_megablast
    ch_megablast_contigs
    ch_blast_db_files
    ch_gettax_sqlite_path

    main:
    BLASTN_CLASSIFY(
        ch_megablast_contigs.combine(ch_blast_db_files)
    )

    def nonEmptyBLASTNResults = BLASTN_CLASSIFY.out.filter { _id, hits_file ->
        hits_file.size() > 0 && hits_file.readLines().size() > 0
    }

    ANNOTATE_BLASTN_RESULTS(
        nonEmptyBLASTNResults.combine(ch_gettax_sqlite_path)
    )

    FILTER_NON_VIRUS_BLASTN_NODES(
        ANNOTATE_BLASTN_RESULTS.out
    )

    MERGE_FILTERED_BLAST_RESULTS(
        ch_filtered_megablast,
        FILTER_NON_VIRUS_BLASTN_NODES.out.filter { _id, hits_file ->
            hits_file.size() > 0 && hits_file.readLines().size() > 0
        }
    )

    // TODO: Eventually, does the pipeline need to support the scenario where either megablast
    // or blastn (though really just blastn) comes up with no hits

    // UPDATE 20250911: With the changes in the megablast subworkflow and the above, where
    // values in Nextflow queue channels are removed if the hits text files are empty, the
    // TODO above should be handled.

    // FIXME
    ch_merged_blast_results = MERGE_FILTERED_BLAST_RESULTS.out

    emit:
    merged_results = ch_merged_blast_results

}
