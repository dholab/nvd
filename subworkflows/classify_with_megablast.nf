include {
    MEGABLAST ;
    ANNOTATE_MEGABLAST_RESULTS ;
    FILTER_NON_VIRUS_MEGABLAST_NODES ;
    REMOVE_MEGABLAST_MAPPED_CONTIGS ;
    SELECT_TOP_BLAST_HITS
} from "../modules/blast"

workflow CLASSIFY_WITH_MEGABLAST {
    take:
    ch_virus_contigs
    ch_blast_db_files
    ch_gettax_sqlite_path

    main:
    MEGABLAST(
        ch_virus_contigs.combine(ch_blast_db_files)
    )

    def nonEmptyMegablastResults = MEGABLAST.out.filter { _id, hits_file ->
        hits_file.size() > 0 && hits_file.readLines().size() > 0
    }

    SELECT_TOP_BLAST_HITS(nonEmptyMegablastResults)

    ANNOTATE_MEGABLAST_RESULTS(
        SELECT_TOP_BLAST_HITS.out.combine(ch_gettax_sqlite_path)
    )

    // Capture this output for the LabKey table
    FILTER_NON_VIRUS_MEGABLAST_NODES(
        ANNOTATE_MEGABLAST_RESULTS.out.hits.filter { _id, hits_file ->
            hits_file.size() > 0 && hits_file.readLines().size() > 0
        }
    )

    REMOVE_MEGABLAST_MAPPED_CONTIGS(
        nonEmptyMegablastResults.join(ch_virus_contigs, by: 0)
    )

    emit:
    filtered_megablast = FILTER_NON_VIRUS_MEGABLAST_NODES.out
    megablast_contigs  = REMOVE_MEGABLAST_MAPPED_CONTIGS.out
    megablast          = SELECT_TOP_BLAST_HITS.out
}
