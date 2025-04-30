include {
    MEGABLAST ;
    ANNOTATE_MEGABLAST_RESULTS ;
    FILTER_NON_VIRUS_MEGABLAST_NODES ;
    REMOVE_MEGABLAST_MAPPED_CONTIGS ;
    BLASTN_CLASSIFY ;
    ANNOTATE_BLASTN_RESULTS ;
    FILTER_NON_VIRUS_BLASTN_NODES ;
    EXTRACT_UNCLASSIFIED_CONTIGS
} from "../modules/blast"

workflow CLASSIFY_WITH_BLASTN {
    take:
    ch_unmapped_contigs
    ch_megablast
    ch_blast_db_files
    ch_gettax_sqlite_path

    main:
    BLASTN_CLASSIFY(
        ch_unmapped_contigs.combine(ch_blast_db_files)
    )

    ANNOTATE_BLASTN_RESULTS(
        BLASTN_CLASSIFY.out.combine(ch_gettax_sqlite_path)
    )

    FILTER_NON_VIRUS_BLASTN_NODES(
        ANNOTATE_BLASTN_RESULTS.out
    )

    ch_merged_blast_results = FILTER_NON_VIRUS_MEGABLAST_NODES.out
        .mix(
            FILTER_NON_VIRUS_BLASTN_NODES
        )
        .groupTuple(by: 0)
        .collectFile { sample_id, blast_file ->
            ["${sample_id}.txt", file(blast_file).text + '\n']
        }

    EXTRACT_UNCLASSIFIED_CONTIGS(
        ch_unmapped_contigs.mix(
            ch_megablast,
            BLASTN_CLASSIFY.out,
        ).groupTuple(by: 0)
    )

    emit:
    ch_merged_blast_results

}
