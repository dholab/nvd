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

    // TODO: Eventually, does the pipeline need to support the scenario where either megablast
    // or blastn (though really just blastn) comes up with no hits?
    ch_merged_blast_results = FILTER_NON_VIRUS_BLASTN_NODES.out
        .join(ch_megablast, by: 0)
        .collectFile { sample_id, _virus_only_txt, blast_file ->
            ["${sample_id}.txt", file(blast_file).text + '\n']
        }

    EXTRACT_UNCLASSIFIED_CONTIGS(
        ch_unmapped_contigs
        .join(ch_megablast, by: 0)
        .join(BLASTN_CLASSIFY.out, by: 0)
    )

    emit:
    ch_merged_blast_results

}
