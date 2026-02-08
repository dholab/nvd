include {
    MEGABLAST ;
    ANNOTATE_MEGABLAST_RESULTS ;
    FILTER_NON_VIRUS_MEGABLAST_NODES ;
    REMOVE_MEGABLAST_MAPPED_CONTIGS ;
    BLASTN_CLASSIFY ;
    ANNOTATE_BLASTN_RESULTS ;
    FILTER_NON_VIRUS_BLASTN_NODES ;
    EXTRACT_UNCLASSIFIED_CONTIGS ;
    MERGE_FILTERED_BLAST_RESULTS ;
    SELECT_TOP_BLAST_HITS
} from "../modules/blast"
include { ANNOTATE_LEAST_COMMON_ANCESTORS  } from "../modules/utils"

workflow CLASSIFY_WITH_BLASTN {
    take:
    ch_filtered_megablast
    ch_megablast_contigs
    ch_blast_db_files
    ch_state_dir      // value channel: state directory path (may be null in stateless mode)
    ch_taxonomy_dir   // value channel: taxonomy directory path for taxonomy lookups

    main:
    BLASTN_CLASSIFY(
        ch_megablast_contigs.combine(ch_blast_db_files)
    )

    SELECT_TOP_BLAST_HITS(BLASTN_CLASSIFY.out)

    ANNOTATE_BLASTN_RESULTS(
        SELECT_TOP_BLAST_HITS.out,
        ch_state_dir,
        ch_taxonomy_dir
    )

    FILTER_NON_VIRUS_BLASTN_NODES(
        ANNOTATE_BLASTN_RESULTS.out
    )

    ch_merged_input = ch_filtered_megablast
        .join(FILTER_NON_VIRUS_BLASTN_NODES.out, by: 0)

    MERGE_FILTERED_BLAST_RESULTS(ch_merged_input)

    ANNOTATE_LEAST_COMMON_ANCESTORS(MERGE_FILTERED_BLAST_RESULTS.out, ch_state_dir, ch_taxonomy_dir)

    emit:
    merged_results = ANNOTATE_LEAST_COMMON_ANCESTORS.out
}
