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
    ch_state_dir  // value channel: state directory path for taxonomy lookups

    main:
    // Static empty blastn template file - using a stable file prevents cache invalidation.
    // Previously, dynamically creating files in the .map block caused MERGE_FILTERED_BLAST_RESULTS
    // to never resume because the file got a new timestamp on each run.
    def emptyBlastnTemplate = file("${projectDir}/assets/empty_blastn.tsv")

    BLASTN_CLASSIFY(
        ch_megablast_contigs.combine(ch_blast_db_files)
    )

    def nonEmptyBLASTNResults = BLASTN_CLASSIFY.out
        .filter { _id, hits_file ->
            hits_file.size() > 0 && hits_file.readLines().size() > 0
        }

    SELECT_TOP_BLAST_HITS(nonEmptyBLASTNResults)

    ANNOTATE_BLASTN_RESULTS(
        SELECT_TOP_BLAST_HITS.out,
        ch_state_dir
    )

    FILTER_NON_VIRUS_BLASTN_NODES(
        ANNOTATE_BLASTN_RESULTS.out
    )

    // Get blastn results (may not have all samples)
    ch_blastn_virus_hits = FILTER_NON_VIRUS_BLASTN_NODES.out
        .filter { _id, hits_file ->
            hits_file.size() > 0 && hits_file.readLines().size() > 0
        }

    // JOIN by sample_id with remainder:true to keep samples that only have megablast hits.
    // For samples without blastn hits, use the static empty template file instead of
    // dynamically creating one (which would break caching).
    ch_merged_input = ch_filtered_megablast
        .join(ch_blastn_virus_hits, by: 0, remainder: true)
        .map { sample_id, megablast_file, blastn_file ->
            tuple(sample_id, megablast_file, blastn_file ?: emptyBlastnTemplate)
        }

    MERGE_FILTERED_BLAST_RESULTS(ch_merged_input)

    ANNOTATE_LEAST_COMMON_ANCESTORS(MERGE_FILTERED_BLAST_RESULTS.out, ch_state_dir)

    emit:
    merged_results = ANNOTATE_LEAST_COMMON_ANCESTORS.out
}
