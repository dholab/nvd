include {
    MEGABLAST ;
    ANNOTATE_MEGABLAST_RESULTS ;
    FILTER_NON_VIRUS_MEGABLAST_NODES ;
    REMOVE_MEGABLAST_MAPPED_CONTIGS ;
    BLASTN_CLASSIFY ;
    ANNOTATE_BLASTN_RESULTS ;
    FILTER_NON_VIRUS_BLASTN_NODES ;
    MERGE_FILTERED_BLAST_RESULTS ;
    SELECT_TOP_BLAST_HITS
} from "../modules/blast"
include { ANNOTATE_LEAST_COMMON_ANCESTORS  } from "../modules/utils"

workflow CLASSIFY_WITH_BLASTN {
    take:
    ch_filtered_megablast
    ch_megablast_contigs
    ch_blast_db_files
    ch_taxonomy_dir   // value channel: taxonomy directory path for taxonomy lookups

    main:
    // The pruned FASTA is an uncompressed pipeline output. Empty FASTA producers
    // create zero-byte files, so this avoids scheduling BLASTN when MEGABLAST
    // left no query sequences. Do not use this predicate for gzipped FASTA inputs.
    ch_blastn_candidates = ch_megablast_contigs.filter { _sample_id, _classified, pruned_contigs ->
        file(pruned_contigs).size() > 0
    }

    BLASTN_CLASSIFY(
        ch_blastn_candidates.combine(ch_blast_db_files)
    )

    SELECT_TOP_BLAST_HITS(BLASTN_CLASSIFY.out)

    ANNOTATE_BLASTN_RESULTS(
        SELECT_TOP_BLAST_HITS.out,
        ch_taxonomy_dir
    )

    FILTER_NON_VIRUS_BLASTN_NODES(
        ANNOTATE_BLASTN_RESULTS.out
    )

    ch_merged_input = ch_filtered_megablast
        .mix(FILTER_NON_VIRUS_BLASTN_NODES.out)
        .filter { _sample_id, hits -> hasDataRows(hits) }
        .groupTuple(by: 0)

    MERGE_FILTERED_BLAST_RESULTS(ch_merged_input)

    ANNOTATE_LEAST_COMMON_ANCESTORS(MERGE_FILTERED_BLAST_RESULTS.out, ch_taxonomy_dir)

    emit:
    merged_results = ANNOTATE_LEAST_COMMON_ANCESTORS.out
}


def hasDataRows(hits) {
    def hits_file = file(hits)
    if (!hits_file.exists() || hits_file.size() == 0) {
        return false
    }
    hits_file.withReader { reader ->
        reader.readLine() != null && reader.readLine() != null
    }
}
