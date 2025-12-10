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
    ch_gettax_sqlite_path

    main:
    BLASTN_CLASSIFY(
        ch_megablast_contigs.combine(ch_blast_db_files)
    )

    def nonEmptyBLASTNResults = BLASTN_CLASSIFY.out
        .filter { _id, hits_file ->
            hits_file.size() > 0 && hits_file.readLines().size() > 0
        }

    SELECT_TOP_BLAST_HITS(nonEmptyBLASTNResults)

    ANNOTATE_BLASTN_RESULTS(
        SELECT_TOP_BLAST_HITS.out.combine(ch_gettax_sqlite_path)
    )

    FILTER_NON_VIRUS_BLASTN_NODES(
        ANNOTATE_BLASTN_RESULTS.out
    )

    // Get blastn results (may not have all samples)
    ch_blastn_virus_hits = FILTER_NON_VIRUS_BLASTN_NODES.out
        .filter { _id, hits_file ->
            hits_file.size() > 0 && hits_file.readLines().size() > 0
        }

    // JOIN by sample_id with remainder:true to keep samples that only have megablast hits
ch_merged_input = ch_filtered_megablast
    .join(ch_blastn_virus_hits, by: 0, remainder: true)
    .map { sample_id, megablast_file, blastn_file ->
        if (blastn_file) {
            tuple(sample_id, megablast_file, blastn_file)
        } else {
            def emptyFile = file("${workflow.workDir}/empty_blastn_${sample_id}.tsv")
            emptyFile.text = "task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids\trank\n"
            tuple(sample_id, megablast_file, emptyFile)
        }
    }

    MERGE_FILTERED_BLAST_RESULTS(ch_merged_input)

    ANNOTATE_LEAST_COMMON_ANCESTORS(MERGE_FILTERED_BLAST_RESULTS.out)

    emit:
    merged_results = ANNOTATE_LEAST_COMMON_ANCESTORS.out
}
