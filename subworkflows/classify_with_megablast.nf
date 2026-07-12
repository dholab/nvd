include {
    MEGABLAST ;
    ANNOTATE_MEGABLAST_RESULTS ;
    FILTER_NON_VIRUS_MEGABLAST_NODES ;
    PARTITION_MEGABLAST_QUERIES ;
    SELECT_TOP_BLAST_HITS
} from "../modules/blast"

workflow CLASSIFY_WITH_MEGABLAST {
    take:
    ch_query_batches  // tuple(sample_id, platform, query_class, fasta, lookup)
    ch_blast_db_files
    ch_taxonomy_dir   // value channel: taxonomy directory path for taxonomy lookups

    main:
    // These FASTA files are uncompressed pipeline outputs. Upstream empty FASTA
    // producers create zero-byte files, so a byte-size check is a cheap scheduling
    // guard. Do not use this predicate for gzipped FASTA inputs.
    ch_megablast_candidates = ch_query_batches.filter { _sample_id, _platform, _query_class, _query_fasta, _lookup ->
          !params.skip_blast
      }
      .filter { _sample_id, _platform, _query_class, query_fasta, _lookup ->
          file(query_fasta).size() > 0
      }
      .multiMap { sample_id, _platform, query_class, query_fasta, _lookup ->
          for_search: tuple(sample_id, query_class, query_fasta)
          for_partition: tuple(sample_id, query_class, query_fasta)
      }

    MEGABLAST(
        ch_megablast_candidates.for_search.combine(ch_blast_db_files)
    )

    SELECT_TOP_BLAST_HITS(MEGABLAST.out)

    ANNOTATE_MEGABLAST_RESULTS(
        SELECT_TOP_BLAST_HITS.out,
        ch_taxonomy_dir
    )

    // Capture this output for the LabKey table
    FILTER_NON_VIRUS_MEGABLAST_NODES(
        ANNOTATE_MEGABLAST_RESULTS.out.hits
    )

    PARTITION_MEGABLAST_QUERIES(
        MEGABLAST.out.join(ch_megablast_candidates.for_partition, by: [0, 1])
    )

    emit:
    filtered_megablast = FILTER_NON_VIRUS_MEGABLAST_NODES.out.hits
    filter_decisions = FILTER_NON_VIRUS_MEGABLAST_NODES.out.decisions
    megablast_query_partition = PARTITION_MEGABLAST_QUERIES.out
    megablast          = SELECT_TOP_BLAST_HITS.out
}
