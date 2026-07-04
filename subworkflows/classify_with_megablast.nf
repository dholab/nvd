include {
    MEGABLAST ;
    ANNOTATE_MEGABLAST_RESULTS ;
    FILTER_NON_VIRUS_MEGABLAST_NODES ;
    REMOVE_MEGABLAST_MAPPED_CONTIGS ;
    SELECT_TOP_BLAST_HITS
} from "../modules/blast"

workflow CLASSIFY_WITH_MEGABLAST {
    take:
    ch_virus_contigs  // tuple(sample_id, platform, read_structure, fasta, lookup)
    ch_blast_db_files
    ch_taxonomy_dir   // value channel: taxonomy directory path for taxonomy lookups

    main:
    // These FASTA files are uncompressed pipeline outputs. Upstream empty FASTA
    // producers create zero-byte files, so a byte-size check is a cheap scheduling
    // guard. Do not use this predicate for gzipped FASTA inputs.
    ch_megablast_candidates = ch_virus_contigs.filter { _sample_id, _platform, _read_structure, _contigs, _lookup ->
          !params.skip_blast
      }
      .filter { _sample_id, _platform, _read_structure, contigs, _lookup ->
          file(contigs).size() > 0
      }
      .multiMap { sample_id, _platform, _read_structure, contigs, _lookup ->
          for_search: tuple(sample_id, contigs)
          for_prune:  tuple(sample_id, contigs)
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

    REMOVE_MEGABLAST_MAPPED_CONTIGS(
        MEGABLAST.out.join(ch_megablast_candidates.for_prune, by: 0)
    )

    emit:
    filtered_megablast = FILTER_NON_VIRUS_MEGABLAST_NODES.out
    megablast_contigs  = REMOVE_MEGABLAST_MAPPED_CONTIGS.out
    megablast          = SELECT_TOP_BLAST_HITS.out
}
