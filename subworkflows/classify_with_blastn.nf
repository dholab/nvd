include {
  MEGABLAST ;
  ANNOTATE_MEGABLAST_RESULTS ;
  FILTER_NON_VIRUS_MEGABLAST_NODES ;
  REMOVE_MEGABLAST_MAPPED_CONTIGS ;
  BLASTN_CLASSIFY ;
  ANNOTATE_BLASTN_RESULTS ;
  FILTER_NON_VIRUS_BLASTN_NODES ;
  COMBINE_BATCH_SEARCH_HITS ;
  STACK_BATCH_BLAST_RESULTS ;
  SELECT_TOP_BLAST_HITS
} from "../modules/blast"
include { ANNOTATE_LEAST_COMMON_ANCESTORS  } from "../modules/utils"

workflow CLASSIFY_WITH_BLASTN {
  take:
  ch_filtered_megablast
  ch_megablast_contigs
  ch_blast_db_files
  ch_taxonomy_dir // value channel: taxonomy directory path for taxonomy lookups

  main:
  // Keep one explicit channel emission per query-class batch that reaches
  // this subworkflow. BLASTN is still skipped when MEGABLAST leaves no pruned
  // query sequences, but the skipped case remains a channel emission instead
  // of disappearing.
  ch_samples_after_megablast = ch_filtered_megablast
    .join(ch_megablast_contigs, by: [0, 1])
    .map { sample_id, query_class, megablast_hits, classified, pruned_contigs ->
      def needs_blastn = file(pruned_contigs).size() > 0
      def meta = [id: sample_id, query_class: query_class, blastn_needed: needs_blastn]
      tuple(meta, megablast_hits, classified, pruned_contigs)
    }

  ch_samples_requiring_blastn = ch_samples_after_megablast.filter { meta, _megablast_hits, _classified, _pruned_contigs -> meta.blastn_needed }

  ch_blastn_context = ch_samples_requiring_blastn.map { meta, megablast_hits, _classified, _pruned_contigs ->
    tuple(meta.id, meta.query_class, meta, megablast_hits)
  }

  BLASTN_CLASSIFY(
    ch_samples_requiring_blastn.map { meta, _megablast_hits, classified, pruned_contigs ->
      tuple(meta.id, meta.query_class, classified, pruned_contigs)
    }.combine(ch_blast_db_files)
  )

  SELECT_TOP_BLAST_HITS(BLASTN_CLASSIFY.out)

  ANNOTATE_BLASTN_RESULTS(
    SELECT_TOP_BLAST_HITS.out,
    ch_taxonomy_dir,
  )

  FILTER_NON_VIRUS_BLASTN_NODES(
    ANNOTATE_BLASTN_RESULTS.out
  )

  ch_blastn_terminal = FILTER_NON_VIRUS_BLASTN_NODES.out
    .join(ch_blastn_context, by: [0, 1])
    .map { sample_id, query_class, blastn_hits, _meta, megablast_hits ->
      tuple(sample_id, query_class, [megablast_hits, blastn_hits])
    }

  ch_samples_skipping_blastn = ch_samples_after_megablast
    .filter { meta, _megablast_hits, _classified, _pruned_contigs -> !meta.blastn_needed }
    .map { meta, megablast_hits, _classified, _pruned_contigs ->
      tuple(meta.id, meta.query_class, [megablast_hits])
    }

  ch_merged_input = ch_blastn_terminal.mix(ch_samples_skipping_blastn)

  COMBINE_BATCH_SEARCH_HITS(ch_merged_input)

  ANNOTATE_LEAST_COMMON_ANCESTORS(COMBINE_BATCH_SEARCH_HITS.out, ch_taxonomy_dir)

  ch_sample_merged_inputs = ANNOTATE_LEAST_COMMON_ANCESTORS.out
    .map { sample_id, _query_class, batch_tsv -> tuple(sample_id, batch_tsv) }
    .groupTuple()

  STACK_BATCH_BLAST_RESULTS(ch_sample_merged_inputs)

  emit:
  merged_results = STACK_BATCH_BLAST_RESULTS.out
}
