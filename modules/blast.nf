process SELECT_BLAST_QUERIES {
  /* Select screened records into BLAST query batches by query class. */

  tag "${sample_id}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(fasta), path(lookup)

  output:
  tuple val(sample_id), val(platform), path("*.blast_queries.fasta"), path(lookup), emit: query_fastas, optional: true

  script:
  """
  select_blast_queries.py \
      --sample-id '${sample_id}' \
      --input-fasta ${fasta} \
      --query-lookup ${lookup} \
      --output-dir .
  """
}

process NORMALIZE_READ_BLAST_QUERIES {
  /* Normalize mapback-unmapped reads into exact-deduplicated BLAST query batches. */

  tag "${sample_id}, ${query_class}"
  label "medium"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(platform), val(query_class), path(reads)

  output:
  tuple val(sample_id), val(platform), val(query_class), path("*.blast_queries.fasta"), path("*.query_sequences.sqlite"), emit: queries, optional: true

  script:
  def producer = query_class == "overlap_merged_pair" ? "bbmerge" : "source_read"
  def support_count_policy = params.dedup || params.dedup_seq ? "one" : "abundance"
  def read_files = reads instanceof List ? reads : [reads]
  def fastq_args = read_files.collect { read_file -> "--fastq ${read_file}" }.join(" ")
  """
  normalize_read_blast_queries.py \
      --sample-id '${sample_id}' \
      --query-class '${query_class}' \
      --producer '${producer}' \
      ${fastq_args} \
      --support-count-policy '${support_count_policy}' \
      --output-dir .
  """
}

process SUMMARIZE_BLAST_QUERY_BATCHES {
  /* Summarize prepared BLAST query batches, including absent query classes. */

  tag "${sample_id}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(platform), val(query_classes), path(query_fastas, stageAs: "query_fastas??????/*"), path(query_lookups, stageAs: "query_lookups??????/*")

  output:
  tuple val(sample_id), path("${sample_id}.blast_query_batches.tsv")

  script:
  def classes = query_classes instanceof List ? query_classes : [query_classes]
  def fastas = query_fastas instanceof List ? query_fastas : [query_fastas]
  def lookups = query_lookups instanceof List ? query_lookups : [query_lookups]
  assert classes.size() == fastas.size()
  assert classes.size() == lookups.size()
  def batch_args = classes.withIndex().collect { query_class, index ->
    "--batch '${query_class}' ${fastas[index]} ${lookups[index]}"
  }.join(" ")
  """
  summarize_blast_query_batches.py \
      --sample-id '${sample_id}' \
      --platform '${platform}' \
      ${batch_args} \
      --output ${sample_id}.blast_query_batches.tsv
  """
}

process SUMMARIZE_EMPTY_BLAST_QUERY_BATCHES {
  /* Emit the same four-row decision table for samples with no query batches. */

  tag "${sample_id}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(platform)

  output:
  tuple val(sample_id), path("${sample_id}.blast_query_batches.tsv")

  script:
  """
  summarize_blast_query_batches.py \
    --sample-id '${sample_id}' \
    --platform '${platform}' \
    --output ${sample_id}.blast_query_batches.tsv
  """
}

/*
Perform rapid sequence similarity search using MEGABLAST.

This rule aligns query sequences against a comprehensive nucleotide database
to identify similar known sequences. It's crucial for initial, broad_scale
identification of viral sequences.
*/
process MEGABLAST {
  tag "${sample_id}, ${query_class}"
  label "high"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(query_fasta), path(blast_db)

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.megablast.txt")

  script:
  def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
  def blast_task = "megablast"
  def index_exists = file("${blast_db}/${params.blast_db_prefix}.00.idx").exists() || file("${blast_db}/${params.blast_db_prefix}.00.00.idx").exists()
  def use_index = blast_task == "megablast" && index_exists ? "-use_index true" : ""
  """
  export BLASTDB=${blast_db}
  blastn -task ${blast_task} \
  -db ${blast_db}/${params.blast_db_prefix} \
  -query ${query_fasta} \
  -num_threads ${task.cpus} \
  -outfmt "${outfmt}" \
  -max_target_seqs ${params.max_blast_targets} \
  ${use_index} \
  > ${sample_id}.${query_class}.megablast.txt
  """
}

//Create file with megablast results annotated with megablast task and full taxonomic rank.
process ANNOTATE_MEGABLAST_RESULTS {

  tag "${sample_id}, ${query_class}"
  label "medium"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(blast_txt)
  val taxonomy_dir

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.annotated_megablast.txt"), emit: hits

  script:
  def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
  def taxonomy_mode_arg = params.taxonomy_mode ? "--taxonomy-mode '${params.taxonomy_mode}'" : ""
  def taxonomy_max_age_arg = params.taxonomy_max_age_days ? "--taxonomy-max-age-days ${params.taxonomy_max_age_days}" : ""
  """
  annotate_blast_results.py \
  --input_file ${blast_txt} \
  --output_file ${sample_id}.${query_class}.annotated_megablast.txt \
  --sample_name ${sample_id} \
  --task 'megablast' \
  ${taxonomy_dir_arg} \
  ${taxonomy_mode_arg} \
  ${taxonomy_max_age_arg}
  """
}

process SELECT_TOP_BLAST_HITS {

  tag "${sample_id}, ${query_class}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(blast_txt)

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.top${params.blast_retention_count}_hits.txt")

  script:
  """
    select_top_blast_hits.py \\
    --input-file ${blast_txt} \\
    --output-file ${sample_id}.${query_class}.top${params.blast_retention_count}_hits.txt \\
    --blast-retention-count ${params.blast_retention_count}
    """
}
/*

Remove any query groups that do not have at least one hit corresponding to viruses.

This rule handles cases where the input file is empty and filters out
blast query read groups where none of the entries are from superkingdom: Viruses
or where the only viral hits are phages. This will leave in any query seq that has things that are 
non virus as long as one of the qseqid reads has a viral hit.
*/
process FILTER_NON_VIRUS_MEGABLAST_NODES {

  tag "${sample_id}, ${query_class}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(annotated_blast)

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.mb_virus_only.txt"), emit: hits
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.megablast_query_filtering.tsv"), emit: decisions

  script:
  """
    filter_non_virus_blast_nodes.py ${annotated_blast} ${sample_id}.${query_class}.mb_virus_only.txt
    summarize_blast_filtering.py \
      --sample-id '${sample_id}' \
      --query-class '${query_class}' \
      --stage megablast_virus_filter \
      --input ${annotated_blast} \
      --output ${sample_id}.${query_class}.mb_virus_only.txt \
      --summary ${sample_id}.${query_class}.megablast_query_filtering.tsv
    """
}

/*
Partition queries according to whether MEGABLAST produced any hit.

This rule generates two outputs:
1. Query IDs accounted for by MEGABLAST
2. A FASTA containing residual queries that are candidates for BLASTN

The candidate FASTA is used as input for more sensitive BLASTN searching.
*/
process PARTITION_MEGABLAST_QUERIES {

  tag "${sample_id}, ${query_class}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(megablast_results), path(query_fasta)

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.megablast_accounted_query_ids.txt"), path("${sample_id}.${query_class}.blastn_candidate_queries.fasta"), path("${sample_id}.${query_class}.megablast_query_partition.tsv")

  script:
  """
    partition_megablast_queries.py \
    --megablast-results ${megablast_results} \
    --query-fasta ${query_fasta} \
    --sample-id '${sample_id}' \
    --query-class '${query_class}' \
    --accounted-query-ids ${sample_id}.${query_class}.megablast_accounted_query_ids.txt \
    --blastn-candidate-fasta ${sample_id}.${query_class}.blastn_candidate_queries.fasta \
    --summary-tsv ${sample_id}.${query_class}.megablast_query_partition.tsv
    """
}

/*
Perform sequence similarity search using BLASTN.

This rule aligns residual query sequences against a comprehensive nucleotide database
to identify similar known sequences with BLASTN
*/
process BLASTN_CLASSIFY {
  tag "${sample_id}, ${query_class}"
  label "high"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(accounted_query_ids), path(blastn_candidate_fasta), path(blast_db)

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.blastn.txt")

  script:
  def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
  def blast_task = "blastn"
  """
    export BLASTDb=${blast_db}/
    blastn -task ${blast_task} \
    -db ${blast_db}/${params.blast_db_prefix} \
    -query ${blastn_candidate_fasta} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${params.max_blast_targets} \
    > ${sample_id}.${query_class}.blastn.txt
    """
}

// Create file with blastn results annotated with blastn task and full taxonomic rank.
process ANNOTATE_BLASTN_RESULTS {

  tag "${sample_id}, ${query_class}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(blastn_txt)
  val taxonomy_dir

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.annotated_blastn.txt")

  script:
  def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
  def taxonomy_mode_arg = params.taxonomy_mode ? "--taxonomy-mode '${params.taxonomy_mode}'" : ""
  def taxonomy_max_age_arg = params.taxonomy_max_age_days ? "--taxonomy-max-age-days ${params.taxonomy_max_age_days}" : ""
  """
    annotate_blast_results.py \
    --sample_name ${sample_id} \
    --input_file ${blastn_txt} \
    --output_file ${sample_id}.${query_class}.annotated_blastn.txt \
    --task 'blastn' \
    ${taxonomy_dir_arg} \
    ${taxonomy_mode_arg} \
    ${taxonomy_max_age_arg}
    """
}

/*
Remove any query groups that do not have at least one hit corresponding to viruses.

This rule handles cases where the input file is empty and filters out
groups where none of the entries are from superkingdom: Viruses
or where the only viral hits are phages.
*/
process FILTER_NON_VIRUS_BLASTN_NODES {

  tag "${sample_id}, ${query_class}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(annotated_blast)

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.nt_virus_only.txt"), emit: hits
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.blastn_query_filtering.tsv"), emit: decisions

  script:
  """
  filter_non_virus_blast_nodes.py ${annotated_blast} ${sample_id}.${query_class}.nt_virus_only.txt
  summarize_blast_filtering.py \
    --sample-id '${sample_id}' \
    --query-class '${query_class}' \
    --stage blastn_virus_filter \
    --input ${annotated_blast} \
    --output ${sample_id}.${query_class}.nt_virus_only.txt \
    --summary ${sample_id}.${query_class}.blastn_query_filtering.tsv
  """
}

// Combine MEGABLAST and optional BLASTN search hits for one query-class batch.
process COMBINE_BATCH_SEARCH_HITS {

  tag "${sample_id}, ${query_class}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(blast_hit_files)

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.top_blast_hits.tsv")

  script:
  def hit_args = blast_hit_files
    .sort { it.name }
    .collect { "--blast-hits ${it}" }
    .join(" ")
  """
  concat_multiblast_hits.py \\
      ${hit_args} \\
      --output-file ${sample_id}.${query_class}.top_blast_hits.tsv
  """
}

// Stack LCA-annotated query-class batch BLAST results for one sample.
process STACK_BATCH_BLAST_RESULTS {

  tag "${sample_id}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), path(blast_hit_files)

  output:
  tuple val(sample_id), path("${sample_id}_blast.merged_with_lca.tsv")

  script:
  def hit_args = blast_hit_files
    .sort { it.name }
    .collect { "--blast-hits ${it}" }
    .join(" ")
  """
  concat_multiblast_hits.py \\
      ${hit_args} \\
      --output-file ${sample_id}_blast.merged_with_lca.tsv
  """
}
