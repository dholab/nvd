/*
 * Record the invocation and source identity that produced the result tree.
 */
process RECORD_RUN_METADATA {

  label "low"
  cache false

  input:
  val command
  val nvd_version
  val source_revision
  val source_dirty

  output:
  path "run_command.sh", emit: command_file
  path "nvd_version.txt", emit: version_file

  script:
  def encoded_command = command.bytes.encodeBase64().toString()
  """
  write_run_metadata.py \
      --command-base64 '${encoded_command}' \
      --version '${nvd_version}' \
      --revision '${source_revision}' \
      --dirty '${source_dirty}'
  """
}

/*
 * Compute the deterministic run context without writing workflow state.
 */
process COMPUTE_RUN_CONTEXT {

  label "low"
  cache false

  input:
  path samplesheet

  output:
  val true, emit: ready
  stdout emit: run_context

  script:
  """
  compute_run_context.py \\
      --verbose \\
      --samplesheet ${samplesheet}
  """
}

/*
 * Ensure the NCBI taxonomy cache is ready before distributed annotation workers run.
 */
process ENSURE_TAXONOMY {

  label "low"
  cache false

  input:
  val taxonomy_dir

  output:
  stdout emit: taxonomy_dir

  script:
  def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
  def taxonomy_mode_arg = params.taxonomy_mode ? "--taxonomy-mode '${params.taxonomy_mode}'" : ""
  def taxonomy_refresh_arg = params.taxonomy_refresh ? "--taxonomy-refresh '${params.taxonomy_refresh}'" : ""
  def taxonomy_max_age_arg = params.taxonomy_max_age_days ? "--taxonomy-max-age-days ${params.taxonomy_max_age_days}" : ""
  """
  ensure_taxonomy.py \\
      --verbose \\
      ${taxonomy_dir_arg} \\
      ${taxonomy_mode_arg} \\
      ${taxonomy_refresh_arg} \\
      ${taxonomy_max_age_arg}
  """
}

process ANNOTATE_LEAST_COMMON_ANCESTORS {

  tag "${sample_id}, ${query_class}"
  label "low"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  input:
  tuple val(sample_id), val(query_class), path(all_blast_hits)
  val taxonomy_dir

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.blast.merged_with_lca.tsv")

  script:
  def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
  def taxonomy_mode_arg = params.taxonomy_mode ? "--taxonomy-mode '${params.taxonomy_mode}'" : ""
  def taxonomy_max_age_arg = params.taxonomy_max_age_days ? "--taxonomy-max-age-days ${params.taxonomy_max_age_days}" : ""
  """
  annotate_blast_lca.py -i ${all_blast_hits} -o ${sample_id}.${query_class}.blast.merged_with_lca.tsv ${taxonomy_dir_arg} ${taxonomy_mode_arg} ${taxonomy_max_age_arg}
  """
}

process ADD_READ_COUNTS_TO_BLAST {
  /*
   * Enrich merged BLAST results with pipeline metadata columns:
   *   query_class / producer / source_id — collected query sequence metadata
   *   mapped_reads  — per-contig read count (joined by qseqid)
   *   total_reads   — sample-level total input reads
   *   blast_db_version — BLAST database version
   *   virus_index_version — STAT k-mer database version
   *   nextflow_run_id — workflow run identifier
   *
   * This consolidates all metadata enrichment into one step so that
   * the published _blast.final.tsv is complete regardless of whether
   * LabKey is enabled. LIMS integration only needs to add experiment_id.
   */

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), path(blast_tsv), val(total_reads), path(contig_count_files), path(query_lookups)
  val run_id
  val blast_db_version
  val virus_index_version

  output:
  tuple val(sample_id), path("${sample_id}_blast.final.tsv")

  script:
  def count_files = contig_count_files instanceof List ? contig_count_files : [contig_count_files]
  assert count_files.size() <= 1 : "Expected at most one mapped-count file for ${sample_id}, received ${count_files.size()}"
  def count_arg = count_files ? "--contig-counts ${count_files[0]}" : ""
  def lookup_files = query_lookups instanceof List ? query_lookups : [query_lookups]
  def lookup_args = lookup_files.collect { lookup -> "--query-lookup ${lookup}" }.join(" ")
  """
  finalize_blast_results.py \\
      --blast-tsv ${blast_tsv} \\
      ${count_arg} \\
      ${lookup_args} \\
      --output ${sample_id}_blast.final.tsv \\
      --total-reads ${total_reads} \\
      --blast-db-version '${blast_db_version}' \\
      --virus-index-version '${virus_index_version}' \\
      --run-id '${run_id}'
  """
}

process BUILD_QUERY_BIG_TABLE {
  /*
   * Build the featured one-row-per-query Big Table artifact from final BLAST rows.
   */

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), path(blast_tsv), path(crumbs_tsv)

  output:
  tuple val(sample_id), path("${sample_id}.query_big_table.tsv")

  script:
  """
  build_query_big_table.py \
      --blast-tsv ${blast_tsv} \
      --crumbs-tsv ${crumbs_tsv} \
      --output ${sample_id}.query_big_table.tsv
  """
}

process BUILD_TAXON_BIG_TABLE {
  /*
   * Build the per-sample one-row-per-taxon Big Table artifact from query Big Tables.
   */

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), path(query_big_table), path(crumbs_taxa_tsv)

  output:
  tuple val(sample_id), path("${sample_id}.taxon_big_table.tsv")

  script:
  """
  build_taxon_big_table.py \
      --query-big-table ${query_big_table} \
      --crumbs-taxa-tsv ${crumbs_taxa_tsv} \
      --output ${sample_id}.taxon_big_table.tsv
  """
}

process CONCATENATE_EXPERIMENT_BLAST {
  /*
   * Concatenate all per-sample _blast.final.tsv files into a single
   * experiment-level TSV. Runs unconditionally (does not require LabKey).
   * Header is taken from the first file; subsequent files contribute
   * data rows only.
   */

  label "low"

  input:
  path "sample_results/*"

  output:
  path "experiment_blast_results.tsv", emit: concatenated_tsv

  script:
  """
  #!/usr/bin/env python3
  from pathlib import Path

  files = sorted(Path("sample_results").glob("*.tsv"))

  if not files:
      Path("experiment_blast_results.tsv").write_text("")
  else:
      with open("experiment_blast_results.tsv", "w") as out:
          # Write header from first file
          with open(files[0]) as f:
              header = f.readline()
              out.write(header)
              for line in f:
                  out.write(line)
          # Append data rows from remaining files
          for tsv in files[1:]:
              with open(tsv) as f:
                  f.readline()  # skip header
                  for line in f:
                      out.write(line)
  """
}

process CONCATENATE_QUERY_BIG_TABLE {
  /*
   * Concatenate per-sample query Big Tables into the featured all-sample artifact.
   */

  label "low"

  input:
  path "sample_big_tables/*"

  output:
  path "query_big_table.tsv", emit: concatenated_tsv

  script:
  """
  stack_big_tables.py \
      --input-dir sample_big_tables \
      --output query_big_table.tsv
  """
}

process CONCATENATE_TAXON_BIG_TABLE {
  /*
   * Concatenate per-sample taxon Big Tables into the featured all-sample artifact.
   */

  label "low"

  input:
  path "sample_big_tables/*"

  output:
  path "taxon_big_table.tsv", emit: concatenated_tsv

  script:
  """
  stack_big_tables.py \
      --input-dir sample_big_tables \
      --output taxon_big_table.tsv
  """
}

process TARGET_ENRICHMENT_REPORT {
  /*
   * Build run-level visualizations from per-sample target enrichment
   * summaries emitted by deacon filter on the raw input reads.
   */

  label "low"

  input:
  path "target_enrichment_summaries/*"

  output:
  path "target_enrichment_summary.tsv", emit: summary_tsv
  path "target_enriched_bases_ranked.png", emit: enriched_bases_ranked_png
  path "target_enriched_bases_ranked.html", emit: enriched_bases_ranked_html
  path "target_retained_vs_filtered_stacked.png", emit: retained_vs_filtered_stacked_png
  path "target_retained_vs_filtered_stacked.html", emit: retained_vs_filtered_stacked_html
  path "target_reads_vs_bases_scatter.png", emit: reads_vs_bases_scatter_png
  path "target_reads_vs_bases_scatter.html", emit: reads_vs_bases_scatter_html

  script:
  """
  plot_target_enrichment_summaries.py \\
      --summaries target_enrichment_summaries/*.json \\
      --outdir .
  """

  stub:
  """
  touch target_enrichment_summary.tsv
  touch target_enriched_bases_ranked.png
  touch target_enriched_bases_ranked.html
  touch target_retained_vs_filtered_stacked.png
  touch target_retained_vs_filtered_stacked.html
  touch target_reads_vs_bases_scatter.png
  touch target_reads_vs_bases_scatter.html
  """
}

/*
 * Send a minimal Slack notification for run completion without workflow state.
 */
process NOTIFY_SLACK {

  tag "${workflow.runName}"
  label "low"
  cache false

  errorStrategy 'ignore'

  secret 'SLACK_BOT_TOKEN'

  input:
  val ready
  val sample_set_id
  val labkey_url

  output:
  val true, emit: done

  when:
  params.slack_enabled && params.slack_channel

  script:
  """
  notify_slack.py \
      --run-id '${workflow.runName}' \
      --experiment-id '${params.experiment_id}' \
      --channel '${params.slack_channel}' \
      --sample-set-id '${sample_set_id}' \
      --labkey-url '${labkey_url}' \
      -v || echo "Slack notification failed (non-fatal)"
  """
}
