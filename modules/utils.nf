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
  tuple val(sample_id), path(blast_tsv), val(total_reads), path(contig_counts), path(query_lookups)
  val run_id

  output:
  tuple val(sample_id), path("${sample_id}_blast.final.tsv")

  script:
  def virus_index_version = NvdUtils.targetEnrichmentEnabled(params) ? params.virus_index_version : "not_used"
  def lookup_files = query_lookups instanceof List ? query_lookups : [query_lookups]
  def lookup_args = lookup_files.collect { lookup -> "--query-lookup ${lookup}" }.join(" ")
  """
  finalize_blast_results.py \\
      --blast-tsv ${blast_tsv} \\
      --contig-counts ${contig_counts} \\
      ${lookup_args} \\
      --output ${sample_id}_blast.final.tsv \\
      --total-reads ${total_reads} \\
      --blast-db-version '${params.blast_db_version}' \\
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
  tuple val(sample_id), path(blast_tsv)

  output:
  tuple val(sample_id), path("${sample_id}.query_big_table.tsv")

  script:
  """
  build_query_big_table.py \
      --blast-tsv ${blast_tsv} \
      --output ${sample_id}.query_big_table.tsv
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
  #!/usr/bin/env python3
  from pathlib import Path

  files = sorted(Path("sample_big_tables").glob("*.tsv"))

  if not files:
      Path("query_big_table.tsv").write_text("")
  else:
      with open("query_big_table.tsv", "w") as out:
          with open(files[0]) as f:
              header = f.readline()
              out.write(header)
              for line in f:
                  out.write(line)
          for tsv in files[1:]:
              with open(tsv) as f:
                  f.readline()
                  for line in f:
                      out.write(line)
  """
}

process VIRUS_ENRICHMENT_REPORT {
  /*
   * Build run-level visualizations from per-sample human-virus enrichment
   * summaries emitted by deacon filter on the raw input reads.
   */

  label "low"

  input:
  path "virus_enrichment_summaries/*"

  output:
  path "virus_enrichment_summary.tsv", emit: summary_tsv
  path "virus_enriched_bases_ranked.png", emit: enriched_bases_ranked_png
  path "virus_enriched_bases_ranked.html", emit: enriched_bases_ranked_html
  path "virus_retained_vs_filtered_stacked.png", emit: retained_vs_filtered_stacked_png
  path "virus_retained_vs_filtered_stacked.html", emit: retained_vs_filtered_stacked_html
  path "virus_reads_vs_bases_scatter.png", emit: reads_vs_bases_scatter_png
  path "virus_reads_vs_bases_scatter.html", emit: reads_vs_bases_scatter_html

  script:
  """
  plot_virus_enrichment_summaries.py \\
      --summaries virus_enrichment_summaries/*.json \\
      --outdir .
  """

  stub:
  """
  touch virus_enrichment_summary.tsv
  touch virus_enriched_bases_ranked.png
  touch virus_enriched_bases_ranked.html
  touch virus_retained_vs_filtered_stacked.png
  touch virus_retained_vs_filtered_stacked.html
  touch virus_reads_vs_bases_scatter.png
  touch virus_reads_vs_bases_scatter.html
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
