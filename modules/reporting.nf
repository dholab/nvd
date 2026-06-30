process BUILD_TAXONOMIC_PROFILE_NORMALIZATION_MAP {

  label "low"

  input:
  path(tax_reports)

  output:
  path("taxonomic_profile.normalization_map.json"), emit: map

  script:
  def report_files = tax_reports instanceof List ? tax_reports.flatten() : [tax_reports]
  def summarized_csvs = report_files.findAll { report -> report.name.endsWith(".summarized.csv") }
  assert summarized_csvs : "Missing taxonomic profile summary CSVs."
  def input_args = summarized_csvs.collect { report -> "--input \"${report}\"" }.join(" ")
  """
  normalize_taxonomic_profile_summary.py build-map \
      ${input_args} \
      --output taxonomic_profile.normalization_map.json
  """
}

process NORMALIZE_TAXONOMIC_PROFILE_SUMMARY {

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(tax_reports), path(normalization_map)

  output:
  tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.taxonomic_profile.summary.normalized.csv"), emit: summary

  script:
  def report_files = tax_reports instanceof List ? tax_reports : [tax_reports]
  def summarized_csv = report_files.find { report -> report.name.endsWith(".summarized.csv") }
  assert summarized_csv : "Missing taxonomic profile summary CSV for ${sample_id}."
  """
  normalize_taxonomic_profile_summary.py apply-map \
      --map "${normalization_map}" \
      --input "${summarized_csv}" \
      --output "${sample_id}.taxonomic_profile.summary.normalized.csv"
  """
}

process RENDER_MERGED_TAXON_ABUNDANCE_SUNBURST {

  label "low"

  input:
  tuple val(sample_ids), path(profile_summaries)

  output:
  path("sourmash.taxburst.html"), emit: report

  script:
  def pairs = [sample_ids, profile_summaries].transpose().sort { left, right -> left[0] <=> right[0] }
  def summary_args = pairs.collect { sample_id, summary -> "--summary \"${sample_id}=${summary}\"" }.join(" ")
  """
  render_multisample_taxburst.py \
      ${summary_args} \
      --output sourmash.taxburst.html
  """
}

process RENDER_TAXON_ABUNDANCE_SUNBURST {

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(profile_summary)

  output:
  tuple val(sample_id), path("${sample_id}.sourmash.taxburst.html"), path("${sample_id}.sourmash.taxburst.json"), emit: reports

  script:
  """
  ln -s "${profile_summary}" "${sample_id}.csv"

  taxburst \
      -F csv_summary \
      "${sample_id}.csv" \
      -o "${sample_id}.sourmash.taxburst.html" \
      --save-json "${sample_id}.sourmash.taxburst.json"
  """
}

process RENDER_SOURMASH_SANKEY {

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(profile_summary)

  output:
  tuple val(sample_id), path("${sample_id}.sourmash.sankey.html"), emit: report

  script:
  """
  ln -s "${profile_summary}" "${sample_id}.csv"

  sourmash scripts sankey \
      --summary-csv "${sample_id}.csv" \
      -o "${sample_id}.sourmash.sankey.html" \
      --title "${sample_id} sourmash taxonomic profile"
  """
}
