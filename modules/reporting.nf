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
