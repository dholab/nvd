process RENDER_TAXON_ABUNDANCE_SUNBURST {

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(tax_reports)

  output:
  tuple val(sample_id), path("${sample_id}.sourmash.taxburst.html"), path("${sample_id}.sourmash.taxburst.json"), emit: reports

  script:
  def report_files = tax_reports instanceof List ? tax_reports : [tax_reports]
  def summarized_csv = report_files.find { report -> report.name.endsWith(".summarized.csv") }
  assert summarized_csv : "Missing sourmash csv_summary report for ${sample_id}."
  """
  taxburst \
      -F csv_summary \
      "${summarized_csv}" \
      -o "${sample_id}.sourmash.taxburst.html" \
      --save-json "${sample_id}.sourmash.taxburst.json"
  """
}

process RENDER_SOURMASH_SANKEY {

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(tax_reports)

  output:
  tuple val(sample_id), path("${sample_id}.sourmash.sankey.html"), emit: report

  script:
  def report_files = tax_reports instanceof List ? tax_reports : [tax_reports]
  def summarized_csv = report_files.find { report -> report.name.endsWith(".summarized.csv") }
  assert summarized_csv : "Missing sourmash csv_summary report for ${sample_id}."
  """
  sourmash scripts sankey \
      --summary-csv "${summarized_csv}" \
      -o "${sample_id}.sourmash.sankey.html" \
      --title "${sample_id} sourmash taxonomic profile"
  """
}
