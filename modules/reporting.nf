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

process BUILD_SEQUENCE_FLOW {

  label "low"

  input:
  path evidence_files, stageAs: "sequence_flow_evidence??????/*"

  output:
  path "sequence_flow.tsv", emit: sequence_flow

  script:
  def evidence_args = evidence_files.collect { evidence -> "--evidence '${evidence}'" }.join(" ")
  """
  build_sequence_flow.py \
      ${evidence_args} \
      --output sequence_flow.tsv
  """
}

process RENDER_CONTIG_ALIGNMENT_PLOTS {

  tag "${sample_id}"
  label "medium"
  cpus 1
  errorStrategy "ignore"

  input:
  tuple val(sample_id), path(contigs), path(bam), path(bai)

  output:
  tuple val(sample_id), path("${sample_id}.alignoth"), emit: plots

  script:
  """
  mkdir "${sample_id}.alignoth"
  samtools faidx "${contigs}"
  samtools idxstats "${bam}" > contig_alignment_regions.tsv

  while IFS=\$'\\t' read -r contig length mapped _unmapped; do
      if [[ "\${contig}" != "*" && "\${mapped}" -gt 0 ]]; then
          alignoth \
              --bam-path "${bam}" \
              --reference "${contigs}" \
              --region "\${contig}:1-\${length}" \
              --max-read-depth 500 \
              --max-width 1024 \
              --mismatch-display-min-percent 1.0 \
              --html \
              > "${sample_id}.alignoth/\${contig}.alignoth.html"
      fi
  done < contig_alignment_regions.tsv
  """
}
