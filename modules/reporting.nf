process BUILD_SEQUENCE_FLOW {

  label "low"

  input:
  path input_files, stageAs: "sequence_flow_inputs??????/*"

  output:
  path "sequence_flow.tsv", emit: sequence_flow

  script:
  def files = input_files instanceof List ? input_files : [input_files]
  def input_args = files.collect { input -> "--input '${input}'" }.join(" ")
  """
  build_sequence_flow.py \
      ${input_args} \
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
