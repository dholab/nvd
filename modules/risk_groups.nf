process ANNOTATE_SOURMASH_RISK_GROUPS {
  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(tax_summary), path(bioboxes_profile)
  path risk_group_lookup

  output:
  tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.sourmash.tax_metagenome.with_risk_groups.csv")

  script:
  """
  annotate_risk_groups.py sourmash-summary \
      --input ${tax_summary} \
      --lookup ${risk_group_lookup} \
      --bioboxes ${bioboxes_profile} \
      --output ${sample_id}.sourmash.tax_metagenome.with_risk_groups.csv
  """
}


process ANNOTATE_BLAST_RISK_GROUPS {
  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), path(blast_tsv)
  path risk_group_lookup
  val taxonomy_dir

  output:
  tuple val(sample_id), path("${sample_id}.blast.with_risk_groups.tsv")

  script:
  """
  annotate_risk_groups.py blast \
      --input ${blast_tsv} \
      --lookup ${risk_group_lookup} \
      --taxonomy-dir ${taxonomy_dir} \
      --output ${sample_id}.blast.with_risk_groups.tsv
  """
}
