/*
 * sourmash: sketch-based metagenomic analysis
 *
 * Experimental branchwater-backed commands live here while the sourmash branch
 * is still being shaped. Keep the interface aligned with post-preprocessing NVD
 * read tuples so later gather/tax steps can extend the same subworkflow.
 */

process SOURMASH_SKETCH_QUERY_METAGENOME {

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(reads)

  output:
  tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.sourmash.query_metagenome.k${params.sourmash_ksize}.scaled${params.sourmash_scaled}.abund.sig.gz"), emit: query_sketches

  script:
  """
  sourmash scripts singlesketch \
      ${reads} \
      --name "${sample_id}" \
      -p dna,k=${params.sourmash_ksize},scaled=${params.sourmash_scaled},abund \
      -o "${sample_id}.sourmash.query_metagenome.k${params.sourmash_ksize}.scaled${params.sourmash_scaled}.abund.sig.gz"
  """
}

process SOURMASH_FETCH_REF_SKETCH {

  tag { url.tokenize('/').last() }
  label "low"

  input:
  val url

  output:
  path "*.sig.zip", arity: '1', emit: ref_sketch

  script:
  def filename = url.tokenize('/').last()
  """
  curl -fsSL "${url}" -o "${filename}"
  """
}

process SOURMASH_SKETCH_REF_FASTA {

  tag "${ref_fasta.simpleName}"
  label "medium"

  input:
  path ref_fasta

  output:
  path "sourmash_reference.k${params.sourmash_ksize}.scaled${params.sourmash_scaled}.sig.zip", emit: ref_sketch

  script:
  """
  sourmash sketch dna \
      ${ref_fasta} \
      --singleton \
      -p dna,k=${params.sourmash_ksize},scaled=${params.sourmash_scaled} \
      -o "sourmash_reference.k${params.sourmash_ksize}.scaled${params.sourmash_scaled}.sig.zip"
  """
}

process SOURMASH_GATHER_QUERY_METAGENOME {

  tag "${sample_id}"
  label "medium"
  cpus 2

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(query_sketch), path(ref_sketch)

  output:
  tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.sourmash.gather.csv"), emit: gather_csv

  script:
  """
  sourmash scripts fastgather \
      ${query_sketch} \
      ${ref_sketch} \
      -k ${params.sourmash_ksize} \
      --scaled ${params.sourmash_scaled} \
      --moltype DNA \
      --cores ${task.cpus} \
      --threshold-bp ${params.sourmash_threshold_bp} \
      -o "${sample_id}.sourmash.gather.csv"
  touch "${sample_id}.sourmash.gather.csv"
  """
}

process SOURMASH_FETCH_LINEAGES {

  tag { url.tokenize('/').last() }
  label "low"

  input:
  val url

  output:
  path "*.csv", arity: '1', emit: lineages

  script:
  def filename = url.tokenize('/').last()
  """
  curl -fsSL "${url}" -o "${filename}"
  """
}

process SOURMASH_STAGE_REFERENCE {

  label "low"

  input:
  tuple path(ref_sketch), path(lineages)

  output:
  path(ref_sketch), emit: ref_sketch
  path(lineages), emit: lineages

  script:
  """
  true
  """
}

process SOURMASH_TAX_METAGENOME {

  tag "${sample_id}"
  label "low"

  input:
  tuple val(sample_id), val(platform), val(read_structure), path(gather_csv), path(lineages_csv)

  output:
  tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.sourmash.tax_metagenome*"), emit: tax_reports

  script:
  def output_formats = params.sourmash_output_bioboxes
      ? "csv_summary lineage_summary krona kreport bioboxes"
      : "csv_summary lineage_summary krona kreport"
  """
  sourmash tax metagenome \
      --gather-csv ${gather_csv} \
      --taxonomy-csv ${lineages_csv} \
      --keep-identifier-versions \
      --use-abundances \
      --output-format ${output_formats} \
      --rank species \
      --output-base "${sample_id}.sourmash.tax_metagenome"
  """
}

process SOURMASH_COLLECT_QUERY_SKETCHES {

  label "low"

  input:
  path query_sketches

  output:
  path ("query_metagenomes.k${params.sourmash_ksize}.scaled${params.sourmash_scaled}.zip"), emit: collection

  script:
  def sketch_list = query_sketches.collect { sketch -> "\"${sketch}\"" }.join(" ")
  """
  sourmash sig cat \
      ${sketch_list} \
      -o query_metagenomes.k${params.sourmash_ksize}.scaled${params.sourmash_scaled}.zip
  """
}

process SOURMASH_COMPARE_QUERY_SKETCHES {

  tag "${metric}"
  label "medium"

  input:
  tuple val(metric), path(query_sketch_collection)

  output:
  tuple val(metric), path("query_sketch_similarity.${metric}.npy"), path("query_sketch_similarity.${metric}.csv"), emit: matrices

  script:
  def ignore_abundance = metric == "noabund" ? "--ignore-abundance" : ""
  """
  sourmash compare \
      ${query_sketch_collection} \
      ${ignore_abundance} \
      --output query_sketch_similarity.${metric}.npy \
      --csv query_sketch_similarity.${metric}.csv
  """
}
