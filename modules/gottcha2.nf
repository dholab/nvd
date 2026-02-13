process GOTTCHA2_PROFILE_NANOPORE {

    tag "${sample_id}"
    label "ludicrous"

    publishDir params.gottcha2_profiles, mode: 'copy', overwrite: false, pattern: "*.sam"
    publishDir params.gottcha2_profiles, mode: 'copy', overwrite: false, pattern: "*.tsv"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(fastq), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    tuple val(sample_id), path("${sample_id}*.sam"), path(ref_mmi), path(stats), path(tax_tsv), emit: aligned
    tuple val(sample_id), path("${sample_id}*.full.tsv"), path(ref_mmi), path(stats), path(tax_tsv), emit: full_tsv
    path "*.tsv", emit: all_stats

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py  \
    --verbose --debug \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --noCutoff --dbLevel strain --threads ${task.cpus} \
    --nanopore \
    --input ${fastq}
    """
}

process GOTTCHA2_PROFILE_ILLUMINA {

    tag "${sample_id}"
    label "ludicrous"

    publishDir params.gottcha2_profiles, mode: 'copy', overwrite: false, pattern: "*.sam"
    publishDir params.gottcha2_profiles, mode: 'copy', overwrite: false, pattern: "*.tsv"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(fastq), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    tuple val(sample_id), path("${sample_id}*.sam"), path(ref_mmi), path(stats), path(tax_tsv), emit: aligned
    tuple val(sample_id), path("${sample_id}*.full.tsv"), path(ref_mmi), path(stats), path(tax_tsv), emit: full_tsv
    path "*.tsv", emit: all_stats

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    --verbose --debug \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --noCutoff --dbLevel strain --threads ${task.cpus} \
    --input ${fastq}
    """
}

process GENERATE_FASTA {

    tag "${sample_id}"
    label "medium"

    publishDir params.extracted_reads, mode: 'copy', overwrite: false, pattern: "*.fasta"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(sam), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    tuple val(sample_id), path("*.extract.fasta"), path("*.full.tsv"), path("*.gottcha_strain.log"), path("*.lineage.tsv"), optional: true

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    --verbose --debug \
    --noCutoff --dbLevel strain --threads ${task.cpus} -ef \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --sam ${sam}

    if [[ -f ${sample_id}.tsv ]]; then
        echo "Renaming the lineage tsv '${sample_id}.tsv' to include the lineage suffix..."
        mv ${sample_id}.tsv ${sample_id}.lineage.tsv
    elif [[ -f ${sample_id}.lineage.tsv ]]; then
        echo "A valid lineage TSV file was successfully produced"
    else
        echo "[WARN] No valid lineage TSV was produced. Downstream processes that depend on it will not run."
    fi
    """
}
// Gottcha2 -ef meaning
// Extract up to 20 sequences per reference from the SAM file and save them to a FASTA file. Equivalent to using: -e 'all:20:fasta'.

/*
 * Register GOTTCHA2 hits with idempotent keys in the state database.
 *
 * This process computes deterministic hit keys from extracted sequences,
 * enabling cross-run and cross-sample deduplication. Taxonomy is parsed
 * directly from the GOTTCHA2 FASTA headers (LEVEL, NAME, TAXID).
 *
 * Input:
 *   - tuple of (sample_id, fasta, full_tsv, sample_set_id, state_dir, hits_dir)
 *   - state_dir may be empty string in stateless mode
 *   - hits_dir is the base directory for parquet output (state_dir/hits or results/hits)
 * Output:
 *   - tuple of (sample_id, log_file, parquet_file)
 */
process REGISTER_GOTTCHA2_HITS {

    tag "${sample_id}"
    label "low"
    maxForks 1

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    // Publish parquet file - runs on every resume even if task is cached
    // Uses Hive-partitioned structure: {hits_dir}/month=NULL/{sample_set_id}/{sample_id}/data.parquet
    publishDir "${hits_dir}/month=NULL/${sample_set_id}/${sample_id}", mode: 'copy', pattern: "data.parquet"

    input:
    tuple val(sample_id), path(fasta), path(full_tsv), val(sample_set_id), val(state_dir), val(hits_dir)

    output:
    tuple val(sample_id), path("${sample_id}_gottcha2_hits_registered.log"), path("data.parquet")

    script:
    assert hits_dir : "hits_dir cannot be null - this indicates a workflow configuration error"
    def labkey_arg = params.labkey ? "--labkey" : ""
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    """
    set -o pipefail
    register_gottcha2_hits.py \\
        --fasta ${fasta} \\
        --full-tsv ${full_tsv} \\
        ${state_dir_arg} \\
        --sample-set-id '${sample_set_id}' \\
        --sample-id '${sample_id}' \\
        --run-id '${workflow.runName}' \\
        --output data.parquet \\
        ${labkey_arg} \\
        -v 2>&1 | tee ${sample_id}_gottcha2_hits_registered.log
    """
}
