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
