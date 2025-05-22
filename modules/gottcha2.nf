process READ_COMPRESSION_PASSTHROUGH {

    input:
    tuple val(sample_id), path(fastq), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    tuple val(sample_id), path("*.fa*"), path(ref_mmi), path(stats), path(tax_tsv)

    script:
    def extension = file(fastq).getBaseName().replace(".gz", "")
    if (file(fastq).getExtension().contains("gz"))
        """
        cat ${fastq} | gzip -c -d > ${sample_id}.${extension}
        """
    else
        """
        cp ${fastq} ${sample_id}.uncompressed.${extension}
        """
}

process GOTTCHA2_PROFILE_NANOPORE {

    tag "${sample_id}"
    // publishDir params.gottcha_sam, mode: 'copy', overwrite: false, pattern: "*.sam"
    // publishDir params.gottcha_stats, mode: 'copy', overwrite: false, pattern: "*.tsv"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 12

    input:
    tuple val(sample_id), path(fastq), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    tuple val(sample_id), path("${sample_id}*.sam"), path(ref_mmi), path(stats), path(tax_tsv), emit: aligned
    tuple val(sample_id), path("${sample_id}*.full.tsv"), path(ref_mmi), path(stats), path(tax_tsv), emit: full_tsv
    path "*.tsv", emit: all_stats

    when:
    (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha")) || params.all || params.gottcha2

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py  \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --noCutoff --dbLevel strain --threads ${task.cpus} \
    --nanopore \
    --input ${fastq}
    """
}

process GOTTCHA2_PROFILE_ILLUMINA {

    tag "${sample_id}"
    // publishDir params.gottcha_sam, mode: 'copy', overwrite: false, pattern: "*.sam"
    // publishDir params.gottcha_stats, mode: 'copy', overwrite: false, pattern: "*.tsv"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 12

    input:
    tuple val(sample_id), path(fastq), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    tuple val(sample_id), path("${sample_id}*.sam"), path(ref_mmi), path(stats), path(tax_tsv), emit: aligned
    tuple val(sample_id), path("${sample_id}*.full.tsv"), path(ref_mmi), path(stats), path(tax_tsv), emit: full_tsv
    path "*.tsv", emit: all_stats

    when:
    (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha")) || params.all || params.gottcha2

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --noCutoff --dbLevel strain --threads ${task.cpus} \
    --input ${fastq}
    """
}

process GENERATE_FASTA {

    tag "${sample_id}"
    // publishDir params.gottcha_fasta, mode: 'copy', overwrite: false, pattern: "*.fasta"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    input:
    tuple val(sample_id), path(sam), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    path "*"

    when:
    (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha")) || params.all || params.gottcha

    script:
    def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
    """
    gottcha2.py \
    --noCutoff --dbLevel strain --threads ${task.cpus} -ef \
    --database ${ref_prefix} \
    --prefix ${sample_id} \
    --sam ${sam}
    """
}
