process EXTRACT_HUMAN_VIRUS_READS {

    tag "${sample_id}"
    label "high_disk"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(fastq), path(stat_dbss), path(stat_annotation), path(human_virus_taxlist)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("*.f*q.gz")

    when:
    params.tools && (params.tools.contains("nvd") || params.tools.contains("stat") || params.tools.contains("blast") || params.tools.contains("stat_blast") || params.tools.contains("stast") || params.tools.contains("all"))

    script:
    """
    seqkit fq2fa --threads 1 ${fastq} -o - \\
    | aligns_to \\
    -dbss ${stat_dbss} \\
    -num_threads ${task.cpus} \\
    -tax_list ${human_virus_taxlist} \\
    stdin \\
    | cut -f1 \\
    | seqkit grep -f - \\
    ${fastq} \\
    -o ${sample_id}.human_virus.fastq.gz
    """
}

process CLASSIFY_CONTIGS_FIRST_PASS {

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(filtered_contigs), path(stat_index)

    output:
    tuple val(sample_id), path(filtered_contigs), path("${sample_id}.firstpass.txt")

    script:
    """
    aligns_to -dbs ${stat_index} \
    -num_threads ${task.cpus} \
    ${filtered_contigs} \
    > "${sample_id}.firstpass.txt"
    """
}

process GENERATE_CONTIGS_TAXA_LIST {

    tag "${sample_id}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(filtered_contigs), path(first_pass_file)

    output:
    tuple val(sample_id), path(filtered_contigs), path("${sample_id}.report")

    script:
    """ 
    generate_tax_list.py ${first_pass_file} ${sample_id}.report
    """
}

process CLASSIFY_CONTIGS_SECOND_PASS {

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 8

    input:
    tuple val(sample_id), path(filtered_contigs), path(tax_list), path(stat_dbss), path(stat_annotation)

    output:
    tuple val(sample_id), path("${sample_id}.secondpass.txt")

    script:
    """
    aligns_to \
    -tax_list ${tax_list} \
    -dbss ${stat_dbss} \
    -num_threads ${task.cpus} \
    ${filtered_contigs} \
    > "${sample_id}.secondpass.txt"
    """
}


process GENERATE_STAT_CONTIG_REPORT {

    tag "${sample_id}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(secondpass_file)
    val state_dir
    val taxonomy_dir

    output:
    tuple val(sample_id), path("${sample_id}.report"), emit: report

    script:
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
    """
    hits_to_report.py \
    --cutoff-percent ${params.cutoff_percent} \
    ${state_dir_arg} \
    ${taxonomy_dir_arg} \
    ${secondpass_file} \
    ${sample_id}.report
    """
}

process IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS {

    tag "${sample_id}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(secondpass_file)
    val state_dir
    val taxonomy_dir

    output:
    tuple val(sample_id), path("${sample_id}.report")

    script:
    // Handle both list (from config) and comma-separated string (from CLI)
    def virus_families = params.human_virus_families instanceof List
        ? params.human_virus_families.join(" ")
        : params.human_virus_families.toString().split(",").collect { it.trim() }.join(" ")
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    def taxonomy_dir_arg = taxonomy_dir ? "--taxonomy-dir '${taxonomy_dir}'" : ""
    """
    extract_taxa_spots.py \
    --hits_file ${secondpass_file} \
    --output_file ${sample_id}.report \
    --taxa ${virus_families} \
    --stringency ${params.tax_stringency} \
    --include_children \
    ${state_dir_arg} \
    ${taxonomy_dir_arg}
    """
}

process SCRUB_HOST_READS {
    /*
     * Remove host (human) reads using STAT's aligns_to + seqkit grep.
     * This approach preserves original FASTQ headers (including Casava 1.8+ metadata),
     * which is required for downstream pair repair.
     */

    tag "${sample_id}"
    label "high_disk"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(fastq), path(human_db)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.scrubbed.fastq.gz")

    script:
    """
    seqkit fq2fa --threads 1 ${fastq} -o - \\
    | aligns_to \\
        -db ${human_db} \\
        -num_threads ${task.cpus} \\
        stdin \\
    | cut -f1 \\
    | seqkit grep -v -f - \\
        ${fastq} \\
        -o ${sample_id}.scrubbed.fastq.gz
    """
}

process SCRUB_HUMAN_READS {
    /*
     * Simplified host scrubbing for clumpify/SRA submission workflow.
     * Uses same STAT approach but with simpler tuple structure.
     */

    tag "${sample_id}"
    label "high_disk"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(fastq), path(human_db)

    output:
    tuple val(sample_id), path("${sample_id}.scrubbed.fastq.gz")

    when:
    params.tools && (params.tools.contains("clump") || params.tools.contains("all"))

    script:
    """
    seqkit fq2fa --threads 1 ${fastq} -o - \\
    | aligns_to \\
        -db ${human_db} \\
        -num_threads ${task.cpus} \\
        stdin \\
    | cut -f1 \\
    | seqkit grep -v -f - \\
        ${fastq} \\
        -o ${sample_id}.scrubbed.fastq.gz
    """
}
