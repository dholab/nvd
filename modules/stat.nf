process EXTRACT_HUMAN_VIRUS_READS {

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(sample_type), path(fastq), path(stat_dbss), path(stat_annotation), path(human_virus_taxlist)

    output:
    tuple val(sample_id), val(sample_type), path("*.f*q.gz")

    when:
    params.tools && (params.tools.contains("nvd") || params.tools.contains("stat") || params.tools.contains("blast") || params.tools.contains("stat_blast") || params.tools.contains("stast") || params.tools.contains("all"))

    script:
    """
    seqkit fq2fa --threads 1 ${fastq} -o ${sample_id}.all_reads.fasta
    aligns_to \
    -dbss ${stat_dbss} \
    -num_threads ${task.cpus} \
    -tax_list ${human_virus_taxlist} \
    ${sample_id}.all_reads.fasta \
    | cut -f1 \
    | seqkit grep -f - \
    ${fastq} \
    -o ${sample_id}.human_virus.fastq.gz

    # Cleanup intermediate large file
    rm ${sample_id}.all_reads.fasta
    """
}

process CLASSIFY_CONTIGS_FIRST_PASS {

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(platform), path(filtered_contigs), path(stat_index)

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

    output:
    tuple val(sample_id), path("${sample_id}.report"), emit: report

    script:
    """
    hits_to_report.py \
    --cutoff-percent ${params.cutoff_percent} \
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

    output:
    tuple val(sample_id), path("${sample_id}.report")

    script:
    // Handle both list (from config) and comma-separated string (from CLI)
    def virus_families = params.human_virus_families instanceof List
        ? params.human_virus_families.join(" ")
        : params.human_virus_families.toString().split(",").collect { it.trim() }.join(" ")
    """
    extract_taxa_spots.py \
    --hits_file ${secondpass_file} \
    --output_file ${sample_id}.report \
    --taxa ${virus_families} \
    --stringency ${params.tax_stringency} \
    --include_children
    """
}

process SCRUB_HUMAN_READS {
    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(fastq), path(human_reads)

    output:
    tuple val(sample_id), path("${sample_id}.scrubbed.fastq.gz")

    when:
	params.tools && (params.tools.contains("clump") || params.tools.contains("all"))

    script:
    """
    # Convert to FASTA
    seqkit fq2fa --threads 1 ${fastq} -o ${sample_id}.all_reads.fasta
    
    # Get human-aligned reads (tabular output like the virus example)
    aligns_to \
        -db ${human_reads} \
        -num_threads ${task.cpus} \
        ${sample_id}.all_reads.fasta \
        | cut -f1 > human_read_ids.txt
    
    # Remove human reads (inverse grep with -v)
    seqkit grep -v -f human_read_ids.txt ${fastq} \
        -o ${sample_id}.scrubbed.fastq.gz
    """
}
