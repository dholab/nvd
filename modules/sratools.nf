process FETCH_FASTQ {

    tag "${run_accession}, ${platform}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    maxForks params.max_concurrent_downloads

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    input:
    tuple val(id), val(platform), val(run_accession)

    output:
    tuple val(id), val(platform), path("${run_accession}*.fastq.gz")

    script:
    id = id ? id : run_accession
    """
    sracha get \
        --split split-3 \
        --gzip-level 1 \
        --threads ${task.cpus} \
        --prefer-ena \
        --no-progress \
        --yes \
        --output-dir . \
        ${run_accession}
    """
}
