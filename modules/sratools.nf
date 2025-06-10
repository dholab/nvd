process FETCH_FASTQ {

    tag "${run_accession}, ${platform}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    maxForks params.max_concurrent_downloads
    cpus 3

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    input:
    tuple val(id), val(platform), val(run_accession)

    output:
    tuple val(id), val(platform), path("${run_accession}*.fastq")

    script:
    id = id ? id : run_accession
    """
    prefetch ${run_accession}
    fasterq-dump --split-3 ${run_accession}
    """
}
