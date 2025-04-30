process FETCH_FASTQ {

    tag "${run_accession}, ${platform}"

    maxForks params.max_concurrent_downloads
    cpus 3

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    input:
    tuple val(platform), val(run_accession)

    output:
    tuple val(run_accession), path("${run_accession}*.fastq")

    script:
    """
    prefetch ${run_accession}
    fasterq-dump --split-3 ${run_accession}
    """
}
