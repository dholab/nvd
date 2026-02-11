/*
 * Deacon: Fast alignment-free decontamination
 * https://github.com/bede/deacon
 *
 * Key features:
 * - Preserves FASTQ headers (critical for read pairing)
 * - Composable indexes via set algebra (union, diff, intersect)
 * - SIMD-accelerated, ~5GB RAM for panhuman index
 */

process DEACON_BUILD_INDEX {
    /*
     * Build a deacon index from FASTA file(s).
     * Use this for custom contaminant sequences.
     */

    tag "${fasta.simpleName}"
    label "medium"

    input:
    path fasta

    output:
    path "*.idx", emit: index

    script:
    def prefix = fasta.simpleName
    """
    deacon index build \\
        --threads ${task.cpus} \\
        -k ${params.deacon_kmer_size} \\
        -w ${params.deacon_window_size} \\
        ${fasta} > ${prefix}.k${params.deacon_kmer_size}w${params.deacon_window_size}.idx
    """
}

process DEACON_FETCH_INDEX {
    /*
     * Download a prebuilt deacon index from URL.
     * Takes the URL as a channel value so the process only runs when
     * the input channel is non-empty (no `when:` guard needed).
     * Caches in work directory; use storeDir for persistent caching.
     */

    label "low"

    input:
    val url

    output:
    path "*.idx", emit: index

    script:
    def filename = url.tokenize('/').last()
    """
    curl -fsSL "${url}" -o ${filename}
    """
}

process DEACON_UNION_INDEXES {
    /*
     * Combine multiple deacon indexes via set union.
     * Only called when both a base index and custom index are present.
     */

    label "low"

    input:
    path indexes  // Collection of .idx files (always 2+)

    output:
    path "combined.idx", emit: index

    script:
    def idx_list = indexes.collect { it.name }.join(' ')
    """
    deacon index union ${idx_list} > combined.idx
    """
}

process DEACON_DEPLETE {
    /*
     * Remove contaminant reads using deacon filter in deplete mode.
     *
     * Critical: This preserves FASTQ headers verbatim, which is required
     * for repair.sh to re-pair reads after filtering. SPAdes paired-end
     * assembly depends on proper read pairing.
     *
     * Deacon natively handles gzipped input/output (since v0.13.0).
     * When writing .gz output via --output, deacon splits --threads 1:1
     * between filtering and compression automatically.
     */

    tag "${sample_id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads), path(index)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.depleted.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.deacon.json"), emit: stats

    script:
    """
    deacon filter \\
        --deplete \\
        --threads ${task.cpus} \\
        --abs-threshold ${params.deacon_abs_threshold} \\
        --rel-threshold ${params.deacon_rel_threshold} \\
        --summary ${sample_id}.deacon.json \\
        --output ${sample_id}.depleted.fastq.gz \\
        ${index} \\
        ${reads}
    """
}
