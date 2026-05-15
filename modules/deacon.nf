/*
 * Deacon: Fast alignment-free decontamination
 * https://github.com/bede/deacon
 *
 * Key features:
 * - Preserves FASTQ headers (critical for read pairing)
 * - Composable indexes via set algebra (union, diff, intersect)
 * - SIMD-accelerated, ~5GB RAM for panhuman index
 */

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
    def union_or_link = indexes.size() == 1
        ? "ln -s ${idx_list} combined.idx"
        : "deacon index union ${idx_list} > combined.idx"
    """
    ${union_or_link}
    """
}

process DEACON_BUILD_INDEX_FROM_FASTA {
    /*
     * Build a host/contaminant deacon index from a custom FASTA file.
     */

    label "medium"

    input:
    path contaminants_fasta

    output:
    path "custom_contaminants.k${params.host_kmer_size}w${params.host_window_size}.idx", emit: index

    script:
    """
    deacon index build \
        -k ${params.host_kmer_size} \
        -w ${params.host_window_size} \
        -o custom_contaminants.k${params.host_kmer_size}w${params.host_window_size}.idx \
        ${contaminants_fasta}
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
        --abs-threshold ${params.host_abs_threshold} \\
        --rel-threshold ${params.host_rel_threshold} \\
        --summary ${sample_id}.deacon.json \\
        --output ${sample_id}.depleted.fastq.gz \\
        ${index} \\
        ${reads}
    """
}

process DEACON_BUILD_VIRUS_INDEX_FROM_FASTA {
    /*
     * Build a deacon virus index from a reference FASTA file.
     *
     * Use this when providing virus_reference_fasta instead of a
     * pre-built index (virus_index) or URL (virus_index_url).
     * Uses virus_kmer_size/virus_window_size for maximum sensitivity
     * (k=31, w=1 by default) to match the behavior of STAT-derived indexes.
     *
     * This runs ONCE per pipeline invocation (not per-sample).
     */

    label "medium"

    input:
    path virus_fasta

    output:
    path "virus_reference.k${params.virus_kmer_size}w${params.virus_window_size}.idx", emit: index

    script:
    """
    deacon index build \\
        -k ${params.virus_kmer_size} \\
        -w ${params.virus_window_size} \\
        -o virus_reference.k${params.virus_kmer_size}w${params.virus_window_size}.idx \\
        ${virus_fasta}
    """
}

process DEACON_FILTER_HUMAN_VIRUS_READS {
    /*
     * Extract human virus reads using deacon filter in non-deplete mode.
     *
     * Accepts pre-interleave input directly from GATHER_READS:
     * - Paired (Illumina): takes R1/R2 as separate files. Deacon does pair-aware
     *   k-mer counting and outputs interleaved FASTQ in a single pass.
     * - Single-end (ONT): takes a single file, filtered as individual reads.
     *
     * Output is always tuple(sample_id, platform, read_structure, fastq) matching
     * the shape expected by downstream preprocessing and SPAdes.
     */

    tag "${sample_id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), path(reads), path(reads2), path(deacon_idx)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.human_virus.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.deacon_filter.json"), emit: stats

    script:
    read_structure = reads2.name != "NO_R2" ? "interleaved" : "single"
    if (read_structure == "interleaved")
        """
        deacon filter \\
            --threads ${task.cpus} \\
            --abs-threshold ${params.virus_abs_threshold} \\
            --rel-threshold ${params.virus_rel_threshold} \\
            --summary ${sample_id}.deacon_filter.json \\
            --output ${sample_id}.human_virus.fastq.gz \\
            ${deacon_idx} \\
            ${reads} ${reads2}
        """
    else
        """
        deacon filter \\
            --threads ${task.cpus} \\
            --abs-threshold ${params.virus_abs_threshold} \\
            --rel-threshold ${params.virus_rel_threshold} \\
            --summary ${sample_id}.deacon_filter.json \\
            --output ${sample_id}.human_virus.fastq.gz \\
            ${deacon_idx} \\
            ${reads}
        """
}

process DEACON_FILTER_CONTIGS {
    /*
     * Extract human virus contigs using deacon filter on assembled FASTA.
     *
     * Replaces the 5-process STAT contig classification chain:
     * CLASSIFY_CONTIGS_FIRST_PASS → GENERATE_CONTIGS_TAXA_LIST →
     * CLASSIFY_CONTIGS_SECOND_PASS → IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS →
     * EXTRACT_HUMAN_VIRUS_CONTIGS
     *
     * Uses the same virus index as read extraction so no additional reference
     * is needed. Contigs are already ≥200 bp (FILTER_SHORT_CONTIGS), so
     * incidental single k-mer matches are rare — the signal-to-noise ratio
     * improves vs reads because contigs have far more k-mers per sequence.
     */

    tag "${sample_id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(fasta), path(deacon_idx)

    output:
    tuple val(sample_id), path("${sample_id}.human_virus.fasta")

    script:
    """
    deacon filter \\
        --threads ${task.cpus} \\
        --abs-threshold ${params.virus_abs_threshold} \\
        --rel-threshold ${params.virus_rel_threshold} \\
        --output ${sample_id}.human_virus.fasta \\
        ${deacon_idx} \\
        ${fasta}
    """
}
