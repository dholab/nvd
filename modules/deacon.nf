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

process DEACON_BUILD_INDEX_FROM_STAT_K_MERS {
    /*
     * Build a deacon index from k-mers extracted from a STAT .dbss database.
     *
     * Converts STAT's MSB-first 2-bit k-mers to deacon's LSB-first format,
     * filtered by a list of target tax IDs (e.g., human-infecting viruses).
     * Uses k=31, w=1 with dual truncation for maximum sensitivity — captures
     * 100% of STAT's virus signal while running 60-80x faster per sample.
     *
     * This runs ONCE per pipeline invocation (not per-sample) since the
     * reference files are the same for all samples. The resulting .idx is
     * then .combine()'d with per-sample reads for filtering.
     */

    label "medium"

    input:
    path stat_dbss
    path stat_annotation
    path human_virus_taxlist

    output:
    path "human_viruses.k31w1.idx", emit: index

    script:
    """
    rust-script ${projectDir}/bin/stat_to_deacon.rs \\
        --target-k 31 \\
        --window-size 1 \\
        --dbss ${stat_dbss} \\
        --annotation ${stat_annotation} \\
        --taxids ${human_virus_taxlist} \\
        --output human_viruses.k31w1.idx
    """
}

process DEACON_FILTER_HUMAN_VIRUS_READS {
    /*
     * Extract human virus reads using deacon filter in non-deplete mode.
     *
     * Replaces STAT aligns_to + seqkit grep for human virus read extraction.
     * Uses a deacon index built from STAT k-mers (via DEACON_BUILD_INDEX_FROM_STAT_K_MERS).
     * Threshold settings (abs=1, rel=0.0) match STAT's single-kmer-hit behavior.
     *
     * Accepts pre-interleave input directly from GATHER_READS:
     * - Paired (Illumina): takes R1/R2 as separate files. Deacon does pair-aware
     *   k-mer counting and outputs interleaved FASTQ in a single pass, replacing
     *   both INTERLEAVE_PAIRS and the old STAT extraction step.
     * - Single-end (ONT): takes a single file, filtered as individual reads.
     *
     * Output is always tuple(sample_id, platform, read_structure, fastq) matching
     * the shape expected by PREPROCESS_READS and downstream SPAdes.
     */

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), path(reads), path(reads2), path(deacon_idx)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.human_virus.fastq.gz")

    script:
    read_structure = reads2.name != "NO_R2" ? "interleaved" : "single"
    if (read_structure == "interleaved")
        """
        deacon filter \\
            --threads ${task.cpus} \\
            --abs-threshold 1 \\
            --rel-threshold 0.0 \\
            --output ${sample_id}.human_virus.fastq.gz \\
            ${deacon_idx} \\
            ${reads} ${reads2}
        """
    else
        """
        deacon filter \\
            --threads ${task.cpus} \\
            --abs-threshold 1 \\
            --rel-threshold 0.0 \\
            --output ${sample_id}.human_virus.fastq.gz \\
            ${deacon_idx} \\
            ${reads}
        """
}