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

    tag "${meta.id}, ${meta.query_class}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("${meta.id}.${meta.query_class}.depleted.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.${meta.query_class}.deacon.json"), emit: stats

    script:
    """
    deacon filter \
        --deplete \
        --threads ${task.cpus} \
        --abs-threshold ${params.host_abs_threshold} \
        --rel-threshold ${params.host_rel_threshold} \
        --summary ${meta.id}.${meta.query_class}.deacon.json \
        --output ${meta.id}.${meta.query_class}.depleted.fastq.gz \
        ${index} \
        ${reads}
    """
}

process DEACON_BUILD_TARGET_INDEX_FROM_FASTA {
    /*
     * Build a deacon target index from a reference FASTA file.
     *
     * Use this when providing virus_reference_fasta instead of a
     * pre-built index (virus_index) or URL (virus_index_url).
     * Uses virus_kmer_size/virus_window_size for maximum sensitivity
     * (k=31, w=1 by default).
     *
     * This runs ONCE per pipeline invocation (not per-sample).
     */

    label "medium"

    input:
    path target_fasta

    output:
    path "target_reference.k${params.virus_kmer_size}w${params.virus_window_size}.idx", emit: index

    script:
    """
    deacon index build \
        -k ${params.virus_kmer_size} \
        -w ${params.virus_window_size} \
        -o target_reference.k${params.virus_kmer_size}w${params.virus_window_size}.idx \
        ${target_fasta}
    """
}

process DEACON_ENRICH_TARGET_READS {
    /* Extract target reads from resolved read bundles, or pass all reads through when target enrichment is disabled. */

    tag "${meta.id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(meta), path(read_files, stageAs: "reads??????/*"), path(deacon_idx), val(target_enrichment_enabled)

    output:
    tuple val(meta.id), val(meta.platform), val(meta.deacon_read_structure), path("${meta.id}.target_enriched.fastq.gz"), emit: reads
    tuple val(meta.id), path("${meta.id}.deacon_filter.json"), emit: stats

    script:
    def files = read_files instanceof List ? read_files : [read_files]
    def r1_count = meta.r1_count as int
    def r1_files = files.take(r1_count)
    def r2_files = files.drop(r1_count)
    def deplete_arg = target_enrichment_enabled ? "" : "--deplete"
    if (meta.read_mode == "single" && r1_files.size() == 1)
        """
        deacon filter \
            ${deplete_arg} \
            --threads ${task.cpus} \
            --abs-threshold ${params.virus_abs_threshold} \
            --rel-threshold ${params.virus_rel_threshold} \
            --summary ${meta.id}.deacon_filter.json \
            --output ${meta.id}.target_enriched.fastq.gz \
            ${deacon_idx} \
            ${r1_files[0]}
        """
    else if (meta.read_mode == "paired" && r1_files.size() == 1 && r2_files.size() == 1)
        """
        deacon filter \
            ${deplete_arg} \
            --threads ${task.cpus} \
            --abs-threshold ${params.virus_abs_threshold} \
            --rel-threshold ${params.virus_rel_threshold} \
            --summary ${meta.id}.deacon_filter.json \
            --output ${meta.id}.target_enriched.fastq.gz \
            ${deacon_idx} \
            ${r1_files[0]} ${r2_files[0]}
        """
    else if (meta.read_mode == "single")
        """
        staged_reads=(reads??????/*)
        printf '%s\n' "\${staged_reads[@]}" > reads.list

        stream_fastqs_to_deacon.py \
            --sample-id ${meta.id} \
            --index ${deacon_idx} \
            ${deplete_arg} \
            --reads-list reads.list \
            --threads ${task.cpus} \
            --abs-threshold ${params.virus_abs_threshold} \
            --rel-threshold ${params.virus_rel_threshold} \
            --summary ${meta.id}.deacon_filter.json \
            --output ${meta.id}.target_enriched.fastq.gz
        """
    else
        """
        staged_reads=(reads??????/*)
        printf '%s\n' "\${staged_reads[@]:0:${r1_count}}" > r1.list
        printf '%s\n' "\${staged_reads[@]:${r1_count}}" > r2.list

        stream_fastqs_to_deacon.py \
            --sample-id ${meta.id} \
            --index ${deacon_idx} \
            ${deplete_arg} \
            --r1-list r1.list \
            --r2-list r2.list \
            --threads ${task.cpus} \
            --abs-threshold ${params.virus_abs_threshold} \
            --rel-threshold ${params.virus_rel_threshold} \
            --summary ${meta.id}.deacon_filter.json \
            --output ${meta.id}.target_enriched.fastq.gz
        """
}

process DEACON_FILTER_CONTIGS {
    /*
     * Retain target contigs using deacon filter on assembled FASTA.
     *
     * Replaces the 5-process STAT contig classification chain:
     * CLASSIFY_CONTIGS_FIRST_PASS → GENERATE_CONTIGS_TAXA_LIST →
     * CLASSIFY_CONTIGS_SECOND_PASS → IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS →
     * EXTRACT_HUMAN_VIRUS_CONTIGS
     *
     * Uses the same target index as read extraction so no additional reference
     * is needed. Contigs are already ≥200 bp (FILTER_SHORT_CONTIGS), so
     * incidental single k-mer matches are rare — the signal-to-noise ratio
     * improves vs reads because contigs have far more k-mers per sequence.
     */

    tag "${sample_id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(fasta), path(query_lookup), path(deacon_idx), val(contig_filter_policy), path(depletion_idx)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.target_filtered.fasta"), path(query_lookup), emit: contigs
    tuple val(sample_id), path("${sample_id}.contig_target_filtering.tsv"), emit: decisions

    script:
    def target_deplete_arg = contig_filter_policy.target_enrichment_enabled ? "" : "--deplete"
    if (contig_filter_policy.depletion_enabled)
        """
        set -euo pipefail

        deacon filter \
            ${target_deplete_arg} \
            --threads ${task.cpus} \
            --abs-threshold ${contig_filter_policy.target_abs_threshold} \
            --rel-threshold ${contig_filter_policy.target_rel_threshold} \
            --summary ${sample_id}.target_filter.deacon.json \
            ${deacon_idx} \
            ${fasta} \
        | deacon filter \
            --deplete \
            --threads ${task.cpus} \
            --abs-threshold ${contig_filter_policy.depletion_abs_threshold} \
            --rel-threshold ${contig_filter_policy.depletion_rel_threshold} \
            --summary ${sample_id}.host_depletion.deacon.json \
            --output ${sample_id}.target_filtered.fasta \
            ${depletion_idx} \
            -

        summarize_contig_filtering.py \
            --sample-id '${sample_id}' \
            --target-summary ${sample_id}.target_filter.deacon.json \
            --depletion-summary ${sample_id}.host_depletion.deacon.json \
            --output ${sample_id}.contig_target_filtering.tsv
        """
    else
        """
        deacon filter \
            ${target_deplete_arg} \
            --threads ${task.cpus} \
            --abs-threshold ${contig_filter_policy.target_abs_threshold} \
            --rel-threshold ${contig_filter_policy.target_rel_threshold} \
            --summary ${sample_id}.target_filter.deacon.json \
            --output ${sample_id}.target_filtered.fasta \
            ${deacon_idx} \
            ${fasta}

        summarize_contig_filtering.py \
            --sample-id '${sample_id}' \
            --target-summary ${sample_id}.target_filter.deacon.json \
            --output ${sample_id}.contig_target_filtering.tsv
        """
}
