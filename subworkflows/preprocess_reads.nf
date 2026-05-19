include { DEACON_FETCH_INDEX as DEACON_FETCH_VIRUS_INDEX  } from "../modules/deacon"
include { DEACON_FETCH_INDEX as DEACON_FETCH_HOST_INDEX   } from "../modules/deacon"
include { DEACON_BUILD_VIRUS_INDEX_FROM_FASTA             } from "../modules/deacon"
include { DEACON_BUILD_INDEX_FROM_FASTA                   } from "../modules/deacon"
include { DEACON_UNION_INDEXES                            } from "../modules/deacon"
include { DEACON_FILTER_HUMAN_VIRUS_READS                 } from "../modules/deacon"
include { DEACON_DEPLETE                                  } from "../modules/deacon"
include { DEDUP_WITH_CLUMPIFY ; TRIM_ADAPTERS ; FILTER_READS ; REPAIR_PAIRS } from "../modules/bbmap"

workflow PREPROCESS_READS {
    take:
    ch_sample_fastqs  // Queue channel: pre-interleave tuples from GATHER_READS
                      // Paired: tuple(sample_id, platform, R1, R2)
                      // Single: tuple(sample_id, platform, fastq)

    main:
    // Normalize mixed-size tuples from GATHER_READS:
    //   Paired: (id, platform, R1, R2) -> (id, platform, R1, R2)
    //   Single: (id, platform, fastq)  -> (id, platform, fastq, NO_R2)
    // The sentinel file NO_R2 lets DEACON_FILTER_HUMAN_VIRUS_READS distinguish
    // paired from single-end input with a fixed-size tuple.
    ch_normalized_samples = ch_sample_fastqs
        .map { items ->
            if (items.size() == 4)
                tuple(items[0], items[1], file(items[2]), file(items[3]))
            else
                tuple(items[0], items[1], file(items[2]), file("NO_R2"))
        }

    // -------------------------------------------------------------------------
    // Step 1: Resolve virus index and frontloaded extraction
    // -------------------------------------------------------------------------
    // Priority: explicit local path → URL download → build from reference FASTA.
    // Each condition guards against the higher-priority source being set, so at
    // most one channel is non-empty and the mix passes through exactly one index.
    ch_local_virus_index = params.virus_index
        ? Channel.fromPath(params.virus_index)
        : Channel.empty()
    ch_virus_fetch_url = (!params.virus_index && params.virus_index_url)
        ? Channel.of(params.virus_index_url)
        : Channel.empty()
    ch_virus_ref_fasta = (!params.virus_index && !params.virus_index_url && params.virus_reference_fasta)
        ? Channel.fromPath(params.virus_reference_fasta)
        : Channel.empty()

    DEACON_FETCH_VIRUS_INDEX(ch_virus_fetch_url)
    DEACON_BUILD_VIRUS_INDEX_FROM_FASTA(ch_virus_ref_fasta)

    ch_virus_index = ch_local_virus_index
        .mix(DEACON_FETCH_VIRUS_INDEX.out.index)
        .mix(DEACON_BUILD_VIRUS_INDEX_FROM_FASTA.out.index)

    // Extract virus reads — runs BEFORE any preprocessing. For paired reads,
    // deacon takes R1/R2 and outputs interleaved FASTQ in one step.
    DEACON_FILTER_HUMAN_VIRUS_READS(
        ch_normalized_samples.combine(ch_virus_index)
    )

    // -------------------------------------------------------------------------
    // Step 2: Inlined preprocessing on virus-only reads
    // -------------------------------------------------------------------------
    ch_virus_reads = DEACON_FILTER_HUMAN_VIRUS_READS.out.reads

    // Extract total read counts from deacon summary JSON (replaces COUNT_READS).
    // The seqs_in field is the total input read count across R1+R2.
    ch_read_counts = DEACON_FILTER_HUMAN_VIRUS_READS.out.stats
        .map { sample_id, json_file ->
            def summary = new groovy.json.JsonSlurper().parse(json_file.toFile())
            tuple(sample_id, summary.seqs_in.toString())
        }

    // 2a. Dedup
    def should_dedup_seq = params.dedup || params.dedup_seq
    ch_after_dedup = should_dedup_seq
        ? DEDUP_WITH_CLUMPIFY(ch_virus_reads)
        : ch_virus_reads

    // 2b. Adapter trim (Illumina only)
    ch_branched_for_trim = ch_after_dedup.branch { _id, platform, _rs, _reads ->
        illumina: platform == "illumina"
        other: true
    }
    ch_after_trim = params.trim_adapters
        ? TRIM_ADAPTERS(ch_branched_for_trim.illumina).mix(ch_branched_for_trim.other)
        : ch_after_dedup

    // 2c. Host/contaminant depletion with deacon (optional)
    def has_host_config = params.host_index || params.host_index_url || params.host_contaminants_fasta
    if (has_host_config) {
        ch_local_host_index = params.host_index
            ? Channel.fromPath(params.host_index)
            : Channel.empty()
        ch_host_fetch_url = (!params.host_index && params.host_index_url)
            ? Channel.of(params.host_index_url)
            : Channel.empty()
        ch_host_contaminants_fasta = params.host_contaminants_fasta
            ? Channel.fromPath(params.host_contaminants_fasta)
            : Channel.empty()

        DEACON_FETCH_HOST_INDEX(ch_host_fetch_url)
        DEACON_BUILD_INDEX_FROM_FASTA(ch_host_contaminants_fasta)

        ch_host_index_sources = ch_local_host_index
            .mix(DEACON_FETCH_HOST_INDEX.out.index)
            .mix(DEACON_BUILD_INDEX_FROM_FASTA.out.index)
            .collect()

        DEACON_UNION_INDEXES(ch_host_index_sources)
        ch_host_index = DEACON_UNION_INDEXES.out.index

        ch_after_scrub = DEACON_DEPLETE(ch_after_trim.combine(ch_host_index)).reads
    } else {
        ch_after_scrub = ch_after_trim
    }

    // 2d. Quality/length filter
    ch_with_quality = ch_after_scrub.map { sample_id, platform, read_structure, reads ->
        def min_qual = platform == "illumina"
            ? params.min_read_quality_illumina
            : params.min_read_quality_nanopore
        tuple(sample_id, platform, read_structure, reads, min_qual)
    }
    ch_after_filter = params.filter_reads
        ? FILTER_READS(ch_with_quality)
        : ch_after_scrub

    // 2e. Repair pairs (interleaved only)
    ch_branched_for_repair = ch_after_filter.branch { _id, _p, read_structure, _r ->
        interleaved: read_structure == "interleaved"
        other: true
    }
    ch_repaired = REPAIR_PAIRS(ch_branched_for_repair.interleaved)
    ch_preprocessed = ch_repaired.mix(ch_branched_for_repair.other)

    emit:
    reads = ch_preprocessed
    read_counts = ch_read_counts
    virus_index = ch_virus_index
}
