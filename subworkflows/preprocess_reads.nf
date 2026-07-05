include { DEACON_FETCH_INDEX as DEACON_FETCH_VIRUS_INDEX  } from "../modules/deacon"
include { DEACON_FETCH_INDEX as DEACON_FETCH_HOST_INDEX   } from "../modules/deacon"
include { DEACON_BUILD_VIRUS_INDEX_FROM_FASTA             } from "../modules/deacon"
include { DEACON_BUILD_INDEX_FROM_FASTA                   } from "../modules/deacon"
include { DEACON_UNION_INDEXES                            } from "../modules/deacon"
include { DEACON_ENRICH_TARGET_READS                     } from "../modules/deacon"
include { DEACON_DEPLETE                                  } from "../modules/deacon"
include { DEDUP_WITH_CLUMPIFY ; TRIM_ADAPTERS ; FILTER_READS ; REPAIR_PAIRS } from "../modules/bbmap"

workflow PREPROCESS_READS {
    take:
    ch_read_bundles  // tuple(meta, read_files) from GATHER_READS

    main:

    // -------------------------------------------------------------------------
    // Step 1: Resolve target index and frontloaded extraction
    // -------------------------------------------------------------------------
    // Priority when target enrichment is enabled: explicit local path → URL
    // download → build from reference FASTA. When target enrichment is disabled,
    // use a committed empty index in deplete mode so deacon keeps all records
    // while preserving the existing one-FASTQ downstream contract.
    def target_enrichment_enabled = NvdUtils.targetEnrichmentEnabled(params)
    ch_target_enrichment_enabled = Channel.value(target_enrichment_enabled)
    ch_local_virus_index = target_enrichment_enabled && params.virus_index
        ? Channel.fromPath(params.virus_index)
        : Channel.empty()
    ch_virus_fetch_url = (target_enrichment_enabled && !params.virus_index && params.virus_index_url)
        ? Channel.of(params.virus_index_url)
        : Channel.empty()
    ch_virus_ref_fasta = (target_enrichment_enabled && !params.virus_index && !params.virus_index_url && params.virus_reference_fasta)
        ? Channel.fromPath(params.virus_reference_fasta)
        : Channel.empty()
    ch_empty_virus_index = target_enrichment_enabled
        ? Channel.empty()
        : Channel.fromPath("${projectDir}/assets/empty_deacon.k31w1.idx")

    DEACON_FETCH_VIRUS_INDEX(ch_virus_fetch_url)
    DEACON_BUILD_VIRUS_INDEX_FROM_FASTA(ch_virus_ref_fasta)

    ch_virus_index = ch_local_virus_index
        .mix(DEACON_FETCH_VIRUS_INDEX.out.index)
        .mix(DEACON_BUILD_VIRUS_INDEX_FROM_FASTA.out.index)
        .mix(ch_empty_virus_index)

    // Extract target reads — runs BEFORE any preprocessing. For paired reads,
    // deacon takes R1/R2 and outputs interleaved FASTQ in one step. When target
    // enrichment is disabled, the empty-index deplete mode retains all reads.
    DEACON_ENRICH_TARGET_READS(
        ch_read_bundles.combine(ch_virus_index)
            .combine(ch_target_enrichment_enabled)
    )

    // -------------------------------------------------------------------------
    // Step 2: Inlined preprocessing on virus-only reads
    // -------------------------------------------------------------------------
    ch_virus_reads = DEACON_ENRICH_TARGET_READS.out.reads

    // Extract total read counts from deacon summary JSON (replaces COUNT_READS).
    // The seqs_in field is the total input read count across R1+R2.
    ch_read_counts = DEACON_ENRICH_TARGET_READS.out.stats
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

    // 2c. Host/contaminant depletion with deacon (optional). The public
    // parameter names remain host_* for compatibility, but this channel is the
    // general depletion index used by both read and contig filtering.
    def has_depletion_config = params.host_index || params.host_index_url || params.host_contaminants_fasta
    if (has_depletion_config) {
        ch_local_depletion_index = params.host_index
            ? Channel.fromPath(params.host_index)
            : Channel.empty()
        ch_depletion_fetch_url = (!params.host_index && params.host_index_url)
            ? Channel.of(params.host_index_url)
            : Channel.empty()
        ch_depletion_contaminants_fasta = params.host_contaminants_fasta
            ? Channel.fromPath(params.host_contaminants_fasta)
            : Channel.empty()

        DEACON_FETCH_HOST_INDEX(ch_depletion_fetch_url)
        DEACON_BUILD_INDEX_FROM_FASTA(ch_depletion_contaminants_fasta)

        ch_depletion_index_sources = ch_local_depletion_index
            .mix(DEACON_FETCH_HOST_INDEX.out.index)
            .mix(DEACON_BUILD_INDEX_FROM_FASTA.out.index)
            .collect()

        DEACON_UNION_INDEXES(ch_depletion_index_sources)
        ch_depletion_index = DEACON_UNION_INDEXES.out.index
        ch_depletion_index_option = ch_depletion_index.map { idx -> tuple(true, idx) }

        ch_after_scrub = DEACON_DEPLETE(ch_after_trim.combine(ch_depletion_index)).reads
    } else {
        ch_depletion_index_option = Channel.value(tuple(false, file("${projectDir}/assets/README.md")))
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
    virus_enrichment_stats = DEACON_ENRICH_TARGET_READS.out.stats
    virus_index = ch_virus_index
    depletion_index = ch_depletion_index_option
}
