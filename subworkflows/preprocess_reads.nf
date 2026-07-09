include { DEACON_FETCH_INDEX as DEACON_FETCH_VIRUS_INDEX  } from "../modules/deacon"
include { DEACON_FETCH_INDEX as DEACON_FETCH_HOST_INDEX   } from "../modules/deacon"
include { DEACON_BUILD_VIRUS_INDEX_FROM_FASTA             } from "../modules/deacon"
include { DEACON_BUILD_INDEX_FROM_FASTA                   } from "../modules/deacon"
include { DEACON_UNION_INDEXES                            } from "../modules/deacon"
include { DEACON_ENRICH_TARGET_READS                     } from "../modules/deacon"
include { DEACON_DEPLETE                                  } from "../modules/deacon"
include { MERGE_PAIRS ; DEDUP_WITH_CLUMPIFY ; TRIM_ADAPTERS ; FILTER_READS ; REPAIR_PAIRS } from "../modules/bbmap"
include { PROFILE_FASTX as PROFILE_READS } from "../modules/fastx"

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

    ch_virus_reads_by_layout = ch_virus_reads.branch { _id, platform, read_structure, _reads ->
        mergeable: platform == "illumina" && read_structure == "interleaved"
        other: true
    }

    ch_reads_to_merge = params.merge_pairs
        ? ch_virus_reads_by_layout.mergeable
        : channel.empty()

    MERGE_PAIRS(ch_reads_to_merge)

    ch_read_shards = params.merge_pairs
        ? MERGE_PAIRS.out.merged
            .mix(MERGE_PAIRS.out.unmerged)
            .mix(ch_virus_reads_by_layout.other.map { sample_id, platform, read_structure, reads ->
                tuple(sample_id, platform, read_structure, "single_read", reads)
            })
        : ch_virus_reads.map { sample_id, platform, read_structure, reads ->
            tuple(sample_id, platform, read_structure, "single_read", reads)
        }

    // 2a. Dedup
    def should_dedup_seq = params.dedup || params.dedup_seq
    ch_after_dedup = should_dedup_seq
        ? DEDUP_WITH_CLUMPIFY(ch_read_shards)
        : ch_read_shards

    // 2b. Adapter trim (Illumina only)
    ch_branched_for_trim = ch_after_dedup.branch { _id, platform, _rs, _ec, _reads ->
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

    // 2d. Independently optional quality/length and low-complexity filters
    ch_with_quality = ch_after_scrub.map { sample_id, platform, read_structure, query_class, reads ->
        def min_qual = platform == "illumina"
            ? params.min_read_quality_illumina
            : params.min_read_quality_nanopore
        tuple(sample_id, platform, read_structure, query_class, reads, min_qual)
    }
    def should_filter_reads = params.filter_reads || params.filter_low_complexity_reads
    ch_after_filter = should_filter_reads
        ? FILTER_READS(ch_with_quality)
        : ch_after_scrub

    // 2e. Repair pairs (interleaved only)
    ch_branched_for_repair = ch_after_filter.branch { _id, _p, read_structure, _ec, _r ->
        interleaved: read_structure == "interleaved"
        other: true
    }
    ch_repaired = REPAIR_PAIRS(ch_branched_for_repair.interleaved)
    ch_preprocessed_shards = ch_repaired.mix(ch_branched_for_repair.other)

    ch_preprocessed_shards_for_profile = ch_preprocessed_shards.map { sample_id, platform, read_structure, query_class, reads ->
        def min_qual = platform == "illumina"
            ? params.min_read_quality_illumina
            : params.min_read_quality_nanopore
        tuple(
            [
                id: sample_id,
                platform: platform,
                read_structure: read_structure,
                query_class: query_class,
                profile_stage: query_class,
                profile_key: "${sample_id}:${read_structure}:${query_class}",
                thresholds: [
                    [name: "min_read_length", axis: "length", value: params.min_read_length],
                    [name: "min_read_quality", axis: "quality", value: min_qual],
                ],
            ],
            reads,
        )
    }

    PROFILE_READS(ch_preprocessed_shards_for_profile)

    ch_preprocessed_shards_by_profile_key = ch_preprocessed_shards_for_profile.map { meta, reads ->
        tuple(meta.profile_key, reads)
    }

    ch_profiles_by_key = PROFILE_READS.out.profiled.map { meta, profile_json, length_histogram, sequence_count ->
        tuple(
            meta.profile_key,
            meta + [sequence_count: sequence_count.text.trim() as long],
            profile_json,
            length_histogram,
        )
    }

    ch_profiled_preprocessed_shards = ch_preprocessed_shards_by_profile_key
        .join(ch_profiles_by_key, by: 0)
        .map { _profile_key, reads, profiled_meta, profile_json, length_histogram ->
            tuple(
                profiled_meta,
                reads,
                profile_json,
                length_histogram,
            )
        }

    ch_mapback_groups = ch_profiled_preprocessed_shards
        .map { meta, reads, _profile_json, _length_histogram ->
            tuple(meta.id, meta.platform, [meta: meta, reads: reads])
        }
        .groupTuple(by: [0, 1])
        .map { sample_id, platform, rows ->
            tuple([id: sample_id, platform: platform], rows)
        }

    ch_mapback_by_layout = ch_mapback_groups.branch { _meta, rows ->
        paired: rows.any { row -> row.meta.read_structure == "interleaved" || row.meta.query_class == "overlap_merged_pair" }
        single: true
    }

    ch_paired_reads_for_mapback = ch_mapback_by_layout.paired.map { meta, rows ->
        def nonempty_rows = rows.findAll { row -> row.meta.sequence_count > 0 }
        def overlap_merged_pair_reads = nonempty_rows.findAll { row -> row.meta.query_class == "overlap_merged_pair" }.collect { row -> row.reads }
        def single_read_reads = nonempty_rows.findAll { row -> row.meta.query_class == "single_read" }.collect { row -> row.reads }
        tuple(
            meta.id,
            meta.platform,
            overlap_merged_pair_reads,
            single_read_reads,
        )
    }.filter { _sample_id, _platform, overlap_merged_pair_reads, single_read_reads ->
        overlap_merged_pair_reads.size() + single_read_reads.size() > 0
    }

    ch_single_reads_for_mapback = ch_mapback_by_layout.single.flatMap { meta, rows ->
        rows.collect { row -> tuple(meta.id, meta.platform, row.meta.read_structure, row.reads) }
    }

    emit:
    reads = ch_preprocessed_shards.map { sample_id, platform, read_structure, _query_class, reads ->
        tuple(sample_id, platform, read_structure, reads)
    }
    read_shards = ch_preprocessed_shards
    // tuple(meta, reads, profile_json, length_histogram)
    // meta: id, platform, read_structure, query_class, profile_stage, sequence_count
    profiled_reads = ch_profiled_preprocessed_shards.map { meta, reads, profile_json, length_histogram ->
        tuple(meta, reads, profile_json, length_histogram)
    }
    // tuple(meta, reads, profile_json, length_histogram)
    // meta: id, platform, read_structure, query_class, profile_stage, sequence_count
    profiled_read_shards = ch_profiled_preprocessed_shards
    paired_reads_for_mapback = ch_paired_reads_for_mapback
    single_reads_for_mapback = ch_single_reads_for_mapback
    read_counts = ch_read_counts
    virus_enrichment_stats = DEACON_ENRICH_TARGET_READS.out.stats
    virus_index = ch_virus_index
    depletion_index = ch_depletion_index_option
}
