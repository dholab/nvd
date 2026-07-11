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

    // Pair merging is where one sample can split into two read batches. We keep
    // those batches separate because later mapback and BLAST results distinguish
    // reads from successfully merged pairs from reads that did not merge. The
    // query_class name records that distinction. sample_batch_count tells later
    // code how many post-QC batches to wait for; an empty output still counts as
    // long as its task completed successfully. MERGE_PAIRS always creates both
    // an overlap_merged_pair batch and a single_read batch, hence two here.
    ch_merged_read_batches = MERGE_PAIRS.out.merged.map { sample_id, platform, read_structure, query_class, reads ->
        tuple(
            [
                id: sample_id,
                platform: platform,
                read_structure: read_structure,
                query_class: query_class,
                sample_batch_count: 2,
            ],
            reads,
        )
    }
    // The unmerged reads are the second query-class batch from the same sample.
    ch_unmerged_read_batches = MERGE_PAIRS.out.unmerged.map { sample_id, platform, read_structure, query_class, reads ->
        tuple(
            [
                id: sample_id,
                platform: platform,
                read_structure: read_structure,
                query_class: query_class,
                sample_batch_count: 2,
            ],
            reads,
        )
    }
    // Inputs that do not enter MERGE_PAIRS remain one single_read batch.
    ch_nonmergeable_read_batches = ch_virus_reads_by_layout.other.map { sample_id, platform, read_structure, reads ->
        tuple(
            [
                id: sample_id,
                platform: platform,
                read_structure: read_structure,
                query_class: "single_read",
                sample_batch_count: 1,
            ],
            reads,
        )
    }

    // With pair merging disabled, every sample likewise has one single_read batch.
    ch_read_batches = params.merge_pairs
        ? ch_merged_read_batches
            .mix(ch_unmerged_read_batches)
            .mix(ch_nonmergeable_read_batches)
        : ch_virus_reads.map { sample_id, platform, read_structure, reads ->
            tuple(
                [
                    id: sample_id,
                    platform: platform,
                    read_structure: read_structure,
                    query_class: "single_read",
                    sample_batch_count: 1,
                ],
                reads,
            )
        }

    // 2a. Dedup
    def should_dedup_seq = params.dedup || params.dedup_seq
    ch_after_dedup = should_dedup_seq
        ? DEDUP_WITH_CLUMPIFY(ch_read_batches)
        : ch_read_batches

    // 2b. Adapter trim (Illumina only)
    ch_branched_for_trim = ch_after_dedup.branch { meta, _reads ->
        illumina: meta.platform == "illumina"
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
    def should_filter_reads = params.filter_reads || params.filter_low_complexity_reads
    ch_after_filter = should_filter_reads
        ? FILTER_READS(ch_after_scrub)
        : ch_after_scrub

    // 2e. Repair pairs (interleaved only)
    ch_branched_for_repair = ch_after_filter.branch { meta, _reads ->
        interleaved: meta.read_structure == "interleaved"
        other: true
    }
    ch_repaired = REPAIR_PAIRS(ch_branched_for_repair.interleaved)
    ch_preprocessed_batches = ch_repaired.mix(ch_branched_for_repair.other)

    ch_preprocessed_batches_for_profile = ch_preprocessed_batches.map { meta, reads ->
        def min_qual = meta.platform == "illumina"
            ? params.min_read_quality_illumina
            : params.min_read_quality_nanopore
        tuple(
            meta + [
                profile_stage: meta.query_class,
                profile_key: "${meta.id}:${meta.read_structure}:${meta.query_class}",
                thresholds: [
                    [name: "min_read_length", axis: "length", value: params.min_read_length],
                    [name: "min_read_quality", axis: "quality", value: min_qual],
                ],
            ],
            reads,
        )
    }

    PROFILE_READS(ch_preprocessed_batches_for_profile)

    ch_preprocessed_batches_by_profile_key = ch_preprocessed_batches_for_profile.map { meta, reads ->
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

    ch_profiled_preprocessed_batches = ch_preprocessed_batches_by_profile_key
        .join(ch_profiles_by_key, by: 0)
        .map { _profile_key, reads, profiled_meta, profile_json, length_histogram ->
            tuple(
                profiled_meta,
                reads,
                profile_json,
                length_histogram,
            )
        }

    // QC runs independently for each query-class batch, so batches from faster
    // samples should not wait for preprocessing to finish across the whole run.
    // An ordinary groupTuple cannot tell when one sample is complete and waits
    // for its input channel to close. groupKey supplies the known count of one
    // or two batches, allowing each complete sample to continue immediately to
    // assembly, mapback, and sourmash. Keeping QC batch-based preserves its
    // existing parallelism and cache reuse; grouping once here also ensures all
    // three consumers see the same available post-QC reads. If an ignored task
    // drops one batch, remainder keeps the surviving batch for best-effort work
    // after the upstream channel closes.
    ch_profiled_batches_by_sample = ch_profiled_preprocessed_batches
        .map { meta, reads, _profile_json, _length_histogram ->
            tuple(
                groupKey(
                    [meta.id, meta.platform, meta.read_structure],
                    meta.sample_batch_count as int,
                ),
                [meta: meta, reads: reads],
            )
        }
        .groupTuple(remainder: true)
        .map { sample_key, batches ->
            def sample = sample_key.getGroupTarget()

            // A missing batch is expected when an upstream task exhausts its
            // retries and is ignored; that case is handled below. Conflicting
            // counts, duplicate query classes, or extra batches instead mean
            // the workflow wired a sample incorrectly. Continuing would risk
            // double-counting reads, so these programming errors stop the run.
            def expected_batch_counts = batches.collect { batch -> batch.meta.sample_batch_count as int }.unique()
            assert expected_batch_counts.size() == 1 : "Sample ${sample} has inconsistent expected batch counts: ${expected_batch_counts}"

            def query_classes = batches.collect { batch -> batch.meta.query_class }
            assert query_classes.unique().size() == query_classes.size() : "Sample ${sample} has duplicate query-class batches: ${query_classes}"

            def expected_batch_count = expected_batch_counts.first()
            def available_batch_count = batches.size()
            assert available_batch_count <= expected_batch_count : "Sample ${sample} produced ${available_batch_count} batches but expected ${expected_batch_count}"

            if (available_batch_count < expected_batch_count) {
                log.warn "Processing available post-QC batches for sample ${sample}: expected ${expected_batch_count}, received ${available_batch_count} (${query_classes.sort()})"
            }

            // Arrival order depends on task timing. Sort before passing files
            // downstream so concatenation order and task hashes stay stable.
            def ordered_batches = batches.sort { left, right -> left.meta.query_class <=> right.meta.query_class }
            // Downstream eligibility checks use the total reads available for
            // this sample, including all query classes that survived QC.
            def sequence_count = ordered_batches.collect { batch -> batch.meta.sequence_count }.sum()
            // Keep the per-batch metadata and files, but provide shared sample
            // metadata once for assembly, mapback, and sourmash.
            tuple(
                [
                    id: sample[0],
                    platform: sample[1],
                    read_structure: sample[2],
                    sequence_count: sequence_count,
                ],
                ordered_batches,
            )
        }

    ch_mapback_by_layout = ch_profiled_batches_by_sample.branch { _meta, batches ->
        paired: batches.any { batch -> batch.meta.read_structure == "interleaved" || batch.meta.query_class == "overlap_merged_pair" }
        single: true
    }

    ch_paired_reads_for_mapback = ch_mapback_by_layout.paired.map { meta, batches ->
        def nonempty_batches = batches.findAll { batch -> batch.meta.sequence_count > 0 }
        def overlap_merged_pair_reads = nonempty_batches.findAll { batch -> batch.meta.query_class == "overlap_merged_pair" }.collect { batch -> batch.reads }
        def single_read_reads = nonempty_batches.findAll { batch -> batch.meta.query_class == "single_read" }.collect { batch -> batch.reads }
        tuple(
            meta.id,
            meta.platform,
            overlap_merged_pair_reads,
            single_read_reads,
        )
    }.filter { _sample_id, _platform, overlap_merged_pair_reads, single_read_reads ->
        overlap_merged_pair_reads.size() + single_read_reads.size() > 0
    }

    ch_single_reads_for_mapback = ch_mapback_by_layout.single.flatMap { meta, batches ->
        batches
            .findAll { batch -> batch.meta.sequence_count > 0 }
            .collect { batch -> tuple(meta.id, meta.platform, batch.meta.read_structure, batch.reads) }
    }

    emit:
    reads = ch_preprocessed_batches.map { meta, reads ->
        tuple(meta.id, meta.platform, meta.read_structure, reads)
    }
    read_batches = ch_preprocessed_batches.map { meta, reads ->
        tuple(meta.id, meta.platform, meta.read_structure, meta.query_class, reads)
    }
    // tuple(meta, reads, profile_json, length_histogram)
    // meta: id, platform, read_structure, query_class, sample_batch_count, profile_stage, sequence_count
    profiled_read_batches = ch_profiled_preprocessed_batches
    // tuple(sample_meta, batches); each batch contains meta and reads
    profiled_batches_by_sample = ch_profiled_batches_by_sample
    paired_reads_for_mapback = ch_paired_reads_for_mapback
    single_reads_for_mapback = ch_single_reads_for_mapback
    read_counts = ch_read_counts
    virus_enrichment_stats = DEACON_ENRICH_TARGET_READS.out.stats
    virus_index = ch_virus_index
    depletion_index = ch_depletion_index_option
}
