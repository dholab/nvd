include { DEACON_BUILD_INDEX ; DEACON_FETCH_INDEX ; DEACON_UNION_INDEXES ; DEACON_DEPLETE } from "../modules/deacon"

workflow HOST_DEPLETION {
    /*
     * Orchestrates deacon-based host/contaminant depletion.
     *
     * Index resolution priority:
     * 1. params.deacon_index (explicit local path)
     * 2. params.deacon_index_url (download prebuilt, e.g., panhuman-1)
     * 3. params.deacon_contaminants_fasta (build from FASTA at runtime)
     *
     * If both a base index (local or URL) and deacon_contaminants_fasta are provided,
     * the indexes are unioned to combine panhuman with custom contaminants.
     *
     * Uses declarative channel ternaries so the DAG is consistent across runs,
     * matching the pattern in preprocess_reads.nf.
     */

    take:
    ch_reads  // tuple(sample_id, platform, read_structure, reads)

    main:
    // --- Resolve base index ---
    // Explicit local path takes priority over URL download.
    // Both produce a channel of one .idx file; unused path emits nothing.
    ch_local_index = params.deacon_index
        ? Channel.fromPath(params.deacon_index)
        : Channel.empty()

    ch_fetch_url = (!params.deacon_index && params.deacon_index_url)
        ? Channel.of(params.deacon_index_url)
        : Channel.empty()

    DEACON_FETCH_INDEX(ch_fetch_url)

    ch_base_index = ch_local_index.mix(DEACON_FETCH_INDEX.out.index)

    // --- Build custom index from FASTA if provided ---
    ch_custom_fasta = params.deacon_contaminants_fasta
        ? Channel.fromPath(params.deacon_contaminants_fasta)
        : Channel.empty()

    DEACON_BUILD_INDEX(ch_custom_fasta)
    ch_custom_index = DEACON_BUILD_INDEX.out.index

    // --- Combine indexes ---
    // Determine at parse time whether we need to union multiple indexes.
    // Both branches produce a value channel of one path (.idx file).
    def needs_union = params.deacon_contaminants_fasta && (params.deacon_index || params.deacon_index_url)

    DEACON_UNION_INDEXES(
        ch_base_index.mix(ch_custom_index).collect()
    )

    ch_index = needs_union
        ? DEACON_UNION_INDEXES.out.index.first()
        : ch_base_index.mix(ch_custom_index).first()

    // --- Run depletion ---
    DEACON_DEPLETE(ch_reads.combine(ch_index))

    emit:
    reads = DEACON_DEPLETE.out.reads
    stats = DEACON_DEPLETE.out.stats
}
