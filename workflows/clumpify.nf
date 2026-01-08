include { CLUMP_READS       } from "../modules/bbmap"
include { SCRUB_HUMAN_READS } from "../modules/stat"

workflow CLUMPIFY_WORKFLOW {
    take:
    ch_gathered_reads  // tuple(sample_id, platform, read_structure, fastq)
    ch_start_gate

    main:
    def gated_reads = ch_gathered_reads
        .combine(ch_start_gate)
        .map { it[0..-2] }

    // Add human read scrubbing step for public posting of clumpified data to SRA
    CLUMP_READS(
        gated_reads.map { id, _platform, _read_structure, reads -> tuple(id, file(reads)) }
    )

    // Resolve human database path, supporting deprecated human_read_scrub param
    def human_db_path = params.sra_human_db ?: params.human_read_scrub

    if (params.human_read_scrub != null && params.sra_human_db == null) {
        log.warn "DEPRECATION WARNING: --human_read_scrub is deprecated. Please use --sra_human_db instead."
    }

    // Check if human database path is a valid file
    if (human_db_path != null && !file(human_db_path).isFile()) {
        error("Error: Human database file does not exist: ${human_db_path}")
    }

    ch_human_reads = human_db_path ? Channel.fromPath( human_db_path ) : Channel.empty()

    SCRUB_HUMAN_READS(
        CLUMP_READS.out.combine(ch_human_reads)
    )
}
