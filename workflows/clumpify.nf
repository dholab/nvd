include { CLUMP_READS       } from "../modules/bbmap"
include { SCRUB_HUMAN_READS } from "../modules/stat"

workflow CLUMPIFY_WORKFLOW {
    take:
    ch_gathered_reads
    ch_start_gate

    main:
    def gated_reads = ch_gathered_reads
        .combine(ch_start_gate)
        .map { it[0..-2] }

    // Add human read scrubbing step for public posting of clumpified data to SRA
    CLUMP_READS(
        gated_reads.map { id, _platform, reads -> tuple(id, file(reads)) }
    )

    // Check if human_read_scrub is a valid file path
    human_db_file = file(params.human_read_scrub)

    if (!human_db_file.exists()) {
        error("Error: Human database file does not exist: ${params.human_read_scrub}")
    }

    if (!human_db_file.isFile()) {
        error("Error: Human database path is not a file: ${params.human_read_scrub}")
    }

    SCRUB_HUMAN_READS(
        CLUMP_READS.out
    )

    SCRUB_HUMAN_READS.out
    .map { sample_id, fastq -> 
        tuple(sample_id, "illumina", fastq)
    }
    .set { scrubbed_reads_with_platform }

    emit:
    scrubbed_reads_with_platform
}

// next process takes
// Queue channel of sample IDs, platforms, and (interleaved) FASTQ files: tuple val(sample_id), val(platform), path(fastq)