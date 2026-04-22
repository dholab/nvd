include { DEACON_BUILD_INDEX_FROM_STAT_K_MERS ; DEACON_FILTER_HUMAN_VIRUS_READS } from "../modules/deacon"
include { MASK_LOW_COMPLEXITY ; FILTER_SHORT_CONTIGS } from "../modules/bbmap"
include { RUN_SPADES                } from "../modules/spades"

workflow PREPROCESS_CONTIGS {
    take:
    ch_sample_fastqs  // tuple(sample_id, platform, read_structure, fastq)
    ch_stat_dbss
    ch_stat_annotation
    ch_human_virus_taxlist

    main:

    // Build deacon index from STAT k-mers (runs once, not per-sample)
    DEACON_BUILD_INDEX_FROM_STAT_K_MERS(
        ch_stat_dbss,
        ch_stat_annotation,
        ch_human_virus_taxlist
    )

    // Extract human virus reads using deacon filter with the STAT-derived index
    DEACON_FILTER_HUMAN_VIRUS_READS(
        ch_sample_fastqs.combine(DEACON_BUILD_INDEX_FROM_STAT_K_MERS.out.index)
    )

    RUN_SPADES(
        DEACON_FILTER_HUMAN_VIRUS_READS.out
            .map { id, platform, read_structure, fq -> tuple(id, platform, read_structure, fq, file(fq).countFastq()) }
            .filter { _id, _platform, _read_structure, _fq, count -> count >= 100 }
            .map { id, platform, read_structure, fq, _count -> tuple(id, platform, read_structure, file(fq)) }
    )

    MASK_LOW_COMPLEXITY(
        RUN_SPADES.out
    )

    FILTER_SHORT_CONTIGS(
        MASK_LOW_COMPLEXITY.out
    )

    ch_filtered_contigs = FILTER_SHORT_CONTIGS.out
        .filter { _id, _platform, _read_structure, fasta ->
            def recordCount = fasta.countFasta()
            recordCount > 0
        }

    emit:
    contigs = ch_filtered_contigs  // tuple(sample_id, platform, read_structure, fasta)
    viral_reads = DEACON_FILTER_HUMAN_VIRUS_READS.out  // tuple(sample_id, platform, read_structure, fastq)
}
