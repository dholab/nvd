include { EXTRACT_HUMAN_VIRUS_READS } from "../modules/stat"
include { MASK_LOW_COMPLEXITY ; FILTER_SHORT_CONTIGS } from "../modules/bbmap"
include { RUN_SPADES                } from "../modules/spades"

workflow PREPROCESS_CONTIGS {
    take:
    ch_sample_fastqs
    ch_stat_dbss
    ch_human_virus_taxlist

    main:
    EXTRACT_HUMAN_VIRUS_READS(
        ch_sample_fastqs.combine(ch_stat_dbss).combine(ch_human_virus_taxlist)
    )

    RUN_SPADES(
        EXTRACT_HUMAN_VIRUS_READS.out
    )

    MASK_LOW_COMPLEXITY(
        RUN_SPADES.out
    )

    FILTER_SHORT_CONTIGS(
        MASK_LOW_COMPLEXITY.out
    )

    emit:
    FILTER_SHORT_CONTIGS.out
}
