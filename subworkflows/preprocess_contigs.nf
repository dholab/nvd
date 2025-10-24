include { EXTRACT_HUMAN_VIRUS_READS } from "../modules/stat"
include { MASK_LOW_COMPLEXITY ; FILTER_SHORT_CONTIGS } from "../modules/bbmap"
include { RUN_SPADES                } from "../modules/spades"

workflow PREPROCESS_CONTIGS {
    take:
    ch_sample_fastqs
    ch_stat_dbss
    ch_stat_annotation
    ch_human_virus_taxlist

    main:
    
    // Uses STAT to get on fastq records that map to specific virus infecting virueses
    EXTRACT_HUMAN_VIRUS_READS(
        ch_sample_fastqs.combine(ch_stat_dbss).combine(ch_stat_annotation).combine(ch_human_virus_taxlist)
    )

    RUN_SPADES(
        EXTRACT_HUMAN_VIRUS_READS.out
            .map { id, platform, fq -> tuple(id, platform, fq, file(fq).countFastq()) }
            .filter { _id, _platform, _fq, count -> count > 0 }
            .map { id, platform, fq, _count -> tuple(id, platform, file(fq)) }
    )

    MASK_LOW_COMPLEXITY(
        RUN_SPADES.out
    )

    FILTER_SHORT_CONTIGS(
        MASK_LOW_COMPLEXITY.out
    )

    ch_filtered_contigs = FILTER_SHORT_CONTIGS.out
        .filter { _id, _platform, fasta ->
            def recordCount = fasta.countFasta()
            recordCount > 0

        }

    emit:
    contigs = ch_filtered_contigs
    viral_reads = EXTRACT_HUMAN_VIRUS_READS.out
}
