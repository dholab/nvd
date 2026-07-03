include { RUN_SPADES } from "../modules/spades"

workflow SHORT_READ_DENOVO_ASSEMBLY {
    take:
    ch_short_read_fastqs  // tuple(sample_id, platform, read_structure, fastq)

    main:

    RUN_SPADES(ch_short_read_fastqs)

    emit:
    contigs = RUN_SPADES.out  // tuple(sample_id, platform, read_structure, fasta)
}
