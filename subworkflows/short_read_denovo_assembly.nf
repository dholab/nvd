include { ASSEMBLE_WITH_SPADES } from "../modules/spades"

workflow SHORT_READ_DENOVO_ASSEMBLY {
    take:
    ch_short_read_fastqs  // tuple(sample_id, platform, read_structure, fastq)

    main:

    ASSEMBLE_WITH_SPADES(ch_short_read_fastqs)

    emit:
    contigs = ASSEMBLE_WITH_SPADES.out  // tuple(sample_id, platform, read_structure, producer, fasta)
}
