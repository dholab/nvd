include { RUN_SPADES } from "../modules/spades"

workflow LONG_READ_DENOVO_ENSEMBLY {
    take:
    ch_long_read_fastqs  // tuple(sample_id, platform, read_structure, fastq)

    main:

    // Transitional implementation: SPAdes stays behind this long-read seam
    // until the ensemble assembler design replaces it.
    RUN_SPADES(ch_long_read_fastqs)

    emit:
    contigs = RUN_SPADES.out  // tuple(sample_id, platform, read_structure, fasta)
}
