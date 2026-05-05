include { MASK_LOW_COMPLEXITY ; FILTER_SHORT_CONTIGS } from "../modules/bbmap"
include { RUN_SPADES                } from "../modules/spades"

workflow PREPROCESS_CONTIGS {
    take:
    ch_sample_fastqs  // tuple(sample_id, platform, read_structure, fastq) — already virus-filtered and preprocessed

    main:

    // Assemble virus reads into contigs with SPAdes.
    // Filter out samples with fewer than 100 reads (insufficient for assembly).
    RUN_SPADES(
        ch_sample_fastqs
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
    viral_reads = ch_sample_fastqs  // the incoming reads ARE the viral reads (extracted upstream)
}
