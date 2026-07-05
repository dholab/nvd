include { DEACON_FILTER_CONTIGS } from "../modules/deacon"
include { MAP_READS_TO_CONTIGS  } from "../modules/minimap2"
include { COUNT_MAPPED_READS    } from "../modules/samtools"

workflow EXTRACT_HUMAN_VIRUSES {
    take:
    ch_filtered_contigs   // tuple(sample_id, platform, read_structure, fasta) from PREPROCESS_CONTIGS
    ch_viral_reads        // tuple(sample_id, platform, read_structure, fastq) from DEACON_FILTER_HUMAN_VIRUS_READS
    ch_virus_index        // path: pre-built or freshly-built virus deacon index
    ch_depletion_index    // tuple(use_depletion, path): resolved host/contaminant depletion index or sentinel

    main:
    // Filter assembled contigs to those sharing k-mers with the virus index.
    // Collapses the former 5-process STAT classification chain into one step.
    DEACON_FILTER_CONTIGS(
        ch_filtered_contigs
            .combine(ch_virus_index)
            .combine(ch_depletion_index)
    )

    // Align virus reads back to the deacon-filtered virus contigs to get
    // per-contig read depth. Join by sample_id so each sample's reads
    // are mapped to its own contigs.
    ch_reads_with_contigs = ch_viral_reads
        .map { sample_id, platform, _read_structure, reads -> tuple(sample_id, platform, reads) }
        .join(
            DEACON_FILTER_CONTIGS.out.map { sample_id, contigs -> tuple(sample_id, contigs) },
            by: 0
        )
        // Result: [sample_id, platform, reads, contigs]

    MAP_READS_TO_CONTIGS(ch_reads_with_contigs)

    COUNT_MAPPED_READS(MAP_READS_TO_CONTIGS.out)

    emit:
    contigs = DEACON_FILTER_CONTIGS.out
    contig_read_counts = COUNT_MAPPED_READS.out.mapped_counts
}
