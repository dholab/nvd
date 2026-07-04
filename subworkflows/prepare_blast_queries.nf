include { DEACON_FILTER_CONTIGS } from "../modules/deacon"
include { SELECT_BLAST_QUERIES } from "../modules/blast_queries"
include { MAP_READS_TO_CONTIGS } from "./map_reads_to_contigs"

workflow PREPARE_BLAST_QUERIES {
    take:
    ch_contigs          // tuple(sample_id, platform, read_structure, fasta, lookup) from PROCESS_CONTIGS
    ch_viral_reads      // tuple(sample_id, platform, read_structure, fastq) from PREPROCESS_READS
    ch_virus_index      // path: pre-built or freshly-built virus deacon index
    ch_depletion_index  // tuple(use_depletion, path): resolved host/contaminant depletion index or sentinel

    main:
    DEACON_FILTER_CONTIGS(
        ch_contigs
            .combine(ch_virus_index)
            .combine(ch_depletion_index)
    )

    MAP_READS_TO_CONTIGS(
        DEACON_FILTER_CONTIGS.out,
        ch_viral_reads,
    )

    SELECT_BLAST_QUERIES(
        DEACON_FILTER_CONTIGS.out,
    )

    emit:
    queries = SELECT_BLAST_QUERIES.out.queries
    contig_read_counts = MAP_READS_TO_CONTIGS.out.contig_read_counts
    filtered_bam = MAP_READS_TO_CONTIGS.out.filtered_bam
    unmapped_reads = MAP_READS_TO_CONTIGS.out.unmapped_reads
    unmapped_read_counts = MAP_READS_TO_CONTIGS.out.unmapped_read_counts
}
