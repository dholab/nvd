include { DEACON_FILTER_CONTIGS } from "../modules/deacon"
include { SELECT_BLAST_QUERIES } from "../modules/blast_queries"
include { CONTIG_READ_MAPBACK } from "./contig_read_mapback"

workflow PREPARE_BLAST_QUERIES {
    take:
    ch_contigs          // tuple(sample_id, platform, read_structure, fasta, lookup) from PROCESS_CONTIGS
    ch_paired_reads     // tuple(sample_id, platform, overlap_merged_pair_reads, single_read_reads) from PREPROCESS_READS
    ch_single_reads     // tuple(sample_id, platform, read_structure, fastq) from PREPROCESS_READS
    ch_virus_index      // path: pre-built or freshly-built virus deacon index
    ch_depletion_index  // tuple(use_depletion, path): resolved host/contaminant depletion index or sentinel

    main:
    DEACON_FILTER_CONTIGS(
        ch_contigs
            .combine(ch_virus_index)
            .combine(ch_depletion_index)
    )

    CONTIG_READ_MAPBACK(
        DEACON_FILTER_CONTIGS.out,
        ch_paired_reads,
        ch_single_reads,
    )

    SELECT_BLAST_QUERIES(
        DEACON_FILTER_CONTIGS.out,
    )

    emit:
    queries = SELECT_BLAST_QUERIES.out.queries
    contig_read_counts = CONTIG_READ_MAPBACK.out.contig_read_counts
    filtered_bam = CONTIG_READ_MAPBACK.out.filtered_bam
    unmapped_reads = CONTIG_READ_MAPBACK.out.unmapped_reads
    unmapped_read_counts = CONTIG_READ_MAPBACK.out.unmapped_read_counts
}
