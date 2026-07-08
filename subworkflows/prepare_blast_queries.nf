include { DEACON_FILTER_CONTIGS } from "../modules/deacon"
include { SELECT_BLAST_QUERIES } from "../modules/blast"
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

    ch_query_batches = SELECT_BLAST_QUERIES.out.query_fastas
        .flatMap { sample_id, platform, fastas, lookup ->
            def batch_fastas = fastas instanceof List ? fastas : [fastas]
            batch_fastas.collect { batch_fasta ->
                def prefix = "${sample_id}."
                def suffix = ".blast_queries.fasta"
                def filename = batch_fasta.name
                assert filename.startsWith(prefix) && filename.endsWith(suffix)
                def evidence_class = filename.substring(prefix.length(), filename.length() - suffix.length())
                tuple(sample_id, platform, evidence_class, batch_fasta, lookup)
            }
        }

    emit:
    queries = ch_query_batches
    contig_sequences = DEACON_FILTER_CONTIGS.out
    contig_read_counts = CONTIG_READ_MAPBACK.out.contig_read_counts
    filtered_bam = CONTIG_READ_MAPBACK.out.filtered_bam
    unmapped_reads = CONTIG_READ_MAPBACK.out.unmapped_reads
    unmapped_read_counts = CONTIG_READ_MAPBACK.out.unmapped_read_counts
}
