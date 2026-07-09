include { DEACON_FILTER_CONTIGS } from "../modules/deacon"
include { SELECT_BLAST_QUERIES ; NORMALIZE_READ_BLAST_QUERIES ; SUMMARIZE_BLAST_QUERY_BATCHES ; SUMMARIZE_EMPTY_BLAST_QUERY_BATCHES } from "../modules/blast"
include { PROFILE_FASTX as PROFILE_FILTERED_CONTIGS } from "../modules/fastx"
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
                def query_class = filename.substring(prefix.length(), filename.length() - suffix.length())
                tuple(sample_id, platform, query_class, batch_fasta, lookup)
            }
        }

    ch_contig_query_lookups = DEACON_FILTER_CONTIGS.out
        .map { sample_id, _platform, _read_structure, _fasta, lookup -> tuple(sample_id, lookup) }

    if (params.experimental == true) {
        NORMALIZE_READ_BLAST_QUERIES(CONTIG_READ_MAPBACK.out.unmapped_reads)
        ch_read_query_batches = NORMALIZE_READ_BLAST_QUERIES.out.queries
            .flatMap { sample_id, platform, query_class, fastas, lookups ->
                def batch_fastas = fastas ? (fastas instanceof List ? fastas : [fastas]) : []
                def query_lookups = lookups ? (lookups instanceof List ? lookups : [lookups]) : []
                assert batch_fastas.size() == query_lookups.size()
                batch_fastas.withIndex().collect { batch_fasta, index ->
                    tuple(sample_id, platform, query_class, batch_fasta, query_lookups[index])
                }
            }
        ch_query_batches = ch_query_batches.mix(ch_read_query_batches)
        ch_query_lookups = ch_contig_query_lookups
            .mix(ch_read_query_batches.map { sample_id, _platform, _query_class, _fasta, lookup -> tuple(sample_id, lookup) })
            .groupTuple()
        SUMMARIZE_BLAST_QUERY_BATCHES(
            ch_query_batches.groupTuple(by: [0, 1])
        )
        ch_samples_without_query_batches = DEACON_FILTER_CONTIGS.out
            .map { sample_id, platform, _read_structure, _fasta, _lookup -> tuple(sample_id, platform, "sample") }
            .unique()
            .join(
                ch_query_batches
                    .map { sample_id, platform, _query_class, _fasta, _lookup -> tuple(sample_id, platform, "batch") }
                    .unique(),
                by: [0, 1],
                remainder: true,
            )
            .filter { row -> row.size() == 3 }
            .map { sample_id, platform, _sample_marker -> tuple(sample_id, platform) }
        SUMMARIZE_EMPTY_BLAST_QUERY_BATCHES(ch_samples_without_query_batches)
    } else {
        ch_query_lookups = ch_contig_query_lookups.groupTuple()
    }

    emit:
    queries = ch_query_batches
    query_lookups = ch_query_lookups
    contig_sequences = DEACON_FILTER_CONTIGS.out
    contig_read_counts = CONTIG_READ_MAPBACK.out.contig_read_counts
    filtered_bam = CONTIG_READ_MAPBACK.out.filtered_bam
    unmapped_reads = CONTIG_READ_MAPBACK.out.unmapped_reads
    unmapped_read_counts = CONTIG_READ_MAPBACK.out.unmapped_read_counts
}
