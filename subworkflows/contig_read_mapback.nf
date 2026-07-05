include { MAP_PAIRED_READS ; MAP_SINGLE_READS ; EXTRACT_UNMAPPED_READS } from "../modules/minimap2"
include { COUNT_MAPPED_READS } from "../modules/samtools"

workflow CONTIG_READ_MAPBACK {
    take:
    ch_screened_contigs   // tuple(sample_id, platform, read_structure, fasta, lookup) from DEACON_FILTER_CONTIGS
    ch_paired_reads       // tuple(sample_id, platform, paired_read_files) from PREPROCESS_READS
    ch_single_reads       // tuple(sample_id, platform, read_structure, fastq) from PREPROCESS_READS

    main:
    ch_contigs = ch_screened_contigs
        .map { sample_id, platform, _read_structure, contigs, _lookup -> tuple(sample_id, platform, contigs) }

    ch_paired_reads_with_contigs = ch_paired_reads
        .combine(
            ch_contigs,
            by: [0, 1]
        )
        .map { sample_id, platform, paired_read_files, contigs -> tuple(sample_id, platform, contigs, paired_read_files) }

    ch_single_reads_with_contigs = ch_single_reads
        .combine(
            ch_contigs,
            by: [0, 1]
        )
        .map { sample_id, platform, read_structure, reads, contigs -> tuple(sample_id, platform, read_structure, reads, contigs) }

    MAP_PAIRED_READS(ch_paired_reads_with_contigs)
    MAP_SINGLE_READS(ch_single_reads_with_contigs)

    ch_mapback_bams = MAP_PAIRED_READS.out.bam.mix(MAP_SINGLE_READS.out.bam)

    COUNT_MAPPED_READS(ch_mapback_bams)

    if (params.experimental == true) {
        EXTRACT_UNMAPPED_READS(ch_single_reads_with_contigs)
        ch_unmapped_reads = MAP_PAIRED_READS.out.overlap_unmapped_reads
            .mix(MAP_PAIRED_READS.out.single_unmapped_reads)
            .mix(EXTRACT_UNMAPPED_READS.out.reads)
        ch_unmapped_read_counts = MAP_PAIRED_READS.out.overlap_unmapped_counts
            .mix(MAP_PAIRED_READS.out.single_unmapped_counts)
            .mix(EXTRACT_UNMAPPED_READS.out.counts)
    } else {
        ch_unmapped_reads = channel.empty()
        ch_unmapped_read_counts = channel.empty()
    }

    emit:
    contig_read_counts = COUNT_MAPPED_READS.out.mapped_counts
    filtered_bam = COUNT_MAPPED_READS.out.filtered_bam
    unmapped_reads = ch_unmapped_reads
    unmapped_read_counts = ch_unmapped_read_counts
}
