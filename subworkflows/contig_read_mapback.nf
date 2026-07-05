include { MAP_PAIRED_READS ; MAP_SINGLE_READS ; EXTRACT_UNMAPPED_READS } from "../modules/minimap2"
include { COUNT_MAPPED_READS } from "../modules/samtools"

workflow CONTIG_READ_MAPBACK {
    take:
    ch_screened_contigs   // tuple(sample_id, platform, read_structure, fasta, lookup) from DEACON_FILTER_CONTIGS
    ch_viral_reads        // tuple(sample_id, platform, read_structure, evidence_class, read_group_id, fastq) from PREPROCESS_READS

    main:
    ch_reads_with_contigs = ch_viral_reads
        .combine(
            ch_screened_contigs.map { sample_id, platform, _read_structure, contigs, _lookup -> tuple(sample_id, platform, contigs) },
            by: [0, 1]
        )
        .map { sample_id, platform, read_structure, evidence_class, read_group_id, reads, contigs -> tuple(sample_id, platform, read_structure, evidence_class, read_group_id, reads, contigs) }

    ch_reads_by_structure = ch_reads_with_contigs.branch { _sample_id, _platform, read_structure, _evidence_class, _read_group_id, _reads, _contigs ->
        paired: read_structure in ["interleaved", "merged"]
        single: true
    }

    ch_paired_reads = ch_reads_by_structure.paired
        .map { sample_id, platform, _read_structure, evidence_class, read_group_id, reads, contigs -> tuple(sample_id, platform, contigs, evidence_class, read_group_id, reads) }
        .groupTuple(by: [0, 1, 2], size: params.merge_pairs ? 2 : 1)

    ch_single_reads = ch_reads_by_structure.single
        .map { sample_id, platform, read_structure, _evidence_class, _read_group_id, reads, contigs -> tuple(sample_id, platform, read_structure, reads, contigs) }

    MAP_PAIRED_READS(ch_paired_reads)
    MAP_SINGLE_READS(ch_single_reads)

    ch_mapback_bams = MAP_PAIRED_READS.out.bam.mix(MAP_SINGLE_READS.out.bam)

    COUNT_MAPPED_READS(ch_mapback_bams)

    if (params.experimental == true) {
        EXTRACT_UNMAPPED_READS(ch_single_reads)
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
