include { MAP_PAIRED_READS ; MAP_SINGLE_READS ; EXTRACT_UNMAPPED_READS } from "../modules/minimap2"
include { COUNT_MAPPED_READS } from "../modules/samtools"

workflow CONTIG_READ_MAPBACK {
    take:
    ch_screened_contigs   // tuple(sample_id, platform, read_structure, fasta, lookup) from DEACON_FILTER_CONTIGS
    ch_viral_reads        // tuple(sample_id, platform, read_structure, fastq) from PREPROCESS_READS

    main:
    ch_reads_with_contigs = ch_viral_reads
        .map { sample_id, platform, read_structure, reads -> tuple(sample_id, platform, read_structure, reads) }
        .join(
            ch_screened_contigs.map { sample_id, _platform, _read_structure, contigs, _lookup -> tuple(sample_id, contigs) },
            by: 0
        )
        .map { sample_id, platform, read_structure, reads, contigs -> tuple(sample_id, platform, read_structure, reads, contigs) }

    ch_reads_by_structure = ch_reads_with_contigs.branch { _sample_id, _platform, read_structure, _reads, _contigs ->
        paired: read_structure == "interleaved"
        single: true
    }

    MAP_PAIRED_READS(ch_reads_by_structure.paired)
    MAP_SINGLE_READS(ch_reads_by_structure.single)

    ch_mapback_bams = MAP_PAIRED_READS.out.bam.mix(MAP_SINGLE_READS.out.bam)

    COUNT_MAPPED_READS(ch_mapback_bams)

    if (params.experimental == true) {
        EXTRACT_UNMAPPED_READS(ch_reads_with_contigs)
        ch_unmapped_reads = EXTRACT_UNMAPPED_READS.out.reads
        ch_unmapped_read_counts = EXTRACT_UNMAPPED_READS.out.counts
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
