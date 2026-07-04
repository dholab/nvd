include { MAP_READS_TO_CONTIGS } from "../modules/minimap2"
include { COUNT_MAPPED_READS   } from "../modules/samtools"

workflow MAP_READS_TO_SCREENED_CONTIGS {
    take:
    ch_screened_contigs   // tuple(sample_id, platform, read_structure, fasta, lookup) from SCREEN_CONTIGS
    ch_viral_reads        // tuple(sample_id, platform, read_structure, fastq) from PREPROCESS_READS

    main:
    ch_reads_with_contigs = ch_viral_reads
        .map { sample_id, platform, _read_structure, reads -> tuple(sample_id, platform, reads) }
        .join(
            ch_screened_contigs.map { sample_id, _platform, _read_structure, contigs, _lookup -> tuple(sample_id, contigs) },
            by: 0
        )
        .map { sample_id, platform, reads, contigs -> tuple(sample_id, platform, reads, contigs) }

    MAP_READS_TO_CONTIGS(ch_reads_with_contigs)

    COUNT_MAPPED_READS(MAP_READS_TO_CONTIGS.out)

    emit:
    contig_read_counts = COUNT_MAPPED_READS.out.mapped_counts
    filtered_bam = COUNT_MAPPED_READS.out.filtered_bam
}
