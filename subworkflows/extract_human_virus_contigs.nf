include {
    CLASSIFY_CONTIGS_FIRST_PASS ;
    GENERATE_CONTIGS_TAXA_LIST ;
    CLASSIFY_CONTIGS_SECOND_PASS ;
    GENERATE_STAT_CONTIG_REPORT ;
    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS
} from "../modules/stat"
include { EXTRACT_HUMAN_VIRUS_CONTIGS         } from "../modules/seqkit"
include { MAP_READS_TO_CONTIGS                } from "../modules/minimap2"
include { MARK_DUPLICATES ; COUNT_MAPPED_READS } from "../modules/samtools"

workflow EXTRACT_HUMAN_VIRUSES {
    take:
    ch_filtered_contigs  // tuple(sample_id, platform, read_structure, fasta) from PREPROCESS_CONTIGS
    ch_viral_reads       // tuple(sample_id, platform, read_structure, fastq) from EXTRACT_HUMAN_VIRUS_READS
    ch_stat_dbs
    ch_stat_dbss
    ch_stat_annotation
    ch_state_dir  // value channel: state directory path for taxonomy lookups

    main:
    CLASSIFY_CONTIGS_FIRST_PASS(
        ch_filtered_contigs.combine(ch_stat_dbs)
    )

    GENERATE_CONTIGS_TAXA_LIST(
        CLASSIFY_CONTIGS_FIRST_PASS.out
            .filter { _id, _contigs, hits_file ->
                hits_file.size() > 0 && hits_file.readLines().size() > 0
            }
    )

    CLASSIFY_CONTIGS_SECOND_PASS(
        GENERATE_CONTIGS_TAXA_LIST.out
            .combine(ch_stat_dbss)
            .combine(ch_stat_annotation)
    )

    GENERATE_STAT_CONTIG_REPORT(
        CLASSIFY_CONTIGS_SECOND_PASS.out
            .filter { _id, hits_file ->
                hits_file.size() > 0 && hits_file.readLines().size() > 0
            },
        ch_state_dir
    )

    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS(
        CLASSIFY_CONTIGS_SECOND_PASS.out
            .filter { _id, hits_file ->
                hits_file.size() > 0 && hits_file.readLines().size() > 0
            },
        ch_state_dir
    )

    EXTRACT_HUMAN_VIRUS_CONTIGS(
        IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS.out
            .filter { _id, hits_file ->
                hits_file.size() > 0 && hits_file.readLines().size() > 0
            }
            .join(
                ch_filtered_contigs
                .map { id, _platform, _read_structure, contigs -> tuple(id, file(contigs)) },
                by: 0
            )
    )

    // Pass through the extracted virus reads from STAT and align them to the assembled and QC-checked
    // SPAdes contigs that have been identified as human-infecting virus family members.
    // Join viral reads with extracted contigs by sample_id
    ch_reads_with_contigs = ch_viral_reads
        .map { sample_id, platform, _read_structure, reads -> tuple(sample_id, platform, reads) }
        .join(
            EXTRACT_HUMAN_VIRUS_CONTIGS.out.map { sample_id, contigs -> tuple(sample_id, contigs) },
            by: 0
        )
        // Result: [sample_id, platform, reads, contigs]

    MAP_READS_TO_CONTIGS(ch_reads_with_contigs)

    MARK_DUPLICATES(MAP_READS_TO_CONTIGS.out)

    COUNT_MAPPED_READS(MARK_DUPLICATES.out)

    emit:
    contigs = EXTRACT_HUMAN_VIRUS_CONTIGS.out
    contig_read_counts = COUNT_MAPPED_READS.out.mapped_counts

}
