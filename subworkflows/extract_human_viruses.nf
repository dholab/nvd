include {
    EXTRACT_HUMAN_VIRUS_READS ;
    CLASSIFY_CONTIGS_FIRST_PASS ;
    GENERATE_CONTIGS_TAXA_LIST ;
    CLASSIFY_CONTIGS_SECOND_PASS ;
    GENERATE_STAT_CONTIG_REPORT ;
    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS
} from "../modules/stat"
include { EXTRACT_HUMAN_VIRUS_CONTIGS         } from "../modules/seqkit"
include { MAP_READS_TO_CONTIGS                } from "../modules/minimap2"
include { COUNT_MAPPED_READS                  } from "../modules/samtools"

workflow EXTRACT_HUMAN_VIRUSES {
    take:
    ch_filtered_reads
    ch_stat_dbs
    ch_stat_dbss
    ch_stat_annotation
    ch_gettax

    main:
    CLASSIFY_CONTIGS_FIRST_PASS(
        ch_filtered_reads.combine(ch_stat_dbs)
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
            }
            .combine(ch_gettax)
    )

    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS(
        CLASSIFY_CONTIGS_SECOND_PASS.out
            .filter { _id, hits_file ->
                hits_file.size() > 0 && hits_file.readLines().size() > 0
            }
            .combine(ch_gettax)
    )

    EXTRACT_HUMAN_VIRUS_CONTIGS(
        IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS.out
            .filter { _id, hits_file ->
                hits_file.size() > 0 && hits_file.readLines().size() > 0
            }
            .join(
                ch_filtered_reads
                .map { id, _platform, reads -> tuple(id, file(reads)) },
                by: 0
            )
    )

    MAP_READS_TO_CONTIGS(
        ch_filtered_reads,
        EXTRACT_HUMAN_VIRUS_CONTIGS.out
    )

    COUNT_MAPPED_READS(MAP_READS_TO_CONTIGS.out)

    emit:
    contigs = EXTRACT_HUMAN_VIRUS_CONTIGS.out
    sqlite = ch_gettax
    
}
