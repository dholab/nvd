include {
    EXTRACT_HUMAN_VIRUS_READS ;
    CLASSIFY_CONTIGS_FIRST_PASS ;
    GENERATE_CONTIGS_TAXA_LIST ;
    CLASSIFY_CONTIGS_SECOND_PASS ;
    GENERATE_STAT_CONTIG_REPORT ;
    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS
} from "../modules/stat"
include { EXTRACT_HUMAN_VIRUS_CONTIGS         } from "../modules/seqkit"

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
    )

    CLASSIFY_CONTIGS_SECOND_PASS(
        GENERATE_CONTIGS_TAXA_LIST.out.combine(ch_stat_dbss).combine(ch_stat_annotation)
    )

    GENERATE_STAT_CONTIG_REPORT(
        CLASSIFY_CONTIGS_SECOND_PASS.out.combine(ch_gettax)
    )

    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS(
        CLASSIFY_CONTIGS_SECOND_PASS.out.combine(ch_gettax)
    )

    EXTRACT_HUMAN_VIRUS_CONTIGS(
        IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS.out
            .join(
                ch_filtered_reads
                .map { id, _platform, reads -> tuple(id, file(reads)) },
                by: 0
            )
    )

    emit:
    contigs = EXTRACT_HUMAN_VIRUS_CONTIGS.out
    sqlite = ch_gettax
    
}
