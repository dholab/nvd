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

    main:
    CLASSIFY_CONTIGS_FIRST_PASS(
        ch_filtered_reads.combine(ch_stat_dbs)
    )

    GENERATE_CONTIGS_TAXA_LIST(
        CLASSIFY_CONTIGS_FIRST_PASS.out
    )

    CLASSIFY_CONTIGS_SECOND_PASS(
        GENERATE_CONTIGS_TAXA_LIST.report.out.combine(ch_stat_dbss)
    )

    GENERATE_STAT_CONTIG_REPORT(
        CLASSIFY_CONTIGS_SECOND_PASS.out
    )

    ch_gettax_sqlite_path = GENERATE_STAT_CONTIG_REPORT.out.gettax_sqlite

    IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS(
        CLASSIFY_CONTIGS_SECOND_PASS.out.combine(ch_gettax_sqlite_path)
    )

    EXTRACT_HUMAN_VIRUS_CONTIGS(
        IDENTIFY_HUMAN_VIRUS_FAMILY_CONTIGS.out.mix(ch_filtered_reads).groupTuple(by: 0)
    )

    emit:
    contigs = EXTRACT_HUMAN_VIRUS_CONTIGS.out
    sqlite = ch_gettax_sqlite_path
    
}
