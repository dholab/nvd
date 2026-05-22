include { FETCH_FASTQ } from "../modules/sratools"

workflow GATHER_READS {

    take:
    ch_samplesheet_row

    main:

        // DOWNLOAD AND HANDLE SRA ACCESSIONS
        // ***************************************************************************/
        ch_sra_accessions = ch_samplesheet_row
            .filter { _id, srr, _platform, _fastq1, _fastq2 -> srr != null && srr != "" }
            .map { id, srr, platform, _fastq1, _fastq2 ->
                assert platform != null && platform != ""
                tuple(id, platform, srr)
            }

        ch_known_pairs = ch_samplesheet_row
            .filter { _id, srr, _platform, fastq1, fastq2 ->
                (srr == null || srr == "") && fastq1 != null && fastq1 != "" && fastq2 != null && fastq2 != ""
            }
            .map { id, _srr, platform, fastq1, fastq2 -> tuple(id, platform, fastq1, fastq2) }

        ch_known_single = ch_samplesheet_row
            .filter { _id, srr, _platform, fastq1, fastq2 ->
                (srr == null || srr == "") && fastq1 != null && fastq1 != "" && (fastq2 == null || fastq2 == "")
            }
            .map { id, _srr, platform, fastq1, _fastq2 -> tuple(id, platform, file(fastq1)) }

        FETCH_FASTQ(ch_sra_accessions)

        // Split the output FASTQs from the SRA download into two channels, where one
        // contains paired-end libraries with >=2 FASTQs, and the other contains the rest
        ch_sorted_fastqs = FETCH_FASTQ.out
            .branch { sample_id, platform, fastq_files ->

                single: fastq_files.size() == 1
                    return tuple(sample_id, platform, file(fastq_files[0]))

                paired: fastq_files.size() > 1 && file(fastq_files[0]).getName().endsWith("1.fastq") && file(fastq_files[1]).getName().endsWith("2.fastq")
                    return tuple(sample_id, platform, file(fastq_files[0]), file(fastq_files[1]))

                triple1: fastq_files.size() > 2 && file(fastq_files[0]).getName().endsWith("1.fastq") && file(fastq_files[1]).getName().endsWith("2.fastq")
                    return tuple(sample_id, platform, file(fastq_files[0]), file(fastq_files[1]))

                triple2: fastq_files.size() > 2 && file(fastq_files[1]).getName().endsWith("1.fastq") && file(fastq_files[2]).getName().endsWith("2.fastq")
                    return tuple(sample_id, platform, file(fastq_files[1]), file(fastq_files[2]))

                other: true
            }

        ch_to_be_paired = ch_known_pairs.mix(
                ch_sorted_fastqs.paired,
                ch_sorted_fastqs.triple1,
                ch_sorted_fastqs.triple2
            )

        // Emit raw reads — no interleaving. Deacon handles R1/R2 directly
        // and outputs interleaved as a byproduct of virus filtering.
        // Paired: tuple(id, platform, R1, R2)
        // Single: tuple(id, platform, fastq)
        ch_raw_reads = ch_to_be_paired.mix(
            ch_known_single,
            ch_sorted_fastqs.single
        )

        emit:
        ch_raw_reads

}
