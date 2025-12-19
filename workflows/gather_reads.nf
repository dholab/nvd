include { FETCH_FASTQ } from "../modules/sratools"
include { INTERLEAVE_PAIRS ; MERGE_PAIRS } from "../modules/bbmap"

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
            .map { id, _srr, platform, fastq1, fastq2 -> tuple( id, platform, fastq1, fastq2 )}

        ch_known_single = ch_samplesheet_row
            .filter { _id, srr, _platform, fastq1, fastq2 ->
                (srr == null || srr == "") && fastq1 != null && fastq1 != "" && (fastq2 == null || fastq2 == "")
            }
            .map { id, _srr, platform, fastq1, _fastq2 -> tuple( id, platform, file(fastq1))}

        FETCH_FASTQ(ch_sra_accessions)

        // split the output FASTQs from the SRA download into two channels, where one
        // contains paired-end libraries with >=2 FASTQs, and the other contains the rest
        ch_sorted_fastqs = FETCH_FASTQ.out
            .branch { sample_id, platform, fastq_files ->

                // If there's just one FASTQ, "unpack" it from the array returned by the glob
                single: fastq_files.size() == 1
                    return tuple(sample_id, platform, file(fastq_files[0]))

                // If there are two FASTQs, expect that the alphanumeric first will end with ".1.fastq" and the second with ".2.fastq",
                // which is a (mostly) reliable SRA convention
                paired: fastq_files.size() > 1 && file(fastq_files[0]).getName().endsWith("1.fastq") && file(fastq_files[1]).getName().endsWith("2.fastq")
                    return tuple(sample_id, platform, file(fastq_files[0]), file(fastq_files[1]))

                // There are a couple common cases of >2 FASTQs per accession that we can handle. The first is where the first two files
                // end with "1.fastq" and ".2.fastq" and the third ends with "3.fastq". Assuming correct alphanumeric sorting, we handle
                // that in this branch.
                triple1: fastq_files.size() > 2 && file(fastq_files[0]).getName().endsWith("1.fastq") && file(fastq_files[1]).getName().endsWith("2.fastq")
                    return tuple(sample_id, platform, file(fastq_files[0]), file(fastq_files[1]))

                // It's also possible that the third, non-R1/R2 reads are in a FASTQ that doesn't have a numbered suffix, e.g.,
                // SRR33146255.fastq. When that's the case, that third FASTQ will end up first when the files are sorted by name.
                // We handle that case in this branch by indexing out the second and third FASTQ (groovy/nextflow are 0-indexed)
                triple2: fastq_files.size() > 2 && file(fastq_files[1]).getName().endsWith("1.fastq") && file(fastq_files[2]).getName().endsWith("2.fastq")
                    return tuple(sample_id, platform, file(fastq_files[1]), file(fastq_files[2]))

                // Other cases are as-yet unsupported
                other: true
            }


        // SEPARATE OUT SETS OF FASTQS THAT WILL NEED TO BE MERGED OR INTERLEAVED
        // ***************************************************************************/
        ch_to_be_paired = ch_known_pairs.mix(
                ch_sorted_fastqs.paired,
                ch_sorted_fastqs.triple1,
                ch_sorted_fastqs.triple2
            )


        // MERGE OR INTERLEAVE ANY PAIRED READS
        // ***************************************************************************/
        ch_combined_fastqs = params.merge_pairs
            ? MERGE_PAIRS(ch_to_be_paired)
            : INTERLEAVE_PAIRS(ch_to_be_paired)


        // MIX READS FROM DIFFERENT PLATFORMS AND EMIT
        // ***************************************************************************/
        ch_gathered_reads = ch_known_single.mix(
            ch_combined_fastqs,
            ch_sorted_fastqs.single
        )

        emit:
        ch_gathered_reads

}
