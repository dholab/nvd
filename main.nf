#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATHER_READS         } from "./workflows/gather_reads"
include { STAT_BLAST_WORKFLOW  } from "./workflows/stat_blast_workflow"

workflow {

    assert params.samplesheet && file(params.samplesheet).isFile()

    ch_input_samplesheet = channel.fromPath(params.samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row -> tuple(row.sample_id, row.srr, row.platform, row.fastq1, row.fastq2) }
        .filter { it -> !it[0].startsWith("#") }
        .unique()

    GATHER_READS(ch_input_samplesheet)

    def stat_blast_results = STAT_BLAST_WORKFLOW(GATHER_READS.out)

    if (params.labkey) {
        stat_blast_results.labkey_log.collectFile(
            name: 'final_labkey_upload.log',
            storeDir: params.results,
        )
    }
}
