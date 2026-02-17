#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATHER_READS         } from "./workflows/gather_reads"
include { PREPROCESS_READS     } from "./workflows/preprocess_reads"
include { STAT_BLAST_WORKFLOW  } from "./workflows/stat_blast_workflow"
include { GOTTCHA2_WORKFLOW    } from "./workflows/gottcha2_workflow"
include { CLUMPIFY_WORKFLOW    } from "./workflows/clumpify"

workflow {

    assert params.samplesheet && file(params.samplesheet).isFile()

    ch_input_samplesheet = channel.fromPath(params.samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row -> tuple(row.sample_id, row.srr, row.platform, row.fastq1, row.fastq2) }
        .filter { it -> !it[0].startsWith("#") }
        .unique()

    GATHER_READS(ch_input_samplesheet)

    PREPROCESS_READS(GATHER_READS.out)

    // STAT+BLAST workflow (human virus detection via STAT + two-phase BLAST)
    // Aliases: nvd, stat, blast, stat_blast, stast (for backward compatibility)
    def stat_blast_results = STAT_BLAST_WORKFLOW(PREPROCESS_READS.out)
    def stat_blast_token = stat_blast_results.completion

    // Collect all LabKey logs as final process
    if (params.labkey) {
        stat_blast_results.labkey_log.collectFile(
            name: 'final_labkey_upload.log',
            storeDir: params.results,
        )
    }

    // GOTTCHA2 workflow
    def gottcha_token = GOTTCHA2_WORKFLOW(PREPROCESS_READS.out).completion
    def start_clumpify = stat_blast_token.mix(gottcha_token).collect().map { true }

    // CLUMPIFY workflow  
    CLUMPIFY_WORKFLOW(
        GATHER_READS.out,
        start_clumpify,
    )
}
