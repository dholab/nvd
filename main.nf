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

    // STAT+BLAST workflow: receives pre-interleave R1/R2 so deacon can filter
    // and interleave in one step, skipping the ~1hr INTERLEAVE_PAIRS bottleneck.
    // Deacon virus extraction + preprocessing happen inside the workflow.
    def stat_blast_results = STAT_BLAST_WORKFLOW(GATHER_READS.out.ch_pre_interleave)
    def stat_blast_token = stat_blast_results.completion

    // Collect all LabKey logs as final process
    if (params.labkey) {
        stat_blast_results.labkey_log.collectFile(
            name: 'final_labkey_upload.log',
            storeDir: params.results,
        )
    }

    // GOTTCHA2 workflow: full interleaved reads, preprocessed — waits for STAT_BLAST
    // to finish so cluster resources from the fast STAT_BLAST path are freed up first.
    def gottcha_start = stat_blast_token.map { true }
    PREPROCESS_READS(
        GATHER_READS.out.ch_gathered_reads.combine(gottcha_start).map { it[0..-2] }
    )
    def gottcha_token = GOTTCHA2_WORKFLOW(PREPROCESS_READS.out).completion

    // CLUMPIFY workflow: raw interleaved reads, waits for both workflows
    def start_clumpify = stat_blast_token.mix(gottcha_token).collect().map { true }
    CLUMPIFY_WORKFLOW(
        GATHER_READS.out.ch_gathered_reads,
        start_clumpify,
    )
}
