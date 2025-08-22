#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATHER_READS } from "./workflows/gather_reads"
include { NVD2_WORKFLOW } from "./workflows/nvd2_workflow"
include { GOTTCHA2_WORKFLOW } from "./workflows/gottcha2_workflow"
include { CLUMPIFY_WORKFLOW } from "./workflows/clumpify"

workflow {

/*
EXAMPLE COMMAND
----------------------------------
nextflow run . \
--tools nvd \              # Run only NVD
--tools gottcha \          # Run only GOTTCHA2
--tools clumpify \         # Run only CLUMPIFY
--tools all \              # Run all tools
--tools nvd,gottcha \      # Run multiple specific tools
--blast_db ./db \
--blast_db_prefix PP819512-nt \
--stat_index ./db/tree_index.dense.dbs \
--stat_dbss ./db/tree_index.dense.dbss \
--stat_annotation ./db/tree_index.dense.dbss.annotation \
--human_virus_taxlist ./db/human_viruses_taxlist.txt \
--samplesheet test.csv \
-resume
*/

    assert params.samplesheet && file(params.samplesheet).isFile()

    ch_input_samplesheet = Channel.fromPath( params.samplesheet )
        .splitCsv(header: true, strip: true)
        .map { row -> tuple(row.sample_id), row.srr, row.platform, row.fastq1, row.fastq2}
        .filter { it -> !it[0].startsWith("#") }
        .unique()

    GATHER_READS(ch_input_samplesheet)

    // NVD workflow
    def nvd_results = NVD2_WORKFLOW(GATHER_READS.out)
    def nvd_token = nvd_results.completion

    // Collect all LabKey logs as final process
    if (params.labkey) {
        nvd_results.labkey_log.collectFile(
            name: 'final_labkey_upload.log',
            storeDir: params.results)
    }

    // GOTTCHA2 workflow
    def gottcha_token = GOTTCHA2_WORKFLOW(GATHER_READS.out).completion

    def start_clumpify = nvd_token.mix(gottcha_token).collect().map { true }

    // CLUMPIFY workflow  
    CLUMPIFY_WORKFLOW(
        GATHER_READS.out,
        start_clumpify
    )

}
