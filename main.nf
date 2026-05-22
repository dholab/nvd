#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NVD_MAIN } from "./workflows/nvd_main"

workflow {

  assert params.samplesheet && file(params.samplesheet).isFile()

  ch_input_samplesheet = channel.fromPath(params.samplesheet)
    .splitCsv(header: true, strip: true)
    .map { row -> tuple(row.sample_id, row.srr, row.platform, row.fastq1, row.fastq2) }
    .filter { it -> !it[0].startsWith("#") }
    .unique()

  nvd_main_results = NVD_MAIN(ch_input_samplesheet)
}
