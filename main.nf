#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NVD_MAIN } from "./workflows/nvd_main"
include { RECORD_RUN_METADATA } from "./modules/utils"

workflow {

  assert params.samplesheet && file(params.samplesheet).isFile()

  def source_revision = System.getenv('NVD_SOURCE_REVISION') ?: workflow.commitId ?: 'unknown'
  def source_dirty = System.getenv('NVD_SOURCE_DIRTY') ?: 'unknown'
  RECORD_RUN_METADATA(
    workflow.commandLine,
    workflow.manifest.version ?: 'unknown',
    source_revision,
    source_dirty,
  )

  ch_input_samplesheet = channel.fromPath(params.samplesheet)

  nvd_main_results = NVD_MAIN(ch_input_samplesheet, RECORD_RUN_METADATA.out.version_file)
}
