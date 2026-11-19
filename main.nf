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

// Nextflow counts ignored task failures and then reports the run as a success
// anyway. Reading that count at the end of the run is the cheapest way to make
// every `errorStrategy 'ignore'` in the pipeline audible, including the ones
// whose failure loses a run-level output rather than a sample and so never
// shows up in the per-sample completeness check.
workflow.onComplete {
  def ignored_failures = NvdUtils.describeIgnoredFailures(
    workflow.stats.ignoredCount,
    workflow.stats.succeedCount,
    "${params.results}/report.html",
  )
  if (ignored_failures) {
    log.warn(ignored_failures)
  } else {
    log.debug "nvd.ignored_failures no task failures were ignored in this run"
  }
}
