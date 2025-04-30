#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATHER_READS } from "./workflows/gather_reads"
include { NVD2_WORKFLOW } from "./workflows/nvd2_workflow"

workflow {

/*

EXAMPLE COMMAND
----------------------------------

nextflow run . \
--tools nvd \
--blast_db ./db \
--blast_db_prefix ./db/PP819512-nt \
--stat_index ./db/tree_index.dense.dbs \
--stat_dbss ./db/tree_index.dense.dbss \
--stat_annotation ./db/tree_index.dense.dbss.annotation \
--human_virus_taxlist ./db/human_viruses_taxlist.txt \
--samplesheet test.csv \
-resume

*/

    assert params.samplesheet && (params.samplesheet).isFile()

    ch_input_samplesheet = Channel.fromPath( params.samplesheet )
        .splitCsv(header: true, strip: true)
        .filter { it -> !it[0].startsWith("#") }
        .map { row -> tuple(row.sample_id), row.srr, row.platform, file(row.fastq1), file(row.fastq2)}
        .unique()

    GATHER_READS(ch_input_samplesheet)

    if (params.all || params.nvd || params.stat || (params.tool && params.tool.contains("nvd"))) {
        NVD2_WORKFLOW (GATHER_READS.out)
    }

    // if (params.all || params.gottcha2 || (params.tool && params.tool.contains("gottcha"))) {
    //     GOTTCHA2_WORKFLOW(GATHER_READS.out)
    // }
    //
    // if (params.all || params.clumpify || (params.tool && params.tool.contains("clump"))
    //     CLUMPIFY_WORKFLOW(GATHER_READS.out)
    // }

}
