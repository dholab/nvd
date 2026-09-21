process RESOLVE_READ_INPUTS {

    label "low"
    // cache false

    input:
    path samplesheet

    output:
    path "resolved_reads.jsonl", emit: jsonl

    script:
    """
    resolve_read_inputs.py \
        --samplesheet ${samplesheet} \
        --output-jsonl resolved_reads.jsonl
    """
}

workflow GATHER_READS {

    take:
    ch_samplesheet

    main:

        RESOLVE_READ_INPUTS(ch_samplesheet)

        ch_resolved_reads = RESOLVE_READ_INPUTS.out.jsonl
            .splitText()
            .filter { line -> line.trim() }
            .map { line -> new groovy.json.JsonSlurper().parseText(line) }

        ch_sra_accessions = ch_resolved_reads
            .filter { rec -> rec.source == "sra" }
            .map { rec ->
                def accession = rec.srr.toString().toUpperCase(java.util.Locale.ROOT)
                tuple(rec.sample_id, rec.platform, accession)
            }

        ch_local_bundles = ch_resolved_reads
             .filter { rec -> rec.source != "sra" }
             .map { rec ->
                def meta = [
                    id: rec.sample_id,
                    platform: rec.platform,
                    source: rec.source,
                    read_mode: rec.source in ["single_file", "single_glob"] ? "single" : "paired",
                ]
                meta.deacon_read_structure = meta.read_mode == "paired" ? "interleaved" : "single"
                if (meta.read_mode == "single") {
                    def reads = rec.reads.collect { file(it) }
                    meta.r1_count = reads.size()
                    return tuple(meta, reads)
                }
                def r1 = rec.r1.collect { file(it) }
                def r2 = rec.r2.collect { file(it) }
                meta.r1_count = r1.size()
                tuple(meta, r1 + r2)
            }

        emit:
        reads = ch_local_bundles
        sra_accessions = ch_sra_accessions
        resolved_reads = RESOLVE_READ_INPUTS.out.jsonl

}
