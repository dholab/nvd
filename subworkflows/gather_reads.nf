include { FETCH_FASTQ } from "../modules/sratools"

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
            .map { rec -> tuple(rec.sample_id, rec.platform, rec.srr) }

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

        FETCH_FASTQ(ch_sra_accessions)

        ch_sra_bundles = FETCH_FASTQ.out
            .map { sample_id, platform, fastq_files ->
                def files = (fastq_files instanceof List ? fastq_files : [fastq_files])
                    .collect { file(it) }
                    .sort { a, b -> a.getName() <=> b.getName() }
                def read1 = files.find { path -> path.getName().endsWith("1.fastq") || path.getName().endsWith("1.fastq.gz") || path.getName().contains("_R1_") }
                def read2 = files.find { path -> path.getName().endsWith("2.fastq") || path.getName().endsWith("2.fastq.gz") || path.getName().contains("_R2_") }
                if (read1 && read2) {
                    def meta = [id: sample_id, platform: platform, source: "sra", read_mode: "paired", r1_count: 1, deacon_read_structure: "interleaved"]
                    return tuple(meta, [read1, read2])
                }
                if (files.size() == 1) {
                    def meta = [id: sample_id, platform: platform, source: "sra", read_mode: "single", r1_count: 1, deacon_read_structure: "single"]
                    return tuple(meta, files)
                }
                throw new IllegalArgumentException("Could not determine SRA FASTQ pairing for ${sample_id}: ${files*.getName()}")
            }

        ch_raw_reads = ch_local_bundles.mix(ch_sra_bundles)

        emit:
        reads = ch_raw_reads
        resolved_reads = RESOLVE_READ_INPUTS.out.jsonl

}
