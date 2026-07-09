include { SOURMASH_SKETCH_QUERY_METAGENOME ; SOURMASH_FETCH_REF_SKETCH ; SOURMASH_SKETCH_REF_FASTA ; SOURMASH_GATHER_QUERY_METAGENOME ; SOURMASH_FETCH_LINEAGES ; SOURMASH_STAGE_REFERENCE ; SOURMASH_TAX_METAGENOME } from "../modules/sourmash"
include { CONCAT_READS_AS_FASTA } from "../modules/seqkit"

workflow RAPID_SCREENING {
    take:
    ch_preprocessed_read_shards  // tuple(meta, reads, profile_json, length_histogram); meta includes id, platform, read_structure, query_class, sequence_count

    main:
    def is_http_url = { value -> value && value ==~ /(?i)^https?:\/\/.+/ }
    assert !is_http_url(params.sourmash_ref_path) : "sourmash_ref_path must be a local path; use sourmash_ref_url for URLs."
    assert !is_http_url(params.sourmash_ref_fasta) : "sourmash_ref_fasta must be a local path."
    assert !is_http_url(params.sourmash_lineages_path) : "sourmash_lineages_path must be a local path; use sourmash_lineages_url for URLs."
    assert !params.sourmash_ref_url || is_http_url(params.sourmash_ref_url) : "sourmash_ref_url must be an HTTP(S) URL."
    assert !params.sourmash_lineages_url || is_http_url(params.sourmash_lineages_url) : "sourmash_lineages_url must be an HTTP(S) URL."

    // Priority: explicit local sketch → URL download → sketch local FASTA.
    // Each condition guards against the higher-priority source being set, so at
    // most one channel is non-empty and the mix passes through exactly one ref.
    ch_local_ref_sketch = params.sourmash_ref_path
        ? channel.fromPath(params.sourmash_ref_path)
        : channel.empty()
    ch_ref_sketch_url = (!params.sourmash_ref_path && params.sourmash_ref_url)
        ? channel.of(params.sourmash_ref_url)
        : channel.empty()
    ch_ref_fasta = (!params.sourmash_ref_path && !params.sourmash_ref_url && params.sourmash_ref_fasta)
        ? channel.fromPath(params.sourmash_ref_fasta)
        : channel.empty()

    SOURMASH_FETCH_REF_SKETCH(ch_ref_sketch_url)
    SOURMASH_SKETCH_REF_FASTA(ch_ref_fasta)

    ch_ref_sketch = ch_local_ref_sketch
        .mix(SOURMASH_FETCH_REF_SKETCH.out.ref_sketch)
        .mix(SOURMASH_SKETCH_REF_FASTA.out.ref_sketch)

    ch_local_lineages = params.sourmash_lineages_path
        ? channel.fromPath(params.sourmash_lineages_path)
        : channel.empty()
    ch_lineages_url = (!params.sourmash_lineages_path && params.sourmash_lineages_url)
        ? channel.of(params.sourmash_lineages_url)
        : channel.empty()

    SOURMASH_FETCH_LINEAGES(ch_lineages_url)

    ch_lineages = ch_local_lineages.mix(SOURMASH_FETCH_LINEAGES.out.lineages)
    SOURMASH_STAGE_REFERENCE(ch_ref_sketch.combine(ch_lineages))

    // Avoid scheduling sourmash for read files that cannot produce a query
    // sketch. Use cached profiler counts rather than scanning FASTQ/FASTA here.
    ch_grouped_read_shards = ch_preprocessed_read_shards
        .map { meta, reads, _profile_json, _length_histogram ->
            tuple(meta.id, meta.platform, meta.read_structure, [meta: meta, reads: reads])
        }
        .groupTuple(by: [0, 1, 2])
        .map { sample_id, platform, read_structure, rows ->
            def count = rows.collect { row -> row.meta.sequence_count }.sum()
            def reads = rows.collect { row -> row.reads }
            tuple([id: sample_id, platform: platform, read_structure: read_structure, sequence_count: count], reads)
        }
        .filter { meta, _reads -> meta.sequence_count > 0 }
        .map { meta, reads -> tuple(meta.id, meta.platform, meta.read_structure, reads) }

    CONCAT_READS_AS_FASTA(ch_grouped_read_shards)
    SOURMASH_SKETCH_QUERY_METAGENOME(CONCAT_READS_AS_FASTA.out)

    ch_gather_inputs = SOURMASH_SKETCH_QUERY_METAGENOME.out.query_sketches
        .combine(SOURMASH_STAGE_REFERENCE.out.ref_sketch)

    SOURMASH_GATHER_QUERY_METAGENOME(ch_gather_inputs)

    ch_taxable_gather_csv = SOURMASH_GATHER_QUERY_METAGENOME.out.gather_csv.filter { _sample_id, _platform, _read_structure, gather_csv ->
        file(gather_csv).size() > 0
    }

    ch_tax_inputs = ch_taxable_gather_csv
        .combine(SOURMASH_STAGE_REFERENCE.out.lineages)

    SOURMASH_TAX_METAGENOME(ch_tax_inputs)

    emit:
    query_sketches = SOURMASH_SKETCH_QUERY_METAGENOME.out.query_sketches
    ref_sketch     = SOURMASH_STAGE_REFERENCE.out.ref_sketch
    lineages       = SOURMASH_STAGE_REFERENCE.out.lineages
    gather_csv     = SOURMASH_GATHER_QUERY_METAGENOME.out.gather_csv
    tax_reports    = SOURMASH_TAX_METAGENOME.out.tax_reports
}
