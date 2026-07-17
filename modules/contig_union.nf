process NORMALIZE_CONTIGS {

    tag "${sample_id}"
    label "medium"

    input:
    tuple val(sample_id), val(platform), val(read_structure), val(assemblers), path(contigs)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("union.representatives.fasta"), path("union.sqlite"), path("union_meta.json"), emit: prepared

    script:
    def manifestRows = assemblers.indices
        .collect { index -> [assembler: assemblers[index], fasta: contigs[index]] }
        .sort { left, right -> left.assembler <=> right.assembler ?: left.fasta.toString() <=> right.fasta.toString() }
        .collect { row -> "printf '%s\\t%s\\n' '${row.assembler}' '${row.fasta}'" }
        .join("\n")
    """
    {
        printf '%s\t%s\n' 'assembler' 'source_fasta'
        ${manifestRows}
    } > assembler_fastas.tsv

    union_contigs.py prepare \
        --sample-id "${sample_id}" \
        --inputs-manifest assembler_fastas.tsv \
        --representatives-fasta union.representatives.fasta \
        --sqlite union.sqlite \
        --meta-json union_meta.json
    """
}

process FIND_CONTAINMENT_DUPLICATES {

    tag "${sample_id}"
    label "medium"

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(representatives), path(sqlite), path(meta_json)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("union.representatives.fasta"), path("union.sqlite"), path("containment_candidates.tsv"), emit: candidates

    script:
    """
    REPRESENTATIVE_COUNT=\$(python -c 'import json, sys; print(json.load(open(sys.argv[1]))["representative_count"])' "${meta_json}")
    MAX_CONTIG_LENGTH=\$(python -c 'import json, sys; print(json.load(open(sys.argv[1]))["max_contig_length"])' "${meta_json}")

    if [ "\$REPRESENTATIVE_COUNT" -lt 2 ]; then
        touch containment_candidates.tsv
    else
        mmseqs easy-search \
            "${representatives}" \
            "${representatives}" \
            containment_candidates.tsv \
            mmseqs_tmp \
            --search-type 3 \
            --threads ${task.cpus} \
            --strand 2 \
            --min-seq-id 1.0 \
            -c 1.0 \
            --cov-mode 2 \
            --alignment-mode 3 \
            --mask 0 \
            --max-seq-len "\$MAX_CONTIG_LENGTH" \
            --max-seqs 1000000 \
            --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,qlen,tlen,qcov,tcov,qframe,tframe"
    fi
    """
}

process EXTRACT_UNIQUE_CONTIGS {

    tag "${sample_id}"
    label "medium"

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(representatives), path(sqlite), path(candidates)

    output:
    tuple val(sample_id), val(platform), val(read_structure), val("long_read_union"), path("${sample_id}.long_read_union.contigs.fasta"), emit: contigs
    tuple val(sample_id), path("${sample_id}.long_read_union.provenance.tsv"), emit: provenance
    tuple val(sample_id), path("${sample_id}.long_read_union.summary.json"), emit: summary

    script:
    """
    union_contigs.py finalize \
        --sample-id "${sample_id}" \
        --representatives-fasta "${representatives}" \
        --candidates-tsv "${candidates}" \
        --sqlite "${sqlite}" \
        --output-fasta ${sample_id}.long_read_union.contigs.fasta \
        --provenance-tsv ${sample_id}.long_read_union.provenance.tsv \
        --summary-json ${sample_id}.long_read_union.summary.json
    """
}
