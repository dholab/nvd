process ASSESS_LONG_READ_ASSEMBLY_ELIGIBILITY {

    tag "${meta.id}"
    label "low"

    input:
    tuple val(meta), path(reads), path(profile_json)

    output:
    tuple val(meta), path(reads), path("${meta.id}.myloasm.run"), optional: true, emit: myloasm
    tuple val(meta), path(reads), path("${meta.id}.metamdbg.run"), optional: true, emit: metamdbg
    tuple val(meta), path(reads), path("${meta.id}.metaflye.run"), optional: true, emit: metaflye
    tuple val(meta.id), path("${meta.id}.long_read_assembly_eligibility.tsv"), emit: report
    tuple val(meta.id), path("${meta.id}.long_read_assembly_eligibility.json"), emit: report_summary

    script:
    """
    assess_long_read_assembly.py assess \
        --profile ${profile_json} \
        --output ${meta.id}.long_read_assembly_eligibility.tsv \
        --summary-json ${meta.id}.long_read_assembly_eligibility.json \
        --myloasm-marker ${meta.id}.myloasm.run \
        --metamdbg-marker ${meta.id}.metamdbg.run \
        --metaflye-marker ${meta.id}.metaflye.run
    """
}

process REPORT_LONG_READ_ASSEMBLY_ELIGIBILITY {

    label "low"

    input:
    path reports

    output:
    path "long_read_assembly_eligibility.tsv", emit: report

    script:
    """
    assess_long_read_assembly.py report \
        --inputs ${reports} \
        --output long_read_assembly_eligibility.tsv
    """
}

process ASSEMBLE_WITH_MYLOASM {

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), val("myloasm"), path("${sample_id}.myloasm.contigs.fasta"), emit: contigs

    script:
    """
    myloasm ${reads} -o output -t ${task.cpus}

    if [ -s output/assembly_primary.fa ]; then
        cp output/assembly_primary.fa ${sample_id}.myloasm.contigs.fasta
    else
        touch ${sample_id}.myloasm.contigs.fasta
    fi
    """
}

process ASSEMBLE_WITH_METAMDBG {

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), val("metamdbg"), path("${sample_id}.metamdbg.contigs.fasta"), emit: contigs

    script:
    // Match metaMDBG's correction overlap floor to the downstream minimum
    // contiguous contig length so short-fragment samples are not excluded by
    // the assembler's more conservative default.
    """
    metaMDBG asm \
        --out-dir output \
        --in-ont ${reads} \
        --threads ${task.cpus} \
        --min-read-overlap ${params.min_consecutive_bases}

    if [ -s output/contigs.fasta.gz ]; then
        gzip -cd output/contigs.fasta.gz > ${sample_id}.metamdbg.contigs.fasta
    elif [ -s output/contigs.fasta ]; then
        cp output/contigs.fasta ${sample_id}.metamdbg.contigs.fasta
    else
        touch ${sample_id}.metamdbg.contigs.fasta
    fi
    """
}

process ASSEMBLE_WITH_METAFLYE {

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), val("metaflye"), path("${sample_id}.metaflye.contigs.fasta"), emit: contigs

    script:
    """
    flye --nano-hq ${reads} --meta --out-dir output --threads ${task.cpus}

    if [ -s output/assembly.fasta ]; then
        cp output/assembly.fasta ${sample_id}.metaflye.contigs.fasta
    else
        touch ${sample_id}.metaflye.contigs.fasta
    fi
    """
}
