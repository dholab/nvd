process ASSEMBLE_WITH_MYLOASM {

    tag "${sample_id}"
    label "high"

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), val("myloasm"), path("${sample_id}.myloasm.contigs.fasta"), emit: contigs

    script:
    """
    set +e
    myloasm ${reads} -o output -t ${task.cpus} > ${sample_id}.myloasm.log 2>&1
    status=\$?
    set -e

    if [ "\$status" -eq 0 ] && [ -s output/assembly_primary.fa ]; then
        cp output/assembly_primary.fa ${sample_id}.myloasm.contigs.fasta
    else
        touch ${sample_id}.myloasm.contigs.fasta
    fi
    """
}

process ASSEMBLE_WITH_METAMDBG {

    tag "${sample_id}"
    label "high"

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), val("metamdbg"), path("${sample_id}.metamdbg.contigs.fasta"), emit: contigs

    script:
    """
    set +e
    metaMDBG asm --out-dir output --in-ont ${reads} --threads ${task.cpus} > ${sample_id}.metamdbg.log 2>&1
    status=\$?
    set -e

    if [ "\$status" -eq 0 ] && [ -s output/contigs.fasta.gz ]; then
        gzip -cd output/contigs.fasta.gz > ${sample_id}.metamdbg.contigs.fasta
    elif [ "\$status" -eq 0 ] && [ -s output/contigs.fasta ]; then
        cp output/contigs.fasta ${sample_id}.metamdbg.contigs.fasta
    else
        touch ${sample_id}.metamdbg.contigs.fasta
    fi
    """
}

process ASSEMBLE_WITH_METAFLYE {

    tag "${sample_id}"
    label "high"

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), val("metaflye"), path("${sample_id}.metaflye.contigs.fasta"), emit: contigs

    script:
    """
    set +e
    flye --nano-hq ${reads} --meta --out-dir output --threads ${task.cpus} > ${sample_id}.metaflye.log 2>&1
    status=\$?
    set -e

    if [ "\$status" -eq 0 ] && [ -s output/assembly.fasta ]; then
        cp output/assembly.fasta ${sample_id}.metaflye.contigs.fasta
    else
        touch ${sample_id}.metaflye.contigs.fasta
    fi
    """
}
