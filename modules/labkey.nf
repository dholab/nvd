/*
 * LabKey validation, formatting, WebDAV upload, and row-upload processes.
 */

process VALIDATE_LABKEY_BLAST_HITS_LIST {
    label 'low'
    secret 'LABKEY_API_KEY'

    input:
    val trigger

    output:
    val true, emit: validated

    script:
    """
    validate_labkey.py \
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_blast_meta_hits_list}' \
        --api_key \$LABKEY_API_KEY \
        --experiment_id ${params.experiment_id} \
        --type blast > blast_validation_report.txt 2>&1
    """
}

process VALIDATE_LABKEY_BLAST_FASTA_LIST {
    label 'low'
    secret 'LABKEY_API_KEY'

    input:
    val trigger

    output:
    val true, emit: validated

    script:
    """
    validate_labkey.py \
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_blast_fasta_list}' \
        --api_key \$LABKEY_API_KEY \
        --experiment_id ${params.experiment_id} \
        --type blast_fasta > fasta_validation_report.txt 2>&1
    """
}

process VALIDATE_LABKEY_EXPERIMENT_FRESH {
    cache false
    secret 'LABKEY_API_KEY'
    tag "exp_${params.experiment_id}"
    label "low"

    input:
    val trigger

    output:
    val true, emit: validated

    script:
    """
    labkey_check_guard_list.py \
        --mode check \
        --server ${params.labkey_server} \
        --container ${params.labkey_project_name} \
        --guard_list ${params.labkey_exp_id_guard_list} \
        --api_key \$LABKEY_API_KEY \
        --experiment_id ${params.experiment_id}
    """
}

process REGISTER_LABKEY_EXPERIMENT {
    cache false
    secret 'LABKEY_API_KEY'
    tag "exp_${params.experiment_id}"
    label 'low'

    input:
    val upload_complete

    output:
    val true, emit: registered

    script:
    """
    labkey_check_guard_list.py \
        --mode register \
        --server ${params.labkey_server} \
        --container ${params.labkey_project_name} \
        --guard_list ${params.labkey_exp_id_guard_list} \
        --api_key \$LABKEY_API_KEY \
        --experiment_id ${params.experiment_id}
    """
}

process WEBDAV_UPLOAD_BLAST {
    tag "${sample_id}"
    label 'low'
    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(blast_csv), path(fasta)
    val validation_complete

    output:
    val true, emit: done

    script:
    """
    gzip -c ${fasta} > ${fasta}.gz
    gzip -c ${blast_csv} > ${blast_csv}.gz

    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${blast_csv}.gz ${params.experiment_id}/${sample_id}/nvd/${blast_csv}.gz

    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${fasta}.gz ${params.experiment_id}/${sample_id}/nvd/${fasta}.gz
    """
}

process WEBDAV_UPLOAD_CONCATENATED {
    label 'low'
    secret 'LABKEY_API_KEY'

    input:
    path concatenated_csv
    val validation_complete

    output:
    val true, emit: done

    script:
    """
    gzip -c ${concatenated_csv} > ${concatenated_csv}.gz

    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${concatenated_csv}.gz ${params.experiment_id}/${concatenated_csv}.gz
    """
}

process PREPARE_BLAST_LABKEY {
    tag "$meta"
    label 'low'

    input:
    tuple val(meta), path(blast_csv), val(total_reads), path(contig_mapped_read_counts), val(output_name)
    val experiment_id
    val run_id
    val validation_complete

    output:
    tuple val(meta), path(output_name), emit: csv

    script:
    """
    prepare_blast_labkey.py \
        --blast-csv ${blast_csv} \
        --contig-counts ${contig_mapped_read_counts} \
        --output ${output_name} \
        --meta '${meta}' \
        --experiment-id ${experiment_id} \
        --run-id '${run_id}' \
        --total-reads ${total_reads} \
        --blast-db-version '${params.blast_db_version}'
    """
}

process CONCAT_ALL_SAMPLE_BLAST_RESULTS {
    label 'low'

    input:
    path "blast_results/*.csv"
    val experiment_id
    val validation_complete

    output:
    path "${experiment_id}_blast_concatenated.csv", emit: concatenated_csv

    script:
    """
    #!/usr/bin/env python3
    import polars as pl
    from pathlib import Path

    files = sorted(Path("blast_results").glob("*.csv"))
    pl.concat([pl.scan_csv(f) for f in files]).sink_csv("${experiment_id}_blast_concatenated.csv")
    """
}

process PREPARE_FASTA_LABKEY {
    tag "$meta"
    label 'low'

    input:
    tuple val(meta), path(fasta), val(output_name)
    val experiment_id
    val run_id
    val validation_complete

    output:
    tuple val(meta), path("${output_name}"), emit: csv

    script:
    """
    #!/usr/bin/env python3

    import csv
    from Bio import SeqIO

    fasta_data = []

    for record in SeqIO.parse('${fasta}', 'fasta'):
        labkey_row = {
            'experiment': ${experiment_id},
            'sample_id': '${meta}',
            'contig_id': record.id,
            'contig_sequence': str(record.seq),
            'notes': '',
            'snakemake_run_id': '${run_id}'
        }
        fasta_data.append(labkey_row)

    if fasta_data:
        with open('${output_name}', 'w') as f:
            fieldnames = ['experiment', 'sample_id', 'contig_id',
                         'contig_sequence', 'notes', 'snakemake_run_id']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(fasta_data)
    else:
        open('${output_name}', 'w').close()
    """
}

process LABKEY_UPLOAD_BLAST {
    tag "${sample_id}"
    label 'low'
    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(csv_file)
    val experiment_id
    val run_id

    output:
    path "blast_labkey_upload.log", emit: log

    script:
    """
    labkey_upload_blast_results.py \
        --experiment-id '${experiment_id}' \
        --run-id '${run_id}' \
        --labkey-server '${params.labkey_server}' \
        --labkey-project-name '${params.labkey_project_name}' \
        --labkey-api-key \$LABKEY_API_KEY \
        --labkey-schema '${params.labkey_schema}' \
        --table-name '${params.labkey_blast_meta_hits_list}'
    """
}

process LABKEY_UPLOAD_FASTA {
    tag "${sample_id}"
    label 'low'
    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(csv_file)
    val experiment_id
    val run_id

    output:
    path "fasta_labkey_upload.log", emit: log

    script:
    """
    labkey_upload_blast_fasta.py \
        --experiment-id '${experiment_id}' \
        --run-id '${run_id}' \
        --labkey-server '${params.labkey_server}' \
        --labkey-project-name '${params.labkey_project_name}' \
        --labkey-api-key \$LABKEY_API_KEY \
        --labkey-schema '${params.labkey_schema}' \
        --table-name '${params.labkey_blast_fasta_list}'
    """
}
