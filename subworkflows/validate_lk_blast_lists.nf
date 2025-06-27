workflow VALIDATE_LK_BLAST {
    VALIDATE_BLAST_HITS_LIST()
    VALIDATE_BLAST_FASTA_LIST()
}

/*
This workflow validates your LabKey configuration for nvd2 pipeline.

It checks:
- API key authentication and permissions
- Server connectivity and project access  
- List existence and schema compatibility
- Insert/delete permissions on target lists

The validation inserts and immediately deletes test records to verify
full read/write access without affecting real data.

Requirements:
- nvd2 secret configured with valid LabKey API key
- API key must have INSERT/DELETE permissions on target lists
- Network access to LabKey server
- Correct nextflow.config parameters
*/

process VALIDATE_BLAST_HITS_LIST {
    label 'process_low'
    secret 'nvd2'

    output:
    path 'blast_validation_report.txt'

    script:
    """
    validate_labkey.py \\
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_blast_meta_hits_list}' \
        --api_key \$nvd2 \
        --experiment_id ${params.experiment_id} \
        --type blast > blast_validation_report.txt 2>&1
    """
}

process VALIDATE_BLAST_FASTA_LIST {
    label 'process_low'
    secret 'nvd2'

    output:
    path 'fasta_validation_report.txt'

    script:
    """
    # Run validation
    validate_labkey.py \
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_blast_fasta_list}' \
        --api_key \$nvd2 \
        --experiment_id ${params.experiment_id} \
        --type fasta > fasta_validation_report.txt 2>&1
    """
}