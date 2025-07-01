workflow BUNDLE_GOTTCHA2_FOR_LABKEY {
    take:
    gottcha2_results     // channel: [ meta, csv ]
    gottcha2_extracted_fastas  // channel: [ meta, fasta ]

    main:

    //gottcha2_results.view()
    //ottcha2_extracted_fastas.view()

    LABKEY_UPLOAD_GOTTCHA2_FULL(gottcha2_results)

    LABKEY_UPLOAD_GOTTCHA2_FASTA(gottcha2_extracted_fastas)

    // WIP this will upload to the final LK list for the verified via blast contigs
    //LABKEY_UPLOAD_GOTTCHA2_HITS_VERIFIED()

}

process LABKEY_UPLOAD_GOTTCHA2_FULL {
    label 'process_low'

    tag "${sample_id}"

    secret 'nvd2'

    input:
    tuple val(sample_id), path(full_tsv), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    path "fasta_labkey_upload.log", emit: log

    script:
    """
    labkey_upload_gottcha2_full.py \
    --input-tsv ${full_tsv} \
    --sample ${sample_id} \
    --experiment '${params.experiment_id}' \
    --labkey-server '${params.labkey_server}' \
    --labkey-project-name '${params.labkey_project_name}' \
    --labkey-api-key \$nvd2 \
    --db-version ${params.gottcha2_db_version} \
    --labkey-schema '${params.labkey_schema}' \
    --table-name '${params.labkey_gottcha_full_list}' \
    """
}

//[illumina_test.extract.fasta, illumina_test.full.tsv, illumina_test.gottcha_strain.log, illumina_test.lineage.tsv, illumina_test.tsv]

process LABKEY_UPLOAD_GOTTCHA2_FASTA {

    tag "${sample_id}"

    secret 'nvd2'

    input:
    tuple val(sample_id), path(fasta), path(full_tsv), path(strain_log), path(lineage_tsv)

    output:
    path("${sample_id}_df.tsv")

    script:
    """
    awk 'BEGIN {
        FS = "\\n"; RS = ">";
        OFS = "\\t";
        print "Experiment", "Sample Id", "Header", "Sequence", "Notes", "Snakemake Run Id", "Upload Type"
    }
    NF > 1 {
        header = \$1;
        sub(/^>/, "", header);
        gsub(/\\n/, "", \$2);  # remove newlines in sequence
        sequence = \$2;
        if (length(sequence) > 0) {
            print ${params.experiment_id}, "${sample_id}", header, sequence, "gottcha2", "GOTTCHA2_UPLOAD", "fasta";
        }
    }' ${fasta} > ${sample_id}_df.tsv

    labkey_upload_gottcha2_fasta.py \
        --server ${params.labkey_server} \
        --container ${params.labkey_project_name} \
        --list ${params.labkey_gottcha_fasta_list} \
        --api_key \$nvd2 \
        --input_tsv ${sample_id}_df.tsv

    """
}


// process GENERATE_FASTA {

//     tag "${sample_id}"
//     // publishDir params.gottcha_fasta, mode: 'copy', overwrite: false, pattern: "*.fasta"

//     errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
// 	maxRetries 2

//     input:
//     tuple val(sample_id), path(sam), path(ref_mmi), path(stats), path(tax_tsv)

//     output:
//     path "*"

//     when:
//     (params.tools && params.tools.contains("gottcha2") || params.tools.contains("gottcha")) || params.all || params.gottcha

//     script:
//     def ref_prefix = file(ref_mmi).getBaseName().toString().replace(".mmi", "")
//     """
//     gottcha2.py \
//     --noCutoff --dbLevel strain --threads ${task.cpus} -ef \
//     --database ${ref_prefix} \
//     --prefix ${sample_id} \
//     --sam ${sam}
//     """
// }
