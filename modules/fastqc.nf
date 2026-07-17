process FASTQC_RAW {
    tag "${unit.alias}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(unit), path(read, stageAs: 'source_fastq')

    output:
    path "${unit.package_name}", arity: '1', emit: packages
    path "${unit.package_name}/${unit.zip_alias}", emit: zips
    path "${unit.package_name}/${unit.html_alias}", emit: htmls

    script:
    """
    mkdir -- fastqc_out "./${unit.package_name}"
    ln -s -- source_fastq "./${unit.staged_filename}"
    fastqc --quiet --outdir fastqc_out "./${unit.staged_filename}"
    cp -- "./fastqc_out/${unit.zip_alias}" "./${unit.package_name}/${unit.zip_alias}"
    cp -- "./fastqc_out/${unit.html_alias}" "./${unit.package_name}/${unit.html_alias}"
    write_nvd_fastqc_receipt.py \
        --output='./${unit.package_name}/unit.json' \
        --sample-id='${unit.sample_id}' \
        --platform='${unit.platform}' \
        --source='${unit.source}' \
        --read-end='${unit.read_end}' \
        --input-ordinal='${unit.input_ordinal}' \
        --alias='${unit.alias}' \
        --zip-alias='${unit.zip_alias}' \
        --html-alias='${unit.html_alias}'
    """
}
