include { FASTQC_RAW } from "../modules/fastqc"
include { BUILD_NVD_MULTIQC_INPUTS; MULTIQC_NVD } from "../modules/multiqc"

workflow NVD_MULTIQC_REPORT {
    take:
    ch_reads
    ch_resolved_reads
    ch_nvd_version
    ch_experimental_enabled
    ch_multiqc_config

    main:
    ch_fastqc_units = ch_reads.flatMap { meta, reads -> NvdReporting.processReadyFastqcTuples(meta, reads) }

    FASTQC_RAW(ch_fastqc_units)

    ch_fastqc_packages = FASTQC_RAW.out.packages.toList()
    ch_fastqc_zips = FASTQC_RAW.out.zips.toList()

    BUILD_NVD_MULTIQC_INPUTS(
        ch_resolved_reads,
        ch_nvd_version,
        ch_fastqc_packages,
        ch_experimental_enabled,
    )

    MULTIQC_NVD(
        ch_fastqc_zips,
        BUILD_NVD_MULTIQC_INPUTS.out.inputs,
        ch_multiqc_config,
    )

    emit:
    report = MULTIQC_NVD.out.report
    data = MULTIQC_NVD.out.data
}
