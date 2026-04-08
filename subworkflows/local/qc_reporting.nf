/*
    QC REPORTING SUBWORKFLOW
    MultiQC report generation
*/

include { MULTIQC } from '../../modules/nf-core/multiqc/main'

workflow QC_REPORTING {
    
    take:
    multiqc_files
    versions
    
    main:
    
    MULTIQC(
        multiqc_files,
        versions
    )
    
    emit:
    report   = MULTIQC.out.report
    versions = MULTIQC.out.versions
}
