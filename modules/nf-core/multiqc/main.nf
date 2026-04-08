/*
    MULTIQC MODULE
    Aggregate QC reports
    
    Output: Batch-specific reports in qc/multiqc/{batch_name}/
*/

process MULTIQC {
    tag "multiqc_${params.batch_name ?: 'all'}"
    label 'process_low'
    container "${params.containers.multiqc}"
    
    // Batch-specific output directory
    publishDir "${params.outdir}/qc/multiqc/${params.batch_name ?: 'default'}", mode: params.publish_dir_mode, failOnError: false
    
    input:
    path multiqc_files
    path versions
    
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data
    path "versions.yml", emit: versions
    
    script:
    def config_arg = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    def batch_name = params.batch_name ?: 'default'
    """
    echo "Generating MultiQC report for batch: ${batch_name}"
    echo "Input files:"
    ls -la
    
    multiqc . ${config_arg} \\
        --filename multiqc_report \\
        --title "WES Pipeline QC - ${batch_name}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
        batch: "${batch_name}"
    END_VERSIONS
    """
}
