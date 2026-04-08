/*
    FASTQC MODULE
    Quality control for FASTQ files
*/

process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    container "${params.containers.fastqc}"
    
    publishDir "${params.outdir}/qc/fastqc", mode: params.publish_dir_mode, failOnError: false
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path "versions.yml", emit: versions
    
    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
