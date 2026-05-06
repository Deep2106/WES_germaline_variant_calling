/*
    FASTP MODULE
    Read trimming and quality control
*/

process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.containers.fastp}"
    
    publishDir "${params.outdir}/qc/fastp", mode: params.publish_dir_mode, failOnError: false, pattern: "*.{html,json}"
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def r1 = reads[0]
    def r2 = reads[1]
    """
    fastp \\
        -i ${r1} \\
        -I ${r2} \\
        -o ${prefix}_R1.trimmed.fastq.gz \\
        -O ${prefix}_R2.trimmed.fastq.gz \\
        --html ${prefix}.fastp.html \\
        --json ${prefix}.fastp.json \\
        --thread ${task.cpus} \\
        --qualified_quality_phred ${params.fastp_qualified_quality} \\
        --length_required ${params.fastp_length_required} \\
        --detect_adapter_for_pe
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """
}
