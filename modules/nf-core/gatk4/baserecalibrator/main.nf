/*
    GATK BASERECALIBRATOR MODULE
    Generates recalibration table for BQSR
*/

process GATK_BASERECALIBRATOR {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.containers.gatk4}"
    
    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path fasta_dict
    path dbsnp
    path dbsnp_tbi
    path known_indels
    path known_indels_tbi
    path intervals
    
    output:
    tuple val(meta), path("*.recal.table"), emit: table
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def interval_arg = intervals.name != 'NO_FILE' ? "-L ${intervals}" : ""
    def memory_gb = task.memory.toGiga() - 4
    """
    export JAVA_OPTS="-Xmx${memory_gb}g"
    
    gatk --java-options "\${JAVA_OPTS}" BaseRecalibrator \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${prefix}.recal.table \\
        --known-sites ${dbsnp} \\
        --known-sites ${known_indels} \\
        ${interval_arg} \\
        --tmp-dir \${TMPDIR:-/tmp}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -1 | sed 's/.*GATK v//')
    END_VERSIONS
    """
}
