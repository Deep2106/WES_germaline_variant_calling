/*
========================================================================================
    GATK APPLYBQSR MODULE
========================================================================================
    Applies base quality score recalibration
    
    PERFORMANCE: Uses local SSD scratch for temp files
----------------------------------------------------------------------------------------
*/

process GATK_APPLYBQSR {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.containers.gatk4}"
    
    publishDir "${params.outdir}/bam/bqsr", mode: params.publish_dir_mode, failOnError: false
    
    input:
    tuple val(meta), path(bam), path(bai), path(table)
    path fasta
    path fasta_fai
    path fasta_dict
    path intervals
    
    output:
    tuple val(meta), path("*.bqsr.bam"), path("*.bqsr.bam.bai"), emit: bam
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def interval_arg = intervals.name != 'NO_FILE' ? "-L ${intervals}" : ""
    def memory_gb = task.memory.toGiga() - 4
    """
    #!/bin/bash
    set -euo pipefail
    
    # =========================================================================
    # Set up LOCAL SSD scratch for temp files
    # =========================================================================
    if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
        LOCAL_TMP="/local/scratch/\${SLURM_JOB_ID}/bqsr_${prefix}"
    elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
        LOCAL_TMP="\${TMPDIR}/bqsr_${prefix}"
    else
        LOCAL_TMP="\${PWD}/tmp_bqsr_${prefix}"
    fi
    mkdir -p \${LOCAL_TMP}
    trap "rm -rf \${LOCAL_TMP}" EXIT
    
    echo "Using temp directory: \${LOCAL_TMP}"
    
    export JAVA_OPTS="-Xmx${memory_gb}g -Djava.io.tmpdir=\${LOCAL_TMP}"
    
    gatk --java-options "\${JAVA_OPTS}" ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${prefix}.bqsr.bam \\
        --bqsr-recal-file ${table} \\
        ${interval_arg} \\
        --create-output-bam-index true \\
        --tmp-dir \${LOCAL_TMP}
    
    # GATK creates index as .bai, rename to .bam.bai for consistency
    mv ${prefix}.bqsr.bai ${prefix}.bqsr.bam.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -1 | sed 's/.*GATK v//')
    END_VERSIONS
    """
}
