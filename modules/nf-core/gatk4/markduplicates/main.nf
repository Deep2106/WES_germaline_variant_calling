/*
========================================================================================
    GATK MARKDUPLICATES MODULE
========================================================================================
    Marks/removes duplicate reads
    Outputs used by BOTH GATK (after BQSR) and DeepVariant (directly)
    
    PERFORMANCE: Uses local SSD scratch for temp files during sorting
----------------------------------------------------------------------------------------
*/

process GATK_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.containers.gatk4}"
    
    publishDir "${params.outdir}/bam/markdup", mode: params.publish_dir_mode, failOnError: false, pattern: "*.md.{bam,bam.bai}"
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.md.bam"), path("*.md.bam.bai"), emit: bam
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def remove_dups = params.remove_duplicates ? "--REMOVE_DUPLICATES true" : ""
    def memory_gb = task.memory.toGiga() - 4
    """
    #!/bin/bash
    set -euo pipefail
    
    # =========================================================================
    # Set up LOCAL SSD scratch for temp files
    # MarkDuplicates generates large temp files during coordinate sorting
    # =========================================================================
    if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
        LOCAL_TMP="/local/scratch/\${SLURM_JOB_ID}/markdup_${prefix}"
    elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
        LOCAL_TMP="\${TMPDIR}/markdup_${prefix}"
    else
        LOCAL_TMP="\${PWD}/tmp_markdup_${prefix}"
    fi
    mkdir -p \${LOCAL_TMP}
    trap "rm -rf \${LOCAL_TMP}" EXIT
    
    echo "Using temp directory: \${LOCAL_TMP}"
    
    export JAVA_OPTS="-Xmx${memory_gb}g -Djava.io.tmpdir=\${LOCAL_TMP}"
    
    gatk --java-options "\${JAVA_OPTS}" MarkDuplicates \\
        -I ${bam} \\
        -O ${prefix}.md.bam \\
        -M ${prefix}.metrics.txt \\
        ${remove_dups} \\
        --CREATE_INDEX true \\
        --TMP_DIR \${LOCAL_TMP}
    
    # GATK creates index as .bai, rename to .bam.bai for consistency
    mv ${prefix}.md.bai ${prefix}.md.bam.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -1 | sed 's/.*GATK v//')
    END_VERSIONS
    """
}
