/*
========================================================================================
    BWA-MEM2 ALIGNMENT MODULE
========================================================================================
    Fast alignment with piped samtools sort
    
    NOTE: Uses original file paths since Singularity binds /home/data
    Index files must exist alongside FASTA or in bwa_index directory
    
    PERFORMANCE: Uses local SSD scratch for samtools sort temp files
----------------------------------------------------------------------------------------
*/

process BWAMEM2_MEM {
    tag "${meta.id}"
    label 'process_high'
    container "${params.containers.bwamem2}"
    
    input:
    tuple val(meta), path(reads)
    val fasta          // Use val - keep original path accessible via Singularity bind
    val bwa_index      // Use val - keep original path accessible via Singularity bind
    
    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${params.library_id}\\tPL:${params.platform}\\tPU:${params.platform_unit}"
    def sort_memory = (task.memory.toGiga() * 0.3).intValue()
    """
    #!/bin/bash
    set -euo pipefail
    
    # =========================================================================
    # Set up LOCAL SSD scratch for samtools sort temp files
    # Sorting generates many temp files - local SSD is much faster than NFS
    # =========================================================================
    if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
        LOCAL_TMP="/local/scratch/\${SLURM_JOB_ID}/bwa_${prefix}"
    elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
        LOCAL_TMP="\${TMPDIR}/bwa_${prefix}"
    else
        LOCAL_TMP="\${PWD}/tmp_bwa_${prefix}"
    fi
    mkdir -p \${LOCAL_TMP}
    trap "rm -rf \${LOCAL_TMP}" EXIT
    
    echo "Using temp directory for sorting: \${LOCAL_TMP}"
    
    # Determine index prefix
    # BWA-MEM2 index files should be: {prefix}.0123, .amb, .ann, .bwt.2bit.64, .pac
    
    # Check if index exists alongside FASTA
    if [ -f "${fasta}.bwt.2bit.64" ]; then
        INDEX_PREFIX="${fasta}"
        echo "Using index alongside FASTA: \${INDEX_PREFIX}"
    # Check in bwa_index directory
    elif [ -d "${bwa_index}" ]; then
        FASTA_NAME=\$(basename ${fasta})
        if [ -f "${bwa_index}/\${FASTA_NAME}.bwt.2bit.64" ]; then
            INDEX_PREFIX="${bwa_index}/\${FASTA_NAME}"
            echo "Using index from bwa_index dir: \${INDEX_PREFIX}"
        else
            # Try to find any index in the directory
            FOUND_INDEX=\$(find ${bwa_index} -name "*.bwt.2bit.64" 2>/dev/null | head -1)
            if [ -n "\${FOUND_INDEX}" ]; then
                INDEX_PREFIX=\${FOUND_INDEX%.bwt.2bit.64}
                echo "Found index: \${INDEX_PREFIX}"
            else
                echo "ERROR: No BWA-MEM2 index found!"
                echo "Looked for: ${fasta}.bwt.2bit.64"
                echo "And in: ${bwa_index}/"
                exit 1
            fi
        fi
    else
        INDEX_PREFIX="${fasta}"
        echo "Using FASTA path as index: \${INDEX_PREFIX}"
    fi
    
    echo "INDEX_PREFIX: \${INDEX_PREFIX}"
    ls -la \${INDEX_PREFIX}* 2>/dev/null || echo "Warning: Cannot list index files"
    
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -M \\
        -I 200,100 \\
        -B 4 \\
        -A 1 \\
        -w 100 \\
        -k 19 \\
        -R "${rg}" \\
        \${INDEX_PREFIX} \\
        ${reads} \\
    | samtools sort \\
        -@ 2 \\
        -m ${sort_memory}G \\
        -T \${LOCAL_TMP}/sort_tmp \\
        -o ${prefix}.sorted.bam \\
        -
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
