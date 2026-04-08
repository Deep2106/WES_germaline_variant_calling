/*
========================================================================================
    GATK SELECTVARIANTS MODULE
========================================================================================
    Select specific variants from a VCF file
    Used to split SNPs and INDELs for separate filtering
    
    Input:
        vcf         - VCF file
        tbi         - VCF index
        type        - Variant type to select (SNP or INDEL)
    
    Output:
        vcf         - Selected variants VCF
        tbi         - VCF index
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

process GATK_SELECTVARIANTS {
    tag "${type}"
    label 'process_low'
    
    container "${params.containers.gatk4}"
    
    input:
    path vcf
    path tbi
    val type
    
    output:
    path "*.${type.toLowerCase()}.vcf.gz"     , emit: vcf
    path "*.${type.toLowerCase()}.vcf.gz.tbi" , emit: tbi
    path "versions.yml"                       , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = vcf.baseName.replace('.vcf', '')
    
    // Calculate Java heap size
    def avail_mem = task.memory ? task.memory.toGiga() - 1 : 7
    def java_opts = "-Xmx${avail_mem}g -Djava.io.tmpdir=\$TMPDIR"
    
    """
    #!/bin/bash
    set -euo pipefail
    
    # =========================================================================
    # Set up LOCAL SSD scratch for temp files
    # =========================================================================
    if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
        LOCAL_TMP="/local/scratch/\${SLURM_JOB_ID}/select_${type}"
    elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
        LOCAL_TMP="\${TMPDIR}/select_${type}"
    else
        LOCAL_TMP="\${PWD}/tmp_select_${type}"
    fi
    mkdir -p \${LOCAL_TMP}
    trap "rm -rf \${LOCAL_TMP}" EXIT
    
    export JAVA_OPTS="-Xmx${avail_mem}g -Djava.io.tmpdir=\${LOCAL_TMP}"
    
    echo "=============================================="
    echo "SelectVariants - ${type}"
    echo "=============================================="
    echo "Input VCF: ${vcf}"
    
    gatk --java-options "\${JAVA_OPTS}" SelectVariants \\
        --variant ${vcf} \\
        --select-type-to-include ${type} \\
        --output ${prefix}.${type.toLowerCase()}.vcf.gz \\
        --tmp-dir \${LOCAL_TMP} \\
        ${args}
    
    # Index the output VCF (GATK doesn't create index automatically)
    tabix -f -p vcf ${prefix}.${type.toLowerCase()}.vcf.gz
    
    # Count selected variants
    VARIANT_COUNT=\$(zcat ${prefix}.${type.toLowerCase()}.vcf.gz | grep -v "^#" | wc -l)
    echo "Selected ${type} count: \${VARIANT_COUNT}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/^.*(GATK) v//' | sed 's/ .*\$//')
        variant_type: "${type}"
        variant_count: \${VARIANT_COUNT}
    END_VERSIONS
    """
    
    stub:
    def prefix = vcf.baseName.replace('.vcf', '')
    """
    touch ${prefix}.${type.toLowerCase()}.vcf.gz
    touch ${prefix}.${type.toLowerCase()}.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
        variant_type: "${type}"
    END_VERSIONS
    """
}
