/*
========================================================================================
    GATK CALCULATEGENOTYPEPOSTERIORS MODULE
========================================================================================
    Refine genotype calls using family pedigree information
    Uses transmission probabilities to improve genotype accuracy
    
    PERFORMANCE: Uses local SSD scratch for temp files
    
    Input:
        vcf         - Joint-called VCF
        tbi         - VCF index
        fasta       - Reference FASTA
        fasta_fai   - Reference FASTA index
        fasta_dict  - Reference sequence dictionary
        pedigree    - PED file for family relationships
    
    Output:
        vcf         - VCF with refined genotype posteriors
        tbi         - VCF index
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

process GATK_CALCULATEGENOTYPEPOSTERIORS {
    tag "calc_posteriors"
    label 'process_medium'
    container "${params.containers.gatk4}"
    
    input:
    path vcf
    path tbi
    path fasta
    path fasta_fai
    path fasta_dict
    path pedigree
    
    output:
    tuple path("*.posteriors.vcf.gz"), path("*.posteriors.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = vcf.baseName.replace('.vcf', '')
    
    // Calculate Java heap size
    def avail_mem = task.memory ? task.memory.toGiga() - 2 : 30
    """
    #!/bin/bash
    set -euo pipefail
    
    # =========================================================================
    # Set up LOCAL SSD scratch for temp files
    # =========================================================================
    if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
        LOCAL_TMP="/local/scratch/\${SLURM_JOB_ID}/cgp_${prefix}"
    elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
        LOCAL_TMP="\${TMPDIR}/cgp_${prefix}"
    else
        LOCAL_TMP="\${PWD}/tmp_cgp_${prefix}"
    fi
    mkdir -p \${LOCAL_TMP}
    trap "rm -rf \${LOCAL_TMP}" EXIT
    
    export JAVA_OPTS="-Xmx${avail_mem}g -Djava.io.tmpdir=\${LOCAL_TMP}"
    
    echo "=============================================="
    echo "CalculateGenotypePosteriors"
    echo "=============================================="
    echo "Input VCF: ${vcf}"
    echo "Pedigree: ${pedigree}"
    echo "Temp dir: \${LOCAL_TMP}"
    
    # Run CalculateGenotypePosteriors
    gatk --java-options "\${JAVA_OPTS}" CalculateGenotypePosteriors \\
        --variant ${vcf} \\
        --reference ${fasta} \\
        --pedigree ${pedigree} \\
        --output ${prefix}.posteriors.vcf.gz \\
        --skip-population-priors \\
        --tmp-dir \${LOCAL_TMP} \\
        ${args}
    
    # Index output
    gatk --java-options "-Xmx4g" IndexFeatureFile \\
        --input ${prefix}.posteriors.vcf.gz
    
    echo "=============================================="
    echo "CalculateGenotypePosteriors complete"
    echo "Output: ${prefix}.posteriors.vcf.gz"
    echo "=============================================="
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/^.*(GATK) v//' | sed 's/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    def prefix = vcf.baseName.replace('.vcf', '')
    """
    touch ${prefix}.posteriors.vcf.gz
    touch ${prefix}.posteriors.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
