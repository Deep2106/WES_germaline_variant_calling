/*
========================================================================================
    GATK VARIANTANNOTATOR MODULE
========================================================================================
    Annotate variants with additional information including de novo status
    Uses PossibleDeNovo annotation for trio analysis
    
    Input:
        vcf         - VCF file (with posteriors)
        tbi         - VCF index
        fasta       - Reference FASTA
        fasta_fai   - Reference FASTA index
        fasta_dict  - Reference sequence dictionary
        pedigree    - PED file for family relationships
    
    Output:
        vcf         - Annotated VCF with de novo flags
        tbi         - VCF index
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

process GATK_VARIANTANNOTATOR {
    tag "annotate_denovo"
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
    tuple path("*.denovo.vcf.gz"), path("*.denovo.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = vcf.baseName.replace('.vcf', '').replace('.posteriors', '')
    
    // Calculate Java heap size
    def avail_mem = task.memory ? task.memory.toGiga() - 2 : 14
    def java_opts = "-Xmx${avail_mem}g -Djava.io.tmpdir=\$TMPDIR"
    
    """
    #!/bin/bash
    set -euo pipefail
    
    # =========================================================================
    # Set up LOCAL SSD scratch for temp files
    # =========================================================================
    if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
        LOCAL_TMP="/local/scratch/\${SLURM_JOB_ID}/annotator_${prefix}"
    elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
        LOCAL_TMP="\${TMPDIR}/annotator_${prefix}"
    else
        LOCAL_TMP="\${PWD}/tmp_annotator_${prefix}"
    fi
    mkdir -p \${LOCAL_TMP}
    trap "rm -rf \${LOCAL_TMP}" EXIT
    
    export JAVA_OPTS="-Xmx${avail_mem}g -Djava.io.tmpdir=\${LOCAL_TMP}"
    
    echo "=============================================="
    echo "VariantAnnotator - De Novo Annotation"
    echo "=============================================="
    echo "Input VCF: ${vcf}"
    echo "Pedigree: ${pedigree}"
    echo "Temp dir: \${LOCAL_TMP}"
    
    # Run VariantAnnotator with PossibleDeNovo annotation
    gatk --java-options "\${JAVA_OPTS}" VariantAnnotator \\
        --variant ${vcf} \\
        --reference ${fasta} \\
        --pedigree ${pedigree} \\
        --annotation PossibleDeNovo \\
        --output ${prefix}.denovo.vcf.gz \\
        --tmp-dir \${LOCAL_TMP} \\
        ${args}
    
    # Index output
    gatk --java-options "-Xmx4g" IndexFeatureFile \\
        --input ${prefix}.denovo.vcf.gz
    
    echo "=============================================="
    echo "De novo annotation complete"
    echo "=============================================="
    
    # Count de novo candidates
    echo "Counting de novo candidates..."
    HI_CONF=\$(zcat ${prefix}.denovo.vcf.gz | grep -v "^#" | grep "hiConfDeNovo" | wc -l || echo "0")
    LO_CONF=\$(zcat ${prefix}.denovo.vcf.gz | grep -v "^#" | grep "loConfDeNovo" | grep -v "hiConfDeNovo" | wc -l || echo "0")
    echo "High-confidence de novo candidates: \${HI_CONF}"
    echo "Low-confidence de novo candidates: \${LO_CONF}"
    echo "Total de novo candidates: \$((\${HI_CONF} + \${LO_CONF}))"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/^.*(GATK) v//' | sed 's/ .*\$//')
        hiConfDeNovo: \${HI_CONF}
        loConfDeNovo: \${LO_CONF}
    END_VERSIONS
    """
    
    stub:
    def prefix = vcf.baseName.replace('.vcf', '').replace('.posteriors', '')
    """
    touch ${prefix}.denovo.vcf.gz
    touch ${prefix}.denovo.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
        hiConfDeNovo: 0
        loConfDeNovo: 0
    END_VERSIONS
    """
}
