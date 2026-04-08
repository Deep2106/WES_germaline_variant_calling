/*
========================================================================================
    GATK MERGE GVCFS
========================================================================================
    Merges per-sample diploid + haploid GVCFs into a single GVCF.
    Used for male samples where HC runs separately on autosomal (ploidy=2)
    and sex chromosome (ploidy=1) regions.
----------------------------------------------------------------------------------------
*/

process GATK_MERGE_GVCFS {
    tag "${meta.id}"
    label 'process_medium'

    container "${params.containers.gatk4}"

    publishDir "${params.outdir}/variants/gvcf", mode: params.publish_dir_mode, failOnError: false, pattern: "*.g.vcf.gz*"

    input:
    tuple val(meta), path(diploid_gvcf), path(diploid_tbi), path(haploid_gvcf), path(haploid_tbi)
    path fasta_dict

    output:
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf
    path "versions.yml", emit: versions

    script:
    def prefix    = "${meta.id}"
    def memory_gb = task.memory.toGiga() - 2
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "MERGE GVCFS - ${prefix}"
    echo "=============================================="
    echo "Diploid GVCF: ${diploid_gvcf}"
    echo "Haploid GVCF: ${haploid_gvcf}"

    export JAVA_OPTS="-Xmx${memory_gb}g"

    gatk --java-options "\${JAVA_OPTS}" MergeVcfs \\
        -I ${diploid_gvcf} \\
        -I ${haploid_gvcf} \\
        -O ${prefix}.g.vcf.gz \\
        -D ${fasta_dict} \\
        --TMP_DIR \${TMPDIR:-/tmp}

    echo "Merged variants: \$(zcat ${prefix}.g.vcf.gz | grep -v '^#' | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -1 | sed 's/.*GATK v//')
    END_VERSIONS
    """
}
