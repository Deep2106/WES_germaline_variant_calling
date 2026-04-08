/*
========================================================================================
    GATK HAPLOTYPECALLER MODULE
========================================================================================
    Per-sample GVCF generation for joint calling.
    
    Supports sex-aware ploidy:
      - ploidy=2 for autosomes + PAR (all samples)
      - ploidy=1 for non-PAR chrX + chrY (males only)
    
    The ploidy parameter controls --ploidy flag (default: 2).
    Different interval BEDs are passed for diploid vs haploid regions.
----------------------------------------------------------------------------------------
*/

process GATK_HAPLOTYPECALLER {
    tag "${meta.id}_ploidy${ploidy}"
    label 'process_high'

    container "${params.containers.gatk4}"

    publishDir "${params.outdir}/variants/gvcf", mode: params.publish_dir_mode, failOnError: false, pattern: "*.g.vcf.gz*"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path fasta_dict
    path dbsnp
    path dbsnp_tbi
    path intervals
    val ploidy

    output:
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf
    path "versions.yml", emit: versions

    script:
    def prefix       = "${meta.id}"
    def suffix       = ploidy == 1 ? ".haploid" : ".diploid"
    def interval_arg = intervals.name != 'NO_FILE' ? "-L ${intervals}" : ""
    def dbsnp_arg    = dbsnp.name != 'NO_FILE' ? "-D ${dbsnp}" : ""
    def memory_gb    = task.memory.toGiga() - 4
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "GATK HaplotypeCaller - ${prefix}"
    echo "=============================================="
    echo "Ploidy: ${ploidy}"
    echo "Intervals: ${intervals}"
    echo "=============================================="

    export JAVA_OPTS="-Xmx${memory_gb}g -XX:+UseParallelGC -XX:ParallelGCThreads=2"

    gatk --java-options "\${JAVA_OPTS}" HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${prefix}${suffix}.g.vcf.gz \\
        -ERC GVCF \\
        --ploidy ${ploidy} \\
        ${interval_arg} \\
        ${dbsnp_arg} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir \${TMPDIR:-/tmp}

    echo "Variants called: \$(zcat ${prefix}${suffix}.g.vcf.gz | grep -v '^#' | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -1 | sed 's/.*GATK v//')
        ploidy: ${ploidy}
        sample: ${prefix}
    END_VERSIONS
    """
}
