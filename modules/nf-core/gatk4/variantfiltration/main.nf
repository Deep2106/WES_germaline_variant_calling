/*
========================================================================================
    GATK VARIANTFILTRATION MODULE
========================================================================================
    Apply hard filters to variants based on annotations
    Different filters for SNPs vs INDELs (GATK best practices, WES-tuned)

    SNP filters:
        QD < 2.0                  - Quality by Depth (depth-normalised confidence)
        FS > 60.0                 - Fisher Strand bias
        MQ < 40.0                 - Mapping Quality
        SOR > 3.0                 - Strand Odds Ratio
        MQRankSum < -12.5         - Mapping Quality Rank Sum
        ReadPosRankSum < -8.0     - Read Position Rank Sum
        QUAL < 30.0               - Raw call confidence (< 99.9%)
        --cluster-window-size 25  - Flag clusters of 3+ SNPs within 25bp window
                                    (WES-tuned: 10 too narrow, 35 too broad)

    INDEL filters:
        QD < 2.0                  - Quality by Depth
        FS > 200.0                - Fisher Strand bias (relaxed for INDELs)
        SOR > 10.0                - Strand Odds Ratio (relaxed for INDELs)
        QUAL < 30.0               - Raw call confidence
        --cluster-window-size 25  - Harmless no-op for INDELs, kept for consistency
        NOTE: ReadPosRankSum intentionally omitted for INDELs per GATK best practice
              (QD + FS + SOR already capture the same artifact signal more reliably)

    Input:
        vcf         - VCF file
        tbi         - VCF index
        fasta       - Reference FASTA
        fasta_fai   - Reference FASTA index
        fasta_dict  - Reference sequence dictionary
        type        - Variant type (SNP or INDEL)

    Output:
        vcf         - Filtered VCF (FILTER field annotated)
        tbi         - VCF index
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

process GATK_VARIANTFILTRATION {
    tag "${type}"
    label 'process_low'

    container "${params.containers.gatk4}"

    input:
    path vcf
    path tbi
    path fasta
    path fasta_fai
    path fasta_dict
    val type

    output:
    path "*.filtered.vcf.gz"     , emit: vcf
    path "*.filtered.vcf.gz.tbi" , emit: tbi
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = vcf.baseName.replace('.vcf', '')

    // Calculate Java heap size
    def avail_mem = task.memory ? task.memory.toGiga() - 1 : 7
    def java_opts = "-Xmx${avail_mem}g -Djava.io.tmpdir=\$TMPDIR"

    // SNP filter thresholds (GATK best practice defaults, WES-tuned)
    def snp_qd              = params.snp_filter_qd              ?: '2.0'
    def snp_fs              = params.snp_filter_fs              ?: '60.0'
    def snp_mq              = params.snp_filter_mq              ?: '40.0'
    def snp_sor             = params.snp_filter_sor             ?: '3.0'
    def snp_mqranksum       = params.snp_filter_mqranksum       ?: '-12.5'
    def snp_readposranksum  = params.snp_filter_readposranksum  ?: '-8.0'
    def snp_qual            = params.snp_filter_qual            ?: '30.0'

    // INDEL filter thresholds
    def indel_qd            = params.indel_filter_qd            ?: '2.0'
    def indel_fs            = params.indel_filter_fs            ?: '200.0'
    def indel_sor           = params.indel_filter_sor           ?: '10.0'
    def indel_qual          = params.indel_filter_qual          ?: '20.0'

    // Cluster window - 25bp WES-tuned (10=too narrow/WGS, 35=full WES, 25=conservative)
    def cluster_window      = params.cluster_window_size        ?: '35'

    """
    #!/bin/bash
    set -euo pipefail

    # =========================================================================
    # Set up LOCAL SSD scratch for temp files
    # =========================================================================
    if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
        LOCAL_TMP="/local/scratch/\${SLURM_JOB_ID}/filter_${type}"
    elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
        LOCAL_TMP="\${TMPDIR}/filter_${type}"
    else
        LOCAL_TMP="\${PWD}/tmp_filter_${type}"
    fi
    mkdir -p \${LOCAL_TMP}
    trap "rm -rf \${LOCAL_TMP}" EXIT

    export JAVA_OPTS="-Xmx${avail_mem}g -Djava.io.tmpdir=\${LOCAL_TMP}"

    echo "=============================================="
    echo "VariantFiltration - ${type}"
    echo "=============================================="
    echo "Input VCF       : ${vcf}"
    echo "Temp dir        : \${LOCAL_TMP}"
    echo "Cluster window  : ${cluster_window}bp"

    if [ "${type}" == "SNP" ]; then
        echo ""
        echo "Applying SNP hard filters (GATK best practices, WES-tuned):"
        echo "  QD             < ${snp_qd}"
        echo "  FS             > ${snp_fs}"
        echo "  MQ             < ${snp_mq}"
        echo "  SOR            > ${snp_sor}"
        echo "  MQRankSum      < ${snp_mqranksum}"
        echo "  ReadPosRankSum < ${snp_readposranksum}"
        echo "  QUAL           < ${snp_qual}"
        echo "  cluster-window-size ${cluster_window}bp"

        gatk --java-options "\${JAVA_OPTS}" VariantFiltration \\
            --variant ${vcf} \\
            --reference ${fasta} \\
            --output ${prefix}.filtered.vcf.gz \\
            --filter-name "QD_filter"             --filter-expression "QD < ${snp_qd}" \\
            --filter-name "FS_filter"             --filter-expression "FS > ${snp_fs}" \\
            --filter-name "MQ_filter"             --filter-expression "MQ < ${snp_mq}" \\
            --filter-name "SOR_filter"            --filter-expression "SOR > ${snp_sor}" \\
            --filter-name "MQRankSum_filter"      --filter-expression "MQRankSum < ${snp_mqranksum}" \\
            --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < ${snp_readposranksum}" \\
            --filter-name "QUAL_filter"           --filter-expression "QUAL < ${snp_qual}" \\
            --cluster-window-size ${cluster_window} \\
            --tmp-dir \${LOCAL_TMP} \\
            ${args}

    else
        echo ""
        echo "Applying INDEL hard filters (GATK best practices, WES-tuned):"
        echo "  QD   < ${indel_qd}"
        echo "  FS   > ${indel_fs}"
        echo "  SOR  > ${indel_sor}"
        echo "  QUAL < ${indel_qual}"
        echo "  cluster-window-size ${cluster_window}bp (no-op for INDELs)"
        echo "  ReadPosRankSum: intentionally omitted per GATK best practice"

        gatk --java-options "\${JAVA_OPTS}" VariantFiltration \\
            --variant ${vcf} \\
            --reference ${fasta} \\
            --output ${prefix}.filtered.vcf.gz \\
            --filter-name "QD_filter"   --filter-expression "QD < ${indel_qd}" \\
            --filter-name "FS_filter"   --filter-expression "FS > ${indel_fs}" \\
            --filter-name "SOR_filter"  --filter-expression "SOR > ${indel_sor}" \\
            --filter-name "QUAL_filter" --filter-expression "QUAL < ${indel_qual}" \\
            --cluster-window-size ${cluster_window} \\
            --tmp-dir \${LOCAL_TMP} \\
            ${args}
    fi

    # Index the output VCF
    tabix -f -p vcf ${prefix}.filtered.vcf.gz

    echo "=============================================="
    echo "VariantFiltration complete for ${type}"
    echo "=============================================="

    # Count PASS and filtered variants
    PASS_COUNT=\$(zcat ${prefix}.filtered.vcf.gz | grep -v "^#" | awk '\$7=="PASS" || \$7=="."' | wc -l)
    FILTERED_COUNT=\$(zcat ${prefix}.filtered.vcf.gz | grep -v "^#" | awk '\$7!="PASS" && \$7!="."' | wc -l)
    TOTAL_COUNT=\$(( PASS_COUNT + FILTERED_COUNT ))
    echo "PASS variants    : \${PASS_COUNT}"
    echo "Filtered variants: \${FILTERED_COUNT}"
    echo "Total variants   : \${TOTAL_COUNT}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/^.*(GATK) v//' | sed 's/ .*\$//')
        variant_type: "${type}"
        cluster_window_size: "${cluster_window}"
        pass_count: \${PASS_COUNT}
        filtered_count: \${FILTERED_COUNT}
        total_count: \${TOTAL_COUNT}
    END_VERSIONS
    """

    stub:
    def prefix = vcf.baseName.replace('.vcf', '')
    """
    touch ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.2.0
        variant_type: "${type}"
    END_VERSIONS
    """
}
