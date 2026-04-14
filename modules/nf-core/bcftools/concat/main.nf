/*
========================================================================================
    BCFTOOLS CONCAT MODULE
========================================================================================
    Concatenates GATK hard-filtered SNP + INDEL VCFs into a single
    GATK-only filtered VCF, prior to merging with GLnexus output.

    Output is named ${batch_name}.gatk_filtered.vcf.gz to clearly distinguish from:
      - variants/merged/   (GATK + GLnexus merged)
      - variants/spliceai/ (SpliceAI-annotated)

    BATCH NAMESPACING: Output files and publishDir are prefixed with
    params.batch_name to prevent overwrite across batches.

    Input:
        vcfs        - VCF files to concatenate (filtered SNPs + filtered INDELs)
        tbis        - VCF index files
    Output:
        vcf         - Concatenated GATK-only filtered VCF
        tbi         - VCF index
        versions    - Software versions
----------------------------------------------------------------------------------------
*/
process BCFTOOLS_CONCAT {
    tag "concat"
    label 'process_low'
    container "${params.containers.bcftools}"

    publishDir "${params.outdir}/variants/gatk_filtered/${params.batch_name}", mode: params.publish_dir_mode, failOnError: false

    input:
    path vcfs
    path tbis

    output:
    path "${params.batch_name}.gatk_filtered.vcf.gz"     , emit: vcf
    path "${params.batch_name}.gatk_filtered.vcf.gz.tbi" , emit: tbi
    path "versions.yml"                                  , emit: versions

    script:
    def batch = params.batch_name
    """
    echo "=============================================="
    echo "BCFtools Concat - GATK-only filtered VCF"
    echo "Batch: ${batch}"
    echo "=============================================="
    echo "Input VCFs:"
    ls -la *.vcf.gz

    # Concatenate filtered SNPs + INDELs, sort, compress, index
    bcftools concat -a ${vcfs} | bcftools sort -Oz -o ${batch}.gatk_filtered.vcf.gz
    tabix -p vcf ${batch}.gatk_filtered.vcf.gz

    # Summary counts
    PASS_COUNT=\$(zcat ${batch}.gatk_filtered.vcf.gz | grep -v "^#" | awk '\$7=="PASS" || \$7=="."' | wc -l)
    FILTERED_COUNT=\$(zcat ${batch}.gatk_filtered.vcf.gz | grep -v "^#" | awk '\$7!="PASS" && \$7!="."' | wc -l)
    TOTAL_COUNT=\$(( PASS_COUNT + FILTERED_COUNT ))
    echo "PASS variants      : \${PASS_COUNT}"
    echo "Filtered (flagged) : \${FILTERED_COUNT}"
    echo "Total variants     : \${TOTAL_COUNT}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        batch: "${batch}"
        pass_count: \${PASS_COUNT}
        filtered_count: \${FILTERED_COUNT}
        total_count: \${TOTAL_COUNT}
    END_VERSIONS
    """
}
