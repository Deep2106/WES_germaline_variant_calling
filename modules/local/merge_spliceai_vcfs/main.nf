/*
========================================================================================
    MERGE SPLICEAI VCFS
========================================================================================
    Merges per-chromosome SpliceAI-annotated VCFs back into a single VCF.
    Accepts UNCOMPRESSED per-chrom VCFs, compresses and indexes the merged output.
    Uses bcftools concat with natural chromosome ordering.

    BATCH NAMESPACING: Output files and publishDir are prefixed with
    params.batch_name to prevent overwrite across batches.
----------------------------------------------------------------------------------------
*/

process MERGE_SPLICEAI_VCFS {
    tag "merge_spliceai"
    label 'process_low'

    container "${params.containers.bcftools}"

    publishDir "${params.outdir}/variants/spliceai/${params.batch_name}", mode: 'copy', pattern: '*.spliceai.merged.*'
    publishDir "${params.outdir}/variants/spliceai/${params.batch_name}", mode: 'copy', pattern: '*.spliceai_summary.txt'

    input:
    path vcfs
    path summaries

    output:
    path "${params.batch_name}.filtered.spliceai.merged.vcf.gz"    , emit: vcf
    path "${params.batch_name}.filtered.spliceai.merged.vcf.gz.tbi", emit: tbi
    path "${params.batch_name}.spliceai_summary.txt"               , emit: summary
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def batch = params.batch_name
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "MERGE SPLICEAI CHROMOSOME VCFS"
    echo "Batch: ${batch}"
    echo "=============================================="

    # Compress and index each per-chromosome VCF
    for vcf_file in *.spliceai.vcf; do
        echo "Compressing \${vcf_file}..."
        bgzip -c "\${vcf_file}" > "\${vcf_file}.gz"
        tabix -p vcf "\${vcf_file}.gz"
    done

    # List all compressed VCFs and sort naturally by chromosome
    ls -1 *.spliceai.vcf.gz | sort -V > vcf_list.txt

    echo "VCFs to merge:"
    cat vcf_list.txt
    echo ""
    echo "Total files: \$(wc -l < vcf_list.txt)"

    # Concatenate in chromosome order
    bcftools concat \\
        --file-list vcf_list.txt \\
        --allow-overlaps \\
        -Oz -o ${batch}.filtered.spliceai.merged.vcf.gz

    tabix -p vcf ${batch}.filtered.spliceai.merged.vcf.gz

    # Verify output
    TOTAL_VARIANTS=\$(bcftools view -H ${batch}.filtered.spliceai.merged.vcf.gz | wc -l)
    echo ""
    echo "Merged VCF: \${TOTAL_VARIANTS} total variants"

    # Count variants with SpliceAI annotation
    ANNOTATED=\$(bcftools view -H ${batch}.filtered.spliceai.merged.vcf.gz | grep -c "SpliceAI=" || echo "0")
    echo "Variants with SpliceAI scores: \${ANNOTATED}"

    # Combine per-chrom summaries into batch-named summary file
    # CRITICAL: use chr*.spliceai_summary.txt NOT *.spliceai_summary.txt
    # The wildcard *.spliceai_summary.txt would match the output file itself
    # (CIINDI_BATCH1.spliceai_summary.txt) causing an infinite read/write loop
    # and a runaway file growing to hundreds of GB.
    {
        echo "=============================================="
        echo "SPLICEAI MERGED SUMMARY"
        echo "Batch: ${batch}"
        echo "=============================================="
        echo "Total variants: \${TOTAL_VARIANTS}"
        echo "Variants with SpliceAI annotation: \${ANNOTATED}"
        echo ""
        echo "--- Per-Chromosome Details ---"
        cat chr*.spliceai_summary.txt
    } > ${batch}.spliceai_summary.txt

    cat ${batch}.spliceai_summary.txt

    # Cleanup
    rm -f vcf_list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        batch: "${batch}"
        total_variants: \${TOTAL_VARIANTS}
        annotated_variants: \${ANNOTATED}
    END_VERSIONS
    """
}
