/*
    BCFTOOLS MERGE CALLERS MODULE
    Merges GATK and DeepVariant VCFs with caller concordance tagging

    IMPORTANT: DeepVariant VCF is filtered to target regions before merging
    (DeepVariant runs without --regions for performance, so we filter here)

    Pre-merge normalization (per-VCF, before concordance comparison):
    1. bcftools norm -m -any   → split multi-allelics into biallelic records
    2. bcftools norm -f ref    → left-align and normalize indels
    3. bcftools view -e '*'    → remove spanning deletion placeholders (ALT=*)

    This fixes ANNOVAR convert2annovar.pl limitations:
    - Multi-allelic sites caused wrong Ref/Alt and wrong zygosity
    - Spanning deletions (*) were converted to Alt=0 (phantom deletions)
    - Inconsistent indel representation caused database lookup mismatches

    Normalization BEFORE merge also improves concordance accuracy — same
    variant won't be represented differently by each caller.

    Output tags:
    - CALLER=BOTH: variant called by both GATK and DeepVariant
    - CALLER=GATK_ONLY: variant called only by GATK
    - CALLER=DV_ONLY: variant called only by DeepVariant

    BATCH NAMESPACING: Output files and publishDir are prefixed with
    params.batch_name to prevent overwrite across batches.
*/

process BCFTOOLS_MERGE_CALLERS {
    tag "merge_callers"
    label 'process_medium'
    container "${params.containers.bcftools}"

    publishDir "${params.outdir}/variants/merged/${params.batch_name}", mode: params.publish_dir_mode, failOnError: false

    input:
    tuple path(gatk_vcf), path(gatk_tbi)
    tuple path(dv_vcf), path(dv_tbi)
    path fasta
    path fasta_fai
    path targets_bed   // Target regions BED file for filtering DeepVariant

    output:
    tuple path("${params.batch_name}.merged.caller_tagged.vcf.gz"), path("${params.batch_name}.merged.caller_tagged.vcf.gz.tbi"), emit: vcf
    path "${params.batch_name}.concordance_stats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    def batch          = params.batch_name
    def filter_regions = targets_bed.name != 'NO_FILE' ? true : false
    """
    echo "=== Caller Concordance Analysis ==="
    echo "Batch: ${batch}"

    # =========================================================================
    # GATK VCF: Normalize → Remove spanning deletions
    # =========================================================================
    echo ""
    echo "=== Normalizing GATK VCF ==="
    GATK_RAW=\$(bcftools view -H ${gatk_vcf} | wc -l)
    GATK_STAR=\$(bcftools view -H ${gatk_vcf} | awk -F'\\t' '\$5 ~ /\\*/' | wc -l)
    GATK_MULTI=\$(bcftools view -H ${gatk_vcf} | awk -F'\\t' '\$5 ~ /,/' | wc -l)

    # Step 1+2: Split multi-allelics and left-align in one pass
    bcftools norm -f ${fasta} -m -any ${gatk_vcf} -Oz -o gatk.split.vcf.gz

    # Step 3: Remove spanning deletion placeholders (ALT=*)
    bcftools view -e 'ALT="*"' gatk.split.vcf.gz -Oz -o gatk.norm.vcf.gz
    tabix -p vcf gatk.norm.vcf.gz

    GATK_FINAL=\$(bcftools view -H gatk.norm.vcf.gz | wc -l)
    echo "GATK: \${GATK_RAW} raw → \${GATK_FINAL} normalized (removed \${GATK_STAR} spanning dels, split \${GATK_MULTI} multi-allelics)"
    rm -f gatk.split.vcf.gz

    # =========================================================================
    # DeepVariant VCF: Filter targets → Normalize → Remove spanning deletions
    # =========================================================================
    echo ""
    echo "=== Normalizing DeepVariant VCF ==="

    if [ "${filter_regions}" == "true" ]; then
        echo "Filtering DeepVariant VCF to target regions: ${targets_bed}"
        BEFORE=\$(bcftools view -H ${dv_vcf} | wc -l)
        bcftools view -R ${targets_bed} ${dv_vcf} -Oz -o dv.filtered.vcf.gz
        tabix -p vcf dv.filtered.vcf.gz
        AFTER=\$(bcftools view -H dv.filtered.vcf.gz | wc -l)
        echo "DeepVariant variants: \${BEFORE} → \${AFTER} (filtered to targets)"

        DV_STAR=\$(bcftools view -H dv.filtered.vcf.gz | awk -F'\\t' '\$5 ~ /\\*/' | wc -l)
        DV_MULTI=\$(bcftools view -H dv.filtered.vcf.gz | awk -F'\\t' '\$5 ~ /,/' | wc -l)

        bcftools norm -f ${fasta} -m -any dv.filtered.vcf.gz -Oz -o dv.split.vcf.gz
    else
        echo "No target BED provided - using all DeepVariant variants"
        DV_STAR=\$(bcftools view -H ${dv_vcf} | awk -F'\\t' '\$5 ~ /\\*/' | wc -l)
        DV_MULTI=\$(bcftools view -H ${dv_vcf} | awk -F'\\t' '\$5 ~ /,/' | wc -l)

        bcftools norm -f ${fasta} -m -any ${dv_vcf} -Oz -o dv.split.vcf.gz
    fi

    # Step 3: Remove spanning deletion placeholders (ALT=*)
    bcftools view -e 'ALT="*"' dv.split.vcf.gz -Oz -o dv.norm.vcf.gz
    tabix -p vcf dv.norm.vcf.gz

    DV_FINAL=\$(bcftools view -H dv.norm.vcf.gz | wc -l)
    echo "DV: \${DV_STAR} spanning dels removed, \${DV_MULTI} multi-allelics split → \${DV_FINAL} final"
    rm -f dv.split.vcf.gz dv.filtered.vcf.gz 2>/dev/null || true

    # =========================================================================
    # Concordance: Find intersection between normalized VCFs
    # =========================================================================
    echo ""
    echo "=== Concordance Analysis ==="

    bcftools isec -p isec_dir gatk.norm.vcf.gz dv.norm.vcf.gz

    GATK_ONLY=\$(bcftools view -H isec_dir/0000.vcf 2>/dev/null | wc -l)
    DV_ONLY=\$(bcftools view -H isec_dir/0001.vcf 2>/dev/null | wc -l)
    BOTH_GATK=\$(bcftools view -H isec_dir/0002.vcf 2>/dev/null | wc -l)

    echo "GATK only: \${GATK_ONLY}"        >  ${batch}.concordance_stats.txt
    echo "DeepVariant only: \${DV_ONLY}"   >> ${batch}.concordance_stats.txt
    echo "Both callers: \${BOTH_GATK}"     >> ${batch}.concordance_stats.txt
    TOTAL=\$((\${GATK_ONLY} + \${DV_ONLY} + \${BOTH_GATK}))
    if [ \${TOTAL} -gt 0 ]; then
        echo "Concordance rate: \$(echo "scale=2; \${BOTH_GATK} * 100 / \${TOTAL}" | bc)%" >> ${batch}.concordance_stats.txt
    fi

    echo ""                                                           >> ${batch}.concordance_stats.txt
    echo "=== Normalization Stats ==="                                >> ${batch}.concordance_stats.txt
    echo "GATK: \${GATK_RAW} raw, \${GATK_MULTI} multi-allelic, \${GATK_STAR} spanning del (*), \${GATK_FINAL} final" >> ${batch}.concordance_stats.txt
    echo "DV: \${DV_MULTI} multi-allelic, \${DV_STAR} spanning del (*), \${DV_FINAL} final"                            >> ${batch}.concordance_stats.txt

    cat ${batch}.concordance_stats.txt

    # =========================================================================
    # Tag and merge
    # =========================================================================
    echo ""
    echo "=== Tagging and Merging ==="

    echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller source">' > header.txt

    bcftools annotate -h header.txt isec_dir/0000.vcf 2>/dev/null | \\
        awk 'BEGIN{OFS="\\t"} /^#/{print} !/^#/{\$8="CALLER=GATK_ONLY;"\$8; print}' | \\
        bgzip -c > gatk_only_tagged.vcf.gz
    tabix -p vcf gatk_only_tagged.vcf.gz

    bcftools annotate -h header.txt isec_dir/0001.vcf 2>/dev/null | \\
        awk 'BEGIN{OFS="\\t"} /^#/{print} !/^#/{\$8="CALLER=DV_ONLY;"\$8; print}' | \\
        bgzip -c > dv_only_tagged.vcf.gz
    tabix -p vcf dv_only_tagged.vcf.gz

    bcftools annotate -h header.txt isec_dir/0002.vcf 2>/dev/null | \\
        awk 'BEGIN{OFS="\\t"} /^#/{print} !/^#/{\$8="CALLER=BOTH;"\$8; print}' | \\
        bgzip -c > both_tagged.vcf.gz
    tabix -p vcf both_tagged.vcf.gz

    bcftools concat -a -Oz -o merged.unsorted.vcf.gz \\
        gatk_only_tagged.vcf.gz \\
        dv_only_tagged.vcf.gz \\
        both_tagged.vcf.gz

    bcftools sort -Oz -o ${batch}.merged.caller_tagged.vcf.gz merged.unsorted.vcf.gz
    tabix -p vcf ${batch}.merged.caller_tagged.vcf.gz

    FINAL_COUNT=\$(bcftools view -H ${batch}.merged.caller_tagged.vcf.gz | wc -l)
    echo ""
    echo "=== Final Output ==="
    echo "Batch         : ${batch}"
    echo "Total merged  : \${FINAL_COUNT}"
    echo "  GATK only   : \${GATK_ONLY}"
    echo "  DV only     : \${DV_ONLY}"
    echo "  Both        : \${BOTH_GATK}"
    ls -lh ${batch}.merged.caller_tagged.vcf.gz

    # Cleanup intermediates
    rm -f merged.unsorted.vcf.gz gatk.norm.vcf.gz* dv.norm.vcf.gz*
    rm -f gatk_only_tagged.vcf.gz* dv_only_tagged.vcf.gz* both_tagged.vcf.gz*
    rm -rf isec_dir header.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        batch: "${batch}"
    END_VERSIONS
    """
}
