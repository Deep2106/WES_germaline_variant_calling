/*
========================================================================================
    BCFTOOLS NORM MODULE
========================================================================================
    Split multi-allelic sites and left-align indels prior to variant filtering.

    WHY THIS IS CRITICAL:
    Multi-allelic records (e.g. C > G,T) cause GATK VariantFiltration to silently
    skip filter expression evaluation in some cases — even when SOR/FS annotations
    clearly exceed thresholds. This was confirmed in production:
        chr3:11018647  SOR=9.883  FS=129.773 → FILTER=. (passed — should have failed)
    The fix is to split all multi-allelics into biallelic records BEFORE SelectVariants
    and VariantFiltration, so each record is evaluated cleanly and independently.

    Operations (in order):
    1. bcftools norm -m -any   → split multi-allelic into biallelic records
    2. bcftools norm -f ref    → left-align and normalize indel representation
    3. bcftools norm -D        → deduplicate any identical records after splitting

    Input:
        vcf         - Joint genotyped VCF (may contain multi-allelics)
        tbi         - VCF index
        fasta       - Reference FASTA (required for left-alignment)
        fasta_fai   - Reference FASTA index

    Output:
        vcf         - Normalised biallelic VCF
        tbi         - VCF index
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

process BCFTOOLS_NORM {
    tag "norm_split_multiallelic"
    label 'process_low'
    container "${params.containers.bcftools}"

    input:
    path vcf
    path tbi
    path fasta
    path fasta_fai

    output:
    path "normalized.vcf.gz"     , emit: vcf
    path "normalized.vcf.gz.tbi" , emit: tbi
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "BCFtools Norm — Split Multi-Allelics"
    echo "=============================================="
    echo "Input VCF : ${vcf}"
    echo "Reference : ${fasta}"

    # Count multi-allelic sites before splitting
    MULTI_BEFORE=\$(bcftools view -H ${vcf} | awk -F'\\t' '\$5 ~ /,/' | wc -l)
    TOTAL_BEFORE=\$(bcftools view -H ${vcf} | wc -l)
    echo "Input     : \${TOTAL_BEFORE} total variants, \${MULTI_BEFORE} multi-allelic"

    # Step 1+2+3: Split multi-allelics, left-align, deduplicate — single pipe
    bcftools norm \\
        -m -any \\
        -f ${fasta} \\
        -D \\
        -Oz \\
        -o normalized.vcf.gz \\
        ${vcf}

    tabix -p vcf normalized.vcf.gz

    # Count after splitting
    TOTAL_AFTER=\$(bcftools view -H normalized.vcf.gz | wc -l)
    MULTI_AFTER=\$(bcftools view -H normalized.vcf.gz | awk -F'\\t' '\$5 ~ /,/' | wc -l)
    SPLIT_COUNT=\$(( TOTAL_AFTER - TOTAL_BEFORE ))

    echo ""
    echo "Output    : \${TOTAL_AFTER} total variants"
    echo "Split     : \${SPLIT_COUNT} new records from multi-allelic sites"
    echo "Remaining multi-allelic: \${MULTI_AFTER} (should be 0)"

    if [ \${MULTI_AFTER} -gt 0 ]; then
        echo "WARNING: \${MULTI_AFTER} multi-allelic records remain after normalisation" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        input_variants: \${TOTAL_BEFORE}
        multi_allelic_split: \${SPLIT_COUNT}
        output_variants: \${TOTAL_AFTER}
    END_VERSIONS
    """

    stub:
    """
    touch normalized.vcf.gz
    touch normalized.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}
