/*
========================================================================================
    MERGE VCF PLOIDY REGIONS
========================================================================================
    Merges diploid (family-analyzed) and haploid (pass-through) VCFs back
    into a single VCF after family analysis.
    
    Diploid VCF has: CGP posteriors + de novo annotations
    Haploid VCF has: original genotypes (no family priors applied)
    
    Uses bcftools concat with --allow-overlaps for safe merging.
----------------------------------------------------------------------------------------
*/

process MERGE_VCF_PLOIDY {
    tag "merge_ploidy"
    label 'process_low'

    container "${params.containers.bcftools}"

    input:
    path diploid_vcf
    path diploid_tbi
    path haploid_vcf
    path haploid_tbi

    output:
    tuple path("merged.family.vcf.gz"), path("merged.family.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                  , emit: versions

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "MERGE PLOIDY REGIONS"
    echo "=============================================="
    echo "Diploid VCF: ${diploid_vcf}"
    echo "Haploid VCF: ${haploid_vcf}"

    DIPLOID_COUNT=\$(bcftools view -H ${diploid_vcf} | wc -l)
    HAPLOID_COUNT=\$(bcftools view -H ${haploid_vcf} | wc -l)

    echo "Diploid variants (with family analysis): \${DIPLOID_COUNT}"
    echo "Haploid variants (pass-through):         \${HAPLOID_COUNT}"

    # Concatenate and sort
    bcftools concat \\
        --allow-overlaps \\
        -Oz -o unsorted.vcf.gz \\
        ${diploid_vcf} ${haploid_vcf}

    # Sort by genomic position
    bcftools sort \\
        -Oz -o merged.family.vcf.gz \\
        unsorted.vcf.gz

    tabix -p vcf merged.family.vcf.gz

    TOTAL=\$(bcftools view -H merged.family.vcf.gz | wc -l)
    echo ""
    echo "Merged total: \${TOTAL} variants"
    echo "Expected:     \$((\${DIPLOID_COUNT} + \${HAPLOID_COUNT})) variants"

    # Cleanup
    rm -f unsorted.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        total_variants: \${TOTAL}
    END_VERSIONS
    """
}
