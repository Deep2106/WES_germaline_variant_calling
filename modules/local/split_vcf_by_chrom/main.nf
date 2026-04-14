/*
========================================================================================
    SPLIT VCF BY CHROMOSOME
========================================================================================
    Splits a multi-sample VCF into per-chromosome VCFs for parallel processing.
    Only emits chromosomes that contain variants.
----------------------------------------------------------------------------------------
*/

process SPLIT_VCF_BY_CHROM {
    tag "split_chroms"
    label 'process_low'

    container "${params.containers.bcftools}"

    input:
    path vcf
    path tbi

    output:
    path "chr*.vcf.gz"    , emit: vcfs
    path "chr*.vcf.gz.tbi", emit: tbis
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "SPLIT VCF BY CHROMOSOME"
    echo "=============================================="
    echo "Input VCF: ${vcf}"

    # Get list of chromosomes that have variants
    CHROMS=\$(tabix -l ${vcf})
    echo "Chromosomes found: \${CHROMS}"

    CHROM_COUNT=0
    for CHROM in \${CHROMS}; do
        # Skip unplaced contigs / alt contigs - keep only standard chroms
        if [[ ! "\${CHROM}" =~ ^chr[0-9XYM]+\$ ]] && [[ ! "\${CHROM}" =~ ^[0-9XYM]+\$ ]]; then
            echo "Skipping non-standard contig: \${CHROM}"
            continue
        fi

        VARIANT_COUNT=\$(bcftools view -H -r "\${CHROM}" ${vcf} | wc -l || true)
        if [ "\${VARIANT_COUNT}" -eq 0 ]; then
            echo "Skipping \${CHROM} (no variants)"
            continue
        fi

        echo "Extracting \${CHROM}..."
        bcftools view -r "\${CHROM}" -Oz -o "\${CHROM}.vcf.gz" ${vcf}
        tabix -p vcf "\${CHROM}.vcf.gz"
        CHROM_COUNT=\$((CHROM_COUNT + 1))

        echo "  \${CHROM}: \$(bcftools view -H \${CHROM}.vcf.gz | wc -l) variants"
    done

    echo "=============================================="
    echo "Split into \${CHROM_COUNT} chromosome VCFs"
    echo "=============================================="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        chromosomes_split: \${CHROM_COUNT}
    END_VERSIONS
    """
}
