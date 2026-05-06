/*
========================================================================================
    SPLIT VCF BY PLOIDY REGIONS
========================================================================================
    Splits a joint VCF into diploid and haploid regions:
    
    - Diploid:  autosomes (chr1-22) + chrM + PAR regions of chrX
    - Haploid:  non-PAR chrX + chrY (contains haploid male genotypes)
    
    This is needed because GATK CalculateGenotypePosteriors and 
    VariantAnnotator (PossibleDeNovo) only support diploid genotypes.
    Male non-PAR chrX/Y variants have 2 PL values which crash these tools.
    
    GRCh38 PAR regions:
      PAR1: chrX:10001-2781479
      PAR2: chrX:155701383-156030895
----------------------------------------------------------------------------------------
*/

process SPLIT_VCF_BY_PLOIDY {
    tag "split_ploidy"
    label 'process_low'

    container "${params.containers.bcftools}"

    input:
    path vcf
    path tbi

    output:
    tuple path("diploid.vcf.gz"), path("diploid.vcf.gz.tbi")           , emit: diploid
    tuple path("haploid.vcf.gz"), path("haploid.vcf.gz.tbi")           , emit: haploid
    path "versions.yml"                                                 , emit: versions

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "SPLIT VCF BY PLOIDY REGIONS"
    echo "=============================================="
    echo "Input VCF: ${vcf}"

    # =========================================================================
    # Create region files
    # =========================================================================

    # GRCh38 PAR regions on chrX (these are diploid even for males)
    cat > par.bed <<'EOF'
chrX	10000	2781479
chrX	155701382	156030895
EOF

    # Get all chromosomes in the VCF
    CHROMS=\$(tabix -l ${vcf})

    # =========================================================================
    # Build diploid regions: autosomes + chrM + PAR
    # =========================================================================
    > diploid_regions.txt
    for CHR in \${CHROMS}; do
        if [[ "\${CHR}" =~ ^chr[0-9]+\$ ]] || [[ "\${CHR}" == "chrM" ]]; then
            echo "\${CHR}" >> diploid_regions.txt
        fi
    done

    # Add PAR regions explicitly
    echo "chrX:10001-2781479" >> diploid_regions.txt
    echo "chrX:155701383-156030895" >> diploid_regions.txt

    echo "Diploid regions:"
    cat diploid_regions.txt

    # =========================================================================
    # Build haploid regions: non-PAR chrX + chrY
    # =========================================================================
    > haploid_regions.txt

    # Non-PAR chrX: between PAR1 end and PAR2 start
    if echo "\${CHROMS}" | grep -q "^chrX\$"; then
        echo "chrX:2781480-155701382" >> haploid_regions.txt
    fi

    # All of chrY
    if echo "\${CHROMS}" | grep -q "^chrY\$"; then
        echo "chrY" >> haploid_regions.txt
    fi

    echo ""
    echo "Haploid regions:"
    cat haploid_regions.txt

    # =========================================================================
    # Split the VCF
    # =========================================================================
    echo ""
    echo "Extracting diploid regions..."
    bcftools view -t \$(paste -sd, diploid_regions.txt) -Oz -o diploid.vcf.gz ${vcf}
    tabix -p vcf diploid.vcf.gz
    DIPLOID_COUNT=\$(bcftools view -H diploid.vcf.gz | wc -l)
    echo "Diploid variants: \${DIPLOID_COUNT}"

    echo "Extracting haploid regions..."
    if [ -s haploid_regions.txt ]; then
        bcftools view -t \$(paste -sd, haploid_regions.txt) -Oz -o haploid.vcf.gz ${vcf}
        tabix -p vcf haploid.vcf.gz
        HAPLOID_COUNT=\$(bcftools view -H haploid.vcf.gz | wc -l)
    else
        # No haploid regions — create empty VCF with header
        bcftools view -h ${vcf} | bgzip > haploid.vcf.gz
        tabix -p vcf haploid.vcf.gz
        HAPLOID_COUNT=0
    fi
    echo "Haploid variants: \${HAPLOID_COUNT}"

    TOTAL=\$(bcftools view -H ${vcf} | wc -l)
    echo ""
    echo "Original: \${TOTAL} variants"
    echo "Diploid:  \${DIPLOID_COUNT} variants"
    echo "Haploid:  \${HAPLOID_COUNT} variants"
    echo "Sum:      \$((\${DIPLOID_COUNT} + \${HAPLOID_COUNT})) variants"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        diploid_variants: \${DIPLOID_COUNT}
        haploid_variants: \${HAPLOID_COUNT}
    END_VERSIONS
    """
}
