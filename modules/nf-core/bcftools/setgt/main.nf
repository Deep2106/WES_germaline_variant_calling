/*
========================================================================================
    BCFTOOLS SETGT  Mask zero-coverage genotypes
========================================================================================
    Converts DP=0 genotypes to missing (./.) after joint calling.

    Why this is needed:
    In multi-kit cohorts, samples captured with kit A have no reads in kit-B-only
    regions. GenotypeGVCFs outputs these as 0/0 (DP=0)  indistinguishable from a
    true homozygous reference. This causes:
      - False de novo calls (parents appear REF with no evidence)
      - Inflated REF allele counts distorting AF calculations
      - Corrupted posterior probability calculations in CGP

    This step runs immediately after GenotypeGVCFs and before any downstream
    analysis. It is a no-op when all samples have coverage (DP>0 everywhere).
========================================================================================
*/

process BCFTOOLS_SETGT {

    tag "${params.batch_name}"
    label 'process_medium'

    container "${params.containers.bcftools}"

    publishDir "${params.outdir}/variants/joint/${params.batch_name}",
        mode: params.publish_dir_mode,
        pattern: "*.masked.vcf.gz*"

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("${params.batch_name}.joint.masked.vcf.gz"),
          path("${params.batch_name}.joint.masked.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                        emit: versions

    script:
    """
    echo "=============================================="
    echo "BCFTOOLS SETGT - Mask DP=0 genotypes"
    echo "Batch: ${params.batch_name}"
    echo "=============================================="

    # Count DP=0 genotypes before masking
    BEFORE=\$(bcftools view -H ${vcf} | \
        bcftools query -f '[%DP\\n]' | \
        awk '\$1==0' | wc -l)
    echo "DP=0 genotypes before masking: \${BEFORE}"

    # Convert any genotype where DP=0 to missing (.)/.)
    # -t q  : query mode  apply to genotypes matching the expression
    # -n .  : set to missing
    # -i    : include filter  only touch GT where DP=0
    bcftools +setGT ${vcf} \
        -- \
        -t q \
        -n . \
        -i 'FMT/DP=0' | \
    bcftools view \
        --output-type z \
        --output ${params.batch_name}.joint.masked.vcf.gz

    bcftools index --tbi ${params.batch_name}.joint.masked.vcf.gz

    # Count after masking  should match BEFORE
    AFTER=\$(bcftools view -H ${params.batch_name}.joint.masked.vcf.gz | \
        bcftools query -f '[%GT\\n]' | \
        grep -c "^\\.\\." || true)
    echo "Missing genotypes after masking: \${AFTER}"
    echo "Expected: \${BEFORE} (should match DP=0 count above)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
