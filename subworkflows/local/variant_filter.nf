/*
========================================================================================
    VARIANT FILTER SUBWORKFLOW
========================================================================================
    Apply GATK hard filters to variants (GATK best practices)

    Flow:
    1. BCFTOOLS_NORM    - Split multi-allelics + left-align (CRITICAL - prevents
                          VariantFiltration silently skipping filter evaluation on
                          multi-allelic records e.g. SOR=9.883 passing SOR>3.0)
    2. SelectVariants   - Split into SNPs and INDELs
    3. VariantFiltration - Apply type-specific hard filters
    4. Concat           - Merge filtered SNPs and INDELs back together

    Input:
        vcf         - VCF with de novo annotations
        tbi         - VCF index
        fasta       - Reference FASTA
        fasta_fai   - Reference FASTA index
        fasta_dict  - Reference sequence dictionary

    Output:
        vcf         - Hard-filtered VCF (SNPs + INDELs)
        tbi         - VCF index
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

include { BCFTOOLS_NORM                                } from '../../modules/local/bcftools_norm/main'
include { GATK_SELECTVARIANTS as SELECT_SNPS           } from '../../modules/nf-core/gatk4/selectvariants/main'
include { GATK_SELECTVARIANTS as SELECT_INDELS         } from '../../modules/nf-core/gatk4/selectvariants/main'
include { GATK_VARIANTFILTRATION as FILTER_SNPS        } from '../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK_VARIANTFILTRATION as FILTER_INDELS      } from '../../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_CONCAT                              } from '../../modules/nf-core/bcftools/concat/main'

workflow VARIANT_FILTER {

    take:
    vcf
    tbi
    fasta
    fasta_fai
    fasta_dict

    main:

    ch_versions = Channel.empty()

    // =========================================================================
    // Step 0: Normalise - split multi-allelics and left-align indels
    // MUST run before SelectVariants to ensure VariantFiltration evaluates
    // each allele independently. Multi-allelic records cause filter expressions
    // to be silently skipped even when SOR/FS clearly exceed thresholds.
    // Confirmed in production: chr3:11018647 SOR=9.883 FS=129.773 → FILTER=.
    // =========================================================================
    BCFTOOLS_NORM(
        vcf,
        tbi,
        fasta,
        fasta_fai
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    // =========================================================================
    // Step 1: Select SNPs and INDELs separately (from normalised VCF)
    // =========================================================================
    SELECT_SNPS(
        BCFTOOLS_NORM.out.vcf,
        BCFTOOLS_NORM.out.tbi,
        'SNP'
    )
    ch_versions = ch_versions.mix(SELECT_SNPS.out.versions)

    SELECT_INDELS(
        BCFTOOLS_NORM.out.vcf,
        BCFTOOLS_NORM.out.tbi,
        'INDEL'
    )
    ch_versions = ch_versions.mix(SELECT_INDELS.out.versions)

    // =========================================================================
    // Step 2: Apply hard filters (GATK best practices, WES-tuned)
    // SNPs:   QD, FS, MQ, SOR, MQRankSum, ReadPosRankSum, QUAL
    // INDELs: QD, FS, SOR, QUAL
    // =========================================================================

    // Modules emit vcf and tbi separately (not as tuple)
    FILTER_SNPS(
        SELECT_SNPS.out.vcf,
        SELECT_SNPS.out.tbi,
        fasta,
        fasta_fai,
        fasta_dict,
        'SNP'
    )
    ch_versions = ch_versions.mix(FILTER_SNPS.out.versions)

    FILTER_INDELS(
        SELECT_INDELS.out.vcf,
        SELECT_INDELS.out.tbi,
        fasta,
        fasta_fai,
        fasta_dict,
        'INDEL'
    )
    ch_versions = ch_versions.mix(FILTER_INDELS.out.versions)

    // =========================================================================
    // Step 3: Concatenate filtered SNPs and INDELs
    // =========================================================================

    // Collect VCFs and TBIs separately for concat
    ch_filtered_vcfs = FILTER_SNPS.out.vcf
        .mix(FILTER_INDELS.out.vcf)
        .collect()

    ch_filtered_tbis = FILTER_SNPS.out.tbi
        .mix(FILTER_INDELS.out.tbi)
        .collect()

    BCFTOOLS_CONCAT(
        ch_filtered_vcfs,
        ch_filtered_tbis
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    emit:
    vcf      = BCFTOOLS_CONCAT.out.vcf
    tbi      = BCFTOOLS_CONCAT.out.tbi
    versions = ch_versions
}
