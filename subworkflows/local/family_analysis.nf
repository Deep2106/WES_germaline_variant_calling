/*
========================================================================================
    FAMILY ANALYSIS SUBWORKFLOW (Ploidy-Aware)
========================================================================================
    Performs family-based analysis with ploidy-aware handling:
    
    1. SPLIT:  Separate diploid (autosomes+PAR) from haploid (non-PAR chrX/Y)
    2. CGP:    CalculateGenotypePosteriors on DIPLOID variants only
    3. DeNovo: VariantAnnotator (PossibleDeNovo) on DIPLOID variants only
    4. MERGE:  Recombine diploid (annotated) + haploid (pass-through)
    5. Peddy:  QC to verify pedigree relationships (optional)
    
    Why split? GATK CGP and VariantAnnotator only support diploid genotypes.
    Male non-PAR chrX/Y variants have haploid PLs (2 values instead of 3)
    which crash these tools. Haploid variants are passed through unchanged.
    
    De novo calling on haploid regions requires specialized hemizygous models
    which are not supported by GATK PossibleDeNovo annotation.

    Input:
        vcf         - Joint-called VCF (may contain haploid genotypes)
        tbi         - VCF index
        fasta       - Reference FASTA
        fasta_fai   - Reference FASTA index
        fasta_dict  - Reference sequence dictionary
        master_ped  - PED file with all samples

    Output:
        vcf         - VCF with de novo annotations (diploid regions)
        tbi         - VCF index
        multiqc_files - Peddy results for MultiQC
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

include { SPLIT_VCF_BY_PLOIDY                 } from '../../modules/local/split_vcf_by_ploidy/main'
include { GATK_CALCULATEGENOTYPEPOSTERIORS     } from '../../modules/nf-core/gatk4/calculategenotypeposteriors/main'
include { GATK_VARIANTANNOTATOR                } from '../../modules/nf-core/gatk4/variantannotator/main'
include { MERGE_VCF_PLOIDY                     } from '../../modules/local/merge_vcf_ploidy/main'
include { PEDDY                                } from '../../modules/nf-core/peddy/main'

workflow FAMILY_ANALYSIS {

    take:
    vcf
    tbi
    fasta
    fasta_fai
    fasta_dict
    master_ped

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // =========================================================================
    // Step 0: Split VCF into diploid and haploid regions
    // Diploid: autosomes + PAR (safe for CGP + VariantAnnotator)
    // Haploid: non-PAR chrX + chrY (skip family priors)
    // =========================================================================

    SPLIT_VCF_BY_PLOIDY(vcf, tbi)

    ch_versions = ch_versions.mix(SPLIT_VCF_BY_PLOIDY.out.versions)

    // =========================================================================
    // Step 1: Calculate Genotype Posteriors (DIPLOID ONLY)
    // Uses pedigree to refine genotype calls using transmission probabilities
    // =========================================================================

    GATK_CALCULATEGENOTYPEPOSTERIORS(
        SPLIT_VCF_BY_PLOIDY.out.diploid.map { vcf, tbi -> vcf },
        SPLIT_VCF_BY_PLOIDY.out.diploid.map { vcf, tbi -> tbi },
        fasta,
        fasta_fai,
        fasta_dict,
        master_ped
    )

    ch_versions = ch_versions.mix(GATK_CALCULATEGENOTYPEPOSTERIORS.out.versions)

    // =========================================================================
    // Step 2: Annotate De Novo Variants (DIPLOID ONLY)
    // Adds hiConfDeNovo and loConfDeNovo annotations
    // =========================================================================

    GATK_VARIANTANNOTATOR(
        GATK_CALCULATEGENOTYPEPOSTERIORS.out.vcf.map { vcf, tbi -> vcf },
        GATK_CALCULATEGENOTYPEPOSTERIORS.out.vcf.map { vcf, tbi -> tbi },
        fasta,
        fasta_fai,
        fasta_dict,
        master_ped
    )

    ch_versions = ch_versions.mix(GATK_VARIANTANNOTATOR.out.versions)

    // =========================================================================
    // Step 3: Merge diploid (annotated) + haploid (pass-through)
    // Recombines into a single VCF preserving chromosome order
    // =========================================================================

    MERGE_VCF_PLOIDY(
        GATK_VARIANTANNOTATOR.out.vcf.map { vcf, tbi -> vcf },
        GATK_VARIANTANNOTATOR.out.vcf.map { vcf, tbi -> tbi },
        SPLIT_VCF_BY_PLOIDY.out.haploid.map { vcf, tbi -> vcf },
        SPLIT_VCF_BY_PLOIDY.out.haploid.map { vcf, tbi -> tbi }
    )

    ch_versions = ch_versions.mix(MERGE_VCF_PLOIDY.out.versions)

    // =========================================================================
    // Step 4: Peddy QC (optional)
    // Verifies pedigree relationships using genetic data
    // Uses the MERGED VCF (full genome) for accurate sex check
    // =========================================================================

    if (params.run_peddy) {
        PEDDY(
            MERGE_VCF_PLOIDY.out.vcf.map { vcf, tbi -> vcf },
            MERGE_VCF_PLOIDY.out.vcf.map { vcf, tbi -> tbi },
            master_ped
        )
        ch_versions = ch_versions.mix(PEDDY.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.results)
    }

    emit:
    vcf           = MERGE_VCF_PLOIDY.out.vcf.map { vcf, tbi -> vcf }
    tbi           = MERGE_VCF_PLOIDY.out.vcf.map { vcf, tbi -> tbi }
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
