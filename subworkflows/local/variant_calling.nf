/*
========================================================================================
    VARIANT CALLING SUBWORKFLOW (Sex-Aware Ploidy)
========================================================================================
    GATK HaplotypeCaller + optional DeepVariant with proper ploidy handling.

    GATK Strategy:
      Males (sex=1):
        - HC run 1: autosomes + PAR regions, ploidy=2
        - HC run 2: non-PAR chrX + chrY, ploidy=1
        - Merge both GVCFs per sample
      Females (sex=2):
        - HC run 1: all targets except chrY, ploidy=2

    DeepVariant Strategy:
      Males:   --haploid_contigs chrX,chrY --par_regions_bed (created inline)
      Females: standard diploid calling

    IMPORTANT:
    - HaplotypeCaller uses BQSR BAM
    - DeepVariant uses MarkDup BAM (no BQSR — DV learns its own quality model)
    - DeepVariant creates PAR BED inline to avoid symlink/permission issues
----------------------------------------------------------------------------------------
*/

include { PREPARE_PLOIDY_BEDS                       } from '../../modules/local/prepare_ploidy_beds/main'
include { GATK_HAPLOTYPECALLER as HC_MALE_DIPLOID   } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK_HAPLOTYPECALLER as HC_MALE_HAPLOID   } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK_HAPLOTYPECALLER as HC_FEMALE         } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK_MERGE_GVCFS                          } from '../../modules/nf-core/gatk4/mergegvcfs/main'
include { DEEPVARIANT                               } from '../../modules/nf-core/deepvariant/main'

workflow VARIANT_CALLING {

    take:
    bam_bqsr            // [meta, bam, bai] — BQSR BAM for GATK
    bam_markdup         // [meta, bam, bai] — MarkDup BAM for DeepVariant
    fasta               // path: reference FASTA
    fasta_fai           // path: reference FASTA index
    fasta_dict          // path: reference dict
    dbsnp               // path: dbSNP VCF
    dbsnp_tbi           // path: dbSNP VCF index
    intervals           // path: original targets BED
    run_deepvariant     // boolean

    main:
    ch_versions = Channel.empty()

    /*
    ========================================================================
        STEP 0: Prepare sex-specific interval BED files (runs ONCE)
    ========================================================================
    */

    // PAR BED file — shipped with pipeline in assets/
    ch_par_bed = Channel.value(file("${projectDir}/assets/grch38_par.bed"))

    PREPARE_PLOIDY_BEDS(intervals, ch_par_bed)

    ch_autosomal_par_bed   = PREPARE_PLOIDY_BEDS.out.autosomal_par
    ch_sexchrom_nonpar_bed = PREPARE_PLOIDY_BEDS.out.sexchrom_nonpar
    ch_no_chrY_bed         = PREPARE_PLOIDY_BEDS.out.no_chrY
    ch_versions = ch_versions.mix(PREPARE_PLOIDY_BEDS.out.versions)

    /*
    ========================================================================
        STEP 1: Branch samples by sex for GATK HaplotypeCaller
    ========================================================================
    */

    bam_bqsr
        .branch {
            male:   it[0].sex == '1' || it[0].sex == 1
            female: true    // sex=2 or unknown defaults to diploid
        }
        .set { ch_bam_by_sex }

    /*
    ========================================================================
        STEP 2a: GATK HC — Males (two runs per sample)
        Run 1: Autosomes + PAR (ploidy=2)
        Run 2: Non-PAR chrX + chrY (ploidy=1)
    ========================================================================
    */

    HC_MALE_DIPLOID(
        ch_bam_by_sex.male,
        fasta,
        fasta_fai,
        fasta_dict,
        dbsnp,
        dbsnp_tbi,
        ch_autosomal_par_bed,
        2   // ploidy
    )
    ch_versions = ch_versions.mix(HC_MALE_DIPLOID.out.versions.first())

    HC_MALE_HAPLOID(
        ch_bam_by_sex.male,
        fasta,
        fasta_fai,
        fasta_dict,
        dbsnp,
        dbsnp_tbi,
        ch_sexchrom_nonpar_bed,
        1   // ploidy
    )
    ch_versions = ch_versions.mix(HC_MALE_HAPLOID.out.versions.first())

    /*
    ========================================================================
        STEP 2b: Merge male diploid + haploid GVCFs per sample
    ========================================================================
    */

    // Join diploid and haploid GVCFs by sample ID
    // HC_MALE_DIPLOID.out.gvcf: [meta, diploid.g.vcf.gz, diploid.g.vcf.gz.tbi]
    // HC_MALE_HAPLOID.out.gvcf: [meta, haploid.g.vcf.gz, haploid.g.vcf.gz.tbi]

    ch_male_to_merge = HC_MALE_DIPLOID.out.gvcf
        .join(HC_MALE_HAPLOID.out.gvcf, by: 0)     // join on meta
        // Result: [meta, diploid.gvcf, diploid.tbi, haploid.gvcf, haploid.tbi]

    GATK_MERGE_GVCFS(
        ch_male_to_merge,
        fasta_dict
    )
    ch_versions = ch_versions.mix(GATK_MERGE_GVCFS.out.versions.first())

    /*
    ========================================================================
        STEP 2c: GATK HC — Females (single run, exclude chrY)
    ========================================================================
    */

    HC_FEMALE(
        ch_bam_by_sex.female,
        fasta,
        fasta_fai,
        fasta_dict,
        dbsnp,
        dbsnp_tbi,
        ch_no_chrY_bed,
        2   // ploidy
    )
    ch_versions = ch_versions.mix(HC_FEMALE.out.versions.first())

    /*
    ========================================================================
        STEP 3: Combine all GVCFs (males merged + females)
    ========================================================================
    */

    ch_gvcf = GATK_MERGE_GVCFS.out.gvcf
        .mix(HC_FEMALE.out.gvcf)

    /*
    ========================================================================
        STEP 4: DeepVariant (optional, sex-aware)
        Males:   --haploid_contigs chrX,chrY --par_regions_bed (inline)
        Females: standard diploid
        PAR BED is created inline in the DV module (avoids symlink issues)
    ========================================================================
    */

    ch_dv_gvcf = Channel.empty()

    if (run_deepvariant) {
        DEEPVARIANT(
            bam_markdup,
            fasta,
            fasta_fai
        )
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())
        ch_dv_gvcf = DEEPVARIANT.out.gvcf
    }

    emit:
    gvcf             = ch_gvcf          // [meta, gvcf, tbi] — merged per sample
    deepvariant_gvcf = ch_dv_gvcf       // [meta, gvcf, tbi] — DV output
    versions         = ch_versions
}
