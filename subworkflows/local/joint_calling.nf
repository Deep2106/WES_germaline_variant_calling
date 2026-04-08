/*
    JOINT CALLING SUBWORKFLOW
    GenomicsDB (incremental) + GenotypeGVCFs + optional GLnexus
    
    IMPORTANT:
    - GATK: Uses GenomicsDB for incremental joint calling (all samples)
    - DeepVariant: Stores GVCFs persistently, GLnexus uses ALL stored GVCFs
    - Both callers produce VCFs with ALL samples for proper merging
    - Backup of GenomicsDB/PED is handled by SLURM script before pipeline runs
*/

include { GATK_GENOMICSDBIMPORT  } from '../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK_GENOTYPEGVCFS     } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GLNEXUS                } from '../../modules/nf-core/glnexus/main'
include { BCFTOOLS_MERGE_CALLERS } from '../../modules/nf-core/bcftools/merge_callers/main'
include { MERGE_PED              } from '../../modules/local/merge_ped/main'
include { COLLECT_DV_GVCFS       } from '../../modules/local/collect_dv_gvcfs/main'

workflow JOINT_CALLING {
    take:
    gvcf             // GATK GVCFs
    dv_gvcf          // DeepVariant GVCFs (optional)
    fasta
    fasta_fai
    fasta_dict
    dbsnp
    dbsnp_tbi
    intervals        // Target BED - used for GenomicsDB AND filtering DV VCF
    master_ped
    genomicsdb_path
    dv_gvcf_path     // Persistent storage for DeepVariant GVCFs
    is_update
    run_deepvariant

    main:
    ch_versions = Channel.empty()

    // =========================================================================
    // GATK Joint Calling (using GenomicsDB)
    // =========================================================================
    
    // Collect all GATK GVCFs
    ch_gvcf_files = gvcf.map { meta, gvcf, tbi -> gvcf }.collect()
    ch_tbi_files = gvcf.map { meta, gvcf, tbi -> tbi }.collect()

    // Prepare existing GenomicsDB path
    ch_existing_db = is_update ?
        Channel.value(file(genomicsdb_path)) :
        Channel.value(file("${workDir}/NO_DB"))

    // GenomicsDBImport (publishDir handles copying to persistent storage)
    GATK_GENOMICSDBIMPORT(
        ch_gvcf_files,
        ch_tbi_files,
        intervals,
        ch_existing_db,
        is_update
    )
    ch_versions = ch_versions.mix(GATK_GENOMICSDBIMPORT.out.versions)

    // GenotypeGVCFs (NO --intervals flag!)
    GATK_GENOTYPEGVCFS(
        GATK_GENOMICSDBIMPORT.out.genomicsdb,
        fasta,
        fasta_fai,
        fasta_dict,
        dbsnp,
        dbsnp_tbi,
        master_ped
    )
    ch_versions = ch_versions.mix(GATK_GENOTYPEGVCFS.out.versions)

    // Output VCF
    ch_joint_vcf = GATK_GENOTYPEGVCFS.out.vcf

    // =========================================================================
    // DeepVariant Joint Calling (using GLnexus with ALL stored GVCFs)
    // =========================================================================
    
    if (run_deepvariant) {
        // Current batch DeepVariant GVCFs
        ch_current_dv_gvcfs = dv_gvcf.map { meta, gvcf, tbi -> gvcf }.collect()
        ch_current_dv_tbis = dv_gvcf.map { meta, gvcf, tbi -> tbi }.collect()

        // Collect ALL DeepVariant GVCFs (existing + current batch)
        COLLECT_DV_GVCFS(
            ch_current_dv_gvcfs,
            ch_current_dv_tbis,
            dv_gvcf_path,
            is_update
        )
        ch_versions = ch_versions.mix(COLLECT_DV_GVCFS.out.versions)

        // GLnexus with ALL DeepVariant GVCFs
        GLNEXUS(COLLECT_DV_GVCFS.out.all_gvcfs)
        ch_versions = ch_versions.mix(GLNEXUS.out.versions)

        // Merge GATK and DeepVariant VCFs with caller tags
        // Both VCFs now have ALL samples
        BCFTOOLS_MERGE_CALLERS(
            GATK_GENOTYPEGVCFS.out.vcf,
            GLNEXUS.out.vcf,
            fasta,
            fasta_fai,
            intervals    // Target BED for filtering DeepVariant
        )
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE_CALLERS.out.versions)

        ch_joint_vcf = BCFTOOLS_MERGE_CALLERS.out.vcf
    }

    emit:
    vcf        = ch_joint_vcf.map { vcf, tbi -> vcf }
    tbi        = ch_joint_vcf.map { vcf, tbi -> tbi }
    genomicsdb = GATK_GENOMICSDBIMPORT.out.genomicsdb
    versions   = ch_versions
}
