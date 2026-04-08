/*
========================================================================================
    ANNOTATION SUBWORKFLOW
========================================================================================
    Annotate variants with functional information - PER SAMPLE
    
    Steps:
        1. Split VCF by sample (only variants where sample has ALT)
        2. Run ANNOVAR for each sample separately
    
    Inputs:
        vcf         - Filtered VCF file (multi-sample)
        tbi         - VCF index
        sample_meta - Channel of sample metadata [meta]
        annovar_db  - ANNOVAR database directory
    
    Outputs:
        txt         - Annotated variants table (per sample)
        vcf         - Annotated VCF (per sample)
        csv         - Annotated CSV (per sample)
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

include { SPLIT_VCF_BY_SAMPLE } from '../../modules/local/split_vcf_by_sample/main'
include { ANNOVAR             } from '../../modules/nf-core/annovar/main'

workflow ANNOTATION {
    
    take:
    vcf          // Path: filtered VCF (multi-sample)
    tbi          // Path: VCF index  
    sample_meta  // Channel: [meta] for each sample to annotate
    annovar_db   // Path: ANNOVAR database directory
    
    main:
    
    ch_versions = Channel.empty()
    
    /*
    ============================================================================
        STEP 1: Split VCF by Sample
    ============================================================================
    */
    
    // Collect sample IDs from meta
    ch_sample_ids = sample_meta
        .map { meta -> meta.sample_id }
        .collect()
    
    // Split VCF into per-sample VCFs
    SPLIT_VCF_BY_SAMPLE(
        vcf,
        tbi,
        ch_sample_ids
    )
    
    ch_versions = ch_versions.mix(SPLIT_VCF_BY_SAMPLE.out.versions)
    
    // Create channel of [meta, vcf, tbi] for each sample
    // Match VCF files back to their metadata
    ch_per_sample_vcf = SPLIT_VCF_BY_SAMPLE.out.vcfs
        .flatten()
        .map { vcf_file ->
            def sample_id = vcf_file.baseName.replace('.per_sample.vcf', '')
            return [sample_id, vcf_file]
        }
    
    ch_per_sample_tbi = SPLIT_VCF_BY_SAMPLE.out.tbis
        .flatten()
        .map { tbi_file ->
            def sample_id = tbi_file.baseName.replace('.per_sample.vcf.gz', '')
            return [sample_id, tbi_file]
        }
    
    // Join VCF and TBI by sample_id
    ch_vcf_tbi = ch_per_sample_vcf
        .join(ch_per_sample_tbi)
        .map { sample_id, vcf_file, tbi_file ->
            return [sample_id, vcf_file, tbi_file]
        }
    
    // Add metadata back
    ch_sample_meta_map = sample_meta
        .map { meta -> [meta.sample_id, meta] }
    
    ch_annotate_input = ch_vcf_tbi
        .join(ch_sample_meta_map)
        .map { sample_id, vcf_file, tbi_file, meta ->
            return [meta, vcf_file, tbi_file]
        }
    
    // Log samples to annotate
    ch_annotate_input
        .map { meta, vcf_file, tbi_file ->
            log.info "Will annotate sample: ${meta.sample_id}"
            return [meta, vcf_file, tbi_file]
        }
    
    /*
    ============================================================================
        STEP 2: ANNOVAR Annotation (Per Sample)
    ============================================================================
    */
    
    ANNOVAR(
        ch_annotate_input,
        annovar_db
    )
    
    ch_versions = ch_versions.mix(ANNOVAR.out.versions.first())
    
    // Log completion
    ANNOVAR.out.txt
        .map { meta, txt_file ->
            log.info "ANNOVAR annotation complete: ${meta.sample_id} -> ${txt_file}"
            return [meta, txt_file]
        }
    
    emit:
    txt      = ANNOVAR.out.txt      // Channel: [meta, multianno.txt]
    vcf      = ANNOVAR.out.vcf      // Channel: [meta, multianno.vcf]
    csv      = ANNOVAR.out.csv      // Channel: [meta, multianno.csv]
    avinput  = ANNOVAR.out.avinput  // Channel: [meta, avinput]
    versions = ch_versions          // Channel: versions.yml
}
