/*
    BAM PROCESSING SUBWORKFLOW
    MarkDuplicates -> BQSR -> Coverage
    
    Outputs BOTH:
    - MarkDup BAM (for DeepVariant)
    - BQSR BAM (for GATK HaplotypeCaller)
*/

include { GATK_MARKDUPLICATES   } from '../../modules/nf-core/gatk4/markduplicates/main'
include { GATK_BASERECALIBRATOR } from '../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK_APPLYBQSR        } from '../../modules/nf-core/gatk4/applybqsr/main'
include { MOSDEPTH              } from '../../modules/nf-core/mosdepth/main'

workflow BAM_PROCESSING {
    
    take:
    bam
    fasta
    fasta_fai
    fasta_dict
    dbsnp
    dbsnp_tbi
    known_indels
    known_indels_tbi
    intervals
    
    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    // MarkDuplicates
    GATK_MARKDUPLICATES(bam)
    ch_versions = ch_versions.mix(GATK_MARKDUPLICATES.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(GATK_MARKDUPLICATES.out.metrics.collect { it[1] })
    
    // BQSR - BaseRecalibrator
    GATK_BASERECALIBRATOR(
        GATK_MARKDUPLICATES.out.bam,
        fasta,
        fasta_fai,
        fasta_dict,
        dbsnp,
        dbsnp_tbi,
        known_indels,
        known_indels_tbi,
        intervals
    )
    ch_versions = ch_versions.mix(GATK_BASERECALIBRATOR.out.versions.first())
    
    // BQSR - ApplyBQSR
    ch_bam_table = GATK_MARKDUPLICATES.out.bam.join(GATK_BASERECALIBRATOR.out.table)
    
    GATK_APPLYBQSR(
        ch_bam_table,
        fasta,
        fasta_fai,
        fasta_dict,
        intervals
    )
    ch_versions = ch_versions.mix(GATK_APPLYBQSR.out.versions.first())
    
    // Coverage calculation
    MOSDEPTH(
        GATK_APPLYBQSR.out.bam,
        intervals
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary.collect { it[1] })
    
    emit:
    bam_markdup   = GATK_MARKDUPLICATES.out.bam  // For DeepVariant
    bam_bqsr      = GATK_APPLYBQSR.out.bam       // For GATK
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
