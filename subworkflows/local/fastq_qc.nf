/*
    FASTQ QC SUBWORKFLOW
    FastQC and Fastp trimming
*/

include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { FASTP  } from '../../modules/nf-core/fastp/main'

workflow FASTQ_QC {
    
    take:
    reads
    skip_fastqc
    skip_fastp
    
    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    // Transform reads channel: [meta, r1, r2] -> [meta, [r1, r2]]
    ch_reads_tuple = reads.map { meta, r1, r2 -> [ meta, [ r1, r2 ] ] }
    
    // FastQC
    if (!skip_fastqc) {
        FASTQC(ch_reads_tuple)
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    }
    
    // Fastp trimming
    if (!skip_fastp) {
        FASTP(ch_reads_tuple)
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] })
        ch_reads_out = FASTP.out.reads
    } else {
        ch_reads_out = ch_reads_tuple
    }
    
    emit:
    reads         = ch_reads_out
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
