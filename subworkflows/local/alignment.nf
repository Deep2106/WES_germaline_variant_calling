/*
    ALIGNMENT SUBWORKFLOW
    BWA-MEM2 alignment
*/

include { BWAMEM2_MEM } from '../../modules/nf-core/bwamem2/mem/main'

workflow ALIGNMENT {
    
    take:
    reads
    fasta      // Path string (not file object)
    bwa_index  // Path string (not file object)
    
    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    // Alignment - pass fasta and bwa_index as values (original paths)
    BWAMEM2_MEM(reads, fasta, bwa_index)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    
    emit:
    bam           = BWAMEM2_MEM.out.bam
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
