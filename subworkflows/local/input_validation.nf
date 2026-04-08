/*
    INPUT VALIDATION SUBWORKFLOW
    Validates CSV, generates PED files, creates sample channels
    
    CSV format: sample_id,family_id,paternal_id,maternal_id,sex,phenotype,fastq_r1,fastq_r2
*/

include { INPUT_CHECK } from '../../modules/local/input_check/main'

workflow INPUT_VALIDATION {
    
    take:
    samplesheet
    master_ped
    
    main:
    
    ch_versions = Channel.empty()
    
    // Run input validation
    INPUT_CHECK(samplesheet, master_ped)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    
    // Parse validated CSV into channel
    // Uses fastq_r1 and fastq_r2 column names
    ch_reads = INPUT_CHECK.out.csv
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                sample_id: row.sample_id,  // For annotation modules
                family_id: row.family_id,
                paternal_id: row.paternal_id,
                maternal_id: row.maternal_id,
                sex: row.sex,
                phenotype: row.phenotype,
                single_end: false
            ]
            def fastq_r1 = file(row.fastq_r1)
            def fastq_r2 = file(row.fastq_r2)
            [ meta, fastq_r1, fastq_r2 ]
        }
    
    emit:
    reads      = ch_reads
    batch_ped  = INPUT_CHECK.out.batch_ped
    master_ped = INPUT_CHECK.out.master_ped
    versions   = ch_versions
}
