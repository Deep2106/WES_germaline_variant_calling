#!/usr/bin/env nextflow

/*
========================================================================================
    WES/WGS GERMLINE VARIANT CALLING PIPELINE v2.0
========================================================================================
    Production pipeline with:
    - Incremental joint calling with GenomicsDB (auto-detection)
    - Automatic PED file generation and management
    - Optional DeepVariant dual-caller (default: OFF)
    - De novo mutation flagging (hiConfDeNovo/loConfDeNovo)
    - Caller concordance tagging (BOTH/GATK_ONLY/DV_ONLY)
    - Batch-specific annotation
    - SpliceAI integration (optional)
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Import for script-level variables accessible in functions
import groovy.transform.Field

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================================================
     WES/WGS Germline Variant Calling Pipeline v${workflow.manifest.version}
    =========================================================================

    Usage:
        nextflow run main.nf -profile slurm --input samplesheet.csv [options]

    Required:
        --input             Path to sample CSV file
        --fasta             Reference genome FASTA
        --fasta_fai         Reference genome FASTA index
        --fasta_dict        Reference sequence dictionary
        --targets_bed       Target regions BED file (for WES)
        --annovar_db        ANNOVAR database directory
        --dbsnp             dbSNP VCF file
        --dbsnp_tbi         dbSNP VCF index
        --known_indels      Known indels VCF
        --known_indels_tbi  Known indels VCF index

    Incremental Joint Calling:
        --genomicsdb_path   Path for GenomicsDB (auto-detects create/update)
        --dv_gvcf_path      Path for DeepVariant GVCFs storage
        --master_ped        Path to existing master PED file

    Optional:
        --run_deepvariant   Run DeepVariant alongside GATK (default: false)
        --run_spliceai      Run SpliceAI annotation (default: false)
        --batch_name        Batch identifier for MultiQC reports
        --outdir            Output directory (default: ./results)

    Sample CSV format (standard PED order):
        family_id,sample_id,paternal_id,maternal_id,sex,phenotype,fastq_r1,fastq_r2

    Example:
        # First batch - creates GenomicsDB and master.ped
        nextflow run main.nf -profile slurm --input batch1.csv \\
            --genomicsdb_path /path/to/genomicsdb --outdir ./results

        # Second batch - auto-detects update mode
        nextflow run main.nf -profile slurm --input batch2.csv \\
            --genomicsdb_path /path/to/genomicsdb \\
            --master_ped ./results/ped/master.ped --outdir ./results
    =========================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    AUTO-DETECTION: Initial vs Incremental Mode
========================================================================================
*/

// Use @Field for script-level variables accessible everywhere
@Field def dbExists = false
@Field def pedExists = false
@Field def isUpdateMode = false
@Field def modeReason = ""

// Initialize mode detection
dbExists = params.genomicsdb_path ? file(params.genomicsdb_path).isDirectory() : false
pedExists = params.master_ped ? file(params.master_ped).exists() : false

if (params.genomicsdb_update) {
    isUpdateMode = true
    modeReason = "Forced by --genomicsdb_update"
} else if (dbExists && pedExists) {
    isUpdateMode = true
    modeReason = "Auto-detected (GenomicsDB + master.ped exist)"
} else if (dbExists && !pedExists) {
    log.warn "GenomicsDB exists but master.ped not found - CREATE mode"
    modeReason = "GenomicsDB exists but no master.ped"
} else {
    modeReason = "First run (no existing GenomicsDB)"
}

/*
========================================================================================
    PRINT BANNER
========================================================================================
*/

log.info """
╔═══════════════════════════════════════════════════════════════════════════════╗
║       WES/WGS GERMLINE VARIANT CALLING PIPELINE v${workflow.manifest.version}                    ║
╠═══════════════════════════════════════════════════════════════════════════════╣
║  Input samplesheet  : ${params.input ?: 'NOT SET'}
║  Batch name         : ${params.batch_name}
║  Reference genome   : ${params.fasta ?: 'NOT SET'}
║  Target regions     : ${params.targets_bed ?: 'WGS mode'}
║  GenomicsDB path    : ${params.genomicsdb_path ?: 'Will create'}
║  DV GVCF path       : ${params.dv_gvcf_path ?: 'Will create'}
║  Mode               : ${isUpdateMode ? 'UPDATE (incremental)' : 'CREATE (initial)'}
║  Mode reason        : ${modeReason}
║  Master PED         : ${params.master_ped ?: 'Will create'}
║  Run DeepVariant    : ${params.run_deepvariant}
║  Run SpliceAI       : ${params.run_spliceai}
║  Output directory   : ${params.outdir}
╚═══════════════════════════════════════════════════════════════════════════════╝
""".stripIndent()

/*
========================================================================================
    VALIDATE PARAMETERS
========================================================================================
*/

def validateParameters() {
    def errors = []

    if (!params.input) errors << "ERROR: --input required"
    else if (!file(params.input).exists()) errors << "ERROR: Input not found: ${params.input}"

    if (!params.fasta) errors << "ERROR: --fasta required"
    if (!params.fasta_fai) errors << "ERROR: --fasta_fai required"
    if (!params.fasta_dict) errors << "ERROR: --fasta_dict required"
    if (!params.dbsnp) errors << "ERROR: --dbsnp required"
    if (!params.known_indels) errors << "ERROR: --known_indels required"
    if (!params.annovar_db) errors << "ERROR: --annovar_db required"

    if (isUpdateMode && !dbExists) {
        errors << "ERROR: Update mode but GenomicsDB not found"
    }

    if (errors) {
        log.error "\nPARAMETER VALIDATION FAILED"
        errors.each { log.error it }
        exit 1
    }

    log.info "Parameters validated"
}

validateParameters()

/*
========================================================================================
    INCLUDE SUBWORKFLOWS
========================================================================================
*/

include { INPUT_VALIDATION } from './subworkflows/local/input_validation'
include { FASTQ_QC         } from './subworkflows/local/fastq_qc'
include { ALIGNMENT        } from './subworkflows/local/alignment'
include { BAM_PROCESSING   } from './subworkflows/local/bam_processing'
include { VARIANT_CALLING  } from './subworkflows/local/variant_calling'
include { JOINT_CALLING    } from './subworkflows/local/joint_calling'
include { FAMILY_ANALYSIS  } from './subworkflows/local/family_analysis'
include { VARIANT_FILTER   } from './subworkflows/local/variant_filter'
include { ANNOTATION       } from './subworkflows/local/annotation'
include { QC_REPORTING     } from './subworkflows/local/qc_reporting'

// SpliceAI module (optional - uncomment when ready to use)
include { SPLICEAI_ANNOTATION } from './subworkflows/local/spliceai_annotation'


/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Ensure assets directory exists
    file("${projectDir}/assets").mkdirs()
    if (!file("${projectDir}/assets/NO_FILE").exists()) {
        file("${projectDir}/assets/NO_FILE").text = ""
    }

    /*
    ============================================================================
        STAGE 1: INPUT VALIDATION AND PED GENERATION
    ============================================================================
    */

    ch_master_ped_input = params.master_ped ?
        Channel.value(file(params.master_ped)) :
        Channel.value(file("${projectDir}/assets/NO_FILE"))

    INPUT_VALIDATION(
        file(params.input),
        ch_master_ped_input
    )

    ch_reads = INPUT_VALIDATION.out.reads
    ch_batch_ped = INPUT_VALIDATION.out.batch_ped
    ch_master_ped = INPUT_VALIDATION.out.master_ped
    ch_versions = ch_versions.mix(INPUT_VALIDATION.out.versions)

    /*
    ============================================================================
        STAGE 2: FASTQ QC AND TRIMMING
    ============================================================================
    */

    FASTQ_QC(
        ch_reads,
        params.skip_fastqc,
        params.skip_fastp
    )

    ch_reads_processed = FASTQ_QC.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC.out.multiqc_files)
    ch_versions = ch_versions.mix(FASTQ_QC.out.versions)

    /*
    ============================================================================
        STAGE 3: ALIGNMENT
    ============================================================================
    */

    // File channels for processes that need staged files
    ch_fasta = Channel.value(file(params.fasta))
    ch_fasta_fai = Channel.value(file(params.fasta_fai))
    ch_fasta_dict = Channel.value(file(params.fasta_dict))
    ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : Channel.value(file("${projectDir}/assets/NO_FILE"))

    // String path channels for BWA-MEM2 (uses original paths via Singularity bind)
    ch_fasta_path = Channel.value(params.fasta)
    ch_bwa_index_path = params.bwa_index ? Channel.value(params.bwa_index) : Channel.value("NO_FILE")

    ALIGNMENT(
        ch_reads_processed,
        ch_fasta_path,
        ch_bwa_index_path
    )

    ch_bam_sorted = ALIGNMENT.out.bam
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.multiqc_files)
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    /*
    ============================================================================
        STAGE 4: BAM PROCESSING
        Outputs BOTH MarkDup BAM (for DeepVariant) and BQSR BAM (for GATK)
    ============================================================================
    */

    ch_dbsnp = Channel.value(file(params.dbsnp))
    ch_dbsnp_tbi = Channel.value(file(params.dbsnp_tbi))
    ch_known_indels = Channel.value(file(params.known_indels))
    ch_known_indels_tbi = Channel.value(file(params.known_indels_tbi))
    ch_targets = params.targets_bed ? Channel.value(file(params.targets_bed)) : Channel.value(file("${projectDir}/assets/NO_FILE"))

    BAM_PROCESSING(
        ch_bam_sorted,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_known_indels,
        ch_known_indels_tbi,
        ch_targets
    )

    ch_bam_markdup = BAM_PROCESSING.out.bam_markdup
    ch_bam_bqsr = BAM_PROCESSING.out.bam_bqsr
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PROCESSING.out.multiqc_files)
    ch_versions = ch_versions.mix(BAM_PROCESSING.out.versions)

    /*
    ============================================================================
        STAGE 5: VARIANT CALLING
        GATK HaplotypeCaller (BQSR BAM) + optional DeepVariant (MarkDup BAM)
    ============================================================================
    */

    VARIANT_CALLING(
        ch_bam_bqsr,
        ch_bam_markdup,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_targets,
        params.run_deepvariant
    )

    ch_gvcf = VARIANT_CALLING.out.gvcf
    ch_dv_gvcf = VARIANT_CALLING.out.deepvariant_gvcf
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

    /*
    ============================================================================
        STAGE 6: JOINT CALLING (Incremental GenomicsDB)
    ============================================================================
    */

    JOINT_CALLING(
        ch_gvcf,
        ch_dv_gvcf,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_targets,
        ch_master_ped,
        params.genomicsdb_path ?: "${params.outdir}/genomicsdb/genomicsdb",
        params.dv_gvcf_path ?: "${params.outdir}/variants/deepvariant",
        isUpdateMode,
        params.run_deepvariant
    )

    ch_joint_vcf = JOINT_CALLING.out.vcf
    ch_joint_tbi = JOINT_CALLING.out.tbi
    ch_versions = ch_versions.mix(JOINT_CALLING.out.versions)

    /*
    ============================================================================
        STAGE 7: FAMILY-AWARE ANALYSIS
    ============================================================================
    */

    FAMILY_ANALYSIS(
        ch_joint_vcf,
        ch_joint_tbi,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_master_ped
    )

    ch_family_vcf = FAMILY_ANALYSIS.out.vcf
    ch_family_tbi = FAMILY_ANALYSIS.out.tbi
    ch_multiqc_files = ch_multiqc_files.mix(FAMILY_ANALYSIS.out.multiqc_files)
    ch_versions = ch_versions.mix(FAMILY_ANALYSIS.out.versions)

    /*
    ============================================================================
        STAGE 8: VARIANT FILTERING
    ============================================================================
    */

    VARIANT_FILTER(
        ch_family_vcf,
        ch_family_tbi,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict
    )

    ch_filtered_vcf = VARIANT_FILTER.out.vcf
    ch_filtered_tbi = VARIANT_FILTER.out.tbi
    ch_versions = ch_versions.mix(VARIANT_FILTER.out.versions)

    /*
    ============================================================================
        STAGE 8.5: SPLICEAI ANNOTATION (OPTIONAL)
        Runs on filtered VCF before per-sample annotation
        Scores unique variants once - scores flow to per-sample VCFs
    ============================================================================
    */

    // Uncomment when SpliceAI module is ready:
    
  if (params.run_spliceai) {
    SPLICEAI_ANNOTATION(
        ch_filtered_vcf,
        ch_filtered_tbi,
        ch_fasta,
        Channel.value(file(params.fasta_fai))
    )
    ch_vcf_for_annotation = SPLICEAI_ANNOTATION.out.vcf
    ch_tbi_for_annotation = SPLICEAI_ANNOTATION.out.tbi
    ch_versions = ch_versions.mix(SPLICEAI_ANNOTATION.out.versions)
} else {
    ch_vcf_for_annotation = ch_filtered_vcf
    ch_tbi_for_annotation = ch_filtered_tbi
} 


    /*
    ============================================================================
        STAGE 9: ANNOTATION (Per-sample, current batch only)
    ============================================================================
    */

    ch_sample_meta = ch_reads.map { meta, r1, r2 -> meta }

    ANNOTATION(
        ch_vcf_for_annotation,
        ch_tbi_for_annotation,
        ch_sample_meta,
        file(params.annovar_db)
    )

    ch_versions = ch_versions.mix(ANNOTATION.out.versions)

    /*
    ============================================================================
        STAGE 10: QC REPORTING (Batch-specific)
    ============================================================================
    */

    ch_versions_filtered = ch_versions.flatten().filter { it != null }

    QC_REPORTING(
        ch_multiqc_files.collect().ifEmpty([]),
        ch_versions_filtered.unique().collectFile(name: 'versions.yml').ifEmpty(file("${projectDir}/assets/NO_FILE"))
    )

    ch_versions_filtered.unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}

/*
========================================================================================
    COMPLETION HANDLERS
========================================================================================
*/

workflow.onComplete {
    log.info """
╔═══════════════════════════════════════════════════════════════════════════════╗
║                         PIPELINE COMPLETED                                    ║
╠═══════════════════════════════════════════════════════════════════════════════╣
║  Status      : ${workflow.success ? 'SUCCESS' : 'FAILED'}
║  Duration    : ${workflow.duration}
║  Batch       : ${params.batch_name}
║  Mode        : ${isUpdateMode ? 'UPDATE' : 'CREATE'}
║  Output      : ${params.outdir}
╚═══════════════════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow.onError {
    log.error """
╔═══════════════════════════════════════════════════════════════════════════════╗
║                           PIPELINE FAILED                                     ║
║  Error: ${workflow.errorMessage}
╚═══════════════════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
