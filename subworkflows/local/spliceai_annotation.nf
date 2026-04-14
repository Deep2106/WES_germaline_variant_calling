/*
========================================================================================
    SPLICEAI ANNOTATION SUBWORKFLOW (Scatter-Gather)
========================================================================================
    Parallelizes SpliceAI across chromosomes:
    
    1. SPLIT:  filtered.vcf.gz → chr1.vcf.gz, chr2.vcf.gz, ..., chrY.vcf.gz
    2. SCORE:  SpliceAI per chromosome (1 CPU, 16GB each - ~20 run in parallel)
               Outputs UNCOMPRESSED VCF (no bgzip/tabix in spliceai container)
    3. MERGE:  bgzip + tabix + bcftools concat → filtered.spliceai.merged.vcf.gz
    
    This achieves ~20x speedup over running SpliceAI on the full VCF since
    SpliceAI is single-threaded and cannot be parallelized internally.
----------------------------------------------------------------------------------------
*/

include { SPLIT_VCF_BY_CHROM } from '../../modules/local/split_vcf_by_chrom/main'
include { SPLICEAI            } from '../../modules/local/spliceai/main'
include { MERGE_SPLICEAI_VCFS } from '../../modules/local/merge_spliceai_vcfs/main'

workflow SPLICEAI_ANNOTATION {

    take:
    ch_vcf          // path: filtered VCF
    ch_tbi          // path: filtered VCF index
    ch_fasta        // path: reference FASTA
    ch_fasta_fai    // path: reference FASTA index

    main:
    ch_versions = Channel.empty()

    /*
    ========================================================================
        STEP 1: Split VCF by chromosome
    ========================================================================
    */
    SPLIT_VCF_BY_CHROM(ch_vcf, ch_tbi)

    ch_versions = ch_versions.mix(SPLIT_VCF_BY_CHROM.out.versions)

    /*
    ========================================================================
        STEP 2: Run SpliceAI per chromosome (scatter)
        .flatten() turns [chr1.vcf.gz, chr2.vcf.gz, ...] into individual items
        Each gets scheduled as a separate SLURM job (1 cpu, 16GB)
    ========================================================================
    */
    ch_per_chrom = SPLIT_VCF_BY_CHROM.out.vcfs
        .flatten()
        .map { vcf ->
            def tbi_path = file("${vcf}.tbi")
            return [vcf, tbi_path]
        }

    SPLICEAI(
        ch_per_chrom.map { vcf, tbi -> vcf },
        ch_per_chrom.map { vcf, tbi -> tbi },
        ch_fasta,
        ch_fasta_fai
    )

    ch_versions = ch_versions.mix(SPLICEAI.out.versions.first())

    /*
    ========================================================================
        STEP 3: Merge per-chromosome results (gather)
        MERGE handles bgzip + tabix + concat (bcftools container)
    ========================================================================
    */
    MERGE_SPLICEAI_VCFS(
        SPLICEAI.out.vcf.collect(),
        SPLICEAI.out.summary.collect()
    )

    ch_versions = ch_versions.mix(MERGE_SPLICEAI_VCFS.out.versions)

    emit:
    vcf      = MERGE_SPLICEAI_VCFS.out.vcf       // path: merged SpliceAI VCF
    tbi      = MERGE_SPLICEAI_VCFS.out.tbi        // path: merged SpliceAI VCF index
    summary  = MERGE_SPLICEAI_VCFS.out.summary    // path: combined summary
    versions = ch_versions                         // channel: versions
}
