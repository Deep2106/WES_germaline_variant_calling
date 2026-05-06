/*
    ALIGNMENT SUBWORKFLOW
    Lane-aware BWA-MEM2 alignment → per-sample BAM merge

    Flow:
      1. Align each lane independently (meta includes lane field)
      2. Group per-lane BAMs by sample_id (strip lane from meta)
      3. Single-lane samples → pass through directly
      4. Multi-lane samples  → SAMTOOLS_MERGE → merged BAM
      5. Emit unified BAM channel (same meta structure for all downstream)
*/

include { BWAMEM2_MEM    } from '../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_MERGE } from '../../modules/local/samtools_merge'

workflow ALIGNMENT {

    take:
    reads      // tuple val(meta), path(reads)  — meta MUST contain 'lane'
    fasta      // val (path string, not file object)
    bwa_index  // val (path string, not file object)

    main:

    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    // -------------------------------------------------------------------
    // 1. Align each lane (meta.id + meta.lane are both set here)
    // -------------------------------------------------------------------
    BWAMEM2_MEM(reads, fasta, bwa_index)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    // -------------------------------------------------------------------
    // 2. Strip lane from meta and group by sample
    //    Keys kept: id, family_id, sex, phenotype (add others as needed)
    // -------------------------------------------------------------------
    ch_grouped = BWAMEM2_MEM.out.bam
        .map { meta, bam ->
            def sample_meta = [
                id:         meta.id,
                family_id:  meta.family_id,
                sex:        meta.sex,
                phenotype:  meta.phenotype
            ]
            [ sample_meta, bam ]
        }
        .groupTuple()           // groups all lane BAMs under the same sample_meta
        .branch {
            single: it[1].size() == 1
            multi:  it[1].size() >  1
        }

    // -------------------------------------------------------------------
    // 3. Single-lane: unwrap the list so downstream sees path, not list
    // -------------------------------------------------------------------
    ch_single = ch_grouped.single.map { meta, bams -> [ meta, bams[0] ] }

    // -------------------------------------------------------------------
    // 4. Multi-lane: merge → indexed merged BAM
    // -------------------------------------------------------------------
    SAMTOOLS_MERGE(ch_grouped.multi)
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    // -------------------------------------------------------------------
    // 5. Unified BAM channel
    // -------------------------------------------------------------------
    ch_bam = ch_single.mix(SAMTOOLS_MERGE.out.bam)

    emit:
    bam           = ch_bam
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
