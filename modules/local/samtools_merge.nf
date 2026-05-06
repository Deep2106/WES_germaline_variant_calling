/*
========================================================================================
    SAMTOOLS MERGE MODULE
========================================================================================
    Merges per-lane BAMs into a single per-sample BAM.

    Input:
        tuple val(meta), path(bams)  - meta.lane is absent/unused at this point;
                                       meta.id is the sample identifier
    Output:
        tuple val(meta), path("*.merged.bam")  - indexed merged BAM
        path "versions.yml"

    Notes:
        - Called only when a sample has > 1 lane (alignment.nf handles branching)
        - Output published to params.merged_bam_path for persistent storage
----------------------------------------------------------------------------------------
*/

process SAMTOOLS_MERGE {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.containers.samtools}"

    publishDir "${params.merged_bam_path ?: "${params.outdir}/bam/merged"}", mode: params.publish_dir_mode ?: 'copy'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.merged.bam"),     emit: bam
    tuple val(meta), path("*.merged.bam.bai"), emit: bai
    path "versions.yml",                        emit: versions

    script:
    def prefix = "${meta.id}"
    """
    samtools merge \\
        -@ ${task.cpus} \\
        -f \\
        ${prefix}.merged.bam \\
        ${bams.sort().join(' ')}

    samtools index \\
        -@ ${task.cpus} \\
        ${prefix}.merged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.merged.bam ${meta.id}.merged.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: stub
    END_VERSIONS
    """
}
