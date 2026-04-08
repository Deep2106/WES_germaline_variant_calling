/*
    MOSDEPTH MODULE
    Fast coverage calculation
*/

process MOSDEPTH {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.containers.mosdepth}"
    
    publishDir "${params.outdir}/qc/coverage", mode: params.publish_dir_mode, failOnError: false
    
    input:
    tuple val(meta), path(bam), path(bai)
    path intervals
    
    output:
    tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(meta), path("*.mosdepth.summary.txt"), emit: summary
    tuple val(meta), path("*.regions.bed.gz"), emit: regions, optional: true
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def interval_arg = intervals.name != 'NO_FILE' ? "--by ${intervals}" : ""
    """
    mosdepth \\
        --threads ${task.cpus} \\
        ${interval_arg} \\
        ${prefix} \\
        ${bam}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
}
