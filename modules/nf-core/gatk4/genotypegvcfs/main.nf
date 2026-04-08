/*
========================================================================================
    GATK GENOTYPEGVCFS MODULE
========================================================================================
    Perform joint genotyping on GenomicsDB workspace
    Supports pedigree file for family-aware calling

    PERFORMANCE NOTE: This module does NOT use --intervals flag.
    When GenomicsDB is created with --merge-input-intervals true, it stores
    chromosome-wide sparse data. Running GenotypeGVCFs WITHOUT intervals
    allows sequential scanning which is ~700x faster than interval-by-interval
    seeking (11 minutes vs 17+ hours for WES with 215K intervals).

    OPTIMIZATION: Copies GenomicsDB to local scratch for faster I/O

    BATCH NAMESPACING: Output files and publishDir are prefixed with
    params.batch_name to prevent overwrite across batches.
----------------------------------------------------------------------------------------
*/

process GATK_GENOTYPEGVCFS {
    tag "joint_genotype"
    label 'process_high'
    container "${params.containers.gatk4}"

    publishDir "${params.outdir}/variants/joint/${params.batch_name}", mode: 'copy'

    input:
    path genomicsdb
    path fasta
    path fasta_fai
    path fasta_dict
    path dbsnp
    path dbsnp_tbi
    path pedigree

    output:
    tuple path("${params.batch_name}.joint.vcf.gz"), path("${params.batch_name}.joint.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def batch       = params.batch_name
    def avail_mem   = task.memory ? task.memory.toGiga() - 8 : 56
    def java_opts   = "-Xmx${avail_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Djava.io.tmpdir=\$TMPDIR -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
    def pedigree_arg = pedigree.name != 'NO_FILE' ? "--pedigree ${pedigree}" : ""
    def dbsnp_arg   = dbsnp.name != 'NO_FILE'    ? "--dbsnp ${dbsnp}"       : ""

    """
    #!/bin/bash
    set -euo pipefail

    # Set TMPDIR - use SLURM job scratch if available
    if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
        export TMPDIR="/local/scratch/\${SLURM_JOB_ID}"
    elif [ -z "\${TMPDIR:-}" ] || [ ! -d "\${TMPDIR}" ]; then
        export TMPDIR=\$(mktemp -d -p \${PWD})
        trap "rm -rf \$TMPDIR" EXIT
    fi

    echo "=============================================="
    echo "GenotypeGVCFs Configuration"
    echo "=============================================="
    echo "Batch      : ${batch}"
    echo "TMPDIR     : \${TMPDIR}"
    echo "GenomicsDB : ${genomicsdb}"
    echo "Pedigree   : ${pedigree_arg ? pedigree : 'Not provided'}"
    echo "NOTE: Running WITHOUT --intervals for optimal performance"

    # =========================================================================
    # PERFORMANCE OPTIMIZATION: Copy GenomicsDB to local scratch
    # =========================================================================
    LOCAL_GENOMICSDB="\${TMPDIR}/genomicsdb_local"

    echo "=============================================="
    echo "Copying GenomicsDB to local scratch..."
    echo "Source      : ${genomicsdb}"
    echo "Destination : \${LOCAL_GENOMICSDB}"
    echo "=============================================="

    START_COPY=\$(date +%s)
    cp -rL "${genomicsdb}" "\${LOCAL_GENOMICSDB}"
    END_COPY=\$(date +%s)
    COPY_TIME=\$((END_COPY - START_COPY))

    echo "GenomicsDB copy completed in \${COPY_TIME} seconds"
    echo "Local GenomicsDB size: \$(du -sh "\${LOCAL_GENOMICSDB}" | cut -f1)"

    # =========================================================================
    # Run GenotypeGVCFs on LOCAL copy - NO --intervals for best performance
    # =========================================================================
    echo "=============================================="
    echo "Starting GenotypeGVCFs on local scratch..."
    echo "=============================================="

    START_GENO=\$(date +%s)

    gatk --java-options "${java_opts}" GenotypeGVCFs \\
        --reference ${fasta} \\
        --variant "gendb://\${LOCAL_GENOMICSDB}" \\
        --output ${batch}.joint.vcf.gz \\
        ${pedigree_arg} \\
        ${dbsnp_arg} \\
        --create-output-variant-index true \\
        --tmp-dir \$TMPDIR \\
        ${args}

    END_GENO=\$(date +%s)
    GENO_TIME=\$((END_GENO - START_GENO))

    echo "=============================================="
    echo "GenotypeGVCFs completed in \${GENO_TIME} seconds (\$((GENO_TIME / 60)) minutes)"
    echo "=============================================="

    # Cleanup local GenomicsDB copy
    rm -rf "\${LOCAL_GENOMICSDB}"

    # Log pedigree usage
    if [ -n "${pedigree_arg}" ]; then
        echo "Pedigree file used: ${pedigree}"
        echo "Number of samples in pedigree: \$(grep -v '^#' ${pedigree} | wc -l)"
    fi

    echo "Output VCF    : ${batch}.joint.vcf.gz"
    echo "Output size   : \$(ls -lh ${batch}.joint.vcf.gz | awk '{print \$5}')"
    echo "Variant count : \$(zcat ${batch}.joint.vcf.gz | grep -v '^#' | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/^.*(GATK) v//' | sed 's/ .*\$//')
        batch: "${batch}"
        pedigree_used: ${pedigree.name != 'NO_FILE' ? 'true' : 'false'}
        intervals_used: false
        local_scratch_used: true
    END_VERSIONS
    """

    stub:
    def batch = params.batch_name
    """
    touch ${batch}.joint.vcf.gz
    touch ${batch}.joint.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.2.0
        batch: "${batch}"
    END_VERSIONS
    """
}
