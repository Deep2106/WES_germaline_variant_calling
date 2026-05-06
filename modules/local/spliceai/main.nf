/*
========================================================================================
    SPLICEAI - PER CHROMOSOME
========================================================================================
    Runs SpliceAI on a single chromosome VCF.
    Designed for scatter-gather parallelism: 24 chroms × 1 CPU each.
    
    SpliceAI is single-threaded — parallelism comes from running
    multiple chromosomes simultaneously via Nextflow scheduling.
    
    Outputs UNCOMPRESSED VCF — compression handled by MERGE_SPLICEAI_VCFS
    (bcftools container) to avoid needing htslib in spliceai container.
----------------------------------------------------------------------------------------
*/

process SPLICEAI {
    tag "spliceai_${vcf.baseName.replaceAll('.vcf', '')}"
    label 'process_spliceai'

    container "${params.containers.spliceai}"

    // SpliceAI can be slow on large chromosomes — allow retries
    errorStrategy 'retry'
    maxRetries 1

    input:
    path vcf
    path tbi
    path fasta
    path fasta_fai

    output:
    path "*.spliceai.vcf"          , emit: vcf
    path "*.spliceai_summary.txt"  , emit: summary
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def chrom      = vcf.baseName.replaceAll('.vcf', '')
    def annotation = params.spliceai_annotation ?: 'grch38'
    def distance   = params.spliceai_distance ?: 500
    def prefix     = "${chrom}.spliceai"
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "SPLICEAI ANNOTATION - ${chrom}"
    echo "=============================================="
    echo "Input VCF: ${vcf}"
    echo "Reference: ${fasta}"
    echo "Annotation: ${annotation}"
    echo "Max distance: ${distance}"
    echo "=============================================="

    # Count input variants
    INPUT_COUNT=\$(zcat ${vcf} | grep -v "^#" | wc -l)
    echo "Input variants: \${INPUT_COUNT}"

    if [ \${INPUT_COUNT} -eq 0 ]; then
        echo "WARNING: No variants in ${chrom} — copying input as-is"
        zcat ${vcf} > ${prefix}.vcf
        echo "No variants to annotate in ${chrom}" > ${chrom}.spliceai_summary.txt
    else
        echo ""
        echo "Running SpliceAI on ${chrom} (\${INPUT_COUNT} variants)..."
        START_TIME=\$(date +%s)

        # Decompress for SpliceAI (requires uncompressed input)
        zcat ${vcf} > input.vcf

        # Run SpliceAI
        spliceai \\
            -I input.vcf \\
            -O ${prefix}.vcf \\
            -R ${fasta} \\
            -A ${annotation} \\
            -D ${distance}

        END_TIME=\$(date +%s)
        ELAPSED=\$((END_TIME - START_TIME))
        echo "SpliceAI completed for ${chrom} in \${ELAPSED} seconds"

        # Output count
        OUTPUT_COUNT=\$(grep -v "^#" ${prefix}.vcf | wc -l)

        # Summary
        {
            echo "SpliceAI Summary - ${chrom}"
            echo "==========================="
            echo "Chromosome: ${chrom}"
            echo "Input variants: \${INPUT_COUNT}"
            echo "Output variants: \${OUTPUT_COUNT}"
            echo "Runtime: \${ELAPSED} seconds"
            echo "Annotation: ${annotation}"
            echo "Max distance: ${distance}"
        } > ${chrom}.spliceai_summary.txt

        cat ${chrom}.spliceai_summary.txt

        # Cleanup temp files
        rm -f input.vcf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spliceai: \$(spliceai --version 2>&1 | head -n1 || echo "unknown")
        chromosome: ${chrom}
        variants_processed: \${INPUT_COUNT}
    END_VERSIONS
    """
}
