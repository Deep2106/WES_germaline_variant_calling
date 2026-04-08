/*
========================================================================================
    COLLECT DEEPVARIANT GVCFS MODULE
========================================================================================
    Manages persistent storage of DeepVariant GVCFs for incremental joint calling
    
    - Copies current batch GVCFs to persistent storage
    - Collects ALL GVCFs (previous + current) for GLnexus
    - Similar concept to GenomicsDB but simpler (just file storage)
========================================================================================
*/

process COLLECT_DV_GVCFS {
    tag "collect_dv_gvcfs"
    label 'process_medium'
    
    container "${params.containers.bcftools}"
    
    // Publish current batch GVCFs to persistent storage
    publishDir "${params.outdir}/deepvariant_gvcfs", mode: 'copy', pattern: "*.g.vcf.gz*"
    
    input:
    path current_gvcfs       // Current batch GVCFs
    path current_tbis        // Current batch TBIs
    val  dv_gvcf_path        // Path to existing GVCFs (or empty for first batch)
    val  is_update           // True if incremental run
    
    output:
    path "all_gvcfs/*.g.vcf.gz", emit: all_gvcfs
    path "versions.yml"        , emit: versions
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=============================================="
    echo "COLLECT DEEPVARIANT GVCFS"
    echo "=============================================="
    echo "Is update mode: ${is_update}"
    echo "Existing GVCF path: ${dv_gvcf_path}"
    
    mkdir -p all_gvcfs
    
    # Count current batch
    CURRENT_COUNT=0
    for gvcf in *.g.vcf.gz; do
        if [ -f "\$gvcf" ] && [ "\$gvcf" != "*.g.vcf.gz" ]; then
            ((CURRENT_COUNT++)) || true
        fi
    done
    echo "Current batch GVCFs: \${CURRENT_COUNT}"
    
    # Collect existing GVCFs if in update mode
    EXISTING_COUNT=0
    if [ "${is_update}" = "true" ] && [ -d "${dv_gvcf_path}" ]; then
        echo ""
        echo "Collecting existing GVCFs from ${dv_gvcf_path}..."
        for existing_gvcf in ${dv_gvcf_path}/*.g.vcf.gz; do
            if [ -f "\$existing_gvcf" ]; then
                basename_gvcf=\$(basename "\$existing_gvcf")
                # Only copy if not in current batch (avoid duplicates)
                if [ ! -f "\$basename_gvcf" ]; then
                    cp -L "\$existing_gvcf" all_gvcfs/
                    ((EXISTING_COUNT++)) || true
                    echo "  Added existing: \$basename_gvcf"
                else
                    echo "  Skip duplicate: \$basename_gvcf"
                fi
            fi
        done
        echo "Existing GVCFs collected: \${EXISTING_COUNT}"
    fi
    
    # Add current batch GVCFs
    echo ""
    echo "Adding current batch GVCFs..."
    for gvcf in *.g.vcf.gz; do
        if [ -f "\$gvcf" ] && [ "\$gvcf" != "*.g.vcf.gz" ]; then
            cp -L "\$gvcf" all_gvcfs/
            echo "  Added: \$gvcf"
        fi
    done
    
    # Summary
    TOTAL_COUNT=\$(ls all_gvcfs/*.g.vcf.gz 2>/dev/null | wc -l)
    echo ""
    echo "=============================================="
    echo "SUMMARY"
    echo "=============================================="
    echo "Existing GVCFs: \${EXISTING_COUNT}"
    echo "Current batch:  \${CURRENT_COUNT}"
    echo "Total GVCFs:    \${TOTAL_COUNT}"
    echo ""
    echo "All GVCFs for GLnexus:"
    ls -1 all_gvcfs/
    echo "=============================================="
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        existing_gvcfs: \${EXISTING_COUNT}
        current_batch: \${CURRENT_COUNT}
        total_gvcfs: \${TOTAL_COUNT}
    END_VERSIONS
    """
}
