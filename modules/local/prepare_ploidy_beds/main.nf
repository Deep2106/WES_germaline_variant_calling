/*
========================================================================================
    PREPARE PLOIDY BEDS
========================================================================================
    Creates sex-specific interval BED files from the original WES targets BED:
    
    1. autosomal_par.bed     - Autosomes + PAR regions (diploid for males)
    2. sexchrom_nonpar.bed   - Non-PAR chrX + chrY (haploid for males)
    3. no_chrY.bed           - Everything except chrY (for females)
    
    GRCh38 PAR regions:
      PAR1: chrX:10001-2781479
      PAR2: chrX:155701383-156030895
    
    Run ONCE per pipeline invocation (not per sample).
----------------------------------------------------------------------------------------
*/

process PREPARE_PLOIDY_BEDS {
    tag "prepare_ploidy_beds"
    label 'process_low'

    container "${params.containers.bcftools}"

    input:
    path targets_bed
    path par_bed

    output:
    path "autosomal_par.bed"    , emit: autosomal_par
    path "sexchrom_nonpar.bed"  , emit: sexchrom_nonpar
    path "no_chrY.bed"          , emit: no_chrY
    path "versions.yml"         , emit: versions

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "PREPARING PLOIDY-AWARE BED FILES"
    echo "=============================================="
    echo "Input targets: ${targets_bed}"
    echo "PAR regions:   ${par_bed}"

    # =========================================================================
    # 1. autosomal_par.bed: Autosomes + PAR regions only
    #    Used for males ploidy=2 run
    # =========================================================================

    # Extract autosomal intervals (chr1-chr22)
    grep -E "^chr[0-9]+\\b" ${targets_bed} > autosomes.bed || true

    # Extract PAR regions from targets (intersect targets with PAR BED)
    # First sort both files
    sort -k1,1 -k2,2n ${targets_bed} > targets_sorted.bed
    sort -k1,1 -k2,2n ${par_bed} > par_sorted.bed

    # Get target intervals that overlap PAR regions
    awk 'NR==FNR{chr[\$1]=1; start[\$1]=\$2; end[\$1]=\$3; next}
         \$1 in chr && \$2 < end[\$1] && \$3 > start[\$1] {
            # Clip to PAR boundaries
            s = (\$2 > start[\$1]) ? \$2 : start[\$1]
            e = (\$3 < end[\$1]) ? \$3 : end[\$1]
            if (s < e) print \$1"\\t"s"\\t"e
         }' par_sorted.bed targets_sorted.bed > par_targets.bed || true

    # Combine autosomes + PAR targets
    cat autosomes.bed par_targets.bed | sort -k1,1 -k2,2n > autosomal_par.bed

    # =========================================================================
    # 2. sexchrom_nonpar.bed: Non-PAR chrX + all chrY
    #    Used for males ploidy=1 run
    # =========================================================================

    # Extract chrX targets
    grep -E "^chrX\\b" ${targets_bed} > chrX_all.bed || true

    # Remove PAR regions from chrX (subtract PAR from chrX targets)
    # Simple approach: keep chrX intervals that DON'T overlap PAR
    if [ -s chrX_all.bed ] && [ -s par_sorted.bed ]; then
        awk 'NR==FNR{
                n++; chr[n]=\$1; start[n]=\$2; end[n]=\$3; next
             }
             {
                is_par=0
                for(i=1;i<=n;i++){
                    if(\$1==chr[i] && \$2 < end[i] && \$3 > start[i]){
                        is_par=1; break
                    }
                }
                if(!is_par) print
             }' par_sorted.bed chrX_all.bed > chrX_nonpar.bed
    else
        cp chrX_all.bed chrX_nonpar.bed 2>/dev/null || touch chrX_nonpar.bed
    fi

    # Extract chrY targets
    grep -E "^chrY\\b" ${targets_bed} > chrY.bed || true

    # Combine non-PAR chrX + chrY
    cat chrX_nonpar.bed chrY.bed | sort -k1,1 -k2,2n > sexchrom_nonpar.bed

    # =========================================================================
    # 3. no_chrY.bed: Everything except chrY
    #    Used for females (diploid everywhere, skip chrY)
    # =========================================================================
    grep -v "^chrY\\b" ${targets_bed} > no_chrY.bed || cp ${targets_bed} no_chrY.bed

    # =========================================================================
    # Summary
    # =========================================================================
    echo ""
    echo "Results:"
    echo "  autosomal_par.bed:    \$(wc -l < autosomal_par.bed) intervals (autosomes + PAR)"
    echo "  sexchrom_nonpar.bed:  \$(wc -l < sexchrom_nonpar.bed) intervals (non-PAR chrX + chrY)"
    echo "  no_chrY.bed:          \$(wc -l < no_chrY.bed) intervals (all except chrY)"
    echo ""

    # Sanity check
    TOTAL_ORIG=\$(wc -l < ${targets_bed})
    TOTAL_MALE=\$(( \$(wc -l < autosomal_par.bed) + \$(wc -l < sexchrom_nonpar.bed) ))
    echo "Original intervals:     \${TOTAL_ORIG}"
    echo "Male combined intervals: \${TOTAL_MALE} (autosomal_par + sexchrom_nonpar)"
    echo "Female intervals:        \$(wc -l < no_chrY.bed)"

    # Cleanup
    rm -f autosomes.bed targets_sorted.bed par_sorted.bed par_targets.bed
    rm -f chrX_all.bed chrX_nonpar.bed chrY.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
