/*
========================================================================================
    MERGE PED MODULE
========================================================================================
    Merge new sample PED entries with existing master PED file
    Used for incremental joint calling to maintain a complete pedigree
    
    Input:
        new_ped         - PED file with new samples from current batch
        existing_ped    - Existing master PED file (optional, may be NO_FILE)
        output_path     - Path to save the merged master PED file
    
    Output:
        merged_ped      - Combined PED file with all samples
        versions        - Software versions
    
    PED file format (6 columns, tab-separated):
        FamilyID  IndividualID  PaternalID  MaternalID  Sex  Phenotype
        
    Notes:
        - Removes duplicate entries (keeps first occurrence)
        - Preserves comments (lines starting with #)
        - Sex: 1=male, 2=female, 0=unknown
        - Phenotype: 1=unaffected, 2=affected, 0=unknown
----------------------------------------------------------------------------------------
*/

process MERGE_PED {
    tag "merge_ped"
    label 'process_single'
    
    // Use python container for scripting
    container "${params.containers?.python ?: 'python:3.9-slim'}"
    
    publishDir "${params.outdir}/pedigree", mode: params.publish_dir_mode, failOnError: false ?: 'copy'
    
    input:
    path new_ped
    path existing_ped
    val output_path
    
    output:
    path "master.ped"   , emit: merged_ped
    path "versions.yml" , emit: versions
    
    script:
    def has_existing = existing_ped.name != 'NO_FILE' && existing_ped.name != 'NO_PED'
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=============================================="
    echo "MERGE PED - Pedigree File Management"
    echo "=============================================="
    echo "New PED file: ${new_ped}"
    echo "Existing PED: ${has_existing ? existing_ped : 'None (first run)'}"
    echo "Output path: ${output_path}"
    
    # Initialize merged file
    > master.ped
    
    # Track seen samples to avoid duplicates
    declare -A seen_samples
    
    # Function to process PED file
    process_ped() {
        local ped_file="\$1"
        local source="\$2"
        
        if [ ! -f "\$ped_file" ]; then
            echo "WARNING: PED file not found: \$ped_file"
            return
        fi
        
        local added=0
        local skipped=0
        
        while IFS=\$'\\t' read -r fam_id ind_id pat_id mat_id sex pheno rest || [ -n "\$fam_id" ]; do
            # Skip empty lines
            [ -z "\$fam_id" ] && continue
            
            # Preserve comments
            if [[ "\$fam_id" == "#"* ]]; then
                echo "\$fam_id\$'\\t'\$ind_id\$'\\t'\$pat_id\$'\\t'\$mat_id\$'\\t'\$sex\$'\\t'\$pheno" >> master.ped
                continue
            fi
            
            # Create unique key for sample
            local key="\${fam_id}_\${ind_id}"
            
            # Skip if already seen
            if [ -n "\${seen_samples[\$key]:-}" ]; then
                ((skipped++))
                continue
            fi
            
            # Mark as seen and add to output
            seen_samples[\$key]=1
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "\$fam_id" "\$ind_id" "\$pat_id" "\$mat_id" "\$sex" "\$pheno" >> master.ped
            ((added++))
            
        done < "\$ped_file"
        
        echo "  \$source: Added \$added samples, skipped \$skipped duplicates"
    }
    
    # Process existing PED first (if present)
    if [ "${has_existing}" = "true" ]; then
        echo ""
        echo "Processing existing master PED..."
        process_ped "${existing_ped}" "Existing"
    fi
    
    # Process new PED entries
    echo ""
    echo "Processing new PED entries..."
    process_ped "${new_ped}" "New"
    
    # Summary
    echo ""
    echo "=============================================="
    echo "MERGE SUMMARY"
    echo "=============================================="
    total_samples=\$(grep -v '^#' master.ped | grep -v '^\$' | wc -l)
    total_families=\$(grep -v '^#' master.ped | grep -v '^\$' | cut -f1 | sort -u | wc -l)
    echo "Total samples: \$total_samples"
    echo "Total families: \$total_families"
    
    # Show family breakdown
    echo ""
    echo "Samples per family:"
    grep -v '^#' master.ped | grep -v '^\$' | cut -f1 | sort | uniq -c | sort -rn | head -20
    
    # Validate PED structure
    echo ""
    echo "Validating PED structure..."
    invalid_lines=0
    line_num=0
    while IFS=\$'\\t' read -r fam_id ind_id pat_id mat_id sex pheno rest || [ -n "\$fam_id" ]; do
        ((line_num++))
        [ -z "\$fam_id" ] && continue
        [[ "\$fam_id" == "#"* ]] && continue
        
        # Check for minimum columns
        if [ -z "\$sex" ] || [ -z "\$pheno" ]; then
            echo "  WARNING: Line \$line_num has missing columns"
            ((invalid_lines++))
        fi
        
        # Check sex values
        if [[ ! "\$sex" =~ ^[012]\$ ]]; then
            echo "  WARNING: Line \$line_num has invalid sex value: \$sex"
            ((invalid_lines++))
        fi
        
        # Check phenotype values  
        if [[ ! "\$pheno" =~ ^[012]\$ ]]; then
            echo "  WARNING: Line \$line_num has invalid phenotype value: \$pheno"
            ((invalid_lines++))
        fi
        
    done < master.ped
    
    if [ \$invalid_lines -eq 0 ]; then
        echo "PED validation: PASSED"
    else
        echo "PED validation: \$invalid_lines warnings"
    fi
    
    # Copy to persistent location if specified
    if [ -n "${output_path}" ] && [ "${output_path}" != "null" ]; then
        echo ""
        echo "Copying master PED to persistent location: ${output_path}"
        mkdir -p \$(dirname ${output_path})
        cp master.ped ${output_path}
        echo "Master PED saved to: ${output_path}"
    fi
    
    echo ""
    echo "Output file: master.ped"
    ls -lh master.ped
    
    # Show first few entries
    echo ""
    echo "First 5 entries:"
    head -5 master.ped
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/.*version //' | sed 's/ .*//')
        total_samples: \$total_samples
        total_families: \$total_families
    END_VERSIONS
    """
    
    stub:
    """
    echo -e "FAM001\\tSAMPLE001\\t0\\t0\\t1\\t2" > master.ped
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: stub
    END_VERSIONS
    """
}
