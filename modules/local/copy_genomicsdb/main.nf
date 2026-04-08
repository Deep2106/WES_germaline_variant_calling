/*
========================================================================================
    COPY GENOMICSDB MODULE
========================================================================================
    Copy GenomicsDB from work directory to persistent storage location
    This ensures the database is preserved for incremental updates
    
    Input:
        genomicsdb      - GenomicsDB workspace from work directory
        output_path     - Destination path for persistent storage
    
    Output:
        genomicsdb_path - Path to copied GenomicsDB
        versions        - Software versions
----------------------------------------------------------------------------------------
*/

process COPY_GENOMICSDB {
    tag "copy_genomicsdb"
    label 'process_medium'
    
    container "${params.containers.python}"
    
    input:
    path genomicsdb
    val output_path
    
    output:
    path "genomicsdb_path.txt", emit: path_file
    path "versions.yml"       , emit: versions
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=============================================="
    echo "COPY GENOMICSDB TO PERSISTENT STORAGE"
    echo "=============================================="
    echo "Source: ${genomicsdb}"
    echo "Destination: ${output_path}"
    
    # Ensure output directory parent exists
    mkdir -p "\$(dirname "${output_path}")"
    
    # Remove existing if present (will be replaced)
    if [ -d "${output_path}" ]; then
        echo "Removing existing GenomicsDB at ${output_path}"
        rm -rf "${output_path}"
    fi
    
    # Copy GenomicsDB to persistent location
    echo "Copying GenomicsDB..."
    cp -r "${genomicsdb}" "${output_path}"
    
    # Verify copy
    if [ -d "${output_path}" ]; then
        echo "GenomicsDB successfully copied to ${output_path}"
        echo "${output_path}" > genomicsdb_path.txt
        echo "Size: \$(du -sh "${output_path}" | cut -f1)"
        
        # Verify key files exist
        if [ -f "${output_path}/callset.json" ]; then
            SAMPLE_COUNT=\$(grep -c "row_idx" "${output_path}/callset.json" || echo "0")
            echo "Samples in GenomicsDB: \${SAMPLE_COUNT}"
        fi
    else
        echo "ERROR: Failed to copy GenomicsDB" >&2
        exit 1
    fi
    
    echo "=============================================="
    echo "Copy complete!"
    echo "=============================================="
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        note: "GenomicsDB copied to persistent storage"
        destination: "${output_path}"
    END_VERSIONS
    """
    
    stub:
    """
    echo "${output_path}" > genomicsdb_path.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        note: "stub mode"
    END_VERSIONS
    """
}
