/*
========================================================================================
    PEDDY MODULE
========================================================================================
    Pedigree and ancestry QC
    Validates pedigree relationships using genetic data
    
    Input:
        vcf         - VCF file (with variants)
        tbi         - VCF index
        ped         - PED file with family relationships
    
    Output:
        results     - Peddy output files (html, csv)
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

process PEDDY {
    tag "peddy"
    label 'process_medium'
    container "${params.containers.peddy}"
    
    publishDir "${params.outdir}/qc/peddy", mode: params.publish_dir_mode, failOnError: false
    
    input:
    path vcf
    path tbi
    path ped
    
    output:
    path "peddy.*", emit: results
    path "versions.yml", emit: versions
    
    script:
    """
    echo "=============================================="
    echo "Peddy - Pedigree QC"
    echo "=============================================="
    echo "VCF: ${vcf}"
    echo "PED: ${ped}"
    
    peddy --plot --prefix peddy ${vcf} ${ped}
    
    echo "Peddy analysis complete"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peddy: \$(peddy --version 2>&1 | sed 's/peddy, version //')
    END_VERSIONS
    """
}
