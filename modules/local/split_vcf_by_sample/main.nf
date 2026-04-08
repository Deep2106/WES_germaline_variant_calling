/*
========================================================================================
    SPLIT VCF BY SAMPLE MODULE
========================================================================================
    Split a multi-sample VCF into individual per-sample VCFs
    Only includes PASS variants where the sample has an ALT allele
    FILTER="PASS" and FILTER="." (unfiltered) are both retained
    Variants tagged with any hard filter (FS_filter, SOR_filter etc.) are excluded
    
    Input:
        vcf         - Multi-sample VCF file
        tbi         - VCF index
        samples     - List of sample IDs to extract
    
    Output:
        vcfs        - Channel of [meta, vcf, tbi] per sample
        versions    - Software versions
----------------------------------------------------------------------------------------
*/

process SPLIT_VCF_BY_SAMPLE {
    tag "split_vcf"
    label 'process_low'
    
    container "${params.containers.bcftools}"
    
    input:
    path vcf
    path tbi
    val samples
    
    output:
    path "*.per_sample.vcf.gz"     , emit: vcfs
    path "*.per_sample.vcf.gz.tbi" , emit: tbis
    path "sample_list.txt"         , emit: sample_list
    path "versions.yml"            , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def sample_list = samples.join('\n')
    
    """
    # Create sample list file
    cat <<EOF > sample_list.txt
${sample_list}
EOF
    
    echo "Splitting VCF for ${samples.size()} samples:"
    cat sample_list.txt
    
    # Split VCF by sample - only PASS variants where sample has ALT
    # -f "PASS,." includes:
    #   PASS = explicitly passed all filters
    #   .    = unfiltered sites (e.g. from upstream callers that don't apply filters)
    # Excludes: FS_filter, SOR_filter, QD_filter, SnpCluster etc.
    while read sample_id; do
        echo "Extracting PASS variants for: \${sample_id}"
        
        bcftools view \\
            --samples \${sample_id} \\
            --apply-filters "PASS,." \\
            --output-type u \\
            ${vcf} \\
        | bcftools view \\
            --include 'GT="alt"' \\
            --output-type z \\
            --output \${sample_id}.per_sample.vcf.gz
        
        # Index
        bcftools index --tbi \${sample_id}.per_sample.vcf.gz
        
        # Count variants
        count=\$(bcftools view -H \${sample_id}.per_sample.vcf.gz | wc -l)
        echo "  PASS variants for \${sample_id}: \${count}"
        
    done < sample_list.txt
    
    echo "VCF splitting complete"
    
    # Generate versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
        samples_processed: ${samples.size()}
    END_VERSIONS
    """
    
    stub:
    """
    cat <<EOF > sample_list.txt
${samples.join('\n')}
EOF
    
    for sample in ${samples.join(' ')}; do
        touch \${sample}.per_sample.vcf.gz
        touch \${sample}.per_sample.vcf.gz.tbi
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.21
        samples_processed: ${samples.size()}
    END_VERSIONS
    """
}
