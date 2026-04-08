/*
    GLNEXUS MODULE
    Joint genotyping for DeepVariant GVCFs
    Much faster than GATK's GenotypeGVCFs for DeepVariant output
*/

process GLNEXUS {
    tag "glnexus_joint"
    label 'process_high'
    container "${params.containers.glnexus}"
    
    publishDir "${params.outdir}/variants/glnexus", mode: params.publish_dir_mode, failOnError: false
    
    input:
    path gvcfs
    
    output:
    tuple path("deepvariant.joint.vcf.gz"), path("deepvariant.joint.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    def config = params.deepvariant_model == 'WGS' ? 'DeepVariantWGS' : 'DeepVariantWES'
    def memory_gb = task.memory.toGiga() - 4
    """
    # Create list of GVCFs
    ls *.dv.g.vcf.gz > gvcf_list.txt
    
    echo "=== GLnexus Joint Calling ==="
    echo "Config: ${config}"
    echo "GVCFs: \$(wc -l < gvcf_list.txt)"
    
    # Run GLnexus
    glnexus_cli \\
        --config ${config} \\
        --threads ${task.cpus} \\
        --mem-gbytes ${memory_gb} \\
        --dir \${TMPDIR:-/tmp}/glnexus_db \\
        *.dv.g.vcf.gz \\
    | bcftools view - \\
    | bgzip -@ ${task.cpus} -c > deepvariant.joint.vcf.gz
    
    tabix -p vcf deepvariant.joint.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: \$(glnexus_cli --version 2>&1 | head -1 || echo "1.2.7")
    END_VERSIONS
    """
}
