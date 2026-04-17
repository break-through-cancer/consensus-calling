process CONSENSUS_SNVS {
    tag "snv_consensus"
    container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"
    
    input:
    path vcfs
    path tbis

    output:
    tuple path("snv_consensus.vcf.gz"), path("snv_consensus.vcf.gz.tbi")

    script:
    """
    ls -1 *.vcf.gz > snv_vcfs.list

    n=\$(wc -l < snv_vcfs.list)

    if [ "$params.snv_min_callers" = "null" ] || [ -z "$params.snv_min_callers" ]; then
        min_callers=\$((n-1))
    else
        min_callers=$params.snv_min_callers
    fi

    bcftools merge --file-list snv_vcfs.list -m none -Oz -o merged_snvs.vcf.gz
    bcftools index -t merged_snvs.vcf.gz

    python ${projectDir}/bin/snv_consensus.py merged_snvs.vcf.gz \$min_callers snv_consensus.vcf.gz

    bcftools index -t snv_consensus.vcf.gz
    """
}