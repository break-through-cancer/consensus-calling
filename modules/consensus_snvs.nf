process CONSENSUS_SNVS {
    tag "snv_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

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

    bcftools merge --force-samples --file-list snv_vcfs.list -m none -Oz -o merged_snvs.vcf.gz
    bcftools index -t merged_snvs.vcf.gz

    bcftools view \
      -i "N_PASS(GT!=\\"mis\\")>=\$min_callers" \
      -Oz \
      -o snv_consensus.vcf.gz \
      merged_snvs.vcf.gz

    bcftools index -t snv_consensus.vcf.gz
    """
}