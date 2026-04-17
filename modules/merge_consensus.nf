process MERGE_CONSENSUS {
    tag "merge_consensus"

    input:
    tuple path(snv_vcf), path(snv_tbi)
    tuple path(indel_vcf), path(indel_tbi)

    output:
    tuple path("final_consensus.vcf.gz"), path("final_consensus.vcf.gz.tbi")

    script:
    """
    bcftools concat -a \
      ${snv_vcf} \
      ${indel_vcf} \
      -Oz -o final_consensus.vcf.gz

    bcftools index -t final_consensus.vcf.gz
    """
}