process SPLIT_SNVS_INDELS {
    tag "${caller_id}"
    container "${params.gatk_docker ?: 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'}"
    input:
    tuple val(caller_id), path(vcf), path(tbi)

    output:
    tuple val(caller_id), path("${caller_id}.snvs.vcf.gz"), path("${caller_id}.snvs.vcf.gz.tbi"), emit: snvs
    tuple val(caller_id), path("${caller_id}.indels.vcf.gz"), path("${caller_id}.indels.vcf.gz.tbi"), emit: indels

    script:
    """
    bcftools norm -m -any ${vcf} -Oz -o ${caller_id}.norm.vcf.gz
    bcftools index -t ${caller_id}.norm.vcf.gz

    bcftools view -v snps ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.snvs.vcf.gz
    bcftools index -t ${caller_id}.snvs.vcf.gz

    bcftools view -v indels ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.indels.vcf.gz
    bcftools index -t ${caller_id}.indels.vcf.gz
    """
}