process LIGHT_FILTER {
    tag "${caller_id}"

    input:
    tuple val(caller_id), path(vcf), path(tbi)

    output:
    tuple val(caller_id), path("${caller_id}.filtered.vcf.gz"), path("${caller_id}.filtered.vcf.gz.tbi")

    script:
    """
    bcftools view \
      -f 'PASS,.' \
      -i 'FORMAT/AD[1] > 3 && FORMAT/AF[0] > 0.02' \
      -Oz \
      -o ${caller_id}.filtered.vcf.gz \
      ${vcf}

    bcftools index -t ${caller_id}.filtered.vcf.gz
    """
}