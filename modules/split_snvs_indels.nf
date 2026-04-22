process SPLIT_SNVS_INDELS {
    tag "${caller_id}"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(caller_id), path(vcf), path(tbi)
    path ref_fasta
    path ref_fai
    val tumor_sample_name

    output:
    tuple val(caller_id), path("${caller_id}.snvs.vcf.gz"), path("${caller_id}.snvs.vcf.gz.tbi"), emit: snvs
    tuple val(caller_id), path("${caller_id}.indels.vcf.gz"), path("${caller_id}.indels.vcf.gz.tbi"), emit: indels

    script:
    """
    echo "=== QC ${caller_id} ==="

    echo -n "RAW total: "
    bcftools view -H ${vcf} | wc -l

    bcftools query -l ${vcf} > samples.txt

    echo "=== GETTING ONLY TUMOR COLUMN ==="
    if grep -Fxq "${tumor_sample_name}" samples.txt; then
        bcftools view -s "${tumor_sample_name}" ${vcf} -Oz -o ${caller_id}.tumor_only.vcf.gz
    else
        echo "ERROR: tumor sample ${tumor_sample_name} not found in ${vcf}" >&2
        cat samples.txt >&2
        exit 1
    fi

    bcftools index -t ${caller_id}.tumor_only.vcf.gz

    echo "=== NORMALIZING ${caller_id} ==="
    bcftools norm -f ${ref_fasta} -m -any ${caller_id}.tumor_only.vcf.gz -Oz -o ${caller_id}.norm.vcf.gz
    bcftools index -t ${caller_id}.norm.vcf.gz

    echo -n "NORMALIZED total: "
    bcftools view -H ${caller_id}.norm.vcf.gz | wc -l

    echo "=== SPLITTING SNVS ==="
    bcftools view -v snps ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.snvs.vcf.gz
    bcftools index -t ${caller_id}.snvs.vcf.gz

    echo "=== SPLITTING INDELS ==="
    bcftools view -v indels ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.indels.vcf.gz
    bcftools index -t ${caller_id}.indels.vcf.gz

    echo -n "SNVs: "
    bcftools view -H ${caller_id}.snvs.vcf.gz | wc -l

    echo -n "INDELs: "
    bcftools view -H ${caller_id}.indels.vcf.gz | wc -l

    echo -n "PASS (SNVs): "
    bcftools view -H -f PASS,. ${caller_id}.snvs.vcf.gz | wc -l

    echo -n "PASS (INDELs): "
    bcftools view -H -f PASS,. ${caller_id}.indels.vcf.gz | wc -l
    """
}