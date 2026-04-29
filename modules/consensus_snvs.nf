process CONSENSUS_SNVS {
    tag "snv_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    path vcfs
    path tbis

    output:
    tuple path("snv_consensus.vcf.gz"), path("snv_consensus.vcf.gz.tbi"), emit: consensus
    path("snv_support_histogram.tsv"), emit: histogram

    script:
    """
    set -euo pipefail

    echo "=== INPUT FILES ==="
    ls -lh *.vcf.gz

    echo "=== VARIANT COUNTS PER INPUT ==="
    for f in *.vcf.gz; do
        echo -n "\$f: "
        bcftools view -H "\$f" | wc -l
    done

    ls -1 *.vcf.gz > snv_vcfs.list
    n=\$(wc -l < snv_vcfs.list)

    if [ "${params.snv_min_callers}" = "null" ] || [ -z "${params.snv_min_callers}" ]; then
        min_callers=\$((n-1))
    else
        min_callers=${params.snv_min_callers}
    fi

    echo "SNV min_callers=\$min_callers"

    # DO NOT use --missing-to-ref
    bcftools merge \
      --force-samples \
      --file-list snv_vcfs.list \
      -m none \
      -Oz \
      -o merged_snvs.vcf.gz

    bcftools index -t merged_snvs.vcf.gz

    echo "=== ALT-support histogram from merged SNVs ==="

    bcftools query -f '[%GT\t]\n' merged_snvs.vcf.gz | \
    awk '
    {
        c=0
        for(i=1;i<=NF;i++) {
            gt=\$i

            # Count only genotypes containing an ALT allele.
            # Do NOT count 0/0, 0|0, ./., ., or missing.
            if (
                gt != "." &&
                gt != "./." &&
                gt != ".|." &&
                gt != "0/0" &&
                gt != "0|0" &&
                gt ~ /[1-9]/
            ) {
                c++
            }
        }
        print c
    }' | sort -n | uniq -c | \
    awk '{print \$2 "\\t" \$1}' \
    > snv_support_histogram.tsv

    echo "=== BUILD SNV CONSENSUS: require ALT support >= \$min_callers ==="

    bcftools view \
      -i "N_PASS(GT!='mis' && GT!='ref')>=\$min_callers" \
      -Oz \
      -o snv_consensus.vcf.gz \
      merged_snvs.vcf.gz

    bcftools index -t snv_consensus.vcf.gz

    echo -n "SNV consensus count: "
    bcftools view -H snv_consensus.vcf.gz | wc -l

    echo "=== CONSENSUS VALIDATION (SNVS) ==="

    expected_consensus_count=\$(awk -v m="\$min_callers" '
    \$1 >= m {s += \$2}
    END {print s+0}
    ' snv_support_histogram.tsv)

    actual_consensus_count=\$(bcftools view -H snv_consensus.vcf.gz | wc -l)

    echo "Expected from histogram: \$expected_consensus_count"
    echo "Actual VCF count:        \$actual_consensus_count"
    """
}