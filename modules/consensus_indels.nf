process CONSENSUS_INDELS {
    tag "indel_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    path vcfs
    path tbis

    output:
    tuple path("indel_consensus.vcf.gz"), path("indel_consensus.vcf.gz.tbi"), emit: consensus
    path("indel_support_histogram.tsv"), emit: histogram

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

    ls -1 *.vcf.gz > indel_vcfs.list

    echo "INDEL min_callers=${params.indel_min_callers}"

    # DO NOT use --missing-to-ref
    bcftools merge \
      --force-samples \
      --file-list indel_vcfs.list \
      -m none \
      -Oz \
      -o merged_indels.vcf.gz

    bcftools index -t merged_indels.vcf.gz

    echo "=== ALT-support histogram from merged INDELs ==="

    bcftools query -f '[%GT\t]\n' merged_indels.vcf.gz | \
    awk '
    {
        c=0
        for(i=1;i<=NF;i++) {
            gt=\$i

            # Count only ALT genotypes.
            # 0/0 is reference, not caller support.
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
    > indel_support_histogram.tsv

    echo "=== BUILD INDEL CONSENSUS: require ALT support >= ${params.indel_min_callers} ==="

    bcftools view \
      -i "N_PASS(GT!='mis' && GT!='ref')>=${params.indel_min_callers}" \
      -Oz \
      -o indel_consensus.vcf.gz \
      merged_indels.vcf.gz

    bcftools index -t indel_consensus.vcf.gz

    echo -n "INDEL consensus count: "
    bcftools view -H indel_consensus.vcf.gz | wc -l

    echo "=== CONSENSUS VALIDATION (INDELS) ==="

    expected_consensus_count=\$(awk -v m="${params.indel_min_callers}" '
    \$1 >= m {s += \$2}
    END {print s+0}
    ' indel_support_histogram.tsv)

    actual_consensus_count=\$(bcftools view -H indel_consensus.vcf.gz | wc -l)

    echo "Expected from histogram: \$expected_consensus_count"
    echo "Actual VCF count:        \$actual_consensus_count"
    """
}