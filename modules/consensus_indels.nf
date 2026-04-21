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
    echo "=== INPUT FILES ==="
    ls -lh *.vcf.gz

    echo "=== VARIANT COUNTS PER INPUT ==="
    for f in *.vcf.gz; do
        echo -n "$f: "
        bcftools view -H "$f" | wc -l
    done
    ls -1 *.vcf.gz > indel_vcfs.list

    bcftools merge --force-samples --file-list indel_vcfs.list -m none -Oz -o merged_indels.vcf.gz
    bcftools index -t merged_indels.vcf.gz

    bcftools view \
      -i "N_PASS(GT!=\\"mis\\")>=$params.indel_min_callers" \
      -Oz \
      -o indel_consensus.vcf.gz \
      merged_indels.vcf.gz

    bcftools index -t indel_consensus.vcf.gz

    echo "=== SUPPORT HISTOGRAM TSV ==="

    bcftools query -f '[%GT\t]\n' merged_indels.vcf.gz | \
    awk '
    {
        c=0
        for(i=1;i<=NF;i++) {
            if(\$i != "./." && \$i != ".") c++
        }
        print c
    }' | sort -n | uniq -c | \
    awk '{print \$2 "\t" \$1}' \
    > indel_support_histogram.tsv
    """
}