process CONSENSUS_SNVS {
    tag "snv_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    path vcfs
    path tbis

    output:
    tuple path("snv_consensus.vcf.gz"), path("snv_consensus.vcf.gz.tbi")
    path("snv_support_histogram.tsv")
    script:
    """
    echo "=== INPUT FILES ==="
    ls -lh *.vcf.gz

    echo "=== VARIANT COUNTS PER INPUT ==="
    for f in *.vcf.gz; do
        echo -n "$f: "
        bcftools view -H "$f" | wc -l
    done

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

    echo "=== SUPPORT HISTOGRAM TSV ==="

    bcftools query -f '[%GT\t]\n' merged_snvs.vcf.gz | \
    awk '
    {
        c=0
        for(i=1;i<=NF;i++) {
            if($i != "./." && $i != ".") c++
        }
        print c
    }' | sort -n | uniq -c | \
    awk '{print $2 "\t" $1}' \
    > snv_support_histogram.tsv.tsv
    """
}