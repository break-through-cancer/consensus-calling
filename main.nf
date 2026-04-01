process consensus_bcftools {
  label 'process_medium'
  container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

  input:
    path vcfs
    val  min_callers

  output:
    path "consensus.vcf.gz",               emit: vcf
    path "consensus.vcf.gz.tbi",           emit: tbi
    path "consensus_summary.tsv",          emit: summary
    path "consensus_support_counts.tsv",   emit: support_counts

  script:
  """
  set -euo pipefail

  echo "=== consensus_bcftools: START ==="
  echo "Minimum caller support required: ${min_callers}"
  echo "Input VCFs:"
  ls -lh *.vcf.gz

  : > consensus_summary.tsv
  echo -e "metric\tcount" > consensus_summary.tsv

  # ── Step 1: normalize each VCF ────────────────────────────────────
  i=0
  for vcf in *.vcf.gz; do
    caller="vcf\$((i+1))"
    echo "=== Normalizing \$caller : \$vcf ==="

    bcftools norm -m -any "\$vcf" -O z -o "norm_\${caller}.vcf.gz"
    bcftools index -t "norm_\${caller}.vcf.gz"

    raw_total=\$(bcftools view -H "\$vcf" | wc -l)
    norm_total=\$(bcftools view -H "norm_\${caller}.vcf.gz" | wc -l)
    snps=\$(bcftools view -H -v snps "norm_\${caller}.vcf.gz" | wc -l)
    indels=\$(bcftools view -H -v indels "norm_\${caller}.vcf.gz" | wc -l)

    echo "\$caller raw total variants:       \$raw_total"
    echo "\$caller normalized total variants: \$norm_total"
    echo "\$caller SNPs:                     \$snps"
    echo "\$caller indels:                   \$indels"

    echo -e "\${caller}_raw_total\t\${raw_total}" >> consensus_summary.tsv
    echo -e "\${caller}_normalized_total\t\${norm_total}" >> consensus_summary.tsv
    echo -e "\${caller}_snps\t\${snps}" >> consensus_summary.tsv
    echo -e "\${caller}_indels\t\${indels}" >> consensus_summary.tsv

    i=\$((i+1))
  done

  echo "=== Finished normalization ==="
  ls -lh norm_*.vcf.gz

  # After normalization, before merge — strip FORMAT to GT only
  for caller in vcf1 vcf2 vcf3; do
    bcftools annotate \
      -x FORMAT \
      norm_\${caller}.vcf.gz \
      -O z -o stripped_\${caller}.vcf.gz
    bcftools index -t stripped_\${caller}.vcf.gz
  done

  ls stripped_*.vcf.gz > norm.list

  # ── Step 2: merge callers into one multi-sample VCF ───────────────
  echo "=== Merging normalized VCFs ==="
  cat norm.list

  bcftools merge \
    --file-list norm.list \
    --force-samples \
    -O z -o merged.vcf.gz

  bcftools index -t merged.vcf.gz

  merged_total=\$(bcftools view -H merged.vcf.gz | wc -l)
  echo "Merged union total variants: \$merged_total"
  echo -e "merged_union_total\t\${merged_total}" >> consensus_summary.tsv

  # ── Step 3: annotate support count per site ───────────────────────
  echo "=== Filling NUM_CALLERS tag ==="
  bcftools +fill-tags merged.vcf.gz -O z -o tagged.vcf.gz -- -t 'NUM_CALLERS=COUNT(GT!="mis")'
  bcftools index -t tagged.vcf.gz

  tagged_total=\$(bcftools view -H tagged.vcf.gz | wc -l)
  echo "Tagged total variants: \$tagged_total"
  echo -e "tagged_total\t\${tagged_total}" >> consensus_summary.tsv

  # ── Step 4: summarize overlap/support distribution ────────────────
  echo "=== Support distribution across callers ==="
  {
    echo -e "support_n\tcount"
    for n in 1 2 3; do
      c=\$(bcftools view -i "INFO/NUM_CALLERS=\${n}" -H tagged.vcf.gz | wc -l)
      echo -e "\${n}\t\${c}"
      echo "Variants supported by \${n} caller(s): \${c}"
      echo -e "support_\${n}\t\${c}" >> consensus_summary.tsv
    done
  } > consensus_support_counts.tsv

  # optional exact intersection counts for 3 callers
  c1_only=\$(bcftools view -i 'INFO/NUM_CALLERS=1' -H tagged.vcf.gz | wc -l)
  c2_overlap=\$(bcftools view -i 'INFO/NUM_CALLERS=2' -H tagged.vcf.gz | wc -l)
  c3_overlap=\$(bcftools view -i 'INFO/NUM_CALLERS=3' -H tagged.vcf.gz | wc -l)

  echo "Singleton calls (only one caller): \$c1_only"
  echo "Shared by exactly two callers:     \$c2_overlap"
  echo "Shared by all three callers:       \$c3_overlap"

  # ── Step 5: filter to final consensus ─────────────────────────────
  echo "=== Filtering final consensus: NUM_CALLERS >= ${min_callers} ==="
  bcftools view -i "INFO/NUM_CALLERS >= ${min_callers}" tagged.vcf.gz \
    -O z -o consensus.vcf.gz

  bcftools index -t consensus.vcf.gz

  final_total=\$(bcftools view -H consensus.vcf.gz | wc -l)
  final_snps=\$(bcftools view -H -v snps consensus.vcf.gz | wc -l)
  final_indels=\$(bcftools view -H -v indels consensus.vcf.gz | wc -l)
  removed=\$((tagged_total - final_total))

  echo "=== Final consensus summary ==="
  echo "Final consensus total: \$final_total"
  echo "Final consensus SNPs:  \$final_snps"
  echo "Final consensus indels:\$final_indels"
  echo "Removed by support filter: \$removed"

  echo -e "final_consensus_total\t\${final_total}" >> consensus_summary.tsv
  echo -e "final_consensus_snps\t\${final_snps}" >> consensus_summary.tsv
  echo -e "final_consensus_indels\t\${final_indels}" >> consensus_summary.tsv
  echo -e "removed_by_support_filter\t\${removed}" >> consensus_summary.tsv

  echo "=== consensus_summary.tsv ==="
  cat consensus_summary.tsv

  echo "=== consensus_support_counts.tsv ==="
  cat consensus_support_counts.tsv

  echo "=== consensus_bcftools: END ==="
  """
}
workflow {

  // if( !params.vcf1 || !params.vcf2 || !params.vcf3 ) {
  //   error "You must provide --vcf1, --vcf2, and --vcf3"
  // }

  // vcf_ch = Channel.of(
  //   file(params.vcf1, checkIfExists: true),
  //   file(params.vcf2, checkIfExists: true),
  //   file(params.vcf3, checkIfExists: true)
  // ).collect()

  // consensus_bcftools(
  //   vcf_ch,
  //   params.snv_min_callers ?: 2   // default: 2-of-3 support
  // )
}