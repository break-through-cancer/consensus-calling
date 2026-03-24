process consensus_bcftools {
  label 'process_medium'
  container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'

  input:
    path vcfs
    val  min_callers

  output:
    path "consensus.vcf.gz",     emit: vcf
    path "consensus.vcf.gz.tbi", emit: tbi

  script:
  """
  set -euo pipefail

  # ── Step 1: normalize each VCF ────────────────────────────────────
  i=0
  for vcf in *.vcf.gz; do
    bcftools norm -m -any "\$vcf" -O z -o "norm_\${i}.vcf.gz"
    bcftools index -t "norm_\${i}.vcf.gz"
    i=\$((i+1))
  done

  # ── Step 2: merge all callers into one VCF ────────────────────────
  ls norm_*.vcf.gz > norm.list

  bcftools merge \
    --file-list norm.list \
    --missing-to-ref \
    --force-samples \
    -O z -o merged.vcf.gz

  bcftools index -t merged.vcf.gz

  # ── Step 3: count caller support and filter by n-1 ───────────────
  bcftools +fill-tags merged.vcf.gz -O z -o tagged.vcf.gz -- -t 'NUM_CALLERS=COUNT(GT!="mis")'
  bcftools index -t tagged.vcf.gz

  bcftools filter -i 'INFO/NUM_CALLERS >= ${min_callers}' tagged.vcf.gz \
    -O z -o consensus.vcf.gz

  bcftools index -t consensus.vcf.gz
  """
}

workflow {

  // Require explicit inputs
  if( !params.vcf1 || !params.vcf2 ) {
    error "You must provide --vcf1 and --vcf2"
  }

  vcf_ch = Channel.of(
    file(params.vcf1, checkIfExists: true),
    file(params.vcf2, checkIfExists: true)
  ).collect()

  consensus_bcftools(
    vcf_ch,
    params.snv_min_callers ?: 1   // default to 1 (since 2 callers → n-1 = 1)
  )
}