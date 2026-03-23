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

  def n = params.n_callers as int

  if (n < 2) {
    error "n_callers must be >= 2 for consensus calling, got ${n}"
  }

  vcf_ch = Channel.fromList(
    (0..<n).collect { i ->
      file(params["caller_vcf_${i}"], checkIfExists: true)
    }
  ).collect()

  consensus_bcftools(
    vcf_ch,
    params.snv_min_callers as int   // n-1, computed in preprocess
  )
}