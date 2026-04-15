nextflow.enable.dsl=2

include { LIGHT_FILTER } from './modules/light_filter'
include { SPLIT_SNVS_INDELS } from './modules/split_snvs_indels'
include { CONSENSUS_SNVS } from './modules/consensus_snvs'
include { CONSENSUS_INDELS } from './modules/consensus_indels'
include { MERGE_CONSENSUS } from './modules/merge_consensus'

workflow {

    // Collect provided VCFs
    vcfs = [
        params.vcf1, params.vcf2, params.vcf3,
        params.vcf4, params.vcf5, params.vcf6, params.vcf7
    ].findAll { it != null }

    if (vcfs.size() < 3) {
        error "Need at least 3 VCFs for consensus"
    }

    // Build channel with paired TBIs
    vcf_ch = Channel
        .fromList(vcfs)
        .map { vcf ->
            def tbi = file(vcf + ".tbi")
            tuple(vcf.tokenize('/')[-1].replace('.vcf.gz',''), file(vcf), tbi)
        }

    filtered_ch = LIGHT_FILTER(vcf_ch)

    split_ch = SPLIT_SNVS_INDELS(filtered_ch)

    snv_consensus = CONSENSUS_SNVS(split_ch.snvs.collect())
    indel_consensus = CONSENSUS_INDELS(split_ch.indels.collect())

    MERGE_CONSENSUS(snv_consensus, indel_consensus)
}

// process consensus_bcftools {
//   label 'process_medium'
//   container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

//   input:
//     path vcfs
//     val  min_callers

//   output:
//     path "consensus.vcf.gz",               emit: vcf
//     path "consensus.vcf.gz.tbi",           emit: tbi
//     path "consensus_summary.tsv",          emit: summary
//     path "consensus_support_counts.tsv",   emit: support_counts

//   script:
//   """
//   set -euo pipefail

//   echo "=== consensus_bcftools: START ==="
//   echo "Minimum caller support required: ${min_callers}"
//   echo "Input VCFs:"
//   ls -lh *.vcf.gz

//   : > consensus_summary.tsv
//   echo -e "metric\tcount" > consensus_summary.tsv

//   # ── Step 1: normalize each VCF ────────────────────────────────────
//   i=0
//   for vcf in *.vcf.gz; do
//     caller="vcf\$((i+1))"
//     echo "=== Normalizing \$caller : \$vcf ==="

//     bcftools norm -m -any "\$vcf" -O z -o "norm_\${caller}.vcf.gz"
//     bcftools index -t "norm_\${caller}.vcf.gz"

//     raw_total=\$(bcftools view -H "\$vcf" | wc -l)
//     norm_total=\$(bcftools view -H "norm_\${caller}.vcf.gz" | wc -l)
//     snps=\$(bcftools view -H -v snps "norm_\${caller}.vcf.gz" | wc -l)
//     indels=\$(bcftools view -H -v indels "norm_\${caller}.vcf.gz" | wc -l)

//     echo "\$caller raw total variants:       \$raw_total"
//     echo "\$caller normalized total variants: \$norm_total"
//     echo "\$caller SNPs:                     \$snps"
//     echo "\$caller indels:                   \$indels"

//     echo -e "\${caller}_raw_total\t\${raw_total}" >> consensus_summary.tsv
//     echo -e "\${caller}_normalized_total\t\${norm_total}" >> consensus_summary.tsv
//     echo -e "\${caller}_snps\t\${snps}" >> consensus_summary.tsv
//     echo -e "\${caller}_indels\t\${indels}" >> consensus_summary.tsv

//     i=\$((i+1))
//   done

//   echo "=== Finished normalization ==="
//   ls -lh norm_*.vcf.gz

//   # After normalization, before merge — strip FORMAT to GT only
//   for caller in vcf1 vcf2 vcf3; do
//     bcftools annotate \
//       -x FORMAT \
//       norm_\${caller}.vcf.gz \
//       -O z -o stripped_\${caller}.vcf.gz
//     bcftools index -t stripped_\${caller}.vcf.gz
//   done

//   ls stripped_*.vcf.gz > norm.list

//   # ── Step 2: merge callers into one multi-sample VCF ───────────────
//   echo "=== Merging normalized VCFs ==="
//   cat norm.list

//   bcftools merge \
//     --file-list norm.list \
//     --force-samples \
//     -O z -o merged.vcf.gz

//   bcftools index -t merged.vcf.gz

//   merged_total=\$(bcftools view -H merged.vcf.gz | wc -l)
//   echo "Merged union total variants: \$merged_total"
//   echo -e "merged_union_total\t\${merged_total}" >> consensus_summary.tsv

//   # ── Step 3: annotate support count per site ───────────────────────
//   echo "=== Filling NUM_CALLERS tag ==="
//   bcftools +fill-tags merged.vcf.gz -O z -o tagged.vcf.gz -- -t 'NUM_CALLERS=COUNT(GT!="mis")'
//   bcftools index -t tagged.vcf.gz

//   tagged_total=\$(bcftools view -H tagged.vcf.gz | wc -l)
//   echo "Tagged total variants: \$tagged_total"
//   echo -e "tagged_total\t\${tagged_total}" >> consensus_summary.tsv

//   # ── Step 4: summarize overlap/support distribution ────────────────
//   echo "=== Support distribution across callers ==="
//   {
//     echo -e "support_n\tcount"
//     for n in 1 2 3; do
//       c=\$(bcftools view -i "INFO/NUM_CALLERS=\${n}" -H tagged.vcf.gz | wc -l)
//       echo -e "\${n}\t\${c}"
//       echo "Variants supported by \${n} caller(s): \${c}"
//       echo -e "support_\${n}\t\${c}" >> consensus_summary.tsv
//     done
//   } > consensus_support_counts.tsv

//   # optional exact intersection counts for 3 callers
//   c1_only=\$(bcftools view -i 'INFO/NUM_CALLERS=1' -H tagged.vcf.gz | wc -l)
//   c2_overlap=\$(bcftools view -i 'INFO/NUM_CALLERS=2' -H tagged.vcf.gz | wc -l)
//   c3_overlap=\$(bcftools view -i 'INFO/NUM_CALLERS=3' -H tagged.vcf.gz | wc -l)

//   echo "Singleton calls (only one caller): \$c1_only"
//   echo "Shared by exactly two callers:     \$c2_overlap"
//   echo "Shared by all three callers:       \$c3_overlap"

//   # ── Step 5: filter to final consensus ─────────────────────────────
//   echo "=== Filtering final consensus: NUM_CALLERS >= ${min_callers} ==="
//   bcftools view -i "INFO/NUM_CALLERS >= ${min_callers}" tagged.vcf.gz \
//     -O z -o consensus.vcf.gz

//   bcftools index -t consensus.vcf.gz

//   final_total=\$(bcftools view -H consensus.vcf.gz | wc -l)
//   final_snps=\$(bcftools view -H -v snps consensus.vcf.gz | wc -l)
//   final_indels=\$(bcftools view -H -v indels consensus.vcf.gz | wc -l)
//   removed=\$((tagged_total - final_total))

//   echo "=== Final consensus summary ==="
//   echo "Final consensus total: \$final_total"
//   echo "Final consensus SNPs:  \$final_snps"
//   echo "Final consensus indels:\$final_indels"
//   echo "Removed by support filter: \$removed"

//   echo -e "final_consensus_total\t\${final_total}" >> consensus_summary.tsv
//   echo -e "final_consensus_snps\t\${final_snps}" >> consensus_summary.tsv
//   echo -e "final_consensus_indels\t\${final_indels}" >> consensus_summary.tsv
//   echo -e "removed_by_support_filter\t\${removed}" >> consensus_summary.tsv

//   echo "=== consensus_summary.tsv ==="
//   cat consensus_summary.tsv

//   echo "=== consensus_support_counts.tsv ==="
//   cat consensus_support_counts.tsv

//   echo "=== consensus_bcftools: END ==="
//   """
// }

//4/13/26 retire
// process consensus_bcftools {
//   label 'process_medium'
//   container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

//   input:
//     path vcfs
//     val  min_callers

//   output:
//     path "consensus.vcf.gz",               emit: vcf
//     path "consensus.vcf.gz.tbi",           emit: tbi
//     path "consensus_summary.tsv",          emit: summary
//     path "consensus_support_counts.tsv",   emit: support_counts

//   script:
//   """
//   set -euo pipefail

//   echo "=== consensus_bcftools: START ==="
//   echo "Minimum caller support required: ${min_callers}"
//   echo "Input VCFs:"
//   ls -lh *.vcf.gz

//   : > consensus_summary.tsv
//   echo -e "metric\tcount" > consensus_summary.tsv

//   # ── Step 1: normalize each VCF ────────────────────────────────────
//   i=0
//   for vcf in *.vcf.gz; do
//     caller="vcf\$((i+1))"
//     echo "=== Normalizing \$caller : \$vcf ==="

//     bcftools norm -m -any "\$vcf" -O z -o "norm_\${caller}.vcf.gz"
//     bcftools index -t "norm_\${caller}.vcf.gz"

//     raw_total=\$(bcftools view -H "\$vcf" | wc -l)
//     norm_total=\$(bcftools view -H "norm_\${caller}.vcf.gz" | wc -l)
//     snps=\$(bcftools view -H -v snps "norm_\${caller}.vcf.gz" | wc -l)
//     indels=\$(bcftools view -H -v indels "norm_\${caller}.vcf.gz" | wc -l)

//     echo "\$caller raw total variants:       \$raw_total"
//     echo "\$caller normalized total variants: \$norm_total"
//     echo "\$caller SNPs:                     \$snps"
//     echo "\$caller indels:                   \$indels"

//     echo -e "\${caller}_raw_total\t\${raw_total}" >> consensus_summary.tsv
//     echo -e "\${caller}_normalized_total\t\${norm_total}" >> consensus_summary.tsv
//     echo -e "\${caller}_snps\t\${snps}" >> consensus_summary.tsv
//     echo -e "\${caller}_indels\t\${indels}" >> consensus_summary.tsv

//     i=\$((i+1))
//   done

//   echo "=== Finished normalization ==="
//   ls -lh norm_*.vcf.gz

//   # ── Step 1b: check sample names before stripping ──────────────────
//   echo "=== Sample names in normalized VCFs (before reheader) ==="
//   for caller in vcf1 vcf2 vcf3; do
//     echo -n "  \${caller} raw header sample col: "
//     bcftools view -h norm_\${caller}.vcf.gz | tail -1
//   done

//   # ── Step 1c: reheader to caller name + strip FORMAT ───────────────
//   echo "=== Reheadering and stripping FORMAT ==="
//   for caller in vcf1 vcf2 vcf3; do
//     # Get the current sample name (last column of header)
//     sample_name=\$(bcftools view -h norm_\${caller}.vcf.gz | tail -1 | awk '{print \$NF}')
//     echo "  \${caller}: renaming sample '\${sample_name}' -> '\${caller}'"

//     bcftools reheader \
//       --samples <(echo "\${sample_name} \${caller}") \
//       norm_\${caller}.vcf.gz \
//     | bcftools annotate \
//       -x FORMAT \
//       -O z -o stripped_\${caller}.vcf.gz

//     bcftools index -t stripped_\${caller}.vcf.gz
//   done

//   # ── Step 1d: verify stripping and reheader worked ─────────────────
//   echo "=== Sample names after reheader (should be vcf1 vcf2 vcf3) ==="
//   for caller in vcf1 vcf2 vcf3; do
//     echo -n "  \${caller} header sample: "
//     bcftools view -h stripped_\${caller}.vcf.gz | tail -1
//   done

//   echo "=== FORMAT fields after stripping (should all be 0) ==="
//   for caller in vcf1 vcf2 vcf3; do
//     format_lines=\$(bcftools view -h stripped_\${caller}.vcf.gz | grep "^##FORMAT" | wc -l)
//     echo "  \${caller} FORMAT lines in header: \${format_lines}  (should be 0)"
//   done

//   echo "=== First variant row per caller (columns 1-10, FORMAT should be absent or GT only) ==="
//   for caller in vcf1 vcf2 vcf3; do
//     echo -n "  \${caller}: "
//     bcftools view -H stripped_\${caller}.vcf.gz | head -1 | cut -f1-10
//   done

//   ls stripped_*.vcf.gz > norm.list

//   echo "=== norm.list contents (should be 3 stripped files) ==="
//   cat norm.list
//   echo "  line count: \$(wc -l < norm.list)  (should be 3)"

//   # ── Step 2: merge callers into one multi-sample VCF ───────────────
//   echo "=== Merging normalized VCFs ==="

//   bcftools merge \
//     --file-list norm.list \
//     --force-samples \
//     -O z -o merged.vcf.gz

//   bcftools index -t merged.vcf.gz

//   # ── Step 2b: verify merge produced 3 distinct sample columns ──────
//   echo "=== Merged VCF header sample columns (should be vcf1 vcf2 vcf3) ==="
//   bcftools view -h merged.vcf.gz | tail -1

//   echo "=== First 5 merged variant rows (cols 9-12) ==="
//   bcftools view -H merged.vcf.gz | head -5 | cut -f9-12

//   merged_total=\$(bcftools view -H merged.vcf.gz | wc -l)
//   echo "Merged union total variants: \$merged_total"
//   echo "  (expect roughly union of all 3 callers, NOT sum)"
//   echo -e "merged_union_total\t\${merged_total}" >> consensus_summary.tsv

//   # ── Step 3: annotate support count per site ───────────────────────
//   echo "=== Filling NUM_CALLERS tag ==="
//   bcftools +fill-tags merged.vcf.gz -O z -o tagged.vcf.gz -- -t 'NUM_CALLERS=COUNT(GT!="mis")'
//   bcftools index -t tagged.vcf.gz

//   tagged_total=\$(bcftools view -H tagged.vcf.gz | wc -l)
//   echo "Tagged total variants: \$tagged_total"
//   echo -e "tagged_total\t\${tagged_total}" >> consensus_summary.tsv

//   # ── Step 3b: spot check NUM_CALLERS tag on a few rows ─────────────
//   echo "=== First 5 tagged rows — check NUM_CALLERS in INFO ==="
//   bcftools view -H tagged.vcf.gz | head -5 | cut -f1-2,8 | grep -o 'NUM_CALLERS=[0-9]*' || echo "  (NUM_CALLERS tag not found — check fill-tags step)"

//   # ── Step 4: summarize overlap/support distribution ────────────────
//   echo "=== Support distribution across callers ==="
//   {
//     echo -e "support_n\tcount"
//     for n in 1 2 3; do
//       c=\$(bcftools view -i "INFO/NUM_CALLERS=\${n}" -H tagged.vcf.gz | wc -l)
//       echo -e "\${n}\t\${c}"
//       echo "Variants supported by \${n} caller(s): \${c}"
//       echo -e "support_\${n}\t\${c}" >> consensus_summary.tsv
//     done
//   } > consensus_support_counts.tsv

//   c1_only=\$(bcftools view -i 'INFO/NUM_CALLERS=1' -H tagged.vcf.gz | wc -l)
//   c2_overlap=\$(bcftools view -i 'INFO/NUM_CALLERS=2' -H tagged.vcf.gz | wc -l)
//   c3_overlap=\$(bcftools view -i 'INFO/NUM_CALLERS=3' -H tagged.vcf.gz | wc -l)

//   echo "Singleton calls (only one caller): \$c1_only"
//   echo "Shared by exactly two callers:     \$c2_overlap"
//   echo "Shared by all three callers:       \$c3_overlap"

//   # sanity check: support_1 + support_2 + support_3 should == tagged_total
//   total_check=\$(( c1_only + c2_overlap + c3_overlap ))
//   echo "=== Sanity check: 1+2+3 caller counts should sum to tagged_total ==="
//   echo "  \${c1_only} + \${c2_overlap} + \${c3_overlap} = \${total_check}  (tagged_total=\${tagged_total})"
//   if [ "\${total_check}" -ne "\${tagged_total}" ]; then
//     echo "  !! WARNING: counts do not sum correctly — NUM_CALLERS tag may be wrong"
//   else
//     echo "  OK: counts sum correctly"
//   fi

//   # ── Step 5: filter to final consensus ─────────────────────────────
//   echo "=== Filtering final consensus: NUM_CALLERS >= ${min_callers} ==="
//   bcftools view -i "INFO/NUM_CALLERS >= ${min_callers}" tagged.vcf.gz \
//     -O z -o consensus.vcf.gz

//   bcftools index -t consensus.vcf.gz

//   final_total=\$(bcftools view -H consensus.vcf.gz | wc -l)
//   final_snps=\$(bcftools view -H -v snps consensus.vcf.gz | wc -l)
//   final_indels=\$(bcftools view -H -v indels consensus.vcf.gz | wc -l)
//   removed=\$((tagged_total - final_total))

//   # sanity check: final should equal c2_overlap + c3_overlap when min_callers=2
//   expected_final=\$(( c2_overlap + c3_overlap ))
//   echo "=== Sanity check: final_total should equal support_2 + support_3 ==="
//   echo "  \${c2_overlap} + \${c3_overlap} = \${expected_final}  (final_total=\${final_total})"
//   if [ "\${final_total}" -ne "\${expected_final}" ]; then
//     echo "  !! WARNING: final count does not match expected — check filter expression"
//   else
//     echo "  OK: final count matches expected"
//   fi

//   echo "=== Final consensus summary ==="
//   echo "Final consensus total: \$final_total"
//   echo "Final consensus SNPs:  \$final_snps"
//   echo "Final consensus indels:\$final_indels"
//   echo "Removed by support filter: \$removed"

//   echo -e "final_consensus_total\t\${final_total}" >> consensus_summary.tsv
//   echo -e "final_consensus_snps\t\${final_snps}" >> consensus_summary.tsv
//   echo -e "final_consensus_indels\t\${final_indels}" >> consensus_summary.tsv
//   echo -e "removed_by_support_filter\t\${removed}" >> consensus_summary.tsv

//   echo "=== consensus_summary.tsv ==="
//   cat consensus_summary.tsv

//   echo "=== consensus_support_counts.tsv ==="
//   cat consensus_support_counts.tsv

//   echo "=== consensus_bcftools: END ==="
//   """
// }

// workflow {

//   if( !params.vcf1 || !params.vcf2 || !params.vcf3 ) {
//     error "You must provide --vcf1, --vcf2, and --vcf3"
//   }

//   vcf_ch = Channel.of(
//     file(params.vcf1, checkIfExists: true),
//     file(params.vcf2, checkIfExists: true),
//     file(params.vcf3, checkIfExists: true)
//   ).collect()

//   consensus_bcftools(
//     vcf_ch,
//     params.snv_min_callers ?: 2   // default: 2-of-3 support
//   )
// }