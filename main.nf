nextflow.enable.dsl=2

params.samples         = params.samples         ?: null   // JSON string from preprocess.py
params.outdir          = params.outdir          ?: "results"
params.ref_fasta       = params.ref_fasta       ?: null
params.ref_fai         = params.ref_fai         ?: null
params.callers         = params.callers         ?: 'mutect,deepsomatic,strelka'   // callers to use (set from the UI checkboxes by preprocess.py)
// Keep only PASS (and '.') records from each caller. DeepSomatic VCFs contain non-PASS records
// (e.g. RefCall) that would otherwise be counted as caller support.
params.pass_only       = params.containsKey('pass_only') ? params.pass_only : true


// Callers that only report SNVs; they take no part in indel consensus.
SNV_ONLY_CALLERS = ['mutect']

// How callers are combined, from the number k of callers in play:
//   k >= 3  n-1         a variant must be called by at least k-1 callers
//   k == 2  intersection  both callers must call it
//   k == 1  passthrough   the single caller's calls are returned as they are
def consensusMethod(int k) {
    k >= 3 ? [name: 'n-1', min: k - 1]
           : k == 2 ? [name: 'intersection', min: 2]
                    : [name: 'passthrough', min: 1]
}


// ── Helper: parse the JSON sample manifest emitted by preprocess.py ──────────
def parseSamples(json_str) {
    return new groovy.json.JsonSlurper().parseText(json_str)
}


// Strelka (nf-core/sarek) writes SNVs and indels as separate VCFs; join them so Strelka
// looks like every other caller (one VCF per sample).
process CONCAT_STRELKA {
    tag "${sample_id}"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(sample_id), path(snvs), path(indels)

    output:
    tuple val("${sample_id}__strelka"), path("${sample_id}.strelka.vcf.gz")

    script:
    """
    set -euo pipefail
    bcftools index -f -t ${snvs}
    bcftools index -f -t ${indels}
    bcftools concat -a ${snvs} ${indels} -Oz -o ${sample_id}.strelka.vcf.gz
    """
}


// Normalise every caller's VCF to a sorted, bgzipped, indexed file. Handles plain .vcf
// (e.g. MuTect) as well as .vcf.gz inputs.
process ENSURE_INDEX {
    tag "${caller_id}"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(caller_id), path(vcf)

    output:
    tuple val(caller_id), path("${caller_id}.input.vcf.gz"), path("${caller_id}.input.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    bcftools sort -Oz -o ${caller_id}.input.vcf.gz ${vcf}
    bcftools index -f -t ${caller_id}.input.vcf.gz
    """
}


process SPLIT_SNVS_INDELS {
    tag "${caller_id}"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(caller_id), path(vcf), path(tbi)
    path ref_fasta
    path ref_fai

    output:
    tuple val(caller_id), path("${caller_id}.snvs.vcf.gz"), path("${caller_id}.snvs.vcf.gz.tbi"), emit: snvs
    tuple val(caller_id), path("${caller_id}.indels.vcf.gz"), path("${caller_id}.indels.vcf.gz.tbi"), emit: indels

    script:
    def pass_filter = params.pass_only.toString() == 'false' ? '' : '-f PASS,.'
    """
    set -euo pipefail

    echo "[SPLIT_SNVS_INDELS] Caller: ${caller_id}"

    bcftools query -l ${vcf} > samples.txt || true
    grep -viE 'PBMC|NORMAL|BLOOD' samples.txt > tumor_samples.txt || true

    # Strelka names its sample columns NORMAL/TUMOR; otherwise take the first
    # non-normal sample. A VCF with no sample columns is used as is.
    selected_sample=""
    if grep -qx "TUMOR" samples.txt; then
        selected_sample="TUMOR"
    elif [ -s samples.txt ]; then
        selected_sample=\$(head -n 1 tumor_samples.txt)
        if [ -z "\$selected_sample" ]; then
            echo "ERROR: Could not identify tumor sample for ${caller_id}" >&2
            exit 1
        fi
    fi

    echo "[SPLIT_SNVS_INDELS] Selected tumor sample: \${selected_sample:-<none: no sample columns>}"

    if [ -n "\$selected_sample" ]; then
        bcftools view ${pass_filter} -s "\$selected_sample" ${vcf} -Oz -o ${caller_id}.tumor_only.vcf.gz
    else
        bcftools view ${pass_filter} ${vcf} -Oz -o ${caller_id}.tumor_only.vcf.gz
    fi
    bcftools index -f -t ${caller_id}.tumor_only.vcf.gz

    bcftools norm \
        -f ${ref_fasta} \
        -m -any \
        ${caller_id}.tumor_only.vcf.gz \
        -Oz -o ${caller_id}.norm.vcf.gz
    bcftools index -f -t ${caller_id}.norm.vcf.gz

    bcftools view -v snps   ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.snvs.vcf.gz
    bcftools index -f -t ${caller_id}.snvs.vcf.gz

    bcftools view -v indels ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.indels.vcf.gz
    bcftools index -f -t ${caller_id}.indels.vcf.gz
    """
}


process CONSENSUS {
    tag "${sample_id}_${vtype}_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), val(vtype), val(method), val(min_callers), path(vcfs), path(tbis)

    output:
    tuple val(sample_id), path("${vtype}_consensus.vcf.gz"), path("${vtype}_consensus.vcf.gz.tbi"), path("${vtype}_consensus_summary.txt"), emit: consensus
    path("${vtype}_caller_support.tsv"),    emit: support
    path("${vtype}_support_histogram.tsv"), emit: histogram

    script:
    """
    set -euo pipefail

    echo "[CONSENSUS] Sample: ${sample_id}  type: ${vtype}  method: ${method}"

    ls -1 *.${vtype}s.vcf.gz > vcfs.list
    n=\$(wc -l < vcfs.list)
    min_callers=${min_callers}

    echo "[CONSENSUS] Input callers: \$n, min_callers: \$min_callers"

    first_vcf=\$(head -n 1 vcfs.list)

    set +o pipefail
    zcat -f "\$first_vcf" | awk '
        BEGIN { OFS="\t" }
        /^##/ { print; next }
        /^#CHROM/ { print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8; exit }
    ' > consensus.header
    set -o pipefail

    > all_sites.tsv

    # Staged files are named <sample>__<caller>.<type>s.vcf.gz; keep just the caller name.
    while read f; do
        caller="\${f%.${vtype}s.vcf.gz}"
        caller="\${caller##*__}"
        zcat -f "\$f" | awk -v caller="\$caller" '
        BEGIN { OFS="\t" }
        !/^#/ {
            key = \$1 ":" \$2 ":" \$4 ":" \$5
            format = (NF >= 9) ? \$9 : "."
            sample = (NF >= 10) ? \$10 : "."
            print key, caller, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, format, sample
        }' >> all_sites.tsv
    done < vcfs.list

    awk -v m="\$min_callers" '
    BEGIN { OFS="\t" }
    {
        key = \$1; caller = \$2
        if (!seen[key, caller]) {
            seen[key, caller] = 1
            count[key]++
            callers[key] = (callers[key] == "") ? caller : callers[key] "," caller
        }
    }
    END {
        for (key in count) {
            print key, count[key], callers[key] > "support.all.tsv"
            print count[key] > "raw_counts.txt"
            if (count[key] >= m) print key > "passing_keys.txt"
        }
    }' all_sites.tsv

    touch passing_keys.txt raw_counts.txt support.all.tsv

    # support (number of callers) <tab> number of distinct variants with that support
    sort -n raw_counts.txt | uniq -c | awk 'BEGIN{OFS="\t"} {print \$2,\$1}' > ${vtype}_support_histogram.tsv

    awk '
    BEGIN { OFS="\t" }
    FILENAME == "passing_keys.txt" { passing[\$1] = 1; next }
    FILENAME == "support.all.tsv"  { support[\$1] = \$2; callers[\$1] = \$3; next }
    FILENAME == "all_sites.tsv" {
        key = \$1
        if (key in passing) {
            print key, support[key], callers[key], \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12 >> "${vtype}_caller_support.body.tsv"
            if (!(key in emitted)) {
                emitted[key] = 1
                print \$3, \$4, \$5, \$6, \$7, \$8, \$9, "." >> "consensus.body"
            }
        }
    }' passing_keys.txt support.all.tsv all_sites.tsv

    { echo -e "variant_id\tsupport\tcallers\tcaller\tchrom\tpos\tid\tref\talt\tqual\tfilter\tinfo\tformat\tsample"
      cat ${vtype}_caller_support.body.tsv 2>/dev/null || true
    } > ${vtype}_caller_support.tsv

    touch consensus.body
    cat consensus.header consensus.body > unsorted_consensus.vcf
    bcftools sort unsorted_consensus.vcf -Oz -o ${vtype}_consensus.vcf.gz

    # One caller: return its (PASS-filtered, normalised) calls with all their fields.
    if [ "${method}" = "passthrough" ]; then
        bcftools view "\$first_vcf" -Oz -o ${vtype}_consensus.vcf.gz
    fi
    bcftools index -f -t ${vtype}_consensus.vcf.gz

    n_consensus=\$(bcftools view -H ${vtype}_consensus.vcf.gz | wc -l)
    echo "[CONSENSUS] ${sample_id} ${vtype} count: \$n_consensus"

    {
        echo "${vtype} consensus - sample ${sample_id}"
        echo "  callers (\$n): \$(sed -E 's/\\.${vtype}s\\.vcf\\.gz\$//; s/^.*__//' vcfs.list | paste -sd, - | sed 's/,/, /g')"
        echo "  method: ${method} - a variant is kept if called by >= \$min_callers of \$n callers"
        echo "  variants per caller (input):"
        awk -F'\t' '{c[\$2]++} END {for (k in c) printf "    %s: %d\\n", k, c[k]}' all_sites.tsv | sort
        echo "  distinct variants by number of callers supporting them:"
        for s in \$(seq \$n -1 1); do
            cnt=\$(awk -v s=\$s '\$1 == s {print \$2}' ${vtype}_support_histogram.tsv)
            tag=""
            [ "\$s" -ge "\$min_callers" ] && tag="  <- in consensus"
            echo "    \$s of \$n callers: \${cnt:-0}\$tag"
        done
        echo "  variants in consensus: \$n_consensus"
    } > ${vtype}_consensus_summary.txt

    cat ${vtype}_consensus_summary.txt
    """
}


process MERGE_CONSENSUS {
    tag "${sample_id}_merge"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcfs), path(tbis), path(summaries)
    val notes

    output:
    tuple val(sample_id), path("${sample_id}.final_consensus.vcf.gz"), path("${sample_id}.final_consensus.vcf.gz.tbi"), emit: merged
    path("${sample_id}.consensus_log.txt"), emit: consensus_log

    script:
    """
    set -euo pipefail

    ls -1 *_consensus.vcf.gz > parts.list

    if [ \$(wc -l < parts.list) -eq 1 ]; then
        # SNV-only run (only SNV-only callers selected): the SNV consensus is the result.
        bcftools view \$(cat parts.list) -Oz -o ${sample_id}.final_consensus.vcf.gz
    else
        # A passthrough part still carries its caller's sample column while the others are
        # sites-only, and concat needs matching samples, so drop genotypes from both.
        for t in snv indel; do
            bcftools view -G \${t}_consensus.vcf.gz -Oz -o \${t}.sites.vcf.gz
            bcftools index -f -t \${t}.sites.vcf.gz
        done
        bcftools concat -a snv.sites.vcf.gz indel.sites.vcf.gz \
            -Oz -o ${sample_id}.final_consensus.vcf.gz
    fi

    bcftools index -f -t ${sample_id}.final_consensus.vcf.gz

    n_final=\$(bcftools view -H ${sample_id}.final_consensus.vcf.gz | wc -l)
    echo "[MERGE_CONSENSUS] ${sample_id} final count: \$n_final"

    {
        echo "Consensus log - sample ${sample_id}"
        echo
        printf '%s\\n' '${notes}' | tr '|' '\\n' | sed 's/^ */  /'
        echo
        cat snv_consensus_summary.txt
        echo
        if [ -f indel_consensus_summary.txt ]; then
            cat indel_consensus_summary.txt
        else
            echo "indel consensus - not run (no selected caller reports indels)"
        fi
        echo
        echo "final consensus variants (SNVs + indels): \$n_final"
    } > ${sample_id}.consensus_log.txt
    """
}


workflow {
    log.info "========================================="
    log.info "   Multi-Sample Somatic Consensus        "
    log.info "========================================="

    if (!params.samples) {
        error "params.samples is required (JSON string from preprocess.py)"
    }
    if (!params.ref_fasta) {
        error "params.ref_fasta is required"
    }

    def sample_list = parseSamples(params.samples)
    log.info "Processing ${sample_list.size()} samples"

    // Callers chosen in the UI. Every caller reports SNVs; SNV-only callers (MuTect) are
    // left out of the indel consensus.
    def snv_callers   = params.callers.toString().tokenize(',').collect { it.trim() }.findAll { it }
    def indel_callers = snv_callers.findAll { !SNV_ONLY_CALLERS.contains(it) }
    if (!snv_callers) {
        error "params.callers is empty; select at least one caller"
    }

    def snv_plan   = consensusMethod(snv_callers.size())
    def indel_plan = indel_callers ? consensusMethod(indel_callers.size()) : null

    def snv_note = "SNV consensus: ${snv_plan.name} of ${snv_callers.join(', ')} (>= ${snv_plan.min} of ${snv_callers.size()} callers)"
    def indel_note = indel_plan
        ? "indel consensus: ${indel_plan.name} of ${indel_callers.join(', ')} (>= ${indel_plan.min} of ${indel_callers.size()} callers)"
        : "indel consensus: not run, ${snv_callers.join(', ')} only call SNVs"
    def notes = "callers selected: ${snv_callers.join(', ')} | ${snv_note} | ${indel_note}".toString()
    log.info notes.replace(' | ', '\n')

    // Build channel: one entry per (caller_id, vcf_path), where caller_id is
    // "<sample_id>__<caller>" to keep process tags and filenames unique.
    // MuTect and DeepSomatic each provide a single VCF per sample.
    single_vcf_ch = Channel.fromList(
        sample_list.collectMany { s ->
            def rows = []
            if (s.vcf_mutect)      rows << [ "${s.sample_id}__mutect",      file(s.vcf_mutect)      ]
            if (s.vcf_deepsomatic) rows << [ "${s.sample_id}__deepsomatic", file(s.vcf_deepsomatic) ]
            rows
        }
    )

    // Strelka provides separate SNV and indel VCFs, joined into one per sample
    strelka_ch = CONCAT_STRELKA(
        Channel.fromList(
            sample_list
                .findAll { s -> s.vcf_strelka_snvs && s.vcf_strelka_indels }
                .collect { s -> tuple(s.sample_id, file(s.vcf_strelka_snvs), file(s.vcf_strelka_indels)) }
        )
    )

    vcf_ch = single_vcf_ch.mix(strelka_ch)

    indexed_ch = ENSURE_INDEX(vcf_ch)

    split_ch = SPLIT_SNVS_INDELS(
        indexed_ch,
        file(params.ref_fasta),
        file(params.ref_fai)
    )

    // Group each sample's per-caller VCFs (sample_id is everything before __).
    // A sample is released as soon as all of its own callers have finished.
    n_callers_by_sample = sample_list.collectEntries { s -> [ (s.sample_id): s.callers.size() ] }
    snv_in = split_ch.snvs
        .map { caller_id, vcf, tbi ->
            def sample_id = caller_id.split('__')[0]
            tuple(groupKey(sample_id, n_callers_by_sample[sample_id]), vcf, tbi)
        }
        .groupTuple(by: 0)
        .map { sample_id, vcfs, tbis ->
            tuple(sample_id.toString(), 'snv', snv_plan.name, snv_plan.min, vcfs, tbis)
        }

    // Indels: only callers that report indels, and only if any were selected.
    indel_in = Channel.empty()
    if (indel_plan) {
        indel_in = split_ch.indels
            .filter { caller_id, vcf, tbi -> indel_callers.contains(caller_id.split('__')[1]) }
            .map { caller_id, vcf, tbi ->
                tuple(groupKey(caller_id.split('__')[0], indel_callers.size()), vcf, tbi)
            }
            .groupTuple(by: 0)
            .map { sample_id, vcfs, tbis ->
                tuple(sample_id.toString(), 'indel', indel_plan.name, indel_plan.min, vcfs, tbis)
            }
    }

    consensus_ch = CONSENSUS(snv_in.mix(indel_in))

    // Per sample: the SNV consensus plus the indel consensus when there is one.
    n_parts = indel_plan ? 2 : 1
    merge_in = consensus_ch.consensus
        .map { sample_id, vcf, tbi, summary -> tuple(groupKey(sample_id, n_parts), vcf, tbi, summary) }
        .groupTuple(by: 0)
        .map { sample_id, vcfs, tbis, summaries -> tuple(sample_id.toString(), vcfs, tbis, summaries) }

    MERGE_CONSENSUS(merge_in, notes)
}
