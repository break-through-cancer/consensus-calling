nextflow.enable.dsl=2

params.samples         = params.samples         ?: null   // JSON string from preprocess.py
params.outdir          = params.outdir          ?: "results"
params.ref_fasta       = params.ref_fasta       ?: null
params.ref_fai         = params.ref_fai         ?: null
params.snv_min_callers = params.snv_min_callers ?: null   // null = n_callers - 1
params.indel_min_callers = params.indel_min_callers ?: 2
params.n_callers       = params.n_callers       ?: 4      // callers expected per sample (from preprocess.py)
params.indel_callers   = params.indel_callers   ?: 'mutect1,deepsomatic,strelka,varscan'   // callers counted for indels
// Keep only PASS (and '.') records from each caller. DeepSomatic VCFs contain non-PASS records
// (e.g. RefCall) that would otherwise be counted as caller support.
params.pass_only       = params.containsKey('pass_only') ? params.pass_only : true


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
// (e.g. MuTect1) as well as .vcf.gz inputs.
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

    # Strelka and VarScan name their sample columns NORMAL/TUMOR; otherwise take the first
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


process CONSENSUS_SNVS {
    tag "${sample_id}_snv_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcfs), path(tbis)

    output:
    tuple val(sample_id), path("snv_consensus.vcf.gz"), path("snv_consensus.vcf.gz.tbi"), emit: consensus
    path("snv_caller_support.tsv"),    emit: support
    path("snv_support_histogram.tsv"), emit: histogram

    script:
    def n_callers = vcfs instanceof List ? vcfs.size() : 1
    def min_callers = params.snv_min_callers ? params.snv_min_callers as int : n_callers - 1
    """
    set -euo pipefail

    echo "[CONSENSUS_SNVS] Sample: ${sample_id}"

    ls -1 *.snvs.vcf.gz > snv_vcfs.list
    n=\$(wc -l < snv_vcfs.list)
    min_callers=${min_callers}

    echo "[CONSENSUS_SNVS] Input callers: \$n, min_callers: \$min_callers"

    first_vcf=\$(head -n 1 snv_vcfs.list)

    set +o pipefail
    zcat -f "\$first_vcf" | awk '
        BEGIN { OFS="\t" }
        /^##/ { print; next }
        /^#CHROM/ { print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8; exit }
    ' > consensus.header
    set -o pipefail

    > all_sites.tsv

    while read f; do
        caller="\${f%.snvs.vcf.gz}"
        zcat -f "\$f" | awk -v caller="\$caller" '
        BEGIN { OFS="\t" }
        !/^#/ {
            key = \$1 ":" \$2 ":" \$4 ":" \$5
            format = (NF >= 9) ? \$9 : "."
            sample = (NF >= 10) ? \$10 : "."
            print key, caller, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, format, sample
        }' >> all_sites.tsv
    done < snv_vcfs.list

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

    sort -n raw_counts.txt | uniq -c | awk 'BEGIN{OFS="\t"} {print \$2,\$1}' > snv_support_histogram.tsv

    awk '
    BEGIN { OFS="\t" }
    FILENAME == "passing_keys.txt" { passing[\$1] = 1; next }
    FILENAME == "support.all.tsv"  { support[\$1] = \$2; callers[\$1] = \$3; next }
    FILENAME == "all_sites.tsv" {
        key = \$1
        if (key in passing) {
            print key, support[key], callers[key], \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12 >> "snv_caller_support.body.tsv"
            if (!(key in emitted)) {
                emitted[key] = 1
                print \$3, \$4, \$5, \$6, \$7, \$8, \$9, "." >> "consensus.body"
            }
        }
    }' passing_keys.txt support.all.tsv all_sites.tsv

    { echo -e "variant_id\tsupport\tcallers\tcaller\tchrom\tpos\tid\tref\talt\tqual\tfilter\tinfo\tformat\tsample"
      cat snv_caller_support.body.tsv 2>/dev/null || true
    } > snv_caller_support.tsv

    touch consensus.body
    cat consensus.header consensus.body > unsorted_consensus.vcf
    bcftools sort unsorted_consensus.vcf -Oz -o snv_consensus.vcf.gz
    bcftools index -f -t snv_consensus.vcf.gz

    echo -n "[CONSENSUS_SNVS] ${sample_id} SNV count: "
    bcftools view -H snv_consensus.vcf.gz | wc -l
    """
}


process CONSENSUS_INDELS {
    tag "${sample_id}_indel_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcfs), path(tbis)

    output:
    tuple val(sample_id), path("indel_consensus.vcf.gz"), path("indel_consensus.vcf.gz.tbi"), emit: consensus
    path("indel_caller_support.tsv"),    emit: support
    path("indel_support_histogram.tsv"), emit: histogram

    script:
    """
    set -euo pipefail

    echo "[CONSENSUS_INDELS] Sample: ${sample_id}"

    # The callers counted for indels are chosen in the workflow (params.indel_callers);
    # use every indel VCF staged here.
    ls -1 *.indels.vcf.gz > filtered_indel_vcfs.list

    n=\$(wc -l < filtered_indel_vcfs.list)
    min_callers=${params.indel_min_callers}

    echo "[CONSENSUS_INDELS] Indel callers: \$n, min_callers: \$min_callers"

    first_vcf=\$(head -n 1 filtered_indel_vcfs.list)

    set +o pipefail
    zcat -f "\$first_vcf" | awk '
        BEGIN { OFS="\t" }
        /^##/ { print; next }
        /^#CHROM/ { print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8; exit }
    ' > consensus.header
    set -o pipefail

    > all_sites.tsv

    while read f; do
        caller="\${f%.indels.vcf.gz}"
        zcat -f "\$f" | awk -v caller="\$caller" '
        BEGIN { OFS="\t" }
        !/^#/ {
            key = \$1 ":" \$2 ":" \$4 ":" \$5
            format = (NF >= 9) ? \$9 : "."
            sample = (NF >= 10) ? \$10 : "."
            print key, caller, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, format, sample
        }' >> all_sites.tsv
    done < filtered_indel_vcfs.list

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

    sort -n raw_counts.txt | uniq -c | awk 'BEGIN{OFS="\t"} {print \$2,\$1}' > indel_support_histogram.tsv

    awk '
    BEGIN { OFS="\t" }
    FILENAME == "passing_keys.txt" { passing[\$1] = 1; next }
    FILENAME == "support.all.tsv"  { support[\$1] = \$2; callers[\$1] = \$3; next }
    FILENAME == "all_sites.tsv" {
        key = \$1
        if (key in passing) {
            print key, support[key], callers[key], \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12 >> "indel_caller_support.body.tsv"
            if (!(key in emitted)) {
                emitted[key] = 1
                print \$3, \$4, \$5, \$6, \$7, \$8, \$9, "." >> "consensus.body"
            }
        }
    }' passing_keys.txt support.all.tsv all_sites.tsv

    { echo -e "variant_id\tsupport\tcallers\tcaller\tchrom\tpos\tid\tref\talt\tqual\tfilter\tinfo\tformat\tsample"
      cat indel_caller_support.body.tsv 2>/dev/null || true
    } > indel_caller_support.tsv

    touch consensus.body
    cat consensus.header consensus.body > unsorted_consensus.vcf
    bcftools sort unsorted_consensus.vcf -Oz -o indel_consensus.vcf.gz
    bcftools index -f -t indel_consensus.vcf.gz

    echo -n "[CONSENSUS_INDELS] ${sample_id} INDEL count: "
    bcftools view -H indel_consensus.vcf.gz | wc -l
    """
}


process MERGE_CONSENSUS {
    tag "${sample_id}_merge"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id),
          path(snv_vcf),   path(snv_tbi),
          path(indel_vcf), path(indel_tbi)

    output:
    tuple val(sample_id), path("${sample_id}.final_consensus.vcf.gz"), path("${sample_id}.final_consensus.vcf.gz.tbi")

    script:
    """
    set -euo pipefail

    bcftools concat -a \
        ${snv_vcf} \
        ${indel_vcf} \
        -Oz -o ${sample_id}.final_consensus.vcf.gz

    bcftools index -f -t ${sample_id}.final_consensus.vcf.gz

    echo -n "[MERGE_CONSENSUS] ${sample_id} final count: "
    bcftools view -H ${sample_id}.final_consensus.vcf.gz | wc -l
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

    // Build channel: one entry per (caller_id, vcf_path), where caller_id is
    // "<sample_id>__<caller>" to keep process tags and filenames unique.
    // MuTect1, DeepSomatic and VarScan each provide a single VCF per sample.
    single_vcf_ch = Channel.fromList(
        sample_list.collectMany { s ->
            [
                [ "${s.sample_id}__mutect1",     file(s.vcf_mutect1)     ],
                [ "${s.sample_id}__deepsomatic", file(s.vcf_deepsomatic) ],
                [ "${s.sample_id}__varscan",     file(s.vcf_varscan)     ],
            ]
        }
    )

    // Strelka provides separate SNV and indel VCFs, joined into one per sample
    strelka_ch = CONCAT_STRELKA(
        Channel.fromList(
            sample_list.collect { s ->
                tuple(s.sample_id, file(s.vcf_strelka_snvs), file(s.vcf_strelka_indels))
            }
        )
    )

    vcf_ch = single_vcf_ch.mix(strelka_ch)

    indexed_ch = ENSURE_INDEX(vcf_ch)

    split_ch = SPLIT_SNVS_INDELS(
        indexed_ch,
        file(params.ref_fasta),
        file(params.ref_fai)
    )

    // Extract sample_id from caller_id (everything before __)
    // Then group the SNV VCFs of all callers per sample
    n_callers = params.n_callers as int
    snv_ch = split_ch.snvs
        .map { caller_id, vcf, tbi ->
            def sample_id = caller_id.split('__')[0]
            tuple(sample_id, vcf, tbi)
        }
        .groupTuple(by: 0, size: n_callers)
        .map { sample_id, vcfs, tbis ->
            tuple(sample_id, vcfs, tbis)
        }

    // For indels: keep only the callers listed in params.indel_callers per sample
    indel_callers = params.indel_callers.toString().tokenize(',').collect { it.trim() }
    log.info "Indel callers: ${indel_callers.join(', ')}"
    indel_ch = split_ch.indels
        .filter { caller_id, vcf, tbi ->
            indel_callers.contains(caller_id.split('__')[1])
        }
        .map { caller_id, vcf, tbi ->
            def sample_id = caller_id.split('__')[0]
            tuple(sample_id, vcf, tbi)
        }
        .groupTuple(by: 0)
        .map { sample_id, vcfs, tbis ->
            tuple(sample_id, vcfs, tbis)
        }

    snv_consensus_ch   = CONSENSUS_SNVS(snv_ch)
    indel_consensus_ch = CONSENSUS_INDELS(indel_ch)

    // Join SNV and INDEL consensus by sample_id for final merge
    merge_input_ch = snv_consensus_ch.consensus
        .map { sample_id, vcf, tbi -> tuple(sample_id, vcf, tbi) }
        .join(
            indel_consensus_ch.consensus
                .map { sample_id, vcf, tbi -> tuple(sample_id, vcf, tbi) },
            by: 0
        )
        .map { sample_id, snv_vcf, snv_tbi, indel_vcf, indel_tbi ->
            tuple(sample_id, snv_vcf, snv_tbi, indel_vcf, indel_tbi)
        }

    MERGE_CONSENSUS(merge_input_ch)
}
