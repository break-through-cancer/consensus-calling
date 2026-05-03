nextflow.enable.dsl=2

params.vcf1 = params.vcf1 ?: null
params.vcf2 = params.vcf2 ?: null
params.vcf3 = params.vcf3 ?: null
params.vcf4 = params.vcf4 ?: null
params.vcf5 = params.vcf5 ?: null
params.vcf6 = params.vcf6 ?: null
params.vcf7 = params.vcf7 ?: null
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
    bcftools merge \\
      --force-samples \\
      --file-list indel_vcfs.list \\
      -m none \\
      -Oz \\
      -o merged_indels.vcf.gz

    bcftools index -f -t merged_indels.vcf.gz

    echo "=== DEBUG: INDEL records passing liberal mis/ref support but not strict ALT support ==="

    bcftools view \\
      -i "N_PASS(GT!='mis' && GT!='ref')>=${params.indel_min_callers} && N_PASS(GT~\\"[1-9]\\")<${params.indel_min_callers}" \\
      merged_indels.vcf.gz | \\
      bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' | \\
      head -50

    echo "=== ALT-support histogram from merged INDELs ==="

    bcftools query -f '[%GT\\t]\\n' merged_indels.vcf.gz | \\
    awk '
    {
        c=0
        for(i=1;i<=NF;i++) {
            gt=\$i
            if (gt ~ /[1-9]/) c++
        }
        print c
    }' | sort -n | uniq -c | \\
    awk '{print \$2 "\\t" \$1}' \\
    > indel_support_histogram.tsv

    echo "=== BUILD INDEL LIBERAL CONSENSUS: require non-missing/non-reference GT support >= ${params.indel_min_callers} ==="

    bcftools view \\
      -i "N_PASS(GT!='mis' && GT!='ref')>=${params.indel_min_callers}" \\
      -Oz \\
      -o indel_consensus.vcf.gz \\
      merged_indels.vcf.gz

    bcftools index -f -t indel_consensus.vcf.gz

    echo -n "INDEL liberal consensus count: "
    bcftools view -H indel_consensus.vcf.gz | wc -l

    echo "=== CONSENSUS COMPARISON (INDELS) ==="

    strict_alt_count=\$(awk -v m="${params.indel_min_callers}" '
    \$1 >= m {s += \$2}
    END {print s+0}
    ' indel_support_histogram.tsv)

    liberal_consensus_count=\$(bcftools view -H indel_consensus.vcf.gz | wc -l)

    echo "Strict ALT-support count from histogram: \$strict_alt_count"
    echo "Actual liberal consensus count:         \$liberal_consensus_count"
    """
}


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

    # ===================== DEBUG BLOCK =====================
    echo "=== DEBUG: mis/ref vs ALT-support disagreement ==="

    bcftools view \
      -i "N_PASS(GT!='mis' && GT!='ref')>=\$min_callers && N_PASS(GT~\"[1-9]\")<\$min_callers" \
      merged_snvs.vcf.gz | \
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | \
      head -50
    # ======================================================

    echo "=== ALT-support histogram from merged SNVs ==="

    bcftools query -f '[%GT\t]\n' merged_snvs.vcf.gz | \
    awk '
    {
        c=0
        for(i=1;i<=NF;i++) {
            gt=\$i
            if (gt ~ /[1-9]/) c++
        }
        print c
    }' | sort -n | uniq -c | \
    awk '{print \$2 "\t" \$1}' \
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


process LIGHT_FILTER {
    tag "${caller_id}"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(caller_id), path(vcf), path(tbi)

    output:
    tuple val(caller_id), path("${caller_id}.filtered.vcf.gz"), path("${caller_id}.filtered.vcf.gz.tbi")

    script:
    """
    set -euo pipefail

    bcftools view \
      -f 'PASS,.' \
      -i 'GT!="mis" && GT!="ref" && FORMAT/AD[1] > 3 && FORMAT/AF[0] > 0.02' \
      -Oz \
      -o ${caller_id}.filtered.vcf.gz \
      ${vcf}

    bcftools index -t ${caller_id}.filtered.vcf.gz

    echo -n "${caller_id} filtered count: "
    bcftools view -H ${caller_id}.filtered.vcf.gz | wc -l
    """
}


process MERGE_CONSENSUS {
    tag "merge_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    input:
    tuple path(snv_vcf), path(snv_tbi)
    tuple path(indel_vcf), path(indel_tbi)

    output:
    tuple path("final_consensus.vcf.gz"), path("final_consensus.vcf.gz.tbi")

    script:
    """
    bcftools concat -a \
      ${snv_vcf} \
      ${indel_vcf} \
      -Oz -o final_consensus.vcf.gz

    bcftools index -t final_consensus.vcf.gz
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
    """
    echo "=== ${caller_id}: DETERMINE TUMOR SAMPLE ==="

    bcftools query -l ${vcf} > samples.txt
    echo "Samples:"
    cat samples.txt

    if echo "${caller_id}" | grep -qi "strelka"; then
        selected_sample="TUMOR"
    else
        selected_sample=\$(grep -viE 'PBMC|NORMAL|BLOOD' samples.txt | head -n 1)
    fi

    if [ -z "\$selected_sample" ]; then
        echo "ERROR: could not identify tumor sample for ${caller_id}" >&2
        cat samples.txt >&2
        exit 1
    fi

    if ! grep -Fxq "\$selected_sample" samples.txt; then
        echo "ERROR: selected sample \$selected_sample not found in ${caller_id}" >&2
        cat samples.txt >&2
        exit 1
    fi

    echo "Selected tumor sample: \$selected_sample"

    bcftools view -s "\$selected_sample" ${vcf} -Oz -o ${caller_id}.tumor_only.vcf.gz
    bcftools index -t ${caller_id}.tumor_only.vcf.gz

    echo "=== QC ${caller_id} ==="
    echo -n "RAW total: "
    bcftools view -H ${vcf} | wc -l

    echo "=== NORMALIZING ${caller_id} ==="
    bcftools norm -f ${ref_fasta} -m -any ${caller_id}.tumor_only.vcf.gz -Oz -o ${caller_id}.norm.vcf.gz
    bcftools index -t ${caller_id}.norm.vcf.gz

    echo -n "NORMALIZED total: "
    bcftools view -H ${caller_id}.norm.vcf.gz | wc -l

    echo "=== SPLITTING SNVS ==="
    bcftools view -v snps ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.snvs.vcf.gz
    bcftools index -t ${caller_id}.snvs.vcf.gz

    echo "=== SPLITTING INDELS ==="
    bcftools view -v indels ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.indels.vcf.gz
    bcftools index -t ${caller_id}.indels.vcf.gz

    echo -n "SNVs: "
    bcftools view -H ${caller_id}.snvs.vcf.gz | wc -l

    echo -n "INDELs: "
    bcftools view -H ${caller_id}.indels.vcf.gz | wc -l

    echo -n "PASS (SNVs): "
    bcftools view -H -f PASS,. ${caller_id}.snvs.vcf.gz | wc -l

    echo -n "PASS (INDELs): "
    bcftools view -H -f PASS,. ${caller_id}.indels.vcf.gz | wc -l
    """
}


//TODO: you could make it such that you input the tbis as well as the file, but maybe this is somethign we can solve later instead of creating the tbi everyt ime
process ENSURE_INDEX {
    tag "${caller_id}"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(caller_id), path(vcf)

    output:
    tuple val(caller_id), path(vcf), path("${vcf}.tbi")

    script:
    """
    bcftools index -t ${vcf}
    """
}


workflow {

    vcfs = [
        params.vcf1, params.vcf2, params.vcf3,
        params.vcf4, params.vcf5, params.vcf6, params.vcf7
    ].findAll { it }

    if (vcfs.size() < 3) {
        error "Need at least 3 VCFs for consensus"
    }

    println "=== RAW PARAM INPUTS ==="
    println vcfs
    println "Total VCFs: ${vcfs.size()}"

    vcf_ch = Channel
    .fromList(vcfs)
    .map { vcf ->
        def base = vcf.tokenize('/')[-1]              // filename
        def caller_id = base.replaceAll(/\.vcf\.gz$/, '')  // strip suffix

        tuple(caller_id, file(vcf))
    }

    vcf_ch.view { "VCF_CH → caller_id=${it[0]}, file=${it[1]}" }

    indexed_ch = ENSURE_INDEX(vcf_ch)

    //dont light filter for now
    filtered_ch =indexed_ch

    split_ch = SPLIT_SNVS_INDELS(
        filtered_ch,
        file(params.ref_fasta),
        file(params.ref_fai)
    )

    snv_vcfs = split_ch.snvs.map { caller_id, vcf, tbi -> vcf }.collect()
    snv_tbis = split_ch.snvs.map { caller_id, vcf, tbi -> tbi }.collect()

    indel_vcfs = split_ch.indels.map { caller_id, vcf, tbi -> vcf }.collect()
    indel_tbis = split_ch.indels.map { caller_id, vcf, tbi -> tbi }.collect()

    snv_consensus = CONSENSUS_SNVS(snv_vcfs, snv_tbis)
    indel_consensus = CONSENSUS_INDELS(indel_vcfs, indel_tbis)

    MERGE_CONSENSUS(snv_consensus.consensus, indel_consensus.consensus)
}

