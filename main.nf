nextflow.enable.dsl=2

// Define pipeline parameters with default fallback values
params.vcf1 = params.vcf1 ?: null
params.vcf2 = params.vcf2 ?: null
params.vcf3 = params.vcf3 ?: null
params.vcf4 = params.vcf4 ?: null
params.vcf5 = params.vcf5 ?: null
params.vcf6 = params.vcf6 ?: null
params.vcf7 = params.vcf7 ?: null
params.outdir = params.outdir ?: "results"
params.qc_outdir = params.qc_outdir ?: "${params.outdir}/qc"

// ==============================================================================
// PROCESS 1: ENSURE_INDEX
// Generates CSI/TBI index files for input VCFs if they are missing or outdated.
// ==============================================================================
process ENSURE_INDEX {
    tag "${caller_id}"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(caller_id), path(vcf)

    output:
    tuple val(caller_id), path(vcf), path("${vcf}.tbi")

    script:
    """
    set -euo pipefail
    echo "[ENSURE_INDEX] Processing caller: ${caller_id}"
    echo "[ENSURE_INDEX] Target file: ${vcf}"
    
    bcftools index -f -t ${vcf}
    
    echo "[ENSURE_INDEX] Successfully indexed ${caller_id}"
    """
}

// ==============================================================================
// PROCESS 2: SPLIT_SNVS_INDELS
// Isolates the tumor sample, drops genotypes to make it sites-only, 
// normalizes variants, and splits them into distinct SNV and INDEL files.
// ==============================================================================
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
    set -euo pipefail
    echo "========================================================="
    echo "[SPLIT_SNVS_INDELS] Starting processing for: ${caller_id}"
    echo "========================================================="
    
    # Extract sample names present in the VCF file
    bcftools query -l ${vcf} > samples.txt
    echo "[SPLIT_SNVS_INDELS] Found samples in VCF:"
    cat samples.txt
    
    # Filter out common control/normal sample names to isolate the tumor sample
    grep -viE 'PBMC|NORMAL|BLOOD' samples.txt > tumor_samples.txt || true
    
    # Strelka outputs rely on standard naming conventions ("TUMOR")
    if echo "${caller_id}" | grep -qi "strelka"; then
        selected_sample="TUMOR"
    else
        selected_sample=\$(head -n 1 tumor_samples.txt)
    fi
    
    if [ -z "\$selected_sample" ]; then
        echo "ERROR: could not identify tumor sample for ${caller_id}" >&2
        cat samples.txt >&2
        exit 1
    fi
    echo "[SPLIT_SNVS_INDELS] Selected tumor sample: \$selected_sample"
    
    # Step 1: Strip out normal samples, keeping only the selected tumor sample
    echo "[SPLIT_SNVS_INDELS] Subsetting to tumor sample only..."
    bcftools view -s "\$selected_sample" ${vcf} -Oz -o ${caller_id}.tumor_only.vcf.gz
    bcftools index -f -t ${caller_id}.tumor_only.vcf.gz
    
    # Step 2: Remove genotype columns (-G) to enable site-level consensus matching
    echo "[SPLIT_SNVS_INDELS] Stripping genotype metrics (-G) for site-only analysis..."
    bcftools view -G ${caller_id}.tumor_only.vcf.gz -Oz -o ${caller_id}.sites_only.vcf.gz
    bcftools index -f -t ${caller_id}.sites_only.vcf.gz
    
    # Step 3: Normalizing variants against the reference and split multiallelics into separate rows
    echo "[SPLIT_SNVS_INDELS] Normalizing variants and splitting multiallelic sites..."
    bcftools norm -f ${ref_fasta} -m -any ${caller_id}.sites_only.vcf.gz -Oz -o ${caller_id}.norm.vcf.gz
    bcftools index -f -t ${caller_id}.norm.vcf.gz
    
    # Step 4: Extract SNVs exclusively
    echo "[SPLIT_SNVS_INDELS] Filtering SNVs..."
    bcftools view -v snps ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.snvs.vcf.gz
    bcftools index -f -t ${caller_id}.snvs.vcf.gz
    
    # Step 5: Extract INDELs exclusively
    echo "[SPLIT_SNVS_INDELS] Filtering INDELs..."
    bcftools view -v indels ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.indels.vcf.gz
    bcftools index -f -t ${caller_id}.indels.vcf.gz
    
    echo -n "[SPLIT_SNVS_INDELS] ${caller_id} Finished. Total SNVs: "
    bcftools view -H ${caller_id}.snvs.vcf.gz | wc -l
    echo -n "[SPLIT_SNVS_INDELS] ${caller_id} Finished. Total INDELs: "
    bcftools view -H ${caller_id}.indels.vcf.gz | wc -l
    """
}

process CONSENSUS_SNVS {
    tag "snv_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}", mode: "copy"

    input:
    path vcfs
    path tbis

    output:
    tuple path("snv_consensus.vcf.gz"), path("snv_consensus.vcf.gz.tbi"), emit: consensus
    path("snv_caller_support.tsv"), emit: support
    path("snv_support_histogram.tsv"), emit: histogram

    script:
    """
    set -euo pipefail
    echo "========================================================="
    echo "[CONSENSUS_SNVS] Starting Low-RAM SNV Consensus"
    echo "========================================================="
    
    ls -1 *.vcf.gz > snv_vcfs.list
    n=\$(wc -l < snv_vcfs.list)
    
    if [ "${params.snv_min_callers}" = "null" ] || [ -z "${params.snv_min_callers}" ]; then
        min_callers=\$((n-1))
    else
        min_callers=${params.snv_min_callers}
    fi
    echo "[CONSENSUS_SNVS] Total input VCF callers: \$n"
    echo "[CONSENSUS_SNVS] Minimum caller threshold (min_callers): \$min_callers"
    
    echo "[CONSENSUS_SNVS] Merging caller headers..."
    echo "[CONSENSUS_SNVS] Getting header from first VCF..."
    first_vcf=\$(ls -1 *.vcf.gz | head -n 1)
    gzip -dc "\$first_vcf" | awk '/^#/ {print} !/^#/ {exit}' > consensus.header
        
    > all_sites.tsv
    
    echo "[CONSENSUS_SNVS] Flattening variant files..."
    for f in *.vcf.gz; do
        caller="\${f%.snvs.vcf.gz}"
        bcftools view -H "\$f" | awk -v caller="\$caller" '
        BEGIN {OFS="\\t"}
        {
            key=\$1 ":" \$2 ":" \$4 ":" \$5
            print key, caller, \$0
        }' >> all_sites.tsv
    done

    echo "[CONSENSUS_SNVS] Running Pass 1 (Counting frequencies in RAM)..."
    awk -v m="\$min_callers" '
    BEGIN { OFS="\\t" }
    {
        key=\$1; caller=\$2
        count[key]++
        
        if (callers[key] == "") {
            callers[key] = caller
        } else if (index(callers[key], caller) == 0) {
            callers[key] = callers[key] "," caller
        }
    }
    END {
        print "variant_id", "support", "callers" > "snv_caller_support.tsv"
        for (key in count) {
            print count[key] > "raw_counts.txt"
            if (count[key] >= m) {
                print key, count[key], callers[key] > "snv_caller_support.tsv"
                print key > "passing_keys.txt" # Offload target keys to disk
            }
        }
    }' all_sites.tsv

    echo "[CONSENSUS_SNVS] Running Pass 2 (Extracting passing rows from disk)..."
    touch passing_keys.txt
    awk '
    BEGIN { OFS="\\t" }
    NR==FNR { passing[\$1]=1; next } # Load just the verified string keys
    {
        key=\$1
        if (key in passing && !(key in seen)) {
            seen[key]=1
            # Reconstruct and output the VCF text line directly to file stream
            vcf_row=\$3
            for(i=4; i<=NF; i++) vcf_row = vcf_row OFS \$i
            print vcf_row > "consensus.body"
        }
    }' passing_keys.txt all_sites.tsv

    echo "[CONSENSUS_SNVS] Compiling metrics..."
    if [ -f raw_counts.txt ]; then
        sort -n raw_counts.txt | uniq -c | awk 'BEGIN{OFS="\\t"} {print \$2,\$1}' > snv_support_histogram.tsv
    else
        touch snv_support_histogram.tsv
    fi
    
    touch consensus.body
    cat consensus.header consensus.body > unsorted_consensus.vcf
    
    bcftools sort unsorted_consensus.vcf -Oz -o snv_consensus.vcf.gz
    bcftools index -f -t snv_consensus.vcf.gz
    
    echo -n "[CONSENSUS_SNVS] Done! Total valid consensus SNVs: "
    bcftools view -H snv_consensus.vcf.gz | wc -l
    """
}

// ==============================================================================
// PROCESS 4: CONSENSUS_INDELS
// Mirror implementation of the SNV processing block tailored for insertion/deletion signatures.
// ==============================================================================
process CONSENSUS_INDELS {
    tag "indel_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}", mode: "copy"

    input:
    path vcfs
    path tbis

    output:
    tuple path("indel_consensus.vcf.gz"), path("indel_consensus.vcf.gz.tbi"), emit: consensus
    path("indel_caller_support.tsv"), emit: support
    path("indel_support_histogram.tsv"), emit: histogram

    script:
    """
    set -euo pipefail
    echo "========================================================="
    echo "[CONSENSUS_INDELS] Starting Low-RAM INDEL Consensus"
    echo "========================================================="
    
    ls -1 *.vcf.gz > indel_vcfs.list
    min_callers=${params.indel_min_callers}
    echo "[CONSENSUS_INDELS] Minimum caller threshold (min_callers): \$min_callers"
    
    echo "[CONSENSUS_INDELS] Merging caller headers..."
    first_vcf=\$(ls -1 *.vcf.gz | head -n 1)
    gzip -dc "\$first_vcf" | awk '/^#/ {print} !/^#/ {exit}' > consensus.header
    
    
    > all_sites.tsv
    
    echo "[CONSENSUS_INDELS] Flattening variant files..."
    for f in *.vcf.gz; do
        caller="\${f%.indels.vcf.gz}"
        bcftools view -H "\$f" | awk -v caller="\$caller" '
        BEGIN {OFS="\\t"}
        {
            key=\$1 ":" \$2 ":" \$4 ":" \$5
            print key, caller, \$0
        }' >> all_sites.tsv
    done

    echo "[CONSENSUS_INDELS] Running Pass 1 (Counting frequencies in RAM)..."
    awk -v m="\$min_callers" '
    BEGIN { OFS="\\t" }
    {
        key=\$1; caller=\$2
        count[key]++
        
        if (callers[key] == "") {
            callers[key] = caller
        } else if (index(callers[key], caller) == 0) {
            callers[key] = callers[key] "," caller
        }
    }
    END {
        print "variant_id", "support", "callers" > "indel_caller_support.tsv"
        for (key in count) {
            print count[key] > "raw_counts.txt"
            if (count[key] >= m) {
                print key, count[key], callers[key] > "indel_caller_support.tsv"
                print key > "passing_keys.txt"
            }
        }
    }' all_sites.tsv

    echo "[CONSENSUS_INDELS] Running Pass 2 (Extracting passing rows from disk)..."
    touch passing_keys.txt
    awk '
    BEGIN { OFS="\\t" }
    NR==FNR { passing[\$1]=1; next }
    {
        key=\$1
        if (key in passing && !(key in seen)) {
            seen[key]=1
            vcf_row=\$3
            for(i=4; i<=NF; i++) vcf_row = vcf_row OFS \$i
            print vcf_row > "consensus.body"
        }
    }' passing_keys.txt all_sites.tsv

    echo "[CONSENSUS_INDELS] Compiling metrics..."
    if [ -f raw_counts.txt ]; then
        sort -n raw_counts.txt | uniq -c | awk 'BEGIN{OFS="\\t"} {print \$2,\$1}' > indel_support_histogram.tsv
    else
        touch indel_support_histogram.tsv
    fi
    
    touch consensus.body
    cat consensus.header consensus.body > unsorted_consensus.vcf
    
    bcftools sort unsorted_consensus.vcf -Oz -o indel_consensus.vcf.gz
    bcftools index -f -t indel_consensus.vcf.gz
    
    echo -n "[CONSENSUS_INDELS] Done! Total valid consensus INDELs: "
    bcftools view -H indel_consensus.vcf.gz | wc -l
    """
}

// ==============================================================================
// PROCESS 5: MERGE_CONSENSUS
// Concatenates the validated consensus SNV and INDEL datasets back into a final VCF.
// ==============================================================================
process MERGE_CONSENSUS {
    tag "merge_consensus"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple path(snv_vcf), path(snv_tbi)
    tuple path(indel_vcf), path(indel_tbi)

    output:
    tuple path("final_consensus.vcf.gz"), path("final_consensus.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    echo "[MERGE_CONSENSUS] Concatenating SNVs and INDELs into unified output file..."
    bcftools concat -a \
      ${snv_vcf} \
      ${indel_vcf} \
      -Oz \
      -o final_consensus.vcf.gz
    bcftools index -f -t final_consensus.vcf.gz
    
    echo -n "[MERGE_CONSENSUS] Workflow finished successfully! Total combined variants: "
    bcftools view -H final_consensus.vcf.gz | wc -l
    """
}

// ==============================================================================
// PIPELINE WORKFLOW ENTRYPOINT
// Evaluates variable configurations, constructs channels, and schedules jobs.
// ==============================================================================
workflow {
    log.info "========================================="
    log.info "       Somatic Consensus Pipeline        "
    log.info "========================================="
    log.info "Outputs Directory : ${params.outdir}"
    log.info "QC Outputs        : ${params.qc_outdir}"
    log.info "Reference FASTA   : ${params.ref_fasta}"
    log.info "-----------------------------------------"

    // Group up arguments and filter out unassigned elements
    vcfs = [
        params.vcf1, params.vcf2, params.vcf3,
        params.vcf4, params.vcf5, params.vcf6, params.vcf7
    ].findAll { it }
    
    if (vcfs.size() < 3) {
        error "Fatal Pipeline Error: Need at least 3 VCFs for consensus filtering routines."
    }
    
    log.info "[WORKFLOW] Detected ${vcfs.size()} input file paths for consolidation."
    
    // Transform file path assignments into (Caller_ID, File) target configurations
    vcf_ch = Channel
        .fromList(vcfs)
        .map { vcf ->
            def base = vcf.tokenize('/')[-1]
            def caller_id = base.replaceAll(/\.vcf\.gz$/, '')
            log.info "[WORKFLOW] Queueing Caller ID -> [ ${caller_id} ]"
            tuple(caller_id, file(vcf))
        }
        
    log.info "[WORKFLOW] Phase 1: Checking indices..."
    indexed_ch = ENSURE_INDEX(vcf_ch)
    
    log.info "[WORKFLOW] Phase 2: Splitting caller profiles into distinct SNV and INDEL vectors..."
    split_ch = SPLIT_SNVS_INDELS(
        indexed_ch,
        file(params.ref_fasta),
        file(params.ref_fai)
    )
    
    // Consolidate output properties across downstream tasks
    snv_vcfs = split_ch.snvs.map { caller_id, vcf, tbi -> vcf }.collect()
    snv_tbis = split_ch.snvs.map { caller_id, vcf, tbi -> tbi }.collect()
    
    indel_vcfs = split_ch.indels.map { caller_id, vcf, tbi -> vcf }.collect()
    indel_tbis = split_ch.indels.map { caller_id, vcf, tbi -> tbi }.collect()
    
    log.info "[WORKFLOW] Phase 3: Executing multi-caller matrix parsing routines..."
    snv_consensus = CONSENSUS_SNVS(snv_vcfs, snv_tbis)
    indel_consensus = CONSENSUS_INDELS(indel_vcfs, indel_tbis)
    
    log.info "[WORKFLOW] Phase 4: Merging final consensus files..."
    MERGE_CONSENSUS(snv_consensus.consensus, indel_consensus.consensus)
}

// nextflow.enable.dsl=2

// params.vcf1 = params.vcf1 ?: null
// params.vcf2 = params.vcf2 ?: null
// params.vcf3 = params.vcf3 ?: null
// params.vcf4 = params.vcf4 ?: null
// params.vcf5 = params.vcf5 ?: null
// params.vcf6 = params.vcf6 ?: null
// params.vcf7 = params.vcf7 ?: null

// params.outdir = params.outdir ?: "results"
// params.qc_outdir = params.qc_outdir ?: "${params.outdir}/qc"


// process CONSENSUS_INDELS {
//     tag "indel_consensus"
//     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

//     input:
//     path vcfs
//     path tbis

//     output:
//     tuple path("indel_consensus.vcf.gz"), path("indel_consensus.vcf.gz.tbi"), emit: consensus
//     tuple path("merged_indels.vcf.gz"), path("merged_indels.vcf.gz.tbi"), emit: merged
//     path("indel_support_histogram.tsv"), emit: histogram

//     script:
//     """
//     set -euo pipefail

//     echo "=== INPUT FILES ==="
//     ls -lh *.vcf.gz

//     echo "=== VARIANT COUNTS PER INPUT ==="
//     for f in *.vcf.gz; do
//         echo -n "\$f: "
//         bcftools view -H "\$f" | wc -l
//     done

//     ls -1 *.vcf.gz > indel_vcfs.list

//     echo "INDEL min_callers=${params.indel_min_callers}"

//     bcftools merge \
//       --force-samples \
//       --file-list indel_vcfs.list \
//       -m none \
//       -Oz \
//       -o merged_indels.vcf.gz

//     bcftools index -f -t merged_indels.vcf.gz

//     echo "=== ALT-support histogram from merged INDELs; no bcftools query ==="
//     bcftools view -H merged_indels.vcf.gz | \
//     awk '
//     {
//         c=0
//         for(i=10;i<=NF;i++) {
//             split(\$i, fields, ":")
//             gt=fields[1]
//             if (gt ~ /[1-9]/) c++
//         }
//         print c
//     }' | sort -n | uniq -c | \
//     awk '{print \$2 "\\t" \$1}' > indel_support_histogram.tsv

//     echo "=== BUILD INDEL CONSENSUS: require ALT support >= ${params.indel_min_callers} ==="
//     bcftools view \
//       -i "N_PASS(GT~\"[1-9]\")>=${params.indel_min_callers}" \
//       -Oz \
//       -o indel_consensus.vcf.gz \
//       merged_indels.vcf.gz

//     bcftools index -f -t indel_consensus.vcf.gz

//     echo -n "INDEL consensus count: "
//     bcftools view -H indel_consensus.vcf.gz | wc -l

//     echo "=== CONSENSUS VALIDATION (INDELS) ==="
//     expected_consensus_count=\$(awk -v m="${params.indel_min_callers}" '\$1 >= m {s += \$2} END {print s+0}' indel_support_histogram.tsv)
//     actual_consensus_count=\$(bcftools view -H indel_consensus.vcf.gz | wc -l)

//     echo "Expected from histogram: \$expected_consensus_count"
//     echo "Actual VCF count:        \$actual_consensus_count"
//     """
// }

// process CONSENSUS_SNVS {
//     tag "snv_consensus"
//     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

//     input:
//     path vcfs
//     path tbis

//     output:
//     tuple path("snv_consensus.vcf.gz"), path("snv_consensus.vcf.gz.tbi"), emit: consensus
//     tuple path("merged_snvs.vcf.gz"), path("merged_snvs.vcf.gz.tbi"), emit: merged
//     path("snv_support_histogram.tsv"), emit: histogram

//     script:
//     """
//     set -euo pipefail

//     echo "=== INPUT FILES ==="
//     ls -lh *.vcf.gz

//     echo "=== VARIANT COUNTS PER INPUT ==="
//     for f in *.vcf.gz; do
//         echo -n "\$f: "
//         bcftools view -H "\$f" | wc -l
//     done

//     ls -1 *.vcf.gz > snv_vcfs.list
//     n=\$(wc -l < snv_vcfs.list)

//     if [ "${params.snv_min_callers}" = "null" ] || [ -z "${params.snv_min_callers}" ]; then
//         min_callers=\$((n-1))
//     else
//         min_callers=${params.snv_min_callers}
//     fi

//     echo "SNV min_callers=\$min_callers"

//     bcftools merge \
//       --force-samples \
//       --file-list snv_vcfs.list \
//       -m none \
//       -Oz \
//       -o merged_snvs.vcf.gz

//     bcftools index -f -t merged_snvs.vcf.gz

//     echo "=== ALT-support histogram from merged SNVs; no bcftools query ==="
//     bcftools view -H merged_snvs.vcf.gz | \
//     awk '
//     {
//         c=0
//         for(i=10;i<=NF;i++) {
//             split(\$i, fields, ":")
//             gt=fields[1]
//             if (gt ~ /[1-9]/) c++
//         }
//         print c
//     }' | sort -n | uniq -c | \
//     awk '{print \$2 "\\t" \$1}' > snv_support_histogram.tsv

//     echo "=== BUILD SNV CONSENSUS: require ALT support >= \$min_callers ==="
//     bcftools view \
//       -i "N_PASS(GT~\"[1-9]\")>=\$min_callers" \
//       -Oz \
//       -o snv_consensus.vcf.gz \
//       merged_snvs.vcf.gz

//     bcftools index -f -t snv_consensus.vcf.gz

//     echo -n "SNV consensus count: "
//     bcftools view -H snv_consensus.vcf.gz | wc -l

//     echo "=== CONSENSUS VALIDATION (SNVS) ==="
//     expected_consensus_count=\$(awk -v m="\$min_callers" '\$1 >= m {s += \$2} END {print s+0}' snv_support_histogram.tsv)
//     actual_consensus_count=\$(bcftools view -H snv_consensus.vcf.gz | wc -l)

//     echo "Expected from histogram: \$expected_consensus_count"
//     echo "Actual VCF count:        \$actual_consensus_count"
//     """
// }

// process LIGHT_FILTER {
//     tag "${caller_id}"
//     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

//     input:
//     tuple val(caller_id), path(vcf), path(tbi)

//     output:
//     tuple val(caller_id), path("${caller_id}.filtered.vcf.gz"), path("${caller_id}.filtered.vcf.gz.tbi")

//     script:
//     """
//     set -euo pipefail

//     bcftools view \
//       -f 'PASS,.' \
//       -i 'GT!="mis" && GT!="ref" && FORMAT/AD[1] > 3 && FORMAT/AF[0] > 0.02' \
//       -Oz \
//       -o ${caller_id}.filtered.vcf.gz \
//       ${vcf}

//     bcftools index -t ${caller_id}.filtered.vcf.gz

//     echo -n "${caller_id} filtered count: "
//     bcftools view -H ${caller_id}.filtered.vcf.gz | wc -l
//     """
// }


// process MERGE_CONSENSUS {
//     tag "merge_consensus"
//     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
//     input:
//     tuple path(snv_vcf), path(snv_tbi)
//     tuple path(indel_vcf), path(indel_tbi)

//     output:
//     tuple path("final_consensus.vcf.gz"), path("final_consensus.vcf.gz.tbi")

//     script:
//     """
//     bcftools concat -a \
//       ${snv_vcf} \
//       ${indel_vcf} \
//       -Oz -o final_consensus.vcf.gz

//     bcftools index -t final_consensus.vcf.gz
//     """
// }


// process SPLIT_SNVS_INDELS {
//     tag "${caller_id}"
//     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

//     input:
//     tuple val(caller_id), path(vcf), path(tbi)
//     path ref_fasta
//     path ref_fai

//     output:
//     tuple val(caller_id), path("${caller_id}.snvs.vcf.gz"), path("${caller_id}.snvs.vcf.gz.tbi"), emit: snvs
//     tuple val(caller_id), path("${caller_id}.indels.vcf.gz"), path("${caller_id}.indels.vcf.gz.tbi"), emit: indels

//     script:
//     """
//     echo "=== ${caller_id}: DETERMINE TUMOR SAMPLE ==="

//     bcftools query -l ${vcf} > samples.txt
//     echo "Samples:"
//     cat samples.txt

//     if echo "${caller_id}" | grep -qi "strelka"; then
//         selected_sample="TUMOR"
//     else
//         selected_sample=\$(grep -viE 'PBMC|NORMAL|BLOOD' samples.txt | head -n 1)
//     fi

//     if [ -z "\$selected_sample" ]; then
//         echo "ERROR: could not identify tumor sample for ${caller_id}" >&2
//         cat samples.txt >&2
//         exit 1
//     fi

//     if ! grep -Fxq "\$selected_sample" samples.txt; then
//         echo "ERROR: selected sample \$selected_sample not found in ${caller_id}" >&2
//         cat samples.txt >&2
//         exit 1
//     fi

//     echo "Selected tumor sample: \$selected_sample"

//     bcftools view -s "\$selected_sample" ${vcf} -Oz -o ${caller_id}.tumor_only.vcf.gz
//     bcftools index -t ${caller_id}.tumor_only.vcf.gz

//     echo "=== QC ${caller_id} ==="
//     echo -n "RAW total: "
//     bcftools view -H ${vcf} | wc -l

//     echo "=== NORMALIZING ${caller_id} ==="
//     bcftools norm -f ${ref_fasta} -m -any ${caller_id}.tumor_only.vcf.gz -Oz -o ${caller_id}.norm.vcf.gz
//     bcftools index -t ${caller_id}.norm.vcf.gz

//     echo -n "NORMALIZED total: "
//     bcftools view -H ${caller_id}.norm.vcf.gz | wc -l

//     echo "=== SPLITTING SNVS ==="
//     bcftools view -v snps ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.snvs.vcf.gz
//     bcftools index -t ${caller_id}.snvs.vcf.gz

//     echo "=== SPLITTING INDELS ==="
//     bcftools view -v indels ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.indels.vcf.gz
//     bcftools index -t ${caller_id}.indels.vcf.gz

//     echo -n "SNVs: "
//     bcftools view -H ${caller_id}.snvs.vcf.gz | wc -l

//     echo -n "INDELs: "
//     bcftools view -H ${caller_id}.indels.vcf.gz | wc -l

//     echo -n "PASS (SNVs): "
//     bcftools view -H -f PASS,. ${caller_id}.snvs.vcf.gz | wc -l

//     echo -n "PASS (INDELs): "
//     bcftools view -H -f PASS,. ${caller_id}.indels.vcf.gz | wc -l
//     """
// }


// //TODO: you could make it such that you input the tbis as well as the file, but maybe this is somethign we can solve later instead of creating the tbi everyt ime
// process ENSURE_INDEX {
//     tag "${caller_id}"
//     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

//     input:
//     tuple val(caller_id), path(vcf)

//     output:
//     tuple val(caller_id), path(vcf), path("${vcf}.tbi")

//     script:
//     """
//     bcftools index -t ${vcf}
//     """
// }



// process MAKE_CONSENSUS_QC_TSVS {
//     tag "${kind}_qc_tsvs"
//     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
//     publishDir "${params.qc_outdir}", mode: 'copy'

//     input:
//     tuple val(kind), path(merged_vcf), path(merged_tbi), path(consensus_vcf), path(consensus_tbi)

//     output:
//     tuple val(kind),
//           path("${kind}_caller_overlap.tsv"),
//           path("${kind}_caller_overlap_counts.tsv"),
//           path("${kind}_caller_counts.tsv"),
//           path("${kind}_consensus_af_depth.tsv"),
//           path("${kind}_chrom_counts.tsv"),
//           emit: qc_tsvs

//     script:
//     """
//     set -euo pipefail

//     echo "=== QC TSVS for ${kind} ==="

//     bcftools view -h ${merged_vcf} \
//       | awk '/^#CHROM/ { for (i=10; i<=NF; i++) print \$i }' \
//       > callers.txt

//     {
//         echo -e "variant_id\tchrom\tpos\tref\talt\tsupport\tcallers"
//         bcftools view -H ${merged_vcf} | awk '
//         BEGIN {
//             while ((getline line < "callers.txt") > 0) caller[++n] = line
//             OFS="\t"
//         }
//         {
//             chrom=\$1; pos=\$2; ref=\$4; alt=\$5
//             variant_id = chrom ":" pos ":" ref ":" alt
//             support = 0; callers = ""
//             for (i=10; i<=NF; i++) {
//                 split(\$i, fields, ":")
//                 gt = fields[1]
//                 if (gt ~ /[1-9]/) {
//                     support++
//                     cidx = i - 9
//                     callers = callers == "" ? caller[cidx] : callers "," caller[cidx]
//                 }
//             }
//             if (callers == "") callers = "none"
//             print variant_id, chrom, pos, ref, alt, support, callers
//         }'
//     } > ${kind}_caller_overlap.tsv

//     tail -n +2 ${kind}_caller_overlap.tsv \
//       | cut -f7 \
//       | sort \
//       | uniq -c \
//       | sort -nr \
//       | awk 'BEGIN{OFS="\t"; print "caller_combo","count"} {count=\$1; \$1=""; sub(/^ /,""); print \$0,count}' \
//       > ${kind}_caller_overlap_counts.tsv

//     {
//         echo -e "caller\talt_variant_count"
//         awk '{print NR "\t" \$0}' callers.txt | while read idx caller; do
//             sample_col=\$((idx + 9))
//             count=\$(bcftools view -H ${merged_vcf} | awk -v sample_col="\$sample_col" '{split(\$sample_col, fields, ":"); gt=fields[1]; if (gt ~ /[1-9]/) c++} END {print c+0}')
//             echo -e "\${caller}\t\${count}"
//         done
//     } > ${kind}_caller_counts.tsv

//     bcftools view -H ${consensus_vcf} \
//       | awk 'BEGIN{OFS="\t"; print "chrom","count"} {c[\$1]++} END{for (chrom in c) print chrom,c[chrom]}' \
//       | sort -V \
//       > ${kind}_chrom_counts.tsv

//     bcftools view -H ${consensus_vcf} | awk '
//     BEGIN {
//         OFS="\t"
//         print "variant_id","chrom","pos","ref","alt","caller","gt","af","dp","ad_ref","ad_alt","raw_sample","format"
//     }
//     {
//         chrom=\$1; pos=\$2; ref=\$4; alt=\$5; fmt=\$9
//         variant_id=chrom ":" pos ":" ref ":" alt
//         split(fmt, keys, ":")
//         delete keyidx
//         for (k=1; k<=length(keys); k++) keyidx[keys[k]]=k
//         for (i=10; i<=NF; i++) {
//             caller_idx=i-9
//             caller="caller_" caller_idx
//             split(\$i, vals, ":")
//             gt=("GT" in keyidx ? vals[keyidx["GT"]] : ".")
//             af=("AF" in keyidx ? vals[keyidx["AF"]] : ".")
//             dp=("DP" in keyidx ? vals[keyidx["DP"]] : ".")
//             ad_ref="."; ad_alt="."
//             if ("AD" in keyidx) {
//                 split(vals[keyidx["AD"]], ads, ",")
//                 ad_ref=ads[1]; ad_alt=ads[2]
//             }
//             if (gt ~ /[1-9]/) print variant_id,chrom,pos,ref,alt,caller,gt,af,dp,ad_ref,ad_alt,\$i,fmt
//         }
//     }' > ${kind}_consensus_af_depth.tsv

//     ls -lh ${kind}_*.tsv
//     """
// }

// process PLOT_CONSENSUS_QC {
//     tag "${kind}_qc_plots"
//     container 'python:3.11-slim'
//     publishDir "${params.qc_outdir}", mode: 'copy'

//     input:
//     tuple val(kind),
//           path(caller_overlap),
//           path(caller_overlap_counts),
//           path(caller_counts),
//           path(af_depth),
//           path(chrom_counts)

//     output:
//     path("${kind}_caller_overlap_barplot.png"), emit: caller_overlap_plot
//     path("${kind}_caller_counts_barplot.png"), emit: caller_counts_plot
//     path("${kind}_consensus_af_distribution.png"), emit: af_plot
//     path("${kind}_consensus_depth_distribution.png"), emit: depth_plot
//     path("${kind}_chrom_counts_barplot.png"), emit: chrom_plot

//     script:
//     """
//     set -euo pipefail

//     python - <<'PYSCRIPT'
// import csv, math
// import matplotlib.pyplot as plt
// kind = "${kind}"

// def read_tsv(path):
//     with open(path, newline='') as f:
//         return list(csv.DictReader(f, delimiter='\t'))

// def safe_float(x):
//     if x is None: return None
//     x = str(x).strip()
//     if x in {"", ".", "NA", "nan"}: return None
//     if "," in x: x = x.split(",")[0]
//     try:
//         v = float(x)
//         return None if math.isnan(v) else v
//     except Exception:
//         return None

// def save_bar(rows, label_col, value_col, outfile, title, xlabel, ylabel, top_n=None):
//     vals=[]
//     for r in rows:
//         label=str(r.get(label_col,""))
//         try: value=int(float(r.get(value_col,0)))
//         except Exception: value=0
//         vals.append((label,value))
//     vals=sorted(vals, key=lambda x:x[1], reverse=True)
//     if top_n: vals=vals[:top_n]
//     if not vals: vals=[("none",0)]
//     labels=[x[0] for x in vals]
//     values=[x[1] for x in vals]
//     plt.figure(figsize=(10, max(4, 0.35*len(labels)+1.5)))
//     plt.barh(labels[::-1], values[::-1])
//     plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title(title)
//     plt.tight_layout(); plt.savefig(outfile, dpi=200); plt.close()

// save_bar(read_tsv("${caller_overlap_counts}"), "caller_combo", "count", f"{kind}_caller_overlap_barplot.png", f"{kind.upper()} caller overlap combinations", "variant count", "caller combo", top_n=30)
// save_bar(read_tsv("${caller_counts}"), "caller", "alt_variant_count", f"{kind}_caller_counts_barplot.png", f"{kind.upper()} pre-consensus ALT counts per caller", "ALT variant count", "caller")
// save_bar(read_tsv("${chrom_counts}"), "chrom", "count", f"{kind}_chrom_counts_barplot.png", f"{kind.upper()} consensus variants per chromosome", "consensus variant count", "chromosome")

// rows=read_tsv("${af_depth}")
// afs=[safe_float(r.get("af")) for r in rows]
// afs=[x for x in afs if x is not None]
// depths=[]
// for r in rows:
//     dp=safe_float(r.get("dp"))
//     if dp is None:
//         ref=safe_float(r.get("ad_ref")); alt=safe_float(r.get("ad_alt"))
//         if ref is not None and alt is not None: dp=ref+alt
//     if dp is not None: depths.append(dp)

// plt.figure(figsize=(8,5))
// if afs: plt.hist(afs, bins=50)
// else: plt.text(0.5,0.5,"No parseable AF values found",ha="center",va="center")
// plt.xlabel("AF"); plt.ylabel("number of ALT-supporting caller genotypes"); plt.title(f"{kind.upper()} consensus AF distribution")
// plt.tight_layout(); plt.savefig(f"{kind}_consensus_af_distribution.png", dpi=200); plt.close()

// plt.figure(figsize=(8,5))
// if depths: plt.hist(depths, bins=50)
// else: plt.text(0.5,0.5,"No parseable DP/AD depth values found",ha="center",va="center")
// plt.xlabel("depth"); plt.ylabel("number of ALT-supporting caller genotypes"); plt.title(f"{kind.upper()} consensus depth distribution")
// plt.tight_layout(); plt.savefig(f"{kind}_consensus_depth_distribution.png", dpi=200); plt.close()
// PYSCRIPT

//     ls -lh ${kind}_*.png
//     """
// }

// workflow {

//     vcfs = [
//         params.vcf1, params.vcf2, params.vcf3,
//         params.vcf4, params.vcf5, params.vcf6, params.vcf7
//     ].findAll { it }

//     if (vcfs.size() < 3) {
//         error "Need at least 3 VCFs for consensus"
//     }

//     println "=== RAW PARAM INPUTS ==="
//     println vcfs
//     println "Total VCFs: ${vcfs.size()}"

//     vcf_ch = Channel
//     .fromList(vcfs)
//     .map { vcf ->
//         def base = vcf.tokenize('/')[-1]              // filename
//         def caller_id = base.replaceAll(/\.vcf\.gz$/, '')  // strip suffix

//         tuple(caller_id, file(vcf))
//     }

//     vcf_ch.view { "VCF_CH → caller_id=${it[0]}, file=${it[1]}" }

//     indexed_ch = ENSURE_INDEX(vcf_ch)

//     //dont light filter for now
//     filtered_ch =indexed_ch

//     split_ch = SPLIT_SNVS_INDELS(
//         filtered_ch,
//         file(params.ref_fasta),
//         file(params.ref_fai)
//     )

//     snv_vcfs = split_ch.snvs.map { caller_id, vcf, tbi -> vcf }.collect()
//     snv_tbis = split_ch.snvs.map { caller_id, vcf, tbi -> tbi }.collect()

//     indel_vcfs = split_ch.indels.map { caller_id, vcf, tbi -> vcf }.collect()
//     indel_tbis = split_ch.indels.map { caller_id, vcf, tbi -> tbi }.collect()

//     snv_consensus = CONSENSUS_SNVS(snv_vcfs, snv_tbis)
//     indel_consensus = CONSENSUS_INDELS(indel_vcfs, indel_tbis)

//     snv_qc_input = snv_consensus.merged
//         .combine(snv_consensus.consensus)
//         .map { merged_tuple, consensus_tuple ->
//             tuple("snv", merged_tuple[0], merged_tuple[1], consensus_tuple[0], consensus_tuple[1])
//         }

//     indel_qc_input = indel_consensus.merged
//         .combine(indel_consensus.consensus)
//         .map { merged_tuple, consensus_tuple ->
//             tuple("indel", merged_tuple[0], merged_tuple[1], consensus_tuple[0], consensus_tuple[1])
//         }

//     consensus_qc_input = snv_qc_input.mix(indel_qc_input)

//     consensus_qc_tsvs = MAKE_CONSENSUS_QC_TSVS(consensus_qc_input)

//     PLOT_CONSENSUS_QC(consensus_qc_tsvs.qc_tsvs)

//     MERGE_CONSENSUS(snv_consensus.consensus, indel_consensus.consensus)
// }


// // nextflow.enable.dsl=2

// // params.vcf1 = params.vcf1 ?: null
// // params.vcf2 = params.vcf2 ?: null
// // params.vcf3 = params.vcf3 ?: null
// // params.vcf4 = params.vcf4 ?: null
// // params.vcf5 = params.vcf5 ?: null
// // params.vcf6 = params.vcf6 ?: null
// // params.vcf7 = params.vcf7 ?: null


// // process CONSENSUS_INDELS {
// //     tag "indel_consensus"
// //     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

// //     input:
// //     path vcfs
// //     path tbis

// //     output:
// //     tuple path("indel_consensus.vcf.gz"), path("indel_consensus.vcf.gz.tbi"), emit: consensus
// //     path("indel_support_histogram.tsv"), emit: histogram

// //     script:
// //     """
// //     set -euo pipefail

// //     echo "=== INPUT FILES ==="
// //     ls -lh *.vcf.gz

// //     echo "=== VARIANT COUNTS PER INPUT ==="
// //     for f in *.vcf.gz; do
// //         echo -n "\$f: "
// //         bcftools view -H "\$f" | wc -l
// //     done

// //     ls -1 *.vcf.gz > indel_vcfs.list

// //     echo "INDEL min_callers=${params.indel_min_callers}"

// //     bcftools merge \\
// //       --force-samples \\
// //       --file-list indel_vcfs.list \\
// //       -m none \\
// //       -Oz \\
// //       -o merged_indels.vcf.gz

// //     bcftools index -f -t merged_indels.vcf.gz

// //     echo "=== ALT-support histogram from merged INDELs ==="

// //     bcftools query -f '[%GT\\t]\\n' merged_indels.vcf.gz | \\
// //     awk '
// //     {
// //         c=0
// //         for(i=1;i<=NF;i++) {
// //             gt=\$i
// //             if (gt ~ /[1-9]/) c++
// //         }
// //         print c
// //     }' | sort -n | uniq -c | \\
// //     awk '{print \$2 "\\t" \$1}' \\
// //     > indel_support_histogram.tsv

// //     echo "=== BUILD INDEL CONSENSUS: require ALT support >= ${params.indel_min_callers} ==="

// //     bcftools view \\
// //       -i "N_PASS(GT~\\"[1-9]\\")>=${params.indel_min_callers}" \\
// //       -Oz \\
// //       -o indel_consensus.vcf.gz \\
// //       merged_indels.vcf.gz

// //     bcftools index -f -t indel_consensus.vcf.gz

// //     echo -n "INDEL consensus count: "
// //     bcftools view -H indel_consensus.vcf.gz | wc -l

// //     echo "=== CONSENSUS VALIDATION (INDELS) ==="

// //     expected_consensus_count=\$(awk -v m="${params.indel_min_callers}" '
// //     \$1 >= m {s += \$2}
// //     END {print s+0}
// //     ' indel_support_histogram.tsv)

// //     actual_consensus_count=\$(bcftools view -H indel_consensus.vcf.gz | wc -l)

// //     echo "Expected from histogram: \$expected_consensus_count"
// //     echo "Actual VCF count:        \$actual_consensus_count"
// //     """
// // }

// // process CONSENSUS_SNVS {
// //     tag "snv_consensus"
// //     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

// //     input:
// //     path vcfs
// //     path tbis

// //     output:
// //     tuple path("snv_consensus.vcf.gz"), path("snv_consensus.vcf.gz.tbi"), emit: consensus
// //     path("snv_support_histogram.tsv"), emit: histogram

// //     script:
// //     """
// //     set -euo pipefail

// //     echo "=== INPUT FILES ==="
// //     ls -lh *.vcf.gz

// //     echo "=== VARIANT COUNTS PER INPUT ==="
// //     for f in *.vcf.gz; do
// //         echo -n "\$f: "
// //         bcftools view -H "\$f" | wc -l
// //     done

// //     ls -1 *.vcf.gz > snv_vcfs.list
// //     n=\$(wc -l < snv_vcfs.list)

// //     if [ "${params.snv_min_callers}" = "null" ] || [ -z "${params.snv_min_callers}" ]; then
// //         min_callers=\$((n-1))
// //     else
// //         min_callers=${params.snv_min_callers}
// //     fi

// //     echo "SNV min_callers=\$min_callers"

// //     bcftools merge \\
// //       --force-samples \\
// //       --file-list snv_vcfs.list \\
// //       -m none \\
// //       -Oz \\
// //       -o merged_snvs.vcf.gz

// //     bcftools index -f -t merged_snvs.vcf.gz

// //     echo "=== ALT-support histogram from merged SNVs ==="

// //     bcftools query -f '[%GT\\t]\\n' merged_snvs.vcf.gz | \\
// //     awk '
// //     {
// //         c=0
// //         for(i=1;i<=NF;i++) {
// //             gt=\$i
// //             if (gt ~ /[1-9]/) c++
// //         }
// //         print c
// //     }' | sort -n | uniq -c | \\
// //     awk '{print \$2 "\\t" \$1}' \\
// //     > snv_support_histogram.tsv

// //     echo "=== BUILD SNV CONSENSUS: require ALT support >= \$min_callers ==="

// //     bcftools view \\
// //       -i "N_PASS(GT~\\"[1-9]\\")>=\$min_callers" \\
// //       -Oz \\
// //       -o snv_consensus.vcf.gz \\
// //       merged_snvs.vcf.gz

// //     bcftools index -f -t snv_consensus.vcf.gz

// //     echo -n "SNV consensus count: "
// //     bcftools view -H snv_consensus.vcf.gz | wc -l

// //     echo "=== CONSENSUS VALIDATION (SNVS) ==="

// //     expected_consensus_count=\$(awk -v m="\$min_callers" '
// //     \$1 >= m {s += \$2}
// //     END {print s+0}
// //     ' snv_support_histogram.tsv)

// //     actual_consensus_count=\$(bcftools view -H snv_consensus.vcf.gz | wc -l)

// //     echo "Expected from histogram: \$expected_consensus_count"
// //     echo "Actual VCF count:        \$actual_consensus_count"
// //     """
// // }

// // process LIGHT_FILTER {
// //     tag "${caller_id}"
// //     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

// //     input:
// //     tuple val(caller_id), path(vcf), path(tbi)

// //     output:
// //     tuple val(caller_id), path("${caller_id}.filtered.vcf.gz"), path("${caller_id}.filtered.vcf.gz.tbi")

// //     script:
// //     """
// //     set -euo pipefail

// //     bcftools view \
// //       -f 'PASS,.' \
// //       -i 'GT!="mis" && GT!="ref" && FORMAT/AD[1] > 3 && FORMAT/AF[0] > 0.02' \
// //       -Oz \
// //       -o ${caller_id}.filtered.vcf.gz \
// //       ${vcf}

// //     bcftools index -t ${caller_id}.filtered.vcf.gz

// //     echo -n "${caller_id} filtered count: "
// //     bcftools view -H ${caller_id}.filtered.vcf.gz | wc -l
// //     """
// // }


// // process MERGE_CONSENSUS {
// //     tag "merge_consensus"
// //     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
// //     input:
// //     tuple path(snv_vcf), path(snv_tbi)
// //     tuple path(indel_vcf), path(indel_tbi)

// //     output:
// //     tuple path("final_consensus.vcf.gz"), path("final_consensus.vcf.gz.tbi")

// //     script:
// //     """
// //     bcftools concat -a \
// //       ${snv_vcf} \
// //       ${indel_vcf} \
// //       -Oz -o final_consensus.vcf.gz

// //     bcftools index -t final_consensus.vcf.gz
// //     """
// // }


// // process SPLIT_SNVS_INDELS {
// //     tag "${caller_id}"
// //     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

// //     input:
// //     tuple val(caller_id), path(vcf), path(tbi)
// //     path ref_fasta
// //     path ref_fai

// //     output:
// //     tuple val(caller_id), path("${caller_id}.snvs.vcf.gz"), path("${caller_id}.snvs.vcf.gz.tbi"), emit: snvs
// //     tuple val(caller_id), path("${caller_id}.indels.vcf.gz"), path("${caller_id}.indels.vcf.gz.tbi"), emit: indels

// //     script:
// //     """
// //     echo "=== ${caller_id}: DETERMINE TUMOR SAMPLE ==="

// //     bcftools query -l ${vcf} > samples.txt
// //     echo "Samples:"
// //     cat samples.txt

// //     if echo "${caller_id}" | grep -qi "strelka"; then
// //         selected_sample="TUMOR"
// //     else
// //         selected_sample=\$(grep -viE 'PBMC|NORMAL|BLOOD' samples.txt | head -n 1)
// //     fi

// //     if [ -z "\$selected_sample" ]; then
// //         echo "ERROR: could not identify tumor sample for ${caller_id}" >&2
// //         cat samples.txt >&2
// //         exit 1
// //     fi

// //     if ! grep -Fxq "\$selected_sample" samples.txt; then
// //         echo "ERROR: selected sample \$selected_sample not found in ${caller_id}" >&2
// //         cat samples.txt >&2
// //         exit 1
// //     fi

// //     echo "Selected tumor sample: \$selected_sample"

// //     bcftools view -s "\$selected_sample" ${vcf} -Oz -o ${caller_id}.tumor_only.vcf.gz
// //     bcftools index -t ${caller_id}.tumor_only.vcf.gz

// //     echo "=== QC ${caller_id} ==="
// //     echo -n "RAW total: "
// //     bcftools view -H ${vcf} | wc -l

// //     echo "=== NORMALIZING ${caller_id} ==="
// //     bcftools norm -f ${ref_fasta} -m -any ${caller_id}.tumor_only.vcf.gz -Oz -o ${caller_id}.norm.vcf.gz
// //     bcftools index -t ${caller_id}.norm.vcf.gz

// //     echo -n "NORMALIZED total: "
// //     bcftools view -H ${caller_id}.norm.vcf.gz | wc -l

// //     echo "=== SPLITTING SNVS ==="
// //     bcftools view -v snps ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.snvs.vcf.gz
// //     bcftools index -t ${caller_id}.snvs.vcf.gz

// //     echo "=== SPLITTING INDELS ==="
// //     bcftools view -v indels ${caller_id}.norm.vcf.gz -Oz -o ${caller_id}.indels.vcf.gz
// //     bcftools index -t ${caller_id}.indels.vcf.gz

// //     echo -n "SNVs: "
// //     bcftools view -H ${caller_id}.snvs.vcf.gz | wc -l

// //     echo -n "INDELs: "
// //     bcftools view -H ${caller_id}.indels.vcf.gz | wc -l

// //     echo -n "PASS (SNVs): "
// //     bcftools view -H -f PASS,. ${caller_id}.snvs.vcf.gz | wc -l

// //     echo -n "PASS (INDELs): "
// //     bcftools view -H -f PASS,. ${caller_id}.indels.vcf.gz | wc -l
// //     """
// // }


// // //TODO: you could make it such that you input the tbis as well as the file, but maybe this is somethign we can solve later instead of creating the tbi everyt ime
// // process ENSURE_INDEX {
// //     tag "${caller_id}"
// //     container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

// //     input:
// //     tuple val(caller_id), path(vcf)

// //     output:
// //     tuple val(caller_id), path(vcf), path("${vcf}.tbi")

// //     script:
// //     """
// //     bcftools index -t ${vcf}
// //     """
// // }


// // workflow {

// //     vcfs = [
// //         params.vcf1, params.vcf2, params.vcf3,
// //         params.vcf4, params.vcf5, params.vcf6, params.vcf7
// //     ].findAll { it }

// //     if (vcfs.size() < 3) {
// //         error "Need at least 3 VCFs for consensus"
// //     }

// //     println "=== RAW PARAM INPUTS ==="
// //     println vcfs
// //     println "Total VCFs: ${vcfs.size()}"

// //     vcf_ch = Channel
// //     .fromList(vcfs)
// //     .map { vcf ->
// //         def base = vcf.tokenize('/')[-1]              // filename
// //         def caller_id = base.replaceAll(/\.vcf\.gz$/, '')  // strip suffix

// //         tuple(caller_id, file(vcf))
// //     }

// //     vcf_ch.view { "VCF_CH → caller_id=${it[0]}, file=${it[1]}" }

// //     indexed_ch = ENSURE_INDEX(vcf_ch)

// //     //dont light filter for now
// //     filtered_ch =indexed_ch

// //     split_ch = SPLIT_SNVS_INDELS(
// //         filtered_ch,
// //         file(params.ref_fasta),
// //         file(params.ref_fai)
// //     )

// //     snv_vcfs = split_ch.snvs.map { caller_id, vcf, tbi -> vcf }.collect()
// //     snv_tbis = split_ch.snvs.map { caller_id, vcf, tbi -> tbi }.collect()

// //     indel_vcfs = split_ch.indels.map { caller_id, vcf, tbi -> vcf }.collect()
// //     indel_tbis = split_ch.indels.map { caller_id, vcf, tbi -> tbi }.collect()

// //     snv_consensus = CONSENSUS_SNVS(snv_vcfs, snv_tbis)
// //     indel_consensus = CONSENSUS_INDELS(indel_vcfs, indel_tbis)

// //     MERGE_CONSENSUS(snv_consensus.consensus, indel_consensus.consensus)
// // }

