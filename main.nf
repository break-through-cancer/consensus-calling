nextflow.enable.dsl=2

params.vcf1 = params.vcf1 ?: null
params.vcf2 = params.vcf2 ?: null
params.vcf3 = params.vcf3 ?: null
params.vcf4 = params.vcf4 ?: null
params.vcf5 = params.vcf5 ?: null
params.vcf6 = params.vcf6 ?: null
params.vcf7 = params.vcf7 ?: null

include { LIGHT_FILTER } from './modules/light_filter'
include { SPLIT_SNVS_INDELS } from './modules/split_snvs_indels'
include { CONSENSUS_SNVS } from './modules/consensus_snvs'
include { CONSENSUS_INDELS } from './modules/consensus_indels'
include { MERGE_CONSENSUS } from './modules/merge_consensus'

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

    vcf_ch = Channel
    .fromList(vcfs)
    .map { vcf ->
        tuple(
            vcf.tokenize('/')[-1].replace('.vcf.gz',''),
            file(vcf)
        )
    }

    indexed_ch = ENSURE_INDEX(vcf_ch)
    filtered_ch = LIGHT_FILTER(indexed_ch)

    split_ch = SPLIT_SNVS_INDELS(filtered_ch)

    snv_vcfs = split_ch.snvs.map { caller_id, vcf, tbi -> vcf }.collect()
    snv_tbis = split_ch.snvs.map { caller_id, vcf, tbi -> tbi }.collect()

    indel_vcfs = split_ch.indels.map { caller_id, vcf, tbi -> vcf }.collect()
    indel_tbis = split_ch.indels.map { caller_id, vcf, tbi -> tbi }.collect()

    snv_consensus = CONSENSUS_SNVS(snv_vcfs, snv_tbis)
    indel_consensus = CONSENSUS_INDELS(indel_vcfs, indel_tbis)

    MERGE_CONSENSUS(snv_consensus, indel_consensus)
}

