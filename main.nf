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

  interval_res = split_intervals(
    file(params.ref_fasta,  checkIfExists: true),
    file(params.ref_fai,    checkIfExists: true),
    file(params.ref_dict,   checkIfExists: true),
    file(params.intervals,  checkIfExists: true),
    params.scatter_count as int
  )

  println "mutect_runs size = ${params.mutect_runs?.size()}"

  Channel.fromList(params.mutect_runs)
    .view { "RUN_RAW: ${it.output_prefix} :: ${it.tumor_reads}" }

  runs_ch = Channel.fromList(params.mutect_runs)
    .map { run ->
        tuple(
            [id: run.output_prefix],
            file(run.tumor_reads),
            file(run.tumor_reads_index),
            run.normal_reads       ?: null,
            run.normal_reads_index ?: null,
            run.tumor_sample_name
        )
    }

  runs_ch.view { "RUNS_CH: $it" }

  intervals_ready = interval_res.interval_shards.collect()

  subset_res = subset_tumor_per_shard(
    runs_ch.map { meta, tbam, tbai, nbam, nbai, tsample ->
        tuple(meta, tbam, tbai)
    },
    intervals_ready
  )

  subset_res.shard_bams.view { "SHARD_BAMS_RAW: $it" }
  subset_res.shard_bais.view { "SHARD_BAIS_RAW: $it" }
  subset_res.shard_intervals.view { "SHARD_INTERVALS_RAW: $it" }

  shard_bams_ch = subset_res.shard_bams
    .transpose()
    .map { sid, f -> tuple(sid, f.name.replaceFirst(/\.bam$/, ''), f) }

  shard_bais_ch = subset_res.shard_bais
    .transpose()
    .map { sid, f -> tuple(sid, f.name.replaceFirst(/\.bam\.bai$/, ''), f) }

  shard_intervals_ch = subset_res.shard_intervals
    .transpose()
    .map { sid, f -> tuple(sid, f.name.replaceFirst(/\.intervals$/, ''), f) }

  shard_bams_ch.view { "SHARD_BAMS: $it" }
  shard_bais_ch.view { "SHARD_BAIS: $it" }
  shard_intervals_ch.view { "SHARD_INTERVALS: $it" }

  mutect_inputs_ch = shard_bams_ch
    .join(shard_bais_ch,      by: [0, 1])
    .join(shard_intervals_ch, by: [0, 1])
    .map { sid, base, bam, bai, interval ->
      tuple(
        sid,
        interval, bam, bai,
        file(params.ref_fasta,         checkIfExists: true),
        file(params.ref_fai,           checkIfExists: true),
        file(params.ref_dict,          checkIfExists: true),
        file(params.germline_resource, checkIfExists: true)
      )
    }

  mutect_inputs_ch.view { "MUTECT_INPUT: $it" }

  normals_ch = runs_ch.map { meta, tbam, tbai, nbam, nbai, tsample ->
    def nbam_file = nbam ? file(nbam) : file(NO_NORMAL_BAM_PATH)
    def nbai_file = nbam ? file(nbai) : file(NO_NORMAL_BAI_PATH)
    tuple(meta.id, nbam_file, nbai_file, tsample)
  }

  mutect_inputs_ch
    .join(normals_ch, by: 0)
    .map { sid, interval, bam, bai, ref, fai, dict, germ, nbam, nbai, tsample ->
      tuple(
        tuple(interval, bam, bai, ref, fai, dict, germ),
        nbam, nbai,
        file(NO_ALLELES_VCF_PATH),
        file(NO_ALLELES_TBI_PATH),
        tsample,
        params.m2_extra_args ?: ''
      )
    }
    .multiMap { main_tuple, nbam, nbai, alleles, alleles_tbi, tsample, extra ->
        main:        main_tuple
        nbam:        nbam
        nbai:        nbai
        alleles:     alleles
        alleles_tbi: alleles_tbi
        tsample:     tsample
        extra:       extra
    }
    .set { mutect_split_ch }

  mutect_res = mutect_wrapper(
    mutect_split_ch.main,
    mutect_split_ch.nbam,
    mutect_split_ch.nbai,
    mutect_split_ch.alleles,
    mutect_split_ch.alleles_tbi,
    mutect_split_ch.tsample,
    mutect_split_ch.extra
  )

  mutect_res.vcf
    .groupTuple(size: params.scatter_count)
    .set { grouped_vcfs_ch }

  gather_vcfs(grouped_vcfs_ch)
}