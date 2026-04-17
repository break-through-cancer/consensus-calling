process CONSENSUS_INDELS {
    tag "indel_consensus"

    input:
    path vcfs
    path tbis

    output:
    tuple path("indel_consensus.vcf.gz"), path("indel_consensus.vcf.gz.tbi")

    script:
    """
    ls -1 *.vcf.gz > indel_vcfs.list

    bcftools merge --file-list indel_vcfs.list -m none -Oz -o merged_indels.vcf.gz
    bcftools index -t merged_indels.vcf.gz

    python ${projectDir}/bin/indel_consensus.py merged_indels.vcf.gz $params.indel_min_callers indel_consensus.vcf.gz

    bcftools index -t indel_consensus.vcf.gz
    """
}