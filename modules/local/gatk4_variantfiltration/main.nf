// modules/local/gatk4_variantfiltration/main.nf
process GATK4_HARDFILTER {
    tag "$population"
    label 'process_low'
    // Local project-provided container (GATK 4.5.0.0).
    // Source: docker://broadinstitute/gatk:4.5.0.0, pulled via
    // depot.galaxyproject.org/singularity/broadinstitute/gatk:4.5.0.0
    // Referenced here as a direct local .sif/.img file path (pre-pulled to the
    // shared cache) so no docker:// pull is ever attempted at pipeline runtime
    // (same rationale as GATK4_HAPLOTYPECALLER, see Task 9's report).
    container '/bigdata/stajichlab/shared/lib/singularity_cache/depot.galaxyproject.org-singularity-broadinstitute-gatk-4.5.0.0.img'

    input:
    tuple val(population), path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(population), path("${population}.filtered.vcf.gz"), path("${population}.filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk SelectVariants -R ${reference} -V ${vcf} --select-type-to-include SNP -O ${population}.snp.vcf.gz
    gatk SelectVariants -R ${reference} -V ${vcf} --select-type-to-include INDEL -O ${population}.indel.vcf.gz

    gatk VariantFiltration -R ${reference} -V ${population}.snp.vcf.gz \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "FS > 60.0" --filter-name "FS60" \\
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \\
        -O ${population}.snp.filtered.vcf.gz

    gatk VariantFiltration -R ${reference} -V ${population}.indel.vcf.gz \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "FS > 200.0" --filter-name "FS200" \\
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        --filter-expression "SOR > 10.0" --filter-name "SOR10" \\
        -O ${population}.indel.filtered.vcf.gz

    gatk MergeVcfs \\
        -I ${population}.snp.filtered.vcf.gz \\
        -I ${population}.indel.filtered.vcf.gz \\
        -O ${population}.filtered.vcf.gz
    """
}
