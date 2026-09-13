// modules/local/gatk4_genotypegvcfs/main.nf
process GATK4_GENOTYPEGVCFS {
    tag "$population"
    label 'process_medium'
    // Local project-provided container (GATK 4.5.0.0). See
    // GATK4_GENOMICSDBIMPORT / GATK4_HAPLOTYPECALLER for the pre-pull rationale.
    container '/bigdata/stajichlab/shared/lib/singularity_cache/depot.galaxyproject.org-singularity-broadinstitute-gatk-4.5.0.0.img'

    input:
    tuple val(population), path(genomicsdb)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(population), path("${population}.vcf.gz"), path("${population}.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk GenotypeGVCFs \\
        -R ${reference} \\
        -V gendb://${genomicsdb} \\
        -O ${population}.vcf.gz
    """
}
