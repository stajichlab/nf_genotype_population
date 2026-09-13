// modules/local/gatk4_genomicsdbimport/main.nf
process GATK4_GENOMICSDBIMPORT {
    tag "$population"
    label 'process_medium'
    // Local project-provided container (GATK 4.5.0.0).
    // Source: docker://broadinstitute/gatk:4.5.0.0, pulled via
    // depot.galaxyproject.org/singularity/broadinstitute/gatk:4.5.0.0
    // Referenced here as a direct local .sif/.img file path (pre-pulled to the
    // shared cache) so no docker:// pull is ever attempted at pipeline runtime
    // (same rationale as GATK4_HAPLOTYPECALLER, see Task 9's report).
    container '/bigdata/stajichlab/shared/lib/singularity_cache/depot.galaxyproject.org-singularity-broadinstitute-gatk-4.5.0.0.img'

    input:
    tuple val(population), path(gvcfs), path(tbis)
    path reference
    path reference_fai
    path reference_dict
    val intervals

    output:
    tuple val(population), path("${population}_gdb"), emit: genomicsdb

    script:
    def gvcf_args = gvcfs.collect { "-V ${it}" }.join(' ')
    """
    gatk GenomicsDBImport \\
        ${gvcf_args} \\
        --genomicsdb-workspace-path ${population}_gdb \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --bypass-feature-reader \\
        -L ${intervals}
    """
}
