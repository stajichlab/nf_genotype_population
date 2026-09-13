// modules/local/gatk4_haplotypecaller/main.nf
process GATK4_HAPLOTYPECALLER {
    tag "$strain"
    label 'process_medium'
    // Local project-provided container (GATK 4.5.0.0).
    // Source: docker://broadinstitute/gatk:4.5.0.0, pulled via
    // depot.galaxyproject.org/singularity/broadinstitute/gatk:4.5.0.0
    // Referenced here as a direct local .sif/.img file path (pre-pulled to the
    // shared cache) so no docker:// pull is ever attempted at pipeline runtime
    // (same rationale as CUSTOM_HET_PLOIDY, see Task 8's report).
    container '/bigdata/stajichlab/shared/lib/singularity_cache/depot.galaxyproject.org-singularity-broadinstitute-gatk-4.5.0.0.img'

    input:
    tuple val(strain), val(ploidy), path(cram), path(crai)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(strain), path("${strain}.g.vcf.gz"), path("${strain}.g.vcf.gz.tbi"), emit: gvcf

    script:
    // Verified (manual `singularity exec ... bash -lc 'command -v gatk'` and via
    // actual Nextflow-driven execution under this project's login-shell
    // process.shell) that this container's own /etc/profile does NOT reset PATH
    // in login-shell mode: `gatk` resolves to /gatk/gatk under both `bash -c` and
    // `bash -lc`, unlike the bcftools/samtools container (see CUSTOM_HET_PLOIDY).
    // No explicit PATH export is needed here.
    """
    gatk HaplotypeCaller \\
        -R ${reference} \\
        -I ${cram} \\
        -O ${strain}.g.vcf.gz \\
        --sample-ploidy ${ploidy} \\
        -ERC GVCF
    """
}
