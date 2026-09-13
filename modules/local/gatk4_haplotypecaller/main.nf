// modules/local/gatk4_haplotypecaller/main.nf
process GATK4_HAPLOTYPECALLER {
    tag "$strain"
    label 'process_medium'
    // FINAL-REVIEW I2: the container (GATK 4.5.0.0) is now declared centrally in
    // conf/modules.config, per spec §5, together with its upstream source
    // (docker://broadinstitute/gatk:4.5.0.0) and the reason a pre-pulled local
    // image is used instead of a runtime docker:// pull.

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
    //
    // FINAL-REVIEW I5: the `gatk` launcher sizes the JVM heap from the NODE's
    // physical RAM, not from the SLURM allocation, so on a large shared node it
    // picks a heap far above this task's cgroup limit and gets OOM killed with
    // an opaque failure. Cap it at ~80% of the memory this task was actually
    // granted. Derived from task.memory rather than hardcoded so it tracks the
    // task.attempt retry escalation in nextflow.config.
    def xmx = (task.memory.toGiga() * 0.8) as int
    """
    set -euo pipefail
    gatk --java-options "-Xmx${xmx}g" HaplotypeCaller \\
        -R "${reference}" \\
        -I "${cram}" \\
        -O "${strain}.g.vcf.gz" \\
        --sample-ploidy ${ploidy} \\
        --native-pair-hmm-threads ${task.cpus} \\
        -ERC GVCF
    """
}
