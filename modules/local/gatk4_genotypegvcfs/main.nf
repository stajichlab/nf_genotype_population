// modules/local/gatk4_genotypegvcfs/main.nf
process GATK4_GENOTYPEGVCFS {
    tag "$population"
    label 'process_medium'
    // FINAL-REVIEW I2: container (GATK 4.5.0.0) declared centrally in
    // conf/modules.config (spec §5). See GATK4_GENOMICSDBIMPORT /
    // GATK4_HAPLOTYPECALLER for the pre-pull rationale.

    input:
    tuple val(population), path(genomicsdb)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(population), path("${population}.vcf.gz"), path("${population}.vcf.gz.tbi"), emit: vcf

    script:
    // FINAL-REVIEW I5: cap the JVM heap at ~80% of this task's granted memory
    // (see GATK4_HAPLOTYPECALLER for the full rationale).
    def xmx = (task.memory.toGiga() * 0.8) as int
    """
    set -euo pipefail
    gatk --java-options "-Xmx${xmx}g" GenotypeGVCFs \\
        -R "${reference}" \\
        -V "gendb://${genomicsdb}" \\
        -O "${population}.vcf.gz"
    """
}
