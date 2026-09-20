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
    // GATK4_GENOTYPEGVCFS: the crash ("Genotype has no likelihoods" in
    // ExcessHet) traces to the GenomicsDB *read* step, not the genotyping
    // step - confirmed by the "too many genotypes" warning still naming 50
    // alleles even after --max-alternate-alleles 24 was set (that flag only
    // trims alleles GATK actually genotypes/outputs, applied too late to
    // prevent the read-time combination). The controlling flag for the read
    // step is --genomicsdb-max-alternate-alleles (GATK 4.5.0.0 --help
    // default: 50), which is what combined 50 alleles into 1275 diploid
    // genotypes - over --max-genotype-count's default cap of 1024 - causing
    // GATK to drop that sample's PL field and then crash annotating it.
    // Capping both at 24 keeps genotype count for a diploid site at
    // 25*26/2=325, safely under the 1024 cap.
    """
    set -euo pipefail
    gatk --java-options "-Xmx${xmx}g" GenotypeGVCFs \\
        -R "${reference}" \\
        -V "gendb://${genomicsdb}" \\
        --genomicsdb-max-alternate-alleles 24 \\
        --max-alternate-alleles 24 \\
        -O "${population}.vcf.gz"
    """
}
