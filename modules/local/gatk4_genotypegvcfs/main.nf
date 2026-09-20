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
    // GATK4_GENOTYPEGVCFS: --max-alternate-alleles default (6, GATK 4.5.0.0
    // docs) undercounts real diversity in this population, but the site that
    // crashed the run had 50 alleles (GenomicsDBImport's own
    // --genomicsdb-max-alternate-alleles default) genotyped at ploidy 2 =
    // 1275 genotypes, over GenotypeGVCFs's --max-genotype-count default of
    // 1024 - GATK then drops that sample's PL field, which crashes the
    // ExcessHet annotator with "Genotype has no likelihoods"
    // (IllegalStateException) rather than handling it gracefully. Capping at
    // 24 alt alleles keeps genotype count for a diploid site at 25*26/2=325,
    // safely under the 1024 cap, while still allowing far more alleles than
    // the GATK default.
    """
    set -euo pipefail
    gatk --java-options "-Xmx${xmx}g" GenotypeGVCFs \\
        -R "${reference}" \\
        -V "gendb://${genomicsdb}" \\
        --max-alternate-alleles 24 \\
        -O "${population}.vcf.gz"
    """
}
