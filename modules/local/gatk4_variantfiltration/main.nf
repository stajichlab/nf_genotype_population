// modules/local/gatk4_variantfiltration/main.nf
process GATK4_HARDFILTER {
    tag "$population"
    label 'process_low'
    // FINAL-REVIEW I2: container (GATK 4.5.0.0) declared centrally in
    // conf/modules.config (spec §5), which records the upstream
    // docker://broadinstitute/gatk:4.5.0.0 source and the pre-pull rationale.

    input:
    tuple val(population), path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict

    output:
    tuple val(population), path("${population}.filtered.vcf.gz"), path("${population}.filtered.vcf.gz.tbi"), emit: vcf

    script:
    // FINAL-REVIEW I3: `set -euo pipefail` is REQUIRED here, and its absence was
    // the most dangerous instance of that omission in the branch. This project's
    // nextflow.config sets process.shell = ['/bin/bash', '-l'], which suppresses
    // the -ue flags Nextflow normally puts in its task shebang, so errexit was
    // OFF. This script runs FIVE sequential gatk commands; without errexit a
    // failure in an early SelectVariants or VariantFiltration would not stop the
    // script, and if MergeVcfs still emitted ${population}.filtered.vcf.gz,
    // Nextflow would see the declared output file and report SUCCESS - a
    // silently incomplete callset. (The same class of bug was already observed
    // for real in this repo in SNPEFF_ANNOTATE, commit 9e44236.)
    //
    // FINAL-REVIEW I5: cap the JVM heap at ~80% of this task's granted memory
    // (see GATK4_HAPLOTYPECALLER for the full rationale).
    def xmx = (task.memory.toGiga() * 0.8) as int
    """
    set -euo pipefail
    gatk --java-options "-Xmx${xmx}g" SelectVariants -R "${reference}" -V "${vcf}" --select-type-to-include SNP -O "${population}.snp.vcf.gz"
    gatk --java-options "-Xmx${xmx}g" SelectVariants -R "${reference}" -V "${vcf}" --select-type-to-include INDEL -O "${population}.indel.vcf.gz"

    gatk --java-options "-Xmx${xmx}g" VariantFiltration -R "${reference}" -V "${population}.snp.vcf.gz" \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "FS > 60.0" --filter-name "FS60" \\
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \\
        -O "${population}.snp.filtered.vcf.gz"

    gatk --java-options "-Xmx${xmx}g" VariantFiltration -R "${reference}" -V "${population}.indel.vcf.gz" \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "FS > 200.0" --filter-name "FS200" \\
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        --filter-expression "SOR > 10.0" --filter-name "SOR10" \\
        -O "${population}.indel.filtered.vcf.gz"

    gatk --java-options "-Xmx${xmx}g" MergeVcfs \\
        -I "${population}.snp.filtered.vcf.gz" \\
        -I "${population}.indel.filtered.vcf.gz" \\
        -O "${population}.filtered.vcf.gz"
    """
}
