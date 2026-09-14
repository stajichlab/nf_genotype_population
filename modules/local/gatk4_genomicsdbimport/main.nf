// modules/local/gatk4_genomicsdbimport/main.nf
process GATK4_GENOMICSDBIMPORT {
    tag "$population"
    label 'process_medium'
    // FINAL-REVIEW I2: container (GATK 4.5.0.0) declared centrally in
    // conf/modules.config (spec §5), which records the upstream
    // docker://broadinstitute/gatk:4.5.0.0 source and the pre-pull rationale.

    input:
    tuple val(population), path(gvcfs), path(tbis)
    path reference
    path reference_fai
    path reference_dict
    val intervals // List<String> of contig names - one -L per contig (multi-contig support)

    output:
    tuple val(population), path("${population}_gdb"), emit: genomicsdb

    script:
    // FINAL-REVIEW C1 (defence in depth): JOINT_GENOTYPING already refuses to
    // build a population with no GVCFs, so this should be unreachable. Kept as
    // a guard because GATK's own error for a zero-`-V` command line is an
    // opaque usage dump.
    if (!gvcfs) {
        error "GATK4_GENOMICSDBIMPORT received zero GVCFs for population '${population}'; refusing to invoke GATK with an empty -V list."
    }
    def gvcf_args = gvcfs.collect { "-V \"${it}\"" }.join(' ')
    // MULTI-CONTIG SUPPORT: `intervals` is now the full list of contig names
    // from the reference .fai (see workflows/genotype_population.nf), not a
    // single contig string. GenomicsDBImport takes one -L per interval, not a
    // comma/space-joined value, so build one flag per contig.
    def interval_args = intervals.collect { "-L \"${it}\"" }.join(' ')
    // FINAL-REVIEW I5: cap the JVM heap at ~80% of this task's granted memory
    // (see GATK4_HAPLOTYPECALLER for the full rationale).
    def xmx = (task.memory.toGiga() * 0.8) as int
    """
    set -euo pipefail
    gatk --java-options "-Xmx${xmx}g" GenomicsDBImport \\
        ${gvcf_args} \\
        --genomicsdb-workspace-path "${population}_gdb" \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --bypass-feature-reader \\
        --merge-input-intervals \\
        ${interval_args}
    """
}
