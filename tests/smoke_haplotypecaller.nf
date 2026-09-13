// tests/smoke_haplotypecaller.nf
nextflow.enable.dsl = 2
include { GATK4_HAPLOTYPECALLER } from '../modules/local/gatk4_haplotypecaller/main.nf'

def ploidyCodeFor(label) {
    if (label == 'haploid') { return 1 }
    if (label == 'diploid') { return 2 }
    error "Unknown ploidy label '${label}' - only 'haploid' or 'diploid' are valid in ploidy_overrides.csv"
}

workflow {
    ref      = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai  = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    ref_dict = file("${projectDir}/../tests/fixtures/ref.dict")

    overrides = Channel.of(
        ['haploid_strain', 'haploid'],
        ['diploid_strain', 'diploid'],
    )
    crams = Channel.of(
        ['haploid_strain', file("${projectDir}/../tests/fixtures/haploid_strain.cram"), file("${projectDir}/../tests/fixtures/haploid_strain.cram.crai")],
        ['diploid_strain', file("${projectDir}/../tests/fixtures/diploid_strain.cram"), file("${projectDir}/../tests/fixtures/diploid_strain.cram.crai")],
    )

    ch = overrides
        .join(crams)
        .map { strain, ploidy_label, cram, crai -> [strain, ploidyCodeFor(ploidy_label), cram, crai] }

    GATK4_HAPLOTYPECALLER(ch, ref, ref_fai, ref_dict)
    GATK4_HAPLOTYPECALLER.out.gvcf.view()
}
