// tests/smoke_joint_genotyping.nf
nextflow.enable.dsl = 2
include { GATK4_HAPLOTYPECALLER } from '../modules/local/gatk4_haplotypecaller/main.nf'
include { JOINT_GENOTYPING }      from '../subworkflows/local/joint_genotyping/main.nf'
include { GATK4_HARDFILTER }      from '../modules/local/gatk4_variantfiltration/main.nf'

def ploidyCodeFor(label) {
    if (label == 'haploid') { return 1 }
    if (label == 'diploid') { return 2 }
    error "Unknown ploidy label '${label}' - only 'haploid' or 'diploid' are valid in ploidy_overrides.csv"
}

workflow {
    ref      = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai  = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    ref_dict = file("${projectDir}/../tests/fixtures/ref.dict")
    population_sets_yaml = file("${projectDir}/../tests/fixtures/population_sets_fixture.yaml")

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

    JOINT_GENOTYPING(
        GATK4_HAPLOTYPECALLER.out.gvcf,
        population_sets_yaml,
        ref,
        ref_fai,
        ref_dict,
        'synth_contig1',
    )

    JOINT_GENOTYPING.out.vcf_by_population.view()

    GATK4_HARDFILTER(
        JOINT_GENOTYPING.out.vcf_by_population,
        ref,
        ref_fai,
        ref_dict,
    )

    GATK4_HARDFILTER.out.vcf.view()
}
