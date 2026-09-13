// tests/smoke_custom_het_ploidy.nf
nextflow.enable.dsl = 2
include { CUSTOM_HET_PLOIDY } from '../modules/local/custom_het_ploidy/main.nf'

workflow {
    ref     = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    ch = Channel.of(
        ['haploid_strain', file("${projectDir}/../tests/fixtures/haploid_strain.cram"), file("${projectDir}/../tests/fixtures/haploid_strain.cram.crai")],
        ['diploid_strain', file("${projectDir}/../tests/fixtures/diploid_strain.cram"), file("${projectDir}/../tests/fixtures/diploid_strain.cram.crai")],
    )
    CUSTOM_HET_PLOIDY(ch, ref, ref_fai)
    CUSTOM_HET_PLOIDY.out.csv.view()
}
