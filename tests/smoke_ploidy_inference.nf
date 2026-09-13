// tests/smoke_ploidy_inference.nf
nextflow.enable.dsl = 2
include { PLOIDY_INFERENCE } from '../subworkflows/local/ploidy_inference/main.nf'

workflow {
    ref     = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    meta    = file("${projectDir}/../tests/fixtures/metadata_fixture.txt")
    ch = Channel.of(
        ['haploid_strain', file("${projectDir}/../tests/fixtures/haploid_strain.cram"), file("${projectDir}/../tests/fixtures/haploid_strain.cram.crai")],
        ['diploid_strain', file("${projectDir}/../tests/fixtures/diploid_strain.cram"), file("${projectDir}/../tests/fixtures/diploid_strain.cram.crai")],
    )
    PLOIDY_INFERENCE(ch, ref, ref_fai, meta)
    PLOIDY_INFERENCE.out.review.view { "review: " + it.text }
}
