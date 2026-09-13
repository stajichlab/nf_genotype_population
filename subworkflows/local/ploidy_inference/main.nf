// subworkflows/local/ploidy_inference/main.nf
include { CUSTOM_HET_PLOIDY } from '../../../modules/local/custom_het_ploidy/main.nf'
include { PLOIDY_CROSSCHECK } from '../../../modules/local/ploidy_crosscheck/main.nf'

workflow PLOIDY_INFERENCE {
    take:
    cram_ch
    reference
    reference_fai
    metadata_txt

    main:
    CUSTOM_HET_PLOIDY(cram_ch, reference, reference_fai)

    CUSTOM_HET_PLOIDY.out.csv
        .map { strain, csv -> csv }
        .collectFile(name: 'ploidy_inference_all.csv', keepHeader: true, skip: 1, storeDir: "${params.outdir}/ploidy")
        .set { inferred_csv }

    PLOIDY_CROSSCHECK(inferred_csv, metadata_txt)

    emit:
    inferred_csv = inferred_csv
    review       = PLOIDY_CROSSCHECK.out.review
}
