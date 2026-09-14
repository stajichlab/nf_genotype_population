// subworkflows/local/ploidy_inference/main.nf
include { CUSTOM_HET_PLOIDY } from '../../../modules/local/custom_het_ploidy/main.nf'
include { NQUIRE            } from '../../../modules/local/nquire/main.nf'
include { PLOIDY_CROSSCHECK } from '../../../modules/local/ploidy_crosscheck/main.nf'

workflow PLOIDY_INFERENCE {
    take:
    cram_ch
    reference
    reference_fai
    metadata_txt

    main:
    // params.ploidy_method selects the pluggable inference method (design spec
    // Phase 2). Both branches emit the same [strain, csv] shape with a
    // 'strain,inferred_ploidy,...' header, so PLOIDY_CROSSCHECK downstream
    // (bin/ploidy_crosscheck.py, which reads only the 'strain' and
    // 'inferred_ploidy' columns by name via csv.DictReader) needs no changes
    // regardless of which method ran.
    //
    // Default 'custom_het_script' is unchanged from Phase 1. 'nquire' is new:
    // it has no haploid model (see bin/nquire_ploidy.py's docstring) and
    // deliberately reports 'unknown' rather than an invented threshold for any
    // sample that isn't a strong, unambiguous diploid fit - PLOIDY_CROSSCHECK
    // already maps 'unknown' to status UNKNOWN_INFERENCE, so no wiring change
    // was needed there either.
    def ploidy_method = params.ploidy_method ?: 'custom_het_script'
    if (ploidy_method == 'custom_het_script') {
        CUSTOM_HET_PLOIDY(cram_ch, reference, reference_fai)
        ploidy_csv_ch = CUSTOM_HET_PLOIDY.out.csv
    } else if (ploidy_method == 'nquire') {
        NQUIRE(cram_ch, reference, reference_fai)
        ploidy_csv_ch = NQUIRE.out.csv
    } else {
        error "Unknown params.ploidy_method '${ploidy_method}' - valid values are 'custom_het_script' or 'nquire'."
    }

    ploidy_csv_ch
        .map { strain, csv -> csv }
        .collectFile(name: 'ploidy_inference_all.csv', keepHeader: true, skip: 1, storeDir: "${params.outdir}/ploidy")
        .set { inferred_csv }

    PLOIDY_CROSSCHECK(inferred_csv, metadata_txt)

    emit:
    inferred_csv = inferred_csv
    review       = PLOIDY_CROSSCHECK.out.review
}
