// tests/smoke_ploidy_inference.nf
//
// Exercises the PLOIDY_INFERENCE subworkflow end to end, which is also the only
// coverage of PLOIDY_CROSSCHECK.
//
// FINAL-REVIEW I8(b): previously only .view()'d the review CSV. It now asserts
// on its content.
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

    // metadata_fixture.txt declares haploid_strain=haploid and
    // diploid_strain=diploid, and the fixtures were built to match, so every row
    // must come back AGREE. A DISAGREE here means either the classifier or the
    // crosscheck regressed.
    PLOIDY_INFERENCE.out.review.subscribe { review_csv ->
        def lines = review_csv.text.readLines()
        if (lines[0] != 'strain,inferred_ploidy,metadata_ploidy,status') {
            error "unexpected review header '${lines[0]}'"
        }
        def rows = lines[1..-1].collect { it.split(',') }
        if (rows.size() != 2) {
            error "expected 2 review rows, got ${rows.size()}"
        }
        def byStrain = rows.collectEntries { [(it[0]): it] }
        ['haploid_strain': 'haploid', 'diploid_strain': 'diploid'].each { strain, ploidy ->
            def row = byStrain[strain]
            if (!row) {
                error "review CSV has no row for ${strain}"
            }
            if (row[1] != ploidy) {
                error "${strain}: inferred_ploidy '${row[1]}', expected '${ploidy}'"
            }
            if (row[3] != 'AGREE') {
                error "${strain}: status '${row[3]}', expected 'AGREE' (inferred=${row[1]}, metadata=${row[2]})"
            }
        }
        println "ASSERTIONS PASSED: both strains AGREE between inference and metadata"
    }

    // The collected per-strain CSV must carry one header plus one row per strain
    // (collectFile with keepHeader/skip:1), not a header repeated per strain.
    PLOIDY_INFERENCE.out.inferred_csv.subscribe { inferred ->
        def lines = inferred.text.readLines().findAll { it.trim() }
        if (lines.size() != 3) {
            error "expected 1 header + 2 data rows in the collected inference CSV, got ${lines.size()} lines: ${lines}"
        }
        if (lines.count { it.startsWith('strain,') } != 1) {
            error "collected inference CSV has a repeated header: ${lines}"
        }
        println "ASSERTIONS PASSED: collected inference CSV has one header and one row per strain"
    }
}
