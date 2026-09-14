// tests/smoke_nquire.nf
//
// Exercises the NQUIRE module end to end against the real container
// (built from assets/nquire.def - no upstream package exists for this tool).
//
// Unlike smoke_custom_het_ploidy.nf, this test does NOT assert a specific
// inferred_ploidy value. nQuire has no haploid model (see
// bin/nquire_ploidy.py's docstring) and its diploid-confidence threshold here
// is a documented, uncalibrated placeholder - there is no validated ground
// truth to assert against for these tiny synthetic fixtures, and asserting
// one would be an invented result, not a real one. This test instead checks
// that the module actually runs the full CRAM->BAM->create->denoise->lrdmodel
// pipeline inside the container and produces a well-formed CSV.
nextflow.enable.dsl = 2
include { NQUIRE } from '../modules/local/nquire/main.nf'

workflow {
    ref     = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    ch = Channel.of(
        ['haploid_strain', file("${projectDir}/../tests/fixtures/haploid_strain.cram"), file("${projectDir}/../tests/fixtures/haploid_strain.cram.crai")],
        ['diploid_strain', file("${projectDir}/../tests/fixtures/diploid_strain.cram"), file("${projectDir}/../tests/fixtures/diploid_strain.cram.crai")],
    )
    NQUIRE(ch, ref, ref_fai)

    NQUIRE.out.csv.view()

    NQUIRE.out.csv
        .toList()
        .subscribe { entries ->
            if (entries.size() != 2) {
                error "expected 2 ploidy CSVs, got ${entries.size()}"
            }
            entries.each { strain, csv ->
                def lines = csv.text.readLines()
                if (lines[0] != 'strain,inferred_ploidy,best_fit_model,diploid_delta,triploid_delta,tetraploid_delta,method') {
                    error "${strain}: unexpected CSV header '${lines[0]}'"
                }
                if (lines.size() != 2) {
                    error "${strain}: expected exactly one data row, got ${lines.size() - 1}"
                }
                def fields = lines[1].split(',')
                if (fields[0] != strain) {
                    error "${strain}: CSV strain column says '${fields[0]}'"
                }
                if (!(fields[1] in ['diploid', 'unknown'])) {
                    error "${strain}: inferred_ploidy '${fields[1]}' is neither 'diploid' nor 'unknown' - call_ploidy_from_deltas() only returns those two values"
                }
                if (fields[6] != 'nquire') {
                    error "${strain}: method column says '${fields[6]}', expected 'nquire'"
                }
            }
            println "ASSERTIONS PASSED: NQUIRE module ran end to end in-container and produced well-formed CSVs for both strains (values NOT asserted - see file header)"
        }
}
