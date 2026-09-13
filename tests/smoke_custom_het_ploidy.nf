// tests/smoke_custom_het_ploidy.nf
//
// FINAL-REVIEW I8(b): this test previously only called .view(), so pass/fail was
// a human reading stdout and nothing could fail automatically. It now asserts on
// the actual CSV content, and errors (non-zero exit) if the content is wrong.
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

    // The fixtures are constructed so the haploid strain has real SNP calls but
    // zero heterozygosity, and the diploid strain has a genuine ~50/50 read
    // mixture at each hapB SNP site (see tests/fixtures/make_fixtures.sh). The
    // classifier must therefore recover the ploidy each fixture was built to
    // have - that is the behaviour under test, not merely "the process exited 0".
    CUSTOM_HET_PLOIDY.out.csv
        .toList()
        .subscribe { entries ->
            if (entries.size() != 2) {
                error "expected 2 ploidy CSVs, got ${entries.size()}"
            }
            entries.each { strain, csv ->
                def lines = csv.text.readLines()
                if (lines[0] != 'strain,inferred_ploidy,het_fraction,method') {
                    error "${strain}: unexpected CSV header '${lines[0]}'"
                }
                if (lines.size() != 2) {
                    error "${strain}: expected exactly one data row, got ${lines.size() - 1}"
                }
                def fields = lines[1].split(',')
                if (fields[0] != strain) {
                    error "${strain}: CSV strain column says '${fields[0]}'"
                }
                def expected = (strain == 'haploid_strain') ? 'haploid' : 'diploid'
                if (fields[1] != expected) {
                    error "${strain}: inferred ploidy '${fields[1]}', expected '${expected}' (het_fraction=${fields[2]})"
                }
                if (fields[2] == 'NA') {
                    error "${strain}: het_fraction is NA - the classifier had no sites to work from"
                }
            }
            println "ASSERTIONS PASSED: haploid_strain=haploid, diploid_strain=diploid, both with real het_fraction values"
        }
}
