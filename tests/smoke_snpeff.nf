// tests/smoke_snpeff.nf
//
// Smoke test for SNPEFF_ANNOTATE (Task 12). Uses a hard-filtered VCF and a
// throwaway SnpEff database built from tests/fixtures/ref.fa +
// tests/fixtures/synth_ref.gtf. tests/run_integration_test.sh builds both and
// can run this test for you; running it by hand requires --vcf, --tbi,
// --snpeff_db_dir and --snpeff_genome_name.
//
// FINAL-REVIEW I8(b): previously only .view()'d the output. It now asserts that
// annotation actually happened.
nextflow.enable.dsl = 2
include { SNPEFF_ANNOTATE } from '../modules/local/snpeff/main.nf'

def gunzip(path) {
    return new java.util.zip.GZIPInputStream(new java.io.FileInputStream(path.toString())).getText('UTF-8')
}

workflow {
    ['vcf', 'tbi', 'snpeff_db_dir', 'snpeff_genome_name'].each { p ->
        if (!params[p]) {
            error "tests/smoke_snpeff.nf requires --${p} (see tests/run_integration_test.sh, which builds the throwaway SnpEff DB this test needs)"
        }
    }

    vcf = file(params.vcf)
    tbi = file(params.tbi)
    db_dir = file(params.snpeff_db_dir)

    ch = Channel.of(['all', vcf, tbi])

    SNPEFF_ANNOTATE(ch, db_dir, params.snpeff_genome_name)

    SNPEFF_ANNOTATE.out.vcf.view()

    SNPEFF_ANNOTATE.out.vcf.subscribe { population, annotated ->
        def lines = gunzip(annotated).readLines()
        def records = lines.findAll { !it.startsWith('#') }
        if (records.isEmpty()) {
            // This is the exact failure the set -euo pipefail fix guards: snpEff
            // failing while bgzip still exits 0 on empty stdin.
            error "annotated VCF contains zero records - snpEff produced no output"
        }
        def unannotated = records.findAll { !it.split('\t')[7].contains('ANN=') }
        if (unannotated) {
            error "${unannotated.size()} of ${records.size()} record(s) have no ANN= field"
        }
        println "ASSERTIONS PASSED: all ${records.size()} records carry an ANN= annotation"
    }
}
