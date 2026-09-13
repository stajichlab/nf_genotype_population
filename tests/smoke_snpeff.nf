// tests/smoke_snpeff.nf
//
// Smoke test for SNPEFF_ANNOTATE (Task 12). Uses Task 11's synthetic
// all.filtered.vcf.gz output and a throwaway SnpEff database built from
// tests/fixtures/ref.fa + tests/fixtures/synth_ref.gtf (see task-12-report
// for the build recipe). This only confirms the module runs end-to-end
// and produces an ANN= field - it is not testing real annotation accuracy.
nextflow.enable.dsl = 2
include { SNPEFF_ANNOTATE } from '../modules/local/snpeff/main.nf'

workflow {
    vcf = file(params.vcf)
    tbi = file(params.tbi)
    db_dir = file(params.snpeff_db_dir)

    ch = Channel.of(['all', vcf, tbi])

    SNPEFF_ANNOTATE(ch, db_dir, params.snpeff_genome_name)

    SNPEFF_ANNOTATE.out.vcf.view()
}
