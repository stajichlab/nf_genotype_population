#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENOTYPE_POPULATION } from './workflows/genotype_population.nf'
include { PLOIDY_INFERENCE }    from './subworkflows/local/ploidy_inference/main.nf'

// CORRECTION 2 (Task 13 binding ruling): Channel.fromFilePairs("*.{cram,cram.crai}", ...)
// does not correctly pair X.cram with X.cram.crai (the {a,b} form in
// fromFilePairs matches alternative FULL suffixes, not an optional trailing
// suffix). Use the bash-brace-expansion-style optional-suffix glob instead:
// "*.cram{,.crai}" matches both "X.cram" and "X.cram.crai" and pairs them
// correctly.
def buildCramCh() {
    return Channel
        .fromFilePairs("${params.cram_dir}/*.cram{,.crai}", checkIfExists: true) { file -> file.baseName.replaceFirst(/\.cram$/, '') }
        .map { strain, files -> [strain, files[0], files[1]] }
}

// CORRECTION 1 (Task 13 binding ruling): the brief only wires up the full
// GENOTYPE_POPULATION workflow, which requires params.ploidy_overrides to
// already exist as a human-reviewed CSV derived from a prior run's
// ploidy_review.csv. That leaves no way to produce the FIRST
// ploidy_review.csv. This entry workflow runs PLOIDY_INFERENCE alone, so a
// human can review its ploidy_review.csv output and hand-build
// ploidy_overrides.csv for the next (default) run.
//
// Run with:
//   nextflow run main.nf -entry PLOIDY_ONLY --cram_dir ... --reference ... --reference_fai ... --metadata ...
// NOTE: on this cluster's installed Nextflow (26.04.6), the default "strict
// parser" rejects the `-entry` CLI flag outright ("The `-entry` option is not
// supported with the strict parser"), regardless of how many named workflow
// blocks this file has. Prefix the run with NXF_SYNTAX_PARSER=v1 to restore
// the legacy parser, which does support -entry, e.g.:
//   NXF_SYNTAX_PARSER=v1 nextflow run main.nf -entry PLOIDY_ONLY --cram_dir ... --reference ... --reference_fai ... --metadata ...
workflow PLOIDY_ONLY {
    cram_ch = buildCramCh()

    PLOIDY_INFERENCE(
        cram_ch,
        file(params.reference),
        file(params.reference_fai),
        file(params.metadata),
    )

    // PLOIDY_INFERENCE.out.review is PLOIDY_CROSSCHECK's ploidy_review.csv -
    // the human-reviewable AGREE/DISAGREE report this entry workflow exists
    // to produce. Publish it explicitly (no publishDir anywhere in this
    // pipeline - see the annotated_vcf publish note below for why).
    PLOIDY_INFERENCE.out.review.subscribe { review_csv ->
        def outDir = file(params.outdir)
        outDir.mkdirs()
        review_csv.copyTo(outDir.resolve('ploidy_review.csv'))
    }
}

workflow {
    cram_ch = buildCramCh()

    GENOTYPE_POPULATION(
        cram_ch,
        file(params.reference),
        file(params.reference_fai),
        file(params.reference_dict),
        file(params.metadata),
        file(params.ploidy_overrides),
        file(params.population_sets),
        file(params.snpeff_db_dir),
        params.snpeff_genome_name,
    )

    // No module in this pipeline declares a publishDir for its final output
    // (established pattern - see Task 9's GVCFs, which are deliberately left
    // as ephemeral work/ artifacts; only PLOIDY_INFERENCE's ploidy CSV uses
    // an explicit collectFile(storeDir: ...)). The final annotated VCF is the
    // end deliverable of the whole pipeline, so publish it explicitly here
    // rather than modifying the already-committed SNPEFF_ANNOTATE module.
    // .subscribe() is safe here (unlike subworkflows/local/joint_genotyping's
    // documented preference for .toList() over .subscribe()): this is a
    // terminal sink with no downstream Nextflow channel consumer ordering it
    // needs to respect, not a side-effect populating a lookup another
    // channel operation reads back from.
    GENOTYPE_POPULATION.out.annotated_vcf.subscribe { population, vcf ->
        def outDir = file(params.outdir)
        outDir.mkdirs()
        vcf.copyTo(outDir.resolve(vcf.name))
    }
}
