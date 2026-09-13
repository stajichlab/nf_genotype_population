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

// FINAL-REVIEW C2: there was no parameter validation anywhere in the branch -
// no nf-schema, no null checks. A missing --cram_dir produced a glob of
// "null/*.cram" rather than a readable error, and a missing --outdir made
// file(params.outdir) throw on a null. Check every required parameter up front
// and report ALL missing ones at once, so a user fixes one command line
// instead of discovering the requirements one failed run at a time.
def requireParams(Map required) {
    def missing = required.findAll { name, description -> !params[name] }
    if (missing) {
        def lines = missing.collect { name, description -> "  --${name}\t${description}" }.join('\n')
        error "Missing required parameter(s):\n${lines}\n\nSee README.md for the full command for this entry point."
    }
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
    requireParams([
        cram_dir     : 'FLAT directory of <strain>.cram / <strain>.cram.crai pairs',
        reference    : 'reference genome FASTA',
        reference_fai: 'reference genome .fai index',
        metadata     : 'metadata.txt (strain -> ploidy, for the AGREE/DISAGREE crosscheck)',
    ])

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
    requireParams([
        cram_dir          : 'FLAT directory of <strain>.cram / <strain>.cram.crai pairs',
        reference         : 'reference genome FASTA',
        reference_fai     : 'reference genome .fai index',
        reference_dict    : 'reference genome GATK sequence dictionary (.dict)',
        metadata          : 'metadata.txt (strain -> ploidy)',
        ploidy_overrides  : 'human-reviewed strain,ploidy CSV - ONE ROW PER STRAIN (see README)',
        population_sets   : 'population_sets YAML (defaults to assets/population_sets.yaml)',
        snpeff_db_dir     : 'prebuilt SnpEff database directory',
        snpeff_genome_name: 'SnpEff genome name inside that database',
    ])

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

    // FINAL-REVIEW I6: PLOIDY_INFERENCE always runs in the default workflow too
    // (a deliberate QC side effect, kept as-is), but its outputs were previously
    // discarded here - so a production run paid for every CUSTOM_HET_PLOIDY job
    // and kept no diagnostic record, and ploidy drift against the frozen
    // ploidy_overrides.csv was invisible. Publish both, matching what
    // PLOIDY_ONLY already does, so every production run leaves the same record.
    GENOTYPE_POPULATION.out.ploidy_review.subscribe { review_csv ->
        def outDir = file(params.outdir)
        outDir.mkdirs()
        review_csv.copyTo(outDir.resolve('ploidy_review.csv'))
    }
    GENOTYPE_POPULATION.out.inferred_csv.subscribe { inferred ->
        def outDir = file(params.outdir)
        outDir.mkdirs()
        inferred.copyTo(outDir.resolve('ploidy_inference_all.csv'))
    }

    // FINAL-REVIEW I10: spec §5 requires per-strain GVCFs to be retained for
    // re-genotyping as new strains are added. They previously existed only in
    // work/, so the first scratch cleanup destroyed them. Publish each strain's
    // .g.vcf.gz and its .tbi index into ${params.outdir}/gvcfs/.
    GENOTYPE_POPULATION.out.gvcf.subscribe { strain, gvcf, tbi ->
        def gvcfDir = file("${params.outdir}/gvcfs")
        gvcfDir.mkdirs()
        gvcf.copyTo(gvcfDir.resolve(gvcf.name))
        tbi.copyTo(gvcfDir.resolve(tbi.name))
    }
}
