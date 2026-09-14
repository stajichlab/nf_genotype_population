// tests/smoke_joint_genotyping.nf
//
// FINAL-REVIEW I8(b): previously only .view()'d both outputs. It now asserts on
// the joint VCF and on the hard-filtered VCF.
// FINAL-REVIEW I9: also asserts the mixed-ploidy claim at the subworkflow level
// (the haploid sample emits one-allele genotypes, the diploid sample two) so a
// regression is caught here and not only by the full integration test.
nextflow.enable.dsl = 2
include { GATK4_HAPLOTYPECALLER } from '../modules/local/gatk4_haplotypecaller/main.nf'
include { JOINT_GENOTYPING }      from '../subworkflows/local/joint_genotyping/main.nf'
include { GATK4_HARDFILTER }      from '../modules/local/gatk4_variantfiltration/main.nf'

def ploidyCodeFor(label) {
    if (label == 'haploid') { return 1 }
    if (label == 'diploid') { return 2 }
    error "Unknown ploidy label '${label}' - only 'haploid' or 'diploid' are valid in ploidy_overrides.csv"
}

def gunzip(path) {
    return new java.util.zip.GZIPInputStream(new java.io.FileInputStream(path.toString())).getText('UTF-8')
}

workflow {
    ref      = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai  = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    ref_dict = file("${projectDir}/../tests/fixtures/ref.dict")
    population_sets_yaml = file("${projectDir}/../tests/fixtures/population_sets_fixture.yaml")

    overrides = Channel.of(
        ['haploid_strain', 'haploid'],
        ['diploid_strain', 'diploid'],
    )
    crams = Channel.of(
        ['haploid_strain', file("${projectDir}/../tests/fixtures/haploid_strain.cram"), file("${projectDir}/../tests/fixtures/haploid_strain.cram.crai")],
        ['diploid_strain', file("${projectDir}/../tests/fixtures/diploid_strain.cram"), file("${projectDir}/../tests/fixtures/diploid_strain.cram.crai")],
    )

    ch = overrides
        .join(crams)
        .map { strain, ploidy_label, cram, crai -> [strain, ploidyCodeFor(ploidy_label), cram, crai] }

    GATK4_HAPLOTYPECALLER(ch, ref, ref_fai, ref_dict)

    JOINT_GENOTYPING(
        GATK4_HAPLOTYPECALLER.out.gvcf,
        population_sets_yaml,
        ref,
        ref_fai,
        ref_dict,
        ['synth_contig1'], // multi-contig support: intervals is now a List<String>
    )

    JOINT_GENOTYPING.out.vcf_by_population.view()

    JOINT_GENOTYPING.out.vcf_by_population.subscribe { population, vcf, tbi ->
        // FINAL-REVIEW C1: `all` is derived from the GVCF channel now, so it
        // must appear even though population_sets_fixture.yaml's `all` key is
        // ignored.
        if (population != 'all') {
            error "expected population 'all', got '${population}'"
        }
        def lines = gunzip(vcf).readLines()
        def chromLine = lines.find { it.startsWith('#CHROM') }
        def samples = chromLine.split('\t')[9..-1]
        if (samples.sort() != ['diploid_strain', 'haploid_strain']) {
            error "joint VCF samples are ${samples}, expected both fixture strains"
        }
        def records = lines.findAll { !it.startsWith('#') }
        if (records.isEmpty()) {
            error "joint VCF contains zero records"
        }

        // FINAL-REVIEW I9: the pipeline's central scientific claim is ploidy-aware
        // per-strain calling surviving a mixed-ploidy GenomicsDB merge. Assert the
        // genotypes themselves, at a site where BOTH samples were called.
        def hapIdx = samples.indexOf('haploid_strain')
        def dipIdx = samples.indexOf('diploid_strain')
        def checked = 0
        records.each { rec ->
            def f = rec.split('\t')
            def hapGt = f[9 + hapIdx].split(':')[0]
            def dipGt = f[9 + dipIdx].split(':')[0]
            if (hapGt == './.' || hapGt == '.' || dipGt == './.' || dipGt == '.') {
                return
            }
            if (!(hapGt ==~ /^[0-9]+$/)) {
                error "haploid_strain GT '${hapGt}' at ${f[0]}:${f[1]} has more than one allele - ploidy-aware calling did not survive joint genotyping"
            }
            if (!(dipGt ==~ /^[0-9]+[\/|][0-9]+$/)) {
                error "diploid_strain GT '${dipGt}' at ${f[0]}:${f[1]} is not a two-allele genotype"
            }
            checked = checked + 1
        }
        if (checked == 0) {
            error "no site had both samples genotyped, so the mixed-ploidy assertion never ran"
        }
        println "ASSERTIONS PASSED: joint VCF has both samples and correct per-strain ploidy at ${checked} shared site(s)"
    }

    GATK4_HARDFILTER(
        JOINT_GENOTYPING.out.vcf_by_population,
        ref,
        ref_fai,
        ref_dict,
    )

    GATK4_HARDFILTER.out.vcf.view()

    GATK4_HARDFILTER.out.vcf.subscribe { population, vcf, tbi ->
        def lines = gunzip(vcf).readLines()
        // VariantFiltration must have declared its filters in the header, and
        // every record must carry a FILTER verdict rather than '.'.
        if (!lines.any { it.startsWith('##FILTER=<ID=SOR3') }) {
            error "hard-filtered VCF header has no SOR3 FILTER declaration"
        }
        def records = lines.findAll { !it.startsWith('#') }
        if (records.isEmpty()) {
            error "hard-filtered VCF contains zero records"
        }
        def unfiltered = records.findAll { it.split('\t')[6] == '.' }
        if (unfiltered) {
            error "${unfiltered.size()} record(s) have an unset FILTER column - VariantFiltration did not run over them"
        }
        println "ASSERTIONS PASSED: hard-filtered VCF has ${records.size()} records, all with a FILTER verdict"
    }
}
