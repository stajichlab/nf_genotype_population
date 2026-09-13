// tests/smoke_haplotypecaller.nf
//
// FINAL-REVIEW I8(b): previously only .view()'d the output tuple, so the test
// could not fail on a bad GVCF. It now decompresses each GVCF and asserts on its
// header and records.
nextflow.enable.dsl = 2
include { GATK4_HAPLOTYPECALLER } from '../modules/local/gatk4_haplotypecaller/main.nf'

def ploidyCodeFor(label) {
    if (label == 'haploid') { return 1 }
    if (label == 'diploid') { return 2 }
    error "Unknown ploidy label '${label}' - only 'haploid' or 'diploid' are valid in ploidy_overrides.csv"
}

// bgzip output is a series of concatenated gzip members; GZIPInputStream reads
// them all, so this yields the whole file as text.
def gunzip(path) {
    return new java.util.zip.GZIPInputStream(new java.io.FileInputStream(path.toString())).getText('UTF-8')
}

workflow {
    ref      = file("${projectDir}/../tests/fixtures/ref.fa")
    ref_fai  = file("${projectDir}/../tests/fixtures/ref.fa.fai")
    ref_dict = file("${projectDir}/../tests/fixtures/ref.dict")

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

    GATK4_HAPLOTYPECALLER.out.gvcf.view()

    GATK4_HAPLOTYPECALLER.out.gvcf
        .toList()
        .subscribe { entries ->
            if (entries.size() != 2) {
                error "expected 2 GVCFs, got ${entries.size()}"
            }
            entries.each { strain, gvcf, tbi ->
                if (!tbi.exists() || tbi.size() == 0) {
                    error "${strain}: missing or empty GVCF index"
                }
                def lines = gunzip(gvcf).readLines()
                def chromLine = lines.find { it.startsWith('#CHROM') }
                if (!chromLine) {
                    error "${strain}: GVCF has no #CHROM header line"
                }
                // The sample column must carry the strain's own name, which is
                // what GATK resolves from the CRAM's @RG SM tag. A mismatch here
                // would silently mislabel that strain in the joint callset.
                def samples = chromLine.split('\t')[9..-1]
                if (samples != [strain]) {
                    error "${strain}: GVCF sample columns are ${samples}, expected exactly ['${strain}']"
                }
                // -ERC GVCF must produce non-variant reference blocks.
                def records = lines.findAll { !it.startsWith('#') }
                if (records.isEmpty()) {
                    error "${strain}: GVCF contains zero records"
                }
                if (!records.any { it.contains('<NON_REF>') }) {
                    error "${strain}: GVCF has no <NON_REF> allele - -ERC GVCF did not take effect"
                }
            }
            println "ASSERTIONS PASSED: both GVCFs carry the right sample name, real records, and <NON_REF> blocks"
        }
}
