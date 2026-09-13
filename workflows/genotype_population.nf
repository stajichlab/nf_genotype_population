// workflows/genotype_population.nf
include { PLOIDY_INFERENCE }       from '../subworkflows/local/ploidy_inference/main.nf'
include { GATK4_HAPLOTYPECALLER }  from '../modules/local/gatk4_haplotypecaller/main.nf'
include { JOINT_GENOTYPING }       from '../subworkflows/local/joint_genotyping/main.nf'
include { GATK4_HARDFILTER }       from '../modules/local/gatk4_variantfiltration/main.nf'
include { SNPEFF_ANNOTATE }        from '../modules/local/snpeff/main.nf'

def ploidyCodeFor(label) {
    if (label == 'haploid') { return 1 }
    if (label == 'diploid') { return 2 }
    error "Unknown ploidy label '${label}' in ploidy_overrides.csv - expected 'haploid' or 'diploid'"
}

workflow GENOTYPE_POPULATION {
    take:
    cram_ch           // channel: [strain, cram, crai]
    reference
    reference_fai
    reference_dict
    metadata_txt
    ploidy_overrides  // path to human-reviewed strain,ploidy CSV
    population_sets_yaml
    snpeff_db_dir
    snpeff_genome_name

    main:
    PLOIDY_INFERENCE(cram_ch, reference, reference_fai, metadata_txt)

    overrides_ch = Channel
        .fromPath(ploidy_overrides)
        .splitCsv(header: true)
        .map { row -> [row.strain, row.ploidy] }

    // failOnMismatch: true guards against a strain silently vanishing from
    // calling because it has a CRAM but no matching row in
    // ploidy_overrides.csv (e.g. a typo) - a plain inner join would otherwise
    // drop it with zero warning, a real correctness hazard at population
    // scale (~278 strains). This also catches the reverse case (an
    // overrides.csv row with no matching CRAM).
    ploidy_cram_ch = overrides_ch
        .join(cram_ch, failOnMismatch: true)
        .map { strain, ploidy_label, cram, crai -> [strain, ploidyCodeFor(ploidy_label), cram, crai] }

    GATK4_HAPLOTYPECALLER(ploidy_cram_ch, reference, reference_fai, reference_dict)

    // CORRECTION 3 (Task 13 binding ruling): Groovy's File.baseName only strips
    // the LAST extension, so for "ref.fa.fai" it yields "ref.fa" - not a valid
    // GATK contig/interval name. Derive the interval from the actual first
    // column of the first line of the .fai file instead (the contig name).
    // NOTE: this only covers a single contig/interval - a known Phase 1
    // limitation for multi-contig real genomes (out of scope here, matches
    // Task 10's own already-documented limitation). Do NOT expand this to
    // multi-interval support - warn at runtime instead so a multi-contig
    // reference doesn't silently only genotype its first contig.
    def faiLines = reference_fai.text.readLines()
    if (faiLines.size() > 1) {
        log.warn "reference_fai has ${faiLines.size()} contigs, but this Phase 1 pipeline only genotypes the FIRST one ('${faiLines[0].split('\t')[0]}'). Multi-contig/multi-interval support is out of scope (see Task 10's own documented limitation)."
    }
    def intervals = faiLines[0].split('\t')[0]

    JOINT_GENOTYPING(
        GATK4_HAPLOTYPECALLER.out.gvcf,
        population_sets_yaml,
        reference,
        reference_fai,
        reference_dict,
        intervals,
    )

    GATK4_HARDFILTER(JOINT_GENOTYPING.out.vcf_by_population, reference, reference_fai, reference_dict)
    SNPEFF_ANNOTATE(GATK4_HARDFILTER.out.vcf, snpeff_db_dir, snpeff_genome_name)

    emit:
    annotated_vcf = SNPEFF_ANNOTATE.out.vcf
}
