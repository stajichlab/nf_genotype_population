// workflows/genotype_population.nf
include { PLOIDY_INFERENCE }       from '../subworkflows/local/ploidy_inference/main.nf'
include { GATK4_HAPLOTYPECALLER }  from '../modules/local/gatk4_haplotypecaller/main.nf'
include { JOINT_GENOTYPING }       from '../subworkflows/local/joint_genotyping/main.nf'
include { GATK4_HARDFILTER }       from '../modules/local/gatk4_variantfiltration/main.nf'
include { VARIANT_QC_FILTER }      from '../modules/local/variant_qc_filter/main.nf'
include { SNPEFF_ANNOTATE }        from '../modules/local/snpeff/main.nf'
include { SNP_ALIGNMENT; IQTREE }  from '../modules/local/snp_tree/main.nf'

def ploidyCodeFor(label) {
    if (label == 'haploid') { return 1 }
    if (label == 'diploid') { return 2 }
    error "Unknown ploidy label '${label}' in ploidy_overrides.csv - expected 'haploid' or 'diploid'"
}

// FINAL-REVIEW M4: spec §2 requires computationally-derived `_hap1`/`_hap2`
// haplotype rows to be "excluded from calling entirely". That rule previously
// had NO implementation - it held only incidentally, because
// reconciliation_table.csv happens to contain zero `_hap` rows today, and the
// population YAML was generated from that table. Any `_hap` CRAM appearing in
// --cram_dir with a matching overrides row would have been called. Make the
// rule explicit so it holds by construction rather than by accident.
def isHaplotypeRow(strain) {
    return strain ==~ /.*_hap\d+$/
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
        .filter { strain, ploidy ->
            if (isHaplotypeRow(strain)) {
                log.warn "excluding computationally-derived haplotype row '${strain}' from calling (spec §2)"
                return false
            }
            return true
        }

    // The same exclusion must be applied to the CRAM side of the join, not just
    // the overrides side: the join below uses failOnMismatch, so filtering only
    // one side would turn a correctly-excluded `_hap` CRAM into a hard abort
    // instead of a clean exclusion.
    callable_cram_ch = cram_ch.filter { strain, cram, crai -> !isHaplotypeRow(strain) }

    // failOnMismatch: true guards against a strain silently vanishing from
    // calling because it has a CRAM but no matching row in
    // ploidy_overrides.csv (e.g. a typo) - a plain inner join would otherwise
    // drop it with zero warning, a real correctness hazard at population
    // scale (~278 strains). This also catches the reverse case (an
    // overrides.csv row with no matching CRAM).
    ploidy_cram_ch = overrides_ch
        .join(callable_cram_ch, failOnMismatch: true)
        .map { strain, ploidy_label, cram, crai -> [strain, ploidyCodeFor(ploidy_label), cram, crai] }

    GATK4_HAPLOTYPECALLER(ploidy_cram_ch, reference, reference_fai, reference_dict)

    // MULTI-CONTIG SUPPORT: derive the FULL interval list from every contig
    // in the .fai (one -L per contig passed through to GenomicsDBImport), not
    // just the first line. This lifts the earlier Phase 1 single-contig
    // limitation (CORRECTION 3's underlying fix - deriving contig names from
    // the .fai's first column rather than reference_fai.baseName's wrong
    // filename-derived value - is unchanged and still required).
    def faiLines = reference_fai.text.readLines().findAll { it.trim() }
    def intervals = faiLines.collect { it.split('\t')[0] }
    log.info "Joint genotyping will cover all ${intervals.size()} contig(s) in ${reference_fai.name}."

    JOINT_GENOTYPING(
        GATK4_HAPLOTYPECALLER.out.gvcf,
        population_sets_yaml,
        reference,
        reference_fai,
        reference_dict,
        intervals,
    )

    GATK4_HARDFILTER(JOINT_GENOTYPING.out.vcf_by_population, reference, reference_fai, reference_dict)

    // Genotype (GQ/DP/allele balance) and site (mask, depth ceiling,
    // missingness, MAF) QC before annotation. Hard filters only flag sites;
    // without this step every hard-filtered record, and every low-quality
    // genotype, reached SnpEff and the published VCF.
    if (!params.mask_bed) {
        log.warn "--mask_bed not set: VARIANT_QC_FILTER will not drop repeat / low-complexity sites"
    }
    def mask_bed = params.mask_bed ? file(params.mask_bed, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    // Population groups. --population_mode subset (default): every named group
    // in the population YAML is cut from the `all` callset here (bcftools -S),
    // then filtered with its own depth ceiling, missingness and MAF. GATK site
    // annotations (QD, FS, SOR, MQ) stay those computed over all strains.
    // --population_mode regenotype: JOINT_GENOTYPING calls each group
    // separately from its GVCFs (slow), and those callsets arrive here.
    // --output_prefix P names every output P.<pop>.* as in the bash pipeline.
    def pop_raw  = (new org.yaml.snakeyaml.Yaml().load(population_sets_yaml.text) ?: [:])
    def pop_sets = (pop_raw.containsKey('Populations') ? (pop_raw.Populations ?: [:]) : pop_raw)
        .findAll { name, strains -> name != 'all' }
    def pfx = params.output_prefix ? "${params.output_prefix}." : ''
    def qc_in_ch = GATK4_HARDFILTER.out.vcf.map { pop, vcf, tbi -> ["${pfx}${pop}".toString(), vcf, tbi, []] }
    if (params.population_mode == 'subset' && pop_sets) {
        pop_sets.each { name, strains ->
            if (!strains) { error "population '${name}' in ${population_sets_yaml.name} lists no strains" }
        }
        log.info "VARIANT_QC_FILTER: ${pop_sets.size()} population group(s) cut from the 'all' callset: ${pop_sets.collect { n, s -> "${n} (${s.size()})" }.join(', ')}"
        qc_in_ch = qc_in_ch.mix(
            GATK4_HARDFILTER.out.vcf
                .filter { pop, vcf, tbi -> pop == 'all' }
                .flatMap { pop, vcf, tbi -> pop_sets.collect { name, strains -> ["${pfx}${name}".toString(), vcf, tbi, strains.collect { it.toString() }] } }
        )
    }
    VARIANT_QC_FILTER(qc_in_ch, mask_bed)

    // Annotate both QC outputs: <pop>.qc (all variant types, no MAF floor) and
    // <pop>.snps.maf (biallelic SNPs, MAF floor).
    SNPEFF_ANNOTATE(VARIANT_QC_FILTER.out.qc.mix(VARIANT_QC_FILTER.out.snps), snpeff_db_dir, snpeff_genome_name)

    // Strain tree from the biallelic SNP set (port of the bash pipeline's
    // 06_make_SNP_tree.sh + 07_iqtree.sh). --skip_tree turns it off.
    def tree_ch = Channel.empty()
    def aln_ch = Channel.empty()
    def aln_stats_ch = Channel.empty()
    if (!params.skip_tree) {
        SNP_ALIGNMENT(VARIANT_QC_FILTER.out.snps)
        IQTREE(SNP_ALIGNMENT.out.alignment)
        tree_ch = IQTREE.out.tree
        aln_ch = SNP_ALIGNMENT.out.alignment
        aln_stats_ch = SNP_ALIGNMENT.out.stats
    }

    emit:
    annotated_vcf = SNPEFF_ANNOTATE.out.vcf
    filter_stats  = VARIANT_QC_FILTER.out.stats
    tree          = tree_ch
    alignment     = aln_ch
    alignment_stats = aln_stats_ch
    // FINAL-REVIEW I6: PLOIDY_INFERENCE runs on every strain in the production
    // workflow (a deliberate QC side effect), but neither of its outputs was
    // emitted, so ~278 CUSTOM_HET_PLOIDY jobs per run produced nothing any
    // consumer could see, and a run showing ploidy drift against the frozen
    // ploidy_overrides.csv left no record of it. Emit both so main.nf can
    // publish them, exactly as PLOIDY_ONLY already does.
    inferred_csv  = PLOIDY_INFERENCE.out.inferred_csv
    ploidy_review = PLOIDY_INFERENCE.out.review
    // FINAL-REVIEW I10: spec §5 requires per-strain GVCFs to be retained
    // ("kept ... for re-genotyping as new strains are added"), but they lived
    // only in work/, which .gitignore excludes and which is routine scratch -
    // the first cleanup destroyed them and the documented incremental
    // re-genotyping story became a full recall from CRAM. Emit them so main.nf
    // publishes them to ${params.outdir}/gvcfs/.
    gvcf          = GATK4_HAPLOTYPECALLER.out.gvcf
}
