// subworkflows/local/joint_genotyping/main.nf
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/local/gatk4_genomicsdbimport/main.nf'
include { GATK4_GENOTYPEGVCFS }   from '../../../modules/local/gatk4_genotypegvcfs/main.nf'

workflow JOINT_GENOTYPING {
    take:
    gvcf_ch              // channel: [strain, gvcf, tbi]
    population_sets_yaml // path to assets/population_sets.yaml
    reference
    reference_fai
    reference_dict
    intervals            // e.g. the whole-genome interval list/contig name

    main:
    // snakeyaml verified working under this project's Nextflow version
    // (26.04.6) inside a DSL2 workflow block: see task-10-report.md.
    def yaml_groups = (new org.yaml.snakeyaml.Yaml().load(population_sets_yaml.text) ?: [:])

    // FINAL-REVIEW C1: the `all` population is NO LONGER read from the YAML.
    //
    // Previously `all` was a frozen 279-entry strain list in
    // assets/population_sets.yaml, intersected against the GVCFs actually
    // produced. That silently dropped (a) any strain that got a GVCF but was
    // absent from the file - the expected steady state as new strains are
    // added, which is precisely the incremental-re-genotyping story spec §5 is
    // built around - and (b) any file-listed strain with no GVCF. Spec §2
    // requires `all` to contain "every included strain ... regardless of
    // population-file completeness", so it must be derived from the run
    // itself, not from a hand-maintained file.
    //
    // The YAML is retained ONLY for named sub-population slices (a Phase 2
    // concept; today the file holds just the legacy `all` placeholder key,
    // which is deliberately ignored here).
    def sub_populations = yaml_groups.findAll { name, strains -> name != 'all' }
    if (yaml_groups.containsKey('all')) {
        log.info "population_sets YAML contains an 'all' key; it is ignored. The 'all' population is derived at runtime from every strain that produced a GVCF in this run (spec §2)."
    }

    // Build the strain -> [gvcf, tbi] lookup via .toList() rather than
    // .subscribe(): a subscribe-based side effect into a plain Groovy map is
    // not guaranteed-ordered relative to downstream channel consumption, and
    // that risk grows with channel size (see Task 13's ~278-strain run). Using
    // .toList() forces full channel materialization before we read it, so the
    // lookup is always complete before population_ch is built.
    gvcf_map_ch = gvcf_ch
        .toList()
        .map { entries ->
            def m = [:]
            entries.each { strain, gvcf, tbi -> m[strain] = [gvcf, tbi] }
            return m
        }

    population_ch = gvcf_map_ch
        .flatMap { gvcf_map ->
            // FINAL-REVIEW C1: never hand GenomicsDBImport an empty -V list;
            // GATK's own failure for that case is an opaque usage error.
            if (gvcf_map.isEmpty()) {
                error "JOINT_GENOTYPING received zero GVCFs, so no population can be genotyped. Check that --cram_dir matched CRAMs and that GATK4_HAPLOTYPECALLER produced output."
            }

            def groups = []

            // The built-in `all` group: every strain with a GVCF in THIS run.
            def all_strains = gvcf_map.keySet().sort()
            log.info "JOINT_GENOTYPING: population 'all' derived from this run's GVCFs (${all_strains.size()} strains)."
            groups << ['all', all_strains.collect { gvcf_map[it][0] }, all_strains.collect { gvcf_map[it][1] }]

            // Named sub-population slices from the YAML (Phase 2).
            sub_populations.each { name, strains ->
                def present = strains.findAll { gvcf_map.containsKey(it) }
                def missing = strains.findAll { !gvcf_map.containsKey(it) }
                if (missing) {
                    log.warn "population '${name}' lists ${missing.size()} strain(s) with no GVCF in this run; they are excluded from that callset: ${missing.join(', ')}"
                }
                if (!present) {
                    error "population '${name}' from the population_sets YAML has no strain with a GVCF in this run, so its GenomicsDBImport -V list would be empty. Remove the group or correct its strain names."
                }
                groups << [name, present.collect { gvcf_map[it][0] }, present.collect { gvcf_map[it][1] }]
            }

            return groups
        }

    GATK4_GENOMICSDBIMPORT(population_ch, reference, reference_fai, reference_dict, intervals)
    GATK4_GENOTYPEGVCFS(GATK4_GENOMICSDBIMPORT.out.genomicsdb, reference, reference_fai, reference_dict)

    emit:
    vcf_by_population = GATK4_GENOTYPEGVCFS.out.vcf
}
