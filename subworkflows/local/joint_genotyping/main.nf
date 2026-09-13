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
    def populations = new org.yaml.snakeyaml.Yaml().load(population_sets_yaml.text)

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
            populations.collect { name, strains ->
                def present = strains.findAll { gvcf_map.containsKey(it) }
                [name, present.collect { gvcf_map[it][0] }, present.collect { gvcf_map[it][1] }]
            }
        }

    GATK4_GENOMICSDBIMPORT(population_ch, reference, reference_fai, reference_dict, intervals)
    GATK4_GENOTYPEGVCFS(GATK4_GENOMICSDBIMPORT.out.genomicsdb, reference, reference_fai, reference_dict)

    emit:
    vcf_by_population = GATK4_GENOTYPEGVCFS.out.vcf
}
