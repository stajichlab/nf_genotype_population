// modules/local/variant_qc_filter/main.nf
process VARIANT_QC_FILTER {
    tag "$population"
    label 'process_medium'
    // Container (bcftools 1.24 with plugins) is declared in conf/modules.config.

    input:
    tuple val(population), path(vcf), path(tbi), val(strains)  // strains: [] = every sample in the VCF
    path mask_bed   // repeat + low-complexity BED, or assets/NO_FILE to skip masking

    output:
    tuple val("${population}.qc"),       path("${population}.qc.vcf.gz"),       path("${population}.qc.vcf.gz.tbi"),       emit: qc
    tuple val("${population}.snps.maf"), path("${population}.snps.maf.vcf.gz"), path("${population}.snps.maf.vcf.gz.tbi"), emit: snps
    tuple val(population), path("${population}.filter_stats.tsv"), path("${population}.depth_threshold.txt"),             emit: stats

    script:
    // Genotype- and site-level QC after GATK hard filtering; see
    // bin/filter_population_vcf.sh for the full step list. Thresholds are
    // params so a population can be re-filtered without re-genotyping.
    def mask_arg = mask_bed.name == 'NO_FILE' ? '' : "-m ${mask_bed}"
    def strains_cmd = strains ? "printf '%s\\n' ${strains.collect { "'${it}'" }.join(' ')} > strains.txt" : ''
    def strains_arg = strains ? '--strains strains.txt' : ''
    // bin_dir and the PATH export follow modules/local/custom_het_ploidy: the
    // login shell (process.shell = bash -l) drops the image's conda env from
    // PATH, and projectDir has no bin/ under the tests/ entry points.
    def bin_dir = file("${projectDir}/bin").exists() ? "${projectDir}/bin" : "${moduleDir}/../../../bin"
    """
    set -euo pipefail
    export PATH="/opt/conda/envs/bcftools_samtools/bin:\$PATH"
    ${strains_cmd}
    ${bin_dir}/filter_population_vcf.sh \\
        -i "${vcf}" \\
        -o "${population}" \\
        ${mask_arg} \\
        ${strains_arg} \\
        -t ${task.cpus} \\
        --min-gq ${params.qc_min_gq} \\
        --min-dp ${params.qc_min_dp} \\
        --min-ab ${params.qc_min_ab} \\
        --hap-min-af ${params.qc_hap_min_af} \\
        --dp-max-factor ${params.qc_dp_max_factor} \\
        --max-missing ${params.qc_max_missing} \\
        --min-maf ${params.qc_min_maf}
    """
}
