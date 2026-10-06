// modules/local/snp_tree/main.nf
//
// Strain tree from the biallelic SNP set (<pop>.snps.maf). Port of
// 06_make_SNP_tree.sh / 07_iqtree.sh from the bash PopGenomics pipeline:
// SNP alignment with the reference as one row (diploid hets as IUPAC codes),
// then IQ-TREE with an ascertainment-bias model (default GTR+ASC), ultrafast
// bootstrap and SH-aLRT.

process SNP_ALIGNMENT {
    tag "$population"
    label 'process_low'
    // Container (bcftools 1.24) is declared in conf/modules.config.

    input:
    tuple val(population), path(vcf), path(tbi)

    output:
    tuple val(population), path("${population}.mfa.gz"), optional: true, emit: alignment  // absent if < 2 samples or 0 sites
    tuple val(population), path("${population}.alignment_stats.tsv"),    emit: stats

    script:
    // bin_dir and the PATH export follow modules/local/variant_qc_filter.
    def bin_dir = file("${projectDir}/bin").exists() ? "${projectDir}/bin" : "${moduleDir}/../../../bin"
    """
    set -euo pipefail
    export PATH="/opt/conda/envs/bcftools_samtools/bin:\$PATH"
    ${bin_dir}/vcf_to_snp_alignment.sh \\
        -i "${vcf}" \\
        -o "${population}" \\
        -r "${params.tree_ref_name}" \\
        -t ${task.cpus} \\
        --max-missing ${params.tree_max_missing}
    """
}

process IQTREE {
    tag "$population"
    label 'process_medium'
    // Container (IQ-TREE 3) is declared in conf/modules.config.

    input:
    tuple val(population), path(aln)

    output:
    tuple val(population), path("${population}.treefile"), path("${population}.iqtree"), path("${population}.log"), emit: tree

    script:
    // tree_bootstraps = 0 turns off UFBoot and SH-aLRT (both need >= 4 sequences).
    def boot = params.tree_bootstraps as int > 0 ? "-B ${params.tree_bootstraps} -alrt ${params.tree_bootstraps}" : ''
    """
    set -euo pipefail
    gzip -dc "${aln}" > "${population}.mfa"
    run_iq() {
        iqtree3 -s "\$1" -st DNA --prefix "${population}" -m "${params.tree_model}" \\
            ${boot} -keep-ident \\
            -T ${task.cpus} --seed 12345 -redo
    }
    # +ASC stops if a column is invariant. IQ-TREE treats an IUPAC code as
    # compatible with its bases, so a column that varies only by diploid het
    # codes (for example A and R) counts as invariant. IQ-TREE then writes
    # <prefix>.varsites.phy with the variable columns only; rerun on that file.
    # -st DNA is required: IQ-TREE cannot detect the type of that PHYLIP file.
    # -keep-ident is required: without it IQ-TREE collapses identical sequences
    # before it writes varsites.phy, and the rerun drops those strains from the tree.
    if ! run_iq "${population}.mfa"; then
        [[ -s "${population}.varsites.phy" ]] || exit 1
        run_iq "${population}.varsites.phy"
    fi
    rm -f "${population}.mfa"
    """
}
