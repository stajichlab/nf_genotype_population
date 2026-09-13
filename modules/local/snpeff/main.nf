// modules/local/snpeff/main.nf
process SNPEFF_ANNOTATE {
    tag "$population"
    label 'process_low'
    // Local project-provided container: mulled combo of snpeff=5.2,tabix=1.11
    // from Galaxy's depot (this is NOT the plain snpeff:5.2--hdfd78af_1
    // biocontainer - that image was checked and has neither bgzip nor tabix,
    // and this module's script pipes snpEff output through bgzip then runs
    // tabix). Pre-pulled to the shared cache and confirmed to have snpEff
    // 5.2, bgzip, and tabix all present and working under this project's
    // login-shell process.shell.
    // Source: docker://quay.io/biocontainers/mulled-v2-2fe536b56916bd1d61a6a1889eb2987d9ea0cd2f:c51b2e46bf63786b2d9a7a7d23680791163ab39a-0
    // (Galaxy depot equivalent: depot.galaxyproject.org/singularity/mulled-v2-2fe536b56916bd1d61a6a1889eb2987d9ea0cd2f:c51b2e46bf63786b2d9a7a7d23680791163ab39a-0)
    container '/bigdata/stajichlab/shared/lib/singularity_cache/depot.galaxyproject.org-singularity-mulled-v2-2fe536b56916bd1d61a6a1889eb2987d9ea0cd2f-c51b2e46bf63786b2d9a7a7d23680791163ab39a-0.img'

    input:
    tuple val(population), path(vcf), path(tbi)
    path snpeff_db_dir
    val snpeff_genome_name

    output:
    tuple val(population), path("${population}.annotated.vcf.gz"), emit: vcf

    script:
    // snpeff_db_dir is the database root as built per conf/snpeff.config's
    // recipe: it contains a snpEff.config (base config + this genome's
    // custom entry, see conf/snpeff.config) at its root, and a data/
    // subdirectory holding data/<snpeff_genome_name>/{sequences.fa,genes.gtf,...}.
    // -c is required: without it snpEff falls back to the container's own
    // bundled default config, which has no entry for our custom genome name
    // and fails with "Property: '<genome>.genome' not found" (confirmed
    // while building the Task 12 smoke test).
    //
    // Nextflow stages `snpeff_db_dir` into the work dir as a symlink named
    // after the input directory's basename (pointing at the real path
    // elsewhere on shared storage). snpEff's database lookup does not
    // tolerate being pointed at that path via a relative symlink - it
    // reports "FATAL ERROR: Failed to download database ..." instead of
    // using the local files, even though the exact same files work fine
    // when addressed by their canonical absolute path (confirmed while
    // debugging the Task 12 smoke test). Resolve to an absolute, symlink-free
    // path with `readlink -f` before invoking snpEff.
    """
    db_dir=\$(readlink -f ${snpeff_db_dir})
    snpEff -c \${db_dir}/snpEff.config -dataDir \${db_dir}/data ${snpeff_genome_name} ${vcf} | bgzip > ${population}.annotated.vcf.gz
    tabix -p vcf ${population}.annotated.vcf.gz
    """
}
