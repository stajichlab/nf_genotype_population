// modules/local/snpeff/main.nf
process SNPEFF_ANNOTATE {
    tag "$population"
    label 'process_low'
    // FINAL-REVIEW I2: the container is now declared centrally in
    // conf/modules.config (spec §5). It is a mulled combo of snpeff=5.2 and
    // tabix=1.11 from Galaxy's depot - NOT the plain snpeff:5.2--hdfd78af_1
    // biocontainer, which was checked and has neither bgzip nor tabix, both of
    // which this module's script needs (it pipes snpEff output through bgzip
    // and then runs tabix). conf/modules.config records the upstream source.

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
    set -euo pipefail
    db_dir=\$(readlink -f "${snpeff_db_dir}")
    snpEff -c "\${db_dir}/snpEff.config" -dataDir "\${db_dir}/data" "${snpeff_genome_name}" "${vcf}" | bgzip > "${population}.annotated.vcf.gz"
    tabix -p vcf "${population}.annotated.vcf.gz"
    """
}
