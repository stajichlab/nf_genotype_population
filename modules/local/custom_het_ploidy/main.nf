// modules/local/custom_het_ploidy/main.nf
process CUSTOM_HET_PLOIDY {
    tag "$strain"
    label 'process_low'
    // FINAL-REVIEW I2: container (bcftools 1.24, samtools 1.24, python 3.14.7)
    // declared centrally in conf/modules.config (spec §5), which records the
    // upstream source
    // docker://ghcr.io/hyphaltip/container_bcftools_samtools/bcftools_samtools:1.24
    // and why a pre-pulled local image is used instead of a runtime docker://
    // pull (sidesteps a singularity-ce 3.9.3 pull bug hit in Task 8).

    input:
    tuple val(strain), path(cram), path(crai)
    path reference
    path reference_fai

    output:
    tuple val(strain), path("${strain}.ploidy_inference.csv"), emit: csv

    script:
    // This project's process.shell = ['/bin/bash', '-l'] (needed elsewhere for HPCC
    // module-loading) invokes a LOGIN shell inside the container too. This container's
    // own /etc/profile resets PATH to a bare system default in login-shell mode,
    // clobbering the conda env (bcftools/samtools/python 3.14) that Singularity's own
    // env-file sourcing had already set up correctly for non-login invocations
    // (verified by direct reproduction: PATH differs between `bash -c` and `bash -lc`
    // inside this exact container). Prepend the conda env's bin dir explicitly so this
    // process doesn't depend on shell-login ordering at all.
    //
    // FINAL-REVIEW I7: prefer ${projectDir}/bin, fall back to the module-relative
    // path.
    //
    // projectDir is the right anchor for a real pipeline run (it resolves to the
    // directory of the entry script, i.e. the repo root for main.nf) and does not
    // depend on how deep this module sits in the tree. But projectDir is
    // entrypoint-dependent: it was MEASURED during this fix round that running
    // `nextflow run tests/smoke_custom_het_ploidy.nf` sets projectDir to
    // <repo>/tests, where no bin/ exists - which is also why the smoke tests
    // address their fixtures as ${projectDir}/../tests/fixtures. A bare
    // ${projectDir}/bin form therefore breaks every module-level smoke test.
    //
    // So: use projectDir when it actually contains bin/ (the production case),
    // and otherwise fall back to this module's own known location in the tree
    // (the test-entrypoint case). Both branches were exercised and verified.
    def bin_dir = file("${projectDir}/bin").exists() ? "${projectDir}/bin" : "${moduleDir}/../../../bin"
    """
    set -euo pipefail
    export PATH="/opt/conda/envs/bcftools_samtools/bin:\$PATH"
    ${bin_dir}/custom_het_ploidy.py \\
        --cram "${cram}" \\
        --reference "${reference}" \\
        --strain "${strain}" \\
        --out "${strain}.ploidy_inference.csv"
    """
}
