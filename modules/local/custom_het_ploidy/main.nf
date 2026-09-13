// modules/local/custom_het_ploidy/main.nf
process CUSTOM_HET_PLOIDY {
    tag "$strain"
    label 'process_low'
    // Local project-provided container (bcftools 1.24, samtools 1.24, python 3.14.7).
    // Source: ghcr.io/hyphaltip/container_bcftools_samtools/bcftools_samtools:1.24
    // Referenced here as a direct local .sif file path (pre-pulled to the shared cache)
    // so no docker:// pull is ever attempted at pipeline runtime (sidesteps a
    // singularity-ce 3.9.3 pull bug hit while testing an earlier mulled-container
    // reference, see Task 8's report).
    container '/bigdata/stajichlab/shared/singularity_cache/bcftools_samtools-1.24.sif'

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
    """
    export PATH="/opt/conda/envs/bcftools_samtools/bin:\$PATH"
    ${moduleDir}/../../../bin/custom_het_ploidy.py \\
        --cram ${cram} \\
        --reference ${reference} \\
        --strain ${strain} \\
        --out ${strain}.ploidy_inference.csv
    """
}
