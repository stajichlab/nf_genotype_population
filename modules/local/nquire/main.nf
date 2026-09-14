// modules/local/nquire/main.nf
process NQUIRE {
    tag "$strain"
    label 'process_low'
    // Container built from source (see .build_tmp/nquire.def in this repo's
    // history / assets/nquire.def) - no upstream container or bioconda package
    // exists for nQuire (github.com/clwgg/nQuire). Built with apptainer
    // --fakeroot from git commit a990a88ef14b38f257f1a0d368ba8be1bd3d7e4b of
    // https://github.com/clwgg/nQuire, on an ubuntu:22.04 base, plus a system
    // samtools package (needed for the CRAM->BAM conversion bin/nquire_ploidy.py
    // does before invoking nQuire - see that script's docstring for why: nQuire's
    // own bundled htslib build explicitly disables full CRAM support).
    // Declared centrally in conf/modules.config (spec §5), matching every other
    // process in this pipeline.

    input:
    tuple val(strain), path(cram), path(crai)
    path reference
    path reference_fai

    output:
    tuple val(strain), path("${strain}.ploidy_inference.csv"), emit: csv

    script:
    // Same projectDir-first, moduleDir-fallback bin/ resolution as
    // CUSTOM_HET_PLOIDY (see that module for the full rationale - this is
    // required so tests/smoke_*.nf entrypoints, whose projectDir has no bin/,
    // still find the script).
    def bin_dir = file("${projectDir}/bin").exists() ? "${projectDir}/bin" : "${moduleDir}/../../../bin"
    """
    set -euo pipefail
    ${bin_dir}/nquire_ploidy.py \\
        --cram "${cram}" \\
        --reference "${reference}" \\
        --strain "${strain}" \\
        --nquire-bin nQuire \\
        --out "${strain}.ploidy_inference.csv"
    """
}
