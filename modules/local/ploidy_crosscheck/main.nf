// modules/local/ploidy_crosscheck/main.nf
process PLOIDY_CROSSCHECK {
    // FINAL-REVIEW M5: this process previously declared no `tag`, unlike every
    // other process in the branch, and lacked the explanatory container comment
    // header the others share. It consumes the single collected
    // ploidy_inference_all.csv (not one strain), so there is no per-strain value
    // to tag with; tag on the population-wide artifact it produces instead, so
    // traces and logs name it consistently with the rest of the pipeline.
    tag "ploidy_review"
    label 'process_low'
    // FINAL-REVIEW I1/I2: this was the one process in the branch still pinned to
    // a mutable floating Docker tag ('python:3.11-slim'), which under
    // singularity.enabled = true resolves to a docker:// pull at runtime - the
    // exact failure mode every other module carries comments explaining it
    // avoids (a singularity-ce 3.9.3 pull bug, Task 8), and an unpinned tag whose
    // Python minor version could change under the pipeline with no commit.
    //
    // It now uses the same pre-pulled local image as CUSTOM_HET_PLOIDY
    // (bcftools 1.24 / samtools 1.24 / python 3.14.7; upstream
    // docker://ghcr.io/hyphaltip/container_bcftools_samtools/bcftools_samtools:1.24),
    // declared centrally in conf/modules.config per spec §5.
    // bin/ploidy_crosscheck.py imports only the Python standard library
    // (argparse, csv), so that image's python3 is sufficient and no new image
    // has to be built or maintained.

    input:
    path inferred_csv
    path metadata_txt

    output:
    path 'ploidy_review.csv', emit: review

    script:
    // Same login-shell PATH caveat as CUSTOM_HET_PLOIDY: this project's
    // process.shell = ['/bin/bash', '-l'] makes the container run its own
    // /etc/profile, which resets PATH and hides the conda env's python3.
    // Prepend the env's bin dir explicitly rather than depending on shell-login
    // ordering.
    //
    // FINAL-REVIEW I7: see CUSTOM_HET_PLOIDY for why the bin directory is
    // resolved with a projectDir-first, moduleDir-fallback form.
    def bin_dir = file("${projectDir}/bin").exists() ? "${projectDir}/bin" : "${moduleDir}/../../../bin"
    """
    set -euo pipefail
    export PATH="/opt/conda/envs/bcftools_samtools/bin:\$PATH"
    ${bin_dir}/ploidy_crosscheck.py \\
        --inferred "${inferred_csv}" \\
        --metadata "${metadata_txt}" \\
        --out ploidy_review.csv
    """
}
