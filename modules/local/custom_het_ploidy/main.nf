// modules/local/custom_het_ploidy/main.nf
process CUSTOM_HET_PLOIDY {
    tag "$strain"
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-b24619740583784f45c6a755710ce30d3e160639:1f6b1b8b7a2df34bea1c43f12a5c4db9c9f10fc3-0'

    input:
    tuple val(strain), path(cram), path(crai)
    path reference
    path reference_fai

    output:
    tuple val(strain), path("${strain}.ploidy_inference.csv"), emit: csv

    script:
    """
    ${moduleDir}/../../../bin/custom_het_ploidy.py \\
        --cram ${cram} \\
        --reference ${reference} \\
        --strain ${strain} \\
        --out ${strain}.ploidy_inference.csv
    """
}
