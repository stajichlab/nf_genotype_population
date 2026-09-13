process PLOIDY_CROSSCHECK {
    label 'process_low'
    container 'python:3.11-slim'

    input:
    path inferred_csv
    path metadata_txt

    output:
    path 'ploidy_review.csv', emit: review

    script:
    """
    ${moduleDir}/../../../bin/ploidy_crosscheck.py --inferred ${inferred_csv} --metadata ${metadata_txt} --out ploidy_review.csv
    """
}
