// modules/local/cram_to_v30/main.nf
//
// Rewrite a CRAM as CRAM 3.0. GATK 4.5 (htsjdk) cannot read CRAM 3.1
// ("CRAM version 3.1 is not supported"), which samtools >= 1.15 writes by
// default (for example `samtools merge` / `samtools view -C`). Without this
// step such a CRAM fails only when HaplotypeCaller reaches it, hours into a run.
// PREPARE_CRAMS sends only non-3.0 CRAMs here.
process CRAM_TO_V30 {
    tag "$strain"
    label 'process_low'
    // Container (bcftools 1.24 / samtools 1.24) is declared in conf/modules.config.

    input:
    tuple val(strain), path(cram), path(crai)
    path reference
    path reference_fai

    output:
    tuple val(strain), path("v30/${strain}.cram"), path("v30/${strain}.cram.crai"), emit: cram

    script:
    """
    set -euo pipefail
    export PATH="/opt/conda/envs/bcftools_samtools/bin:\$PATH"
    mkdir -p v30
    samtools view -@ ${task.cpus} -C --output-fmt-option version=3.0 -T "${reference}" -o "v30/${strain}.cram" "${cram}"
    samtools index "v30/${strain}.cram"
    a=\$(samtools view -c -@ ${task.cpus} --reference "${reference}" "${cram}")
    b=\$(samtools view -c -@ ${task.cpus} --reference "${reference}" "v30/${strain}.cram")
    if [[ "\$a" -ne "\$b" ]]; then
        echo "CRAM_TO_V30 ${strain}: record count changed \$a -> \$b" >&2
        exit 1
    fi
    """
}
