#!/usr/bin/env bash
# tests/run_integration_test.sh
#
# FINAL-REVIEW I8(a) + I9.
#
# End-to-end integration test for the assembled pipeline. Before this script,
# NOTHING tested workflows/genotype_population.nf or either main.nf entry
# workflow: commit a58a51b claimed a full synthetic integration run, but no test
# file was committed, so the assembled system - buildCramCh()'s glob, the
# failOnMismatch join, the .fai interval derivation, the publish subscribes -
# was untested and unreproducible. This script does everything that run did by
# hand, and then asserts on the result instead of asking a human to read stdout.
#
# What it does:
#   1. builds a throwaway flat CRAM dir from the committed synthetic fixtures
#   2. builds a throwaway SnpEff database from those fixtures
#   3. runs `-entry PLOIDY_ONLY`
#   4. derives ploidy_overrides.csv from its output using EXACTLY the one-liner
#      the README documents (so the README's instructions are themselves tested)
#   5. runs the full default workflow
#   6. asserts on results/all.qc.annotated.vcf.gz and
#      results/all.snps.maf.annotated.vcf.gz (VARIANT_QC_FILTER outputs), the
#      filter stats, the published GVCFs, and the published ploidy QC record
#
# Every assertion failure exits non-zero. Expect roughly 10-20 minutes, mostly
# SLURM queueing behind other jobs on this account.
#
# Run from the repo root:
#   tests/run_integration_test.sh
#
# NOTE: no BASH_SOURCE-based path resolution anywhere here - it resolves to the
# wrong thing under SLURM work directories. The repo root comes from git, or
# from NF_GENOTYPE_POPULATION_ROOT if you set it.

set -euo pipefail

REPO_ROOT="${NF_GENOTYPE_POPULATION_ROOT:-$(git rev-parse --show-toplevel)}"
cd "${REPO_ROOT}"

FIXTURES="${REPO_ROOT}/tests/fixtures"
SCRATCH_DIR="${REPO_ROOT}/.integration_test"
CRAM_DIR="${SCRATCH_DIR}/cram"
DB_DIR="${SCRATCH_DIR}/snpeff_db"
OUTDIR="${SCRATCH_DIR}/results"
OVERRIDES="${SCRATCH_DIR}/ploidy_overrides.csv"
GENOME="synth_genome"

SNPEFF_CONTAINER="${SNPEFF_CONTAINER:-/bigdata/stajichlab/shared/singularity_cache/depot.galaxyproject.org-singularity-mulled-v2-2fe536b56916bd1d61a6a1889eb2987d9ea0cd2f-c51b2e46bf63786b2d9a7a7d23680791163ab39a-0.img}"

FAILURES=0

pass() { echo "  PASS: $*"; }
fail() { echo "  FAIL: $*"; FAILURES=$((FAILURES + 1)); }

section() { echo; echo "=== $* ==="; }

cleanup_previous() {
    rm -rf "${SCRATCH_DIR}"
    # A stale work/ from an earlier run would let -resume serve cached results
    # and hide a real regression, so start clean.
    rm -rf "${REPO_ROOT}/work" "${REPO_ROOT}/.nextflow" "${REPO_ROOT}"/.nextflow.log*
    mkdir -p "${CRAM_DIR}" "${SCRATCH_DIR}"
}

build_cram_dir() {
    section "Building throwaway flat CRAM directory"
    for strain in haploid_strain diploid_strain; do
        cp "${FIXTURES}/${strain}.cram"      "${CRAM_DIR}/${strain}.cram"
        cp "${FIXTURES}/${strain}.cram.crai" "${CRAM_DIR}/${strain}.cram.crai"
    done
    ls -1 "${CRAM_DIR}"
}

build_snpeff_db() {
    section "Building throwaway SnpEff database (${GENOME})"
    module load singularity 2>/dev/null || true
    mkdir -p "${DB_DIR}/data/${GENOME}"
    cp "${FIXTURES}/ref.fa"         "${DB_DIR}/data/${GENOME}/sequences.fa"
    cp "${FIXTURES}/synth_ref.gtf"  "${DB_DIR}/data/${GENOME}/genes.gtf"
    singularity exec "${SNPEFF_CONTAINER}" cat /usr/local/share/snpeff-5.2-0/snpEff.config > "${DB_DIR}/snpEff.config"
    echo "${GENOME}.genome : synthetic integration test genome" >> "${DB_DIR}/snpEff.config"
    singularity exec "${SNPEFF_CONTAINER}" snpEff build \
        -gtf22 -dataDir "${DB_DIR}/data" -c "${DB_DIR}/snpEff.config" \
        -noCheckCds -noCheckProtein -v "${GENOME}"
    test -s "${DB_DIR}/data/${GENOME}/snpEffectPredictor.bin"
    echo "SnpEff DB built."
}

run_ploidy_only() {
    section "Phase 1: nextflow run main.nf -entry PLOIDY_ONLY"
    # NXF_SYNTAX_PARSER=v1 is required on this cluster's Nextflow (26.04.6):
    # its default strict parser rejects -entry outright.
    NXF_SYNTAX_PARSER=v1 nextflow run main.nf -entry PLOIDY_ONLY -profile hpcc \
        --cram_dir      "${CRAM_DIR}" \
        --reference     "${FIXTURES}/ref.fa" \
        --reference_fai "${FIXTURES}/ref.fa.fai" \
        --metadata      "${FIXTURES}/metadata_fixture.txt" \
        --outdir        "${OUTDIR}"
}

build_overrides() {
    section "Deriving ploidy_overrides.csv from the review CSV"
    # This is EXACTLY the command the README tells a user to run. Keeping it
    # identical means a drift between the docs and the pipeline breaks this test.
    { echo "strain,ploidy"; tail -n +2 "${OUTDIR}/ploidy_review.csv" | cut -d, -f1,2; } > "${OVERRIDES}"
    cat "${OVERRIDES}"

    # The README's own rule: one row per strain with a CRAM.
    local n_crams n_rows
    n_crams=$(find "${CRAM_DIR}" -name '*.cram' | wc -l)
    n_rows=$(($(wc -l < "${OVERRIDES}") - 1))
    if [[ "${n_crams}" -ne "${n_rows}" ]]; then
        echo "FATAL: ${n_crams} CRAMs but ${n_rows} overrides rows; the README derivation is wrong." >&2
        exit 1
    fi
    if grep -q ',unknown$' "${OVERRIDES}"; then
        echo "FATAL: overrides file contains an 'unknown' ploidy; hand-edit required (see README)." >&2
        exit 1
    fi
}

run_full_workflow() {
    section "Phase 2: nextflow run main.nf (default workflow)"
    nextflow run main.nf -profile hpcc \
        --cram_dir           "${CRAM_DIR}" \
        --reference          "${FIXTURES}/ref.fa" \
        --reference_fai      "${FIXTURES}/ref.fa.fai" \
        --reference_dict     "${FIXTURES}/ref.dict" \
        --metadata           "${FIXTURES}/metadata_fixture.txt" \
        --ploidy_overrides   "${OVERRIDES}" \
        --population_sets    "${FIXTURES}/population_sets_fixture.yaml" \
        --snpeff_db_dir      "${DB_DIR}" \
        --snpeff_genome_name "${GENOME}" \
        --mask_bed           "${FIXTURES}/mask_fixture.bed" \
        --outdir             "${OUTDIR}"
}

assert_outputs() {
    local vcf="${OUTDIR}/all.qc.annotated.vcf.gz"

    section "Assertions on ${vcf}"

    if [[ ! -s "${vcf}" ]]; then
        echo "  FAIL: ${vcf} missing or empty" >&2
        exit 1
    fi
    pass "final annotated VCF exists and is non-empty"

    if [[ ! -s "${vcf}.tbi" ]]; then
        echo "  FAIL: ${vcf}.tbi missing or empty" >&2
        exit 1
    fi
    pass "final annotated VCF index (.tbi) published"

    local body
    body=$(zcat "${vcf}" | grep -v '^#')

    # ---- record count -----------------------------------------------------
    # 34 joint-called records -> 27 PASS -> 24 outside mask_fixture.bed
    # (522, 565, 590 masked) -> 22 after missingness: at 75 and 150 the haploid
    # has FORMAT/DP=3 < qc_min_dp, so its GT is set missing and F_MISSING=0.5.
    local n_records
    n_records=$(printf '%s\n' "${body}" | grep -c . || true)
    if [[ "${n_records}" -eq 22 ]]; then
        pass "record count is 22"
    else
        fail "record count is ${n_records}, expected 22"
    fi

    # ---- ANN= on every record --------------------------------------------
    local n_ann
    n_ann=$(printf '%s\n' "${body}" | awk -F'\t' '$8 ~ /ANN=/' | grep -c . || true)
    if [[ "${n_ann}" -eq "${n_records}" ]]; then
        pass "all ${n_records} records carry an ANN= annotation"
    else
        fail "only ${n_ann} of ${n_records} records carry ANN="
    fi

    # ---- QC filter effects ------------------------------------------------
    local n_notpass
    n_notpass=$(printf '%s\n' "${body}" | awk -F'\t' '$7 != "PASS"' | grep -c . || true)
    if [[ "${n_notpass}" -eq 0 ]]; then
        pass "every record is FILTER=PASS"
    else
        fail "${n_notpass} record(s) are not PASS after VARIANT_QC_FILTER"
    fi
    local n_masked
    n_masked=$(printf '%s\n' "${body}" | awk -F'\t' '$2 > 500 && $2 <= 600' | grep -c . || true)
    if [[ "${n_masked}" -eq 0 ]]; then
        pass "no record inside mask_fixture.bed"
    else
        fail "${n_masked} record(s) inside mask_fixture.bed"
    fi
    local n_lowgt
    n_lowgt=$(printf '%s\n' "${body}" | awk -F'\t' '
        { split($9, f, ":"); for (i in f) { if (f[i]=="GQ") gq=i; if (f[i]=="DP") dp=i }
          for (c = 10; c <= NF; c++) { split($c, g, ":");
            if (g[1] !~ /\./ && ((g[gq] != "." && g[gq] < 20) || (g[dp] != "." && g[dp] < 5))) print } }' | grep -c . || true)
    if [[ "${n_lowgt}" -eq 0 ]]; then
        pass "no called genotype has GQ < 20 or DP < 5"
    else
        fail "${n_lowgt} called genotype(s) have GQ < 20 or DP < 5"
    fi

    local stats="${OUTDIR}/all.filter_stats.tsv"
    if [[ -s "${stats}" ]] && awk -F'\t' '$1=="input"{i=$2} $1=="pass_unmasked"{p=$2} END{exit !(i==34 && p<i)}' "${stats}"; then
        pass "filter stats published; hard filters flagged records (input 34 > pass_unmasked)"
    else
        fail "filter stats missing or hard filters flagged nothing: $(tr '\n' ' ' < "${stats}" 2>/dev/null)"
    fi

    local snps="${OUTDIR}/all.snps.maf.annotated.vcf.gz"
    local n_snps
    n_snps=$(zcat "${snps}" 2>/dev/null | grep -vc '^#' || true)
    if [[ -s "${snps}.tbi" && "${n_snps}" -eq 22 ]]; then
        pass "biallelic SNP set published with index, 22 records"
    else
        fail "biallelic SNP set: ${n_snps} records (expected 22) or index missing"
    fi

    # ---- both sample columns present --------------------------------------
    local header samples
    header=$(zcat "${vcf}" | grep '^#CHROM')
    samples=$(printf '%s\n' "${header}" | cut -f10- | tr '\t' ' ')
    local hap_col dip_col
    hap_col=$(printf '%s\n' "${header}" | tr '\t' '\n' | grep -n '^haploid_strain$' | cut -d: -f1 || true)
    dip_col=$(printf '%s\n' "${header}" | tr '\t' '\n' | grep -n '^diploid_strain$' | cut -d: -f1 || true)
    if [[ -n "${hap_col}" && -n "${dip_col}" ]]; then
        pass "both sample columns present: ${samples}"
    else
        echo "  FAIL: sample columns are '${samples}', expected both fixture strains" >&2
        exit 1
    fi

    # ---- THE central assertion (FINAL-REVIEW I9) --------------------------
    # Ploidy-aware per-strain calling surviving a mixed-ploidy GenomicsDB merge
    # is this pipeline's central scientific claim (spec §2), and nothing verified
    # it at the output level. At a site where BOTH samples are genotyped, the
    # haploid sample's GT must have exactly ONE allele and the diploid sample's
    # GT must have TWO.
    section "Assertion: mixed-ploidy genotypes (FINAL-REVIEW I9)"
    local gt_result
    gt_result=$(printf '%s\n' "${body}" | awk -F'\t' -v hap="${hap_col}" -v dip="${dip_col}" '
        {
            split($hap, h, ":"); split($dip, d, ":");
            hgt = h[1]; dgt = d[1];
            if (hgt == "./." || hgt == "." || dgt == "./." || dgt == ".") next;
            checked++;
            if (hgt !~ /^[0-9]+$/) { print "BAD_HAP\t" $1 ":" $2 "\t" hgt; bad++ }
            else if (dgt !~ /^[0-9]+[\/|][0-9]+$/) { print "BAD_DIP\t" $1 ":" $2 "\t" dgt; bad++ }
        }
        END { print "CHECKED\t" checked+0 "\tBAD\t" bad+0 }
    ')
    printf '%s\n' "${gt_result}" | grep -E '^BAD_' || true

    local n_checked n_bad
    n_checked=$(printf '%s\n' "${gt_result}" | awk -F'\t' '$1=="CHECKED" {print $2}')
    n_bad=$(printf '%s\n' "${gt_result}" | awk -F'\t' '$1=="CHECKED" {print $4}')

    if [[ "${n_checked}" -eq 0 ]]; then
        fail "no site had both samples genotyped - the mixed-ploidy assertion never ran"
    elif [[ "${n_bad}" -ne 0 ]]; then
        fail "${n_bad} of ${n_checked} shared site(s) have the wrong genotype arity (see BAD_ lines above)"
    else
        pass "at all ${n_checked} shared site(s): haploid_strain GT has exactly 1 allele, diploid_strain GT has 2"
        printf '%s\n' "${body}" | awk -F'\t' -v hap="${hap_col}" -v dip="${dip_col}" '
            { split($hap,h,":"); split($dip,d,":");
              if (h[1] != "./." && h[1] != "." && d[1] != "./." && d[1] != ".") {
                  print "        " $1 ":" $2 "  haploid_strain GT=" h[1] "   diploid_strain GT=" d[1]; n++ }
              if (n >= 5) exit }'
    fi

    # ---- published GVCFs (FINAL-REVIEW I10) -------------------------------
    section "Assertions on published GVCFs (FINAL-REVIEW I10)"
    for strain in haploid_strain diploid_strain; do
        if [[ -s "${OUTDIR}/gvcfs/${strain}.g.vcf.gz" && -s "${OUTDIR}/gvcfs/${strain}.g.vcf.gz.tbi" ]]; then
            pass "${strain}: GVCF and index retained in results/gvcfs/"
        else
            fail "${strain}: GVCF or index not published to ${OUTDIR}/gvcfs/"
        fi
    done

    # ---- published ploidy QC record (FINAL-REVIEW I6) ---------------------
    section "Assertions on the ploidy QC record (FINAL-REVIEW I6)"
    if [[ -s "${OUTDIR}/ploidy_review.csv" ]]; then
        pass "ploidy_review.csv published by the DEFAULT workflow"
    else
        fail "ploidy_review.csv not published by the default workflow"
    fi
    if [[ -s "${OUTDIR}/ploidy_inference_all.csv" ]]; then
        pass "ploidy_inference_all.csv published by the DEFAULT workflow"
    else
        fail "ploidy_inference_all.csv not published by the default workflow"
    fi
}

main() {
    cleanup_previous
    build_cram_dir
    build_snpeff_db
    run_ploidy_only
    build_overrides
    run_full_workflow
    assert_outputs

    section "RESULT"
    if [[ "${FAILURES}" -eq 0 ]]; then
        echo "INTEGRATION TEST PASSED - all assertions held."
        exit 0
    fi
    echo "INTEGRATION TEST FAILED - ${FAILURES} assertion(s) failed." >&2
    exit 1
}

main "$@"
