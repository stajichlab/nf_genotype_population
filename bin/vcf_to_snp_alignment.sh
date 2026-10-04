#!/usr/bin/env bash
# bin/vcf_to_snp_alignment.sh
#
# Build a SNP alignment (FASTA, one row per sample plus the reference) from a
# biallelic SNP VCF, for a strain tree. Port of 06_make_SNP_tree.sh from the
# bash PopGenomics pipeline.
#   - keeps sites with F_MISSING <= MAX_MISSING (computed over the VCF samples)
#   - drops sites where every called allele is ALT (the reference row would be
#     the only difference)
#   - haploid calls give one base; diploid hets give an IUPAC ambiguity code;
#     missing calls give N
#   - writes no alignment (exit 0, with a warning) if no site is left or the
#     VCF has fewer than 2 samples: a tree needs at least 3 sequences
#
# Usage:
#   vcf_to_snp_alignment.sh -i snps.vcf.gz -o prefix [-r reference_name] [--max-missing 0] [-t threads]
# Writes <prefix>.mfa.gz and <prefix>.alignment_stats.tsv
set -euo pipefail

MAX_MISSING=0
REF_NAME=reference
THREADS=1
IN=""
PREFIX=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) IN=$2; shift 2 ;;
        -o) PREFIX=$2; shift 2 ;;
        -r) REF_NAME=$2; shift 2 ;;
        -t) THREADS=$2; shift 2 ;;
        --max-missing) MAX_MISSING=$2; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$IN" && -n "$PREFIX" ]] || { echo "need -i and -o" >&2; exit 2; }

TMP="${SCRATCH:-$PWD}/aln_tmp.$$"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

bcftools view --threads "$THREADS" -m2 -M2 -v snps \
    -i "F_MISSING <= ${MAX_MISSING}" -Ou "$IN" \
  | bcftools view --threads "$THREADS" -e 'AC==AN' -Ob -o "$TMP/sites.bcf"

NSITES=$(bcftools view -H "$TMP/sites.bcf" | wc -l)
NSAMP=$(bcftools query -l "$TMP/sites.bcf" | wc -l)
printf "max_missing\tsamples\tsites\n%s\t%s\t%s\n" "$MAX_MISSING" "$NSAMP" "$NSITES" > "${PREFIX}.alignment_stats.tsv"
if [[ "$NSITES" -eq 0 || "$NSAMP" -lt 2 ]]; then
    echo "WARNING: ${NSAMP} sample(s), ${NSITES} site(s) after F_MISSING <= ${MAX_MISSING}: no alignment written" >&2
    exit 0
fi

# One pass over the VCF; transpose site rows into sequence rows in awk.
{ printf "%s" "$REF_NAME"; bcftools query -l "$TMP/sites.bcf" | awk '{printf "\t%s", $1}'; printf "\n"
  bcftools query -f '%REF[\t%IUPACGT]\n' "$TMP/sites.bcf"
} | awk -F'\t' '
    NR==1 { n=NF; for (i=1;i<=n;i++) name[i]=$i; next }
    { for (i=1;i<=n;i++) { b=$i; if (b=="." || b=="./." || b==".|.") b="N"; seq[i]=seq[i] b } }
    END { for (i=1;i<=n;i++) printf ">%s\n%s\n", name[i], seq[i] }
  ' | gzip -c > "${PREFIX}.mfa.gz"

# Every row must have NSITES columns; IUPACGT can emit more than one character
# only for a malformed genotype, so check rather than assume.
BAD=$(gzip -dc "${PREFIX}.mfa.gz" | awk -v L="$NSITES" '!/^>/ && length($0)!=L' | wc -l)
if [[ "$BAD" -ne 0 ]]; then
    echo "${BAD} alignment rows do not have ${NSITES} columns" >&2
    exit 1
fi
