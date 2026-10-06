#!/usr/bin/env bash
# bin/filter_population_vcf.sh
#
# Genotype- and site-level QC for a hard-filtered, mixed-ploidy population VCF.
# Runs after GATK4_HARDFILTER and before SNPEFF_ANNOTATE.
#
# Steps (counts of records after every step go to <prefix>.filter_stats.tsv):
#   1. optional sample subset (-S VCF sample names, or --strains strain names;
#      a strain matches sample STRAIN or STRAIN_STRAIN); keep FILTER=PASS sites (GATK hard
#      filters); recompute INFO/DP as the sum of FORMAT/DP over the kept samples
#   2. drop sites whose REF span overlaps the repeat / low-complexity mask BED
#   3. drop sites with INFO/DP > DP_MAX_FACTOR x median INFO/DP (collapsed
#      repeats / CNV)
#   4. set genotypes to missing where FORMAT/GQ < MIN_GQ, or FORMAT/DP < MIN_DP
#      for a genotype that carries an ALT allele, or FORMAT/DP < MIN_DP_HOMREF
#      for a hom-ref genotype (GATK genotype-level filters, then
#      --set-filtered-gt-to-nocall equivalent). Hom-ref genotypes have their own
#      DP threshold (default 0 = no DP mask) because GenotypeGVCFs takes their
#      FORMAT/DP from the gVCF reference block, which is the block MIN_DP, not
#      the depth at the site. A DP mask on that value removes most hom-ref calls
#      of strains close to the reference. GQ < MIN_GQ still applies to them.
#   5. set heterozygous genotypes to missing where one allele carries more than
#      1-MIN_AB of the allelic reads (for a biallelic site: AB < MIN_AB or
#      AB > 1-MIN_AB). smpl_max/smpl_sum over FORMAT/AD works for any number of
#      alleles, unlike a fixed AD[:1] index.
#   6. set haploid and homozygous genotypes to missing where no single allele
#      carries at least HAP_MIN_AF of the allelic reads (mixed reads in a call
#      that should show one allele)
#   7. drop ALT alleles that no remaining genotype carries, drop sites with no
#      ALT allele left, recompute AC/AN/AF/MAF/F_MISSING (ploidy aware)
#   8. drop sites with F_MISSING > MAX_MISSING
#   -> <prefix>.qc.vcf.gz          all variant types, no MAF filter
#   9. biallelic SNPs only (no '*' allele), MAF >= MIN_MAF
#   -> <prefix>.snps.maf.vcf.gz    population-genetics SNP set
#
# Usage:
#   filter_population_vcf.sh -i in.vcf.gz -o prefix [-m mask.bed] [-S samples.txt | --strains strains.txt] [-t threads]
#     [--min-gq 20] [--min-dp 5] [--min-dp-homref 0] [--min-ab 0.2] [--hap-min-af 0.8]
#     [--dp-max-factor 2] [--max-missing 0.1] [--min-maf 0.05]
set -euo pipefail

MIN_GQ=20
MIN_DP=5
MIN_DP_HOMREF=0
MIN_AB=0.2
HAP_MIN_AF=0.8
DP_MAX_FACTOR=2
MAX_MISSING=0.1
MIN_MAF=0.05
THREADS=1
MASK=""
SAMPLES=""
STRAINS=""
IN=""
PREFIX=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) IN=$2; shift 2 ;;
        -o) PREFIX=$2; shift 2 ;;
        -m) MASK=$2; shift 2 ;;
        -S) SAMPLES=$2; shift 2 ;;
        --strains) STRAINS=$2; shift 2 ;;
        -t) THREADS=$2; shift 2 ;;
        --min-gq) MIN_GQ=$2; shift 2 ;;
        --min-dp) MIN_DP=$2; shift 2 ;;
        --min-dp-homref) MIN_DP_HOMREF=$2; shift 2 ;;
        --min-ab) MIN_AB=$2; shift 2 ;;
        --hap-min-af) HAP_MIN_AF=$2; shift 2 ;;
        --dp-max-factor) DP_MAX_FACTOR=$2; shift 2 ;;
        --max-missing) MAX_MISSING=$2; shift 2 ;;
        --min-maf) MIN_MAF=$2; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$IN" && -n "$PREFIX" ]] || { echo "need -i and -o" >&2; exit 2; }
[[ -n "$SAMPLES" && -n "$STRAINS" ]] && { echo "use -S or --strains, not both" >&2; exit 2; }

TMP="${SCRATCH:-$PWD}/filter_tmp.$$"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
STATS="${PREFIX}.filter_stats.tsv"
printf "step\trecords\n" > "$STATS"
count() { printf "%s\t%s\n" "$1" "$(bcftools index -n "$2")" >> "$STATS"; }

# Strain names -> VCF sample names. Sarek-derived CRAMs give sample STRAIN_STRAIN.
if [[ -n "$STRAINS" ]]; then
    bcftools query -l "$IN" > "$TMP/vcf_samples.txt"
    awk -v miss="$TMP/missing.txt" 'NR==FNR {s[$1]=1; next}
        NF { if ($1 in s) print $1; else if (($1 "_" $1) in s) print $1 "_" $1; else print $1 > miss }' \
        "$TMP/vcf_samples.txt" "$STRAINS" | sort -u > "$TMP/samples.txt"
    if [[ -s "$TMP/missing.txt" ]]; then
        echo "WARNING: $(wc -l < "$TMP/missing.txt") strain(s) not in the VCF, skipped: $(tr '\n' ' ' < "$TMP/missing.txt")" >&2
    fi
    [[ -s "$TMP/samples.txt" ]] || { echo "no strain from ${STRAINS} is in the VCF" >&2; exit 1; }
    printf "strains_requested\t%s\nstrains_in_vcf\t%s\n" "$(grep -c . "$STRAINS")" "$(wc -l < "$TMP/samples.txt")" >> "$STATS"
    SAMPLES="$TMP/samples.txt"
fi

printf "input\t%s\n" "$(bcftools view -H "$IN" | wc -l)" >> "$STATS"

# 1-2. PASS, outside mask
MASK_ARGS=()
if [[ -n "$MASK" ]]; then
    MASK_ARGS=(-T "^${MASK}" --targets-overlap 1)
fi
SAMPLE_ARGS=()
if [[ -n "$SAMPLES" ]]; then
    SAMPLE_ARGS=(-S "$SAMPLES")
fi
bcftools view --threads "$THREADS" -f PASS "${MASK_ARGS[@]}" "${SAMPLE_ARGS[@]}" -Ou "$IN" \
  | bcftools +fill-tags -Ob -o "$TMP/s1.bcf" -- -t 'INFO/DP:1=int(sum(FMT/DP))'
bcftools index "$TMP/s1.bcf"
count "pass_unmasked" "$TMP/s1.bcf"

# 3. site depth ceiling from the median INFO/DP of the remaining sites
MEDIAN_DP=$(bcftools query -f '%INFO/DP\n' "$TMP/s1.bcf" | sort -n | awk '{a[NR]=$1} END{print (NR%2 ? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2)}')
DP_MAX=$(awk -v m="$MEDIAN_DP" -v f="$DP_MAX_FACTOR" 'BEGIN{printf "%d", m*f}')
printf "median_INFO_DP=%s\tDP_MAX=%s\n" "$MEDIAN_DP" "$DP_MAX" > "${PREFIX}.depth_threshold.txt"

# 3-6. depth ceiling, then genotype masking
bcftools view --threads "$THREADS" -e "INFO/DP > ${DP_MAX}" -Ou "$TMP/s1.bcf" \
  | bcftools +setGT -Ou -- -t q -n . \
        -i "FMT/GQ<${MIN_GQ} | (FMT/DP<${MIN_DP} & GT!=\"ref\") | (FMT/DP<${MIN_DP_HOMREF} & GT=\"ref\")" \
  | bcftools +setGT -Ou -- -t q -n . \
        -i "GT=\"het\" & smpl_sum(FMT/AD)>0 & smpl_max(FMT/AD)/smpl_sum(FMT/AD) > (1-${MIN_AB})" \
  | bcftools +setGT -Ou -- -t q -n . \
        -i "(GT=\"hap\" | GT=\"hom\") & smpl_sum(FMT/AD)>0 & smpl_max(FMT/AD)/smpl_sum(FMT/AD) < ${HAP_MIN_AF}" \
  | bcftools view --threads "$THREADS" -a -Ou \
  | bcftools view --threads "$THREADS" -c 1 -Ou \
  | bcftools +fill-tags -Ob -o "$TMP/s2.bcf" -- -t AN,AC,AF,MAF,F_MISSING
bcftools index "$TMP/s2.bcf"
count "depth_ok_gt_masked_polymorphic" "$TMP/s2.bcf"

# 8. site missingness
bcftools view --threads "$THREADS" -i "INFO/F_MISSING <= ${MAX_MISSING}" -Oz -o "${PREFIX}.qc.vcf.gz" "$TMP/s2.bcf"
bcftools index -t "${PREFIX}.qc.vcf.gz"
count "missing_le_${MAX_MISSING}" "${PREFIX}.qc.vcf.gz"

# 9. biallelic SNPs, MAF
bcftools view --threads "$THREADS" -m2 -M2 -v snps -i "ALT!=\"*\" && INFO/MAF >= ${MIN_MAF}" \
    -Oz -o "${PREFIX}.snps.maf.vcf.gz" "${PREFIX}.qc.vcf.gz"
bcftools index -t "${PREFIX}.snps.maf.vcf.gz"
count "biallelic_snps_maf_ge_${MIN_MAF}" "${PREFIX}.snps.maf.vcf.gz"
