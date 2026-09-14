import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from custom_het_ploidy import het_fraction_from_vcf, call_ploidy

HAPLOID_LIKE_VCF = [
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\n",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t1/1\n",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t0/0\n",
]

DIPLOID_LIKE_VCF = [
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n",
    "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t1/1\n",
]


def test_het_fraction_haploid_like_is_zero():
    assert het_fraction_from_vcf(HAPLOID_LIKE_VCF) == 0.0


def test_het_fraction_diploid_like_is_high():
    frac = het_fraction_from_vcf(DIPLOID_LIKE_VCF)
    assert abs(frac - (2 / 3)) < 1e-9


def test_het_fraction_none_when_no_sites():
    assert het_fraction_from_vcf(["#comment\n"]) is None


def test_het_fraction_skips_multiallelic_and_indels():
    lines = [
        "chr1\t100\t.\tA\tG,T\t50\tPASS\t.\tGT\t1/2\n",   # multiallelic, skip
        "chr1\t200\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1\n",     # indel, skip
        "chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n",       # counts
    ]
    assert het_fraction_from_vcf(lines) == 1.0


def test_call_ploidy_thresholds():
    assert call_ploidy(0.0) == "haploid"
    assert call_ploidy(0.005) == "haploid"
    assert call_ploidy(0.02) == "diploid"
    assert call_ploidy(None) == "unknown"


# Real-strain pattern found debugging EXF_12768 (2026-09-13): a het GT call is
# not proof of real heterozygosity. Genuine diploid het sites (confirmed
# against strain DBVPG_3857, metadata-diploid) carry ~50% allele balance;
# EXF_12768 (metadata-haploid) called "het" at 99.8% of biallelic sites
# genome-wide but with allele balance skewed to ~0.17-0.35 (median 0.26) at
# every one of them - a uniform low-level minor-allele signal, not real
# heterozygosity. A GT=0/1 call with skewed AD must not count as heterozygous.

REAL_DIPLOID_HET_VCF = [
    # AD ref=11,alt=9 -> alt fraction 0.45, within the real-het band
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:11,9\n",
]

SKEWED_CONTAMINATION_LIKE_VCF = [
    # AD ref=37,alt=13 -> alt fraction 0.26, matches the EXF_12768 pattern
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:37,13\n",
]


def test_het_fraction_counts_balanced_het_site():
    assert het_fraction_from_vcf(REAL_DIPLOID_HET_VCF) == 1.0


def test_het_fraction_excludes_skewed_allele_balance_site():
    assert het_fraction_from_vcf(SKEWED_CONTAMINATION_LIKE_VCF) == 0.0


def test_het_fraction_falls_back_to_gt_only_when_ad_absent():
    # No AD in FORMAT - can't check allele balance, so a het GT call still
    # counts (existing behavior preserved for callers that don't emit AD).
    lines = ["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"]
    assert het_fraction_from_vcf(lines) == 1.0
