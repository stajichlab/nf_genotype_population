import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from custom_het_ploidy import (het_fraction_from_vcf, call_ploidy, count_het_sites,
                               count_callable, het_per_mb, load_mask)

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


def test_call_ploidy_thresholds_per_mb():
    assert call_ploidy(0.0, 1000) == "haploid"
    assert call_ploidy(999.9, 1000) == "haploid"
    assert call_ploidy(1000.1, 1000) == "diploid"
    assert call_ploidy(None, 1000) == "unknown"


def test_call_ploidy_unknown_below_min_callable():
    assert call_ploidy(5000.0, 1000, callable_bp=1_500_000, min_callable_bp=5_000_000) == "unknown"
    assert call_ploidy(5000.0, 1000, callable_bp=19_000_000, min_callable_bp=5_000_000) == "diploid"


def test_het_per_mb():
    assert het_per_mb(31, 1_000_000) == 31.0
    assert het_per_mb(5, 0) is None


def test_count_het_sites_returns_sites_and_balanced_hets():
    assert count_het_sites(DIPLOID_LIKE_VCF) == (3, 2)
    assert count_het_sites(HAPLOID_LIKE_VCF) == (3, 0)


# Regression for the 2026-10-03 defect: a strain close to the reference has
# few hom-alt sites, so the old ratio (balanced hets / variant sites) is large
# from a handful of noise hets. On contig CM179498.1, TFCN_25-332D-2 had 37
# hom-alt and 3 balanced hets over ~1 Mb. Per callable Mb the same strain is
# clearly haploid.
NEAR_REFERENCE_VCF = (
    ["chr1\t%d\t.\tA\tG\t50\tPASS\t.\tGT:AD\t1/1:0,20\n" % (100 + i) for i in range(4)]
    + ["chr1\t%d\t.\tC\tT\t50\tPASS\t.\tGT:AD\t0/1:10,10\n" % (900 + i) for i in range(1)]
)


def test_near_reference_strain_old_ratio_vs_per_mb():
    n_sites, n_het = count_het_sites(NEAR_REFERENCE_VCF)
    assert (n_sites, n_het) == (5, 1)
    assert n_het / n_sites > 0.01                       # old metric: "diploid"
    assert call_ploidy(het_per_mb(n_het, 1_000_000), 1000) == "haploid"


def test_count_callable_applies_depth_floor_and_mask(tmp_path):
    bed = tmp_path / "m.bed"
    bed.write_text("chr1\t1\t3\n")                     # masks 1-based positions 2 and 3
    mask = load_mask(str(bed), {"chr1": 10})
    depth = ["chr1\t%d\t%d\n" % (p, d) for p, d in [(1, 12), (2, 30), (3, 30), (4, 9), (5, 10)]]
    # pos1 kept, pos2-3 masked, pos4 below floor, pos5 kept
    assert count_callable(depth, 10, mask) == 2
    assert count_callable(depth, 10, {}) == 4


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
