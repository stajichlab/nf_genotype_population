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
