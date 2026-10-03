#!/usr/bin/env python3
"""Estimate ploidy (haploid vs diploid) from heterozygous sites per callable Mb.

Calls variants in forced diploid mode via `bcftools mpileup | bcftools call`,
counts biallelic SNPs called heterozygous with a balanced allele ratio, and
divides by the number of callable bases (read depth >= MIN_DP at the same
mapping and base quality floors, outside the optional repeat mask).

A true haploid genome shows near-zero balanced heterozygosity (mapping and
base-calling noise only); a diploid or hybrid genome shows many balanced het
sites per Mb.

Why per callable Mb (2026-10-03): the earlier metric divided by the number of
variant sites (het + hom-alt). `bcftools call -v` writes only variant sites,
so that denominator is mostly the hom-alt count, i.e. the distance from the
reference. For a strain close to the reference it is tiny, and a few noise
hets in repeats gave het fractions of 0.1-0.2. In the DH4148 population, 114
of 115 near-reference haploids were called diploid that way (see
Rhodotorula_mucilaginosa_DH4148_ref/docs/ploidy_assessment_2026-10-03.md).
The old ratio is still reported as `het_fraction` for comparison only.
"""

import argparse
import subprocess
import sys


def count_het_sites(vcf_lines, min_allele_balance=0.35, max_allele_balance=0.65):
    """Count biallelic-SNP records and the balanced heterozygous ones.

    Returns (n_sites, n_het). n_sites counts every biallelic SNP with a called
    diploid GT. A GT=0/1 call counts as heterozygous only if its allele
    balance (alt fraction of AD) is within [min_allele_balance,
    max_allele_balance]. Genuine diploid het sites carry ~50% allele balance
    (confirmed against strain DBVPG_3857, metadata-diploid); a uniform skewed
    signal (~0.25-0.3, strain EXF_12768, metadata-haploid) must not count. If
    AD is absent from FORMAT, the GT call is trusted as-is.
    """
    n_sites = 0
    n_het = 0
    for line in vcf_lines:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        ref, alt = fields[3], fields[4]
        if len(ref) != 1 or len(alt) != 1:
            continue  # biallelic SNPs only - skip indels and multiallelic sites
        format_keys = fields[8].split(":")
        sample_values = fields[9].split(":")
        if "GT" not in format_keys:
            continue
        gt = sample_values[format_keys.index("GT")]
        alleles = gt.replace("|", "/").split("/")
        if len(alleles) != 2 or "." in alleles:
            continue
        n_sites += 1
        if alleles[0] == alleles[1]:
            continue
        if "AD" in format_keys:
            ad = sample_values[format_keys.index("AD")].split(",")
            if len(ad) == 2 and "." not in ad:
                ref_reads, alt_reads = int(ad[0]), int(ad[1])
                total = ref_reads + alt_reads
                if total > 0:
                    alt_fraction = alt_reads / total
                    if not (min_allele_balance <= alt_fraction <= max_allele_balance):
                        continue
        n_het += 1
    return n_sites, n_het


def het_fraction_from_vcf(vcf_lines, min_allele_balance=0.35, max_allele_balance=0.65):
    """Old metric: balanced hets / biallelic variant sites. None if no sites.

    Kept for comparison only. Do not use it to call ploidy: its denominator
    shrinks with distance from the reference (see the module docstring).
    """
    n_sites, n_het = count_het_sites(vcf_lines, min_allele_balance, max_allele_balance)
    if n_sites == 0:
        return None
    return n_het / n_sites


def load_mask(bed_path, contig_lengths):
    """BED -> {contig: bytearray}, 1 at masked 0-based positions."""
    mask = {}
    if not bed_path:
        return mask
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.rstrip("\n").split("\t")[:3]
            if c not in contig_lengths:
                continue
            if c not in mask:
                mask[c] = bytearray(contig_lengths[c])
            s, e = max(0, int(s)), min(int(e), contig_lengths[c])
            if e > s:
                mask[c][s:e] = b"\x01" * (e - s)
    return mask


def count_callable(depth_lines, min_dp, mask):
    """Count positions with depth >= min_dp outside the mask.

    depth_lines: `samtools depth` output (contig, 1-based pos, depth).
    """
    n = 0
    for line in depth_lines:
        c, p, d = line.rstrip("\n").split("\t")[:3]
        if int(d) < min_dp:
            continue
        m = mask.get(c)
        if m is not None and m[int(p) - 1]:
            continue
        n += 1
    return n


def het_per_mb(n_het, callable_bp):
    if not callable_bp:
        return None
    return n_het / callable_bp * 1e6


def call_ploidy(value, threshold, callable_bp=None, min_callable_bp=0):
    """value: balanced hets per callable Mb.

    'unknown' if value is None or callable_bp is below min_callable_bp: with
    little callable sequence the rate is unstable (DH4148 strains with
    0.9-2.8 Mb callable gave 655-10,336 hets/Mb whatever their ploidy).
    """
    if value is None:
        return "unknown"
    if callable_bp is not None and callable_bp < min_callable_bp:
        return "unknown"
    return "diploid" if value > threshold else "haploid"


def read_fai(fai):
    lengths = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            lengths[f[0]] = int(f[1])
    return lengths


def run_variant_calls(cram, reference, min_mapq, min_baseq, min_qual, min_dp, mask_bed, region=None):
    mpileup = ["bcftools", "mpileup", "-Ou", "-f", reference, "-a", "FORMAT/AD,FORMAT/DP",
               "-q", str(min_mapq), "-Q", str(min_baseq)]
    if region:
        mpileup += ["-r", region]
    mpileup.append(cram)
    call = ["bcftools", "call", "-mv", "--ploidy", "2", "-Ou"]
    view = ["bcftools", "view", "-i", f"QUAL>={min_qual} && FMT/DP>={min_dp}"]
    if mask_bed:
        view += ["-T", f"^{mask_bed}", "--targets-overlap", "1"]
    p1 = subprocess.Popen(mpileup, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(call, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    result = subprocess.run(view, stdin=p2.stdout, capture_output=True, text=True, check=True)
    p2.stdout.close()
    for p in (p1, p2):
        if p.wait() != 0:
            raise subprocess.CalledProcessError(p.returncode, p.args)
    return result.stdout.splitlines()


def run_callable(cram, reference, min_mapq, min_baseq, min_dp, mask, region=None):
    cmd = ["samtools", "depth", "--reference", reference, "-q", str(min_baseq), "-Q", str(min_mapq)]
    if region:
        cmd += ["-r", region]
    cmd.append(cram)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    n = count_callable(p.stdout, min_dp, mask)
    if p.wait() != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd)
    return n


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cram", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--reference-fai", default=None, help="default: <reference>.fai")
    parser.add_argument("--strain", required=True)
    parser.add_argument("--region", default=None)
    parser.add_argument("--mask", default=None, help="repeat / low-complexity BED to exclude")
    parser.add_argument("--threshold", type=float, default=1000.0,
                        help="diploid if balanced hets per callable Mb is above this")
    parser.add_argument("--min-callable-bp", type=int, default=5_000_000,
                        help="'unknown' if fewer callable bases than this")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--min-qual", type=float, default=30)
    parser.add_argument("--min-dp", type=int, default=10)
    parser.add_argument("--min-allele-balance", type=float, default=0.35)
    parser.add_argument("--max-allele-balance", type=float, default=0.65)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    lengths = read_fai(args.reference_fai or args.reference + ".fai")
    mask = load_mask(args.mask, lengths)
    vcf_lines = run_variant_calls(args.cram, args.reference, args.min_mapq, args.min_baseq,
                                  args.min_qual, args.min_dp, args.mask, args.region)
    n_sites, n_het = count_het_sites(vcf_lines, args.min_allele_balance, args.max_allele_balance)
    callable_bp = run_callable(args.cram, args.reference, args.min_mapq, args.min_baseq,
                               args.min_dp, mask, args.region)
    value = het_per_mb(n_het, callable_bp)
    ploidy = call_ploidy(value, args.threshold, callable_bp, args.min_callable_bp)

    with open(args.out, "w") as fh:
        fh.write("strain,inferred_ploidy,het_fraction,method,het_per_mb,n_het,n_snp_sites,callable_bp\n")
        hf = f"{n_het / n_sites:.6f}" if n_sites else "NA"
        hpm = f"{value:.3f}" if value is not None else "NA"
        fh.write(f"{args.strain},{ploidy},{hf},custom_het_per_mb,{hpm},{n_het},{n_sites},{callable_bp}\n")


if __name__ == "__main__":
    sys.exit(main())
