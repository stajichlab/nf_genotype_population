#!/usr/bin/env python3
"""Estimate ploidy (haploid vs diploid) from CRAM heterozygous-site rate.

Calls variants at forced diploid genotype-likelihood mode via
`bcftools mpileup | bcftools call`, then computes the fraction of covered,
biallelic SNP sites called heterozygous. A true haploid genome shows
near-zero heterozygosity (mapping/base-calling noise only); a diploid or
hybrid genome shows a substantial het fraction. This is the fast,
transparent Phase 1 ploidy check - nQuire/nQuack are added in Phase 2 as
additional cross-validated methods, not replacements for this one.
"""
import argparse
import subprocess


def het_fraction_from_vcf(vcf_lines):
    """Fraction of biallelic-SNP records called heterozygous.

    vcf_lines: iterable of VCF data/header lines (header lines, i.e. those
    starting with '#', are ignored). Returns None if no qualifying sites
    were found (e.g. no coverage at all).
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
        if alleles[0] != alleles[1]:
            n_het += 1
    if n_sites == 0:
        return None
    return n_het / n_sites


def call_ploidy(het_fraction, threshold=0.01):
    if het_fraction is None:
        return "unknown"
    return "diploid" if het_fraction > threshold else "haploid"


def run_bcftools_forced_diploid_call(cram, reference, region=None):
    mpileup_cmd = ["bcftools", "mpileup", "-Ou", "-f", reference]
    if region:
        mpileup_cmd += ["-r", region]
    mpileup_cmd.append(cram)
    call_cmd = ["bcftools", "call", "-mv", "--ploidy", "2"]
    mpileup = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
    result = subprocess.run(call_cmd, stdin=mpileup.stdout, capture_output=True, text=True, check=True)
    mpileup.stdout.close()
    mpileup.wait()
    return result.stdout.splitlines()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cram", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--strain", required=True)
    parser.add_argument("--region", default=None)
    parser.add_argument("--threshold", type=float, default=0.01)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    vcf_lines = run_bcftools_forced_diploid_call(args.cram, args.reference, args.region)
    het_fraction = het_fraction_from_vcf(vcf_lines)
    ploidy = call_ploidy(het_fraction, args.threshold)

    with open(args.out, "w") as fh:
        fh.write("strain,inferred_ploidy,het_fraction,method\n")
        hf = f"{het_fraction:.6f}" if het_fraction is not None else "NA"
        fh.write(f"{args.strain},{ploidy},{hf},custom_het_script\n")


if __name__ == "__main__":
    main()
