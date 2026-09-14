#!/usr/bin/env python3
"""Estimate ploidy signal from a CRAM using nQuire (github.com/clwgg/nQuire).

Pipeline: CRAM -> BAM (nQuire's own bundled htslib build explicitly disables
full CRAM support - see its Makefile's `htsconf` target comment - so CRAM
must be converted first) -> `nQuire create` -> `nQuire denoise` ->
`nQuire lrdmodel`.

IMPORTANT METHODOLOGICAL CAVEAT, not yet resolved - do not treat this
script's ploidy call as authoritative without validation:

nQuire's `lrdmodel` only tests three FIXED models: diploid, triploid,
tetraploid (it reports a delta log-likelihood per model - lower delta means
that model fits better; see its own README). It was designed for detecting
polyploidy (e.g. in plants), and has **no built-in haploid model at all**.
There is no published, verified threshold (in this project) for converting
"diploid fits best, with delta X" into a confident haploid-vs-diploid call
for a fungal organism, which is what this project actually needs.

REVISION HISTORY, both changes driven by real data, not guesses:

1. An earlier version of this function used a fixed absolute cutoff
   (MAX_CONFIDENT_DIPLOID_DELTA = 0.05) on the diploid delta to decide
   "confident enough to call diploid". Running this against 6 REAL DH4148
   CRAMs (not the synthetic test fixtures) showed the delta magnitude is
   NOT a normalized confidence score - it scales with the number of
   informative sites, and came back in the hundreds to tens-of-thousands
   for every real strain (e.g. 453.6, 6025.1, 12768.8). An absolute
   threshold of 0.05 therefore rejected every single real strain as
   "unknown" - it was unusable on real data. That threshold has been
   removed.

2. In its place, this function now uses ONLY `best_fit_model` (which of
   diploid/triploid/tetraploid nQuire's lrdmodel scored lowest,
   i.e. best) - not any delta magnitude. In a 6-strain real-data pilot
   (3 metadata-haploid, 3 metadata-diploid DH4148 strains), best_fit_model
   perfectly separated the two groups: all 3 diploid-by-metadata strains
   had best_fit_model == "diploid"; all 3 haploid-by-metadata strains had
   best_fit_model == "tetraploid" (nQuire, lacking a haploid model,
   apparently defaults to its most extreme available model for haploid
   input). n=6 is FAR too small to validate this as a general rule - it is
   reported here as an observed pattern from real data, not a proven one.

Given (2), this script maps best_fit_model == "diploid" to inferred_ploidy
"diploid", and best_fit_model in {"triploid", "tetraploid"} to
"non_diploid" - deliberately NOT "haploid": nQuire never tests a haploid
model, so calling a tetraploid-best-fit result "haploid" would be an
unsupported leap dressed up as a direct measurement. "non_diploid" says
exactly what was actually observed (this sample did not fit the diploid
model) and leaves the haploid interpretation to whoever reviews the result
alongside other evidence (e.g. custom_het_script's independent inference,
metadata). Degenerate/uninformative lrdmodel output (see
call_ploidy_from_deltas below) still maps to "unknown".
"""
import argparse
import os
import shutil
import subprocess
import tempfile


def convert_cram_to_bam(cram, reference, out_bam):
    with open(out_bam, "wb") as fh:
        subprocess.run(
            ["samtools", "view", "-b", "-T", reference, cram],
            stdout=fh,
            check=True,
        )
    subprocess.run(["samtools", "index", out_bam], check=True)


def run_nquire_pipeline(nquire_bin, bam, workdir):
    base = os.path.join(workdir, "base")
    denoised = os.path.join(workdir, "base_denoised")
    subprocess.run([nquire_bin, "create", "-b", bam, "-o", base], check=True)
    subprocess.run(
        [nquire_bin, "denoise", base + ".bin", "-o", denoised], check=True
    )
    result = subprocess.run(
        [nquire_bin, "lrdmodel", denoised + ".bin"],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def parse_lrdmodel_output(lrdmodel_stdout):
    """Parse nQuire lrdmodel's tab-separated stdout.

    Columns (per nQuire's own documentation): filename, free-model
    log-likelihood, diploid/triploid/tetraploid fixed-model log-likelihoods,
    then diploid/triploid/tetraploid delta log-likelihoods (8 columns total).
    Returns a dict with the three delta values, or None if no data line
    was found (e.g. nQuire produced no output for an empty/low-coverage
    input).
    """
    for line in lrdmodel_stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("file"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        return {
            "diploid_delta": float(fields[5]),
            "triploid_delta": float(fields[6]),
            "tetraploid_delta": float(fields[7]),
        }
    return None


def call_ploidy_from_deltas(deltas):
    """Return (inferred_ploidy, best_fit_model) from nQuire's raw deltas.

    See the module docstring for the two real-data-driven revisions this
    function has been through. Current rule: best_fit_model == "diploid" ->
    "diploid"; best_fit_model in {"triploid", "tetraploid"} -> "non_diploid"
    (deliberately not "haploid" - nQuire never tests a haploid model, so that
    would be an unvalidated leap). Based on only a 6-strain real-data pilot -
    treat as an observed pattern, not a proven rule.

    CONFIRMED FAILURE MODE (found by running the real haploid test fixture,
    not hypothesised): a sample with too few informative sites for `create`/
    `denoise` to work with - which is exactly what a true haploid sample
    looks like to nQuire, since it has no het-site signal at all and no
    haploid model - makes `lrdmodel` report all three fixed-model deltas as
    identical (observed: 0.000000/0.000000/0000000 for the haploid fixture).
    That is degenerate output (lrdmodel could not distinguish any fixed model
    from the free model, i.e. zero information), not a valid best fit, so it
    is rejected as "unknown" rather than mapped to any model. NaN values
    (also observed, on the tiny diploid fixture) are handled the same way
    for the same reason: no real information to call from.
    """
    if deltas is None:
        return "unknown", None
    values = (deltas["diploid_delta"], deltas["triploid_delta"], deltas["tetraploid_delta"])
    if any(v != v for v in values):  # NaN != NaN
        return "unknown", None
    if len(set(values)) == 1:
        return "unknown", None
    best_fit_model = min(
        ("diploid", "triploid", "tetraploid"),
        key=lambda m: deltas[f"{m}_delta"],
    )
    inferred_ploidy = "diploid" if best_fit_model == "diploid" else "non_diploid"
    return inferred_ploidy, best_fit_model


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cram", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--strain", required=True)
    parser.add_argument("--nquire-bin", default="nQuire")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    with tempfile.TemporaryDirectory() as workdir:
        bam = os.path.join(workdir, f"{args.strain}.bam")
        convert_cram_to_bam(args.cram, args.reference, bam)
        stdout = run_nquire_pipeline(args.nquire_bin, bam, workdir)

    deltas = parse_lrdmodel_output(stdout)
    inferred_ploidy, best_fit_model = call_ploidy_from_deltas(deltas)

    with open(args.out, "w") as fh:
        fh.write(
            "strain,inferred_ploidy,best_fit_model,diploid_delta,"
            "triploid_delta,tetraploid_delta,method\n"
        )
        if deltas is None:
            fh.write(f"{args.strain},unknown,NA,NA,NA,NA,nquire\n")
        else:
            fh.write(
                f"{args.strain},{inferred_ploidy},{best_fit_model or 'NA'},"
                f"{deltas['diploid_delta']:.6f},"
                f"{deltas['triploid_delta']:.6f},"
                f"{deltas['tetraploid_delta']:.6f},nquire\n"
            )


if __name__ == "__main__":
    main()
