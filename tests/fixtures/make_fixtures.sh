#!/usr/bin/env bash
# tests/fixtures/make_fixtures.sh
#
# Generates a tiny, deterministic synthetic reference + haploid/diploid CRAM
# fixture pair for Phase 1 smoke tests (Tasks 5, 7, 9, 11, 13).
#
# Requires (load as environment modules on UCR HPCC, run interactively, not
# via SLURM):
#   module load bwa-mem2/2.3
#   module load samtools/1.19.2   # bundles the `wgsim` read simulator
#   module load gatk/4.6.2.0      # for CreateSequenceDictionary (or picard/2.26.11)
#   /usr/bin/python3.12 for the reference/haplotype generation snippets
#
# NOTE on design: the original brief's approach was to use dwgsim to mutate a
# single template and simulate reads from it. That produces 100%
# homozygous-alt calls at every mutated site (every read at a site comes from
# the one mutated copy) -- indistinguishable from a haploid fixture.
#
# The first fix attempted here was `wgsim -h` (documented as "haplotype
# mode"). Empirically verified via `samtools mpileup` on the resulting CRAM:
# -h does NOT mix reference- and alternate-allele reads in the samtools
# 1.19.2-bundled wgsim build -- every read at every mutated site still
# carried only the ALT allele (0 ref-supporting reads at each site checked).
# So -h was rejected after direct verification, not assumed to work.
#
# Working approach used instead: build two haplotype FASTAs (haplotype A =
# the unmutated reference itself, haplotype B = the reference with ~1% SNPs
# introduced at deterministic seeded positions), simulate an equal number of
# read pairs independently from each haplotype with `wgsim -r 0` (no further
# mutation during simulation), concatenate the two read sets, and map the
# combined pool back to the single reference. This produces a genuine ~50/50
# read depth split at each mutated site -- true heterozygous (0/1) signal
# when force-called at ploidy 2.
set -euo pipefail
cd "$(dirname "$0")"

# 2kb synthetic single-contig reference, deterministic (seeded), plus a
# second haplotype (hapB) carrying ~1% SNPs at deterministic positions.
/usr/bin/python3.12 - <<'PY'
import random

random.seed(42)
seq = "".join(random.choice("ACGT") for _ in range(2000))
with open("ref.fa", "w") as fh:
    fh.write(">synth_contig1\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i + 60] + "\n")

# Haplotype B: independent seed, ~1% SNP rate, always substitutes to a
# different base than the reference at each chosen site.
random.seed(43)
bases = "ACGT"
seqB = list(seq)
for i in range(len(seqB)):
    if random.random() < 0.01:
        ref_b = seqB[i]
        seqB[i] = random.choice([b for b in bases if b != ref_b])
seqB = "".join(seqB)
with open("hapB.fa", "w") as fh:
    fh.write(">synth_contig1\n")
    for i in range(0, len(seqB), 60):
        fh.write(seqB[i:i + 60] + "\n")
PY
samtools faidx ref.fa
rm -f ref.dict
gatk CreateSequenceDictionary -R ref.fa -O ref.dict 2>/dev/null || \
    picard CreateSequenceDictionary R=ref.fa O=ref.dict
bwa-mem2 index ref.fa

# Haploid strain: simulate reads from the reference itself (-r 0: no
# injected mutations at all), so heterozygosity should be ~0.
wgsim -N 2000 -1 100 -2 100 -r 0 -e 0.001 -S 42 \
    ref.fa haploid_strain.read1.fq haploid_strain.read2.fq
bwa-mem2 mem ref.fa haploid_strain.read1.fq haploid_strain.read2.fq \
    | samtools sort -O cram --reference ref.fa -o haploid_strain.cram -
samtools index haploid_strain.cram

# Diploid strain: simulate half the read pairs from haplotype A (the
# reference, unmutated) and half from haplotype B (~1% SNPs), each with
# -r 0 so wgsim injects no further mutation during simulation, then pool
# both read sets before mapping to the single reference. This gives a
# genuine ~50/50 read mixture at each hapB SNP site -- true heterozygous
# depth, not single-haplotype homozygous-alt.
bwa-mem2 index hapB.fa
wgsim -N 1000 -1 100 -2 100 -r 0 -e 0.001 -S 42 \
    ref.fa diploid_hapA.read1.fq diploid_hapA.read2.fq
wgsim -N 1000 -1 100 -2 100 -r 0 -e 0.001 -S 43 \
    hapB.fa diploid_hapB.read1.fq diploid_hapB.read2.fq
cat diploid_hapA.read1.fq diploid_hapB.read1.fq > diploid_strain.read1.fq
cat diploid_hapA.read2.fq diploid_hapB.read2.fq > diploid_strain.read2.fq
bwa-mem2 mem ref.fa diploid_strain.read1.fq diploid_strain.read2.fq \
    | samtools sort -O cram --reference ref.fa -o diploid_strain.cram -
samtools index diploid_strain.cram

rm -f *.read1.fq *.read2.fq *.mutations.txt hapB.fa hapB.fa.* \
    ref.fa.amb ref.fa.ann ref.fa.bwt.2bit.64 ref.fa.pac ref.fa.0123
echo "Fixtures written to $(pwd)"
