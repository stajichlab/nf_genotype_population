# nf_genotype_population

Ploidy-aware per-strain GATK4 calling, joint genotyping, and hard-filtering
for Rhodotorula population genomics on UCR HPCC. Consumes CRAM output from
an nf-core/sarek run (Sarek itself only used through alignment/dedup - see
this project's design spec).

Design spec: see `Rhodotorula_mucilaginosa_DH4148_ref/docs/superpowers/specs/2026-09-12-nf-genotype-population-design.md`.

## Usage

This pipeline runs as a two-phase workflow:

### Phase 1: Ploidy Inference Review

Run the PLOIDY_ONLY entry point to infer ploidy per strain and generate a human-reviewable report:

```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf -entry PLOIDY_ONLY -profile hpcc \
    --cram_dir /path/to/flat_cram_dir \
    --reference /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna \
    --reference_fai /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna.fai \
    --metadata /path/to/metadata.txt
```

**Important:** This cluster's Nextflow (26.04.6) requires the environment variable `NXF_SYNTAX_PARSER=v1` set before the command. Without it, Nextflow errors with "The `-entry` option is not supported with the strict parser".

This produces `results/ploidy_review.csv`, with columns `strain`, `inferred_ploidy` (haploid/diploid), `metadata_ploidy`, and `status` (AGREE/DISAGREE/NO_METADATA/UNKNOWN_INFERENCE).

Add `--mask_bed` here too: the default method (`custom_het_script`, `bin/custom_het_ploidy.py`) excludes the masked repeats.

How `custom_het_script` decides:
1. It calls variants with `bcftools mpileup -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 | bcftools call -mv --ploidy 2`.
2. It keeps sites with QUAL >= 30 and DP >= 10 that are outside the mask.
3. It counts the biallelic SNPs called heterozygous with allele balance 0.35-0.65.
4. It divides that count by the callable bases: `samtools depth` at the same quality floors, depth >= 10, outside the mask.
5. A strain is diploid if it has more than `--ploidy_het_per_mb_threshold` (default 1000) balanced het SNPs per callable Mb. It is `unknown` if its callable bases are less than `--ploidy_min_callable_fraction` (default 0.25) of the reference length.

These defaults come from a run on all 319 DH4148 strains (2026-10-03, genome 20.4 Mb):
- Among strains with at least 10 Mb callable, haploids scored at most 453 hets/Mb and diploids at least 1,785.
- One exception: EXF_12768, a haploid with mixed reads, scored 2,295.
- Strains with 0.9-2.8 Mb callable (4-14% of the genome) gave unstable values (655-10,336).
- A strain with 5.2 Mb callable (25%) gave a normal haploid value.

`results/ploidy/ploidy_inference_all.csv` reports `het_per_mb`, `n_het`, `n_snp_sites` and `callable_bp`. It also reports the old ratio `het_fraction` (balanced hets / variant sites) for comparison only.

The old ratio must not be used to call ploidy. Its denominator is mostly the hom-alt count, so it shrinks for strains close to the reference. In the DH4148 population it called 114 of 115 near-reference haploids diploid.

### Phase 2: Manual Review and Full Genotyping

#### Building `ploidy_overrides.csv`

**`ploidy_overrides.csv` must contain ONE ROW PER STRAIN that has a CRAM in `--cram_dir` — not just the rows you want to change.**

This file is the *only* source of calling ploidy. The inferred ploidy is never
consulted for calling under any code path, and a strain with a CRAM but no row
here **aborts the run** (the join in `workflows/genotype_population.nf` uses
`failOnMismatch: true`, deliberately, so a strain cannot silently vanish from
the callset). The reverse also aborts: a row here with no matching CRAM.

1. **Review** `results/ploidy_review.csv`, paying attention to the `status`
   column. `DISAGREE` means the inference and `metadata.txt` disagree;
   `UNKNOWN_INFERENCE` means the inference could not decide.
2. **Derive** the full-coverage file from the review CSV:

   ```bash
   { echo "strain,ploidy"; tail -n +2 results/ploidy_review.csv | cut -d, -f1,2; } > ploidy_overrides.csv
   ```

   That takes the `strain` and `inferred_ploidy` columns for **every** strain
   and writes them under the header `strain,ploidy` that the pipeline expects.
3. **Hand-edit** the result before using it. Set the correct value for every
   `DISAGREE` and `UNKNOWN_INFERENCE` row — an `UNKNOWN_INFERENCE` row carries
   the literal value `unknown`, which is not a valid ploidy and will stop the
   run with an "Unknown ploidy label" error. Only `haploid` and `diploid` are
   accepted.
4. **Run** the full pipeline:

```bash
nextflow run main.nf -profile hpcc \
    --cram_dir /path/to/flat_cram_dir \
    --reference /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna \
    --reference_fai /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna.fai \
    --reference_dict /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.dict \
    --metadata /path/to/metadata.txt \
    --ploidy_overrides /path/to/ploidy_overrides.csv \
    --population_sets assets/population_sets.yaml \
    --snpeff_db_dir /path/to/snpeff_db \
    --snpeff_genome_name RmucDH4148 \
    --mask_bed /path/to/mask.bed
```

#### Variant QC

`VARIANT_QC_FILTER` (`bin/filter_population_vcf.sh`) runs between GATK hard
filtering and SnpEff. GATK hard filters only flag sites in FILTER; this step
removes them and filters genotypes. In order:

1. keep `FILTER=PASS`; drop sites that overlap `--mask_bed` (repeats, low complexity)
2. drop sites with summed `FORMAT/DP` > `qc_dp_max_factor` x the median site depth
3. set genotypes to missing where `GQ < qc_min_gq` or `DP < qc_min_dp`
4. set het genotypes to missing where one allele has more than `1 - qc_min_ab` of the allelic reads
5. set haploid / homozygous genotypes to missing where no allele has at least `qc_hap_min_af` of the reads
6. drop ALT alleles no genotype still carries and sites left monomorphic; recompute AC/AN/AF/MAF/F_MISSING
7. drop sites with `F_MISSING > qc_max_missing` → `<pop>.qc`
8. biallelic SNPs with `MAF >= qc_min_maf` → `<pop>.snps.maf`

The record count after each step is written to `<pop>.filter_stats.tsv`.

#### Strain tree

SNP_ALIGNMENT and IQTREE port `06_make_SNP_tree.sh` and `07_iqtree.sh` from the
bash PopGenomics pipeline. They run on `<pop>.snps.maf`:

1. keep sites with a missing-genotype fraction `<= tree_max_missing` (default 0, as in the bash pipeline)
2. write one FASTA row per strain plus one row for the reference; haploid calls give one base, diploid het calls give an IUPAC code, missing calls give `N`
3. IQ-TREE 3 (`-st DNA`) with `tree_model` (default `GTR+ASC`), 1000 ultrafast bootstraps and 1000 SH-aLRT replicates. If `+ASC` stops on invariant columns, IQ-TREE writes `<pop>.snps.maf.varsites.phy` and the task reruns on it. IQ-TREE treats an IUPAC code as compatible with its bases, so a column that varies only by diploid het codes counts as invariant and is dropped.

`--skip_tree` turns the step off.

#### Outputs

- `results/<pop>.qc.annotated.vcf.gz` (+ `.tbi`) — QC-filtered, annotated callset: all variant types, no MAF floor.
- `results/<pop>.snps.maf.annotated.vcf.gz` (+ `.tbi`) — biallelic SNPs with the MAF floor, annotated (population structure, PCA, GWAS-type analyses).
- `results/<pop>.filter_stats.tsv`, `results/<pop>.depth_threshold.txt` — records kept at each QC step, and the site-depth ceiling used.
- `results/strain_tree/<pop>.snps.maf.mfa.gz`, `<pop>.snps.maf.alignment_stats.tsv` — SNP alignment and its site count.
- `results/strain_tree/<pop>.snps.maf.treefile`, `.iqtree`, `.log` — IQ-TREE tree (Newick, UFBoot/SH-aLRT support), report, and log.
- `results/gvcfs/<strain>.g.vcf.gz` (+ `.tbi`) — **per-strain GVCFs, retained.** Design spec §5 keeps these so new strains can be added and the population re-genotyped without recalling every existing strain from CRAM. They are copied out of `work/`, so they survive a scratch cleanup.
- `results/ploidy_review.csv` and `results/ploidy_inference_all.csv` — the ploidy QC record. Ploidy inference runs on every strain in the full workflow too, so each production run leaves a diagnostic trail. Compare `ploidy_review.csv` against the `ploidy_overrides.csv` you supplied to spot ploidy drift.

#### Populations

The built-in `all` population is derived at runtime from **every strain that
produced a GVCF in the run**, not from `assets/population_sets.yaml`. Adding a
strain to `--cram_dir` and `ploidy_overrides.csv` is therefore sufficient to get
it into the `all` callset; no file needs to be kept in sync. The YAML is
reserved for named sub-population slices (Phase 2); any `all` key in it is
ignored.

Named groups use the bash PopGenomics pipeline format. A flat map without the
`Populations:` key also works:

```yaml
Populations:
  rmuc_core:
    - DBVPG_3044
    - EXF_1510
  hybrid_diploids:
    - DBVPG_3045
```

Each group gets its own QC filter, SnpEff annotation, alignment and tree, named
`<pop>.*`, or `<output_prefix>.<pop>.*` when `--output_prefix` is set (the bash
pipeline's `$PREFIX.$POPNAME` naming). `--population_mode` sets how a group is
built:

- `subset` (default): VARIANT_QC_FILTER cuts the group from the `all` callset.
  It recomputes depth, genotype masks, AC/AN/AF, missingness and MAF over the
  group. GATK site annotations (QD, FS, SOR, MQ) and the hard-filter flags stay
  as computed over all strains. This needs no extra joint genotyping.
- `regenotype`: JOINT_GENOTYPING runs GenomicsDBImport and GenotypeGVCFs again
  for each group from its GVCFs. All site annotations then come from the
  group only, but each group costs a full joint-genotyping run.

Strain names in the YAML match VCF sample `STRAIN` or `STRAIN_STRAIN` (the
form Sarek CRAMs give). A strain with no match is skipped with a warning, and
`<pop>.filter_stats.tsv` records `strains_requested` and `strains_in_vcf`.

### Parameters

- **cram_dir** (required): A FLAT directory containing `<strain>.cram`/`<strain>.cram.crai` pairs directly (no subdirectories — `buildCramCh()` in `main.nf` globs `${cram_dir}/*.cram{,.crai}`, one level only). **Known Phase 1 limitation**: Sarek's own output is nested per-sample (`results/preprocessing/markduplicates/<strain>/<strain>.md.cram`), so it cannot be pointed at directly yet. Flatten it into a single directory first, e.g.:
  ```bash
  mkdir -p flat_cram_dir
  find /path/to/sarek/results/preprocessing/markduplicates -name '*.md.cram*' \
      -exec sh -c 'ln -s "$1" "flat_cram_dir/$(basename "$1" | sed "s/\.md//")"' _ {} \;
  ```
  (renames `<strain>.md.cram` → `<strain>.cram` and `<strain>.md.cram.crai` → `<strain>.cram.crai` via symlinks, matching the strain-name key `buildCramCh()` derives from the filename). Supporting Sarek's nested layout directly (e.g. a recursive glob) is a good candidate for a small follow-up fix, not yet done.
- **reference** (required for both phases): Reference genome FASTA file.
- **reference_fai** (required for both phases): Reference genome FAI index file.
- **reference_dict** (required for Phase 2 only): Reference genome GATK sequence dictionary file (`.dict`).
- **metadata** (required for both phases): Metadata file mapping strain names to attributes (e.g., `metadata.txt`).
- **ploidy_overrides** (required for Phase 2 only): CSV with columns `strain,ploidy`, **one row for every strain that has a CRAM**. Values must be `haploid` or `diploid`. See "Building `ploidy_overrides.csv`" above — a partial file aborts the run.
- **population_sets** (default: `assets/population_sets.yaml`): YAML defining named sub-population groups (Phase 2). The built-in `all` group is **not** read from this file; it is derived at runtime from every strain that produced a GVCF, and an `all` key in the YAML is ignored.
- **snpeff_db_dir** (required for Phase 2 only): Directory containing SnpEff database files.
- **snpeff_genome_name** (required for Phase 2 only): SnpEff genome database name (e.g., `RmucDH4148`).
- **mask_bed** (optional, strongly recommended): BED of repeat and low-complexity intervals to exclude (for example RepeatMasker + dustmasker + TRF, merged). If unset the run warns and does no masking.
- **ploidy_het_per_mb_threshold** (default 1000): `custom_het_script` calls a strain diploid above this many balanced het SNPs per callable Mb.
- **ploidy_min_callable_fraction** (default 0.25): if callable bases are below this fraction of the reference length, `custom_het_script` reports `unknown`.
- **qc_min_gq** (default 20), **qc_min_dp** (default 5), **qc_min_ab** (default 0.2), **qc_hap_min_af** (default 0.8), **qc_dp_max_factor** (default 2), **qc_max_missing** (default 0.1), **qc_min_maf** (default 0.05): VARIANT_QC_FILTER thresholds; see "Variant QC" above.
- **population_mode** (default `subset`): `subset` or `regenotype`; see "Populations" above.
- **output_prefix** (default none): if set, outputs are named `<output_prefix>.<pop>.*`.
- **skip_tree** (default false), **tree_max_missing** (default 0), **tree_model** (default `GTR+ASC`; must include `+ASC`, because the alignment has only variable sites), **tree_bootstraps** (default 1000; 0 turns off UFBoot and SH-aLRT, which need at least 4 sequences), **tree_ref_name** (default `reference`): strain tree settings; see "Strain tree" above.
- **outdir** (default: `${launchDir}/results`, i.e. `results` relative to wherever you run `nextflow run` from): Output directory for results.

### SnpEff Database

A pre-built SnpEff database for R. mucilaginosa DH4148 is available at:
```
/bigdata/stajichlab/shared/lib/snpeff_db/RmucDH4148
```

Use `--snpeff_db_dir /bigdata/stajichlab/shared/lib/snpeff_db/RmucDH4148` and `--snpeff_genome_name RmucDH4148` for real production runs.

## Containers

Every process's container is declared centrally in `conf/modules.config`, not
inside the module files. By default each resolves to a pre-pulled local image
under `/bigdata/stajichlab/shared/singularity_cache` (override the root with
`--singularity_cache_dir`). Pulling `docker://` references at runtime trips a
singularity-ce 3.9.3 bug on this cluster, which is why local images are the
default here.

On a machine without that shared image store, run with `-profile container_pull`
to pull each image from its upstream registry instead.

`assets/container_checksums.sha256` records sha256 checksums of the exact local
images the pipeline executes. Verify them with:

```bash
cd /bigdata/stajichlab/shared/singularity_cache
sha256sum -c /path/to/nf_genotype_population/assets/container_checksums.sha256
```

## Tests

**Integration test (end to end, real SLURM jobs):**

```bash
tests/run_integration_test.sh
```

This builds throwaway synthetic fixtures and a throwaway SnpEff database, runs
`-entry PLOIDY_ONLY`, derives `ploidy_overrides.csv` from its output exactly as
the Phase 2 instructions above describe, runs the full default workflow, and
then asserts on `results/all.qc.annotated.vcf.gz` and
`results/all.snps.maf.annotated.vcf.gz` — record counts, `ANN=` on every
record, only `PASS` records, no record in the fixture mask, no called genotype
below the GQ/DP floors, the filter stats, the SNP alignment and the IQ-TREE
treefile (bootstraps off: the fixture has only 3 sequences), both sample
columns, the `test_subset` group cut from `all` (2 of its 3 listed strains are
in the VCF), and
(most importantly) that the haploid strain emits single-allele genotypes while
the diploid strain emits two-allele genotypes at a shared site. Any failed
assertion exits non-zero. It takes roughly 10-20 minutes, mostly SLURM queueing.

**Unit tests** (pure Python, no cluster needed):

```bash
python3 -m pytest tests/ -v
```

**Module smoke tests** (one SLURM run each; each asserts on real output content):

```bash
nextflow run tests/smoke_custom_het_ploidy.nf -profile hpcc
nextflow run tests/smoke_ploidy_inference.nf  -profile hpcc
nextflow run tests/smoke_haplotypecaller.nf   -profile hpcc
nextflow run tests/smoke_joint_genotyping.nf  -profile hpcc
# smoke_snpeff.nf needs a prebuilt SnpEff DB; run_integration_test.sh builds one.
```

## Phase 2+ (not yet implemented)

See `Rhodotorula_mucilaginosa_DH4148_ref/docs/superpowers/specs/2026-09-12-nf-genotype-population-design.md`
for nQuire/nQuack ploidy methods, population-sliced joint genotyping beyond
`all`, mosdepth/CNV visualization, empirical population-structure proposal,
and full CNV calling.
