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
    --snpeff_genome_name RmucDH4148
```

#### Outputs

- `results/all.annotated.vcf.gz` — the final annotated, hard-filtered joint callset. (Note: its `.tbi` index is not currently published to `results/` — regenerate with `tabix -p vcf results/all.annotated.vcf.gz` if needed.)
- `results/gvcfs/<strain>.g.vcf.gz` (+ `.tbi`) — **per-strain GVCFs, retained.** Design spec §5 keeps these so new strains can be added and the population re-genotyped without recalling every existing strain from CRAM. They are copied out of `work/`, so they survive a scratch cleanup.
- `results/ploidy_review.csv` and `results/ploidy_inference_all.csv` — the ploidy QC record. Ploidy inference runs on every strain in the full workflow too, so each production run leaves a diagnostic trail. Compare `ploidy_review.csv` against the `ploidy_overrides.csv` you supplied to spot ploidy drift.

#### Populations

The built-in `all` population is derived at runtime from **every strain that
produced a GVCF in the run**, not from `assets/population_sets.yaml`. Adding a
strain to `--cram_dir` and `ploidy_overrides.csv` is therefore sufficient to get
it into the `all` callset; no file needs to be kept in sync. The YAML is
reserved for named sub-population slices (Phase 2); any `all` key in it is
ignored.

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
then asserts on the final `results/all.annotated.vcf.gz` — record count, `ANN=`
on every record, both `PASS` and named FILTER values, both sample columns, and
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
