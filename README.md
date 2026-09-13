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
NXF_SYNTAX_PARSER=v1 nextflow run main.nf -entry PLOIDY_ONLY \
    --cram_dir /path/to/sarek/results/preprocessing/*/*.cram \
    --reference /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna \
    --reference_fai /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna.fai \
    --metadata /path/to/metadata.txt
```

**Important:** This cluster's Nextflow (26.04.6) requires the environment variable `NXF_SYNTAX_PARSER=v1` set before the command. Without it, Nextflow errors with "The `-entry` option is not supported with the strict parser".

This produces `results/ploidy_review.csv`, with columns `strain`, `inferred_ploidy` (haploid/diploid), and `inference_status` (AGREE/DISAGREE/NO_METADATA/UNKNOWN_INFERENCE).

### Phase 2: Manual Review and Full Genotyping

1. **Review** the `ploidy_review.csv` file.
2. **Create** a `ploidy_overrides.csv` file with header `strain,ploidy` and values `haploid` or `diploid` for any strains you want to override the inference. For strains not listed, the inferred ploidy is used.
3. **Run** the full pipeline:

```bash
nextflow run main.nf -profile hpcc \
    --cram_dir /path/to/sarek/results/preprocessing/*/*.cram \
    --reference /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna \
    --reference_fai /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna.fai \
    --reference_dict /path/to/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.dict \
    --metadata /path/to/metadata.txt \
    --ploidy_overrides /path/to/ploidy_overrides.csv \
    --population_sets assets/population_sets.yaml \
    --snpeff_db_dir /path/to/snpeff_db \
    --snpeff_genome_name RmucDH4148
```

This produces `results/all.annotated.vcf.gz` (and its `.tbi` index).

### Parameters

- **cram_dir** (required): Directory containing CRAM files and index files from Sarek preprocessing (e.g., `/path/to/sarek/results/preprocessing/*/*.cram`).
- **reference** (required for both phases): Reference genome FASTA file.
- **reference_fai** (required for both phases): Reference genome FAI index file.
- **reference_dict** (required for Phase 2 only): Reference genome GATK sequence dictionary file (`.dict`).
- **metadata** (required for both phases): Metadata file mapping strain names to attributes (e.g., `metadata.txt`).
- **ploidy_overrides** (required for Phase 2 only): CSV file with columns `strain,ploidy` containing manual ploidy assignments for strains. Generated from `ploidy_review.csv` after human review.
- **population_sets** (default: `assets/population_sets.yaml`): YAML file defining population groups. The default includes a built-in `all` group covering every confirmed R. mucilaginosa strain.
- **snpeff_db_dir** (required for Phase 2 only): Directory containing SnpEff database files.
- **snpeff_genome_name** (required for Phase 2 only): SnpEff genome database name (e.g., `RmucDH4148`).
- **outdir** (default: `results`): Output directory for results.

### SnpEff Database

A pre-built SnpEff database for R. mucilaginosa DH4148 is available at:
```
/bigdata/stajichlab/shared/lib/snpeff_db/RmucDH4148
```

Use `--snpeff_db_dir /bigdata/stajichlab/shared/lib/snpeff_db` and `--snpeff_genome_name RmucDH4148` for real production runs.

## Phase 2+ (not yet implemented)

See `Rhodotorula_mucilaginosa_DH4148_ref/docs/superpowers/specs/2026-09-12-nf-genotype-population-design.md`
for nQuire/nQuack ploidy methods, population-sliced joint genotyping beyond
`all`, mosdepth/CNV visualization, empirical population-structure proposal,
and full CNV calling.
