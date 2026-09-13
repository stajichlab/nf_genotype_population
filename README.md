# nf_genotype_population

Ploidy-aware per-strain GATK4 calling, joint genotyping, and hard-filtering
for Rhodotorula population genomics on UCR HPCC. Consumes CRAM output from
an nf-core/sarek run (Sarek itself only used through alignment/dedup - see
this project's design spec).

Design spec: see `Rhodotorula_mucilaginosa_DH4148_ref/docs/superpowers/specs/2026-09-12-nf-genotype-population-design.md`.
