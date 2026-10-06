# tests/test_variant_qc_options.py
# Decisions of 2026-10-05: hom-ref genotypes have their own DP threshold
# (gVCF block MIN_DP), and IQ-TREE keeps identical sequences.


def read(path):
    with open(path) as fh:
        return fh.read()


def test_homref_dp_threshold_is_separate():
    text = read("bin/filter_population_vcf.sh")
    assert 'FMT/DP<${MIN_DP} & GT!=\\"ref\\"' in text
    assert 'FMT/DP<${MIN_DP_HOMREF} & GT=\\"ref\\"' in text
    assert "MIN_DP_HOMREF=0" in text
    assert "--min-dp-homref ${params.qc_min_dp_homref}" in read("modules/local/variant_qc_filter/main.nf")
    assert "qc_min_dp_homref   = 0" in read("nextflow.config")


def test_iqtree_keeps_identical_sequences():
    assert "-keep-ident" in read("modules/local/snp_tree/main.nf")
