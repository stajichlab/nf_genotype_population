# tests/test_filter_expressions.py
import re

MODULE_PATH = "modules/local/gatk4_variantfiltration/main.nf"


def read_module():
    with open(MODULE_PATH) as fh:
        return fh.read()


def test_snp_filter_thresholds_present():
    text = read_module()
    for expr in ["QD < 2.0", "FS > 60.0", "MQ < 40.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0", "SOR > 3.0"]:
        assert expr in text, f"missing SNP filter expression: {expr}"


def test_indel_filter_thresholds_present():
    text = read_module()
    for expr in ["QD < 2.0", "FS > 200.0", "ReadPosRankSum < -20.0", "SOR > 10.0"]:
        assert expr in text, f"missing indel filter expression: {expr}"
