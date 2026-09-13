import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from ploidy_crosscheck import crosscheck


def test_agree():
    rows = [{"strain": "S1", "inferred_ploidy": "haploid"}]
    meta = {"S1": "haploid"}
    assert crosscheck(rows, meta)[0]["status"] == "AGREE"


def test_disagree():
    rows = [{"strain": "S1", "inferred_ploidy": "diploid"}]
    meta = {"S1": "haploid"}
    assert crosscheck(rows, meta)[0]["status"] == "DISAGREE"


def test_no_metadata():
    rows = [{"strain": "S1", "inferred_ploidy": "haploid"}]
    meta = {}
    assert crosscheck(rows, meta)[0]["status"] == "NO_METADATA"


def test_haploid_from_hybrid_normalizes_to_haploid():
    rows = [{"strain": "S1", "inferred_ploidy": "haploid"}]
    meta = {"S1": "haploid_from_hybrid"}
    assert crosscheck(rows, meta)[0]["status"] == "AGREE"


def test_unknown_inference():
    rows = [{"strain": "S1", "inferred_ploidy": "unknown"}]
    meta = {"S1": "diploid"}
    assert crosscheck(rows, meta)[0]["status"] == "UNKNOWN_INFERENCE"
