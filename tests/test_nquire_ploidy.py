import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from nquire_ploidy import parse_lrdmodel_output, call_ploidy_from_deltas


LRDMODEL_STDOUT_CONFIDENT_DIPLOID = (
    "file\tfree\tdip\ttri\ttet\tddip\tdtri\tdtet\n"
    "base.bin\t-100.0\t-100.02\t-105.0\t-110.0\t0.02\t5.0\t10.0\n"
)

LRDMODEL_STDOUT_DEGENERATE_ZERO = (
    "file\tfree\tdip\ttri\ttet\tddip\tdtri\tdtet\n"
    "base.bin\t-50.0\t-50.0\t-50.0\t-50.0\t0.000000\t0.000000\t0.000000\n"
)

LRDMODEL_STDOUT_NAN = (
    "file\tfree\tdip\ttri\ttet\tddip\tdtri\tdtet\n"
    "base.bin\tnan\tnan\tnan\tnan\tnan\tnan\tnan\n"
)


def test_parse_lrdmodel_output_extracts_deltas():
    deltas = parse_lrdmodel_output(LRDMODEL_STDOUT_CONFIDENT_DIPLOID)
    assert deltas == {"diploid_delta": 0.02, "triploid_delta": 5.0, "tetraploid_delta": 10.0}


def test_parse_lrdmodel_output_none_when_no_data_line():
    assert parse_lrdmodel_output("file\tfree\tdip\ttri\ttet\tddip\tdtri\tdtet\n") is None


def test_call_ploidy_confident_diploid_fit():
    deltas = parse_lrdmodel_output(LRDMODEL_STDOUT_CONFIDENT_DIPLOID)
    ploidy, best_fit = call_ploidy_from_deltas(deltas)
    assert ploidy == "diploid"
    assert best_fit == "diploid"


def test_call_ploidy_none_deltas_is_unknown():
    ploidy, best_fit = call_ploidy_from_deltas(None)
    assert ploidy == "unknown"
    assert best_fit is None


def test_call_ploidy_degenerate_all_equal_deltas_is_unknown():
    # CONFIRMED FAILURE MODE (found by running the real haploid test fixture):
    # a sample with too few informative sites - exactly what a true haploid
    # sample looks like to nQuire, which has no haploid model - makes
    # lrdmodel report identical deltas across all three fixed models. That
    # must be rejected as uninformative, not accepted as a confident diploid
    # fit (0.0 <= MAX_CONFIDENT_DIPLOID_DELTA would otherwise wrongly pass).
    deltas = parse_lrdmodel_output(LRDMODEL_STDOUT_DEGENERATE_ZERO)
    ploidy, best_fit = call_ploidy_from_deltas(deltas)
    assert ploidy == "unknown"


def test_call_ploidy_nan_deltas_is_unknown():
    deltas = parse_lrdmodel_output(LRDMODEL_STDOUT_NAN)
    ploidy, best_fit = call_ploidy_from_deltas(deltas)
    assert ploidy == "unknown"
