#!/usr/bin/env python3
import exact_tests
import math
import pytest

# First four test cases are from the R HardyWeinberg package's vignettes.
# Tests #5-8 are from https://github.com/jeremymcrae/snphwe .
# Entries are of the form (hom1, hets, hom2, alternative, midp, p_want).
r_HWE_cases = [
    (298, 489, 213, "two-sided", True, 0.6330965),
    (230, 314, 107, "two-sided", True, 0.9676001),
    (298, 489, 213, "two-sided", False, 0.6556635),
    (230, 314, 107, "two-sided", False, 1.0),
    (10, 500, 5000, "two-sided", False, 0.65157189991),
    (20, 1000, 5000, "two-sided", False, 1.26598491e-5),
    (200000, 200000, 200000, "two-sided", False, 0.0),
    (495000, 4990, 10, "two-sided", False, 0.570223198305),
    ]


def test_HWE():
    for test_case in r_HWE_cases:
        pval = exact_tests.HWE_test(test_case[0], test_case[1], test_case[2], alternative=test_case[3], midp=test_case[4])
        logp = exact_tests.HWE_test(test_case[0], test_case[1], test_case[2], alternative=test_case[3], midp=test_case[4], logp=True)
        assert pval == pytest.approx(test_case[5], rel=1e-6, abs=0), str(test_case)
        assert math.exp(logp) == pytest.approx(pval, rel=1e-13, abs=0), str("logp")
    # todo: test exception-throwing cases
