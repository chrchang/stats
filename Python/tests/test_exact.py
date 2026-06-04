#!/usr/bin/env python3
import exact_tests
import math
import pytest
import sys

DBL_MIN = sys.float_info.min

# These test cases were scraped from
#   scipy/stats/tests/test_morestats.py
# which were in turn mostly(?) gathered from R 3.6.2 (Dec 2019).
# Entries are of the form (k, n, p, alternative, p_want, rtol).
#
# SciPy has the BSD-3-Clause license.
#
# Possible todo: write a function which performs the entire calculation with
# double-double arithmetic, so we have easy access to ground truth (within 1
# ULP, anyway) for these and other test cases.
scipy_binom_cases = [
    # R is obviously wrong here, k=10080000 has higher likelihood!
    # (10079999, 21000000, 0.48, "two-sided", 1.0, 1e-10),
    (10079990, 21000000, 0.48, "two-sided", 0.9966892187965, 1e-10),
    # R is wrong here as well: it treats k=10079991 as
    # equal-or-lesser-likelihood than k=10080009 when it actually has
    # ~0.00000689% greater likelihood.
    # (10080009, 21000000, 0.48, "two-sided", 0.9970377203856, 1e-10),
    (10080017, 21000000, 0.48, "two-sided", 0.9940754817328, 1e-9),
    (9, 21, 0.48, "two-sided", 0.6689672431939, 1e-10),
    (4, 21, 0.48, "two-sided", 0.0081395634521, 1e-10),
    (11, 21, 0.48, "two-sided", 0.8278629664608, 1e-10),
    (7, 21, 0.48, "two-sided", 0.1966772901718, 1e-10),
    (3, 10, 0.5, "two-sided", 0.34375, 1e-10),
    (2, 2, 0.4, "two-sided", 0.16, 1e-10),
    (2, 4, 0.3, "two-sided", 0.5884, 1e-10),
    (484, 967, 0.5, "two-sided", 1.0, 1e-10),
    (3, 47, 3/47, "two-sided", 1.0, 1e-10),
    (13, 46, 13/46, "two-sided", 1.0, 1e-10),
    (15, 44, 15/44, "two-sided", 1.0, 1e-10),
    (7, 13, 0.5, "two-sided", 1.0, 1e-10),
    (6, 11, 0.5, "two-sided", 1.0, 1e-10),
    (20, 100, 0.25, "less", 0.148831050443, 1e-12),
    (20, 100, 0.25, "greater", 0.9004695898947, 1e-12),
    (20, 100, 0.25, "two-sided", 0.2983720970096, 1e-12),
    (3, 50, 0.2, "less", 0.005656361, 1e-6),
    (3, 50, 0.2, "greater", 0.9987146, 1e-6),
    (3, 50, 0.2, "two-sided", 0.01191714, 1e-6),
    (0, 10, 0.25, "less", 0.05631351, 1e-6),
    (0, 10, 0.25, "greater", 1.0, 1e-6),
    (0, 10, 0.25, "two-sided", 0.07604122, 1e-6),
    (10, 10, 0.25, "less", 1.0, 1e-6),
    (10, 10, 0.25, "greater", 9.536743e-07, 1e-6),
    (10, 10, 0.25, "two-sided", 9.536743e-07, 1e-6),
    (4, 16, 0.25, "two-sided", 1.0, 1e-10),
    (450, 501, 0.1, "two-sided", 0.0, 1e-10),
    (450, 501, 0.125, "two-sided", 0.0, 1e-10),
    (450, 501, 0.15, "two-sided", 1.0159969301994141e-304, 1e-10),
    (450, 501, 0.175, "two-sided", 2.9752418572150531e-275, 1e-10),
    (450, 501, 0.2, "two-sided", 7.7668382922535275e-250, 1e-10),
    (450, 501, 0.45, "two-sided", 2.3381250925167094e-099, 1e-10),
    (450, 501, 0.5, "two-sided", 7.8284591587323951e-081, 1e-10),
    (450, 501, 0.55, "two-sided", 9.9155947819961383e-065, 1e-10),
    (450, 501, 0.6, "two-sided", 2.8729390725176308e-050, 1e-10),
    (450, 501, 0.65, "two-sided", 1.7175066298388421e-037, 1e-10),
    (450, 501, 0.85, "two-sided", 0.0021070691951093692, 1e-10),
    (450, 501, 0.875, "two-sided", 0.12044570587262322, 1e-10),
    (450, 501, 0.9, "two-sided", 0.88154763174802508, 1e-10),
    (450, 501, 0.925, "two-sided", 0.027120993063129286, 1e-10),
    (450, 501, 0.95, "two-sided", 2.6102587134694721e-006, 1e-10),
    (50, 100, 0.1, "two-sided", 5.8320387857343647e-024, 1e-10),
    ]

# These test cases were scraped from
#   scipy/stats/tests/test_{discrete_distns,multivariate}.py
# (obvious todo: add more tests)
# Entries are of the form (k, n, p, p_want)
scipy_dbinom_cases = [
    (996, 1000, 0.01, 0.0),
    (3, 7, 0.3, 0.2268944999999999),
    (6, 14, 0.1, 0.0012926930316300004),
    ]

# These test cases were scraped from
#   scipy/stats/tests/data/fisher_exact_results_from_r.py
# which were in turn gathered from R 3.6.2.
# Entries are of the form (table, alternative, p_want).
scipy_fisher22_cases = [
    ([[100, 2], [1000, 5]], "two-sided", 0.1300759363430016),
    ([[2, 7], [8, 2]], "two-sided", 0.02301413756522116),
    ([[5, 1], [10, 10]], "two-sided", 0.1973244147157191),
    ([[5, 15], [20, 20]], "two-sided", 0.09580440012477633),
    ([[5, 16], [16, 25]], "two-sided", 0.2697004098849359),
    ([[10, 5], [10, 1]], "two-sided", 0.1973244147157192),
    ([[10, 5], [10, 0]], "two-sided", 0.06126482213438735),
    ([[5, 0], [1, 4]], "two-sided", 0.04761904761904762),
    ([[0, 5], [1, 4]], "two-sided", 1),
    ([[5, 1], [0, 4]], "two-sided", 0.04761904761904761),
    ([[0, 1], [3, 2]], "two-sided", 1),
    ([[200, 7], [8, 300]], "two-sided", 2.005657880389071e-122),
    ([[28, 21], [6, 1957]], "two-sided", 5.728437460831947e-44),
    ([[190, 800], [200, 900]], "two-sided", 0.574111858126088),
    ([[100, 2], [1000, 5]], "less", 0.1300759363430016),
    ([[2, 7], [8, 2]], "less", 0.0185217259520665),
    ([[5, 1], [10, 10]], "less", 0.9782608695652173),
    ([[5, 15], [20, 20]], "less", 0.05625775074399956),
    ([[5, 16], [16, 25]], "less", 0.1808979350599346),
    ([[10, 5], [10, 1]], "less", 0.1652173913043479),
    ([[10, 5], [10, 0]], "less", 0.0565217391304348),
    ([[5, 0], [1, 4]], "less", 1),
    ([[0, 5], [1, 4]], "less", 0.5),
    ([[5, 1], [0, 4]], "less", 1),
    ([[0, 1], [3, 2]], "less", 0.4999999999999999),
    ([[200, 7], [8, 300]], "less", 1),
    ([[28, 21], [6, 1957]], "less", 1),
    ([[190, 800], [200, 900]], "less", 0.7416227010368963),
    ([[100, 2], [1000, 5]], "greater", 0.979790445314723),
    ([[2, 7], [8, 2]], "greater", 0.9990149169715733),
    ([[5, 1], [10, 10]], "greater", 0.1652173913043478),
    ([[5, 15], [20, 20]], "greater", 0.9849086665340765),
    ([[5, 16], [16, 25]], "greater", 0.9330176609214881),
    ([[10, 5], [10, 1]], "greater", 0.9782608695652174),
    ([[10, 5], [10, 0]], "greater", 1),
    ([[5, 0], [1, 4]], "greater", 0.02380952380952382),
    ([[0, 5], [1, 4]], "greater", 1),
    ([[5, 1], [0, 4]], "greater", 0.0238095238095238),
    ([[0, 1], [3, 2]], "greater", 1),
    ([[200, 7], [8, 300]], "greater", 2.005657880388915e-122),
    ([[28, 21], [6, 1957]], "greater", 5.728437460831983e-44),
    ([[190, 800], [200, 900]], "greater", 0.2959825901308897),
    ]

# These test cases were scraped from
#   scipy/stats/tests/test_distributions.py
# Entries are of the form (k, M, n, N, p_want, rtol)
scipy_dhyper_cases = [
    (2, 2500, 50, 500, 0.0010114963068932233, 1e-11),
    (0, 2, 1, 0, 1.0, 1e-11),
    # This corresponds to negative m21, which we don't support.
    # (1, 2, 1, 0, 0.0, 1e-11),
    (0, 2, 0, 2, 1.0, 1e-11),
    ]

# These test cases were scraped from
#   scipy/stats/tests/test_{discrete_distns,distributions}.py
# Entries are of the form (k, M, n, N, complement, logp, p_want, rtol)
scipy_phyper_cases = [
    (3, 10, 4, 5, False, False, 0.9761904761904762, 1e-15),
    (107, 10000, 3000, 215, False, False, 0.9999999997226765, 1e-15),
    (10, 10000, 3000, 215, False, False, 2.681682217692179e-21, 5e-11),
    (25, 10000, 3000, 215, True, False, 0.9999999999052958, 1e-15),
    (125, 10000, 3000, 215, True, False, 1.4416781705752128e-18, 5e-11),
    # cdf value from R; scipy test only confirms this is in [0, 1]
    (30, 13397950, 4363, 12390, False, True, -1.3299429823230965e-17, 1e-11),
    (2e4, 9.9e4+1.1e5, 9.9e4, 3e4, True, False, 0, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 3.8e4, True, False, 1.904153e-114, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 3.9e4, True, False, 2.752693e-66, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 4e4, True, False, 4.931217e-32, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 4.1e4, True, False, 8.265601e-11, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 4.2e4, True, False, 0.1237904, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 5e4, True, False, 1, 5e-7),
    (1e4, 1e7, 1e6, 5e4, True, True, -2239.771249920399, 1e-11),
    (1, 1600, 600, 300, True, True, -2.5665672697340924e-68, 1e-11),
    (1, 1e7, 1e6, 5e4, False, True, -5273.3351516552739, 1e-11),
    (40, 1600, 50, 300, False, True, -7.5651488792285293e-23, 1e-11),
    (125, 1600, 250, 500, False, True, -4.2426884594741954e-12, 1e-11),
    # This corrresponds to negative m22.
    # (9, 1e5, 1e5, 10, True, False, 1, 1e-11),
    (9, 1e6, 1e5, 10, True, False, 9.9959506789896877e-11, 1e-11),
    (9, 1e7, 1e5, 10, True, False, 9.9955458497748593e-21, 1e-11),
    (9, 1e8, 1e5, 10, True, False, 9.9955053678820622e-31, 1e-11),
    (9, 1e9, 1e5, 10, True, False, 9.9955013197030621e-41, 1e-11),
    (9, 1e10, 1e5, 10, True, False, 9.9955009148852711e-51, 1e-11),
    (9, 1e11, 1e5, 10, True, False, 9.9955008744034664e-61, 1e-11),
    (9, 1e12, 1e5, 10, True, False, 9.9955008703553125e-71, 1e-11),
    (9, 1e13, 1e5, 10, True, False, 9.9955008699504913e-81, 1e-11),
    (9, 1e14, 1e5, 10, True, False, 9.9955008699100105e-91, 1e-11),
    (9, 1e15, 1e5, 10, True, False, 9.9955008699059437e-101, 1e-11),
    ]

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


def test_binom():
    for test_case in scipy_binom_cases:
        pval = exact_tests.binom(test_case[0], test_case[1], test_case[2], alternative=test_case[3])
        assert pval == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN), str(test_case)
    # possible todo: port test_binomtest2, test_binomtest3

    assert exact_tests.binom(3999, 10000, 0.4, alternative="less", midp=True) == pytest.approx(0.49212893433522414, rel=1e-13, abs=0)
    assert exact_tests.binom(4000, 10000, 0.4, alternative="less", midp=True) == pytest.approx(0.5002714161582071, rel=1e-13, abs=0)
    assert exact_tests.binom(4001, 10000, 0.4, alternative="less", midp=True) == pytest.approx(0.5084135588241734, rel=1e-13, abs=0)
    assert exact_tests.binom(3999, 10000, 0.4, alternative="greater", midp=True) == pytest.approx(0.5078710656647758, rel=1e-13, abs=0)
    assert exact_tests.binom(4000, 10000, 0.4, alternative="greater", midp=True) == pytest.approx(0.4997285838417929, rel=1e-13, abs=0)
    assert exact_tests.binom(4001, 10000, 0.4, alternative="greater", midp=True) == pytest.approx(0.4915864411758265, rel=1e-13, abs=0)
    assert exact_tests.binom(0, 8, midp=True) == pytest.approx(1/256, rel=1e-13, abs=0)
    assert exact_tests.binom(0, 0) == 1.0
    assert exact_tests.binom(0, 0, midp=True) == 0.5
    assert exact_tests.binom(0, 2, 0.0) == 1.0
    assert exact_tests.binom(0, 2, 0.0, midp=True) == 0.5
    assert exact_tests.binom(0, 2, 1.0) == 0.0
    assert exact_tests.binom(0, 2, 1.0, midp=True) == 0.0
    assert math.isnan(exact_tests.binom(0, 2, 1.0, logp=True))
    assert math.isnan(exact_tests.binom(0, 2, 1.0, midp=True, logp=True))
    assert exact_tests.binom(1, 2, 0.0) == 0.0
    assert exact_tests.binom(1, 2, 0.0, midp=True) == 0.0
    assert math.isnan(exact_tests.binom(1, 2, 0.0, logp=True))
    assert math.isnan(exact_tests.binom(1, 2, 0.0, midp=True, logp=True))
    assert exact_tests.binom(1, 2, 1.0) == 0.0
    assert exact_tests.binom(1, 2, 1.0, midp=True) == 0.0
    assert math.isnan(exact_tests.binom(1, 2, 1.0, logp=True))
    assert math.isnan(exact_tests.binom(1, 2, 1.0, midp=True, logp=True))
    assert exact_tests.binom(2, 2, 0.0) == 0.0
    assert exact_tests.binom(2, 2, 0.0, midp=True) == 0.0
    assert math.isnan(exact_tests.binom(2, 2, 0.0, logp=True))
    assert math.isnan(exact_tests.binom(2, 2, 0.0, midp=True, logp=True))
    assert exact_tests.binom(2, 2, 1.0) == 1.0
    assert exact_tests.binom(2, 2, 1.0, midp=True) == 0.5

    assert exact_tests.binom(0, 0, alternative="less") == 1.0
    assert exact_tests.binom(0, 0, alternative="less", midp=True) == 0.5
    assert exact_tests.binom(0, 2, 0.0, alternative="less") == 1.0
    assert exact_tests.binom(0, 2, 0.0, alternative="less", midp=True) == 0.5
    assert exact_tests.binom(0, 2, 1.0, alternative="less") == 0.0
    assert exact_tests.binom(0, 2, 1.0, alternative="less", midp=True) == 0.0
    assert math.isnan(exact_tests.binom(0, 2, 1.0, alternative="less", logp=True))
    assert math.isnan(exact_tests.binom(0, 2, 1.0, alternative="less", midp=True, logp=True))
    assert exact_tests.binom(1, 2, 0.0, alternative="less") == 1.0
    assert exact_tests.binom(1, 2, 0.0, alternative="less", midp=True) == 1.0
    assert exact_tests.binom(1, 2, 1.0, alternative="less") == 0.0
    assert exact_tests.binom(1, 2, 1.0, alternative="less", midp=True) == 0.0
    assert math.isnan(exact_tests.binom(1, 2, 1.0, alternative="less", logp=True))
    assert math.isnan(exact_tests.binom(1, 2, 1.0, alternative="less", midp=True, logp=True))
    assert exact_tests.binom(2, 2, 0.0, alternative="less") == 1.0
    assert exact_tests.binom(2, 2, 0.0, alternative="less", midp=True) == 1.0
    assert exact_tests.binom(2, 2, 1.0, alternative="less") == 1.0
    assert exact_tests.binom(2, 2, 1.0, alternative="less", midp=True) == 0.5

    assert exact_tests.binom(0, 0, alternative="greater") == 1.0
    assert exact_tests.binom(0, 0, alternative="greater", midp=True) == 0.5
    assert exact_tests.binom(0, 2, 0.0, alternative="greater") == 1.0
    assert exact_tests.binom(0, 2, 0.0, alternative="greater", midp=True) == 0.5
    assert exact_tests.binom(0, 2, 1.0, alternative="greater") == 1.0
    assert exact_tests.binom(0, 2, 1.0, alternative="greater", midp=True) == 1.0
    assert exact_tests.binom(1, 2, 0.0, alternative="greater") == 0.0
    assert exact_tests.binom(1, 2, 0.0, alternative="greater", midp=True) == 0.0
    assert math.isnan(exact_tests.binom(1, 2, 0.0, alternative="greater", logp=True))
    assert math.isnan(exact_tests.binom(1, 2, 0.0, alternative="greater", midp=True, logp=True))
    assert exact_tests.binom(1, 2, 1.0, alternative="greater") == 1.0
    assert exact_tests.binom(1, 2, 1.0, alternative="greater", midp=True) == 1.0
    assert exact_tests.binom(2, 2, 0.0, alternative="greater") == 0.0
    assert exact_tests.binom(2, 2, 0.0, alternative="greater", midp=True) == 0.0
    assert math.isnan(exact_tests.binom(2, 2, 0.0, alternative="greater", logp=True))
    assert math.isnan(exact_tests.binom(2, 2, 0.0, alternative="greater", midp=True, logp=True))
    assert exact_tests.binom(2, 2, 1.0, alternative="greater") == 1.0
    assert exact_tests.binom(2, 2, 1.0, alternative="greater", midp=True) == 0.5

    assert exact_tests.binom(0, 1022) == pytest.approx(0.5 ** 1021, rel=1e-13, abs=0)
    # Tests are intended to be agnostic to denormal-flushing behavior.
    #
    # In this case, true value of exact_tests.binom(0, 1023) is DBL_MIN, we
    # return 0 in my testing since computed log-pvalue is a tad smaller than
    # log(DBL_MIN) and we flush the would-be denormal to zero, that's
    # acceptable.
    #
    # Returning a number within a small relative error of DBL_MIN would also be
    # acceptable.  But I'd rather not allow returning something like 1.75 *
    # DBL_MIN, so I've set up this test accordingly.
    pval = exact_tests.binom(0, 1023)
    if pval != 0.0:
        assert pval == pytest.approx(0.5 ** 1022, rel=1e-13, abs=0)

    # True value is small enough that zero should be returned even if we don't
    # flush denormals.
    assert exact_tests.binom(0, 1077) == 0.0

    # huge-magnitude log
    assert exact_tests.binom(0, 999999999, logp=True) == pytest.approx(-999999998 * math.log(2), rel=1e-13, abs=0)
    assert exact_tests.binom(1, 999999999, logp=True) == pytest.approx(-999999998 * math.log(2) + math.log(1000000000), rel=1e-13, abs=0)
    # tiny-magnitude log, unimportant case but may as well capture that we get
    # it right
    assert exact_tests.binom(5548, 9999, 0.37, alternative="less", logp=True) == pytest.approx(-8.201532972018594e-308, rel=1e-13, abs=0)
    assert exact_tests.binom(5549, 9999, 0.37, alternative="less", logp=True) == pytest.approx(-3.8607083741381037e-308, rel=1e-13, abs=0)
    # accept either denormal or flush-to-zero
    assert exact_tests.binom(5550, 9999, 0.37, alternative="less", logp=True) == pytest.approx(0.0, abs=DBL_MIN)

    assert exact_tests.binom(9998, 9999, alternative="less", logp=True) == 0.0
    # todo: test exception-throwing cases


def test_dbinom():
    for test_case in scipy_dbinom_cases:
        pval = exact_tests.dbinom(test_case[0], test_case[1], test_case[2])
        assert pval == pytest.approx(test_case[3], rel=1e-10, abs=DBL_MIN), str(test_case)
    assert exact_tests.dbinom(0, 0) == 1.0
    assert exact_tests.dbinom(-1, 0) == 0.0
    assert math.isnan(exact_tests.dbinom(-1, 0, logp=True))
    assert exact_tests.dbinom(1, 0) == 0.0
    assert math.isnan(exact_tests.dbinom(1, 0, logp=True))
    assert exact_tests.dbinom(0, 2, 0.0) == 1.0
    assert exact_tests.dbinom(1, 2, 0.0) == 0.0
    assert math.isnan(exact_tests.dbinom(1, 2, 0.0, logp=True))
    assert exact_tests.dbinom(2, 2, 0.0) == 0.0
    assert math.isnan(exact_tests.dbinom(2, 2, 0.0, logp=True))
    assert exact_tests.dbinom(0, 2, 1.0) == 0.0
    assert math.isnan(exact_tests.dbinom(0, 2, 1.0, logp=True))
    assert exact_tests.dbinom(1, 2, 1.0) == 0.0
    assert math.isnan(exact_tests.dbinom(1, 2, 1.0, logp=True))
    assert exact_tests.dbinom(2, 2, 1.0) == 1.0
    assert exact_tests.dbinom(0, 999999999, logp=True) == pytest.approx(-999999999 * math.log(2), rel=1e-15, abs=0)
    assert exact_tests.dbinom(1, 999999999, logp=True) == pytest.approx(-999999999 * math.log(2) + math.log(999999999), rel=1e-15, abs=0)
    # todo: test exception-throwing cases


def test_pbinom():
    for test_case in scipy_binom_cases:
        if test_case[3] == "less":
            pval = exact_tests.pbinom(test_case[0], test_case[1], test_case[2])
            assert pval == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN)
    assert exact_tests.pbinom(0, 0) == 1.0
    assert exact_tests.pbinom(-1, 0) == 0.0
    assert math.isnan(exact_tests.pbinom(-1, 0, logp=True))
    assert exact_tests.pbinom(1, 0) == 1.0
    assert exact_tests.pbinom(0, 2, 0.0) == 1.0
    assert exact_tests.pbinom(0, 2, 1.0) == 0.0
    assert math.isnan(exact_tests.pbinom(0, 2, 1.0, logp=True))
    assert exact_tests.pbinom(1, 2, 0.0) == 1.0
    assert exact_tests.pbinom(1, 2, 1.0) == 0.0
    assert math.isnan(exact_tests.pbinom(1, 2, 1.0, logp=True))
    assert exact_tests.pbinom(2, 2, 0.0) == 1.0
    assert exact_tests.pbinom(2, 2, 1.0) == 1.0
    assert exact_tests.pbinom(0, 999999999, logp=True) == pytest.approx(-999999999 * math.log(2), rel=1e-15, abs=0)
    assert exact_tests.pbinom(1, 999999999, logp=True) == pytest.approx(-999999999 * math.log(2) + math.log(1000000000), rel=1e-15, abs=0)
    assert exact_tests.pbinom(5548, 9999, 0.37, logp=True) == pytest.approx(-8.201532972018594e-308, rel=1e-15, abs=0)
    assert exact_tests.pbinom(5549, 9999, 0.37, logp=True) == pytest.approx(-3.8607083741381037e-308, rel=1e-15, abs=0)
    assert exact_tests.pbinom(5550, 9999, 0.37, logp=True) == pytest.approx(0.0, abs=DBL_MIN)
    # todo: test exception-throwing cases


def test_qbinom():
    for test_case in scipy_binom_cases:
        if test_case[3] == "less":
            # Precise inversion of pbinom().
            pval = exact_tests.pbinom(test_case[0], test_case[1], test_case[2])
            assert test_case[0] == exact_tests.qbinom(pval, test_case[1], test_case[2])
            if pval < 1.0:
                assert test_case[0] + 1 == exact_tests.qbinom(pval * (1 + 0.5 ** 52), test_case[1], test_case[2])
    assert exact_tests.qbinom(0, 0) == 0
    assert exact_tests.qbinom(1, 0) == 0
    assert exact_tests.qbinom(0, 2, 0.0) == 0
    assert exact_tests.qbinom(0.5, 2, 0.0) == 0
    assert exact_tests.qbinom(1, 2, 0.0) == 0
    assert exact_tests.qbinom(0, 2, 1.0) == 0
    assert exact_tests.qbinom(0.5, 2, 1.0) == 2
    assert exact_tests.qbinom(1, 2, 1.0) == 2
    assert exact_tests.qbinom(-1000000000 * math.log(2), 999999999, logTarget=True) == 0
    assert exact_tests.qbinom(-999999998 * math.log(2), 999999999, logTarget=True) == 1
    # todo: test exception-throwing cases


def test_fisher():
    # 1e-13 is several hundred ULPs; good implementations should agree to that
    # relative tolerance for small cases (such as the scipy cases where no
    # table sum > 3000).
    for test_case in scipy_fisher22_cases:
        pval = exact_tests.fisher(test_case[0], alternative=test_case[1])
        assert pval == pytest.approx(test_case[2], rel=1e-13, abs=DBL_MIN), str(test_case)
    pval = exact_tests.fisher([[4, 0], [0, 4]], midp=True)
    logp = exact_tests.fisher([[4, 0], [0, 4]], midp=True, logp=True)
    assert pval == pytest.approx(1/70, rel=1e-13, abs=0), "midp"
    assert logp == pytest.approx(math.log(1/70), rel=1e-13, abs=0), "logp"
    # R print(phyper(1e9, 4e9, 7.999e9, 3e9), digits=17)
    # (scipy.stats.hypergeom.cdf(1e9, 11.999e9, 3e9, 4e9) is substantially
    # less accurate.)
    assert exact_tests.fisher([[1e9, 2e9], [3e9, 5.999e9]], alternative="less") == pytest.approx(9.6863818919688989e-05, rel=1e-8, abs=0)
    # todo: test exception-throwing cases


def test_dhyper():
    for test_case in scipy_dhyper_cases:
        # scipy.stats.hypergeom.pmf(m11, m11+m12+m21+m22, m11+m12, m11+m21) is
        # equivalent to R dhyper(m11, m11+m21, m12+m22, m11+m12)
        # -> ([0], [3], [1]-[3], [2])
        a = test_case[0]
        b = test_case[3]
        c = test_case[1] - b
        d = test_case[2]
        p_got = exact_tests.dhyper(a, b, c, d)
        lnp_got = exact_tests.dhyper(a, b, c, d, logp=True)
        assert p_got == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN)
        if math.isnan(lnp_got):
            assert 0.0 == test_case[4]
        else:
            assert math.exp(lnp_got) == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN)
    assert exact_tests.dhyper(1e9, 4e9, 7.999e9, 3e9) == pytest.approx(1.7179298149371888e-08, rel=1e-8, abs=0)
    assert exact_tests.dhyper(0, 0, 0, 0) == 1.0
    assert exact_tests.dhyper(2, 2, 2, 4) == 1.0
    # todo: test exception-throwing cases


def test_phyper():
    for test_case in scipy_phyper_cases:
        a = test_case[0]
        b = test_case[3]
        c = test_case[1] - b
        d = test_case[2]
        p_got_approx = exact_tests.phyper(a, b, c, d, complement=test_case[4], approx=True, logp=test_case[5])
        p_got = exact_tests.phyper(a, b, c, d, complement=test_case[4], logp=test_case[5])
        assert p_got_approx == pytest.approx(test_case[6], rel=test_case[7], abs=DBL_MIN)
        assert p_got == pytest.approx(test_case[6], rel=test_case[7], abs=DBL_MIN)
    assert exact_tests.phyper(1e9, 4e9, 7.999e9, 3e9) == pytest.approx(9.6863818919688989e-05, rel=1e-8, abs=0)
    assert exact_tests.phyper(1e9, 4e9, 7.999e9, 3e9, approx=True) == pytest.approx(9.6863818919688989e-05, rel=1e-8, abs=0)
    assert exact_tests.phyper(0, 0, 0, 0) == 1.0
    assert exact_tests.phyper(2, 2, 2, 4) == 1.0
    # todo: test exception-throwing cases


def test_qhyper():
    for test_case in scipy_phyper_cases:
        if test_case[4] == False:
            # Precise inversion of phyper().
            a = test_case[0]
            b = test_case[3]
            c = test_case[1] - b
            d = test_case[2]
            logp = test_case[5]
            pval = exact_tests.phyper(a, b, c, d, logp=logp)
            assert test_case[0] == exact_tests.qhyper(pval, b, c, d, logp=logp)
            if not logp:
                if pval < 1.0:
                    assert test_case[0] + 1 == exact_tests.qhyper(pval * (1 + 0.5 ** 52), b, c, d)
            else:
                if pval < 0.0:
                    assert test_case[0] + 1 == exact_tests.qhyper(pval * (1 - 0.5 ** 52), b, c, d, logp=True)
    assert exact_tests.qhyper(0, 0, 0, 0) == 0
    assert exact_tests.qhyper(1, 0, 0, 0) == 0
    assert exact_tests.qhyper(0, 2, 2, 4) == 2
    assert exact_tests.qhyper(1, 2, 2, 4) == 2
    assert exact_tests.qhyper(0, 1, 1, 1, logp=True) == 1


def test_HWE():
    for test_case in r_HWE_cases:
        pval = exact_tests.HWE(test_case[0], test_case[1], test_case[2], alternative=test_case[3], midp=test_case[4])
        logp = exact_tests.HWE(test_case[0], test_case[1], test_case[2], alternative=test_case[3], midp=test_case[4], logp=True)
        assert pval == pytest.approx(test_case[5], rel=1e-6, abs=0), str(test_case)
        assert math.exp(logp) == pytest.approx(pval, rel=1e-13, abs=0), str("logp")
    # todo: test exception-throwing cases
