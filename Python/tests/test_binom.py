#!/usr/bin/env python3
import exact_tests
import math
import numpy as np
import pytest
import sys

DBL_MIN = sys.float_info.min

# See also ../utils/pbinom_accuracy.py and ../utils/pbinom_benchmark.py .

# These test cases were scraped from
#   scipy/stats/tests/test_{discrete_distns,multivariate}.py
# (obvious todo: add more tests)
# Entries are of the form (k, n, p, p_want)
#
# SciPy has the BSD-3-Clause license.
scipy_dbinom_cases = [
    (996, 1000, 0.01, 0.0),
    (3, 7, 0.3, 0.2268945),
    (6, 14, 0.1, 0.0012926930316300004),
    ]

# These test cases were scraped from
#   scipy/stats/tests/test_morestats.py
# which were in turn mostly(?) gathered from R 3.6.2 (Dec 2019).
# Entries are of the form (k, n, p, alternative, p_want, rtol).
# rtol=1e-13 iff p_want was recalculated with MPFR.
scipy_binomtest_cases = [
    # R is obviously wrong here, k=10080000 has higher likelihood!
    # (10079999, 21000000, 0.48, "two-sided", 1.0, 1e-10),
    (10079999, 21000000, 0.48, "two-sided", 0.9998257478590489, 1e-13),
    (10079990, 21000000, 0.48, "two-sided", 0.9966892187964642, 1e-13),
    # R is wrong here as well: it treats k=10079991 as
    # equal-or-lesser-likelihood than k=10080009 when it actually has
    # ~0.00000689% greater likelihood.
    # (10080009, 21000000, 0.48, "two-sided", 0.9970377203856, 1e-10),
    (10080009, 21000000, 0.48, "two-sided", 0.9968634695850518, 1e-13),
    (10080017, 21000000, 0.48, "two-sided", 0.9940754817328065, 1e-13),
    (9, 21, 0.48, "two-sided", 0.6689672431938855, 1e-13),
    (4, 21, 0.48, "two-sided", 0.008139563452105921, 1e-13),
    (11, 21, 0.48, "two-sided", 0.8278629664608204, 1e-13),
    (7, 21, 0.48, "two-sided", 0.19667729017182273, 1e-13),
    (3, 10, 0.5, "two-sided", 0.34375, 1e-13),
    (2, 2, 0.4, "two-sided", 0.16000000000000003, 1e-13),
    (2, 4, 0.3, "two-sided", 0.5884, 1e-13),
    (484, 967, 0.5, "two-sided", 1.0, 1e-13),
    (3, 47, 3/47, "two-sided", 1.0, 1e-13),
    (13, 46, 13/46, "two-sided", 1.0, 1e-13),
    (15, 44, 15/44, "two-sided", 1.0, 1e-13),
    (7, 13, 0.5, "two-sided", 1.0, 1e-13),
    (6, 11, 0.5, "two-sided", 1.0, 1e-13),
    (20, 100, 0.25, "less", 0.14883105044299125, 1e-13),
    (20, 100, 0.25, "greater", 0.900469589894686, 1e-13),
    (20, 100, 0.25, "two-sided", 0.2983720970095647, 1e-13),
    (3, 50, 0.2, "less", 0.005656361012155304, 1e-13),
    (3, 50, 0.2, "greater", 0.9987145850467567, 1e-13),
    (3, 50, 0.2, "two-sided", 0.011917135597017913, 1e-13),
    (0, 10, 0.25, "less", 0.056313514709472656, 1e-13),
    (0, 10, 0.25, "greater", 1.0, 1e-13),
    (0, 10, 0.25, "two-sided", 0.07604122161865234, 1e-13),
    (10, 10, 0.25, "less", 1.0, 1e-13),
    (10, 10, 0.25, "greater", 9.5367431640625e-07, 1e-13),
    (10, 10, 0.25, "two-sided", 9.5367431640625e-07, 1e-13),
    (4, 16, 0.25, "two-sided", 1.0, 1e-13),
    (450, 501, 0.1, "two-sided", 0.0, 1e-13),
    (450, 501, 0.125, "two-sided", 0.0, 1e-13),
    (450, 501, 0.15, "two-sided", 1.0159969301992047e-304, 1e-13),
    (450, 501, 0.175, "two-sided", 2.975241857214944e-275, 1e-13),
    (450, 501, 0.2, "two-sided", 7.766838292252921e-250, 1e-13),
    (450, 501, 0.45, "two-sided", 2.338125092516658e-99, 1e-13),
    (450, 501, 0.5, "two-sided", 7.828459158733017e-81, 1e-13),
    (450, 501, 0.55, "two-sided", 9.915594781996304e-65, 1e-13),
    (450, 501, 0.6, "two-sided", 2.8729390725175505e-50, 1e-13),
    (450, 501, 0.65, "two-sided", 1.7175066298389118e-37, 1e-13),
    (450, 501, 0.85, "two-sided", 0.0021070691951093705, 1e-13),
    (450, 501, 0.875, "two-sided", 0.12044570587262306, 1e-13),
    (450, 501, 0.9, "two-sided", 0.8815476317480244, 1e-13),
    (450, 501, 0.925, "two-sided", 0.02712099306312912, 1e-13),
    (450, 501, 0.95, "two-sided", 2.6102587134694674e-06, 1e-13),
    (50, 100, 0.1, "two-sided", 5.832038785734363e-24, 1e-13),
    ]


def test_dbinom():
    for test_case in scipy_dbinom_cases:
        pval = exact_tests.dbinom(test_case[0], test_case[1], test_case[2])
        assert pval == pytest.approx(test_case[3], rel=1e-15, abs=DBL_MIN), str(test_case)
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
    # broadcast tests
    assert exact_tests.dbinom([[[0, 1, 2], [3, 4, 5]]], 10) == pytest.approx(np.array([[[0.0009765625, 0.009765625, 0.0439453125], [0.1171875, 0.205078125, 0.24609375]]]), rel=1e-15, abs=0)
    assert exact_tests.binom.pmf(1, 3, [[[0, 0.2, 0.4], [0.6, 0.8, 1]]]) == pytest.approx(np.array([[[0, 0.384, 0.432], [0.288, 0.096, 0]]]), rel=1e-15, abs=0)
    assert exact_tests.binom.pmf([0, 1, 2], 2, [[0], [0.25], [0.5], [1]]) == pytest.approx(np.array([[1, 0, 0], [0.5625, 0.375, 0.0625], [0.25, 0.5, 0.25], [0, 0, 1]]), rel=1e-15, abs=0)
    assert exact_tests.binom.pmf([0, 1, 2], [2, 3, 4], [[0], [0.25], [0.5], [1]]) == pytest.approx(np.array([[1, 0, 0], [0.5625, 0.421875, 0.2109375], [0.25, 0.375, 0.375], [0, 0, 0]]), rel=1e-15, abs=0)
    # todo: test exception-throwing cases


def test_pbinom():
    for test_case in scipy_binomtest_cases:
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
    assert exact_tests.pbinom(2**50 - 555, 2**51, 0.499999) == pytest.approx(1, rel=1e-15, abs=0)
    # broadcast tests
    assert exact_tests.pbinom([[[0, 1, 2], [3, 4, 5]]], 5) == pytest.approx(np.array([[[0.03125, 0.1875, 0.5], [0.8125, 0.96875, 1]]]), rel=1e-15, abs=0)
    assert exact_tests.binom.cdf([0, 1, 2], 2, [[0], [0.25], [0.5], [1]]) == pytest.approx(np.array([[1, 1, 1], [0.5625, 0.9375, 1], [0.25, 0.75, 1], [0, 0, 1]]), rel=1e-15, abs=0)
    assert exact_tests.binom.cdf([0, 1, 2], [2, 3, 4], [[0], [0.25], [0.5], [1]]) == pytest.approx(np.array([[1, 1, 1], [0.5625, 0.84375, 0.94921875], [0.25, 0.5, 0.6875], [0, 0, 0]]), rel=1e-15, abs=0)
    # todo: test exception-throwing cases


def test_qbinom():
    for test_case in scipy_binomtest_cases:
        if test_case[3] == "less":
            # Precise inversion of pbinom().
            pval = exact_tests.pbinom(test_case[0], test_case[1], test_case[2])
            assert test_case[0] == exact_tests.qbinom(pval, test_case[1], test_case[2])
            if pval < 1.0:
                assert test_case[0] + 1 == exact_tests.qbinom(pval * (1 + 0.5 ** 52), test_case[1], test_case[2])
    # Note that both R qbinom() and scipy.stats.binom.ppf() handle some of
    # edge cases differently.
    # We consistently return the minimum k in {0, 1, ..., n} for which cdf(k,
    # n, succP) >= targetP.
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
    # broadcast tests
    np.testing.assert_array_equal(exact_tests.qbinom([[[0, 0.2, 0.4], [0.6, 0.8, 1]]], 7), np.array([[[0, 2, 3], [4, 5, 7]]]))
    np.testing.assert_array_equal(exact_tests.binom.ppf([0, 0.5, 1], 4, [[0], [0.25], [0.5], [1]]), np.array([[0, 0, 0], [0, 1, 4], [0, 2, 4], [0, 4, 4]]))
    np.testing.assert_array_equal(exact_tests.binom.ppf([0, 0.5, 1], [4, 5, 6], [[0], [0.25], [0.5], [1]]), np.array([[0, 0, 0], [0, 1, 6], [0, 2, 6], [0, 5, 6]]))
    # todo: test exception-throwing cases


def test_binomtest():
    for test_case in scipy_binomtest_cases:
        pval = exact_tests.binomtest(test_case[0], test_case[1], test_case[2], alternative=test_case[3])
        assert pval == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN), str(test_case)
    # possible todo: port test_binomtest2, test_binomtest3

    assert exact_tests.binomtest(3999, 10000, 0.4, alternative="less", midp=True) == pytest.approx(0.49212893433522414, rel=1e-13, abs=0)
    assert exact_tests.binomtest(4000, 10000, 0.4, alternative="less", midp=True) == pytest.approx(0.5002714161582071, rel=1e-13, abs=0)
    assert exact_tests.binomtest(4001, 10000, 0.4, alternative="less", midp=True) == pytest.approx(0.5084135588241734, rel=1e-13, abs=0)
    assert exact_tests.binomtest(3999, 10000, 0.4, alternative="greater", midp=True) == pytest.approx(0.5078710656647758, rel=1e-13, abs=0)
    assert exact_tests.binomtest(4000, 10000, 0.4, alternative="greater", midp=True) == pytest.approx(0.4997285838417929, rel=1e-13, abs=0)
    assert exact_tests.binomtest(4001, 10000, 0.4, alternative="greater", midp=True) == pytest.approx(0.4915864411758265, rel=1e-13, abs=0)
    assert exact_tests.binomtest(0, 8, midp=True) == pytest.approx(1/256, rel=1e-13, abs=0)
    assert exact_tests.binomtest(0, 0) == 1.0
    assert exact_tests.binomtest(0, 0, midp=True) == 0.5
    assert exact_tests.binomtest(0, 2, 0.0) == 1.0
    assert exact_tests.binomtest(0, 2, 0.0, midp=True) == 0.5
    assert exact_tests.binomtest(0, 2, 1.0) == 0.0
    assert exact_tests.binomtest(0, 2, 1.0, midp=True) == 0.0
    assert math.isnan(exact_tests.binomtest(0, 2, 1.0, logp=True))
    assert math.isnan(exact_tests.binomtest(0, 2, 1.0, midp=True, logp=True))
    assert exact_tests.binomtest(1, 2, 0.0) == 0.0
    assert exact_tests.binomtest(1, 2, 0.0, midp=True) == 0.0
    assert math.isnan(exact_tests.binomtest(1, 2, 0.0, logp=True))
    assert math.isnan(exact_tests.binomtest(1, 2, 0.0, midp=True, logp=True))
    assert exact_tests.binomtest(1, 2, 1.0) == 0.0
    assert exact_tests.binomtest(1, 2, 1.0, midp=True) == 0.0
    assert math.isnan(exact_tests.binomtest(1, 2, 1.0, logp=True))
    assert math.isnan(exact_tests.binomtest(1, 2, 1.0, midp=True, logp=True))
    assert exact_tests.binomtest(2, 2, 0.0) == 0.0
    assert exact_tests.binomtest(2, 2, 0.0, midp=True) == 0.0
    assert math.isnan(exact_tests.binomtest(2, 2, 0.0, logp=True))
    assert math.isnan(exact_tests.binomtest(2, 2, 0.0, midp=True, logp=True))
    assert exact_tests.binomtest(2, 2, 1.0) == 1.0
    assert exact_tests.binomtest(2, 2, 1.0, midp=True) == 0.5

    assert exact_tests.binomtest(0, 0, alternative="less") == 1.0
    assert exact_tests.binomtest(0, 0, alternative="less", midp=True) == 0.5
    assert exact_tests.binomtest(0, 2, 0.0, alternative="less") == 1.0
    assert exact_tests.binomtest(0, 2, 0.0, alternative="less", midp=True) == 0.5
    assert exact_tests.binomtest(0, 2, 1.0, alternative="less") == 0.0
    assert exact_tests.binomtest(0, 2, 1.0, alternative="less", midp=True) == 0.0
    assert math.isnan(exact_tests.binomtest(0, 2, 1.0, alternative="less", logp=True))
    assert math.isnan(exact_tests.binomtest(0, 2, 1.0, alternative="less", midp=True, logp=True))
    assert exact_tests.binomtest(1, 2, 0.0, alternative="less") == 1.0
    assert exact_tests.binomtest(1, 2, 0.0, alternative="less", midp=True) == 1.0
    assert exact_tests.binomtest(1, 2, 1.0, alternative="less") == 0.0
    assert exact_tests.binomtest(1, 2, 1.0, alternative="less", midp=True) == 0.0
    assert math.isnan(exact_tests.binomtest(1, 2, 1.0, alternative="less", logp=True))
    assert math.isnan(exact_tests.binomtest(1, 2, 1.0, alternative="less", midp=True, logp=True))
    assert exact_tests.binomtest(2, 2, 0.0, alternative="less") == 1.0
    assert exact_tests.binomtest(2, 2, 0.0, alternative="less", midp=True) == 1.0
    assert exact_tests.binomtest(2, 2, 1.0, alternative="less") == 1.0
    assert exact_tests.binomtest(2, 2, 1.0, alternative="less", midp=True) == 0.5

    assert exact_tests.binomtest(0, 0, alternative="greater") == 1.0
    assert exact_tests.binomtest(0, 0, alternative="greater", midp=True) == 0.5
    assert exact_tests.binomtest(0, 2, 0.0, alternative="greater") == 1.0
    assert exact_tests.binomtest(0, 2, 0.0, alternative="greater", midp=True) == 0.5
    assert exact_tests.binomtest(0, 2, 1.0, alternative="greater") == 1.0
    assert exact_tests.binomtest(0, 2, 1.0, alternative="greater", midp=True) == 1.0
    assert exact_tests.binomtest(1, 2, 0.0, alternative="greater") == 0.0
    assert exact_tests.binomtest(1, 2, 0.0, alternative="greater", midp=True) == 0.0
    assert math.isnan(exact_tests.binomtest(1, 2, 0.0, alternative="greater", logp=True))
    assert math.isnan(exact_tests.binomtest(1, 2, 0.0, alternative="greater", midp=True, logp=True))
    assert exact_tests.binomtest(1, 2, 1.0, alternative="greater") == 1.0
    assert exact_tests.binomtest(1, 2, 1.0, alternative="greater", midp=True) == 1.0
    assert exact_tests.binomtest(2, 2, 0.0, alternative="greater") == 0.0
    assert exact_tests.binomtest(2, 2, 0.0, alternative="greater", midp=True) == 0.0
    assert math.isnan(exact_tests.binomtest(2, 2, 0.0, alternative="greater", logp=True))
    assert math.isnan(exact_tests.binomtest(2, 2, 0.0, alternative="greater", midp=True, logp=True))
    assert exact_tests.binomtest(2, 2, 1.0, alternative="greater") == 1.0
    assert exact_tests.binomtest(2, 2, 1.0, alternative="greater", midp=True) == 0.5

    assert exact_tests.binomtest(0, 1022) == pytest.approx(0.5 ** 1021, rel=1e-13, abs=0)
    # Tests are intended to be agnostic to denormal-flushing behavior.
    #
    # In this case, true value of exact_tests.binomtest(0, 1023) is DBL_MIN,
    # we return 0 in my testing since computed log-pvalue is a tad smaller than
    # log(DBL_MIN) and we flush the would-be denormal to zero, that's
    # acceptable.
    #
    # Returning a number within a small relative error of DBL_MIN would also be
    # acceptable.  But I'd rather not allow returning something like 1.75 *
    # DBL_MIN, so I've set up this test accordingly.
    pval = exact_tests.binomtest(0, 1023)
    if pval != 0.0:
        assert pval == pytest.approx(0.5 ** 1022, rel=1e-13, abs=0)

    # True value is small enough that zero should be returned even if we don't
    # flush denormals.
    assert exact_tests.binomtest(0, 1077) == 0.0

    # huge-magnitude log
    assert exact_tests.binomtest(0, 999999999, logp=True) == pytest.approx(-999999998 * math.log(2), rel=1e-13, abs=0)
    assert exact_tests.binomtest(1, 999999999, logp=True) == pytest.approx(-999999998 * math.log(2) + math.log(1000000000), rel=1e-13, abs=0)
    assert exact_tests.binomtest(2.9e12, 5e12, 0.52, logp=True) == pytest.approx(-36262022199.691765, rel=1e-13, abs=0)
    # tiny-magnitude log, unimportant case but may as well capture that we get
    # it right
    assert exact_tests.binomtest(5548, 9999, 0.37, alternative="less", logp=True) == pytest.approx(-8.201532972018594e-308, rel=1e-13, abs=0)
    assert exact_tests.binomtest(5549, 9999, 0.37, alternative="less", logp=True) == pytest.approx(-3.8607083741381037e-308, rel=1e-13, abs=0)
    # accept either denormal or flush-to-zero
    assert exact_tests.binomtest(5550, 9999, 0.37, alternative="less", logp=True) == pytest.approx(0.0, abs=DBL_MIN)

    assert exact_tests.binomtest(9998, 9999, alternative="less", logp=True) == 0.0
    # todo: test exception-throwing cases
