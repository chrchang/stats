#!/usr/bin/env python3
import exact_tests
import pytest

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
    (10079999, 21000000, 0.48, "two-sided", 1.0, 1e-10),
    (10079990, 21000000, 0.48, "two-sided", 0.9966892187965, 1e-10),
    (10080009, 21000000, 0.48, "two-sided", 0.9970377203856, 1e-10),
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


def test_binom():
    for test_case in scipy_binom_cases:
        pval = exact_tests.binom(test_case[0], test_case[1], test_case[2], alternative=test_case[3])
        assert pval == pytest.approx(test_case[4], rel=test_case[5], abs=0), str(test_case)


def test_fisher():
    # 1e-13 is several hundred ULPs; good implementations should agree to that
    # relative tolerance for small cases (such as the scipy cases where no
    # table sum > 3000).
    for test_case in scipy_fisher22_cases:
        pval = exact_tests.fisher(test_case[0], alternative=test_case[1])
        assert pval == pytest.approx(test_case[2], rel=1e-13, abs=0), str(test_case)
    pval = exact_tests.fisher([[4, 0], [0, 4]], midp=True)
    assert pval == pytest.approx(1 / 70, rel=1e-13, abs=0), "midp"
