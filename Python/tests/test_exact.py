#!/usr/bin/env python3
import exact_tests
import pytest

# These test cases were scraped from
# scipy/stats/tests/data/fisher_exact_results_from_r.py , which were in turn
# gathered from R 3.6.2 (Dec 2019).
# SciPy has the BSD-3-Clause license.
#
# Possible todo: write a function which performs the entire calculation with
# double-double arithmetic, so we have easy access to ground truth (within 1
# ULP, anyway) for these and other test cases.
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

def test_fisher():
    # 1e-13 is several hundred ULPs; good implementations should agree to that
    # relative tolerance for small cases (such as the scipy cases where no
    # table sum > 3000).
    for test_case in scipy_fisher22_cases:
        pval = exact_tests.fisher(test_case[0], alternative=test_case[1])
        assert pval == pytest.approx(test_case[2], rel=1e-13, abs=0)
    pval = exact_tests.fisher([[4, 0], [0, 4]], midp=True)
    assert pval == pytest.approx(1 / 70, rel=1e-13, abs=0)
