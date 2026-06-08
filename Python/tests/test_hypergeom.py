#!/usr/bin/env python3
import exact_tests
import math
import pytest
import sys

DBL_MIN = sys.float_info.min

# These test cases were scraped from
#   scipy/stats/tests/test_distributions.py
# Entries are of the form (k, M, n, N, p_want, rtol)
scipy_hypergeom_pmf_cases = [
    (2, 2500, 50, 500, 0.0010114963068932233, 1e-11),
    (0, 2, 1, 0, 1.0, 1e-11),
    # This corresponds to negative m21, which we don't support.
    # (1, 2, 1, 0, 0.0, 1e-11),
    (0, 2, 0, 2, 1.0, 1e-11),
    ]

# These test cases were scraped from
#   scipy/stats/tests/test_{discrete_distns,distributions}.py
# Entries are of the form (k, M, n, N, p_want, rtol)
scipy_hypergeom_cdf_cases = [
    (3, 10, 4, 5, 0.9761904761904762, 1e-15),
    (107, 10000, 3000, 215, 0.9999999997226765, 1e-15),
    (10, 10000, 3000, 215, 2.681682217692179e-21, 5e-11),
    ]
scipy_hypergeom_logcdf_cases = [
    # cdf value from R; scipy test only confirms this is in [0, 1]
    (30, 13397950, 4363, 12390, -1.3299429823230965e-17, 1e-11),
    (1, 1e7, 1e6, 5e4, -5273.3351516552739, 1e-11),
    (40, 1600, 50, 300, -7.5651488792285293e-23, 1e-11),
    (125, 1600, 250, 500, -4.2426884594741954e-12, 1e-11),
    ]
scipy_hypergeom_sf_cases = [
    (25, 10000, 3000, 215, 0.9999999999052958, 1e-15),
    (125, 10000, 3000, 215, 1.4416781705752128e-18, 5e-11),
    (2e4, 9.9e4+1.1e5, 9.9e4, 3e4, 0, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 3.8e4, 1.904153e-114, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 3.9e4, 2.752693e-66, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 4e4, 4.931217e-32, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 4.1e4, 8.265601e-11, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 4.2e4, 0.1237904, 5e-7),
    (2e4, 9.9e4+1.1e5, 9.9e4, 5e4, 1, 5e-7),
    # This corrresponds to negative m22.
    # (9, 1e5, 1e5, 10, 1, 1e-11),
    (9, 1e6, 1e5, 10, 9.9959506789896877e-11, 1e-11),
    (9, 1e7, 1e5, 10, 9.9955458497748593e-21, 1e-11),
    (9, 1e8, 1e5, 10, 9.9955053678820622e-31, 1e-11),
    (9, 1e9, 1e5, 10, 9.9955013197030621e-41, 1e-11),
    (9, 1e10, 1e5, 10, 9.9955009148852711e-51, 1e-11),
    (9, 1e11, 1e5, 10, 9.9955008744034664e-61, 1e-11),
    (9, 1e12, 1e5, 10, 9.9955008703553125e-71, 1e-11),
    (9, 1e13, 1e5, 10, 9.9955008699504913e-81, 1e-11),
    (9, 1e14, 1e5, 10, 9.9955008699100105e-91, 1e-11),
    (9, 1e15, 1e5, 10, 9.9955008699059437e-101, 1e-11),
    ]
scipy_hypergeom_logsf_cases = [
    (1e4, 1e7, 1e6, 5e4, -2239.771249920399, 1e-11),
    (1, 1600, 600, 300, -2.5665672697340924e-68, 1e-11),
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


def test_dhyper():
    for test_case in scipy_hypergeom_pmf_cases:
        p_got = exact_tests.hypergeom.pmf(test_case[0], test_case[1], test_case[2], test_case[3])
        lnp_got = exact_tests.hypergeom.logpmf(test_case[0], test_case[1], test_case[2], test_case[3])
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
    for test_case in scipy_hypergeom_cdf_cases:
        p_got = exact_tests.hypergeom.cdf(test_case[0], test_case[1], test_case[2], test_case[3])
        assert p_got == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN)
    for test_case in scipy_hypergeom_logcdf_cases:
        lnp_got = exact_tests.hypergeom.logcdf(test_case[0], test_case[1], test_case[2], test_case[3])
        assert lnp_got == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN)
    for test_case in scipy_hypergeom_sf_cases:
        p_got = exact_tests.hypergeom.sf(test_case[0], test_case[1], test_case[2], test_case[3])
        assert p_got == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN)
    for test_case in scipy_hypergeom_logsf_cases:
        lnp_got = exact_tests.hypergeom.logsf(test_case[0], test_case[1], test_case[2], test_case[3])
        assert lnp_got == pytest.approx(test_case[4], rel=test_case[5], abs=DBL_MIN)
    for approx in [False, True]:
        assert exact_tests.phyper(1e9, 4e9, 7.999e9, 3e9, approx=approx) == pytest.approx(9.6863818919688989e-05, rel=1e-8, abs=0)
        assert exact_tests.phyper(0, 0, 0, 0, approx=approx) == 1.0
        assert exact_tests.phyper(2, 2, 2, 4, approx=approx) == 1.0
        # todo: test exception-throwing cases


def test_qhyper():
    for test_case in scipy_hypergeom_cdf_cases:
        # Precise inversion of hypergeom.cdf().
        # Ok for a future scipy_hypergeom_cdf_cases entry to fail this test and
        # be exempted (e.g. we could define
        # scipy_hypergeom_invertible_cdf_cases to be subslice of
        # scipy_hypergeom_cdf_cases in the future), as long as there is a good
        # reason for precise inversion to fail and that reason is documented
        # (e.g. there are multiple values of x for which cdf(x, ...) yields the
        # same float64 pval, and the test case doesn't cover the smallest such
        # x).
        pval = exact_tests.hypergeom.cdf(test_case[0], test_case[1], test_case[2], test_case[3])
        assert test_case[0] == exact_tests.hypergeom.ppf(pval, test_case[1], test_case[2], test_case[3])
        if pval < 1.0:
            assert test_case[0] + 1 == exact_tests.hypergeom.ppf(pval * (1 + 0.5 ** 52), test_case[1], test_case[2], test_case[3])
    for test_case in scipy_hypergeom_logcdf_cases:
        ln_pval = exact_tests.hypergeom.logcdf(test_case[0], test_case[1], test_case[2], test_case[3])
        # scipy interface doesn't include inverse-logcdf, so we explicitly
        # convert to R-style parameters here.
        r_m = test_case[3]
        r_n = test_case[1] - test_case[3]
        r_k = test_case[2]
        assert test_case[0] == exact_tests.qhyper(ln_pval, r_m, r_n, r_k, logp=True)
        if ln_pval < 0.0:
            assert test_case[0] + 1 == exact_tests.qhyper(ln_pval * (1 - 0.5 ** 52), r_m, r_n, r_k, logp=True)
    assert exact_tests.qhyper(0, 0, 0, 0) == 0
    assert exact_tests.qhyper(1, 0, 0, 0) == 0
    assert exact_tests.qhyper(0, 2, 2, 4) == 2
    assert exact_tests.qhyper(1, 2, 2, 4) == 2
    assert exact_tests.qhyper(0, 1, 1, 1, logp=True) == 1
    assert exact_tests.qhyper(-0.1, 1, 1, 1, logp=True) == 1


def test_fisher():
    # 1e-13 is several hundred ULPs; preexisting good implementations should
    # agree to that relative tolerance for small cases (such as the scipy cases
    # where no table sum > 3000).
    for test_case in scipy_fisher22_cases:
        pval = exact_tests.fisher_test(test_case[0], alternative=test_case[1])
        assert pval == pytest.approx(test_case[2], rel=1e-13, abs=DBL_MIN), str(test_case)
    pval = exact_tests.fisher_test([[4, 0], [0, 4]], midp=True)
    logp = exact_tests.fisher_test([[4, 0], [0, 4]], midp=True, logp=True)
    assert pval == pytest.approx(1/70, rel=1e-13, abs=0), "midp"
    assert logp == pytest.approx(math.log(1/70), rel=1e-13, abs=0), "logp"
    # R print(phyper(1e9, 4e9, 7.999e9, 3e9), digits=17)
    # (scipy.stats.hypergeom.cdf(1e9, 11.999e9, 3e9, 4e9) is substantially
    # less accurate.)
    assert exact_tests.fisher_test([[1e9, 2e9], [3e9, 5.999e9]], alternative="less") == pytest.approx(9.6863818919688989e-05, rel=1e-8, abs=0)
    # todo: test exception-throwing cases
