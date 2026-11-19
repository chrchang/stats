#!/usr/bin/env python3
import exact_tests
from fractions import Fraction
import math
import numpy as np
import pytest
import sys


def brute_pmfs(mu, kmax):
    return [math.exp(k * math.log(mu) - mu - math.lgamma(k + 1)) for k in range(kmax + 1)]


# Entries are of the form (k, mu, logp, pmf_want); pmf_want is from R 4.7's
# dpois(), which uses Loader's algorithm and is accurate to ~1 ULP.
dpois_cases = [
    (0, 2.5, False, 8.2084998623898800e-02),
    (3, 5.0, False, 1.4037389581428059e-01),
    (1000000000, 1e9, False, 1.2615662609049496e-05),
    (950, 1000.5, False, 3.5395832418587216e-03),
    (123456789, 123400000.0, False, 7.6032678540385891e-11),
    (0, 2.5, True, -2.5),
    (3, 5.0, True, -1.9634457319257539),
    (1000000000, 1e9, True, -11.2805714517612099),
    (950, 1000.5, True, -5.6437462870554693),
    (123456789, 123400000.0, True, -23.2998578872629167),
    ]

def test_dpois():
    for k, mu, logp, want in dpois_cases:
        assert exact_tests.dpois(k, mu, logp=logp) == pytest.approx(want, rel=1e-14, abs=0), (k, mu, logp)
        if logp:
            assert exact_tests.poisson.logpmf(k, mu) == pytest.approx(want, rel=1e-14, abs=0)
        else:
            assert exact_tests.poisson.pmf(k, mu) == pytest.approx(want, rel=1e-14, abs=0)

    assert exact_tests.dpois(-1, 3.0) == 0.0
    assert exact_tests.dpois(2.5, 3.0) == 0.0
    assert exact_tests.dpois(0, 0.0) == 1.0
    assert exact_tests.dpois(1, 0.0) == 0.0
    assert exact_tests.dpois(1, 0.0, logp=True) == -math.inf
    assert exact_tests.dpois(10, 2000.0) == 0.0  # underflows

    pmfs = brute_pmfs(3.7, 30)
    np.testing.assert_allclose(exact_tests.dpois(np.arange(31), 3.7), pmfs, rtol=1e-13)
    # scipy-style broadcasting over both k and mu
    got = exact_tests.poisson.pmf(np.array([[1], [2]]), np.array([0.5, 1.5, 2.5]))
    assert got.shape == (2, 3)
    assert got[1, 2] == exact_tests.dpois(2, 2.5)
    with pytest.raises(RuntimeError):
        exact_tests.dpois(1, -1.0)
    with pytest.raises(RuntimeError):
        exact_tests.dpois(1, math.inf)


# Entries are of the form (k, mu, complement, logp, cdf_want), where cdf_want
# was computed with a 320-bit MPFR implementation.
ppois_cases = [
    (3, 5.0, False, False, 0.26502591529736169),
    (3, 5.0, True, False, 0.73497408470263825),
    (0, 1e-05, False, False, 0.99999000004999983),
    (1000, 1000.5, False, False, 0.50210232422301759),
    (950, 1000.5, False, True, -2.8816172696753406),
    (1003000, 1e6, True, False, 0.0013535901291012492),
    (10, 2000.0, False, True, -1939.0903779584701),
    (5000100000, 5e9, True, True, -2.5427614870316204),
    ]

def test_ppois():
    for k, mu, complement, logp, want in ppois_cases:
        got = exact_tests.ppois(k, mu, complement=complement, logp=logp)
        assert got == pytest.approx(want, rel=1e-14, abs=0), (k, mu, complement, logp)
        if complement:
            got2 = exact_tests.poisson.logsf(k, mu) if logp else exact_tests.poisson.sf(k, mu)
        else:
            got2 = exact_tests.poisson.logcdf(k, mu) if logp else exact_tests.poisson.cdf(k, mu)
        assert got2 == got
        dist = exact_tests.Poisson(mu=mu)
        if complement:
            got3 = dist.logccdf(k) if logp else dist.ccdf(k)
        else:
            got3 = dist.logcdf(k) if logp else dist.cdf(k)
        assert got3 == got

    pmfs = brute_pmfs(37.5, 300)
    np.testing.assert_allclose(exact_tests.ppois(np.arange(101), 37.5), [math.fsum(pmfs[:k+1]) for k in range(101)], rtol=1e-13)
    np.testing.assert_allclose(exact_tests.ppois(np.arange(101), 37.5, complement=True),
                               [math.fsum(pmfs[k+1:]) for k in range(101)], rtol=1e-13)

    assert exact_tests.ppois(-1, 3.0) == 0.0
    assert exact_tests.ppois(-1, 3.0, complement=True) == 1.0
    assert exact_tests.ppois(math.inf, 3.0) == 1.0
    assert exact_tests.ppois(math.inf, 3.0, complement=True) == 0.0
    assert exact_tests.ppois(0, 0.0) == 1.0
    assert exact_tests.ppois(0, 0.0, complement=True) == 0.0
    assert exact_tests.ppois(2.9, 5.0) == exact_tests.ppois(2, 5.0)
    assert math.isnan(exact_tests.ppois(math.nan, 5.0))
    # deep tails which underflow in non-log space
    assert exact_tests.ppois(10, 2000.0) == 0.0
    assert exact_tests.ppois(10**6, 1e4, complement=True, logp=True) < -3e6
    # huge mu
    assert exact_tests.ppois(0, 1e100, logp=True) == -1e100
    with pytest.raises(RuntimeError):
        exact_tests.ppois(2.0**54, 2.0**52)


def ulp_of_exact(x):
    # Spacing between consecutive doubles in the binade containing the exact
    # value x (a Fraction).  Unlike math.ulp(float(x)), this is correct when x
    # is just below a power of 2 and float(x) rounds up to it.
    x = abs(x)
    if x < Fraction(sys.float_info.min):
        return Fraction(1, 2**1074)
    e = x.numerator.bit_length() - x.denominator.bit_length()
    if Fraction(2)**e > x:
        e -= 1
    elif Fraction(2)**(e + 1) <= x:
        e += 1
    return Fraction(2)**(e - 52)

def assert_within_ulps(got, want_str, max_ulps):
    # want_str is the exact (infinite-precision) value to >= 40 significant
    # digits.  The error bound is strict, so with max_ulps=1, only the two
    # doubles bracketing the exact value are accepted (or just one, if the
    # exact value is itself a double).
    want = Fraction(want_str)
    assert abs(Fraction(got) - want) < Fraction(max_ulps) * ulp_of_exact(want), (got, want_str)

# Entries are of the form (k, mu, complement, logp, cdf_want), where cdf_want
# is the exact value to 40 significant digits (320-bit MPFR).  These cover
# direct summation, Temme's expansion, both continued fractions, and
# subnormal results.
ppois_ulp_cases = [
    (141, 147.20167892657778, False, False, "3.2304912873868626971708020570449528795488e-01"),
    (12663, 12686.929247057338, True, False, "5.8183571812835849256984833640330551334087e-01"),
    (1076091, 1081711.834412181, False, False, "3.1792000688403860517851490342466605538187e-08"),
    (1846969520, 1846972037.0, True, False, "5.2334533361158147469107299102850825586071e-01"),
    (67534192, 67535519.68782584, False, False, "4.3585850014403995534920817754114154457697e-01"),
    (168122, 163465.91228615338, True, False, "1.0030570566323526339401009456046814545951e-30"),
    (236634, 244837.89991574993, True, True, "-1.0400949324237104381799170364576936337331e-62"),
    (1197, 3000.75, False, False, "4.8915362820691371222748901856173272531538e-308"),
    (960000, 1000000.5, False, True, "-8.1549327880972367259529500916973386110255e+02"),
    (1040000, 1000000.5, True, True, "-7.9414974396212741140676557029250080890784e+02"),
    (2, 0.0028336472934657963, False, False, "9.9999999621589469309393198721189565326894e-01"),
    ]

def test_ppois_ulps():
    for k, mu, complement, logp, want in ppois_ulp_cases:
        # Ppois() promises <1 ULP error relative to the exact value, so only the
        # two doubles bracketing it are acceptable.  (As of this writing, all
        # of these cases are within 0.5 ULP, i.e. correctly rounded, but we
        # don't want to fail if a future speed optimization gives up
        # unpromised accuracy.)  Several of these exact values are within
        # 0.001 ULP of a rounding boundary.
        assert_within_ulps(exact_tests.ppois(k, mu, complement=complement, logp=logp), want, 1)
        # PpoisApprox() is usually within a few ULPs.
        assert_within_ulps(exact_tests.ppois(k, mu, complement=complement, logp=logp, approx=True), want, 8)
    assert exact_tests.ppois(1, 30000.25) == 0.0
    assert exact_tests.ppois(1, 30000.25, approx=True) == 0.0
    dist = exact_tests.Poisson(mu=1081711.834412181, tol=1e-6)
    assert dist.cdf(1076091) == exact_tests.ppois(1076091, 1081711.834412181, approx=True)
    assert exact_tests.poisson.cdf(1076091, 1081711.834412181, approx=True) == exact_tests.ppois(1076091, 1081711.834412181, approx=True)


# Entries are of the form (k, mu, midp, logp, p_want), where p_want was
# computed with a 320-bit MPFR implementation.
poissontest_cases = [
    (1000123456, 1e9, False, False, 9.4638293091475092e-05),
    (999876544, 1e9, False, False, 9.4582633468018055e-05),
    (10, 2000.0, False, True, -1938.6127473836495),
    (137, 24.19893, False, False, 2.845227264114493e-56),
    (0, 3.7, False, False, 0.059964256942478408),
    (12, 5.0, False, False, 0.0054530919130093593),
    (3, 5.0, False, False, 0.50284245232442304),
    (7, 3.3, True, False, 0.035371368582143993),
    (30000, 31000.25, False, False, 1.1425590657444108e-08),
    ]

def test_poissontest():
    for k, mu, midp, logp, want in poissontest_cases:
        got = exact_tests.poissontest(k, 1.0, mu, midp=midp, logp=logp)
        assert got == pytest.approx(want, rel=1e-14, abs=0), (k, mu, midp, logp)
        # t and r are interchangeable
        assert exact_tests.poissontest(k, mu, 1.0, midp=midp, logp=logp) == got

    # One-sided tests, including midp.
    mu = 3.3
    pmfs = brute_pmfs(mu, 200)
    for k in range(15):
        less = math.fsum(pmfs[:k+1])
        greater = math.fsum(pmfs[k:])
        assert exact_tests.poissontest(k, r=mu, alternative="less") == pytest.approx(less, rel=1e-13)
        assert exact_tests.poissontest(k, r=mu, alternative="greater") == pytest.approx(greater, rel=1e-13)
        assert exact_tests.poissontest(k, r=mu, alternative="less", midp=True) == pytest.approx(less - 0.5 * pmfs[k], rel=1e-13)
        assert exact_tests.poissontest(k, r=mu, alternative="greater", midp=True) == pytest.approx(greater - 0.5 * pmfs[k], rel=1e-13)
        # Two-sided: sum of probabilities no greater than pmf(k).
        twosided = math.fsum(p for p in pmfs if p <= pmfs[k] * (1 + 1e-12))
        assert exact_tests.poissontest(k, r=mu) == pytest.approx(twosided, rel=1e-13)
    assert exact_tests.poissontest(0, r=mu, alternative="greater") == 1.0

    # Ties: with mu=6, pmf(5) == pmf(6), and pmf(k) == pmf(k') never happens
    # otherwise for this mu except between modes.
    assert exact_tests.poissontest(5, r=6.0) == 1.0
    assert exact_tests.poissontest(6, r=6.0) == 1.0
    pmfs = brute_pmfs(6.0, 200)
    assert exact_tests.poissontest(5, r=6.0, midp=True) == pytest.approx(1 - pmfs[5], rel=1e-13)
    # t=3, r=1/3: the float product is slightly less than 1, so x=0 is the
    # unique mode; Fraction(1, 3) makes the product exactly 1, so pmf(0) and
    # pmf(1) tie.
    assert exact_tests.poissontest(0, 3.0, 1/3, midp=True) == pytest.approx(1 - 0.5 * math.exp(-1), rel=1e-13)
    assert exact_tests.poissontest(0, 3.0, Fraction(1, 3), midp=True) == pytest.approx(1 - math.exp(-1), rel=1e-13)
    assert exact_tests.poissontest(0, t=0.0) == 1.0
    assert exact_tests.poissontest(3, t=0.0) == 0.0
    assert exact_tests.poissontest(0, t=0.0, midp=True) == 0.5

    # Broadcasting.
    ks = np.array([[0, 3], [12, 7]])
    got = exact_tests.poissontest(ks, 1.0, 5.0)
    assert got.shape == (2, 2)
    assert got[1, 0] == exact_tests.poissontest(12, 1.0, 5.0)
    got = exact_tests.poissontest(3, np.array([1.0, 2.0]), 2.5)
    assert got[1] == exact_tests.poissontest(3, 2.0, 2.5)

    with pytest.raises(RuntimeError):
        exact_tests.poissontest(1, r=2.0**51)
    with pytest.raises(RuntimeError):
        exact_tests.poissontest(-1, r=1.0)
    with pytest.raises(RuntimeError):
        exact_tests.poissontest(1, r=-1.0)
    with pytest.raises(RuntimeError):
        exact_tests.poissontest(1, alternative="two.sided")


def test_poissontest2():
    # R 4.6.1:
    #   poisson.test(c(11, 6+8+7), c(800, 1083+1050+878))$p.value
    assert exact_tests.poissontest2(11, 6+8+7, 800, 1083+1050+878) == pytest.approx(0.07966863303332952, rel=1e-12)
    # Reduction to binomtest(), with the success probability computed exactly.
    for alt in ("two-sided", "less", "greater"):
        for midp in (False, True):
            got = exact_tests.poissontest2(7, 13, 1.5, 2.25, 0.75, alternative=alt, midp=midp)
            want = exact_tests.binomtest(7, 20, Fraction(9, 8) / (Fraction(9, 8) + Fraction(9, 4)), alternative=alt, midp=midp)
            assert got == want
    got = exact_tests.poissontest2(np.array([3, 7]), 13, 1.5, 2.25, 0.75)
    assert got[1] == exact_tests.poissontest2(7, 13, 1.5, 2.25, 0.75)
    with pytest.raises(RuntimeError):
        exact_tests.poissontest2(1, 2, 0.0, 0.0)
