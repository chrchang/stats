#!/usr/bin/env python3
import exact_tests
import math
import pytest
import sys

DBL_MIN = sys.float_info.min

# These test cases were scraped from
#   scipy/stats/tests/data/fisher_exact_results_from_r.py
# which were in turn gathered from R 3.6.2.
# Entries are of the form (table, alternative, p_want, odds_ratio, ci95_low,
#                          ci95_high, ci99_low, ci99_high).
# Note that R odds-ratio numbers have low accuracy.
scipy_fisher22_cases = [
    ([[100, 2], [1000, 5]], "two-sided", 0.1300759363430016, 0.25055839934223, 0.04035202926536294, 2.662846672960251, 0.02502345007115455, 6.304424772117853),
    ([[2, 7], [8, 2]], "two-sided", 0.02301413756522116, 0.0858623513573622, 0.004668988338943325, 0.895792956493601, 0.001923034001462487, 1.53670836950172),
    ([[5, 1], [10, 10]], "two-sided", 0.1973244147157191, 4.725646047336587, 0.4153910882532168, 259.2593661129417, 0.2397970951413721, 1291.342011095509),
    ([[5, 15], [20, 20]], "two-sided", 0.09580440012477633, 0.3394396617440851, 0.08056337526385809, 1.22704788545557, 0.05127576113762925, 1.717176678806983),
    ([[5, 16], [16, 25]], "two-sided", 0.2697004098849359, 0.4937791394540491, 0.1176691231650079, 1.787463657995973, 0.07498546954483619, 2.506969905199901),
    ([[10, 5], [10, 1]], "two-sided", 0.1973244147157192, 0.2116112781158479, 0.003857141267422399, 2.407369893767229, 0.0007743881879531337, 4.170192301163831),
    ([[10, 5], [10, 0]], "two-sided", 0.06126482213438735, 0, 0, 1.451643573543705, 0, 2.642491011905582),
    ([[5, 0], [1, 4]], "two-sided", 0.04761904761904762, math.inf, 1.024822256141754, math.inf, 0.496935393325443, math.inf),
    ([[0, 5], [1, 4]], "two-sided", 1, 0, 0, 39.00054996869288, 0, 198.019801980198),
    ([[5, 1], [0, 4]], "two-sided", 0.04761904761904761, math.inf, 1.024822256141754, math.inf, 0.496935393325443, math.inf),
    ([[0, 1], [3, 2]], "two-sided", 1, 0, 0, 39.00054996869287, 0, 198.019801980198),
    ([[200, 7], [8, 300]], "two-sided", 2.005657880389071e-122, 977.7866978606228, 349.2595113327733, 3630.382605689872, 270.0334165523604, 5461.333333326708),
    ([[28, 21], [6, 1957]], "two-sided", 5.728437460831947e-44, 425.2403028434684, 152.4166024390096, 1425.700792178893, 116.7944750275836, 1931.995993191814),
    ([[190, 800], [200, 900]], "two-sided", 0.574111858126088, 1.068697577856801, 0.8520462587912048, 1.340148950273938, 0.7949398282935892, 1.436229679394333),
    ([[100, 2], [1000, 5]], "less", 0.1300759363430016, 0.25055839934223, 0, 1.797867027270803, 0, 4.375946050832565),
    ([[2, 7], [8, 2]], "less", 0.0185217259520665, 0.0858623513573622, 0, 0.6785254803404526, 0, 1.235282118191202),
    ([[5, 1], [10, 10]], "less", 0.9782608695652173, 4.725646047336587, 0, 127.8497388102893, 0, 657.2063583945989),
    ([[5, 15], [20, 20]], "less", 0.05625775074399956, 0.3394396617440851, 0, 1.032332939718425, 0, 1.498867660683128),
    ([[5, 16], [16, 25]], "less", 0.1808979350599346, 0.4937791394540491, 0, 1.502407513296985, 0, 2.186159386716762),
    ([[10, 5], [10, 1]], "less", 0.1652173913043479, 0.2116112781158479, 0, 1.820421051562392, 0, 3.335351451901569),
    ([[10, 5], [10, 0]], "less", 0.0565217391304348, 0, 0, 1.06224603077045, 2.075407697450433),
    ([[5, 0], [1, 4]], "less", 1, math.inf, 0, math.inf, 0, math.inf),
    ([[0, 5], [1, 4]], "less", 0.5, 0, 0, 19.00192394479939, 0, 99.00009507969122),
    ([[5, 1], [0, 4]], "less", 1, math.inf, 0, math.inf, 0, math.inf),
    ([[0, 1], [3, 2]], "less", 0.4999999999999999, 0, 0, 19.00192394479939, 0, 99.00009507969123),
    ([[200, 7], [8, 300]], "less", 1, 977.7866978606228, 0, 3045.460216525746, 0, 4503.078257659934),
    ([[28, 21], [6, 1957]], "less", 1, 425.2403028434684, 0, 1186.440170942579, 0, 1811.766127544222),
    ([[190, 800], [200, 900]], "less", 0.7416227010368963, 1.068697577856801, 0, 1.293551891610822, 0, 1.396522811516685),
    ([[100, 2], [1000, 5]], "greater", 0.979790445314723, 0.25055839934223, 0.05119649909830196, math.inf, 0.03045407081240429, math.inf),
    ([[2, 7], [8, 2]], "greater", 0.9990149169715733, 0.0858623513573622, 0.007163749169069961, math.inf, 0.002768053063547901, math.inf),
    ([[5, 1], [10, 10]], "greater", 0.1652173913043478, 4.725646047336587, 0.5493234651081089, math.inf, 0.2998184792279909, math.inf),
    ([[5, 15], [20, 20]], "greater", 0.9849086665340765, 0.3394396617440851, 0.1003538933958604, math.inf, 0.06180414342643172, math.inf),
    ([[5, 16], [16, 25]], "greater", 0.9330176609214881, 0.4937791394540491, 0.146507416280863, math.inf, 0.09037094010066403, math.inf),
    ([[10, 5], [10, 1]], "greater", 0.9782608695652174, 0.2116112781158479, 0.007821681994077808, math.inf, 0.001521592095430679, math.inf),
    ([[10, 5], [10, 0]], "greater", 1, 0, 0, math.inf, 0, math.inf),
    ([[5, 0], [1, 4]], "greater", 0.02380952380952382, math.inf, 1.487678929918272, math.inf, 0.6661157890359722, math.inf),
    ([[0, 5], [1, 4]], "greater", 1, 0, 0, math.inf, 0, math.inf),
    ([[5, 1], [0, 4]], "greater", 0.0238095238095238, math.inf, 1.487678929918272, math.inf, 0.6661157890359725, math.inf),
    ([[0, 1], [3, 2]], "greater", 1, 0, 0, math.inf, 0, math.inf),
    ([[200, 7], [8, 300]], "greater", 2.005657880388915e-122, 977.7866978606228, 397.784359748113, math.inf, 297.9619252357688, math.inf),
    ([[28, 21], [6, 1957]], "greater", 5.728437460831983e-44, 425.2403028434684, 174.7148056880929, math.inf, 130.3213490295859, math.inf),
    ([[190, 800], [200, 900]], "greater", 0.2959825901308897, 1.068697577856801, 0.8828406663967776, math.inf, 0.8176272148267533, math.inf),
    ]


def test_fisher():
    # 1e-13 is several hundred ULPs; preexisting good implementations should
    # agree to that relative tolerance for small cases (such as the scipy cases
    # where no table sum > 3000).
    for test_case in scipy_fisher22_cases:
        pval = exact_tests.fisher_exact(test_case[0], alternative=test_case[1])
        assert pval == pytest.approx(test_case[2], rel=1e-13, abs=DBL_MIN), str(test_case)
    pval = exact_tests.fisher_exact([[4, 0], [0, 4]], midp=True)
    logp = exact_tests.fisher_exact([[4, 0], [0, 4]], midp=True, logp=True)
    assert pval == pytest.approx(1/70, rel=1e-13, abs=0), "midp"
    assert logp == pytest.approx(math.log(1/70), rel=1e-13, abs=0), "logp"
    assert exact_tests.fisher_exact([[1e9, 2e9], [3e9, 5.999e9]], alternative="less") == pytest.approx(9.686381892010733e-05, rel=1e-10, abs=0)
    assert exact_tests.fisher_exact([[2e12, 3e12], [4e12, 5.9999e12]]) == pytest.approx(2.95770177417893e-50, rel=1e-8, abs=0)
    # todo: test exception-throwing cases


def test_cond_odds_ratio():
    for test_case in scipy_fisher22_cases:
        if len(test_case) < 8:
            continue
        m11 = test_case[0][0][0]
        m12 = test_case[0][0][1]
        m21 = test_case[0][1][0]
        m22 = test_case[0][1][1]
        cond_odds_ratio = exact_tests.cond_odds_ratio(m11, m12, m21, m22)
        # Imitate scipy test_odds_ratio.py.
        or_rtol = 5e-4
        ci_rtol = 2e-2
        if cond_odds_ratio >= 400:
            or_rtol = 5e-2
            ci_rtol = 1e-1
        ci95_low = 0.025
        ci95_high = 0.975
        ci99_low = 0.005
        ci99_high = 0.995
        if test_case[1] == "less":
            ci95_low = 0
            ci95_high = 0.95
            ci99_low = 0
            ci99_high = 0.99
        elif test_case[1] == "greater":
            ci95_low = 0.05
            ci95_high = 1
            ci99_low = 0.01
            ci99_high = 1
        ci95 = exact_tests.cond_odds_ratio_ci(m11, m12, m21, m22, ci95_low, ci95_high)
        ci99 = exact_tests.cond_odds_ratio_ci(m11, m12, m21, m22, ci99_low, ci99_high)
        assert cond_odds_ratio == pytest.approx(test_case[3], rel=or_rtol, abs=DBL_MIN), str(test_case)
        assert ci95[0] == pytest.approx(test_case[4], rel=ci_rtol, abs=DBL_MIN), str(test_case)
        assert ci95[1] == pytest.approx(test_case[5], rel=ci_rtol, abs=DBL_MIN), str(test_case)
        assert ci99[0] == pytest.approx(test_case[6], rel=ci_rtol, abs=DBL_MIN), str(test_case)
        assert ci99[1] == pytest.approx(test_case[7], rel=ci_rtol, abs=DBL_MIN), str(test_case)
    # these values from scipy.stats.contingency.odds_ratio()
    assert exact_tests.cond_odds_ratio(1e4, 2e4, 3e4, 5e4) == pytest.approx(0.8333346994547735, rel=1e-6, abs=0)
    ci95 = exact_tests.cond_odds_ratio_ci(1e4, 2e4, 3e4, 5e4)
    assert ci95[0] == pytest.approx(0.8102717994162574, rel=1e-6, abs=0)
    assert ci95[1] == pytest.approx(0.857026539321631, rel=1e-6, abs=0)
    # these must not enter infinite loop
    assert exact_tests.cond_odds_ratio(2**51 - 2, 1, 1, 2**51 - 1) == pytest.approx(3.036688e+30, rel=1e-6, abs=0)
    assert exact_tests.cond_odds_ratio(1, 2**51 - 2, 2**51 - 1, 1) == pytest.approx(3.293061e-31, rel=1e-6, abs=0)
