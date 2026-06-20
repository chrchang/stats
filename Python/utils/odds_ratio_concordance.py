#!/usr/bin/env python3
import argparse
import exact_tests
import math
import random
import scipy
import sys
"""
Checks exact_tests.cond_odds_ratio() and .cond_odds_ratio_ci() for approximate
concordance with scipy.stats.contingency.odds_ratio().
"""


def parse_range_string(input_str: str):
    result = []
    for part in input_str.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            result.extend(range(start, end + 1))
        else:
            result.append(int(part))
    return result


def flush_if_denormal(x: float):
    if abs(x) < sys.float_info.min:
        return 0.0
    return x


def compute_relerr(got: float, want: float):
    if want == 0.0:
        if got == 0.0:
            return 0.0
        return math.inf
    elif math.isnan(want):
        if math.isnan(got):
            return 0.0
        return math.inf
    return (got - want) / want


def odds_ratio_concordance_test(ab: float, ac: float, z: float, pow2s: list[int], num_trials_per_pow2: int, omit_ci: bool):
    ns = []
    for pow2 in pow2s:
        if pow2 >= 31:
            # scipy doesn't support this
            continue
        min_n = 2 ** pow2
        n_limit = min_n * 2
        if num_trials_per_pow2 >= n_limit - min_n:
            ns = range(min_n, n_limit)
        else:
            ns = random.sample(range(min_n, n_limit), k=num_trials_per_pow2)
        relerr_est_ssq = 0.0
        relerr_ci95_low_ssq = 0.0
        relerr_ci95_high_ssq = 0.0
        for n in ns:
            cur_ab = round(ab * n)
            cur_ac = round(ac * n)
            d_minus_a = n - cur_ab - cur_ac
            min_a = max(0, -d_minus_a)
            stdev = 0
            if n > 1:
                variance = cur_ab * cur_ac * (n - cur_ab) * (n - cur_ac) / (n * n * (n - 1))
                stdev = math.sqrt(variance)
            a = round((cur_ab * cur_ac / n) + z * stdev)
            if a < min_a:
                a = min_a
            elif a > min(cur_ab, cur_ac):
                a = min(cur_ab, cur_ac)
            b = cur_ab - a
            c = cur_ac - a
            d = d_minus_a + a
            table = [[a, b], [c, d]]
            # assume for now that scipy values are accurate
            want = scipy.stats.contingency.odds_ratio(table).statistic
            got = exact_tests.cond_odds_ratio(a, b, c, d)
            relerr = compute_relerr(got, want)
            relerr_est_ssq += relerr * relerr
            if not omit_ci:
                want_ci95 = scipy.stats.contingency.odds_ratio(table).confidence_interval()
                got_ci95 = exact_tests.cond_odds_ratio_ci(a, b, c, d)
                relerr_ci95_low = compute_relerr(got_ci95[0], want_ci95.low)
                relerr_ci95_high = compute_relerr(got_ci95[1], want_ci95.high)
                relerr_ci95_low_ssq = relerr_ci95_low * relerr_ci95_low
                relerr_ci95_high_ssq = relerr_ci95_high * relerr_ci95_high
        num_trials = len(ns)
        est_rms = math.sqrt(relerr_est_ssq / num_trials)
        print_str = f"n in [2^{pow2}, 2^{pow2+1}): estErrRMS={est_rms:.3g}"
        if not omit_ci:
            ci95_low_rms = math.sqrt(relerr_ci95_low_ssq / num_trials)
            ci95_high_rms = math.sqrt(relerr_ci95_high_ssq / num_trials)
            print_str += f"  ciLowErrRMS={ci95_low_rms:.3g}  ciHighErrRMS={ci95_high_rms:.3g}"
        print(print_str)


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    requiredarg = parser.add_argument_group('Required Arguments')
    requiredarg.add_argument('-z', '--z-score', type=float, required=True,
                             help="Z-score to test.")
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-r', '--row-prop', type=float, default=0.5,
                             help="Proportion of total in first row.")
    optionalarg.add_argument('-c', '--col-prop', type=float, default=0.5,
                             help="Proportion of total in first column.")
    optionalarg.add_argument('-e', '--exps', type=str, default="5,10,15",
                             help="Test n~=2**<these values>.")
    optionalarg.add_argument('-n', '--number', type=int, default=10,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-s', '--seed', type=int, default=1,
                             help="RNG seed.")
    optionalarg.add_argument('-o', '--omit-ci', action="store_true",
                             help="Skip CI.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    random.seed(cmd_args.seed)
    odds_ratio_concordance_test(cmd_args.row_prop,
                                cmd_args.col_prop,
                                cmd_args.z_score,
                                parse_range_string(cmd_args.exps),
                                cmd_args.number,
                                cmd_args.omit_ci)


if __name__ == '__main__':
    main()
