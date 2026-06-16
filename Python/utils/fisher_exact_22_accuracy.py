#!/usr/bin/env python3
import argparse
import exact_tests
import mpfr_impls
import math
import random
import scipy
import sys
"""
Checks exact_tests.fisher_exact() and (when logp=False and pow2 < 31)
scipy.stats.fisher_exact()'s accuracy against a MPFR-based implementation
(slow!).
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


def fisher_exact_22_accuracy_test(ab: float, ac: float, z: float, pow2s: list[int], num_trials_per_pow2: int, bits: int, logp: bool):
    ns = []
    for pow2 in pow2s:
        min_n = 2 ** pow2
        n_limit = min_n * 2
        if num_trials_per_pow2 >= n_limit - min_n:
            ns = range(min_n, n_limit)
        else:
            ns = random.sample(range(min_n, n_limit), k=num_trials_per_pow2)
        relerr_base_ssq = 0.0
        relerr_scipy_ssq = 0.0
        include_scipy = (pow2 < 31) and not logp
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
            want = flush_if_denormal(mpfr_impls.fisher_exact_22(a, b, c, d, bits=bits, return_log=logp))
            got_base = flush_if_denormal(exact_tests.fisher_exact(table, logp=logp))
            relerr = compute_relerr(got_base, want)
            relerr_base_ssq += relerr * relerr
            if include_scipy:
                got_scipy = scipy.stats.fisher_exact(table).pvalue
                relerr = compute_relerr(got_scipy, want)
                relerr_scipy_ssq += relerr * relerr
        num_trials = len(ns)
        base_rms = math.sqrt(relerr_base_ssq / num_trials)
        print_str = f"n in [2^{pow2}, 2^{pow2+1}): errRMS={base_rms:.3g}"
        if include_scipy:
            scipy_rms = math.sqrt(relerr_scipy_ssq / num_trials)
            print_str += f"  scipyErrRMS={scipy_rms:.3g}"
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
    optionalarg.add_argument('-e', '--exps', type=str, default="5,20,35",
                             help="Test n~=2**<these values>.")
    optionalarg.add_argument('-n', '--number', type=int, default=10,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-b', '--bits', type=int, default=256,
                             help="MPFR precision.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True (disables scipy).")
    optionalarg.add_argument('-s', '--seed', type=int, default=1,
                             help="RNG seed.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    random.seed(cmd_args.seed)
    fisher_exact_22_accuracy_test(cmd_args.row_prop,
                                  cmd_args.col_prop,
                                  cmd_args.z_score,
                                  parse_range_string(cmd_args.exps),
                                  cmd_args.number,
                                  cmd_args.bits,
                                  cmd_args.logp)


if __name__ == '__main__':
    main()
