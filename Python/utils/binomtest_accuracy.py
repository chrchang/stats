#!/usr/bin/env python3
import argparse
import exact_tests
import mpfr_impls
import math
import random
import scipy
import sys
"""
Primary mode checks exact_tests.pbinom(approx=False) and
scipy.stats.binom.cdf()'s accuracy against a simple MPFR-based implementation
(slow!).

"--bits 0" mode checks exact_tests.pbinom(approx=True) and
scipy.stats.binom.{,log}cdf(), treating exact_tests.pbinom(approx=False) as
ground truth.
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


def binomtest_accuracy_test(p: float, z: float, pow2s: list[int], num_trials_per_pow2: int, bits: int, logp: bool):
    n = 0
    pq = p * (1.0 - p)
    want = 0.0
    got_scipy = 0.0
    ns = []
    for pow2 in pow2s:
        min_n = 2 ** pow2
        n_limit = min_n * 2
        relerr_base_ssq = 0.0
        relerr_scipy_ssq = 0.0
        if num_trials_per_pow2 >= n_limit - min_n:
            ns = range(min_n, n_limit)
        else:
            ns = random.sample(range(min_n, n_limit), k=num_trials_per_pow2)
        for n in ns:
            stdev = math.sqrt(n * pq)
            k = round(n * p + z * stdev + 0.5)
            if k < 0:
                k = 0
            elif k > n:
                k = n
            got_base = flush_if_denormal(exact_tests.binomtest(k, n, p, logp=logp))
            want = flush_if_denormal(mpfr_impls.binomtest(k, n, p, bits=bits, return_log=logp))
            relerr = compute_relerr(got_base, want)
            relerr_base_ssq += relerr * relerr
            if not logp:
                got_scipy = scipy.stats.binomtest(k, n, p).pvalue
                relerr = compute_relerr(got_scipy, want)
                relerr_scipy_ssq += relerr * relerr
        num_trials = len(ns)
        base_rms = math.sqrt(relerr_base_ssq / num_trials)
        if not logp:
            scipy_rms = math.sqrt(relerr_scipy_ssq / num_trials)
            print(f"n in [2^{pow2}, 2^{pow2+1}): errRMS={base_rms:.6g}  scipyErrRMS={scipy_rms:.6g}")
        else:
            print(f"n in [2^{pow2}, 2^{pow2+1}): errRMS={base_rms:.6g}")


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-p', '--succ-prob', type=float, default=0.5,
                             help="Binomial distribution success-probability to test.")
    optionalarg.add_argument('-z', '--z-score', type=float, default=0.0,
                             help="Success-count z-score to test.")
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
    binomtest_accuracy_test(cmd_args.succ_prob,
                            cmd_args.z_score,
                            parse_range_string(cmd_args.exps),
                            cmd_args.number,
                            cmd_args.bits,
                            cmd_args.logp)


if __name__ == '__main__':
    main()
