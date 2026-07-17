#!/usr/bin/env python3
import argparse
import exact_tests
import math
import numpy as np
import scipy
import timeit
"""
This collects benchmark results for exact_tests.binom.pmf() on array-like
input.
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


def binom_pmf_benchmark(p: float, z: float, size: int, pow2s: list[int], num_trials_per_pow2: int, logp=False):
    warmup = exact_tests.binom.pmf(1, 2, [0.25, 0.5, 0.75])
    warmup = scipy.stats.binom.pmf(1, 2, [0.25, 0.5, 0.75])
    pq = p * (1.0 - p)
    odds = p / (1.0 - p)
    secs_main = 0.0
    secs_scipy = 0.0
    for pow2 in pow2s:
        n = 2**pow2 - 1
        stdev = math.sqrt(n * pq)
        k = round(n * p + z * stdev)
        if k < 0:
            k = 0
        elif k > n:
            k = n
        # assumes pow2 > 0
        ps = [o/(o+1) for o in np.linspace(odds / (1+1/n), odds * (1+1/n), num=size)]
        if not logp:
            secs_main = timeit.timeit(lambda: exact_tests.binom.pmf(k, n, ps), number=num_trials_per_pow2) / num_trials_per_pow2
            secs_scipy = timeit.timeit(lambda: scipy.stats.binom.pmf(k, n, ps), number=num_trials_per_pow2) / num_trials_per_pow2
        else:
            secs_main = timeit.timeit(lambda: exact_tests.binom.logpmf(k, n, ps), number=num_trials_per_pow2) / num_trials_per_pow2
            secs_scipy = timeit.timeit(lambda: scipy.stats.binom.logpmf(k, n, ps), number=num_trials_per_pow2) / num_trials_per_pow2
        print(f"n=(2^{pow2})-1: base={secs_main:.3g}  scipy={secs_scipy:.3g} sec/iter")


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-p', '--succ-prob', type=float, default=0.5,
                             help="Central binomial distribution success-probability to test.")
    optionalarg.add_argument('-z', '--z-score', type=float, default=0.0,
                             help="Success-count z-score to test.")
    optionalarg.add_argument('-s', '--size', type=int, default=10,
                             help="Size of p-vector to test.")
    optionalarg.add_argument('-e', '--exps', type=str, default="5,20,35",
                             help="Test n=2**<these values> - 1.")
    optionalarg.add_argument('-n', '--number', type=int, default=10,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    binom_pmf_benchmark(cmd_args.succ_prob,
                        cmd_args.z_score,
                        cmd_args.size,
                        parse_range_string(cmd_args.exps),
                        cmd_args.number,
                        cmd_args.logp)


if __name__ == '__main__':
    main()
