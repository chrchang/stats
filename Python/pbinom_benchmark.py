#!/usr/bin/env python3
import argparse
import exact_tests
import math
import scipy
import timeit
"""
This collects benchmark results for exact_tests.pbinom().
"""

def pbinom_benchmark(p: float, z: float, min_pow2: int, max_pow2: int, num_trials_per_pow2: int, approx=False, logp=False, bench_scipy=False):
    pq = p * (1.0 - p)
    if bench_scipy:
        warmup = scipy.stats.binom.cdf(1, 2, 0.5)
    for pow2 in range(min_pow2, max_pow2 + 1):
        n = 2**pow2 - 1
        stdev = math.sqrt(n * pq)
        k = round(n * p + z * stdev)
        if k < 0:
            k = 0
        elif k > n:
            k = n
        secs = 0.0
        if not bench_scipy:
            secs = timeit.timeit(lambda: exact_tests.pbinom(k, n, p, approx=approx, logp=logp), number=num_trials_per_pow2)
        else:
            if not logp:
                secs = timeit.timeit(lambda: scipy.stats.binom.cdf(k, n, p), number=num_trials_per_pow2)
            else:
                secs = timeit.timeit(lambda: scipy.stats.binom.logcdf(k, n, p), number=num_trials_per_pow2)
        print("n=(2^" + str(pow2) + ")-1: " + str(secs))



def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-p', '--succ-prob', type=float, default=0.5,
                             help="Binomial distribution success-probability to test.")
    optionalarg.add_argument('-z', '--z-score', type=float, default=0.0,
                             help="Success-count z-score to test.")
    optionalarg.add_argument('-f', '--from-pow2', type=int, default=15,
                             help="Start testing at n=2**<this value> - 1.")
    optionalarg.add_argument('-t', '--to-pow2', type=int, default=33,
                             help="Continue testing up to n=2**(<this value>+1) - 1.")
    optionalarg.add_argument('-n', '--number', type=int, default=25,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-a', '--approx', action="store_true",
                             help="Test approx=True.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True.")
    optionalarg.add_argument('-s', '--scipy', action="store_true",
                             help="Benchmark scipy.stats.binom instead of exact_tests.pbinom.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    pbinom_benchmark(cmd_args.succ_prob,
                     cmd_args.z_score,
                     cmd_args.from_pow2,
                     cmd_args.to_pow2,
                     cmd_args.number,
                     cmd_args.approx,
                     cmd_args.logp,
                     cmd_args.scipy)


if __name__ == '__main__':
    main()
