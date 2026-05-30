#!/usr/bin/env python3
import argparse
import exact_tests
import math
import random
import scipy
"""
This checks exact_tests.pbinom(approx=True)'s accuracy against that of
scipy.stats.binom.{,log}cdf().
"""

def pbinom_accuracy_test(p: float, z: float, min_pow2: int, max_pow2: int, num_trials_per_pow2: int, logp=False):
    pq = p * (1.0 - p)
    for pow2 in range(min_pow2, max_pow2 + 1):
        min_n = 2 ** pow2
        max_n = min_n * 2 - 1
        relerr_ssq = 0.0
        relerr_scipy_ssq = 0.0
        for trial_idx in range(0, num_trials_per_pow2):
            n = random.randint(min_n, max_n)
            stdev = math.sqrt(n * pq)
            k = round(n * p + z * stdev + 0.5)
            if k < 0:
                k = 0
            elif k > n:
                k = n
            want = exact_tests.pbinom(k, n, p, logp=logp)
            got = exact_tests.pbinom(k, n, p, approx=True, logp=logp)
            relerr = (got - want) / want
            relerr_ssq += relerr * relerr
            if logp:
                got = scipy.stats.binom.logcdf(k, n, p)
            else:
                got = scipy.stats.binom.cdf(k, n, p)
            relerr = (got - want) / want
            relerr_scipy_ssq += relerr * relerr
        approx_rms = math.sqrt(relerr_ssq / num_trials_per_pow2)
        scipy_rms = math.sqrt(relerr_scipy_ssq / num_trials_per_pow2)
        print("[2^" + str(pow2) + ", 2^" + str(pow2+1) + " - 1): approxRMS=" + str(approx_rms) + "  scipyRMS=" + str(scipy_rms))


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-p', '--succ-prob', type=float, default=0.5,
                             help="Binomial distribution success-probability to test.")
    optionalarg.add_argument('-z', '--z-score', type=float, default=0.0,
                             help="Success-count z-score to test.")
    optionalarg.add_argument('-f', '--from-pow2', type=int, default=9,
                             help="Start testing at n=2**<this value>.")
    optionalarg.add_argument('-t', '--to-pow2', type=int, default=33,
                             help="Continue testing up to n=2**(<this value>+1) - 1.")
    optionalarg.add_argument('-n', '--number', type=int, default=1000,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True.")
    optionalarg.add_argument('-s', '--seed', type=int, default=1,
                             help="RNG seed.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    random.seed(cmd_args.seed)
    pbinom_accuracy_test(cmd_args.succ_prob,
                         cmd_args.z_score,
                         cmd_args.from_pow2,
                         cmd_args.to_pow2,
                         cmd_args.number,
                         cmd_args.logp)


if __name__ == '__main__':
    main()
