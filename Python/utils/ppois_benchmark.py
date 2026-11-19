#!/usr/bin/env python3
import argparse
import exact_tests
import math
import scipy
import timeit
"""
This collects benchmark results for exact_tests.ppois().
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


def ppois_benchmark(z: float, pow2s: list[int], num_trials_per_pow2: int, logp=False):
    warmup = exact_tests.ppois(1, 2.5)
    warmup = scipy.stats.poisson.cdf(1, 2.5)
    for pow2 in pow2s:
        mu = 2**pow2 - 0.5
        k = round(mu + z * math.sqrt(mu))
        if k < 0:
            k = 0
        secs_noapprox = timeit.timeit(lambda: exact_tests.ppois(k, mu, logp=logp), number=num_trials_per_pow2) / num_trials_per_pow2
        secs_approx = timeit.timeit(lambda: exact_tests.ppois(k, mu, approx=True, logp=logp), number=num_trials_per_pow2) / num_trials_per_pow2
        if not logp:
            secs_scipy = timeit.timeit(lambda: scipy.stats.poisson.cdf(k, mu), number=num_trials_per_pow2) / num_trials_per_pow2
        else:
            secs_scipy = timeit.timeit(lambda: scipy.stats.poisson.logcdf(k, mu), number=num_trials_per_pow2) / num_trials_per_pow2
        print(f"mu=(2^{pow2})-0.5: base={secs_noapprox:.3g}  approx={secs_approx:.3g}  scipy={secs_scipy:.3g} sec/iter")


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-z', '--z-score', type=float, default=0.0,
                             help="Event-count z-score to test.")
    optionalarg.add_argument('-e', '--exps', type=str, default="5,20,35",
                             help="Test mu=2**<these values> - 0.5.")
    optionalarg.add_argument('-n', '--number', type=int, default=10,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    ppois_benchmark(cmd_args.z_score,
                    parse_range_string(cmd_args.exps),
                    cmd_args.number,
                    cmd_args.logp)


if __name__ == '__main__':
    main()
