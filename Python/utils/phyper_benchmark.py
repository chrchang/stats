#!/usr/bin/env python3
import argparse
import exact_tests
import math
import scipy
import timeit
"""
This collects benchmark results for exact_tests.pbinom().
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


def phyper_benchmark(ab: float, ac: float, z: float, pow2s: list[int], num_trials_per_pow2: int, logp=False, omit_scipy=False):
    if not omit_scipy:
        warmup = scipy.stats.hypergeom.cdf(1, 4, 2, 2)
    secs_scipy = 0.0
    for pow2 in pow2s:
        n = 2**pow2 - 1
        cur_ab = round(ab * n)
        cur_ac = round(ac * n)
        stdev = 0
        if n > 1:
            variance = cur_ab * cur_ac * (n - cur_ab) * (n - cur_ac) / (n * n * (n - 1))
            stdev = math.sqrt(variance)
        k = round(cur_ab * cur_ac / n)
        secs_noapprox = timeit.timeit(lambda: exact_tests.phyper(k, cur_ac, n - cur_ac, cur_ab, logp=logp), number=num_trials_per_pow2) / num_trials_per_pow2
        secs_approx = timeit.timeit(lambda: exact_tests.phyper(k, cur_ac, n - cur_ac, cur_ab, approx=True, logp=logp), number=num_trials_per_pow2) / num_trials_per_pow2
        if not omit_scipy:
            if not logp:
                secs_scipy = timeit.timeit(lambda: scipy.stats.hypergeom.cdf(k, n, cur_ab, cur_ac), number=num_trials_per_pow2) / num_trials_per_pow2
            else:
                secs_scipy = timeit.timeit(lambda: scipy.stats.hypergeom.logcdf(k, n, cur_ab, cur_ac), number=num_trials_per_pow2) / num_trials_per_pow2
        print_str = "n=(2^" + str(pow2) + ")-1: base=" + str(secs_noapprox) + "  approx=" + str(secs_approx)
        if not omit_scipy:
            print_str += "  scipy=" + str(secs_scipy)
        print(print_str + " sec/iter")


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-r', '--row-prop', type=float, default=0.5,
                             help="Proportion of total in first row.")
    optionalarg.add_argument('-c', '--col-prop', type=float, default=0.5,
                             help="Proportion of total in first column.")
    optionalarg.add_argument('-z', '--z-score', type=float, default=0.0,
                             help="Z-score to test.")
    optionalarg.add_argument('-e', '--exps', type=str, default="5,20,35",
                             help="Test n=2**<these values> - 1.")
    optionalarg.add_argument('-n', '--number', type=int, default=3,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True.")
    optionalarg.add_argument('-o', '--omit-scipy', action="store_true",
                             help="Skip scipy.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    phyper_benchmark(cmd_args.row_prop,
                     cmd_args.col_prop,
                     cmd_args.z_score,
                     parse_range_string(cmd_args.exps),
                     cmd_args.number,
                     cmd_args.logp,
                     cmd_args.omit_scipy)


if __name__ == '__main__':
    main()
