#!/usr/bin/env python3
import argparse
import exact_tests
import math
import scipy
import timeit
"""
This collects benchmark results for exact_tests.fisher_exact().
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


def fisher_exact_22_benchmark(ab: float, ac: float, z: float, pow2s: list[int], num_trials_per_pow2: int):
    warmup = exact_tests.fisher_exact([[1, 1], [1, 1]])
    warmup = scipy.stats.fisher_exact([[1, 1], [1, 1]])
    for pow2 in pow2s:
        n = 2**pow2 - 1
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
        table = [[a, cur_ac - a], [cur_ab - a, d_minus_a + a]]
        secs_base = timeit.timeit(lambda: exact_tests.fisher_exact(table), number=num_trials_per_pow2) / num_trials_per_pow2
        if pow2 < 32:
            # scipy implementation uses int32s
            secs_scipy = timeit.timeit(lambda: scipy.stats.fisher_exact(table), number=num_trials_per_pow2) / num_trials_per_pow2
            print(f"n=(2^{pow2})-1: base={secs_base:.3g}  scipy={secs_scipy:.3g} sec/iter")
        else:
            print(f"n=(2^{pow2})-1: base={secs_base:.3g} sec/iter")


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    requiredarg = parser.add_argument_group('Required Arguments')
    # z=0 is not representative.
    requiredarg.add_argument('-z', '--z-score', type=float, required=True,
                             help="Z-score to test.")
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-r', '--row-prop', type=float, default=0.5,
                             help="Proportion of total in first row.")
    optionalarg.add_argument('-p', '--col-prop', type=float, default=0.5,
                             help="Proportion of total in first column.")
    optionalarg.add_argument('-e', '--exps', type=str, default="5,20,35",
                             help="Test n=2**<these values> - 1.")
    optionalarg.add_argument('-n', '--number', type=int, default=10,
                             help="Number of trials per power-of-2 tier.")
    # scipy fisher_exact does not support logp.
    # optionalarg.add_argument('-l', '--logp', action="store_true",
    #                          help="Test logp=True.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    fisher_exact_22_benchmark(cmd_args.row_prop,
                              cmd_args.col_prop,
                              cmd_args.z_score,
                              parse_range_string(cmd_args.exps),
                              cmd_args.number)


if __name__ == '__main__':
    main()
