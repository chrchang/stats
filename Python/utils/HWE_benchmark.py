#!/usr/bin/env python3
import argparse
import exact_tests
import math
import snphwe
import timeit
"""
This collects benchmark results for exact_tests.HWE_exact().

Note that you may need to downgrade setuptools (v80.10.2 is ok) to get snphwe
to work.
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


def HWE_exact_benchmark(maf: float, z: float, pow2s: list[int], num_trials_per_pow2: int):
    warmup = exact_tests.HWE_exact(2, 2, 2)
    warmup = snphwe.snphwe(2, 2, 2)
    for pow2 in pow2s:
        n = 2**pow2 - 1
        n_allele = 2 * n
        n_rare = round(n_allele * maf)
        n_common = n_allele - n_rare
        mean = n_rare * n_common / (n_allele - 1)
        stdev = 0
        if n > 1:
            variance = mean * (1 + (n_rare - 1) * (n_common - 1) / (n_allele - 3) - (n_rare * n_common) / (n_allele - 1))
            stdev = math.sqrt(variance)
        hets_raw = mean + z * stdev
        # Round to nearest integer with correct parity.
        parity = n_rare % 2
        hets = 2 * round((hets_raw - parity) * 0.5) + parity
        if hets < 0:
            hets = parity
        homr = (n_rare - hets) // 2
        if homr < 0:
            homr = 0
            hets = n_rare
        homc = (n_common - hets) // 2

        secs_base = timeit.timeit(lambda: exact_tests.HWE_exact(homr, hets, homc), number=num_trials_per_pow2) / num_trials_per_pow2
        secs_snphwe = timeit.timeit(lambda: snphwe.snphwe(hets, homr, homc), number=num_trials_per_pow2) / num_trials_per_pow2
        print(f"n=(2^{pow2})-1: base={secs_base:.3g}  snphwe={secs_snphwe:.3g} sec/iter")


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    requiredarg = parser.add_argument_group('Required Arguments')
    # z=0 is not representative.
    requiredarg.add_argument('-z', '--z-score', type=float, required=True,
                             help="Z-score to test.")
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-m', '--maf', type=float, default=0.5,
                             help="Minor allele frequency.")
    optionalarg.add_argument('-e', '--exps', type=str, default="10,15,20",
                             help="Test n=2**<these values> - 1.")
    optionalarg.add_argument('-n', '--number', type=int, default=25,
                             help="Number of trials per power-of-2 tier.")
    # scipy fisher_exact does not support logp.
    # optionalarg.add_argument('-l', '--logp', action="store_true",
    #                          help="Test logp=True.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    HWE_exact_benchmark(cmd_args.maf,
                        cmd_args.z_score,
                        parse_range_string(cmd_args.exps),
                        cmd_args.number)


if __name__ == '__main__':
    main()
