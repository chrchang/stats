#!/usr/bin/env python3
import argparse
import exact_tests
import mpfr_impls
import math
import random
import snphwe
import sys
"""
Checks exact_tests.HWE_exact() and (when logp=False) snphwe.snphwe()'s accuracy
against a MPFR-based implementation (slow!).

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


def HWE_accuracy_test(maf: float, z: float, pow2s: list[int], num_trials_per_pow2: int, bits: int, logp: bool):
    ns = []
    for pow2 in pow2s:
        min_n = 2 ** pow2
        n_limit = min_n * 2
        if num_trials_per_pow2 >= n_limit - min_n:
            ns = range(min_n, n_limit)
        else:
            ns = random.sample(range(min_n, n_limit), k=num_trials_per_pow2)
        relerr_base_ssq = 0.0
        relerr_snphwe_ssq = 0.0
        for n in ns:
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

            want = flush_if_denormal(mpfr_impls.snphwe(hets, homr, homc, bits=bits, return_log=logp))
            got_base = flush_if_denormal(exact_tests.HWE_exact(homr, hets, homc, logp=logp))
            relerr = compute_relerr(got_base, want)
            relerr_base_ssq += relerr * relerr
            if not logp:
                got_snphwe = snphwe.snphwe(hets, homr, homc)
                relerr = compute_relerr(got_snphwe, want)
                relerr_snphwe_ssq += relerr * relerr
        num_trials = len(ns)
        base_rms = math.sqrt(relerr_base_ssq / num_trials)
        print_str = f"n in [2^{pow2}, 2^{pow2+1}): errRMS={base_rms:.3g}"
        if not logp:
            snphwe_rms = math.sqrt(relerr_snphwe_ssq / num_trials)
            print_str += f"  snphweErrRMS={snphwe_rms:.3g}"
        print(print_str)


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    requiredarg = parser.add_argument_group('Required Arguments')
    requiredarg.add_argument('-z', '--z-score', type=float, required=True,
                             help="Z-score to test.")
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-m', '--maf', type=float, default=0.5,
                             help="Minor allele frequency.")
    optionalarg.add_argument('-e', '--exps', type=str, default="10,15,20",
                             help="Test n~=2**<these values>.")
    optionalarg.add_argument('-n', '--number', type=int, default=10,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-b', '--bits', type=int, default=256,
                             help="MPFR precision.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True (disables snphwe).")
    optionalarg.add_argument('-s', '--seed', type=int, default=1,
                             help="RNG seed.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    random.seed(cmd_args.seed)
    HWE_accuracy_test(cmd_args.maf,
                      cmd_args.z_score,
                      parse_range_string(cmd_args.exps),
                      cmd_args.number,
                      cmd_args.bits,
                      cmd_args.logp)


if __name__ == '__main__':
    main()
