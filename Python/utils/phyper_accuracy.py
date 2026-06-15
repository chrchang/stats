#!/usr/bin/env python3
import argparse
import exact_tests
import mpfr_impls
import math
import random
import scipy
import sys
"""
Primary mode checks exact_tests.phyper(approx=False) and
scipy.stats.hypergeom.cdf()'s accuracy against a simple MPFR-based
implementation (slow!).

"--bits 0" mode checks exact_tests.phyper(approx=True) and
scipy.stats.hypergeom.{,log}cdf(), treating exact_tests.phyper(approx=False) as
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


def phyper_accuracy_test(ab: float, ac: float, z: float, pow2s: list[int], num_trials_per_pow2: int, bits: int, logp: bool, omit_scipy: bool):
    want = 0.0
    got_scipy = 0.0
    relerr_noapprox_ssq = 0.0
    relerr_approx_ssq = 0.0
    relerr_scipy_ssq = 0.0
    ns = []
    for pow2 in pow2s:
        min_n = 2 ** pow2
        n_limit = min_n * 2
        if num_trials_per_pow2 >= n_limit - min_n:
            ns = range(min_n, n_limit)
        else:
            ns = random.sample(range(min_n, n_limit), k=num_trials_per_pow2)
        for n in ns:
            cur_ab = round(ab * n)
            cur_ac = round(ac * n)
            stdev = 0
            if n > 1:
                variance = cur_ab * cur_ac * (n - cur_ab) * (n - cur_ac) / (n * n * (n - 1))
                stdev = math.sqrt(variance)
            k = round(cur_ab * cur_ac / n)
            got_noapprox = flush_if_denormal(exact_tests.phyper(k, cur_ac, n - cur_ac, cur_ab, logp=logp))
            got_approx = flush_if_denormal(exact_tests.phyper(k, cur_ac, n - cur_ac, cur_ab, approx=True, logp=logp))
            if not omit_scipy:
                if logp:
                    got_scipy = scipy.stats.hypergeom.logcdf(k, n, cur_ab, cur_ac)
                else:
                    got_scipy = scipy.stats.hypergeom.cdf(k, n, cur_ab, cur_ac)
            if bits == 0:
                want = got_noapprox
            else:
                want = flush_if_denormal(mpfr_impls.hypergeom_cdf(k, n, cur_ab, cur_ac, bits=bits, return_log=logp))
                relerr = compute_relerr(got_noapprox, want)
                relerr_noapprox_ssq += relerr * relerr
            relerr = compute_relerr(got_approx, want)
            relerr_approx_ssq += relerr * relerr
            relerr = compute_relerr(got_scipy, want)
            relerr_scipy_ssq += relerr * relerr
            if not omit_scipy:
                if logp:
                    got = scipy.stats.hypergeom.logcdf(k, n, cur_ab, cur_ac)
                else:
                    got = scipy.stats.hypergeom.cdf(k, n, cur_ab, cur_ac)
                relerr = compute_relerr(got, want)
                relerr_scipy_ssq += relerr * relerr
        num_trials = len(ns)
        print_str = f"n in [2^{pow2}, 2^{pow2+1}): "
        if bits > 0:
            noapprox_rms = math.sqrt(relerr_noapprox_ssq / num_trials)
            print_str += f"errRMS={noapprox_rms:.6g}  "
        approx_rms = math.sqrt(relerr_approx_ssq / num_trials)
        print_str += f"approxErrRMS={approx_rms:.6g}"
        if not omit_scipy:
            scipy_rms = math.sqrt(relerr_scipy_ssq / num_trials)
            print_str += f"  scipyErrRMS={scipy_rms:.6g}"
        print(print_str)


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
                             help="Test n~=2**<these values>.")
    optionalarg.add_argument('-n', '--number', type=int, default=10,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-b', '--bits', type=int, default=256,
                             help="MPFR precision; or if 0, test approx=True and scipy against approx=False.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True.")
    optionalarg.add_argument('-o', '--omit-scipy', action="store_true",
                             help="Skip scipy calculations.")
    optionalarg.add_argument('-s', '--seed', type=int, default=1,
                             help="RNG seed.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    random.seed(cmd_args.seed)
    phyper_accuracy_test(cmd_args.row_prop,
                         cmd_args.col_prop,
                         cmd_args.z_score,
                         parse_range_string(cmd_args.exps),
                         cmd_args.number,
                         cmd_args.bits,
                         cmd_args.logp,
                         cmd_args.omit_scipy)


if __name__ == '__main__':
    main()
