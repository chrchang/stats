#!/usr/bin/env python3
import argparse
import exact_tests
import mpfr_impls
import math
import random
import scipy
import sys
"""
Primary mode checks exact_tests.ppois(approx=False), exact_tests.ppois(approx=True),
and scipy.stats.poisson.cdf()'s accuracy against a simple MPFR-based
implementation (slow!).

"--bits 0" mode checks exact_tests.ppois(approx=True) and
scipy.stats.poisson.{,log}cdf(), treating exact_tests.ppois(approx=False) as
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


def ppois_accuracy_test(z: float, pow2s: list[int], num_trials_per_pow2: int, bits: int, logp: bool):
    want = 0.0
    for pow2 in pow2s:
        relerr_noapprox_ssq = 0.0
        relerr_approx_ssq = 0.0
        relerr_scipy_ssq = 0.0
        for _ in range(num_trials_per_pow2):
            # Non-integer mean in [2^pow2, 2^(pow2+1)).
            mu = (2 ** pow2) * (1 + random.random())
            k = round(mu + z * math.sqrt(mu) + 0.5)
            if k < 0:
                k = 0
            got_noapprox = flush_if_denormal(exact_tests.ppois(k, mu, logp=logp))
            got_approx = flush_if_denormal(exact_tests.ppois(k, mu, approx=True, logp=logp))
            if logp:
                got_scipy = scipy.stats.poisson.logcdf(k, mu)
            else:
                got_scipy = scipy.stats.poisson.cdf(k, mu)
            if bits == 0:
                want = got_noapprox
            else:
                want = flush_if_denormal(mpfr_impls.pois_cdf(k, mu, bits=bits, return_log=logp))
                relerr = compute_relerr(got_noapprox, want)
                relerr_noapprox_ssq += relerr * relerr
            relerr = compute_relerr(got_approx, want)
            relerr_approx_ssq += relerr * relerr
            relerr = compute_relerr(got_scipy, want)
            relerr_scipy_ssq += relerr * relerr
        approx_rms = math.sqrt(relerr_approx_ssq / num_trials_per_pow2)
        scipy_rms = math.sqrt(relerr_scipy_ssq / num_trials_per_pow2)
        if bits == 0:
            print(f"mu in [2^{pow2}, 2^{pow2+1}): approxErrRMS={approx_rms:.3g}  scipyErrRMS={scipy_rms:.3g}")
        else:
            noapprox_rms = math.sqrt(relerr_noapprox_ssq / num_trials_per_pow2)
            print(f"mu in [2^{pow2}, 2^{pow2+1}): errRMS={noapprox_rms:.3g}  approxErrRMS={approx_rms:.3g}  scipyErrRMS={scipy_rms:.3g}")


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-z', '--z-score', type=float, default=0.0,
                             help="Event-count z-score to test.")
    optionalarg.add_argument('-e', '--exps', type=str, default="5,20,30",
                             help="Test mu~=2**<these values>.")
    optionalarg.add_argument('-n', '--number', type=int, default=10,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-b', '--bits', type=int, default=256,
                             help="MPFR precision; or if 0, treat approx=False as ground truth.")
    optionalarg.add_argument('-l', '--logp', action="store_true",
                             help="Test logp=True.")
    optionalarg.add_argument('-s', '--seed', type=int, default=1,
                             help="RNG seed.")
    cmd_args = parser.parse_args()
    return cmd_args


def main():
    cmd_args = parse_commandline_args()
    random.seed(cmd_args.seed)
    ppois_accuracy_test(cmd_args.z_score,
                        parse_range_string(cmd_args.exps),
                        cmd_args.number,
                        cmd_args.bits,
                        cmd_args.logp)


if __name__ == '__main__':
    main()
