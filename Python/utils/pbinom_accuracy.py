#!/usr/bin/env python3
import argparse
import exact_tests
import gmpy2
import math
import random
import scipy
"""
Primary mode checks exact_tests.pbinom(approx=False) and
scipy.stats.binom.cdf()'s accuracy against a simple MPFR-based implementation
(slow!).

"--bits 0" mode checks exact_tests.pbinom(approx=True) and
scipy.statas.binom.{,log}cdf(), treating exact_tests.pbinom(approx=False) as
ground truth.
"""

def pbinom_mpfr(k: int, n: int, p: float, bits: int, logp: bool):
    # (Doesn't use continued fraction, since much of the point of this function
    # is to validate our continued fraction implementation.)
    #
    # 1. If we're to the right of the mode, invert.  (Ok if the decision is
    #    imprecise; but 1-p must be represented precisely.)
    # Then, with high precision:
    # 2. Compute logpmf(k).
    # 3. Sum relative likelihoods for k, k-1, k-2, ... until configured
    #    precision limit.
    # 4. If we didn't invert, return exp(logpmf(k)) * lik_sum if logp=False, or
    #    logpmf(k) + log(lik_sum) if logp=True.
    #    If we did invert, return 1 - (exp(logpmf(k)) * lik_sum) if logp=False,
    #    or log1p(-exp(logpmf(k)) * lik_sum) if logp=True.
    if k == n:
        # Special-cased so that inversion can't yield k=-1.
        if logp:
            return 0.0
        return 1.0

    gmpy2.get_context().precision = bits

    p = gmpy2.mpfr(p)
    # must not confuse this with logp parameter...
    lnp = gmpy2.log(p)
    q = 1 - p
    lnq = gmpy2.log1p(-p)

    invert = False
    if k > n * p:
        invert = True
        k = n - k - 1
        p, q = q, p
        lnp, lnq = lnq, lnp

    nmk = n - k
    qdp = q / p
    lnpmf = lnp * k + lnq * nmk + gmpy2.lngamma(n + 1) - gmpy2.lngamma(k + 1) - gmpy2.lngamma(nmk + 1)
    lik = gmpy2.mpfr(1.0)
    lik_sum = lik
    while True:
        nmk += 1
        lik *= (qdp * k) / nmk
        k -= 1
        preadd = lik_sum
        lik_sum += lik
        if preadd == lik_sum:
            break

    if logp and not invert:
        return float(lnpmf + gmpy2.log(lik_sum))

    # For n in [0, 2^52), lik_sum is in [1, 2^52].
    # Add 52 to lnpmf and divide lik_sum by 2^52 so that, for small lnpmf, we
    # don't risk underflow at the exp(lnpmf) step when the final result
    # doesn't underflow.
    lnpmf_shifted = lnpmf + 52 * gmpy2.const_log2()
    lik_sum_shifted = gmpy2.mul_2exp(lik_sum, -52)
    product = gmpy2.exp(lnpmf_shifted) * lik_sum_shifted
    if not invert:
        return float(product)
    if not logp:
        return float(1 - product)
    return float(gmpy2.log1p(-product))


def pbinom_accuracy_test(p: float, z: float, min_pow2: int, max_pow2: int, num_trials_per_pow2: int, bits: int, logp: bool):
    pq = p * (1.0 - p)
    want = 0.0
    got = 0.0
    test_str = "RMS="
    if bits == 0:
        test_str = "approxRMS="
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
            if bits == 0:
                want = exact_tests.pbinom(k, n, p, logp=logp)
                got = exact_tests.pbinom(k, n, p, approx=True, logp=logp)
            else:
                want = pbinom_mpfr(k, n, p, bits=bits, logp=logp)
                got = exact_tests.pbinom(k, n, p, logp=logp)
            relerr = (got - want) / want
            relerr_ssq += relerr * relerr
            if logp:
                got = scipy.stats.binom.logcdf(k, n, p)
            else:
                got = scipy.stats.binom.cdf(k, n, p)
            relerr = (got - want) / want
            relerr_scipy_ssq += relerr * relerr
        test_rms = math.sqrt(relerr_ssq / num_trials_per_pow2)
        scipy_rms = math.sqrt(relerr_scipy_ssq / num_trials_per_pow2)
        print("[2^" + str(pow2) + ", 2^" + str(pow2+1) + " - 1): " + test_str + str(test_rms) + "  scipyRMS=" + str(scipy_rms))


def parse_commandline_args():
    parser = argparse.ArgumentParser(description=__doc__)
    optionalarg = parser.add_argument_group('Optional Arguments')
    optionalarg.add_argument('-p', '--succ-prob', type=float, default=0.5,
                             help="Binomial distribution success-probability to test.")
    optionalarg.add_argument('-z', '--z-score', type=float, default=0.0,
                             help="Success-count z-score to test.")
    optionalarg.add_argument('-f', '--from-pow2', type=int, default=8,
                             help="Start testing at n=2**<this value>.")
    optionalarg.add_argument('-t', '--to-pow2', type=int, default=33,
                             help="Continue testing up to n=2**(<this value>+1) - 1.")
    optionalarg.add_argument('-n', '--number', type=int, default=25,
                             help="Number of trials per power-of-2 tier.")
    optionalarg.add_argument('-b', '--bits', type=int, default=256,
                             help="MPFR precision; or if 0, test approx=True and scipy against approx=False.")
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
                         cmd_args.bits,
                         cmd_args.logp)


if __name__ == '__main__':
    main()
