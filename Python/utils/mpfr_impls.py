#!/usr/bin/env python3
import gmpy2
import math


# MPFR-based implementations of most functions, for verification purposes.
# Clarity is prioritized over speed here.

def to_float(x):
    bits = gmpy2.get_context().precision
    if bits >= 106:
        # If we have enough bits to work with, perform two-step rounding (first
        # to ~3/4 of bits, then to float64).  This ~eliminates the risk of
        # mishandling cases where the true probability is exactly on the
        # boundary between two float64s.
        bits = (bits * 3) // 4
    return float(gmpy2.round2(x, bits))

def float_from_ln(lnp, return_log):
    if return_log:
        return to_float(lnp)
    return to_float(gmpy2.exp(lnp))

def join_log_and_nonlog(log_part, nonlog_part, return_log, invert):
    # Returns exp(log_part) * nonlog_part if invert=False, subtracts that from
    # 1 if invert=True.
    # Assumes nonlog_part <= 2^52.
    if return_log and not invert:
        return to_float(log_part + gmpy2.log(nonlog_part))

    # Add 52 to log_part and divide nonlog_part by 2^52 so that, for small
    # log_part, we don't risk underflow at the exp(log_part) step when the
    # final result doesn't underflow.
    log_part_shifted = log_part + 52 * gmpy2.const_log2()
    nonlog_part_shifted = gmpy2.mul_2exp(nonlog_part, -52)
    product = gmpy2.exp(log_part_shifted) * nonlog_part_shifted
    if not invert:
        return to_float(product)
    if not return_log:
        return to_float(1 - product)
    return to_float(gmpy2.log1p(-product))

def logspace_add(x, y):
    # Returns ln(exp(x) + exp(y)).
    if x < y:
        x, y = y, x
    return x + gmpy2.log1p(gmpy2.exp(y - x))


def binom_pmf(k: int, n: int, p: float, bits: int, return_log: bool):
    gmpy2.get_context().precision = bits
    p = gmpy2.mpfr(p)
    lnp = gmpy2.log(p)
    q = 1 - p
    lnq = gmpy2.log1p(-p)
    nmk = n - k
    lnpmf = lnp * k + lnq * nmk + gmpy2.lngamma(n + 1) - gmpy2.lngamma(k + 1) - gmpy2.lngamma(nmk + 1)
    return float_from_ln(lnpmf, return_log)


def binom_cdf(k: int, n: int, p: float, bits: int, return_log: bool):
    gmpy2.get_context().precision = bits
    # Assumes 0 < p < 1.
    # Implementation doesn't take advantage of the Aroian/DiDonato/Morris
    # continued fraction, since much of the point of this function is to
    # validate our continued fraction implementation.
    #
    # 1. If we're to the right of the mode, invert.  (Ok if the decision is
    #    imprecise; but 1-p must be represented precisely.)
    # Then, with high precision:
    # 2. Compute logpmf(k).
    # 3. Sum relative likelihoods for k, k-1, k-2, ... until configured
    #    precision limit.
    # 4. If we didn't invert, return exp(logpmf(k)) * lik_sum if
    #    return_log=False, or logpmf(k) + log(lik_sum) if return_log=True.
    #    If we did invert, return 1 - (exp(logpmf(k)) * lik_sum) if
    #    return_log=False, or log1p(-exp(logpmf(k)) * lik_sum) if
    #    return_log=True.
    if k == n:
        # Special-cased so that inversion can't yield k=-1.
        return float_from_ln(0.0, return_log)

    p = gmpy2.mpfr(p)
    # We avoid the variable name 'logp' in this function; 'return_log' and
    # 'lnp' are less ambiguous.
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
    lnpmf = lnp * k + lnq * nmk + gmpy2.lngamma(n + 1) - gmpy2.lngamma(k + 1) - gmpy2.lngamma(nmk + 1)
    qdp = q / p
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

    return join_log_and_nonlog(lnpmf, lik_sum, return_log, invert)


def binomtest(obs_k: int, n: int, p: float, bits: int, return_log: bool):
    gmpy2.get_context().precision = bits
    # Assumes 0 < p < 1.
    # Assumes alternative="two-sided", since the "less" and "greater" cases
    # reduce to single binom_cdf() calls.
    #
    # With high precision:
    # 1. Determine if we're at or to the right of the mode(s).  If we're at a
    #    mode, immediately return 1 or log(1).  If we're to the right, invert.
    # 2. Use binom_cdf() to calculate left-tail probability.
    # 3. Perform binary search for smallest k > mode for which pmf(k) <=
    #    pmf(obs_k).
    # 4. Use binom_cdf() to calculate right-tail probability, return sum of
    #    tails.

    eps_bits = (bits * 3) // 4
    eps = gmpy2.exp2(-eps_bits)

    p = gmpy2.mpfr(p)
    premode = gmpy2.round2((n+1) * p, eps_bits)
    mode = gmpy2.floor(premode)
    if obs_k == mode or (premode.is_integer() and obs_k == mode - 1):
        return float_from_ln(0.0, return_log)

    q = 1 - p
    if obs_k > mode:
        obs_k = n - obs_k
        mode = n - mode
        p, q = q, p

    ln_left = binom_cdf(obs_k, n, p, bits, True)

    # _adj to indicate we've premultiplied by (1 - eps).
    ln_starting_pmf_adj = binom_pmf(obs_k, n, p, bits, True) * (1 - eps)

    ln_pmf_n = binom_pmf(n, n, p, bits, True)
    if ln_pmf_n > ln_starting_pmf_adj:
        # All contingency tables to right of mode have higher probability than
        # our starting table.
        return float_from_ln(ln_left, return_log)

    # Binary search.  Invariant: pmf(min_k-1) > pmf(obs_k) >= pmf(max_k).
    min_k = mode + 1
    max_k = n
    while min_k < max_k:
        k = (min_k + max_k) // 2
        ln_pmf_k = binom_pmf(k, n, p, bits, True)
        if ln_pmf_k > ln_starting_pmf_adj:
            min_k = k + 1
        else:
            max_k = k
    ln_right = binom_cdf(n - min_k, n, q, bits, True)
    ln_total = logspace_add(ln_left, ln_right)
    return float_from_ln(ln_total, return_log)


def hypergeom_lnpmf_internal(a: int, b: int, c: int, d: int):
    ln_numer = gmpy2.lngamma(a+b+1) + gmpy2.lngamma(c+d+1) + gmpy2.lngamma(a+c+1) + gmpy2.lngamma(b+d+1)
    ln_denom = gmpy2.lngamma(a+1) + gmpy2.lngamma(b+1) + gmpy2.lngamma(c+1) + gmpy2.lngamma(d+1) + gmpy2.lngamma(a+b+c+d+1)
    return ln_numer - ln_denom


def hypergeom_pmf(k: int, M: int, n: int, N: int, bits: int, return_log: bool):
    gmpy2.get_context().precision = bits
    a = k
    b = n - k
    c = N - k
    d = M - N - b
    return float_from_ln(hypergeom_lnpmf_internal(a, b, c, d), return_log)


def hypergeom_cdf(k: int, M: int, n: int, N: int, bits: int, return_log: bool):
    # Very similar to binom_cdf().
    #
    # 1. If we're to the right of the mode, invert.
    # Then, with high precision:
    # 2. Compute logpmf(k).
    # 3. Sum relative likelihoods for k, k-1, k-2, ... until configured
    #    precision limit.
    # 4. If we didn't invert, return exp(logpmf(k)) * lik_sum if
    #    return_log=False, or logpmf(k) + log(lik_sum) if return_log=True.
    #    If we did invert, return 1 - (exp(logpmf(k)) * lik_sum) if
    #    return_log=False, or log1p(-exp(logpmf(k)) * lik_sum) if
    #    return_log=True.
    gmpy2.get_context().precision = bits
    a = k
    b = n - k
    c = N - k
    d = M - N - b

    if b == 0 or c == 0:
        # Special-cased so that inversion can't yield a=-1.
        return float_from_ln(0.0, return_log)

    invert = False
    if a * d > b * c:
        invert = True
        a, b, c, d = b-1, a+1, d+1, c-1

    lnpmf = hypergeom_lnpmf_internal(a, b, c, d)
    lik = gmpy2.mpfr(1.0)
    lik_sum = lik
    while True:
        b += 1
        c += 1
        lik *= a * d
        lik /= b * c
        a -= 1
        d -= 1
        preadd = lik_sum
        lik_sum += lik
        if preadd == lik_sum:
            break

    return join_log_and_nonlog(lnpmf, lik_sum, return_log, invert)


def fisher_exact_22(obs_a: int, obs_b: int, obs_c: int, obs_d: int, bits: int, return_log: bool):
    gmpy2.get_context().precision = bits

    eps_bits = (bits * 3) // 4
    eps = gmpy2.exp2(-eps_bits)

    abcd = obs_a + obs_b + obs_c + obs_d
    ab = obs_a + obs_b
    ac = obs_a + obs_c
    denom = gmpy2.mpfr(abcd + 2)
    premode = gmpy2.round2(gmpy2.mpfr((ab + 1) * (ac + 1)) / denom, eps_bits)
    mode = gmpy2.floor(premode)
    if obs_a == mode or (premode.is_integer() and obs_a == mode - 1):
        return float_from_ln(0.0, return_log)

    if obs_a > mode:
        obs_a, obs_b, obs_c, obs_d = obs_b, obs_a, obs_d, obs_c
        ac = obs_a + obs_c
        premode = gmpy2.round2(gmpy2.mpfr((ab + 1) * (ac + 1)) / denom, eps_bits)
        mode = gmpy2.floor(premode)

    d_minus_a = obs_d - obs_a
    min_a = max(0, -d_minus_a)
    max_a = obs_a + min(obs_b, obs_c)

    ln_left = hypergeom_cdf(obs_a, abcd, ab, ac, bits, True)

    # _adj to indicate we've premultiplied by (1 - eps).
    ln_starting_pmf_adj = hypergeom_lnpmf_internal(obs_a, obs_b, obs_c, obs_d) * (1 - eps)

    ln_pmf_maxa = hypergeom_lnpmf_internal(max_a, ab - max_a, ac - max_a, d_minus_a + max_a)
    if ln_pmf_maxa > ln_starting_pmf_adj:
        # All contingency tables to right of mode have higher probability than
        # our starting table.
        return float_from_ln(ln_left, return_log)

    # Binary search.  Invariant: pmf(min_a-1) > pmf(obs_a) >= pmf(max_a).
    min_a = mode + 1
    while min_a < max_a:
        a = (min_a + max_a) // 2
        ln_pmf_a = hypergeom_lnpmf_internal(a, ab - a, ac - a, d_minus_a + a)
        if ln_pmf_a > ln_starting_pmf_adj:
            min_a = a + 1
        else:
            max_a = a
    # Invert: a, b, c, d = b, a, d, c
    # a+b and a+b+c+d remain the same
    ln_right = hypergeom_cdf(ab - min_a, abcd, ab, abcd - ac, bits, True)
    ln_total = logspace_add(ln_left, ln_right)
    return float_from_ln(ln_total, return_log)


# For speed, we use roughly the same Newton-iteration logic as the code we're
# verifying.  This is low-risk since we're dealing with smooth strictly
# increasing/decreasing functions and simple termination conditions; if there's
# a shared Newton-iteration bug, we fail to terminate or obviously
# overflow/underflow, rather than subtly returning a shared wrong result.
def nchypergeom_fisher_approx_mode(m1, m2, n, odds):
    if odds == 1:
        # Avoid division by zero.
        return gmpy2.floor(gmpy2.mpfr((m1 + 1) * (n + 1)) / gmpy2.mpfr(m1 + m2 + 2))
    # Solve for m11 satisfying m11 * m22 == m12 * m21 * odds.
    aa = 1 - odds
    bb = (m1 + n) * odds + (m2 - n)
    neg_c = odds * (m1 * n)
    discrim = bb * bb + 4 * aa * neg_c
    sqrt_discrim = gmpy2.mpfr(0.0)
    if discrim > 0:
        sqrt_discrim = gmpy2.sqrt(discrim)
    return int((sqrt_discrim - bb) / (2 * aa))

def nchypergeom_fisher_mean(m1, m2, n, odds):
    mode = nchypergeom_fisher_approx_mode(m1, m2, n, odds)
    lik = gmpy2.mpfr(1.0)
    m11 = mode
    m12 = m1 - mode
    m21 = n - mode
    m22 = m2 - m21
    # Iterate rightward until convergence.
    rnumer = gmpy2.mpfr(mode)
    rdenom = gmpy2.mpfr(1.0)
    while True:
        m11 += 1
        m22 += 1
        lik *= (m12 * m21) * odds
        lik /= m11 * m22
        m12 -= 1
        m21 -= 1
        preadd = rnumer
        rnumer += lik * m11
        if preadd == rnumer:
            break
        rdenom += lik
    # Jump back to mode, and then iterate leftward until left-sums converge.
    lik = gmpy2.mpfr(1.0)
    m11 = mode
    m12 = m1 - mode
    m21 = n - mode
    m22 = m2 - m21
    lnumer = gmpy2.mpfr(0.0)
    ldenom = gmpy2.mpfr(0.0)
    while True:
        m12 += 1
        m21 += 1
        lik *= m11 * m22
        lik /= (m12 * m21) * odds
        m11 -= 1
        m22 -= 1
        preadd = ldenom
        ldenom += lik
        if preadd == ldenom:
            break
        lnumer += lik * m11
    return (lnumer + rnumer) / (ldenom + rdenom)

def nchypergeom_fisher_variance_from_mean(m1, m2, n, odds, mean):
    total = m1 + m2
    if odds == 1:
        # Avoid division by zero.
        if total < 2:
            return gmpy2.mpfr(0.0)
        return gmpy2.mpfr(m1 * n * (total - m1) * (total - n)) / gmpy2.mpfr(total * total * (total - 1))
    first_term = ((m1 * n) * odds - total * mean) / (1 - odds)
    second_term = mean * (mean - (m1 + n))
    return first_term - second_term

def nchypergeom_fisher_mean_deriv_odds(m1, m2, n, odds, mean):
    return nchypergeom_fisher_variance_from_mean(m1, m2, n, odds, mean) / odds

def nchypergeom_fisher_cdf2(obs_m11, obs_m12, obs_m21, obs_m22, odds1, odds2):
    # Assumes odds1 <= odds2, and they're close to each other.
    lik1 = gmpy2.mpfr(1.0)
    lik2 = gmpy2.mpfr(1.0)
    m11 = obs_m11
    m12 = obs_m12
    m21 = obs_m21
    m22 = obs_m22
    # Iterate rightward until convergence.
    right_sum1 = gmpy2.mpfr(0.0)
    right_sum2 = gmpy2.mpfr(0.0)
    while True:
        m11 += 1
        m22 += 1
        shared_mult = gmpy2.mpfr(m12 * m21) / gmpy2.mpfr(m11 * m22)
        lik1 *= shared_mult * odds1
        lik2 *= shared_mult * odds2
        m12 -= 1
        m21 -= 1
        preadd = right_sum2
        right_sum2 += lik2
        if preadd == right_sum2:
            break
        right_sum1 += lik1
    # Don't need to worry about right_sum2 overflow here.
    lik1 = gmpy2.mpfr(1.0)
    lik2 = gmpy2.mpfr(1.0)
    m11 = obs_m11
    m12 = obs_m12
    m21 = obs_m21
    m22 = obs_m22
    inv_odds1 = 1.0 / odds1
    inv_odds2 = 1.0 / odds2
    left_sum1 = lik1
    left_sum2 = lik2
    while True:
        m12 += 1
        m21 += 1
        shared_mult = gmpy2.mpfr(m11 * m22) / gmpy2.mpfr(m12 * m21)
        lik1 *= shared_mult * inv_odds1
        lik2 *= shared_mult * inv_odds2
        m11 -= 1
        m22 -= 1
        preadd = left_sum1
        left_sum1 += lik1
        if preadd == left_sum1:
            break
        left_sum2 += lik2
    return (left_sum1 / (left_sum1 + right_sum1)), (left_sum2 / (left_sum2 + right_sum2))

def logit(p):
    return gmpy2.log(p / (1-p))

def cond_odds_ratio_quantile_match(m11, m12, m21, m22, target_p):
    # Degenerate cases.
    if (m12 == 0) or (m21 == 0):
        if (m11 == 0) or (m22 == 0):
            return 1.0
        return math.inf
    if target_p <= 0.0:
        return math.inf
    if target_p >= 1.0:
        return 0.0

    bits = gmpy2.get_context().precision
    eps = gmpy2.exp2(-(bits // 2))
    half_eps = gmpy2.exp2(-(1 + (bits // 2)))
    one_plus_half_eps = 1 + half_eps
    one_minus_half_eps = 1 - half_eps
    # Just use initial guess of 1, this forces a different iteration path than
    # code we're verifying
    odds1 = gmpy2.mpfr(1.0)
    logit_target_p = logit(target_p)
    while True:
        odds2 = odds1 * one_plus_half_eps
        odds1 *= one_minus_half_eps
        # As odds ratio increases, p decreases, so we'll get p2 <= p1.
        p1, p2 = nchypergeom_fisher_cdf2(m11, m12, m21, m22, odds1, odds2)
        # Interpolate or take a Newton step, treating logit(p) as a linear
        # function of log-odds.
        logit_p1 = logit(p1)
        logit_p2 = logit(p2)
        if (target_p >= p2) and (target_p <= p1):
            interp = (logit_target_p - logit_p2) / (logit_p1 - logit_p2)
            return odds2 * (1 - eps * interp)
        logitp_deriv_lnodds = (logit_p2 - logit_p1) / eps
        if target_p < p2:
            lnodds_incr = (logit_target_p - logit_p2) / logitp_deriv_lnodds
            odds1 = odds2 * gmpy2.exp(lnodds_incr)
        else:
            lnodds_incr = (logit_target_p - logit_p1) / logitp_deriv_lnodds
            odds1 *= gmpy2.exp(lnodds_incr)


def cond_odds_ratio(m11: int, m12: int, m21: int, m22: int, bits: int):
    gmpy2.get_context().precision = bits
    # Degenerate cases.
    m11_is_maximal = (m12 == 0) or (m21 == 0)
    if (m11 == 0) or (m22 == 0):
        if m11_is_maximal:
            return 1.0
        return 0.0
    if m11_is_maximal:
        return math.inf
    eps = gmpy2.exp2(-(bits // 2))
    m1 = m11 + m12
    m2 = m21 + m22
    n = m11 + m21
    # +0.5 to force different iteration path than code we're verifying
    odds = gmpy2.mpfr(m11 + 0.5) * gmpy2.mpfr(m22 + 0.5) / (gmpy2.mpfr(m12 + 0.5) * gmpy2.mpfr(m21 + 0.5))
    while True:
        mean = nchypergeom_fisher_mean(m1, m2, n, odds)
        mean_delta = m11 - mean
        mean_deriv_odds = nchypergeom_fisher_mean_deriv_odds(m1, m2, n, odds, mean)
        if mean_deriv_odds <= 0:
            raise RuntimeError("non-positive mean_deriv_odds")
        odds_incr = mean_delta / mean_deriv_odds
        odds += odds_incr
        if abs(odds_incr) < odds * eps:
            return odds


# def cond_odds_ratio_quantile_match(m11, m12, m21, m22, target_p):
def cond_odds_ratio_ci(m11: int, m12: int, m21: int, m22: int, low=0.025, high=0.975, bits=256):
    gmpy2.get_context().precision = bits
    low_result = to_float(1.0 / cond_odds_ratio_quantile_match(m12, m11, m22, m21, gmpy2.mpfr(low)))
    high_result = to_float(cond_odds_ratio_quantile_match(m11, m12, m21, m22, 1.0 - gmpy2.mpfr(high)))
    return low_result, high_result


# For 2x3 and larger fisher_exact tables, verification should involve some
# simulation.  An arbitrary-precision calculation either handles larger cases
# too slowly, or will be too likely to share bugs with the code we're trying to
# verify.
#
# However, we'll still want the arbitrary-precision calculation to help us
# identify and plug low-hanging accuracy leaks.


def snphwe(obs_hets: int, hom1: int, hom2: int, bits: int, return_log: bool):
    gmpy2.get_context().precision = bits
    # After one-sided HWE tests are implemented, we'll implement a
    # Levene-Haldane cdf function, and maybe modify this function to call it.
    # But this implementation is efficient enough if we aren't very far from
    # the mode.

    lik = gmpy2.mpfr(1.0)  # Normalize starting table to likelihood 1.
    tail_sum = gmpy2.mpfr(1.0)  # It's in the tail.
    center_sum = gmpy2.mpfr(0.0)
    # Iterate down "left" (low hets) side.
    hets = obs_hets
    homc = max(hom1, hom2)
    homr = min(hom1, hom2)
    while True:
        homc += 1
        homr += 1
        lik *= hets * (hets - 1)
        lik /= 4 * homc * homr
        hets -= 2
        if lik <= gmpy2.mpfr(1.0):
            break
        center_sum += lik
    while True:
        preadd = tail_sum
        tail_sum += lik
        if tail_sum == preadd:
            break
        homc += 1
        homr += 1
        lik *= hets * (hets - 1)
        lik /= 4 * homc * homr
        hets -= 2
    # Jump back to starting point, iterate down right side.
    lik = gmpy2.mpfr(1.0)
    hets = obs_hets
    homc = max(hom1, hom2)
    homr = min(hom1, hom2)
    while True:
        hets += 2
        lik *= 4 * homc * homr
        lik /= hets * (hets - 1)
        homc -= 1
        homr -= 1
        if lik <= gmpy2.mpfr(1.0):
            break
        center_sum += lik
    while True:
        preadd = tail_sum
        tail_sum += lik
        if tail_sum == preadd:
            break
        hets += 2
        lik *= 4 * homc * homr
        lik /= hets * (hets - 1)
        homc -= 1
        homr -= 1
    result = tail_sum / (tail_sum + center_sum)
    if return_log:
        return to_float(gmpy2.log(result))
    return to_float(result)
