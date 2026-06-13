#!/usr/bin/env python3
import gmpy2


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


def dbinom(k: int, n: int, p: float, bits: int, return_log: bool):
    gmpy2.get_context().precision = bits
    p = gmpy2.mpfr(p)
    lnp = gmpy2.log(p)
    q = 1 - p
    lnq = gmpy2.log1p(-p)
    nmk = n - k
    lnpmf = lnp * k + lnq * nmk + gmpy2.lngamma(n + 1) - gmpy2.lngamma(k + 1) - gmpy2.lngamma(nmk + 1)
    return float_from_ln(lnpmf, return_log)


def pbinom(k: int, n: int, p: float, bits: int, return_log: bool):
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

    gmpy2.get_context().precision = bits

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
    # Assumes 0 < p < 1.
    # Assumes alternative="two-sided", since the "less" and "greater" cases
    # reduce to single pbinom() calls.
    #
    # With high precision:
    # 1. Determine if we're at or to the right of the mode(s).  If we're at a
    #    mode, immediately return 1 or log(1).  If we're to the right, invert.
    # 2. Use pbinom() to calculate left-tail probability.
    # 3. Perform binary search for smallest k > mode for which pmf(k) <=
    #    pmf(obs_k).
    # 4. Use pbinom() to calculate right-tail probability, return sum of tails.
    gmpy2.get_context().precision = bits

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

    ln_left = pbinom(obs_k, n, p, bits, True)

    # _adj to indicate we've premultiplied by (1 - eps).
    ln_starting_pmf_adj = dbinom(obs_k, n, p, bits, True) * (1 - eps)

    ln_pmf_n = dbinom(n, n, p, bits, True)
    if ln_pmf_n > ln_starting_pmf_adj:
        # All contingency tables to right of mode have higher probability than
        # our starting table.
        return float_from_ln(ln_left, return_log)

    # Binary search.  Invariant: pmf(min_k-1) > pmf(obs_k) >= pmf(max_k).
    min_k = mode + 1
    max_k = n
    while min_k < max_k:
        k = (min_k + max_k) // 2
        ln_pmf_k = dbinom(k, n, p, bits, True)
        if ln_pmf_k > ln_starting_pmf_adj:
            min_k = k + 1
        else:
            max_k = k
    ln_right = pbinom(n - min_k, n, q, bits, True)
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
    # Very similar to pbinom().
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


# todo: at least 2x2 fisher_exact
