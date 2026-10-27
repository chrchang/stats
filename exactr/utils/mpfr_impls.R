suppressPackageStartupMessages(library(Rmpfr))

# MPFR-based implementations of most functions, for verification purposes.
# This is a port of Python/utils/mpfr_impls.py.  Clarity is prioritized over
# speed, with one exception: long series are summed in vectorized chunks (see
# series_sums()), since elementwise Rmpfr operations are ~50x slower than
# vectorized ones.
#
# Unlike the Python version, intermediate log-probabilities are kept in MPFR
# form instead of being rounded to float64 between steps.

to_float <- function(x, bits) {
  if (bits >= 106) {
    # If we have enough bits to work with, perform two-step rounding (first to
    # ~3/4 of bits, then to float64).  This ~eliminates the risk of mishandling
    # cases where the true probability is exactly on the boundary between two
    # float64s.
    bits <- (bits * 3) %/% 4
  }
  asNumeric(roundMpfr(x, bits))
}

float_from_ln <- function(lnp, bits, return_log) {
  if (return_log) {
    return(to_float(lnp, bits))
  }
  to_float(exp(lnp), bits)
}

join_log_and_nonlog <- function(log_part, nonlog_part, bits, return_log, invert) {
  # Returns exp(log_part) * nonlog_part if invert=FALSE, subtracts that from 1
  # if invert=TRUE.
  # Assumes nonlog_part <= 2^52.
  if (return_log && !invert) {
    return(to_float(log_part + log(nonlog_part), bits))
  }

  # Add 52 to log_part and divide nonlog_part by 2^52 so that, for small
  # log_part, we don't risk underflow at the exp(log_part) step when the final
  # result doesn't underflow.
  log_part_shifted <- log_part + 52 * Const("log2", bits)
  nonlog_part_shifted <- nonlog_part / 2^52
  product <- exp(log_part_shifted) * nonlog_part_shifted
  if (!invert) {
    return(to_float(product, bits))
  }
  if (!return_log) {
    return(to_float(1 - product, bits))
  }
  to_float(log1p(-product), bits)
}

logspace_add <- function(x, y) {
  # Returns ln(exp(x) + exp(y)).
  if (x < y) {
    tmp <- x
    x <- y
    y <- tmp
  }
  x + log1p(exp(y - x))
}

mpfr_vec <- function(x, bits) {
  mpfr(x, precBits = bits)
}

# Let lik_j := prod_{i=1}^j (ratio_fn(i) * mult), for each mult in mults.
# Returns list(sums, wsums), where sums[m] = sum_{j>=1} lik_j and (if weight_fn
# is provided) wsums[m] = sum_{j>=1} weight_fn(j) * lik_j.
#
# ratio_fn() must accept a vector of indexes and return an mpfr vector; it is
# never called with an index > max_terms.  The sequence of terms must be
# unimodal; summation stops once the last term of a chunk is negligible at the
# working precision.
series_sums <- function(ratio_fn, max_terms, bits, mults = list(mpfr(1, bits)), weight_fn = NULL) {
  nmult <- length(mults)
  sums <- rep(list(mpfr(0, bits)), nmult)
  wsums <- rep(list(mpfr(0, bits)), nmult)
  carries <- rep(list(mpfr(1, bits)), nmult)
  converged <- rep(FALSE, nmult)
  start <- 1
  chunk_size <- 256
  while ((start <= max_terms) && !all(converged)) {
    end <- min(start + chunk_size - 1, max_terms)
    idxs <- start:end
    ratios <- ratio_fn(idxs)
    weights <- if (is.null(weight_fn)) NULL else weight_fn(idxs)
    for (m in seq_len(nmult)) {
      if (converged[m]) {
        next
      }
      liks <- cumprod(c(carries[[m]], ratios * mults[[m]]))[-1]
      chunk_sum <- sum(liks)
      last_lik <- liks[length(liks)]
      presum <- sums[[m]]
      sums[[m]] <- presum + chunk_sum
      done <- (sums[[m]] + last_lik == sums[[m]])
      if (!is.null(weights)) {
        wsums[[m]] <- wsums[[m]] + sum(liks * weights)
        last_wlik <- last_lik * weights[length(weights)]
        done <- done && (wsums[[m]] + last_wlik == wsums[[m]])
      }
      converged[m] <- done
      carries[[m]] <- last_lik
    }
    start <- end + 1
    chunk_size <- min(chunk_size * 2, 65536)
  }
  list(sums = sums, wsums = wsums)
}


binom_lnpmf_internal <- function(k, n, p, bits) {
  p <- mpfr(p, bits)
  nmk <- n - k
  log(p) * k + log1p(-p) * nmk + lgamma(mpfr(n + 1, bits)) - lgamma(mpfr(k + 1, bits)) - lgamma(mpfr(nmk + 1, bits))
}

binom_pmf <- function(k, n, p, bits, return_log) {
  float_from_ln(binom_lnpmf_internal(k, n, p, bits), bits, return_log)
}


# Returns list(ln_part, nonlog_part, invert); see join_log_and_nonlog().
binom_cdf_parts <- function(k, n, p, bits) {
  # Assumes 0 < p < 1, k < n.
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
  p <- mpfr(p, bits)
  lnp <- log(p)
  q <- 1 - p
  lnq <- log1p(-p)

  invert <- FALSE
  if (k > n * p) {
    invert <- TRUE
    k <- n - k - 1
    tmp <- p
    p <- q
    q <- tmp
    tmp <- lnp
    lnp <- lnq
    lnq <- tmp
  }

  nmk <- n - k
  lnpmf <- lnp * k + lnq * nmk + lgamma(mpfr(n + 1, bits)) - lgamma(mpfr(k + 1, bits)) - lgamma(mpfr(nmk + 1, bits))
  qdp <- q / p
  # lik(k-j) / lik(k-j+1) = (q/p) * (k-j+1) / (n-k+j)
  tail <- series_sums(function(j) qdp * mpfr_vec(k - j + 1, bits) / mpfr_vec(nmk + j, bits),
                      k, bits)
  list(ln_part = lnpmf, nonlog_part = 1 + tail$sums[[1]], invert = invert)
}

binom_cdf <- function(k, n, p, bits, return_log) {
  if (k == n) {
    # Special-cased so that inversion can't yield k=-1.
    return(float_from_ln(mpfr(0, bits), bits, return_log))
  }
  parts <- binom_cdf_parts(k, n, p, bits)
  # 4. If we didn't invert, return exp(logpmf(k)) * lik_sum if
  #    return_log=FALSE, or logpmf(k) + log(lik_sum) if return_log=TRUE.
  #    If we did invert, return 1 - (exp(logpmf(k)) * lik_sum) if
  #    return_log=FALSE, or log1p(-exp(logpmf(k)) * lik_sum) if
  #    return_log=TRUE.
  join_log_and_nonlog(parts$ln_part, parts$nonlog_part, bits, return_log, parts$invert)
}

# Natural log of left-tail cdf, as an mpfr.  Assumes k is to the left of the
# mode.
binom_ln_left_cdf <- function(k, n, p, bits) {
  parts <- binom_cdf_parts(k, n, p, bits)
  stopifnot(!parts$invert)
  parts$ln_part + log(parts$nonlog_part)
}


binomtest <- function(obs_k, n, p, bits, return_log) {
  # Assumes 0 < p < 1.
  # Assumes alternative="two-sided", since the "less" and "greater" cases
  # reduce to single binom_cdf() calls.
  #
  # With high precision:
  # 1. Determine if we're at or to the right of the mode(s).  If we're at a
  #    mode, immediately return 1 or log(1).  If we're to the right, invert.
  # 2. Calculate left-tail probability.
  # 3. Perform binary search for smallest k > mode for which pmf(k) <=
  #    pmf(obs_k).
  # 4. Calculate right-tail probability, return sum of tails.
  eps_bits <- (bits * 3) %/% 4
  eps <- 2^mpfr(-eps_bits, bits)

  p <- mpfr(p, bits)
  premode <- roundMpfr((n + 1) * p, eps_bits)
  mode <- asNumeric(floor(premode))
  if ((obs_k == mode) || ((floor(premode) == premode) && (obs_k == mode - 1))) {
    return(float_from_ln(mpfr(0, bits), bits, return_log))
  }

  q <- 1 - p
  if (obs_k > mode) {
    obs_k <- n - obs_k
    mode <- n - mode
    tmp <- p
    p <- q
    q <- tmp
  }

  ln_left <- binom_ln_left_cdf(obs_k, n, p, bits)

  # _adj to indicate we've premultiplied by (1 - eps).
  ln_starting_pmf_adj <- binom_lnpmf_internal(obs_k, n, p, bits) * (1 - eps)

  ln_pmf_n <- binom_lnpmf_internal(n, n, p, bits)
  if (ln_pmf_n > ln_starting_pmf_adj) {
    # All outcomes to right of mode have higher probability than our starting
    # outcome.
    return(float_from_ln(ln_left, bits, return_log))
  }

  # Binary search.  Invariant: pmf(min_k-1) > pmf(obs_k) >= pmf(max_k).
  min_k <- mode + 1
  max_k <- n
  while (min_k < max_k) {
    k <- (min_k + max_k) %/% 2
    ln_pmf_k <- binom_lnpmf_internal(k, n, p, bits)
    if (ln_pmf_k > ln_starting_pmf_adj) {
      min_k <- k + 1
    } else {
      max_k <- k
    }
  }
  ln_right <- binom_ln_left_cdf(n - min_k, n, q, bits)
  ln_total <- logspace_add(ln_left, ln_right)
  float_from_ln(ln_total, bits, return_log)
}


hypergeom_lnpmf_internal <- function(a, b, c, d, bits) {
  ln_numer <- lgamma(mpfr(a + b + 1, bits)) + lgamma(mpfr(c + d + 1, bits)) + lgamma(mpfr(a + c + 1, bits)) + lgamma(mpfr(b + d + 1, bits))
  ln_denom <- lgamma(mpfr(a + 1, bits)) + lgamma(mpfr(b + 1, bits)) + lgamma(mpfr(c + 1, bits)) + lgamma(mpfr(d + 1, bits)) + lgamma(mpfr(a + b + c + d + 1, bits))
  ln_numer - ln_denom
}


# scipy-style parameterization: M = total, n = first-row total, N = first-col
# total.
hypergeom_pmf <- function(k, M, n, N, bits, return_log) {
  a <- k
  b <- n - k
  c <- N - k
  d <- M - N - b
  float_from_ln(hypergeom_lnpmf_internal(a, b, c, d, bits), bits, return_log)
}


# Returns list(ln_part, nonlog_part, invert); see join_log_and_nonlog().
hypergeom_cdf_parts <- function(a, b, c, d, bits) {
  # Very similar to binom_cdf_parts().
  invert <- FALSE
  if (mpfr(a, bits) * d > mpfr(b, bits) * c) {
    invert <- TRUE
    new_a <- b - 1
    new_b <- a + 1
    new_c <- d + 1
    new_d <- c - 1
    a <- new_a
    b <- new_b
    c <- new_c
    d <- new_d
  }

  lnpmf <- hypergeom_lnpmf_internal(a, b, c, d, bits)
  # lik(a-j) / lik(a-j+1) = (a-j+1)(d-j+1) / ((b+j)(c+j))
  tail <- series_sums(function(j) {
    (mpfr_vec(a - j + 1, bits) * mpfr_vec(d - j + 1, bits)) / (mpfr_vec(b + j, bits) * mpfr_vec(c + j, bits))
  }, min(a, d), bits)
  list(ln_part = lnpmf, nonlog_part = 1 + tail$sums[[1]], invert = invert)
}

hypergeom_cdf <- function(k, M, n, N, bits, return_log) {
  a <- k
  b <- n - k
  c <- N - k
  d <- M - N - b

  if ((b == 0) || (c == 0)) {
    # Special-cased so that inversion can't yield a=-1.
    return(float_from_ln(mpfr(0, bits), bits, return_log))
  }
  parts <- hypergeom_cdf_parts(a, b, c, d, bits)
  join_log_and_nonlog(parts$ln_part, parts$nonlog_part, bits, return_log, parts$invert)
}

# Natural log of left-tail cdf, as an mpfr.  Assumes we're to the left of the
# mode.
hypergeom_ln_left_cdf <- function(k, M, n, N, bits) {
  parts <- hypergeom_cdf_parts(k, n - k, N - k, M - N - n + k, bits)
  stopifnot(!parts$invert)
  parts$ln_part + log(parts$nonlog_part)
}


fisher_exact_22 <- function(obs_a, obs_b, obs_c, obs_d, bits, return_log) {
  eps_bits <- (bits * 3) %/% 4
  eps <- 2^mpfr(-eps_bits, bits)

  abcd <- obs_a + obs_b + obs_c + obs_d
  ab <- obs_a + obs_b
  ac <- obs_a + obs_c
  denom <- mpfr(abcd + 2, bits)
  premode <- roundMpfr(mpfr(ab + 1, bits) * (ac + 1) / denom, eps_bits)
  mode <- asNumeric(floor(premode))
  if ((obs_a == mode) || ((floor(premode) == premode) && (obs_a == mode - 1))) {
    return(float_from_ln(mpfr(0, bits), bits, return_log))
  }

  if (obs_a > mode) {
    tmp <- obs_a
    obs_a <- obs_b
    obs_b <- tmp
    tmp <- obs_c
    obs_c <- obs_d
    obs_d <- tmp
    ac <- obs_a + obs_c
    premode <- roundMpfr(mpfr(ab + 1, bits) * (ac + 1) / denom, eps_bits)
    mode <- asNumeric(floor(premode))
  }

  d_minus_a <- obs_d - obs_a
  max_a <- obs_a + min(obs_b, obs_c)

  ln_left <- hypergeom_ln_left_cdf(obs_a, abcd, ab, ac, bits)

  # _adj to indicate we've premultiplied by (1 - eps).
  ln_starting_pmf_adj <- hypergeom_lnpmf_internal(obs_a, obs_b, obs_c, obs_d, bits) * (1 - eps)

  ln_pmf_maxa <- hypergeom_lnpmf_internal(max_a, ab - max_a, ac - max_a, d_minus_a + max_a, bits)
  if (ln_pmf_maxa > ln_starting_pmf_adj) {
    # All contingency tables to right of mode have higher probability than our
    # starting table.
    return(float_from_ln(ln_left, bits, return_log))
  }

  # Binary search.  Invariant: pmf(min_a-1) > pmf(obs_a) >= pmf(max_a).
  min_a <- mode + 1
  while (min_a < max_a) {
    a <- (min_a + max_a) %/% 2
    ln_pmf_a <- hypergeom_lnpmf_internal(a, ab - a, ac - a, d_minus_a + a, bits)
    if (ln_pmf_a > ln_starting_pmf_adj) {
      min_a <- a + 1
    } else {
      max_a <- a
    }
  }
  # Invert: a, b, c, d = b, a, d, c
  # a+b and a+b+c+d remain the same
  ln_right <- hypergeom_ln_left_cdf(ab - min_a, abcd, ab, abcd - ac, bits)
  ln_total <- logspace_add(ln_left, ln_right)
  float_from_ln(ln_total, bits, return_log)
}


# For speed, we use roughly the same Newton-iteration logic as the code we're
# verifying.  This is low-risk since we're dealing with smooth strictly
# increasing/decreasing functions and simple termination conditions; if there's
# a shared Newton-iteration bug, we fail to terminate or obviously
# overflow/underflow, rather than subtly returning a shared wrong result.
nchypergeom_fisher_approx_mode <- function(m1, m2, n, odds, bits) {
  if (odds == 1) {
    # Avoid division by zero.
    return(asNumeric(floor(mpfr(m1 + 1, bits) * (n + 1) / mpfr(m1 + m2 + 2, bits))))
  }
  # Solve for m11 satisfying m11 * m22 == m12 * m21 * odds.
  aa <- 1 - odds
  bb <- (m1 + n) * odds + (m2 - n)
  neg_c <- odds * (mpfr(m1, bits) * n)
  discrim <- bb * bb + 4 * aa * neg_c
  sqrt_discrim <- mpfr(0, bits)
  if (discrim > 0) {
    sqrt_discrim <- sqrt(discrim)
  }
  asNumeric(trunc((sqrt_discrim - bb) / (2 * aa)))
}

nchypergeom_fisher_mean <- function(m1, m2, n, odds, bits) {
  mode <- nchypergeom_fisher_approx_mode(m1, m2, n, odds, bits)
  m11 <- mode
  m12 <- m1 - mode
  m21 <- n - mode
  m22 <- m2 - m21
  # Iterate rightward until convergence.
  right <- series_sums(function(j) {
    (mpfr_vec(m12 - j + 1, bits) * mpfr_vec(m21 - j + 1, bits)) / (mpfr_vec(m11 + j, bits) * mpfr_vec(m22 + j, bits))
  }, min(m12, m21), bits, mults = list(odds), weight_fn = function(j) mpfr_vec(m11 + j, bits))
  rnumer <- mode + right$wsums[[1]]
  rdenom <- 1 + right$sums[[1]]
  # Iterate leftward from mode until convergence.
  left <- series_sums(function(j) {
    (mpfr_vec(m11 - j + 1, bits) * mpfr_vec(m22 - j + 1, bits)) / (mpfr_vec(m12 + j, bits) * mpfr_vec(m21 + j, bits))
  }, min(m11, m22), bits, mults = list(1 / odds), weight_fn = function(j) mpfr_vec(m11 - j, bits))
  lnumer <- left$wsums[[1]]
  ldenom <- left$sums[[1]]
  (lnumer + rnumer) / (ldenom + rdenom)
}

nchypergeom_fisher_variance_from_mean <- function(m1, m2, n, odds, mean, bits) {
  total <- m1 + m2
  if (odds == 1) {
    # Avoid division by zero.
    if (total < 2) {
      return(mpfr(0, bits))
    }
    return(mpfr(m1, bits) * n * (total - m1) * (total - n) / (mpfr(total, bits) * total * (total - 1)))
  }
  first_term <- ((mpfr(m1, bits) * n) * odds - total * mean) / (1 - odds)
  second_term <- mean * (mean - (m1 + n))
  first_term - second_term
}

nchypergeom_fisher_mean_deriv_odds <- function(m1, m2, n, odds, mean, bits) {
  nchypergeom_fisher_variance_from_mean(m1, m2, n, odds, mean, bits) / odds
}

nchypergeom_fisher_cdf2 <- function(obs_m11, obs_m12, obs_m21, obs_m22, odds1, odds2, bits) {
  # Assumes odds1 <= odds2, and they're close to each other.
  # Iterate rightward until convergence.
  right <- series_sums(function(j) {
    (mpfr_vec(obs_m12 - j + 1, bits) * mpfr_vec(obs_m21 - j + 1, bits)) / (mpfr_vec(obs_m11 + j, bits) * mpfr_vec(obs_m22 + j, bits))
  }, min(obs_m12, obs_m21), bits, mults = list(odds1, odds2))
  right_sum1 <- right$sums[[1]]
  right_sum2 <- right$sums[[2]]
  # Iterate leftward until convergence.
  left <- series_sums(function(j) {
    (mpfr_vec(obs_m11 - j + 1, bits) * mpfr_vec(obs_m22 - j + 1, bits)) / (mpfr_vec(obs_m12 + j, bits) * mpfr_vec(obs_m21 + j, bits))
  }, min(obs_m11, obs_m22), bits, mults = list(1 / odds1, 1 / odds2))
  left_sum1 <- 1 + left$sums[[1]]
  left_sum2 <- 1 + left$sums[[2]]
  list(left_sum1 / (left_sum1 + right_sum1), left_sum2 / (left_sum2 + right_sum2))
}

logit <- function(p) {
  log(p / (1 - p))
}

cond_odds_ratio_quantile_match <- function(m11, m12, m21, m22, target_p, bits) {
  # Degenerate cases.
  if ((m12 == 0) || (m21 == 0)) {
    if ((m11 == 0) || (m22 == 0)) {
      return(mpfr(1, bits))
    }
    return(mpfr(Inf, bits))
  }
  if (target_p <= 0) {
    return(mpfr(Inf, bits))
  }
  if (target_p >= 1) {
    return(mpfr(0, bits))
  }

  eps <- 2^mpfr(-(bits %/% 2), bits)
  half_eps <- 2^mpfr(-(1 + (bits %/% 2)), bits)
  one_plus_half_eps <- 1 + half_eps
  one_minus_half_eps <- 1 - half_eps
  # Just use initial guess of 1, this forces a different iteration path than
  # code we're verifying
  odds1 <- mpfr(1, bits)
  logit_target_p <- logit(target_p)
  repeat {
    odds2 <- odds1 * one_plus_half_eps
    odds1 <- odds1 * one_minus_half_eps
    # As odds ratio increases, p decreases, so we'll get p2 <= p1.
    ps <- nchypergeom_fisher_cdf2(m11, m12, m21, m22, odds1, odds2, bits)
    p1 <- ps[[1]]
    p2 <- ps[[2]]
    # Interpolate or take a Newton step, treating logit(p) as a linear
    # function of log-odds.
    logit_p1 <- logit(p1)
    logit_p2 <- logit(p2)
    if ((target_p >= p2) && (target_p <= p1)) {
      interp <- (logit_target_p - logit_p2) / (logit_p1 - logit_p2)
      return(odds2 * (1 - eps * interp))
    }
    logitp_deriv_lnodds <- (logit_p2 - logit_p1) / eps
    if (target_p < p2) {
      lnodds_incr <- (logit_target_p - logit_p2) / logitp_deriv_lnodds
      odds1 <- odds2 * exp(lnodds_incr)
    } else {
      lnodds_incr <- (logit_target_p - logit_p1) / logitp_deriv_lnodds
      odds1 <- odds1 * exp(lnodds_incr)
    }
  }
}


cond_odds_ratio <- function(m11, m12, m21, m22, bits) {
  # Degenerate cases.
  m11_is_maximal <- (m12 == 0) || (m21 == 0)
  if ((m11 == 0) || (m22 == 0)) {
    if (m11_is_maximal) {
      return(1)
    }
    return(0)
  }
  if (m11_is_maximal) {
    return(Inf)
  }
  eps <- 2^mpfr(-(bits %/% 2), bits)
  m1 <- m11 + m12
  m2 <- m21 + m22
  n <- m11 + m21
  # +0.5 to force different iteration path than code we're verifying
  odds <- mpfr(m11 + 0.5, bits) * (m22 + 0.5) / (mpfr(m12 + 0.5, bits) * (m21 + 0.5))
  repeat {
    mean <- nchypergeom_fisher_mean(m1, m2, n, odds, bits)
    mean_delta <- m11 - mean
    mean_deriv_odds <- nchypergeom_fisher_mean_deriv_odds(m1, m2, n, odds, mean, bits)
    if (mean_deriv_odds <= 0) {
      stop("non-positive mean_deriv_odds")
    }
    odds_incr <- mean_delta / mean_deriv_odds
    odds <- odds + odds_incr
    if (abs(odds_incr) < odds * eps) {
      return(to_float(odds, bits))
    }
  }
}


cond_odds_ratio_ci <- function(m11, m12, m21, m22, low = 0.025, high = 0.975, bits = 256) {
  low_result <- to_float(1 / cond_odds_ratio_quantile_match(m12, m11, m22, m21, mpfr(low, bits), bits), bits)
  high_result <- to_float(cond_odds_ratio_quantile_match(m11, m12, m21, m22, 1 - mpfr(high, bits), bits), bits)
  c(low_result, high_result)
}
