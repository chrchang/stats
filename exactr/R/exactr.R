is.boolean <- function(x) {
  is.logical(x) && length(x) == 1 && !is.na(x)
}

#' Binomial distribution cmf
#'
#' Cumulative mass function for binomial distribution with parameters
#' `size` and `prob`.
#'
#' @param q vector of success counts.
#' @param size number of trials (zero or more).
#' @param prob probability of success on each trial.
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\],
#'   otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param approx logical; if FALSE (the default), aim for bit-exact results,
#'   otherwise use a faster algorithm (which is still more accurate than R
#'   4.6.1 pbinom() and scipy 1.18 stats.binom.cdf()).
#'
#' @details Based on the BFRAC component of Boost's TOMS 708 implementation,
#' and the QD high-precision library.
#'
#' @references DiDonato AR, Morris AH (1992) Algorithm 708: Significant digit
#'   computation of the incomplete beta function ratios.  ACM Transactions on
#'   Mathematical Software (TOMS), 18.
#'
#' @references Hida Y, Li XS, Bailey DH (2001) Algorithms for quad-double
#'   precision floating point arithmetic.  Proceedings of the 15th IEEE
#'   Symposium on Computer Arithmetic.
#'
#' @return cmf(q).
#' @export
pbinom <- function(q, size, prob=0.5, lower.tail=TRUE, log.p=FALSE, approx=FALSE) {
  if (!is.numeric(q)) {
    stop("'q' must be a numeric vector.")
  }
  if (!is.numeric(size) || length(size) != 1) {
    stop("'size' must be a single numeric value.")
  }
  if (!is.numeric(prob) || length(prob) != 1) {
    stop("'prob' must be a single numeric value.")
  }
  if (!is.boolean(lower.tail)) {
    stop("'lower.tail' must be a single boolean value.")
  }
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value.")
  }
  if (!is.boolean(approx)) {
    stop("'approx' must be a single boolean value.")
  }
  result <- pbinom_cpp(q, size, prob, lower.tail, log.p, FALSE, approx, 1.0)
  return(result)
}

#' Binomial distribution ppf
#'
#' Probability point function for binomial distribution with parameters
#' `size` and `prob`.
#'
#' @param p vector of probabilities to find quantiles for.
#' @param size number of trials (zero or more).
#' @param prob probability of success on each trial.
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\],
#'   otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @details Starts with CLT-based initial guess, follows with Newton iteration
#' to find a close-enough pmf value before evaluating cmf.
#'
#' @return ppf(p).
#' @export
qbinom <- function(p, size, prob=0.5, lower.tail=TRUE, log.p=FALSE) {
  if (!is.numeric(p)) {
    stop("'p' must be a numeric vector.")
  }
  if (!is.numeric(size) || length(size) != 1) {
    stop("'size' must be a single numeric value.")
  }
  if (!is.numeric(prob) || length(prob) != 1) {
    stop("'prob' must be a single numeric value.")
  }
  if (!is.boolean(lower.tail)) {
    stop("'lower.tail' must be a single boolean value.")
  }
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value.")
  }
  result <- qbinom_cpp(p, size, prob, lower.tail, log.p)
  return(result)
}

#' Exact Binomial Test
#'
#' Performs an exact test of a simple null hypothesis about the probability of
#' success in a Bernoulli experiment.
#'
#' @param x number of successes, or a vector of length 2 giving the number of
#'   successes and failures, respectively.
#' @param n number of trials; ignored if `x` has length 2.
#' @param p hypothesized probability of success; or numerator if p.denom
#'   specified.
#' @param alternative indicates the alternative hypothesis and must be one of
#'   "`two.sided`", "`greater`", or "`less`".  You can specify just the initial
#'   letter.
#' @param conf.level confidence level for the returned confidence interval.
#' @param midp logical; if TRUE, midp-value is returned in place of p-value.
#' @param p.denom optional denominator of hypothesized success probability.
#'
#' @details Confidence intervals are obtained by a procedure first given in
#'   Clopper and Pearson (1934).  This guarantees that the confidence level is
#'   at least `conf.level`, but in general does not give the shortest-length
#'   confidence intervals.
#'
#'   For the two-sided test, outcome probabilities are only considered to be
#'   tied if they agree to ~40 significant digits.  This behavior allows the
#'   test to remain highly accurate for very large cases.  However, it comes
#'   with a drawback: in small textbook cases where `p` is an exact fraction
#'   subject to rounding error (e.g. 1/3) and ties are not improbable, it is no
#'   longer sufficient to specify `p`=1/3, because that incurs too much
#'   rounding error.  For reliable tie recognition in that scenario, represent
#'   the fraction exactly by specifying e.g. `p`=1 and `p.denom`=3.
#'
#' @references Clopper CJ, Pearson ES (1934) The use of confidence or fiducial
#'   limits illustrated in the case of the binomial.  Biometrika, 26.
#'   <https://doi.org/10.2307/2331986>.
#'
#' @return A list with class "`htest`" containing the following components:
#'   \item{statistic}{the number of successes.}
#'   \item{parameter}{the number of trials.}
#'   \item{p.value}{the p- or midp-value of the test.}
#'   \item{log.p.value}{the log-p- or log-midp-value of the test.}
#'   \item{conf.int}{a confidence interval for the probability of success.}
#'   \item{estimate}{the estimated probability of success.}
#'   \item{null.value}{the probability of success under the null, `p`.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{method}{the character string "`Exact binomial test`" or "`Exact binomial test (midp)`".}
#'   \item{data.name}{a character string giving the names of the data.}
#' @export
binom.test <- function(x, n, p=0.5, alternative=c("two.sided", "less", "greater"), conf.level=0.95, midp=FALSE, p.denom=1) {
  # Wrapper code is mostly identical to (R 4.6.0) binom.test.R since this is
  # designed to be a drop-in replacement.
  #
  # If p.denom is judged to be an inadequate solution to the small-case
  # tie-recognition problem, an alternative is to define e.g. snap.max=1e6 and
  # snap.eps=1e-14, and snap p to "exactly" (up to the internal ~48 digit
  # td_real precision) Fraction.limit_denominator(snap.max) when that
  # corresponds to relative error <= snap.eps.  This would be friendlier to
  # some students, but I don't like it because it distorts results in a
  # nonobvious way across the board.
  DNAME <- deparse1(substitute(x))
  xr <- round(x)

  if (any(is.na(x) | (x < 0)) || max(abs(x-xr)) > 1e-7) {
    stop("'x' must be nonnegative and integer")
  }
  x <- xr
  if (length(x) == 2L) {
    ## x gives successes and failures
    n <- sum(x)
    x <- x[1L]
  } else if (length(x) == 1L) {
    ## x gives successes, n gives trials
    nr <- round(n)
    if ((length(n) > 1L) || is.na(n) || (n < 1) || abs(n-nr) > 1e-7 || (x > nr)) {
      stop("'n' must be a positive intever >= 'x'")
    }
    DNAME <- paste(DNAME, "and", deparse1(substitute(n)))
    n <- nr
  } else {
    stop("incorrect length of 'x'")
  }

  if (length(p) > 1L || is.na(p) || length(p.denom) > 1L || is.na(p.denom) || p.denom == 0 || sign(p) == -sign(p.denom) || abs(p) > abs(p.denom)) {
    stop("'p'/'p.denom' must be a single number between 0 and 1")
  }

  alternative <- match.arg(alternative)

  if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
        (conf.level > 0) && (conf.level < 1))) {
    stop("'conf.level' must be a single number between 0 and 1")
  }

  if (!is.boolean(midp)) {
    stop("'midp' must be a single boolean value.")
  }

  # todo: fill in two.sided stub
  LNPVAL <- switch(alternative,
                   less = pbinom_cpp(x, n, p, TRUE, TRUE, midp, TRUE, p.denom),
                   greater = pbinom_cpp(x - 1, n, p, FALSE, TRUE, midp, TRUE, p.denom),
                   two.sided = binom_test_lnpval(x, n, p, midp, p.denom))
  PVAL <- exp(LNPVAL)

  p.L <- function(x, alpha) {
    if (x == 0) {  # No solution
      0
    } else {
      # todo: check if this is also worth reimplementing within the package
      qbeta(alpha, x, n - x + 1)
    }
  }
  p.U <- function(x, alpha) {
    if (x == n) {
      1
    } else {
      qbeta(1 - alpha, x + 1, n - x)
    }
  }
  CINT <- switch(alternative,
                 less = c(0, p.U(x, 1 - conf.level)),
                 greater = c(p.L(x, 1 - conf.level), 1),
                 two.sided = {
                   alpha <- (1 - conf.level) / 2
                   c(p.L(x, alpha), p.U(x, alpha))
                 })
  attr(CINT, "conf.level") <- conf.level

  ESTIMATE <- x / n

  names(x) <- "number of successes"
  names(n) <- "number of trials"
  names(ESTIMATE) <-
  names(p) <- "probability of success"

  METHOD <- "Exact binomial test"
  if (midp) {
    METHOD <- "Exact binomial test (midp)"
  }

  structure(list(statistic = x,
                 parameter = n,
                 p.value = PVAL,
                 log.p.value = LNPVAL,
                 conf.int = CINT,
                 estimate = ESTIMATE,
                 null.value = p,
                 alternative = alternative,
                 method = METHOD,
                 data.name = DNAME),
            class = "htest")
}

#' Hypergeometric distribution cmf
#'
#' Cumulative mass function for hypergeometric distribution with parameters
#' `m`, `n`, and `k`.
#'
#' @param q vector which can be interpreted as numbers of white balls drawn
#'   without replacement from an urn which contains both black and white balls.
#' @param m the number of white balls in the urn.
#' @param n the number of black balls in the urn.
#' @param k the number of balls drawn from the urn, must be in \[0, m+n\].
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\],
#'   otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param approx logical; if FALSE (the default), aim for bit-exact results,
#'   otherwise use a faster algorithm (which is still more accurate than R
#'   4.6.1 phyper() and scipy 1.18 stats.hypergeom.cdf()).
#'
#' @details Based on the QD high-precision library.
#'
#' @references Hida Y, Li XS, Bailey DH (2001) Algorithms for quad-double
#'   precision floating point arithmetic.  Proceedings of the 15th IEEE
#'   Symposium on Computer Arithmetic.
#'
#' @return cmf(q).
#' @export
phyper <- function(q, m, n, k, lower.tail=TRUE, log.p=FALSE, approx=FALSE) {
  if (!is.numeric(q)) {
    stop("'q' must be a numeric vector.")
  }
  if (!is.numeric(m) || length(m) != 1) {
    stop("'m' must be a single numeric value.")
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("'n' must be a single numeric value.")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("'k' must be a single numeric value.")
  }
  if (!is.boolean(lower.tail)) {
    stop("'lower.tail' must be a single boolean value.")
  }
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value.")
  }
  if (!is.boolean(approx)) {
    stop("'approx' must be a single boolean value.")
  }
  result <- phyper_cpp(q, m, n, k, lower.tail, log.p, approx)
  return(result)
}

#' Hypergeometric distribution ppf
#'
#' Probability point function for hypergeometric distribution with parameters
#' `m`, `n`, and `k`.
#'
#' @param p vector of probabilities to find quantiles for.
#' @param m the number of white balls in the urn.
#' @param n the number of black balls in the urn.
#' @param k the number of balls drawn from the urn, must be in \[0, m+n\].
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\],
#'   otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @details Starts with CLT-based initial guess, follows with Newton iteration
#' to find a close-enough pmf value before evaluating cmf.
#'
#' @return ppf(p).
#' @export
qhyper <- function(p, m, n, k, lower.tail=TRUE, log.p=FALSE) {
  if (!is.numeric(p)) {
    stop("'p' must be a numeric vector.")
  }
  if (!is.numeric(m) || length(m) != 1) {
    stop("'m' must be a single numeric value.")
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("'n' must be a single numeric value.")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("'k' must be a single numeric value.")
  }
  if (!is.boolean(lower.tail)) {
    stop("'lower.tail' must be a single boolean value.")
  }
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value.")
  }
  result <- qhyper_cpp(p, m, n, k, lower.tail, log.p)
  return(result)
}
