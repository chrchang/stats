is.boolean <- function(x) {
  is.logical(x) && length(x) == 1 && !is.na(x)
}

#' Binomial distribution cmf
#'
#' Cumulative mass function for binomial distribution with parameters
#' `size` and `prob`.  Based on the BFRAC component of Boost's TOMS 708
#' implementation, and the QD high-precision library.  The default (`approx` =
#' FALSE) mode aims for bit-exact results.
#'
#' @references DiDonato AR, Morris AH (1992) Algorithm 708: Significant digit
#'   computation of the incomplete beta function ratios.  ACM Transactions on
#'   Mathematical Software (TOMS), 18.
#'
#' @references Hida Y, Li XS, Bailey DH (2001) Algorithms for quad-double
#'   precision floating point arithmetic.  Proceedings of the 15th IEEE
#'   Symposium on Computer Arithmetic.
#'
#' @param q vector of success counts.
#' @param size number of trials (zero or more).
#' @param prob probability of success on each trial.
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\],
#'   otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param approx logical; if TRUE, faster algorithm which doesn't aim for
#'   bit-exact results is used.  (This fast mode is still more accurate than R
#'   4.6.1 pbinom() and scipy 1.18 stats.binom.cdf().)
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
  result <- pbinom_cpp(q, size, prob, lower.tail, log.p, approx)
  return(result)
}

#' Binomial distribution ppf
#'
#' Probability point function for binomial distribution with parameters
#' `size` and `prob`.  Starts with CLT-based initial guess, follows with Newton
#' iteration to find a close-enough pmf value before evaluating cmf.
#'
#' @param p vector of probabilities to find quantiles for.
#' @param size number of trials (zero or more).
#' @param prob probability of success on each trial.
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\],
#'   otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
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

#' Hypergeometric distribution cmf
#'
#' Cumulative mass function for hypergeometric distribution with parameters
#' `m`, `n`, and `k`.  Based on the QD high-precision library.  The default
#' (`approx` = FALSE) mode aims for bit-exact results.
#'
#' @references Hida Y, Li XS, Bailey DH (2001) Algorithms for quad-double
#'   precision floating point arithmetic.  Proceedings of the 15th IEEE
#'   Symposium on Computer Arithmetic.
#'
#' @param q vector which can be interpreted as numbers of white balls drawn
#'   without replacement from an urn which contains both black and white balls.
#' @param m the number of white balls in the urn.
#' @param n the number of black balls in the urn.
#' @param k the number of balls drawn from the urn, must be in \[0, m+n\].
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\],
#'   otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param approx logical; if TRUE, faster algorithm which doesn't aim for
#'   bit-exact results is used.  (This fast mode is still more accurate than R
#'   4.6.1 phyper() and scipy 1.18 stats.hypergeom.cdf().)
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
#' `m`, `n`, and `k`.  Starts with CLT-based initial guess, follows with Newton
#' iteration to find a close-enough pmf value before evaluating cmf.
#'
#' @param p vector of probabilities to find quantiles for.
#' @param m the number of white balls in the urn.
#' @param n the number of black balls in the urn.
#' @param k the number of balls drawn from the urn, must be in \[0, m+n\].
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\],
#'   otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
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
