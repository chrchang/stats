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
#'   4.6.0 pbinom() and scipy 1.18 stats.binom.cdf().)
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