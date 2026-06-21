#' Binomial distribution cmf
#'
#' Cumulative mass function for binomial distribution with parameters
#' `size` and `prob`.  Based on the BFRAC component of Boost's TOMS 708
#' implementation, and the [QD high-precision
#' library](https://github.com/BL-highprecision/QD).
#'
#' @references https://dl.acm.org/doi/10.1145/131766.131776
#'
#' @param q vector of success counts.
#' @param size number of trials (zero or more).
#' @param prob probability of success on each trial.
#' @param lower.tail logical; if TRUE, probabilities are P\[X <= x\], otherwise, P\[X > x\].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param approx logical; if TRUE, faster algorithm which doesn't aim for bit-exact results is used
#' @return cmf(q).
#' @export
pbinom <- function(q, size, prob=0.5, lower.tail=FALSE, log.p=FALSE, approx=FALSE) {
  if (!is.numeric(q)) {
    stop("'q' must be a numeric vector.")
  }
  if (!is.numeric(size) || length(size) != 1) {
    stop("'size' must be a single numeric value.")
  }
  if (!is.numeric(prob) || length(prob) != 1) {
    stop("'prob' must be a single numeric value.")
  }
  result <- pbinom_cpp(q, size, prob, lower.tail, log.p, approx)
  return(result)
}