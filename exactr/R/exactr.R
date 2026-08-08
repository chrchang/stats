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
    stop("'q' must be a numeric vector")
  }
  if (!is.numeric(size) || length(size) != 1) {
    stop("'size' must be a single numeric value")
  }
  if (!is.numeric(prob) || length(prob) != 1) {
    stop("'prob' must be a single numeric value")
  }
  if (!is.boolean(lower.tail)) {
    stop("'lower.tail' must be a single boolean value")
  }
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value")
  }
  if (!is.boolean(approx)) {
    stop("'approx' must be a single boolean value")
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
    stop("'p' must be a numeric vector")
  }
  if (!is.numeric(size) || length(size) != 1) {
    stop("'size' must be a single numeric value")
  }
  if (!is.numeric(prob) || length(prob) != 1) {
    stop("'prob' must be a single numeric value")
  }
  if (!is.boolean(lower.tail)) {
    stop("'lower.tail' must be a single boolean value")
  }
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value")
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
#' @param log.p logical; if TRUE, log.p.value is calculated and returned.
#'   (p.value is still returned.)  log.p.value remains accurate when p.value
#'   underflows to zero.
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
#'   \item{log.p.value}{the log-p- or log-midp-value of the test.  Only present if argument `log.p = TRUE`.}
#'   \item{conf.int}{a confidence interval for the probability of success.}
#'   \item{estimate}{the estimated probability of success.}
#'   \item{null.value}{the probability of success under the null, `p`.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{method}{the character string "`Exact binomial test`" or "`Exact binomial test (midp)`".}
#'   \item{data.name}{a character string giving the names of the data.}
#' @export
binom.test <- function(x, n, p=0.5,
                       alternative=c("two.sided", "less", "greater"),
                       conf.level=0.95, midp=FALSE, log.p=FALSE, p.denom=1) {
  ## This closely follows the GPL-2.0-or-later binom.test.R from R 4.6.0, since
  ## it is designed to be a drop-in replacement.
  ##
  ## If p.denom is judged to be an inadequate solution to the small-case
  ## tie-recognition problem, an alternative is to define e.g. snap.max=1e6 and
  ## snap.eps=1e-14, and snap p to "exactly" (up to the internal ~48 digit
  ## td_real precision) Fraction.limit_denominator(snap.max) when that
  ## corresponds to relative error <= snap.eps.  This would be friendlier to
  ## some students, but I don't like it because it distorts results in a
  ## nonobvious way across the board.
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
    stop("'midp' must be a single boolean value")
  }

  P_OR_LNPVAL <- switch(alternative,
                        less = pbinom_cpp(x, n, p, TRUE, log.p, midp, TRUE, p.denom),
                        greater = pbinom_cpp(x - 1, n, p, FALSE, log.p, midp, TRUE, p.denom),
                        two.sided = binom_test_pval(x, n, p, midp, log.p, p.denom))

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
    METHOD <- paste(METHOD, "(midp)")
  }

  RVAL <- list(statistic = x, parameter = n)
  if (log.p) {
    RVAL <- c(RVAL, p.value = exp(P_OR_LNPVAL), log.p.value = P_OR_LNPVAL)
  } else {
    RVAL <- c(RVAL, p.value = P_OR_LNPVAL)
  }
  structure(c(RVAL,
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
    stop("'q' must be a numeric vector")
  }
  if (!is.numeric(m) || length(m) != 1) {
    stop("'m' must be a single numeric value")
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("'n' must be a single numeric value")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("'k' must be a single numeric value")
  }
  if (!is.boolean(lower.tail)) {
    stop("'lower.tail' must be a single boolean value")
  }
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value")
  }
  if (!is.boolean(approx)) {
    stop("'approx' must be a single boolean value")
  }
  result <- phyper_cpp(q, m, n, k, lower.tail, log.p, FALSE, approx)
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
    stop("'p' must be a numeric vector")
  }
  if (!is.numeric(m) || length(m) != 1) {
    stop("'m' must be a single numeric value")
  }
  if (!is.numeric(n) || length(n) != 1) {
    stop("'n' must be a single numeric value")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("'k' must be a single numeric value")
  }
  if (!is.boolean(lower.tail)) {
    stop("'lower.tail' must be a single boolean value")
  }
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value")
  }
  result <- qhyper_cpp(p, m, n, k, lower.tail, log.p)
  return(result)
}

#' Fisher's Exact Test for Count Data
#'
#' Performs Fisher's exact test for testing the null of independence of rows
#' and columns in a contingency table with fixed marginals.
#'
#' @param x either a two-dimensional contingency table in matrix form, or a
#'   factor object.
#' @param y a factor object; ignored if `x` is a matrix.
#' @param workspace an integer specifying the size of the workspace used in the
#'   network algorithm.  In units of 4 bytes.  Only used for non-simulated
#'   p-values on larger-than-2x3 tables.  This also increases the internal
#'   stack size which allows larger problems to be solved, sometimes needing
#'   hours.  In such cases, `simulate.p.values = TRUE` may be more reasonable.
#' @param hybrid a logical.  Only used for larger-than-2x3 tables, in which
#'   cases it indicates whether the exact probabilities (default) or a hybrid
#'   approximation thereof should be computed.
#' @param hybridPars a numeric vector of length 3, by default describing
#'   "Cochran's conditions" for the validity of the chi-squared approximation,
#'   see 'Details'.
#' @param control a list with named components for low-level algorithm control.
#'   At present the only one used is "`mult`", a positive integer >= 2 with
#'   default 30 used only for larger-than-2x3 tables.  This says how many times
#'   as much space should be allocated to paths as to keys: see file
#'   '`fexact.c`' in the base R stats source code.
#' @param or the hypothesized odds ratio; or numerator if or.denom specified.
#'   Only used in the 2x2 case.
#' @param alternative indicates the alternative hypothesis and must be one of
#'   "`two.sided`", "`greater`" or "`less`".  You can specify just the initial
#'   letter.  Only used in the 2x2 case.
#' @param conf.int logical indicating if a confidence interval for the odds
#'   ratio in a 2x2 table should be computed (and returned).
#' @param conf.level confidence level for the returned confidence interval.
#'   Only used in the 2x2 case and if `conf.int = TRUE`.
#' @param simulate.p.value a logical indicating whether to compute p-values by
#'   Monte Carlo simulation, in larger-than-2x3 tables.
#' @param B an integer specifying the number of replicates used in the Monte
#'   Carlo test when `simulate.p.value` is true.
#' @param midp logical; if TRUE, midp-value is returned in place of p-value.
#'   Not yet supported for larger-than-2x3 tables.
#' @param log.p logical; if TRUE, log.p.value is calculated and returned.
#'   (p.value is still returned.)  log.p.value usually remains accurate when
#'   p.value underflows to zero; a warning is printed when that is not the
#'   case.
#' @param or.denom optional denominator of hypothesized odds ratio.
#'
#' @details If `x` is a matrix, it is taken as a two-dimensional contingency
#'   table, and hence its entries should be nonnegative integers.  Otherwise,
#'   both `x` and `y` must be vectors or factors of the same length.
#'   Incomplete cases are removed, vectors are coerced into factor objects, and
#'   the contingency table is computed from these.
#'
#'   For `or=1` two-sided tests on 2x2 tables, and 2x3 tables, computations are
#'   based on the algorithm described in Chang et al. (2015).  For other tests
#'   on 2x2 tables, p-values are obtained directly using the (central or
#'   non-central) hypergeometric distribution.
#'
#'   For larger tables, computations are based on a C version of the FORTRAN
#'   subroutine `FEXACT` which implements the network developed by Mehta and
#'   Patel (1983, 1986) and improved by Clarkson, Fan and Joe (1993).  The
#'   FORTRAN code can be obtained from <https://netlib.org/toms/643>.  Note the
#'   latter fails (with an error message) when the entries of the table are too
#'   large.  (It transposes the table if necessary so it has no more rows than
#'   columns.  One constraint is that the product of the row marginals be less
#'   than (2^31 - 1).)
#'
#'   For 2x2 tables, the null of conditional independence is equivalent to the
#'   hypothesis that the odds ratio equals one.  'Exact' inference can be based
#"   on observing that in general, given all marginal totals fixed, the first
#'   element of the contingency table has a non-central hypergeometric
#'   distribution with non-centrality parameter given by the odds ratio
#'   (Fisher, 1935).  The alternative for a one-sided test is based on the odds
#'   ratio, so `alternative = greater` is a test of the odds ratio being bigger
#'   than `or`.  Two-sided tests are based on the probabilities of the tables,
#'   and take as "more extreme" all tables with probabilities less than or
#'   equal to that of the observed table, the p-value being the sum of such
#'   probabilities.
#'
#'   For larger-than-2x3 tables and `hybrid = TRUE`, asymptotic chi-squared
#'   probabilities are only used of the 'Cochran conditions' (or modified
#'   version thereof) specified by
#'   `hybridPars = c(expect = 5, percent = 80, Emin = 1)` are satisfied, that
#'   is if no cell has expected counts less than 1 (`= Emin`) and more than 80%
#'   (`= percent`) of the cells have expected counts at least 5 (`= expect`),
#'   otherwise the exact calculation is used.  A corresponding `if()` decision
#'   is made for all sub-tables considered.
#'
#'   In the `r`x`c` case with `r>2` or `c>2`, internal tables can be too
#'   large for the exact test in which case an error is signalled.  Apart from
#'   increasing `workspace` sufficiently, which then may lead to very long
#'   running times, using `simulate.p.value = TRUE` may then often be
#'   sufficient and hence advisable.
#'
#'   Simulation is done conditional on the row and column marginals, and works
#'   only if the marginals are strictly positive.  (A C translation of the
#'   algorithm of Patefield (1981) is used.)  Note that the default number of
#'   replicates (`B = 2000`) implies a minimum p-value of about 0.0005
#'   (1/(`B`+1)).
#'
#' @references Fisher RA (1935) The logic of inductive inference.  Journal of
#'   the Royal Statistical Society Series A, 98.
#'   <https://doi.org/10.2307/2342435>.
#'
#' @references Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ
#'   (2015) Second-generation PLINK: rising to the challenge of larger and
#'   richer datasets.  GigaScience, 4.
#'   <https://doi.org/10.1186/s13742-015-0047-8>.
#'
#' @references Mehta CR and Patel NR (1983) A network algorithm for performing
#'   Fisher's exact test in `r`x`c` contingency tables.  Journal of the
#'   American Statistical Association, 78.
#'   <https://doi.org/10.1080/01621459.1983.10477989>.
#'
#' @references Mehta CR and Patel NR (1986) Algorithm 643: FEXACT, a FORTRAN
#'   subroutine for Fisher's exact test on unordered `r`x`c` contingency
#'   tables.  ACM Transactions on Mathematical Software, 12.
#'   <https://doi.org/10.1145/6497.214326>.
#'
#' @references Clarkson DB, Fan Y, Joe H (1993) A Remark on Algorithm 643:
#'   FEXACT: An Algorithm for Performing Fisher's Exact Test in `r`x`c`
#'   Contingency Tables.  ACM Transactions on Mathematical Software, 19.
#'   <https://doi.org/10.1145/168173.168412>.
#'
#' @references Patefield WM (1981) Algorithm AS 159: An efficient method of
#'   generating `r`x`c` tables with given row and column totals.  Applied
#'   Statistics, 30.  <https://doi.org/10.2307/2346669>.
#'
#' @return A list with class "`htest`" containing the following components:
#'   \item{p.value}{the p- or midp-value of the test.}
#'   \item{log.p.value}{the log-p- or log-midp-value of the test.  Only present if argument `log.p = TRUE`.}
#'   \item{conf.int}{a confidence interval for the odds ratio.  Only present in the 2x2 case and if argument `conf.int = TRUE`.}
#'   \item{estimate}{an estimate of the odds ratio.  Note that the _conditional_ Maximum Likelihood Estimate (MLE) rather than the unconditional MLE (the sample odds ratio) is used.  Only present in the 2x2 case.}
#'   \item{null.value}{the odds ratio under the null, `or`.  Only present in the 2x2 case.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{method}{character string starting with "`Fisher's Exact Test for Count Data`".}
#'   \item{data.name}{a character string giving the name(s) of the data.}
#' @export
fisher.test <- function(x, y = NULL, workspace = 200000, hybrid = FALSE,
                        hybridPars = c(expect = 5, percent = 80, Emin = 1),
                        control = list(), or = 1, alternative = "two.sided",
                        conf.int = TRUE, conf.level = 0.95,
                        simulate.p.value = FALSE, B = 2000, midp = FALSE,
                        log.p = FALSE, or.denom = 1) {
  ## This closely follows parts of the GPL-2.0-or-later fisher.test.R from R
  ## 4.6.0, since it is designed to be a drop-in replacement.
  ##
  ## Start with check for conditions that currently require us to rely on the
  ## stats::fisher.test() fallback.  Specifically, if either:
  ##   nr + nc > 5
  ##   nr == 2, nc == 2, and or != 1
  ## we just forward to stats::fisher.test().

  expr.x <- substitute(x)
  expr.y <- substitute(y)
  ## preserve original x and y until we know we aren't falling back on
  ## stats::fisher.test()
  xmat <- x
  if (is.data.frame(xmat)) {
    xmat <- as.matrix(xmat)
  }
  was.matrix <- is.matrix(xmat)
  if (!was.matrix) {
    if (is.null(y)) {
      stop("if 'x' is not a matrix, 'y' must be given")
    }
    if (length(xmat) != length(y)) {
      stop("'x' and 'y' must have the same length")
    }
    OK <- complete.cases(xmat, y)
    ## use as.factor rather than factor here to be consistent with
    ## pre-tabulated data
    xf <- as.factor(xmat[OK])
    yf <- as.factor(y[OK])
    if ((nlevels(xf) < 2L) || (nlevels(yf) < 2L)) {
      stop("'x' and 'y' must have at least 2 levels")
    }
    xmat <- table(xf, yf)
  }

  nr <- nrow(xmat)
  nc <- ncol(xmat)
  if (!is.boolean(log.p)) {
    stop("'log.p' must be a single boolean value")
  }
  if (nr + nc > 5) {
    if (midp) {
      stop("'midp = TRUE' not yet supported in this case")
    }
    RET <- eval(substitute(stats::fisher.test(ex, ey, workspace, hybrid,
                                              hybridPars, control, or,
                                              alternative, conf.int,
                                              conf.level, simulate.p.value, B),
                           list(ex = expr.x, ey = expr.y)))
    if (log.p) {
      if (RET$p.value <= 0.0) {
        warning("p.value underflowed to zero, and direct calculation of log.p.value not yet supported in this case; setting log.p.value to -Inf")
        RET$log.p.value <- -Inf
      } else {
        RET$log.p.value <- log(RET$p.value)
      }
    }
    return(RET)
  }

  DNAME <- deparse1(expr.x)
  x <- xmat
  METHOD <- "Fisher's Exact Test for Count Data"
  if (was.matrix) {
    if (nr < 2L || nc < 2L) {
      stop("'x' must have at least 2 rows and columns")
    }
    if (!is.numeric(x) || any(x < 0) || anyNA(x)) {
      stop("all entries of 'x' must be nonnegative and finite")
    }
    if (!is.integer(x)) {
      xo <- x
      x <- round(x)
      if (!identical(TRUE, (ax <- all.equal(xo, x)))) {
        warning(gettextf("'x' has been rounded to integer: %s", ax), domain = NA)
      }
      ## R integers are 32-bit; we support "53-bit" integers in the 2x2 case.
      if (nr > 2L || nc > 2L) {
        if (any(x > .Machine$integer.max)) {
          stop("'x' has entries too large to be integer")
        }
        storage.mode(x) <- "integer"
      }
    }
  } else {
    DNAME <- paste(DNAME, "and", deparse1(expr.y))
  }

  have.2x2 <- (nr == 2) && (nc == 2)
  if (have.2x2) {
    alternative <- char.expand(alternative,
                               c("two.sided", "less", "greater"))
    if (length(alternative) > 1L || is.na(alternative)) {
      stop("alternative must be \"two.sided\", \"less\", or \"greater\"")
    }
    if (!((length(conf.level) == 1L) && (conf.level >= 0.5**23) &&
          (!(conf.level > 1 - 1.1920928955078125e-07)))) {
      stop("'conf.level' must be a single number in [2^{-23}, 1 - 2^{-23}]")
    }
    if ((!missing(or) || !missing(or.denom)) && (length(or) > 1L || is.na(or) || length(or.denom) > 1L || is.na(or.denom) || sign(or) != sign(or.denom) || abs(or) >= abs(or.denom) * 2.5822498780869086e+120 || abs(or.denom) >= abs(or) * 2.5822498780869086e+120)) {
      stop("'or'/'or.denom' must be a single number in (2^{-400}, 2^400)")
    }
  }

  if (hybrid) {
    warning("'hybrid' is ignored for 2x2 and 2x3 tables")
  }
  if (midp) {
    METHOD <- paste(METHOD, "(midp)")
  }

  RVAL <- NULL
  if (!have.2x2) {
    ## 2x3 case, otherwise we would have forwarded to stats::fisher.test().
    P_OR_LNPVAL <- fisher23_test_pval(x, midp, log.p)
    if (log.p) {
      RVAL <- list(p.value = exp(P_OR_LNPVAL),
                   log.p.value = P_OR_LNPVAL)
    } else {
      RVAL <- list(p.value = P_OR_LNPVAL)
    }
  } else {
    x11 <- x[1, 1]
    x12 <- x[1, 2]
    x21 <- x[2, 1]
    x22 <- x[2, 2]
    NVAL <- c("odds ratio" = or)

    P_OR_LNPVAL <- NULL
    if (or == 1) {
      P_OR_LNPVAL <- switch(alternative,
                            less = phyper_cpp(x11, x11+x21, x12+x22, x11+x12, TRUE, log.p, midp, TRUE),
                            greater = phyper_cpp(x11-1, x11+x21, x12+x22, x11+x12, FALSE, log.p, midp, TRUE),
                            two.sided = fisher22_2sided_pval(x11, x12, x21, x22, midp, log.p))
    } else {
      P_OR_LNPVAL <- switch(alternative,
                            less = fisher22_1sided_pval_ex(x11, x12, x21, x22, or/or.denom, TRUE, midp, log.p),
                            greater = fisher22_1sided_pval_ex(x11, x12, x21, x22, or/or.denom, FALSE, midp, log.p),
                            two.sided = fisher22_2sided_pval_ex(x11, x12, x21, x22, or, midp, log.p, or.denom))
    }
    ESTIMATE <- c("odds ratio" = odds_ratio_22(x11, x12, x21, x22))
    if (conf.int) {
      CINT <- switch(alternative,
                     less = odds_ratio_ci_22(x11, x12, x21, x22, 0, conf.level),
                     greater = odds_ratio_ci_22(x11, x12, x21, x22, 1 - conf.level, 1),
                     two.sided = odds_ratio_ci_22(x11, x12, x21, x22, (1 - conf.level) / 2, (1 + conf.level) / 2))
      attr(CINT, "conf.level") <- conf.level
    }
    if (log.p) {
      RVAL <- list(p.value = exp(P_OR_LNPVAL),
                   log.p.value = P_OR_LNPVAL)
    } else {
      RVAL <- list(p.value = P_OR_LNPVAL)
    }
    RVAL <- c(RVAL,
              list(conf.int = if (conf.int) CINT,
                   estimate = ESTIMATE,
                   null.value = NVAL))
  }

  structure(c(RVAL,
              alternative = alternative,
              method = METHOD,
              data.name = DNAME),
            class = "htest")
}