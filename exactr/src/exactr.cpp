#include <math.h>

#include "include/binom.h"
#include "include/fisher.h"
#include "include/hypergeom.h"
#include "include/plink2_base.h"
#include "include/plink2_float.h"
#include "include/plink2_highprec.h"

#include <Rcpp.h>
using namespace Rcpp;

static inline bool nonint_warn(double x, double roundx) {
  // Similar to R_D_nonint_check().
  if (fabs(x - roundx) <= 1e-9 * MAXV(1, fabs(x))) {
    return false;
  }
  warning("non-integer x = %.17g", x);
  return true;
}

plink2::td_real make_prob_tdr(double prob, double prob_denom) {
  plink2::td_real prob_tdr = plink2::tdr_make1(prob);
  if (prob_denom != 1.0) {
    prob_tdr = plink2::tdr_divd(prob_tdr, prob_denom);
  }
  if ((!(prob_tdr.x[0] >= 0.0)) || plink2::tdr_gtd(prob_tdr, 1.0)) {
    stop("(prob / prob_denom) is not in [0, 1]");
  }
  return prob_tdr;
}

//' Binomial distribution pmf
//'
//' Mass function for binomial distribution with parameters `size` and `prob`.
//' Implementation is based on log-factorial functions utilizing the QD
//' high-precision library.
//'
//' @references Hida Y, Li XS, Bailey DH (2001) Algorithms for quad-double
//'   precision floating point arithmetic.  Proceedings of the 15th IEEE
//'   Symposium on Computer Arithmetic.
//'
//' @param x vector of success counts.
//' @param size number of trials (zero or more).
//' @param prob probability of success on each trial.
//' @param log logical; if TRUE, probabilities are returned as logarithms.
//' @return pmf(x).
//' @export
// [[Rcpp::export]]
NumericVector dbinom(NumericVector x, double size, double prob = 0.5, bool log = false) {
  const double size_round = round(size);
  if ((size_round < 0) || (!(size_round < (1LL << 52)))) {
    stop("size is not in [0, 2^52 - 1]");
  }
  const plink2::td_real prob_tdr = make_prob_tdr(prob, 1.0);
  const int64_t n = static_cast<int64_t>(size_round);
  // todo: avoid repeating log(n!), log(p), log(1-p) computations when x_len >
  // 1.  BinomMassShared() can compute those values, then BinomMassFinish() can
  // compute the final value given k and those precomputed values.
  // possible todo: use parallel-for from RcppParallel for large x_len.
  // possible todo: when len(x) > max(x) - min(x) + 1 (so there is at least one
  // repeat value), memoize.  (how often does this come up?)
  const uint32_t x_len = x.size();
  NumericVector results = NumericVector(x_len);
  for (uint32_t idx = 0; idx < x_len; ++idx) {
    const double k_float = x[idx];
    if (isnan(k_float)) {
      results[idx] = k_float;
      continue;
    }
    const double k_round = round(k_float);
    if ((k_round < 0) || (k_round > size_round) || nonint_warn(k_float, k_round)) {
      results[idx] = log? (0.0 / 0.0) : 0.0;
      continue;
    }
    results[idx] = plink2::BinomMass(static_cast<int64_t>(k_round), n, prob_tdr, log);
  }
  return results;
}

//' @title Binomial distribution cmf
//' @description Backend for pbinom(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector pbinom_cpp(NumericVector q, double size, double prob = 0.5, bool lower_tail = true, bool log_p = false, bool approx = false) {
  const double size_round = round(size);
  if ((size_round < 0) || (!(size_round < (1LL << 52)))) {
    stop("size is not in [0, 2^52 - 1]");
  }
  const plink2::td_real prob_tdr = make_prob_tdr(prob, 1.0);
  const int64_t n = static_cast<int64_t>(size_round);
  // Unfortunately, can't take advantage of some vectorization opportunities
  // without potentially changing the last bit of some results, so I won't plan
  // on writing Pbinom{,Approx}Multi() functions for now.
  // possible todo: use parallel-for from RcppParallel for large q_len.
  // possible todo: when len(q) > max(q) - min(q) + 1 (so there is at least one
  // repeat value), memoize.  (how often does this come up?)
  const uint32_t q_len = q.size();
  NumericVector results = NumericVector(q_len);
  for (uint32_t idx = 0; idx < q_len; ++idx) {
    const double k_float = q[idx];
    if (isnan(k_float)) {
      results[idx] = k_float;
      continue;
    }
    // Imitate R pbinom().
    const double k_floor = floor(k_float + 1e-7);
    int64_t k;
    // Avoid underflow/overflow in float64 -> int64 conversion.
    if (k_floor < 0) {
      k = -1;
    } else if (k_floor <= size_round) {
      k = static_cast<int64_t> (k_floor);
    } else {
      k = n + 1;
    }
    double result;
    if (approx) {
      result = plink2::PbinomApprox(k, n, prob_tdr, !lower_tail, 0, log_p);
    } else {
      result = plink2::Pbinom(k, n, prob_tdr, !lower_tail, log_p);
    }
    results[idx] = result;
  }
  return results;
}

//' @title Binomial distribution ppf
//' @description Backend for qbinom(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector qbinom_cpp(NumericVector p, double size, double prob = 0.5, bool lower_tail = true, bool log_p = false) {
  const double size_round = round(size);
  if ((size_round < 0) || (!(size_round < (1LL << 52)))) {
    stop("size is not in [0, 2^52 - 1]");
  }
  const plink2::td_real prob_tdr = make_prob_tdr(prob, 1.0);
  const int64_t n = static_cast<int64_t>(size_round);
  // Unfortunately, with current backend implementation, it's tricky to benefit
  // much from vectorization without destabilizing results when a probability
  // is within epsilon of a cdf value, so I won't plan on writing QbinomMulti()
  // functions for now.
  // possible todo: use parallel-for from RcppParallel for large p_len.
  // possible todo: implement a function which precomputes cmf(k_min),
  // cmf(k_min+1), ..., cmf(k_max) in a manner perfectly consistent with
  // Qbinom() (this will probably involve modifying Qbinom() to reliably invert
  // Pbinom()).  Then, when len(p) is large enough, we can just call that
  // function and perform binary searches on its result.
  const uint32_t p_len = p.size();
  NumericVector results = NumericVector(p_len);
  uint32_t nans_produced = 0;
  for (uint32_t idx = 0; idx < p_len; ++idx) {
    const double p_float = p[idx];
    if (isnan(p_float)) {
      nans_produced = 1;
      results[idx] = p_float;
      continue;
    }
    if (log_p) {
      if (p_float > 0.0) {
        nans_produced = 1;
        results[idx] = 0.0 / 0.0;
        continue;
      }
    } else {
      if ((p_float < 0.0) || (p_float > 1.0)) {
        nans_produced = 1;
        results[idx] = 0.0 / 0.0;
      }
    }
    plink2::dd_real p_ddr = plink2::ddr_maked(p_float);
    if (!lower_tail) {
      if (log_p) {
        // ugh
        p_ddr = plink2::ddr_log(plink2::ddr_negate(plink2::ddr_addd(plink2::ddr_exp(p_ddr), -1)));
      } else {
        p_ddr = plink2::ddr_negate(plink2::ddr_add2d(p_float, -1));
      }
    }
    results[idx] = plink2::QbinomHalfUlp(p_ddr, n, prob_tdr, log_p);
  }
  if (nans_produced) {
    warning("NaNs produced");
  }
  return results;
}

// ...more functions coming...
