#include <math.h>

#include "include/binom.h"
#include "include/binom_detail.h"
#include "include/fisher.h"
#include "include/hypergeom.h"
#include "include/hypergeom_detail.h"
#include "include/nchypergeom_fisher.h"
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
  const uint32_t x_len = x.size();
  NumericVector results = NumericVector(x_len);
  if (size == plink2::INFINITY_D) {
    for (uint32_t idx = 0; idx < x_len; ++idx) {
      const double k_float = x[idx];
      if (isnan(k_float)) {
        results[idx] = k_float;
        continue;
      }
      if ((prob == 0.0) && (round(k_float) == 0)) {
        results[idx] = log? 0.0 : 1.0;
        continue;
      }
      results[idx] = log? (0.0 / 0.0) : 0.0;
    }
  } else {
    const double size_round = round(size);
    // This excludes a bit of stats::dbinom()'s domain.
    if ((size_round < 0) || (!(size_round < (1LL << 52)))) {
      stop("size is not in {[0, 2^52 - 1] U Inf}");
    }
    const plink2::td_real prob_tdr = make_prob_tdr(prob, 1.0);
    const int64_t n = static_cast<int64_t>(size_round);
    uint32_t p_is_half;
    plink2::td_real lfact_n_tdr;
    plink2::td_real lnp_tdr;
    plink2::td_real lnq_tdr;
    plink2::BinomMassMultiKPrecomp(n, prob_tdr, &p_is_half, &lfact_n_tdr, &lnp_tdr, &lnq_tdr);
    // possible todo: use parallel-for from RcppParallel for large x_len.
    // possible todo: when len(x) > max(x) - min(x) + 1 (so there is at least
    // one repeat value), memoize.  (how often does this come up?)
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
      const int64_t k = static_cast<int64_t>(k_round);
      if ((prob == 0.0) || (prob == 1.0)) {
        if (((prob == 0.0) && (k == 0)) || ((k == n) && (prob == 1.0))) {
          results[idx] = log? 0.0 : 1.0;
        } else {
          results[idx] = log? (0.0 / 0.0) : 0.0;
        }
        continue;
      }
      results[idx] = plink2::BinomMassJustK(k, n, p_is_half, lfact_n_tdr, lnp_tdr, lnq_tdr, log);
    }
  }
  results.attr("dim") = x.attr("dim");
  return results;
}

//' @title Binomial distribution cmf
//' @description Backend for pbinom(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector pbinom_cpp(NumericVector q, double size, double prob, bool lower_tail, bool log_p, bool midp, bool approx, double prob_denom) {
  const double size_round = round(size);
  if ((size_round < 0) || (!(size_round < (1LL << 52)))) {
    stop("size is not in [0, 2^52 - 1]");
  }
  const plink2::td_real prob_tdr = make_prob_tdr(prob, prob_denom);
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
      result = plink2::PbinomApprox(k, n, prob_tdr, !lower_tail, midp, log_p);
    } else {
      result = plink2::Pbinom(k, n, prob_tdr, !lower_tail, log_p);
    }
    results[idx] = result;
  }
  results.attr("dim") = q.attr("dim");
  return results;
}

//' @title Binomial distribution ppf
//' @description Backend for qbinom(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector qbinom_cpp(NumericVector p, double size, double prob, bool lower_tail, bool log_p) {
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
    bool cur_log = log_p;
    if (!lower_tail) {
      if (log_p) {
        p_ddr = plink2::ddr_negate(plink2::ddr_addd(plink2::ddr_exp(p_ddr), -1));
        cur_log = false;
      } else {
        p_ddr = plink2::ddr_negate(plink2::ddr_add2d(p_float, -1));
      }
    }
    results[idx] = plink2::QbinomHalfUlp(p_ddr, n, prob_tdr, cur_log);
  }
  if (nans_produced) {
    warning("NaNs produced");
  }
  results.attr("dim") = p.attr("dim");
  return results;
}

//' @title Exact binomial two-sided test p-value
//' @description Implements main p-value calculation for 2-sided binom.test().
//' @noRd
// [[Rcpp::export]]
double binom_2sided_pval(double x_round, double size_round, double prob, bool midp, bool log_p, double prob_denom) {
  // size_round >= 1 checked by caller
  if (!(size_round < (1LL << 52))) {
    stop("size is not in [1, 2^52 - 1]");
  }
  const plink2::td_real prob_tdr = make_prob_tdr(prob, prob_denom);
  const int64_t x = static_cast<int64_t>(x_round);
  const int64_t n = static_cast<int64_t>(size_round);
  return BinomTwoSidedP(x, n, prob_tdr, midp, log_p);
}

//' Hypergeometric distribution pmf
//'
//' Mass function for hypergeometric distribution with parameters `m`, `n`, and
//' `k`.  Implementation is based on log-factorial functions utilizing the QD
//' high-precision library.
//'
//' @references Hida Y, Li XS, Bailey DH (2001) Algorithms for quad-double
//'   precision floating point arithmetic.  Proceedings of the 15th IEEE
//'   Symposium on Computer Arithmetic.
//'
//' @param x vector which can be interpreted as numbers of white balls drawn
//'   without replacement from an urn which contains both black and white
//'   balls.
//' @param m the number of white balls in the urn.
//' @param n the number of black balls in the urn.
//' @param k the number of balls drawn from the urn, must be in \[0, m+n\].
//' @param log logical; if TRUE, probabilities are returned as logarithms.
//' @return pmf(x).
//' @export
// [[Rcpp::export]]
NumericVector dhyper(NumericVector x, double m, double n, double k, bool log = false) {
  const double m_round = round(m);
  const double n_round = round(n);
  const double k_round = round(k);
  if ((m_round < 0) || (!(m_round < (1LL << 52)))) {
    stop("m is not in [0, 2^52 - 1]");
  }
  if ((n_round < 0) || (!(n_round < (1LL << 52)))) {
    stop("n is not in [0, 2^52 - 1]");
  }
  const int64_t mi = static_cast<int64_t>(m_round);
  const int64_t ni = static_cast<int64_t>(n_round);
  const int64_t total = mi + ni;
  if (total >= (1LL << 52)) {
    stop("m+n is not in [0, 2^52 - 1]");
  }

  if ((k_round < 0) || (!(k_round < (1LL << 52)))) {
    stop("k is not in [0, m+n]");
  }
  const int64_t ki = static_cast<int64_t>(k_round);
  if (ki > total) {
    stop("k is not in [0, m+n]");
  }

  const double x_min = MAXV(0, k_round - n_round);
  const double x_max = MINV(m_round, k_round);
  plink2::td_real lfact_m1x_tdr;
  plink2::td_real lfact_m2x_tdr;
  plink2::td_real lfact_mx1_tdr;
  plink2::td_real lfact_mx2_tdr;
  plink2::td_real lfact_mxx_tdr;
  plink2::HypergeomMassMultiKPrecomp(total, ki, mi, &lfact_m1x_tdr, &lfact_m2x_tdr, &lfact_mx1_tdr, &lfact_mx2_tdr, &lfact_mxx_tdr);
  const uint32_t x_len = x.size();
  NumericVector results = NumericVector(x_len);
  for (uint32_t idx = 0; idx < x_len; ++idx) {
    const double x_float = x[idx];
    if (isnan(x_float)) {
      results[idx] = x_float;
      continue;
    }
    const double x_round = round(x_float);
    if ((x_round < x_min) || (!(x_round <= x_max)) || nonint_warn(x_float, x_round)) {
      results[idx] = log? (0.0 / 0.0) : 0.0;
      continue;
    }
    results[idx] = plink2::HypergeomMassJustK(static_cast<int64_t>(x_round), total, ki, mi, lfact_m1x_tdr, lfact_m2x_tdr, lfact_mx1_tdr, lfact_mx2_tdr, lfact_mxx_tdr, log);
  }
  results.attr("dim") = x.attr("dim");
  return results;
}

//' @title Hypergeometric distribution cmf
//' @description Backend for phyper(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector phyper_cpp(NumericVector q, double m, double n, double k, bool lower_tail, bool log_p, bool midp, bool approx) {
  const double m_round = round(m);
  const double n_round = round(n);
  const double k_round = round(k);
  if ((m_round < 0) || (!(m_round < (1LL << 52)))) {
    stop("m is not in [0, 2^52 - 1]");
  }
  if ((n_round < 0) || (!(n_round < (1LL << 52)))) {
    stop("n is not in [0, 2^52 - 1]");
  }
  const int64_t ac = static_cast<int64_t>(m_round);
  const int64_t bd = static_cast<int64_t>(n_round);
  const int64_t total = ac + bd;
  if (total >= (1LL << 52)) {
    stop("m+n is not in [0, 2^52 - 1]");
  }

  if ((k_round < 0) || (!(k_round < (1LL << 52)))) {
    stop("k is not in [0, m+n]");
  }
  const int64_t ab = static_cast<int64_t>(k_round);
  if (ab > total) {
    stop("k is not in [0, m+n]");
  }

  const double a_min = MAXV(0, k_round - n_round);
  const double a_max = MINV(m_round, k_round);

  const uint32_t q_len = q.size();
  NumericVector results = NumericVector(q_len);
  for (uint32_t idx = 0; idx < q_len; ++idx) {
    const double a_float = q[idx];
    if (isnan(a_float)) {
      results[idx] = a_float;
      continue;
    }
    // Imitate R phyper().
    const double a_floor = floor(a_float + 1e-7);
    double result;
    if ((a_floor < a_min) || (a_floor > a_max)) {
      if ((a_floor < a_min) == lower_tail) {
        result = log_p? (0.0 / 0.0) : 0.0;
      } else {
        result = log_p? 0.0 : 1.0;
      }
    } else {
      int64_t a = static_cast<int64_t>(a_floor);
      int64_t b = ab - a;
      int64_t c = ac - a;
      int64_t d = bd - b;
      if (!lower_tail) {
        plink2::swap_i64(&a, &b);
        plink2::swap_i64(&c, &d);
        a -= 1;
        b += 1;
        c += 1;
        d -= 1;
      }
      if ((a < 0) || (d < 0)) {
        result = log_p? (0.0 / 0.0) : 0.0;
      } else if (approx) {
        result = plink2::PhyperApprox(a, b, c, d, 0, midp, log_p);
      } else {
        result = plink2::Phyper(a, b, c, d, log_p);
      }
    }
    results[idx] = result;
  }
  results.attr("dim") = q.attr("dim");
  return results;
}

//' @title Hypergeometric distribution ppf
//' @description Backend for qbinom(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector qhyper_cpp(NumericVector p, double m, double n, double k, bool lower_tail, bool log_p) {
  const double m_round = round(m);
  const double n_round = round(n);
  const double k_round = round(k);
  if ((m_round < 0) || (!(m_round < (1LL << 52)))) {
    stop("m is not in [0, 2^52 - 1]");
  }
  if ((n_round < 0) || (!(n_round < (1LL << 52)))) {
    stop("n is not in [0, 2^52 - 1]");
  }
  const int64_t ac = static_cast<int64_t>(m_round);
  const int64_t bd = static_cast<int64_t>(n_round);
  const int64_t total = ac + bd;
  if (total >= (1LL << 52)) {
    stop("m+n is not in [0, 2^52 - 1]");
  }

  if ((k_round < 0) || (!(k_round < (1LL << 52)))) {
    stop("k is not in [0, m+n]");
  }
  const int64_t ab = static_cast<int64_t>(k_round);
  if (ab > total) {
    stop("k is not in [0, m+n]");
  }

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
    bool cur_log = log_p;
    if (!lower_tail) {
      if (log_p) {
        p_ddr = plink2::ddr_negate(plink2::ddr_addd(plink2::ddr_exp(p_ddr), -1));
        cur_log = false;
      } else {
        p_ddr = plink2::ddr_negate(plink2::ddr_add2d(p_float, -1));
      }
    }
    results[idx] = plink2::QhyperHalfUlp(p_ddr, ac, bd, ab, cur_log);
  }
  if (nans_produced) {
    warning("NaNs produced");
  }
  results.attr("dim") = p.attr("dim");
  return results;
}

//' @title Fisher 2x3 test log-p-value
//' @description Implements main p-value calculation for 2x3 tables.
//' @noRd
// [[Rcpp::export]]
double fisher23_test_pval(IntegerMatrix x, bool midp, bool log_p) {
  const int64_t m11 = x(0, 0);
  int64_t m12;
  int64_t m13;
  int64_t m21;
  const int64_t m22 = x(1, 1);
  int64_t m23;
  if (x.nrow() == 2) {
    m12 = x(0, 1);
    m13 = x(0, 2);
    m21 = x(1, 0);
    m23 = x(1, 2);
  } else {
    m12 = x(1, 0);
    m13 = x(2, 0);
    m21 = x(0, 1);
    m23 = x(2, 1);
  }
  if (m11 + m12 + m13 + m21 + m22 + m23 >= (1LL << 31)) {
    stop("problem instance too large (2x3 table entries must sum to less than 2^31)");
  }
  const double ln_pval = plink2::Fisher23LnP(m11, m12, m13, m21, m22, m23, midp);
  if (log_p) {
    return ln_pval;
  }
  return exp(ln_pval);
}

//' @title Fisher 2x2 two-sided test p-value
//' @description Implements main two-sided p-value calculation for 2x2 tables.
//' @noRd
// [[Rcpp::export]]
double fisher22_2sided_pval(double x11, double x12, double x21, double x22, bool midp, bool log_p) {
  // Assumes {x11,x12,x21,x22} are nonnegative integers.
  if (x11 + x12 + x21 + x22 >= (1LL << 52)) {
    stop("problem instance too large (2x2 table entries must sum to less than 2^52)");
  }
  return plink2::Fisher22TwoSidedP(static_cast<int64_t>(x11), static_cast<int64_t>(x12), static_cast<int64_t>(x21), static_cast<int64_t>(x22), midp, log_p);
}

//' @title Fisher 2x2 one-sided test p-value, nonstandard odds ratio
//' @description This is equivalent to evaluation of the Fisher's noncentral
//'   hypergeometric distribution cdf.
//' @noRd
// [[Rcpp::export]]
double fisher22_1sided_pval_ex(double x11, double x12, double x21, double x22, double odds, bool lower_tail, bool midp, bool log_p) {
  // Assumes {x11,x12,x21,x22} are nonnegative integers, and odds has already
  // been verified to be in (2^{-400}, 2^400).
  if (x11 + x12 + x21 + x22 >= (1LL << 52)) {
    stop("problem instance too large (2x2 table entries must sum to less than 2^52)");
  }
  return plink2::P_FNCHypergeo(static_cast<int64_t>(x11), static_cast<int64_t>(x12), static_cast<int64_t>(x21), static_cast<int64_t>(x22), odds, !lower_tail, midp, log_p);
}

//' @title Fisher 2x2 two-sided test p-value, nonstandard odds ratio
//' @description Implements main two-sided p-value calculation for 2x2 tables
//'   with hypothesized odds!=1.
//' @noRd
// [[Rcpp::export]]
double fisher22_2sided_pval_ex(double x11, double x12, double x21, double x22, double odds, bool midp, bool log_p, double odds_denom) {
  if (x11 + x12 + x21 + x22 >= (1LL << 52)) {
    stop("problem instance too large (2x2 table entries must sum to less than 2^52)");
  }
  plink2::td_real odds_tdr = plink2::tdr_make1(odds);
  if (odds_denom != 1.0) {
    odds_tdr = plink2::tdr_divd(odds_tdr, odds_denom);
  }
  return plink2::Fisher22TwoSidedPEx(static_cast<int64_t>(x11), static_cast<int64_t>(x12), static_cast<int64_t>(x21), static_cast<int64_t>(x22), odds_tdr, midp, log_p);
}

//' @title Odds ratio, 2x2 table
//' @description Implements odds ratio point-estimate for 2x2 tables.
//' @noRd
// [[Rcpp::export]]
double odds_ratio_22(double x11, double x12, double x21, double x22) {
  // Assumes parameters have been validated.
  return plink2::Fisher22OddsRatio(static_cast<int64_t>(x11), static_cast<int64_t>(x12), static_cast<int64_t>(x21), static_cast<int64_t>(x22));
}

//' @title Odds ratio confidence interval, 2x2 table
//' @description Implements odds ratio confidence-interval calculation for 2x2
//'   tables.
//' @noRd
// [[Rcpp::export]]
NumericVector odds_ratio_ci_22(double x11, double x12, double x21, double x22, double low, double high) {
  // Assumes parameters have been validated to be in-range.  (In particular, we
  // currently enforce conf.level in [2^{-23}, 1 - 2^{-23}] to stay within
  // Fisher22OddsRatioCI()'s domain.)
  NumericVector result(2);
  double low_result;
  double high_result;
  plink2::Fisher22OddsRatioCI(static_cast<int64_t>(x11), static_cast<int64_t>(x12), static_cast<int64_t>(x21), static_cast<int64_t>(x22), low, high, &low_result, &high_result);
  result[0] = low_result;
  result[1] = high_result;
  return result;
}
