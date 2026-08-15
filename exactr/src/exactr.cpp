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

static inline bool nonint(double x, double roundx) {
  // R_nonint()
  return (fabs(x - roundx) > 1e-9 * MAXV(1, fabs(x)));
}

static inline bool nonint_warn(double x, double roundx) {
  // R_D_nonint_check()
  if (nonint(x, roundx)) {
    warning("non-integer x = %.17g", x);
    return true;
  }
  return false;
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
NumericVector dbinom(NumericVector x, NumericVector size, NumericVector prob = NumericVector::create(0.5), bool log = false) {
  // Imitate SETUP_Math3 macro in R src/library/stats/src/distn.c .
  const uint32_t x_len = x.size();
  const uint32_t size_len = size.size();
  const uint32_t prob_len = prob.size();
  if ((x_len == 0) || (size_len == 0) || (prob_len == 0)) {
    NumericVector results = NumericVector(0);
    if (x_len == 0) {
      results.attr("dim") = x.attr("dim");
    }
    return results;
  }
  const uint32_t results_len = std::max({x_len, size_len, prob_len});
  NumericVector results = NumericVector(results_len);
  uint32_t nans_produced = 0;

  // Imitate mod_iterate3 macro in R src/library/stats/src/distn.c .
  // possible todo: use parallel-for from RcppParallel for large results_len.
  uint32_t x_idx = 0;
  uint32_t size_idx = 0;
  uint32_t prob_idx = 0;
  for (uint32_t ridx = 0; ridx < results_len; ++ridx) {
    const double k_float = x[x_idx++];
    const double n_float = size[size_idx++];
    const double p = prob[prob_idx++];
    x_idx = (x_idx == x_len)? 0 : x_idx;
    size_idx = (size_idx == size_len)? 0 : size_idx;
    prob_idx = (prob_idx == prob_len)? 0 : prob_idx;
    // if_NA_Math3_set()
    if (traits::is_nan<REALSXP>(k_float) || traits::is_nan<REALSXP>(n_float) || traits::is_nan<REALSXP>(p)) {
      if (NumericVector::is_na(k_float) || NumericVector::is_na(n_float) || NumericVector::is_na(p)) {
        results[ridx] = NA_REAL;
      } else {
        results[ridx] = R_NaN;
      }
      // This counts as a preexisting rather than a 'produced' NaN.
      continue;
    }

    // src/nmath/dbinom.c
    const double n_round = nearbyint(n_float);
    // strangely, R has a test which requires n_float rather than just n_round
    // < 0, when that is not the case for e.g. pbinom().
    if ((p < 0) || (p > 1) || (n_float < 0) || nonint(n_float, n_round)) {
      results[ridx] = R_NaN;
      nans_produced = 1;
      continue;
    }
    const double k_round = nearbyint(k_float);
    if (nonint_warn(k_float, k_round) || (k_float < 0) || (k_float == R_PosInf)) {
      results[ridx] = log? R_NegInf : 0.0;
      continue;
    }
    if ((p == 0.0) || (p == 1.0)) {
      if (((p == 0.0) && (k_round == 0)) || ((k_round == n_round) && (p == 1.0))) {
        results[ridx] = log? 0.0 : 1.0;
      } else {
        results[ridx] = log? R_NegInf : 0.0;
      }
      continue;
    }
    // stats::dbinom() is expected to handle size in [2^52, Inf), and denormal
    // p.  If we need to handle such arguments, I'd prefer to do so by
    // delegating to Rmpfr; that lets us continue to promise <= 1 ULP error.
    // (todo: check whether the 1 ULP error promise is currently broken for p
    // near but not below DBL_MIN, and fix that if so.)
    //
    // In the meantime, plink2::BinomMass()'s domain could be extended up to
    // n=INT64_MAX.
    if (n_round >= (1LL << 52)) {
      if (n_round != R_PosInf) {
        stop("size values in [2^52, Inf) not currently supported");
      }
      results[ridx] = log? R_NegInf : 0.0;
      continue;
    }
    results[ridx] = plink2::BinomMass(static_cast<int64_t>(k_round), static_cast<int64_t>(n_round), plink2::tdr_make1(p), log);
  }
  // Imitate FINISH_Math3 macro in R src/library/stats/src/distn.c .
  if (nans_produced) {
    warning("NaNs produced");
  }
  if (results_len == x_len) {
    results.attr("dim") = x.attr("dim");
  } else if (results_len == size_len) {
    results.attr("dim") = size.attr("dim");
  } else {
    results.attr("dim") = prob.attr("dim");
  }
  return results;
}

//' @title Binomial distribution cdf
//' @description Backend for pbinom(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector pbinom_cpp(NumericVector q, NumericVector size, NumericVector prob, bool lower_tail, bool log_p, bool midp, bool approx, double prob_denom) {
  // Imitate SETUP_Math3 macro in R src/library/stats/src/distn.c .
  const uint32_t q_len = q.size();
  const uint32_t size_len = size.size();
  const uint32_t prob_len = prob.size();
  if ((q_len == 0) || (size_len == 0) || (prob_len == 0)) {
    NumericVector results = NumericVector(0);
    if (q_len == 0) {
      results.attr("dim") = q.attr("dim");
    }
    return results;
  }

  const uint32_t results_len = std::max({q_len, size_len, prob_len});
  NumericVector results = NumericVector(results_len);
  uint32_t nans_produced = 0;

  // Imitate mod_iterate3 macro in R src/library/stats/src/distn.c .
  // possible todo: use parallel-for from RcppParallel for large results_len.
  uint32_t q_idx = 0;
  uint32_t size_idx = 0;
  uint32_t prob_idx = 0;
  for (uint32_t ridx = 0; ridx < results_len; ++ridx) {
    const double k_float = q[q_idx++];
    const double n_float = size[size_idx++];
    const double p = prob[prob_idx++];
    q_idx = (q_idx == q_len)? 0 : q_idx;
    size_idx = (size_idx == size_len)? 0 : size_idx;
    prob_idx = (prob_idx == prob_len)? 0 : prob_idx;
    // if_NA_Math3_set()
    if (traits::is_nan<REALSXP>(k_float) || traits::is_nan<REALSXP>(n_float) || traits::is_nan<REALSXP>(p)) {
      if (NumericVector::is_na(k_float) || NumericVector::is_na(n_float) || NumericVector::is_na(p)) {
        results[ridx] = NA_REAL;
      } else {
        results[ridx] = R_NaN;
      }
      // This counts as a preexisting rather than a 'produced' NaN.
      continue;
    }

    // src/nmath/pbinom.c
    const double n_round = nearbyint(n_float);
    plink2::td_real p_tdr = plink2::tdr_make1(p);
    if (prob_denom != 1.0) {
      p_tdr = plink2::tdr_divd(p_tdr, prob_denom);
    }
    if ((n_float == R_PosInf) || nonint_warn(n_float, n_round) || (n_round < 0) || (!(p_tdr.x[0] >= 0.0)) || plink2::tdr_gtd(p_tdr, 1.0)) {
      results[ridx] = R_NaN;
      nans_produced = 1;
      continue;
    }
    if (k_float < 0) {
      results[ridx] = log_p? R_NegInf : 0.0;
      continue;
    }
    const double k_floor = floor(k_float + 1e-7);
    if (n_round <= k_floor) {
      results[ridx] = log_p? 0.0 : 1.0;
      continue;
    }

    if (n_round >= (1LL << 52)) {
      stop("size values in [2^52, Inf) not currently supported");
    }

    const int64_t k = static_cast<int64_t>(k_floor);
    const int64_t n = static_cast<int64_t>(n_round);
    double result;
    if (approx) {
      result = plink2::PbinomApprox(k, n, p_tdr, !lower_tail, midp, log_p);
    } else {
      result = plink2::Pbinom(k, n, p_tdr, !lower_tail, log_p);
    }
    results[ridx] = result;
  }
  // Imitate FINISH_Math3 macro in R src/library/stats/src/distn.c .
  if (nans_produced) {
    warning("NaNs produced");
  }
  if (results_len == q_len) {
    results.attr("dim") = q.attr("dim");
  } else if (results_len == size_len) {
    results.attr("dim") = size.attr("dim");
  } else {
    results.attr("dim") = prob.attr("dim");
  }
  return results;
}

//' @title Binomial distribution ppf
//' @description Backend for qbinom(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector qbinom_cpp(NumericVector p, NumericVector size, NumericVector prob, bool lower_tail, bool log_p) {
  // Imitate SETUP_Math3 macro in R src/library/stats/src/distn.c .
  const uint32_t p_len = p.size();
  const uint32_t size_len = size.size();
  const uint32_t prob_len = prob.size();
  if ((p_len == 0) || (size_len == 0) || (prob_len == 0)) {
    NumericVector results = NumericVector(0);
    if (p_len == 0) {
      results.attr("dim") = p.attr("dim");
    }
    return results;
  }

  const uint32_t results_len = std::max({p_len, size_len, prob_len});
  NumericVector results = NumericVector(results_len);
  uint32_t nans_produced = 0;

  // Imitate mod_iterate3 macro in R src/library/stats/src/distn.c .
  // possible todo: use parallel-for from RcppParallel for large results_len.
  uint32_t p_idx = 0;
  uint32_t size_idx = 0;
  uint32_t prob_idx = 0;
  for (uint32_t ridx = 0; ridx < results_len; ++ridx) {
    const double p_float = p[p_idx++];
    const double n_float = size[size_idx++];
    const double prob_float = prob[prob_idx++];
    p_idx = (p_idx == p_len)? 0 : p_idx;
    size_idx = (size_idx == size_len)? 0 : size_idx;
    prob_idx = (prob_idx == prob_len)? 0 : prob_idx;
    // if_NA_Math3_set()
    if (traits::is_nan<REALSXP>(p_float) || traits::is_nan<REALSXP>(n_float) || traits::is_nan<REALSXP>(prob_float)) {
      if (NumericVector::is_na(p_float) || NumericVector::is_na(n_float) || NumericVector::is_na(prob_float)) {
        results[ridx] = NA_REAL;
      } else {
        results[ridx] = R_NaN;
      }
      // This counts as a preexisting rather than a 'produced' NaN.
      continue;
    }

    // src/nmath/qbinom.c
    const double n_round = nearbyint(n_float);
    if ((n_float == R_PosInf) || (prob_float < 0.0) || (prob_float > 1.0) || (n_round < 0)) {
      results[ridx] = R_NaN;
      nans_produced = 1;
      continue;
    }
    if (log_p) {
      if (p_float > 0) {
        results[ridx] = R_NaN;
        nans_produced = 1;
        continue;
      }
      if (p_float == 0) {  // upper bound
        results[ridx] = lower_tail? n_round : 0;
        continue;
      }
      if (p_float == R_NegInf) {
        results[ridx] = lower_tail? 0 : n_round;
        continue;
      }
    } else {  // !log_p
      if ((p_float < 0) || (p_float > 1)) {
        results[ridx] = R_NaN;
        nans_produced = 1;
        continue;
      }
      if (p_float == 0) {
        results[ridx] = lower_tail? 0 : n_round;
        continue;
      }
      if (p_float == 1) {
        results[ridx] = lower_tail? n_round : 0;
      }
    }
    if ((prob_float == 0) || (n_round == 0)) {
      results[ridx] = 0;
      continue;
    }
    if (prob_float == 1) {
      results[ridx] = n_round;
      continue;
    }

    if (n_round >= (1LL << 52)) {
      stop("size values in [2^52, Inf) not currently supported");
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
    results[ridx] = plink2::QbinomHalfUlp(p_ddr, static_cast<int64_t>(n_round), plink2::tdr_make1(prob_float), cur_log);
  }
  // Imitate FINISH_Math3 macro in R src/library/stats/src/distn.c .
  if (nans_produced) {
    warning("NaNs produced");
  }
  if (results_len == p_len) {
    results.attr("dim") = p.attr("dim");
  } else if (results_len == size_len) {
    results.attr("dim") = size.attr("dim");
  } else {
    results.attr("dim") = prob.attr("dim");
  }
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
NumericVector dhyper(NumericVector x, NumericVector m, NumericVector n, NumericVector k, bool log = false) {
  // Imitate SETUP_Math4 macro in R src/library/stats/src/distn.c .
  const uint32_t x_len = x.size();
  const uint32_t m_len = m.size();
  const uint32_t n_len = n.size();
  const uint32_t k_len = k.size();
  if ((x_len == 0) || (m_len == 0) || (n_len == 0) || (k_len == 0)) {
    NumericVector results = NumericVector(0);
    if (x_len == 0) {
      results.attr("dim") = x.attr("dim");
    }
    return results;
  }
  const uint32_t results_len = std::max({x_len, m_len, n_len, k_len});
  NumericVector results = NumericVector(results_len);
  uint32_t nans_produced = 0;

  // Imitate mod_iterate4 macro in R src/library/stats/src/distn.c .
  uint32_t x_idx = 0;
  uint32_t m_idx = 0;
  uint32_t n_idx = 0;
  uint32_t k_idx = 0;
  for (uint32_t ridx = 0; ridx < results_len; ++ridx) {
    const double a_float = x[x_idx++];
    const double ac_float = m[m_idx++];
    const double bd_float = n[n_idx++];
    const double ab_float = k[k_idx++];
    x_idx = (x_idx == x_len)? 0 : x_idx;
    m_idx = (m_idx == m_len)? 0 : m_idx;
    n_idx = (n_idx == n_len)? 0 : n_idx;
    k_idx = (k_idx == k_len)? 0 : k_idx;
    // if_NA_Math4_set()
    if (traits::is_nan<REALSXP>(a_float) || traits::is_nan<REALSXP>(ac_float) || traits::is_nan<REALSXP>(bd_float) || traits::is_nan<REALSXP>(ab_float)) {
      if (NumericVector::is_na(a_float) || NumericVector::is_na(ac_float) || NumericVector::is_na(bd_float) || NumericVector::is_na(ab_float)) {
        results[ridx] = NA_REAL;
      } else {
        results[ridx] = R_NaN;
      }
      // This counts as a preexisting rather than a 'produced' NaN.
      continue;
    }

    // src/nmath/dhyper.c
    const double ac_round = nearbyint(ac_float);
    const double bd_round = nearbyint(bd_float);
    const double ab_round = nearbyint(ab_float);
    if ((ac_float < 0) || nonint(ac_float, ac_round) || (bd_float < 0) || nonint(bd_float, bd_round) || (ab_float < 0) || nonint(ab_float, ab_round) || (ab_float > ac_float + bd_float)) {
      results[ridx] = R_NaN;
      nans_produced = 1;
      continue;
    }
    if (a_float < 0) {
      results[ridx] = log? R_NegInf : 0.0;
      continue;
    }
    const double a_round = nearbyint(a_float);
    if (nonint_warn(a_float, a_round)) {
      results[ridx] = R_NaN;
      nans_produced = 1;
      continue;
    }
    if ((ab_round < a_round) || (ac_round < a_round) || (ab_round - a_round > bd_round)) {
      results[ridx] = log? R_NegInf : 0.0;
      continue;
    }
    if (ab_round == 0) {
      if (a_round == 0) {
        results[ridx] = log? 0.0 : 1.0;
      } else {
        results[ridx] = log? R_NegInf : 0.0;
      }
      continue;
    }

    if (ac_round + bd_round >= (1LL << 52)) {
      // stats::dhyper() handles infinities in a fiddly manner which isn't
      // covered by tests, so I won't try to replicate that for now.
      // As with dbinom(), if we need to handle larger finite cases I'd prefer
      // to delegate those to Rmpfr.
      stop("m+n values >= 2^52 not currently supported");
    }

    const int64_t a = static_cast<int64_t>(a_round);
    const int64_t b = static_cast<int64_t>(ab_round) - a;
    const int64_t c = static_cast<int64_t>(ac_round) - a;
    const int64_t d = static_cast<int64_t>(bd_round) - b;
    results[ridx] = plink2::HypergeomMass(a, b, c, d, log);
  }

  // Imitate FINISH_Math4 macro in R src/library/stats/src/distn.c .
  if (nans_produced) {
    warning("NaNs produced");
  }
  if (results_len == x_len) {
    results.attr("dim") = x.attr("dim");
  } else if (results_len == m_len) {
    results.attr("dim") = m.attr("dim");
  } else if (results_len == n_len) {
    results.attr("dim") = n.attr("dim");
  } else {
    results.attr("dim") = k.attr("dim");
  }
  return results;
}

//' @title Hypergeometric distribution cdf
//' @description Backend for phyper(), separated since dots aren't permitted in
//'   C++ parameter names.
//' @noRd
// [[Rcpp::export]]
NumericVector phyper_cpp(NumericVector q, NumericVector m, NumericVector n, NumericVector k, bool lower_tail, bool log_p, bool midp, bool approx) {
  // Imitate SETUP_Math4 macro in R src/library/stats/src/distn.c .
  const uint32_t q_len = q.size();
  const uint32_t m_len = m.size();
  const uint32_t n_len = n.size();
  const uint32_t k_len = k.size();
  if ((q_len == 0) || (m_len == 0) || (n_len == 0) || (k_len == 0)) {
    NumericVector results = NumericVector(0);
    if (q_len == 0) {
      results.attr("dim") = q.attr("dim");
    }
    return results;
  }
  const uint32_t results_len = std::max({q_len, m_len, n_len, k_len});
  NumericVector results = NumericVector(results_len);
  uint32_t nans_produced = 0;

  // Imitate mod_iterate4 macro in R src/library/stats/src/distn.c .
  uint32_t q_idx = 0;
  uint32_t m_idx = 0;
  uint32_t n_idx = 0;
  uint32_t k_idx = 0;
  for (uint32_t ridx = 0; ridx < results_len; ++ridx) {
    const double a_float = q[q_idx++];
    const double ac_float = m[m_idx++];
    const double bd_float = n[n_idx++];
    const double ab_float = k[k_idx++];
    q_idx = (q_idx == q_len)? 0 : q_idx;
    m_idx = (m_idx == m_len)? 0 : m_idx;
    n_idx = (n_idx == n_len)? 0 : n_idx;
    k_idx = (k_idx == k_len)? 0 : k_idx;
    // if_NA_Math4_set()
    if (traits::is_nan<REALSXP>(a_float) || traits::is_nan<REALSXP>(ac_float) || traits::is_nan<REALSXP>(bd_float) || traits::is_nan<REALSXP>(ab_float)) {
      if (NumericVector::is_na(a_float) || NumericVector::is_na(ac_float) || NumericVector::is_na(bd_float) || NumericVector::is_na(ab_float)) {
        results[ridx] = NA_REAL;
      } else {
        results[ridx] = R_NaN;
      }
      // This counts as a preexisting rather than a 'produced' NaN.
      continue;
    }

    // src/nmath/phyper.c boundary checks
    const double ac_round = nearbyint(ac_float);
    const double bd_round = nearbyint(bd_float);
    const double ab_round = nearbyint(ab_float);
    const double total = ac_round + bd_round;
    if ((ac_round < 0) || (bd_round < 0) || (total == R_PosInf) || (ab_round < 0) || (ab_round > total)) {
      results[ridx] = R_NaN;
      nans_produced = 1;
      continue;
    }
    if (total >= (1LL << 52)) {
      stop("m+n values in [2^52, Inf) not currently supported");
    }
    const double a_min = std::max(0.0, ab_round - bd_round);
    const double a_max = std::min(ab_round, ac_round);
    const double a_floor = floor(a_float + 1e-7);
    if ((a_floor < a_min) || (a_floor > a_max)) {
      if ((a_floor < a_min) == lower_tail) {
        results[ridx] = log_p? R_NegInf : 0.0;
      } else {
        results[ridx] = log_p? 0.0 : 1.0;
      }
      continue;
    }
    int64_t a = static_cast<int64_t>(a_floor);
    int64_t b = static_cast<int64_t>(ab_round) - a;
    int64_t c = static_cast<int64_t>(ac_round) - a;
    int64_t d = static_cast<int64_t>(bd_round) - b;
    if (!lower_tail) {
      plink2::swap_i64(&a, &b);
      plink2::swap_i64(&c, &d);
      a -= 1;
      b += 1;
      c += 1;
      d -= 1;
    }
    double result;
    if ((a < 0) || (d < 0)) {
      result = log_p? R_NegInf : 0.0;
    } else if (approx) {
      result = plink2::PhyperApprox(a, b, c, d, 0, midp, log_p);
    } else {
      result = plink2::Phyper(a, b, c, d, log_p);
    }
    results[ridx] = result;
  }

  // Imitate FINISH_Math4 macro in R src/library/stats/src/distn.c .
  if (nans_produced) {
    warning("NaNs produced");
  }
  if (results_len == q_len) {
    results.attr("dim") = q.attr("dim");
  } else if (results_len == m_len) {
    results.attr("dim") = m.attr("dim");
  } else if (results_len == n_len) {
    results.attr("dim") = n.attr("dim");
  } else {
    results.attr("dim") = k.attr("dim");
  }
  return results;
}

// TODO: fix vectorization/recycling in remaining functions.

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
