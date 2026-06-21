#include <math.h>

#include "include/binom.h"
#include "include/fisher.h"
#include "include/hypergeom.h"
#include "include/plink2_base.h"
#include "include/plink2_float.h"
#include "include/plink2_highprec.h"

#include <Rcpp.h>
using namespace Rcpp;

static inline int64_t nonint_warn(double x, double roundx) {
  // Similar to R_D_nonint_check().
  if (fabs(x - roundx) > 1e-9 * MAXV(1, fabs(x))) {
    warning("non-integer x = %.17g", x);
  }
  return roundx;
}

plink2::td_real make_p_tdr(double prob, double prob_denom) {
  plink2::td_real p_tdr = plink2::tdr_make1(prob);
  if (prob_denom != 1.0) {
    p_tdr = plink2::tdr_divd(p_tdr, prob_denom);
  }
  if ((!(p_tdr.x[0] >= 0.0)) || plink2::tdr_gtd(p_tdr, 1.0)) {
    stop("(prob / prob_denom) is not in [0, 1]");
  }
  return p_tdr;
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
  const plink2::td_real p_tdr = make_p_tdr(prob, 1.0);
  const int64_t n = nonint_warn(size, size_round);
  // todo: implement BinomMassMulti() which avoids repeating log(n!), log(p),
  // log(1-p) computations
  const uint32_t x_len = x.size();
  NumericVector results = NumericVector(x_len);
  for (uint32_t idx = 0; idx < x_len; ++idx) {
    const double k_float = x[idx];
    if (isnan(k_float)) {
      results[idx] = k_float;
      continue;
    }
    const double k_round = round(k_float);
    if ((k_round < 0) || (k_round > size_round)) {
      results[idx] = log? (0.0 / 0.0) : 0.0;
      continue;
    }
    const int64_t k = nonint_warn(k_float, k_round);
    results[idx] = plink2::BinomMass(k, n, p_tdr, log);
  }
  return results;
}

//' @title Binomial distribution cmf
//' @description Backend for pbinom()
//' @noRd
// [[Rcpp::export]]
NumericVector pbinom_cpp(NumericVector q, double size, double prob = 0.5, bool lower_tail = true, bool log_p = false, bool approx = false) {
  const double size_round = round(size);
  if ((size_round < 0) || (!(size_round < (1LL << 52)))) {
    stop("size is not in [0, 2^52 - 1]");
  }
  const plink2::td_real p_tdr = make_p_tdr(prob, 1.0);
  const int64_t n = nonint_warn(size, size_round);
  // Unfortunately, can't take advantage of some vectorization opportunities
  // without potentially changing the last bit of some results, so I won't plan
  // on writing Pbinom{,Approx}Multi() functions for now.
  const uint32_t q_len = q.size();
  NumericVector results = NumericVector(q_len);
  for (uint32_t idx = 0; idx < q_len; ++idx) {
    const double k_float = q[idx];
    if (isnan(k_float)) {
      results[idx] = k_float;
      continue;
    }
    const double k_round = round(k_float);
    int64_t k = nonint_warn(k_float, k_round);
    if (k_round < 0) {
      k = -1;
    } else if (k_round > size_round) {
      k = n + 1;
    }
    double result;
    if (approx) {
      result = plink2::PbinomApprox(k, n, p_tdr, !lower_tail, 0, log_p);
    } else {
      result = plink2::Pbinom(k, n, p_tdr, !lower_tail, log_p);
    }
    results[idx] = result;
  }
  return results;
}

// ...more functions coming...
