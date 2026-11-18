// Poisson Exact Test library, copyright (C) 2026 Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "pois.h"

#include <math.h>

#include "plink2_float.h"
#include "plink2_highprec.h"
#include "special_func.h"

#ifdef __cplusplus
namespace plink2 {
#endif

double PoisMass(double k, dd_real lambda_ddr, uint32_t logp) {
  const dd_real ln_prob_ddr = pois_ln_prob_loader(k, lambda_ddr);
  if (logp) {
    return ln_prob_ddr.x[0];
  }
  return ddr_exp(ln_prob_ddr).x[0];
}

// Tail-sum method selection.
//
// When summing adjacent-outcome likelihood ratios outward from distance
// dist := |k - lambda| from the mean, the j-th relative likelihood is roughly
// exp(-(dist * j + j^2 / 2) / lambda); it falls below 2^{-54} at
//   j ~= sqrt(dist^2 + 75 * lambda) - dist.
// Direct summation is used when that is small (or when there are few terms
// left in a lower tail).
//
// Otherwise, the continued fractions are used when dist exceeds
// sqrt(kPoisAsymMaxZ2 * lambda); their iteration count blows up as dist/sqrt(lambda)
// approaches 0.  Temme's expansion (gasym_approx()) is used inside that radius.
CONSTI32(kPoisSimpleMaxIter, 128);
static const double kPoisAsymMaxZ2 = 9.0;

static inline uint32_t pois_simple_sum_ok(double dist, double lambda, double max_terms) {
  const double est_terms = sqrt(dist * dist + 75 * lambda) - dist;
  return (MINV(est_terms, max_terms) <= kPoisSimpleMaxIter);
}

// Returns
//   pmf(k) + pmf(k-1) + ... + pmf(0)
//   --------------------------------
//               pmf(k)
// with 0.5 subtracted when midp is true.  Assumes k < lambda.
static double pois_ltail_lik_simple(double k, double lambda, uint32_t midp) {
  double lik = 1;
  double tail_sum = 1 - 0.5 * midp;
  while (k > 0) {
    lik *= k / lambda;
    k -= 1;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
  }
  return tail_sum;
}

// Returns
//   pmf(k) + pmf(k+1) + pmf(k+2) + ...
//   ----------------------------------
//                pmf(k)
// with 0.5 subtracted when midp is true.  Assumes k >= lambda - 1 (or that
// lambda is small enough that intermediate overflow is impossible).
static double pois_rtail_lik_simple(double k, double lambda, uint32_t midp) {
  double lik = 1;
  double tail_sum = 1 - 0.5 * midp;
  while (1) {
    k += 1;
    lik *= lambda / k;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      return tail_sum;
    }
  }
}

// Same return value as pois_ltail_lik_simple(), computed with Legendre's
// continued fraction for the upper incomplete gamma function:
//   Gamma(a, x) = e^{-x} x^a / (x+1-a + 1(a-1)/(x+3-a + 2(a-2)/(x+5-a + ...)))
// With a = k+1 and x = lambda, the tail-sum is lambda / <denominator>.
// This is the n -> infinity limit of the binomial continued fraction in
// ibeta_continued_fraction_recip_d() when (a, p) <-> (b, q) are swapped.
// Assumes k + 1 < lambda, and is efficient when (lambda - k) is at least a few
// standard deviations.
static double pois_ltail_lik_cf(double k, dd_real lambda_ddr, uint32_t midp) {
  const double aa = k + 1;
  // x + 1 - a
  const double b0 = ddr_subd(lambda_ddr, k).x[0];
  dd_real ff_ddr = ddr_maked(b0);
  double cc = b0;
  double dd = 0.0;
  for (double mm = 1; ; mm += 1) {
    const double cur_a = mm * (aa - mm);
    if (cur_a == 0.0) {
      break;
    }
    const double cur_b = b0 + 2 * mm;
    dd = 1.0 / prefer_fma(cur_a, dd, cur_b);
    cc = cur_b + cur_a / cc;
    const dd_real delta_ddr = ddr_mul2d(cc, dd);
    if (fabs(delta_ddr.x[0] - 1) <= k2m52) {
      break;
    }
    ff_ddr = ddr_mul(ff_ddr, delta_ddr);
  }
  return lambda_ddr.x[0] / ff_ddr.x[0] - 0.5 * midp;
}

// Same return value as pois_rtail_lik_simple(), computed with the
// n -> infinity limit of the Aroian/DiDonato/Morris continued fraction in
// ibeta_continued_fraction_recip_d(): with a = k, x = lambda, and
// lambda' = a - x (the limit of ay - bx), its terms become
//   c_0 = (lambda' + 1) * a / (a+1)
//   cur_a_m = ((a+m-1)/(a+2m-1)) * (x/(a+2m-1)) * x * m
//   cur_b_m = ((a+m)/(a+2m+1)) * (2m + lambda' + 1) + m + m*x/(a+2m-1)
// and the tail-sum is k / <continued fraction>.
// Assumes k > lambda, and is efficient when (k - lambda) is at least a few
// standard deviations.
static double pois_rtail_lik_cf(double k, dd_real lambda_ddr, uint32_t midp) {
  const double aa = k;
  const double xx = lambda_ddr.x[0];
  const double lambda_p1 = ddr_negate(ddr_subd(lambda_ddr, aa)).x[0] + 1.0;
  double cc = (aa / (aa + 1.0)) * lambda_p1;
  dd_real ff_ddr = ddr_maked(cc);
  double dd = 0.0;
  double mm = 1.0;
  while (1) {
    const double denom = aa + 2 * mm - 1;
    const double x_div_denom = xx / denom;
    const double cur_a = ((aa + mm - 1) / denom) * x_div_denom * (xx * mm);
    const double cur_b = prefer_fma((aa + mm) / (aa + 2 * mm + 1), 2 * mm + lambda_p1, prefer_fma(mm, x_div_denom, mm));
    mm += 1.0;
    dd = 1.0 / prefer_fma(cur_a, dd, cur_b);
    cc = cur_b + cur_a / cc;
    const dd_real delta_ddr = ddr_mul2d(cc, dd);
    if (fabs(delta_ddr.x[0] - 1) <= k2m52) {
      break;
    }
    ff_ddr = ddr_mul(ff_ddr, delta_ddr);
  }
  return aa / ff_ddr.x[0] - 0.5 * midp;
}

// Adjusts log(tail-prob) returned by gasym_approx() for complement and midp.
static dd_real pois_gasym_finish(dd_real ln_tail_ddr, uint32_t invert, dd_real lnpmf_ddr, uint32_t midp) {
  if (invert) {
    ln_tail_ddr = ddr_log1p(ddr_negate(ddr_exp(ln_tail_ddr)));
  }
  if (midp) {
    ln_tail_ddr = ddr_logspace_sub(ln_tail_ddr, ddr_add(lnpmf_ddr, _ddr_log05));
  }
  return ln_tail_ddr;
}

// Returns log(P[X <= k] - 0.5 * midp * P[X = k]), given lnpmf_ddr =
// log(P[X = k]).  Assumes k < lambda.
static dd_real pois_ln_ltail(double k, dd_real lambda_ddr, dd_real lnpmf_ddr, uint32_t midp) {
  const double lambda = lambda_ddr.x[0];
  const double dist = ddr_subd(lambda_ddr, k).x[0];
  double tail_lik;
  if (pois_simple_sum_ok(dist, lambda, k + 1)) {
    tail_lik = pois_ltail_lik_simple(k, lambda, midp);
  } else if ((dist * dist < kPoisAsymMaxZ2 * lambda) && (k >= 99)) {
    const double aa = k + 1;
    const uint32_t invert = !ddr_gtd(lambda_ddr, aa);
    return pois_gasym_finish(gasym_approx(aa, lambda_ddr), invert, lnpmf_ddr, midp);
  } else {
    tail_lik = pois_ltail_lik_cf(k, lambda_ddr, midp);
  }
  return ddr_addd(lnpmf_ddr, log(tail_lik));
}

// Returns log(P[X >= k] - 0.5 * midp * P[X = k]), given lnpmf_ddr =
// log(P[X = k]).  Assumes k > lambda.
static dd_real pois_ln_rtail(double k, dd_real lambda_ddr, dd_real lnpmf_ddr, uint32_t midp) {
  const double lambda = lambda_ddr.x[0];
  const double dist = ddr_negate(ddr_subd(lambda_ddr, k)).x[0];
  double tail_lik;
  if (pois_simple_sum_ok(dist, lambda, INFINITY_D)) {
    tail_lik = pois_rtail_lik_simple(k, lambda, midp);
  } else if ((dist * dist < kPoisAsymMaxZ2 * lambda) && (k >= 100)) {
    return pois_gasym_finish(gasym_approx(k, lambda_ddr), 0, lnpmf_ddr, midp);
  } else {
    tail_lik = pois_rtail_lik_cf(k, lambda_ddr, midp);
  }
  return ddr_addd(lnpmf_ddr, log(tail_lik));
}

static inline double pois_finish_small(dd_real ln_small_ddr, uint32_t invert, uint32_t logp) {
  if (!invert) {
    if (ln_small_ddr.x[0] > 0.0) {
      return logp? 0.0 : 1.0;
    }
    return join_log_and_nonlog(ln_small_ddr, 1.0, logp);
  }
  if (ln_small_ddr.x[0] >= 0.0) {
    return logp? -INFINITY_D : 0.0;
  }
  if (logp) {
    return ddr_log1p_negexp(ln_small_ddr);
  }
  return ddr_negate(ddr_subd(ddr_exp(ln_small_ddr), 1.0)).x[0];
}

double PpoisApprox(double obs_k, dd_real lambda_ddr, int32_t complement, int32_t midp, uint32_t logp) {
  if (obs_k < 0) {
    if (complement) {
      return logp? 0.0 : 1.0;
    }
    return logp? -INFINITY_D : 0.0;
  }
  obs_k = floor(obs_k);
  const double lambda = lambda_ddr.x[0];
  if (lambda == 0.0) {
    // Point mass at 0.
    if ((obs_k > 0) || (!midp)) {
      if (complement) {
        return logp? -INFINITY_D : 0.0;
      }
      return logp? 0.0 : 1.0;
    }
    return logp? -kLn2 : 0.5;
  }
  if (obs_k >= S_CAST(double, 1LL << 53)) {
    // Caller is responsible for ensuring lambda < 2^51, in which case
    // P[X > obs_k] is far smaller than DBL_MIN.
    if (complement) {
      return logp? -INFINITY_D : 0.0;
    }
    return logp? 0.0 : 1.0;
  }
  const dd_real lnpmf_ddr = pois_ln_prob_loader(obs_k, lambda_ddr);
  dd_real ln_small_ddr;
  uint32_t small_is_lower;
  if (ddr_gtd(lambda_ddr, obs_k)) {
    ln_small_ddr = pois_ln_ltail(obs_k, lambda_ddr, lnpmf_ddr, midp);
    small_is_lower = 1;
  } else {
    if (midp) {
      ln_small_ddr = pois_ln_rtail(obs_k, lambda_ddr, lnpmf_ddr, 1);
    } else {
      ln_small_ddr = pois_ln_rtail(obs_k + 1, lambda_ddr, pois_ln_prob_loader(obs_k + 1, lambda_ddr), 0);
    }
    small_is_lower = 0;
  }
  return pois_finish_small(ln_small_ddr, (small_is_lower != (!complement)), logp);
}

// Sets *diff_ptr to (approximately) log(pmf(k)) - log(pmf(obs_k)), and
// returns its sign, with near-ties resolved by CompareFactorialProducts().
static intptr_t pois_cmp(int64_t obs_k, int64_t k, dd_real lambda_ddr, dd_real obs_lnpmf_ddr, td_real* starting_lnprobv_tdr_ptr, td_real* ln_lambda_tdr_ptr, double* diff_ptr, dd_real* lnpmf_ddr_ptr) {
  const dd_real lnpmf_ddr = pois_ln_prob_loader(S_CAST(double, k), lambda_ddr);
  *lnpmf_ddr_ptr = lnpmf_ddr;
  const double diff = ddr_sub(lnpmf_ddr, obs_lnpmf_ddr).x[0];
  *diff_ptr = diff;
  // Loader's algorithm with dd_reals should have absolute error well below
  // 2^{-90} times the magnitude of its components.
  const double tol = (1.0 / (k2p50 * (1LL << 35))) * (64 + fabs(lnpmf_ddr.x[0]) + fabs(obs_lnpmf_ddr.x[0]));
  if (diff > tol) {
    return 1;
  }
  if (diff < -tol) {
    return -1;
  }
  // pmf(k) / pmf(obs_k) = lambda^{k - obs_k} * obs_k! / k!
  uint64_t numer_factorial_args[1];
  numer_factorial_args[0] = obs_k;
  uint64_t denom_factorial_args[1];
  denom_factorial_args[0] = k;
  const td_real lambda_tdr = tdr_make(lambda_ddr.x[0], lambda_ddr.x[1], 0.0);
  double dummy;
  return CompareFactorialProducts(1, lambda_tdr, k - obs_k, obs_k, numer_factorial_args, denom_factorial_args, starting_lnprobv_tdr_ptr, ln_lambda_tdr_ptr, &dummy);
}

double PoisTwoSidedP(int64_t obs_k, dd_real lambda_ddr, int32_t midp, uint32_t logp) {
  const double lambda = lambda_ddr.x[0];
  if (lambda == 0.0) {
    if (obs_k) {
      return logp? -INFINITY_D : 0.0;
    }
    if (midp) {
      return logp? -kLn2 : 0.5;
    }
    return logp? 0.0 : 1.0;
  }
  const double xx = obs_k;
  // pmf(x) < pmf(x+1) iff lambda > x+1, and pmf(x) < pmf(x-1) iff lambda < x.
  const uint32_t x_left = ddr_gtd(lambda_ddr, xx + 1);
  const uint32_t x_right = ddr_ltd(lambda_ddr, xx);
  if ((!x_left) && (!x_right)) {
    // obs_k is a mode.
    if (!midp) {
      return logp? 0.0 : 1.0;
    }
    // Other mode with identical probability exists iff lambda is an integer.
    const uint32_t tie = (ddr_is(lambda_ddr, xx) && (obs_k > 0)) || ddr_is(lambda_ddr, xx + 1);
    dd_real ln_sub_ddr = pois_ln_prob_loader(xx, lambda_ddr);
    if (!tie) {
      ln_sub_ddr = ddr_add(ln_sub_ddr, _ddr_log05);
    }
    return pois_finish_small(ln_sub_ddr, 1, logp);
  }
  const dd_real obs_lnpmf_ddr = pois_ln_prob_loader(xx, lambda_ddr);
  const dd_real ln_tail1_ddr = x_left? pois_ln_ltail(xx, lambda_ddr, obs_lnpmf_ddr, midp) : pois_ln_rtail(xx, lambda_ddr, obs_lnpmf_ddr, midp);

  // Now find the furthest-inward outcome on the other side of the mode with
  // probability <= pmf(obs_k).
  td_real starting_lnprobv_tdr = tdr_make1(DBL_MAX);
  td_real ln_lambda_tdr = tdr_make1(DBL_MAX);
  // Largest integer < lambda, and smallest integer > lambda.
  double below_lambda = ceil(lambda) - 1;
  if ((below_lambda + 1 == lambda) && (lambda_ddr.x[1] > 0.0)) {
    below_lambda += 1;
  }
  double above_lambda = floor(lambda) + 1;
  if ((above_lambda - 1 == lambda) && (lambda_ddr.x[1] < 0.0)) {
    above_lambda -= 1;
  }
  double kk = round(2 * lambda - xx);
  double diff;
  dd_real lnpmf_ddr;
  dd_real lnpmf2_ddr;
  intptr_t cmp_result;
  if (x_left) {
    // Find smallest k > lambda with pmf(k) <= pmf(obs_k).  This always
    // exists.  The log-pmf is strictly decreasing over this range, and
    // pmf(above_lambda - 1) > pmf(obs_k).
    if (kk < above_lambda) {
      kk = above_lambda;
    }
    while (1) {
      cmp_result = pois_cmp(obs_k, S_CAST(int64_t, kk), lambda_ddr, obs_lnpmf_ddr, &starting_lnprobv_tdr, &ln_lambda_tdr, &diff, &lnpmf_ddr);
      if (cmp_result > 0) {
        // Successive decrements in log-pmf have increasing magnitude, so this
        // may overshoot.
        const double ll_deriv = log(lambda / (kk + 1));
        double step = ceil(diff / (-ll_deriv));
        if (!(step >= 1)) {
          step = 1;
        }
        kk += step;
        continue;
      }
      if (kk == above_lambda) {
        break;
      }
      cmp_result = pois_cmp(obs_k, S_CAST(int64_t, kk) - 1, lambda_ddr, obs_lnpmf_ddr, &starting_lnprobv_tdr, &ln_lambda_tdr, &diff, &lnpmf2_ddr);
      if (cmp_result > 0) {
        break;
      }
      // pmf(kk-1) <= pmf(obs_k).  Successive increments in log-pmf as we move
      // left have decreasing magnitude, so this never overshoots.
      kk -= 1;
      const double ll_deriv = log(kk / lambda);
      double step = trunc((-diff) / ll_deriv);
      if (step > 0) {
        kk -= step;
        if (kk < above_lambda) {
          kk = above_lambda;
        }
      }
    }
  } else {
    // Find largest k < lambda with pmf(k) <= pmf(obs_k), if it exists.  The
    // log-pmf is strictly increasing over [0, below_lambda], and
    // pmf(below_lambda) > pmf(obs_k).
    if (kk > below_lambda) {
      kk = below_lambda;
    } else if (kk < 0) {
      kk = 0;
    }
    while (1) {
      cmp_result = pois_cmp(obs_k, S_CAST(int64_t, kk), lambda_ddr, obs_lnpmf_ddr, &starting_lnprobv_tdr, &ln_lambda_tdr, &diff, &lnpmf_ddr);
      if (cmp_result > 0) {
        if (kk == 0) {
          // No outcomes on the left side contribute.
          return pois_finish_small(ln_tail1_ddr, 0, logp);
        }
        // May overshoot.
        const double ll_deriv = log(kk / lambda);
        double step = ceil(diff / (-ll_deriv));
        if (!(step >= 1)) {
          step = 1;
        }
        kk -= step;
        if (kk < 0) {
          kk = 0;
        }
        continue;
      }
      if (kk == below_lambda) {
        break;
      }
      cmp_result = pois_cmp(obs_k, S_CAST(int64_t, kk) + 1, lambda_ddr, obs_lnpmf_ddr, &starting_lnprobv_tdr, &ln_lambda_tdr, &diff, &lnpmf2_ddr);
      if (cmp_result > 0) {
        break;
      }
      // Never overshoots.
      kk += 1;
      const double ll_deriv = log(lambda / (kk + 1));
      double step = trunc((-diff) / ll_deriv);
      if (step > 0) {
        kk += step;
        if (kk > below_lambda) {
          kk = below_lambda;
        }
      }
    }
  }
  // Recompute the comparison at the final kk, to determine whether it's an
  // exact tie.  (This is cheap relative to everything else.)
  cmp_result = pois_cmp(obs_k, S_CAST(int64_t, kk), lambda_ddr, obs_lnpmf_ddr, &starting_lnprobv_tdr, &ln_lambda_tdr, &diff, &lnpmf_ddr);
  const uint32_t tail2_midp = midp && (cmp_result == 0);
  const dd_real ln_tail2_ddr = x_left? pois_ln_rtail(kk, lambda_ddr, lnpmf_ddr, tail2_midp) : pois_ln_ltail(kk, lambda_ddr, lnpmf_ddr, tail2_midp);
  return pois_finish_small(ddr_logspace_add(ln_tail1_ddr, ln_tail2_ddr), 0, logp);
}

#ifdef __cplusplus
}
#endif
