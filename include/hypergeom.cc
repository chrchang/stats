// Fisher's Exact Test library, copyright (C) 2013-2026 Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include <math.h>

#include "hypergeom_detail.h"
#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// obs_m11 + obs_m12 + obs_m21 + obs_m22 assumed to be <2^52.  (I'll keep the
// old Fisher's 2x2 exact test variable names since they're relatively
// self-explanatory, and scipy and R don't agree on hypergeometric-distribution
// parameterization.)
// See Phyper() below for a higher-accuracy variant of this function; this one
// mostly avoids use of dd_reals.
// TODO: learn more about continued fractions, and investigate whether
// something like the TOMS 708 algorithm can help here.
double PhyperApprox(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp, uint32_t logp) {
  // Normalize.
  if (obs_m11 < obs_m22) {
    swap_i64(&obs_m11, &obs_m22);
  }
  if (obs_m12 < obs_m21) {
    swap_i64(&obs_m12, &obs_m21);
  }
  // Flipping m11<->m12 and m21<->m22 also flips the direction of the
  // alternative hypothesis.  So we flip on m11-is-greater alternative
  // hypothesis here to allow the rest of the code to assume m11-is-less.
  if (m11_is_greater_alt) {
    swap_i64(&obs_m11, &obs_m12);
    swap_i64(&obs_m21, &obs_m22);
  }
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;
  if (m11 * m22 >= m12 * m21) {
    // We're at or to the right of the mode.
    // Start by computing an upper bound on the right-sum, and then iterating
    // leftward until we either know the p-value > 1 - logp? DBL_MIN : 2^{-54}
    // (at which point we just return log(1) or 1; in the logp case, don't want
    // to risk imposing a surprising denormal-handling performance penalty for
    // no good reason), or remaining left likelihoods are smaller than the
    // prevision limit.
    const double first_right_mult = m12 * m21 / ((m11 + 1) * (m22 + 1));
    // r + r^2 + ... = r / (1-r)
    const double right_upper_bound = 0.5 * midp + first_right_mult / (1 - first_right_mult);
    if (right_upper_bound == 0.0) {
      // p-value is exactly 1 when m12*m21==0 and midp is false
      return logp? 0 : 1;
    }

    // Scale our starting likelihood so that we overflow to INFINITY when we'd
    // want to early-exit and return log(1) or 1; this saves us a comparison in
    // the loop.
    const double start_lik = (DBL_MAX * (logp? DBL_MIN : k2m54)) / right_upper_bound;
    double lik = start_lik;
    double left_sum = start_lik;
    while (1) {
      m12 += 1;
      m21 += 1;
      lik *= m11 * m22 / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      const double preadd = left_sum;
      left_sum += lik;
      if (left_sum == preadd) {
        break;
      }
    }
    if (left_sum == INFINITY_D) {
      return logp? 0 : 1;
    }

    // Now compute the right-sum to the precision limit.
    double right_sum = first_right_mult * start_lik;
    m11 = obs_m11 + 1.0;
    m12 = obs_m12 - 1;
    m21 = obs_m21 - 1;
    m22 = obs_m22 + 1;
    lik = right_sum;
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= m12 * m21 / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      const double preadd = right_sum;
      right_sum += lik;
      if (right_sum == preadd) {
        break;
      }
    }
    // For one-sided test, slightly more convenient to exclude midp term from
    // left_sum and right_sum since it just cancels out in denom
    const double midp_numer = -0.5 * midp * start_lik;
    const double denom = right_sum + left_sum;
    if (!logp) {
      return (left_sum + midp_numer) / denom;
    }
    return log1p((midp_numer - right_sum) / denom);
  }
  // We're to the left of the mode, and are responsible for tiny p-values.
  // If we're close enough to the mode that a simple left_sum / (left_sum +
  // right_sum) calculation doesn't risk overflow with the initial
  // relative-likelihood set to 1, just do that.
  // Otherwise, we evaluate the starting log-likelihood with dd_reals to work
  // around catastrophic cancellation, and then iterate leftward to the
  // precision limit.
  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx2 = obs_m12 + obs_m22;
  const double mxx = m1x + m2x;
  const double modal_m22 = m2x * mx2 / mxx;
  // ((168^168) / 168!)^4 ~= 6.3e285
  // Need a bit more headroom than int32_t case.
  if (modal_m22 <= obs_m22 + 168) {
    double lik = 1;
    double right_sum = 0;
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= m12 * m21 / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      const double preadd = right_sum;
      right_sum += lik;
      if (right_sum == preadd) {
        break;
      }
    }
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21;
    m22 = obs_m22;
    lik = 1;
    double left_sum = 1;
    while (1) {
      m12 += 1;
      m21 += 1;
      lik *= m11 * m22 / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      const double preadd = left_sum;
      left_sum += lik;
      if (left_sum == preadd) {
        break;
      }
    }
    const double pval = (left_sum - 0.5 * midp) / (left_sum + right_sum);
    return logp? log(pval) : pval;
  }
  const dd_real starting_lnprob_ddr = hypergeom_ln_prob_internal(obs_m11, obs_m12, obs_m21, obs_m22);
  // left_sum is the sum of < 2^52 terms, each of which is <= 1, so if
  // starting_lnprob < DBL_MIN / 2^52, final return value should always be 0
  // when logp=false and we're flushing denormals to zero.  DBL_MIN is
  // 2^{-1022}.
  //
  // 2^{-1074} is the smallest positive denormal, and (1 + epsilon) * 2^{-1075}
  // is the smallest number that should be rounded up to it, so -1074 can be
  // replaced with -1127 if we want this function to return denormals.
  //
  // (Yes, a tighter bound could be established for left_sum if it matters.)
  if ((!logp) && (starting_lnprob_ddr.x[0] < -1074 * kLn2)) {
    return 0;
  }
  double lik = 1;
  double left_sum = 1 - 0.5 * midp;
  while (1) {
    m12 += 1;
    m21 += 1;
    lik *= m11 * m22 / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    const double preadd = left_sum;
    left_sum += lik;
    if (left_sum == preadd) {
      break;
    }
  }
  return join_log_and_nonlog(starting_lnprob_ddr, left_sum, logp);
}

// Assumes parameters are nonnegative, and add up to less than 2^52.  Aims for
// <1 ULP relative error.
double Phyper(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, uint32_t logp) {
  // Normalize.
  if (obs_m11 < obs_m22) {
    swap_i64(&obs_m11, &obs_m22);
  }
  if (obs_m12 < obs_m21) {
    swap_i64(&obs_m12, &obs_m21);
  }
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;
  if (m11 * m22 >= m12 * m21) {
    // We're at or to the right of the mode.  Don't have to worry about cmf
    // values < DBL_MIN, so even for e.g. sum > 2^51 our relative error is not
    // inflated by a need to work in log-space.  (todo: we should still replace
    // calculation of left_sum with a log-factorial-based likelihood evaluation
    // when that would be faster and not blow our error budget.  But let's get
    // a minimal implementation working first.)
    const dd_real first_right_mult_ddr = ddr_accurate_div(ddr_mul2d(m12, m21), ddr_mul2d(m11 + 1, m22 + 1));
    if (ddr_is_zero(first_right_mult_ddr)) {
      return logp? 0.0 : 1.0;
    }
    const double right_upper_bound = first_right_mult_ddr.x[0] / ddr_negate(ddr_subd(first_right_mult_ddr, 1.0)).x[0];

    const double start_lik = (DBL_MAX * (logp? DBL_MIN : k2m54)) / right_upper_bound;
    // We want to compute left_sum_ddr to at most 2^{-67} relative error.  See
    // related discussion in Pbinom() implementation.
    const double min_incr_left = (1.0 / (1 << 15)) / (m22 * m22);
    dd_real lik_ddr = ddr_maked(start_lik);
    dd_real left_sum_ddr = lik_ddr;
    do {
      m12 += 1;
      m21 += 1;
      lik_ddr = ddr_mul(lik_ddr, ddr_accurate_div(ddr_mul2d(m11, m22), ddr_mul2d(m12, m21)));
      m11 -= 1;
      m22 -= 1;
      left_sum_ddr = ddr_add(left_sum_ddr, lik_ddr);
    } while (lik_ddr.x[0] > left_sum_ddr.x[0] * min_incr_left);
    if (!(left_sum_ddr.x[0] < INFINITY_D)) {
      return logp? 0.0 : 1.0;
    }
    if (m22 > 0) {
      // Continue the calculation with ordinary precision.
      double lik = lik_ddr.x[0];
      double left_tail_sum = 0.0;
      while (1) {
        m12 += 1;
        m21 += 1;
        lik *= m11 * m22 / (m12 * m21);
        m11 -= 1;
        m22 -= 1;
        const double preadd = left_tail_sum;
        left_tail_sum += lik;
        if (left_tail_sum == preadd) {
          break;
        }
      }
      left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
      if (left_sum_ddr.x[0] == INFINITY_D) {
        return logp? 0.0 : 1.0;
      }
    }

    // Now compute the right-sum to at most 2^{-67} relative error.
    dd_real right_sum_ddr = ddr_muld(first_right_mult_ddr, start_lik);
    m11 = obs_m11 + 1;
    m12 = obs_m12 - 1;
    m21 = obs_m21 - 1;
    m22 = obs_m22 + 1;
    lik_ddr = right_sum_ddr;
    if (m21 > 0) {
      const double min_incr_right = (1.0 / (1 << 15)) / (m21 * m21);
      do {
        m11 += 1;
        m22 += 1;
        lik_ddr = ddr_mul(lik_ddr, ddr_accurate_div(ddr_mul2d(m12, m21), ddr_mul2d(m11, m22)));
        m12 -= 1;
        m21 -= 1;
        right_sum_ddr = ddr_add(right_sum_ddr, lik_ddr);
      } while (lik_ddr.x[0] > right_sum_ddr.x[0] * min_incr_right);
      if (m21 > 0) {
        // Continue the calculation with ordinary precision.
        double lik = lik_ddr.x[0];
        double right_tail_sum = 0.0;
        while (1) {
          m11 += 1;
          m22 += 1;
          lik *= m12 * m21 / (m11 * m22);
          m12 -= 1;
          m21 -= 1;
          const double preadd = right_tail_sum;
          right_tail_sum += lik;
          if (right_tail_sum == preadd) {
            break;
          }
        }
        right_sum_ddr = ddr_addd(right_sum_ddr, right_tail_sum);
      }
    }
    const dd_real denom_ddr = ddr_add(left_sum_ddr, right_sum_ddr);
    if (!(denom_ddr.x[0] < INFINITY_D)) {
      return logp? 0.0 : 1.0;
    }
    const dd_real one_minus_prob_ddr = ddr_accurate_div(right_sum_ddr, denom_ddr);
    if (!logp) {
      return -ddr_subd(one_minus_prob_ddr, 1.0).x[0];
    }
    return ddr_log1p(ddr_negate(one_minus_prob_ddr)).x[0];
  }
  // We're at or to the left of the mode, and are responsible for tiny cmf
  // values.
  // If we're close enough to the mode that a simple left_sum / (left_sum +
  // right_sum) calculation doesn't risk overflow with the initial
  // relative-likelihood set to 1, just do that.
  // Otherwise, we use ddr_lfact() and friends to compute the starting
  // log-likelihood, and accumulate the tail-sum from there.
  // (todo: opportunistically use the _lfact() path when that's within the
  // error budget and rates to be faster than the left_sum / (left_sum +
  // right_sum) approach.)
  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx2 = obs_m12 + obs_m22;
  const double mxx = m1x + m2x;
  const double modal_m22 = m2x * mx2 / mxx;
  if (modal_m22 <= m22 + 168) {
    dd_real lik_ddr = ddr_maked(1.0);
    dd_real right_sum_ddr = ddr_maked(0.0);
    if (m21 > 0) {
      const double min_incr_right = (1.0 / (1 << 15)) / (m21 * m21);
      do {
        m11 += 1;
        m22 += 1;
        lik_ddr = ddr_mul(lik_ddr, ddr_accurate_div(ddr_mul2d(m12, m21), ddr_mul2d(m11, m22)));
        m12 -= 1;
        m21 -= 1;
        right_sum_ddr = ddr_add(right_sum_ddr, lik_ddr);
      } while (lik_ddr.x[0] > right_sum_ddr.x[0] * min_incr_right);
      if (m21 > 0) {
        double lik = lik_ddr.x[0];
        double right_tail_sum = 0.0;
        while (1) {
          m11 += 1;
          m22 += 1;
          lik *= m12 * m21 / (m11 * m22);
          m12 -= 1;
          m21 -= 1;
          const double preadd = right_tail_sum;
          right_tail_sum += lik;
          if (right_tail_sum == preadd) {
            break;
          }
        }
        right_sum_ddr = ddr_addd(right_sum_ddr, right_tail_sum);
      }
    }
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21;
    m22 = obs_m22;
    lik_ddr = ddr_maked(1.0);
    dd_real left_sum_ddr = lik_ddr;
    if (m22 > 0) {
      const double min_incr_left = (1.0 / (1 << 15)) / (m22 * m22);
      do {
        m12 += 1;
        m21 += 1;
        lik_ddr = ddr_mul(lik_ddr, ddr_accurate_div(ddr_mul2d(m11, m22), ddr_mul2d(m12, m21)));
        m11 -= 1;
        m22 -= 1;
        left_sum_ddr = ddr_add(left_sum_ddr, lik_ddr);
      } while (lik_ddr.x[0] > left_sum_ddr.x[0] * min_incr_left);
      if (m22 > 0) {
        double lik = lik_ddr.x[0];
        double left_tail_sum = 0.0;
        while (1) {
          m12 += 1;
          m21 += 1;
          lik *= m11 * m22 / (m12 * m21);
          m11 -= 1;
          m22 -= 1;
          const double preadd = left_tail_sum;
          left_tail_sum += lik;
          if (left_tail_sum == preadd) {
            break;
          }
        }
        left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
      }
    }
    const dd_real prob_ddr = ddr_accurate_div(left_sum_ddr, ddr_add(left_sum_ddr, right_sum_ddr));
    if (!logp) {
      return prob_ddr.x[0];
    }
    return ddr_log(prob_ddr).x[0];
  }
  dd_real ln_prob_ddr = hypergeom_ln_prob_internal(obs_m11, obs_m12, obs_m21, obs_m22);
  dd_real lik_ddr = ddr_maked(1.0);
  dd_real left_sum_ddr = lik_ddr;
  if (m22 > 0) {
    const double min_incr_left = (1.0 / (1 << 15)) / (m22 * m22);
    do {
      m12 += 1;
      m21 += 1;
      lik_ddr = ddr_mul(lik_ddr, ddr_accurate_div(ddr_mul2d(m11, m22), ddr_mul2d(m12, m21)));
      m11 -= 1;
      m22 -= 1;
      left_sum_ddr = ddr_add(left_sum_ddr, lik_ddr);
    } while (lik_ddr.x[0] > left_sum_ddr.x[0] * min_incr_left);
    if (m22 > 0) {
      double lik = lik_ddr.x[0];
      double left_tail_sum = 0.0;
      while (1) {
        m12 += 1;
        m21 += 1;
        lik *= m11 * m22 / (m12 * m21);
        m11 -= 1;
        m22 -= 1;
        const double preadd = left_tail_sum;
        left_tail_sum += lik;
        if (left_tail_sum == preadd) {
          break;
        }
      }
      left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
    }
  }
  ln_prob_ddr = ddr_add(ln_prob_ddr, ddr_log(left_sum_ddr));
  if (logp) {
    return ln_prob_ddr.x[0];
  }
  return ddr_exp(ln_prob_ddr).x[0];
}

// Returns smallest 'a' in the distribution support for which cdf(a) >= p.
// Assumes 0 <= a+b+c+d < 2^52, and should achieve p relative error < 2^{-54}
// unless a+b+c+d is well over 2^31.  Since fisher22 parameterization doesn't
// work here, we mirror R qhyper()'s parameters.
//
// The QhyperHalfUlp() entry point subtracts the natural epsilon value (0.5
// times the value of the least significant bit in p_ddr.x[0]) off of p_ddr,
// for the goal described above, before calling Qhyper().
//
// Probable todo: (except possibly on some very large cases) start with faster
// lower-accuracy interval-math calculation, and fall back on reliable
// high-accuracy calculation only when needed.
int64_t Qhyper(dd_real p_or_lnp_ddr, int64_t ac, int64_t bd, int64_t ab, uint32_t logp) {
  const int64_t abcd = ac + bd;
  // Normalize so that d's support is of the form [0, max_d], to minimize
  // differences from the earlier functions in this file.  We perform most
  // calculations in terms of d, and then add a_minus_d to the return value.
  int64_t a_minus_d = ab - bd;
  int64_t final_return_incr = a_minus_d;
  if (a_minus_d < 0) {
    a_minus_d = -a_minus_d;
    final_return_incr = 0;
    ab = abcd - ab;
    swap_i64(&ac, &bd);
  }
  if (ab < ac) {
    swap_i64(&ab, &ac);
    bd = abcd - ac;
  }
  const int64_t max_d = abcd - ab;
  if ((ddr_is_zero(p_or_lnp_ddr) && (!logp)) || (max_d == 0)) {
    return final_return_incr;
  }
  if ((ddr_is(p_or_lnp_ddr, 1) && (!logp)) || (ddr_is_zero(p_or_lnp_ddr) && logp)) {
    return max_d + final_return_incr;
  }
  // If p > 0.5, work with (1-p) and flipped columns.
  const uint32_t inv = ((!logp) && (p_or_lnp_ddr.x[0] > 0.5)) || (logp && (p_or_lnp_ddr.x[0] > _ddr_log05.x[0]));
  if (inv) {
    swap_i64(&ac, &bd);
    a_minus_d = ab - bd;
    if (!logp) {
      p_or_lnp_ddr = ddr_negate(ddr_subd(p_or_lnp_ddr, 1.0));
    } else {
      p_or_lnp_ddr = ddr_negate(ddr_expm1(p_or_lnp_ddr));
      logp = 0;
    }
  }

  const double m1x = ab;
  const double m2x = abcd - ab;  // = max_d
  const double mx2 = bd;
  const double mxx = abcd;
  // max_d=1 doesn't play well with the current initial-guess algorithm, and is
  // straightforward to handle directly.
  if (max_d == 1) {
    double m11 = m1x - mx2;
    double m12 = mx2;
    double m21 = m2x;
    // double m22 = 0;

    m11 += 1;
    // m22 += 1;
    // min(m12, m21) is guaranteed to be 1, so we don't need to multiply with
    // dd_real precision.
    dd_real lik_ddr = ddr_divd(ddr_maked(m12 * m21), m11);
    // m12 -= 1;
    // m21 -= 1;

    // lik is now equal to pmf(1) / pmf(0).
    // cdf(0) = pmf(0) / (pmf(0) + pmf(1))
    //        = 1 / (1 + lik)
    // p < cdf(0) -> p < 1 / (1 + lik)
    //               p * (1 + lik) < 1
    const dd_real p_ddr = logp? ddr_exp(p_or_lnp_ddr) : p_or_lnp_ddr;
    const int64_t d = ddr_geqd(ddr_mul(p_ddr, ddr_addd(lik_ddr, 1.0)), 1.0);
    return final_return_incr + (inv? (1 - d) : d);
  }
  // We make an initial guess, use Newton's method to refine it (in the same
  // way as Fisher22TwoSidedP() when jumping from one tail to the other),
  // calculate tail probability to sufficient accuracy, and then start moving
  // d inward and updating tail probability until it crosses p.  This avoids
  // repetition of the tail-probability calculation.
  //
  // Initial guess is based on fitting a quadratic to three points of the
  // log-probability function near the mode.  Log-probability is
  //   log(ab! cd! ac! bd! / (a! b! c! d! abcd!))
  //   = <constant> - log(a! b! c! d!)
  // First derivative w.r.t. d is
  //   -digamma(a+1) + digamma(b+1) + digamma(c+1) - digamma(d+1)
  // Second derivative is
  //   -trigamma(a+1) - trigamma(b+1) - trigamma(c+1) - trigamma(d+1)
  //   ~= -(1/(a+1) + 1/(b+1) + 1/(c+1) + 1/(d+1))
  // which is slowly varying near the mode in larger cases.
  //
  // Probable todo: Cornish-Fisher expansion should be applicable here.
  int64_t modal_d = S_CAST(int64_t, 0.5 + (m2x * mx2) / mxx);
  modal_d = modal_d + (modal_d == 0) - (modal_d == abcd);
  const double modal_dd = modal_d;
  // {a, b, c, d} are int64s, {m11, m12, m21, m22} are float64s.
  const double m11_minus_m22 = a_minus_d;
  double m22 = modal_dd;
  double m11 = m11_minus_m22 + m22;
  double m12 = mx2 - m22;
  double m21 = m2x - m22;
  td_real lnprobf_tdr;
  const uint32_t use_tdr = use_tdr_for_hypergeom_lnprob(abcd);
  if (!use_tdr) {
    lnprobf_tdr = tdr_make_dd(ddr_sub(ddr_sort_and_add_4_lfacts(m1x, m2x, mxx - mx2, mx2),
                                      ddr_lfact(mxx)));
  } else {
    td_real numer_tdrs[4];
    numer_tdrs[0] = tdr_lfact(m1x);
    numer_tdrs[1] = tdr_lfact(m2x);
    numer_tdrs[2] = tdr_lfact(mxx - mx2);
    numer_tdrs[3] = tdr_lfact(mx2);
    lnprobf_tdr = tdr_sub(tdr_sort_and_add(4, numer_tdrs),
                          tdr_lfact(mxx));
  }
  const dd_real mode_lnprob_ddr = ddr_sub(ddr_make_td(lnprobf_tdr),
                                          ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
  const dd_real modem1_lnprob_incr_ddr = ddr_log(ddr_accurate_div(ddr_mul2d(m11, m22), ddr_mul2d(m12 + 1, m21 + 1)));
  const dd_real modep1_lnprob_incr_ddr = ddr_log(ddr_accurate_div(ddr_mul2d(m12, m21), ddr_mul2d(m11 + 1, m22 + 1)));

  const double x2_coeff = 0.5 * ddr_add(modem1_lnprob_incr_ddr, modep1_lnprob_incr_ddr).x[0];
  const double x1_coeff = 0.5 * ddr_sub(modep1_lnprob_incr_ddr, modem1_lnprob_incr_ddr).x[0];
  // Better numerical behavior if we let x=0 correspond to the mode.
  const double x0_coeff = mode_lnprob_ddr.x[0];

  // 1. Identify d<mode such that
  //      pmf(d) * (d+1) < p
  //    If no such d>0 exists, initialize d=0.  Also ok to initialize d=0 if
  //    mode is small.
  // 2. Compute pmf(d) to high accuracy.
  // 3. Sum left-tail (<= d) likelihoods.  (Guaranteed to be < p unless d=0, in
  //    which case we can immediately return 0.)  Use float64 instead of
  //    dd_real precision when we can get away with it.
  // 4. Sum inward until the sum >= p, at which point we return d.
  //    Again, use float64 precision as far as we can.
  const dd_real target_lnprob_ddr = logp? p_or_lnp_ddr : ddr_log(p_or_lnp_ddr);
  // Find the x on the left side where the quadratic crosses
  // y=log(p/mode).  If there's no such point, just start at x=mode-1.
  // (Could recalculate target_lnprob with mode replaced with d when
  // log(pmf(d)) is too high.)
  const double search_lnprob = target_lnprob_ddr.x[0] - log(m22);
  // (-b - sqrt(b^2 - 4ac)) / 2a
  const double discrim = x1_coeff * x1_coeff - 4 * x2_coeff * (x0_coeff - search_lnprob);
  if (discrim < 0.0) {
    m22 -= 1;
  } else {
    double sqrt_discrim = sqrt(discrim);
    if (x2_coeff > 0.0) {
      // this shouldn't happen
      sqrt_discrim = -sqrt_discrim;
    }
    m22 = trunc(modal_dd + (sqrt_discrim - x1_coeff) / (2 * x2_coeff));
    if (m22 < 0) {
      m22 = 0;
    }
  }
  // Our relative error budget is usually 2^{-54}.
  // We want to ensure that accumulated error when evaluating the outer part
  // of the tailsum using plain float64 arithmetic < 2^{-55}.  Then the other
  // half of the budget covers dd_real-precision evaluation of log(pmf(d))
  // and the inner part of the tailsum.
  // 2^{-55} corresponds to 0.125 ULPs; Pbinom() comments elaborate on the
  // squared term in the denominator.
  const double tailsum_ddr_end = -log(8 * modal_dd * modal_dd);
  // |log(pmf(0) / pmf(1))| is larger than all the other gaps between
  // adjacent log(pmf()) points to the left of the mode.  So this ensures
  // there's at least one value of d where log(pmf(d)) - target_lnprob is in
  // (lnprob_diff_min, 0], letting us exit the loop; unless log(pmf(0)) >=
  // target_lnprob, in which case we exit the loop at d=0.
  double lnprob_diff_min = log((m11_minus_m22 + 1) / (mx2 * m2x)) * (1 + kSmallEpsilon);
  if (lnprob_diff_min > tailsum_ddr_end) {
    lnprob_diff_min = tailsum_ddr_end;
  }
  dd_real cur_lnprob_ddr;
  while (1) {
    m11 = m11_minus_m22 + m22;
    m12 = mx2 - m22;
    m21 = m2x - m22;
    if (!use_tdr) {
      cur_lnprob_ddr = ddr_sub(ddr_make_td(lnprobf_tdr),
                               ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
    } else {
      td_real tdrs[4];
      tdrs[0] = tdr_lfact(m11);
      tdrs[1] = tdr_lfact(m12);
      tdrs[2] = tdr_lfact(m21);
      tdrs[3] = tdr_lfact(m22);
      cur_lnprob_ddr = ddr_make_td(tdr_sub(lnprobf_tdr,
                                           tdr_sort_and_add(4, tdrs)));
    }
    const double lnprob_diff = cur_lnprob_ddr.x[0] - search_lnprob;
    if (lnprob_diff > 0) {
      if (m22 == 0) {
        break;
      }
      const double ll_deriv = log((m12 + 1) * (m21 + 1) / (m11 * m22));
      m22 -= ceil(lnprob_diff / ll_deriv);
      if (m22 < 0) {
        m22 = 0;
      }
    } else if (lnprob_diff > lnprob_diff_min) {
      break;
    } else {
      const double ll_deriv = log(m12 * m21 / ((m11 + 1) * (m22 + 1)));
      m22 += S_CAST(int64_t, -lnprob_diff / ll_deriv);
    }
  }
  // Express current likelihood as a fraction of p.
  const double tailenter_m22 = m22;
  const dd_real tailenter_lik_ddr = ddr_exp(ddr_sub(cur_lnprob_ddr, target_lnprob_ddr));
  dd_real lik_ddr = tailenter_lik_ddr;
  dd_real tailsum_ddr = tailenter_lik_ddr;
  if (m22 > 0) {
    // Could use geometric-series upper bound on tailsum to raise this
    // threshold.
    const double min_incr_left = (1.0 / (1 << 14)) / (m22 * m22);
    do {
      m12 += 1;
      m21 += 1;
      lik_ddr = ddr_mul(lik_ddr, ddr_accurate_div(ddr_mul2d(m11, m22), ddr_mul2d(m12, m21)));
      m11 -= 1;
      m22 -= 1;
      tailsum_ddr = ddr_add(tailsum_ddr, lik_ddr);
    } while (lik_ddr.x[0] > tailsum_ddr.x[0] * min_incr_left);
    if (m22 > 0) {
      double lik = lik_ddr.x[0];
      double tailsum = 0.0;
      while (1) {
        m12 += 1;
        m21 += 1;
        lik *= m11 * m22 / (m12 * m21);
        m11 -= 1;
        m22 -= 1;
        const double preadd = tailsum;
        tailsum += lik;
        if (tailsum == preadd) {
          break;
        }
      }
      tailsum_ddr = ddr_addd(tailsum_ddr, tailsum);
    }
    lik_ddr = tailenter_lik_ddr;
    m22 = tailenter_m22;
    m11 = m11_minus_m22 + m22;
    m12 = mx2 - m22;
    m21 = m2x - m22;
  }
  while (ddr_ltd(tailsum_ddr, 1.0)) {
    m11 += 1;
    m22 += 1;
    lik_ddr = ddr_mul(lik_ddr, ddr_accurate_div(ddr_mul2d(m12, m21), ddr_mul2d(m11, m22)));
    m12 -= 1;
    m21 -= 1;
    tailsum_ddr = ddr_add(tailsum_ddr, lik_ddr);
  }
  return final_return_incr + S_CAST(int64_t, inv? (m2x - m22) : m22);
}

#ifdef __cplusplus
}
#endif
