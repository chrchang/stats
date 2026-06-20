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

#include "fisher.h"

#include <math.h>

#include "hypergeom_detail.h"
#include "nchypergeom_fisher.h"
#include "plink2_float.h"
#include "plink2_highprec.h"
#include "special_func.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// obs_m11 + obs_m12 + obs_m21 + obs_m22 assumed to be <2^52.
// They're defined as int64_ts instead of uint64_ts since signed int <->
// floating-point conversions are sometimes faster than the same-width unsigned
// int <-> floating-point conversions.
//
// possible todo: add support for R fisher.test's 'or' (null-hypothesis odds
// ratio) parameter.  This is straightforward for logp=False, and there should
// be no practical need for logp=True (which isn't supported by R) there.
double Fisher22TwoSidedP(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, int32_t midp, uint32_t logp) {
  // Normalize: m11 >= m22, m12 >= m21, m11*m22 < m12*m21.
  // Note that the first two are reversed from PLINK 1.9, to get rid of
  // spurious index differences between Fisher22 and Fisher23.
  if (obs_m11 < obs_m22) {
    swap_i64(&obs_m11, &obs_m22);
  }
  if (obs_m12 < obs_m21) {
    swap_i64(&obs_m12, &obs_m21);
  }
  if (ddr_gt(ddr_mul2d(obs_m11, obs_m22), ddr_mul2d(obs_m12, obs_m21))) {
    swap_i64(&obs_m11, &obs_m12);
    swap_i64(&obs_m21, &obs_m22);
  }
  // I experimented with making more int -> double casts explicit, but
  // concluded that it was hurting readability more than it was helping.
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;
  if (!midp) {
    // Fast path for p=1.
    if (ddr_geq(ddr_mul2d(m11 + 1, m22 + 1), ddr_mul2d(m12, m21))) {
      return logp? 0.0 : 1.0;
    }
  }

  // Iterate outward to floating-point precision limit.
  double lik = 1;
  double tail_sum = 1 - 0.5 * midp;
  while (1) {
    m12 += 1;
    m21 += 1;
    lik *= (m11 * m22) / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
  }
  // In the common case, where we're close enough to the mode that float64
  // underflow/overflow isn't an issue, use the original algorithm: sum all
  // center relative-likelihoods, sum far-tail relative-likelihoods to
  // floating-point precision limit, return
  //   log(tail_sum / (tail_sum + center_sum))
  //
  // As with HweLnP(), we note that if we're within 172 steps of the mode and
  // the starting relative-likelihood is normalized to 1, the modal
  // relative-likelihood can be loosely bounded above by
  //   ((172^172) / 172!)^4 ~= 5.3e+292
  // which leaves enough headroom to accumulate the rest of the center-sum and
  // represent intermediate values without overflowing.
  if (ddr_geq(ddr_mul2d(obs_m11 + 172, obs_m22 + 172), ddr_mul2d(obs_m12 - 172, obs_m21 - 172))) {
    lik = 1;
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21;
    m22 = obs_m22;
    double center_sum = 0.5 * midp;
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= (m12 * m21) / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      // Number of center contingency tables is maximized with obs_m22 = 0,
      // modal_m22 = 172, other values large.
      // Since 1 + 1/2 + ... + 1/172 < 1/173 + ... + 1/53000, we're limited to
      // ~53000 tables.  Each lik update involves 4 operations which can each
      // introduce up to 0.5 ULP relative error under the default rounding
      // mode.
      if (lik < 1 + 53000 * 2 * k2m52) {
        if (lik <= 1 - 53000 * 2 * k2m52) {
          tail_sum += lik;
          break;
        }
        // Near-tie.  True value of lik can be greater than, equal to, or
        // less than 1.
        const int64_t m22_incr = S_CAST(int64_t, m22) - obs_m22;
        td_real starting_lnprobv_tdr = tdr_make1(DBL_MAX);
        const intptr_t cmp_result = HypergeomCompare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, &starting_lnprobv_tdr, &lik);
        if (cmp_result <= 0) {
          tail_sum += lik;
          if (midp && (cmp_result == 0)) {
            tail_sum -= 0.5;
            center_sum += 0.5;
          }
          break;
        }
      }
      center_sum += lik;
    }
    // Continue down tail to floating-point precision limit.
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= (m12 * m21) / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        break;
      }
    }
    const double pval = tail_sum / (tail_sum + center_sum);
    return logp? log(pval) : pval;
  }

  // Now we want to jump near the other tail, without evaluating that many
  // contingency table log-likelihoods along the way.
  //
  // Each full log-likelihood evaluation requires 4 ddr_lfact() or tdr_lfact()
  // calls.  Since they are now performed with extra precision, they require
  // hundreds or thousands of floating-point operations, so we want to limit
  // ourselves to 1-2 full evaluations most of the time.  (Possible todo: use
  // lower-accuracy Lfact() to jump around, followed by {d,t}dr_lfact() when
  // exiting the loop.  Should be an easy performance win, but there's a
  // complexity cost so I'll wait until I see a scenario where this branch
  // executes frequently...)
  //
  // The current heuristic starts by reflecting (obs_m21 + m21) * 0.5 across
  // the mode, performing a full log-likelihood check at an adjacent valid
  // point.  (It is convenient to focus on m21 here, since m21=0 corresponds to
  // the outermost table on this tail.)  Hopefully we find that we're in
  // (starting_lnprob - lnprobv_diff_min, starting_lnprob], so we're at or
  // near a table that actually contributes to the tail-sum.  (This window is
  // chosen to be wide enough to guarantee that at least one point falls
  // inside.)
  //
  // If not, we jump again, using Newton's method.
  // If m21 is too high (i.e. current log-likelihood is too high), decreasing
  // m12 by 1 would multiply the likelihood by
  //   (m11 + 1) * (m22 + 1) / (m12 * m21)
  // If m12 is too low, increasing m12 by 1 would multiply the likelihood by
  //   (m12 + 1) * (m21 + 1) / (m11 * m22)
  // We use the negative-log of the first expression as the Newton's method
  // f'(x) when we're jumping to lower m12, and the log of the second
  // expression when we're jumping to higher m12.
  // f''(x) is always negative, so we can aim for starting_lnprob instead of
  // the middle of the interval.

  const double m12_minus_m21 = obs_m12 - obs_m21;
  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx1 = obs_m11 + obs_m21;
  const double mxx = m1x + m2x;
  {
    // x=modal_m21 satisfies
    //    (m2x - x) * (mx1 - x) = x * (x + m12_minus_m21)
    // -> (x - m2x) * (x - mx1) = x * (x + m12_minus_m21)
    // -> x^2 + x*(-m2x - mx1) + m2x*mx1 = x^2 + x*m12_minus_m21
    // -> m2x*mx1 = x*(m12_minus_m21 + m2x + mx1)
    // -> x = m2x*mx1 / mxx
    const double modal_m21 = m2x * mx1 / mxx;
    m21 = 2 * modal_m21 - (m21 + obs_m21) * 0.5;
    // Round down (to guarantee we've actually moved to the other side of the
    // mode) and clamp.
    m21 = S_CAST(int64_t, m21);
    if (m21 < 0) {
      m21 = 0;
    }
  }
  // Extremal case: m11=mx1, m12=m12_minus_m21, m21=0, m22=m2x
  //   log((m12 + 1) * (m21 + 1) / (m11 * m22))
  // = log((m12_minus_m21 + 1) / (mx1 * m2x))
  double lnprobv_diff_min = log((m12_minus_m21 + 1) / (mx1 * m2x)) * (1 + kSmallEpsilon);
  const uint32_t tdr_lnprobv_needed = use_tdr_for_hypergeom_lnprob(obs_m11 + obs_m12 + obs_m21 + obs_m22);
  td_real starting_lnprobv_tdr;
  dd_real starting_lnprob_ddr;
  if (!tdr_lnprobv_needed) {
    const dd_real starting_lnprobv_ddr = ddr_negate(ddr_sort_and_add_4_lfacts(obs_m11, obs_m12, obs_m21, obs_m22));
    starting_lnprobv_tdr = tdr_make_dd(starting_lnprobv_ddr);

    const dd_real lnprobf_ddr =
      ddr_sub(ddr_sort_and_add_4_lfacts(m1x, m2x, mx1, obs_m12 + obs_m22),
              ddr_lfact(mxx));
    starting_lnprob_ddr = ddr_add(lnprobf_ddr, starting_lnprobv_ddr);
  } else {
    td_real tdrs[4];
    tdrs[0] = tdr_lfact(obs_m11);
    tdrs[1] = tdr_lfact(obs_m12);
    tdrs[2] = tdr_lfact(obs_m21);
    tdrs[3] = tdr_lfact(obs_m22);
    starting_lnprobv_tdr = tdr_negate(tdr_sort_and_add(4, tdrs));

    tdrs[0] = tdr_lfact(m1x);
    tdrs[1] = tdr_lfact(m2x);
    tdrs[2] = tdr_lfact(mx1);
    tdrs[3] = tdr_lfact(obs_m12 + obs_m22);
    const td_real lnprobf_tdr = tdr_sub(tdr_sort_and_add(4, tdrs),
                                        tdr_lfact(mxx));
    starting_lnprob_ddr = ddr_make_td(tdr_add(lnprobf_tdr, starting_lnprobv_tdr));
  }
  while (1) {
    m11 = mx1 - m21;
    m12 = m12_minus_m21 + m21;
    m22 = m2x - m21;
    double lnprobv_diff;
    if (!tdr_lnprobv_needed) {
      const dd_real lnprobv_ddr =
        ddr_negate(ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
      lnprobv_diff = ddr_sub(lnprobv_ddr, ddr_make_td(starting_lnprobv_tdr)).x[0];
    } else {
      td_real tdrs[4];
      tdrs[0] = tdr_lfact(m11);
      tdrs[1] = tdr_lfact(m12);
      tdrs[2] = tdr_lfact(m21);
      tdrs[3] = tdr_lfact(m22);
      const td_real lnprobv_tdr = tdr_negate(tdr_sort_and_add(4, tdrs));
      lnprobv_diff = tdr_sub(lnprobv_tdr, starting_lnprobv_tdr).x[0];
    }
    // Could tighten this threshold further.  But code is correct as long as
    // we're guaranteed to enter the "lik < 2 - one_minus_scaled_eps" branch
    // for positive lnprobv_diff.
    if (lnprobv_diff >= k2m53) {
      if (m21 == 0) {
        // All tables on this tail have higher likelihood than the starting
        // table.  Exit.
        return join_log_and_nonlog(starting_lnprob_ddr, tail_sum, logp);
      }
      const double ll_deriv = log(m12 * m21 / ((m11 + 1) * (m22 + 1)));
      // Round up, to guarantee that we make progress.
      // (lnprobv_diff is positive and ll_deriv is negative.)
      // This may overshoot.  But the function is guaranteed to terminate
      // because we never overshoot (and we do always make progress on each
      // step) once we're on the other side.
      m21 -= ceil(-lnprobv_diff / ll_deriv);
      if (m21 < 0) {
        m21 = 0;
      }
    } else {
      double ll_deriv = DBL_MAX;
      if ((lnprobv_diff_min < -53 * kLn2) && (m21 > 0)) {
        // Tighten this threshold if that lets us sum fewer terms later.
        ll_deriv = log((m12 + 1) * (m21 + 1) / (m11 * m22));
        lnprobv_diff_min = ll_deriv * (1 + kSmallEpsilon);
      }
      if (lnprobv_diff > lnprobv_diff_min) {
        lik = exp(lnprobv_diff);
        break;
      } else {
        if (ll_deriv == DBL_MAX) {
          ll_deriv = log((m12 + 1) * (m21 + 1) / (m11 * m22));
        }
        // Round down, to guarantee we don't overshoot.
        // We're guaranteed to make progress because of how lnprobv_diff_min
        // was set.
        // m11 * m22 < exp(-lnprobv_diff_min), and (m12 + 1) * (m21 + 1) >= 1.
        m21 += S_CAST(int64_t, lnprobv_diff / ll_deriv);
      }
    }
  }
  // Sum toward center, until lik >= 1.
  // lik should be accurate to 3 ULP as we enter this loop (max 1.5 ULP
  // observed error from exp, tiny bit over 0.5 from lnprobv_diff), so near-tie
  // detection can use a tight epsilon here.
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  // Save where we're starting on this tail, which isn't necessarily on the
  // boundary.  We sum inward until relative-likelihood > 1, then we jump back
  // to tailenter_{m11,m12,m21,m22} and sum outward.
  const double tailenter_lik = lik;
  const double tailenter_m11 = m11;
  const double tailenter_m12 = m12;
  const double tailenter_m21 = m21;
  const double tailenter_m22 = m22;
  while (1) {
    while (lik <= one_minus_scaled_eps) {
      tail_sum += lik;
      m12 += 1;
      m21 += 1;
      lik *= m11 * m22 / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      one_minus_scaled_eps -= 2 * k2m52;
    }
    if (lik >= 2 - one_minus_scaled_eps) {
      break;
    }
    const int64_t m22_incr = S_CAST(int64_t, m22) - obs_m22;
    if (!tdr_lnprobv_needed) {
      starting_lnprobv_tdr.x[2] = DBL_MAX;
    }
    const intptr_t cmp_result = HypergeomCompare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, &starting_lnprobv_tdr, &lik);
    if (cmp_result >= 0) {
      if (cmp_result == 0) {
        tail_sum += 1 - 0.5 * midp;
      }
      break;
    }
    one_minus_scaled_eps = 1 - 3 * k2m52;
    // In very large cases, additional step(s) may be required?
  }
  // Sum away from center, until sums stop changing.
  lik = tailenter_lik;
  m11 = tailenter_m11;
  m12 = tailenter_m12;
  m21 = tailenter_m21;
  m22 = tailenter_m22;
  while (1) {
    m11 += 1;
    m22 += 1;
    lik *= m12 * m21 / (m11 * m22);
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
    m12 -= 1;
    m21 -= 1;
  }
  return join_log_and_nonlog(starting_lnprob_ddr, tail_sum, logp);
}

double Fisher22OddsRatio(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22) {
  // R's conditional MLE.
  const int32_t m11_is_maximal = (obs_m12 == 0) || (obs_m21 == 0);
  if ((obs_m11 == 0) || (obs_m22 == 0)) {
    // If table is degenerate (a row and/or column sums to 0), I think 1 is a
    // more appropriate return value than 0.
    return S_CAST(double, m11_is_maximal);
  }
  if (m11_is_maximal) {
    return INFINITY_D;
  }
  const double m11d = obs_m11;
  const int64_t m1 = obs_m11 + obs_m12;
  const int64_t m2 = obs_m21 + obs_m22;
  const int64_t n = obs_m11 + obs_m21;
  // Set initial guess to unconditional MLE, use Newton's method until we have
  // float32-level accuracy (relative error < 2^{-24}).
  double odds = S_CAST(double, obs_m11) * S_CAST(double, obs_m22) / (S_CAST(double, obs_m12) * S_CAST(double, obs_m21));
  while (1) {
    const dd_real mean_ddr = MeanFNCHypergeo(m1, m2, n, odds);
    const double mean_delta = -ddr_subd(mean_ddr, m11d).x[0];
    const double mean_deriv_odds = MeanFNCHypergeoDerivOdds(m1, m2, n, odds, mean_ddr);
    assert(mean_deriv_odds > 0.0);
    const double odds_incr = mean_delta / mean_deriv_odds;
    // printf("%.17g: %.17g %.17g %.17g\n", odds, mean_delta, mean_deriv_odds, odds_incr);
    odds += odds_incr;
    if (fabs(odds_incr) < odds * k2m24) {
      return odds;
    }
  }
}

static inline double logit(double p) {
  return log(p / (1-p));
}

double Fisher22OddsRatioQuantileMatch(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double target_p) {
  if ((obs_m12 == 0) || (obs_m21 == 0)) {
    // Degenerate cases.
    if ((obs_m11 == 0) || (obs_m22 == 0)) {
      return 1.0;
    }
    return INFINITY_D;
  }
  if (target_p < k2m24) {
    return INFINITY_D;
  }
  if (target_p > 1 - k2m24) {
    return 0.0;
  }
  // There should be no risk of odds overflow/underflow with these target_p
  // bounds.

  // Initial guess: calculate unconditional MLE, map target_p to z-score, etc.
  const double m11_p05 = S_CAST(double, obs_m11) + 0.5;
  const double m12_p05 = S_CAST(double, obs_m12) + 0.5;
  const double m21_p05 = S_CAST(double, obs_m21) + 0.5;
  const double m22_p05 = S_CAST(double, obs_m22) + 0.5;
  const double ln_or = log(m11_p05 * m22_p05 / (m12_p05 * m21_p05));
  const double se_ln_or = sqrt(1.0/m11_p05 + 1.0/m12_p05 + 1.0/m21_p05 + 1.0/m22_p05);
  const double zscore = QuantileToZscoreD(target_p, 0);
  double odds1 = exp(ln_or - zscore * se_ln_or);
  const double logit_target_p = logit(target_p);
  while (1) {
    const double odds2 = odds1 * (1.0 + k2m25);
    odds1 = odds1 * (1.0 - k2m25);
    // As odds ratio increases, p decreases, so we'll get p2 <= p1.
    double p1;
    double p2;
    P_FNCHypergeoTwoOdds(obs_m11, obs_m12, obs_m21, obs_m22, odds1, odds2, &p1, &p2);
    // Interpolate or take a Newton step, treating logit(p) as a linear
    // function of log-odds.
    const double logit_p1 = logit(p1);
    const double logit_p2 = logit(p2);
    if ((target_p >= p2) && (target_p <= p1)) {
      const double interp = (logit_target_p - logit_p2) / (logit_p1 - logit_p2);
      return odds2 * (1 - k2m24 * interp);
    }
    const double logitp_deriv_lnodds = (logit_p2 - logit_p1) / k2m24;  // negative
    if (target_p < p2) {
      const double lnodds_incr = (logit_target_p - logit_p2) / logitp_deriv_lnodds;
      odds1 = odds2 * exp(lnodds_incr);
    } else {
      // target_p > p1
      const double lnodds_incr = (logit_target_p - logit_p1) / logitp_deriv_lnodds;
      odds1 *= exp(lnodds_incr);
    }
  }
}

// Switch between log- and regular representations at kSwitchThresh.
// 2^890 leaves enough headroom for at least 4 more multiplies by obs_total.
// (Useful to reduce this to e.g. 2^150 when testing correctness, though it
// isn't currently safe to go all the way down to 1 due to algorithmic
// assumptions.)
static const double kSwitchThresh = k2p800 * k2p50 * (1LL << 40);
static const double kLnSwitchThresh = 890.0 * kLn2;
static const double kJumpThresh = 314.0; // chosen to guarantee base_lik < kSwitchThresh in non-jumping case
// static const double kSwitchThresh = k2p100 * k2p50;
// static const double kLnSwitchThresh = 150.0 * kLn2;

// Since each Fisher's exact test contigency table has rows and columns, we use
// 'rank' instead of 'row' to refer to the set of tables with 3rd column held
// constant.
void Fisher23LnStartingRank(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, double* tail_sum_ptr, dd_real* starting_lnprobv_ddr_ptr, int32_t* tie_ct_ptr, double* orig_base_likl_ptr, double* orig_base_lnlikl_ptr, double* orig_base_epsl_ptr, double* orig_base_likr_ptr, double* orig_base_lnlikr_ptr, double* orig_base_epsr_ptr, double* orig_saved_l11_ptr, double* orig_saved_l12_ptr, double* orig_saved_l21_ptr, double* orig_saved_l22_ptr, double* orig_saved_r11_ptr, double* orig_saved_r12_ptr, double* orig_saved_r21_ptr, double* orig_saved_r22_ptr) {
  // possible todo: have this and Fisher22TwoSidedP() call a shared function
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;
  *starting_lnprobv_ddr_ptr =
    ddr_negate(ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
  double lik = 1;
  int32_t tie_ct = 1;
  double tail_sum = 1;
  // No guaranteed normalization beyond m11+m21 >= m12+m22.  In particular,
  // m11 < m22 is possible.

  const double m1x = m11 + m12;
  const double m2x = m21 + m22;
  const double mx2 = m12 + m22;
  const double mxx = m1x + m2x;
  const double modal_m22 = m2x * mx2 / mxx;
  const double delta = m22 - modal_m22;
  // As with PLINK 1.9 fisher23(), "left tail" corresponds to small m22, and
  // "right tail" corresponds to large m22.
  if (delta > 0) {
    // Starting m22 is beginning of right tail.
    *orig_base_likr_ptr = 1;
    *orig_base_epsr_ptr = k2m52;
    *orig_saved_r11_ptr = m11;
    *orig_saved_r12_ptr = m12;
    *orig_saved_r21_ptr = m21;
    *orig_saved_r22_ptr = m22;
    // Iterate outward to floating-point precision limit.
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= m12 * m21 / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        break;
      }
    }
    if (delta < kJumpThresh) {
      // Jump back to starting table, and iterate inward until we find the
      // start of the other tail.
      lik = 1;
      m11 = obs_m11;
      m12 = obs_m12;
      m21 = obs_m21;
      m22 = obs_m22;
      double one_plus_scaled_eps = 1 + k2m52;
      while (1) {
        const double m11_x_m22 = m11 * m22;
        // Don't want to wait until lik becomes 0, since we want m11 >= 0 and
        // m22 >= 0 guaranteed when we save m11 and m22 on loop exit.
        if (m11_x_m22 == 0) {
          break;
        }
        m12 += 1;
        m21 += 1;
        lik *= m11_x_m22 / (m12 * m21);
        m11 -= 1;
        m22 -= 1;
        one_plus_scaled_eps += 2 * k2m52;
        if (lik < one_plus_scaled_eps) {
          if (lik <= 2 - one_plus_scaled_eps) {
            tail_sum += lik;
            break;
          }
          // Near-tie.  True value of lik can be greater than, equal to, or
          // less than 1.
          const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
          const intptr_t cmp_result = HypergeomCompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, *starting_lnprobv_ddr_ptr, &lik);
          one_plus_scaled_eps = 1 + 3 * k2m52;
          if (cmp_result <= 0) {
            tail_sum += lik;
            tie_ct += (cmp_result == 0);
            break;
          }
        }
      }
      *orig_saved_l11_ptr = m11;
      *orig_saved_l12_ptr = m12;
      *orig_saved_l21_ptr = m21;
      *orig_saved_l22_ptr = m22;
      *orig_base_likl_ptr = lik;
      *orig_base_epsl_ptr = one_plus_scaled_eps - 1;
    } else {
      // Jump to other tail.  Round down and clamp.
      const double m22_minus_m11 = m22 - m11; // usually but not always nonpositive
      m22 = S_CAST(int32_t, modal_m22 - delta);
      const double min_m22 = MAXV(0, m22_minus_m11);
      if (m22 < min_m22) {
        m22 = min_m22;
      }
      while (1) {
        m11 = m22 - m22_minus_m11;
        m12 = mx2 - m22;
        m21 = m2x - m22;
        const dd_real lnprobv_ddr =
          ddr_negate(ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
        const dd_real lnprobv_diff_ddr = ddr_sub(lnprobv_ddr, *starting_lnprobv_ddr_ptr);
        const double lnprobv_diff = lnprobv_diff_ddr.x[0];
        if (lnprobv_diff >= k2m53) {
          if (m22 == min_m22) {
            // All tables on this tail have higher likelihood than the starting
            // table.  Exit.
            *tie_ct_ptr = 1;
            *orig_saved_l11_ptr = m11;
            *orig_saved_l12_ptr = m12;
            *orig_saved_l21_ptr = m21;
            *orig_saved_l22_ptr = m22;
            *tail_sum_ptr = tail_sum;
            if (lnprobv_diff < kLnSwitchThresh) {
              *orig_base_likl_ptr = ddr_exp(lnprobv_diff_ddr).x[0];
              *orig_base_epsl_ptr = 2 * k2m52;
            } else {
              *orig_base_likl_ptr = 0;
              *orig_base_lnlikl_ptr = lnprobv_diff;
              *orig_base_epsl_ptr = (1 + ceil(lnprobv_diff)) * k2m52;
            }
            return;
          }
          const double ll_deriv = log((m12 + 1) * (m21 + 1) / (m11 * m22));
          // Round up, to guarantee that we make progress.
          // This may overshoot.  But the function is guaraneed to terminate
          // because we never overshoot (and we do always make progress on each
          // step) once we're on the other side.
          m22 -= ceil(lnprobv_diff / ll_deriv);
          if (m22 < min_m22) {
            m22 = min_m22;
          }
        } else if (lnprobv_diff > -62 * kLn2) {
          lik = exp(lnprobv_diff);
          break;
        } else {
          const double ll_deriv = log(m12 * m21 / ((m11 + 1) * (m22 + 1)));
          // Round down, to guarantee we don't overshoot.
          // (lnprobv_diff is negative and ll_deriv is positive.)
          // We're guaranteed to make progress, since lnprobv_diff <=
          // -62 * log(2) and m12 * m21 >= 1.
          m22 -= S_CAST(int64_t, lnprobv_diff / ll_deriv);
        }
      }
      // Sum toward center, until lik >= 1.
      double one_minus_scaled_eps = 1 - 3 * k2m52;
      const double tailenter_lik = lik;
      const double m11_tail = m11;
      const double m12_tail = m12;
      const double m21_tail = m21;
      const double m22_tail = m22;
      while (lik <= one_minus_scaled_eps) {
        tail_sum += lik;
        m11 += 1;
        m22 += 1;
        lik *= m12 * m21 / (m11 * m22);
        m12 -= 1;
        m21 -= 1;
        one_minus_scaled_eps -= 2 * k2m52;
      }
      if (lik < 2 - one_minus_scaled_eps) {
        const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
        const intptr_t cmp_result = HypergeomCompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, *starting_lnprobv_ddr_ptr, &lik);
        one_minus_scaled_eps = 1 - 3 * k2m52;
        if (cmp_result <= 0) {
          tail_sum += lik;
          tie_ct += (cmp_result == 0);
        }
      }
      *orig_saved_l11_ptr = m11;
      *orig_saved_l12_ptr = m12;
      *orig_saved_l21_ptr = m21;
      *orig_saved_l22_ptr = m22;
      *orig_base_likl_ptr = lik;
      *orig_base_epsl_ptr = 1 - one_minus_scaled_eps;
      lik = tailenter_lik;
      m11 = m11_tail;
      m12 = m12_tail;
      m21 = m21_tail;
      m22 = m22_tail;
    }
    *tie_ct_ptr = tie_ct;
    // Sum away from center, until sums stop changing.
    while (1) {
      m12 += 1;
      m21 += 1;
      lik *= m11 * m22 / (m12 * m21);
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        break;
      }
      m11 -= 1;
      m22 -= 1;
    }
    *tail_sum_ptr = tail_sum;
    return;
  }
  *orig_base_likl_ptr = 1;
  *orig_base_epsl_ptr = k2m52;
  *orig_saved_l11_ptr = obs_m11;
  *orig_saved_l12_ptr = obs_m12;
  *orig_saved_l21_ptr = obs_m21;
  *orig_saved_l22_ptr = obs_m22;
  while (1) {
    m12 += 1;
    m21 += 1;
    lik *= m11 * m22 / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
  }
  if (delta > -kJumpThresh) {
    lik = 1;
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21;
    m22 = obs_m22;
    double one_plus_scaled_eps = 1 + k2m52;
    while (1) {
      const double m12_x_m21 = m12 * m21;
      if (m12_x_m21 == 0) {
        break;
      }
      m11 += 1;
      m22 += 1;
      lik *= m12_x_m21 / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      one_plus_scaled_eps += 2 * k2m52;
      if (lik < one_plus_scaled_eps) {
        if (lik <= 2 - one_plus_scaled_eps) {
          tail_sum += lik;
          break;
        }
        const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
        const intptr_t cmp_result = HypergeomCompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, *starting_lnprobv_ddr_ptr, &lik);
        one_plus_scaled_eps = 1 + 3 * k2m52;
        if (cmp_result <= 0) {
          tail_sum += lik;
          tie_ct += (cmp_result == 0);
          break;
        }
      }
    }
    *orig_saved_r11_ptr = m11;
    *orig_saved_r12_ptr = m12;
    *orig_saved_r21_ptr = m21;
    *orig_saved_r22_ptr = m22;
    *orig_base_likr_ptr = lik;
    *orig_base_epsr_ptr = one_plus_scaled_eps - 1;
  } else {
    // Jump to other tail.
    const double m21_minus_m12 = m21 - m12; // usually but not always nonpositive
    const double mx1 = obs_m11 + obs_m21;
    const double modal_m21 = m2x - modal_m22;
    m21 = S_CAST(int32_t, modal_m21 + delta);
    const double min_m21 = MAXV(0, m21_minus_m12);
    if (m21 < min_m21) {
      m21 = min_m21;
    }
    while (1) {
      m11 = mx1 - m21;
      m12 = m21 - m21_minus_m12;
      m22 = m2x - m21;
      const dd_real lnprobv_ddr =
        ddr_negate(ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
      const dd_real lnprobv_diff_ddr = ddr_sub(lnprobv_ddr, *starting_lnprobv_ddr_ptr);
      const double lnprobv_diff = lnprobv_diff_ddr.x[0];
      if (lnprobv_diff >= k2m53) {
        if (m21 == min_m21) {
          // All tables on this tail have higher likelihood than the starting
          // table.  Exit.
          *tie_ct_ptr = 1;
          *orig_saved_r11_ptr = m11;
          *orig_saved_r12_ptr = m12;
          *orig_saved_r21_ptr = m21;
          *orig_saved_r22_ptr = m22;
          *tail_sum_ptr = tail_sum;
          if (lnprobv_diff < kLnSwitchThresh) {
            *orig_base_likr_ptr = ddr_exp(lnprobv_diff_ddr).x[0];
            *orig_base_epsr_ptr = 2 * k2m52;
          } else {
            *orig_base_likr_ptr = 0;
            *orig_base_lnlikr_ptr = lnprobv_diff;
            *orig_base_epsr_ptr = (1 + ceil(lnprobv_diff)) * k2m52;
          }
          return;
        }
        // Derivative is w.r.t. m21.
        const double ll_deriv = log((m11 + 1) * (m22 + 1) / (m12 * m21));
        m21 -= ceil(lnprobv_diff / ll_deriv);
        if (m21 < min_m21) {
          m21 = min_m21;
        }
      } else if (lnprobv_diff > -62 * kLn2) {
        lik = exp(lnprobv_diff);
        break;
      } else {
        const double ll_deriv = log(m11 * m22 / ((m12 + 1) * (m21 + 1)));
        m21 += S_CAST(int64_t, -lnprobv_diff / ll_deriv);
      }
    }
    // Sum toward center, until lik >= 1.
    double one_minus_scaled_eps = 1 - 3 * k2m52;
    const double tailenter_lik = lik;
    const double m11_tail = m11;
    const double m12_tail = m12;
    const double m21_tail = m21;
    const double m22_tail = m22;
    while (lik <= one_minus_scaled_eps) {
      tail_sum += lik;
      m12 += 1;
      m21 += 1;
      lik *= m11 * m22 / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      one_minus_scaled_eps -= 2 * k2m52;
    }
    if (lik < 2 - one_minus_scaled_eps) {
      const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
      const intptr_t cmp_result = HypergeomCompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, *starting_lnprobv_ddr_ptr, &lik);
      one_minus_scaled_eps = 1 - 3 * k2m52;
      if (cmp_result <= 0) {
        tail_sum += lik;
        tie_ct += (cmp_result == 0);
      }
    }
    *orig_saved_r11_ptr = m11;
    *orig_saved_r12_ptr = m12;
    *orig_saved_r21_ptr = m21;
    *orig_saved_r22_ptr = m22;
    *orig_base_likr_ptr = lik;
    *orig_base_epsr_ptr = 1 - one_minus_scaled_eps;
    lik = tailenter_lik;
    m11 = m11_tail;
    m12 = m12_tail;
    m21 = m21_tail;
    m22 = m22_tail;
  }
  *tie_ct_ptr = tie_ct;
  // Sum away from center, until sums stop changing.
  while (1) {
    m11 += 1;
    m22 += 1;
    lik *= m12 * m21 / (m11 * m22);
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
    m12 -= 1;
    m21 -= 1;
  }
  *tail_sum_ptr = tail_sum;
  return;
}

intptr_t Fisher23Compare(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m13, uint32_t obs_m21, uint32_t obs_m22, uint32_t obs_m23, uint32_t cur_m11, uint32_t cur_m12, dd_real neg_numer_ddr, double* dbl_ptr) {
  // Fisher 2x3 contingency table probability is
  //
  //   (m11+m12+m13)! (m21+m22+m23)! (m11+m21)! (m12+m22)! (m13+m23)!
  //   --------------------------------------------------------------
  //      m11! m12! m13! m21! m22! m23! (m11+m12+m13+m21+m22+m23)!
  //
  // so likelihood ratio of interest is
  //
  //   obs_m11! obs_m12! obs_m13! obs_m21! obs_m22! obs_m23!
  //   -----------------------------------------------------
  //   cur_m11! cur_m12! cur_m13! cur_m21! cur_m22! cur_m23!
  //
  // where we infer the remaining cur_ values from the fixed row/column sums.
  const uint32_t cur_m13 = obs_m11 + obs_m12 + obs_m13 - cur_m11 - cur_m12;
  const uint32_t cur_m21 = obs_m11 + obs_m21 - cur_m11;
  const uint32_t cur_m22 = obs_m12 + obs_m22 - cur_m12;
  const uint32_t cur_m23 = obs_m13 + obs_m23 - cur_m13;
  uint64_t numer_factorial_args[6];
  numer_factorial_args[0] = obs_m11;
  numer_factorial_args[1] = obs_m12;
  numer_factorial_args[2] = obs_m13;
  numer_factorial_args[3] = obs_m21;
  numer_factorial_args[4] = obs_m22;
  numer_factorial_args[5] = obs_m23;
  uint64_t denom_factorial_args[6];
  denom_factorial_args[0] = cur_m11;
  denom_factorial_args[1] = cur_m12;
  denom_factorial_args[2] = cur_m13;
  denom_factorial_args[3] = cur_m21;
  denom_factorial_args[4] = cur_m22;
  denom_factorial_args[5] = cur_m23;
  td_real ln_odds_ratio_tdr = tdr_make1(0.0);
  td_real neg_numer_tdr = {{neg_numer_ddr.x[0], neg_numer_ddr.x[1], DBL_MAX}};
  return CompareFactorialProducts(6, tdr_make1(1.0), 0, 0, numer_factorial_args, denom_factorial_args, &neg_numer_tdr, &ln_odds_ratio_tdr, dbl_ptr);
}

// 'Left' = small m11 and m22.
void Fisher23LnPLeftTailsum(dd_real starting_lnprobv_ddr, uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m13, uint32_t obs_m21, uint32_t obs_m22, uint32_t obs_m23, double* base_likp, double* base_lnlikp, double* base_epsp, double* saved_m11p, double* saved_m12p, double* saved_m21p, double* saved_m22p, int32_t* tie_ctp, double* totalp, uint32_t* center_is_emptyp) {
  double total = 0;
  double lik = *base_likp;
  double cur_eps = *base_epsp;
  double m11 = *saved_m11p;
  double m12 = *saved_m12p;
  double m21 = *saved_m21p;
  double m22 = *saved_m22p;
  // identify beginning (center-facing side) of tail
  if (lik == 0.0) {
    double last_lnp = *base_lnlikp;
    if (last_lnp >= kLnSwitchThresh) {
      while (1) {
        const double prev_numer = m11 * m22;
        if (prev_numer == 0) {
          // lowest-likelihood table on this side is still too probable
          *base_lnlikp = last_lnp;
          *base_epsp = cur_eps;
          *center_is_emptyp = 0;
          *saved_m11p = m11;
          *saved_m12p = m12;
          *saved_m21p = m21;
          *saved_m22p = m22;
          return;
        }
        m12 += 1;
        m21 += 1;
        const double lnlik_incr = log(prev_numer / (m12 * m21));
        last_lnp += lnlik_incr;
        m11 -= 1;
        m21 -= 1;
        cur_eps += 3 * k2m52;
        if (lnlik_incr <= -2) {
          cur_eps += (trunc(-lnlik_incr) - 1) * k2m52;
        }
        // no risk of last_lnp dropping below cur_eps
      }
    }
    lik = exp(last_lnp);
  }
  double m11_tail;
  double m12_tail;
  double m21_tail;
  double m22_tail;
  if (lik >= 1 + cur_eps) {
    while (1) {
      const double prev_numer = m11 * m22;
      if (prev_numer == 0) {
        // lowest-likelihood table on this side is still too probable
        *center_is_emptyp = 0;
        *saved_m11p = m11;
        *saved_m12p = m12;
        *saved_m21p = m21;
        *saved_m22p = m22;
        if (lik < kSwitchThresh) {
          *base_likp = lik;
          *base_epsp = cur_eps;
        } else {
          *base_likp = 0;
          const double last_lnp = log(lik);
          *base_lnlikp = last_lnp;
          *base_epsp = cur_eps + ceil(last_lnp) * k2m52;
        }
        return;
      }
      m12 += 1;
      m21 += 1;
      lik *= prev_numer / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      cur_eps += 2 * k2m52;
      if (lik < 1 + cur_eps) {
        if (lik <= 1 - cur_eps) {
          break;
        }
        const intptr_t cmp_result = Fisher23Compare(obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, S_CAST(int32_t, m11), S_CAST(int32_t, m12), starting_lnprobv_ddr, &lik);
        cur_eps = 3 * k2m52;
        if (cmp_result <= 0) {
          *tie_ctp += (cmp_result == 0);
          break;
        }
      }
    }
    *base_likp = lik;
    total = lik;
    m11_tail = m11;
    m12_tail = m12;
    m21_tail = m21;
    m22_tail = m22;
  } else {
    m11_tail = m11;
    m12_tail = m12;
    m21_tail = m21;
    m22_tail = m22;
    const double tailenter_lik = lik;
    uint32_t tie_ct_incr = 0;
    while (1) {
      if (lik > 1 - cur_eps) {
        if (lik >= 1 + cur_eps) {
          break;
        }
        const intptr_t cmp_result = Fisher23Compare(obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, S_CAST(int32_t, m11), S_CAST(int32_t, m12), starting_lnprobv_ddr, &lik);
        cur_eps = 3 * k2m52;
        if (cmp_result > 0) {
          break;
        }
        tie_ct_incr += (cmp_result == 0);
      }
      total += lik;
      m11 += 1;
      m22 += 1;
      const double lik_mult = m12 * m21 / (m11 * m22);
      if (lik_mult <= 1 + 2 * k2m52) {
        const int64_t m11_i = S_CAST(int32_t, m11);
        const int64_t m12_i = S_CAST(int32_t, m12);
        const int64_t m21_i = S_CAST(int32_t, m21);
        const int64_t m22_i = S_CAST(int32_t, m22);
        const int64_t lik_mult_numer = m12_i * m21_i;
        const int64_t lik_mult_denom = m11_i * m22_i;
        if (lik_mult_numer <= lik_mult_denom) {
          // Interestingly, this case is well-behaved in a way that isn't true
          // for chrX HWE.
          // Given that the current rank's mode has relative-likelihood <= 1,
          // consider the triangle containing the contingency tables accessible
          // by taking m13-=1 m23+=1 m12+=1 m22-=1 and m13-=1 m23+=1 m11+=1
          // m21-=1 steps from the *previous* rank's mode if m13_decreasing is
          // true, or m13+=1 m23-=1 m11-=1 m21+=1 and m13+=1 m23-=1 m12-=1
          // m22+=1 steps if m13_decreasing is false.
          // - Modes on subsequent ranks must lie within this triangle, because
          //   the likelihood multipliers for moving out of the left or right
          //    edges of the triangle in a later rank are bounded above by the
          //   corresponding likelihood multipliers in the current rank.
          // - Likelihood multipliers for moving to the next rank within the
          //   triangle are less than the larger of the two likelihood
          //   multipliers for moving from the previous rank's mode to an
          //   adjacent table on the current rank; by assumption this value <=
          //   1.
          if (tie_ct_incr) {
            *tie_ctp += tie_ct_incr + (lik_mult_numer == lik_mult_denom);
          }
          *center_is_emptyp = 1;
          return;
        }
      }
      lik *= lik_mult;
      m12 -= 1;
      m21 -= 1;
      cur_eps += 2 * k2m52;
    }
    *base_likp = lik;
    *tie_ctp += tie_ct_incr;
    lik = tailenter_lik;
  }
  *base_epsp = cur_eps;
  *saved_m11p = m11;
  *saved_m12p = m12;
  *saved_m21p = m21;
  *saved_m22p = m22;
  // {m11,m12,m21,m22}_tail is now at a table with relative-likelihood <= 1,
  // and 'total' is the sum of all relative-likelihoods for tables at least as
  // close to the center.

  // sum tail to floating-point precision limit
  m11 = m11_tail;
  m12 = m12_tail;
  m21 = m21_tail;
  m22 = m22_tail;
  while (1) {
    m12 += 1;
    m21 += 1;
    lik *= m11 * m22 / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    const double preadd = total;
    total += lik;
    if (total == preadd) {
      break;
    }
  }
  *totalp = total;
  *center_is_emptyp = 0;
}

static inline void Fisher23LnPRightTailsum(dd_real starting_lnprobv_ddr, uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m13, uint32_t obs_m21, uint32_t obs_m22, uint32_t obs_m23, double* base_likp, double* base_lnlikp, double* base_epsp, double* saved_m11p, double* saved_m12p, double* saved_m21p, double* saved_m22p, int32_t* tie_ctp, double* totalp) {
  uint32_t unused;
  Fisher23LnPLeftTailsum(starting_lnprobv_ddr, obs_m12, obs_m11, obs_m13, obs_m22, obs_m21, obs_m23, base_likp, base_lnlikp, base_epsp, saved_m12p, saved_m11p, saved_m22p, saved_m21p, tie_ctp, totalp, &unused);
}

// obs_m11 + obs_m12 + obs_m13 + obs_m21 + obs_m22 + obs_m23 assumed to be
// <2^31.
double Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp) {
  // Normalize: m11 + m21 >= m12 + m22 >= m13 + m23,
  //            m13 * (m21 + m22) <= m23 * (m11 + m12)
  // Note that columns are reversed from PLINK 1.9 fisher23().
  {
    int32_t mx1i = obs_m11 + obs_m21;
    {
      int32_t mx2i = obs_m12 + obs_m22;
      if (mx1i < mx2i) {
        swap_i32(&obs_m11, &obs_m12);
        swap_i32(&obs_m21, &obs_m22);
        swap_i32(&mx1i, &mx2i);
      }
      {
        int32_t mx3i = obs_m13 + obs_m23;
        if (mx2i < mx3i) {
          swap_i32(&obs_m12, &obs_m13);
          swap_i32(&obs_m22, &obs_m23);
          swap_i32(&mx2i, &mx3i);
        }
        if (mx3i == 0) {
          return Fisher22TwoSidedP(obs_m11, obs_m12, obs_m21, obs_m22, midp, 1);
        }
      }
      if (mx1i < mx2i) {
        swap_i32(&obs_m11, &obs_m12);
        swap_i32(&obs_m21, &obs_m22);
        mx1i = mx2i;
      }
    }
    if (S_CAST(int64_t, obs_m13) * (obs_m21 + obs_m22) > S_CAST(int64_t, obs_m23) * (obs_m11 + obs_m12)) {
      swap_i32(&obs_m11, &obs_m21);
      swap_i32(&obs_m12, &obs_m22);
      swap_i32(&obs_m13, &obs_m23);
    }
  }
  double orig_base_lnlikl = 0;
  double orig_base_lnlikr = 0;
  double orig_base_likl;
  double orig_base_epsl;
  double orig_base_likr;
  double orig_base_epsr;
  double outer_sum;
  dd_real starting_lnprobv_ddr;
  int32_t tie_ct;
  double orig_saved_l11;
  double orig_saved_l12;
  double orig_saved_l21;
  double orig_saved_l22;
  double orig_saved_r11;
  double orig_saved_r12;
  double orig_saved_r21;
  double orig_saved_r22;
  Fisher23LnStartingRank(obs_m11, obs_m12, obs_m21, obs_m22, &outer_sum, &starting_lnprobv_ddr, &tie_ct, &orig_base_likl, &orig_base_lnlikl, &orig_base_epsl, &orig_base_likr, &orig_base_lnlikr, &orig_base_epsr, &orig_saved_l11, &orig_saved_l12, &orig_saved_l21, &orig_saved_l22, &orig_saved_r11, &orig_saved_r12, &orig_saved_r21, &orig_saved_r22);

  // Returned starting_lnprobv_ddr is for a single rank, corresponding to
  // 1 / (obs_m11! obs_m12! obs_m21! obs_m22!).
  // Include m13! and m23! before we iterate over other ranks.
  starting_lnprobv_ddr =
    ddr_sub(starting_lnprobv_ddr,
            ddr_add_lfacts(obs_m13, obs_m23));

  // Other log-factorial expressions we want:
  //
  // * lnprobf, so we can add it to starting_lnprobv at the end and convert
  //   outer_sum (which is a multiple of starting_prob) into the final
  //   log-[mid]p-value.
  //
  //     (m11+m12+m13)! (m21+m22+m23)! (m11+m21)! (m12+m22)! (m13+m23)!
  //     --------------------------------------------------------------
  //                       (m11+m12+m13+m21+m22+m23)!
  //
  // * Likelihood sum over a rank is
  //
  //     (m11+m12+m13)! (m21+m22+m23)! (m13+m23)!    (m11+m12+m21+m22)!
  //     ---------------------------------------- * ---------------------
  //       (m11+m12+m13+m21+m22+m23)! m13! m23!     (m11+m12)! (m21+m22)!
  //
  //   i.e. the 2x2 likelihood after merging columns 1 and 2.
  //
  //   Converting to starting_prob units:
  //
  //                    1                           (m11+m12+m21+m22)!
  //     ------------------------------- * ------------------------------------
  //     m13! m23! (m11+m12)! (m21+m22)!   (m11+m21)! (m12+m22)! starting_probv
  const double m1x = obs_m11 + obs_m12 + obs_m13;
  const double m2x = obs_m21 + obs_m22 + obs_m23;
  const double mx1 = obs_m11 + obs_m21;
  const double mx2 = obs_m12 + obs_m22;
  const double mx3 = obs_m13 + obs_m23;
  const double mxx = m1x + m2x;
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_sort_and_add_5_lfacts(m1x, m2x, mx1, mx2, mx3),
            ddr_lfact(mxx));
  const dd_real rank_relative_lnprobf_ddr =
    ddr_sub(ddr_lfact(mx1 + mx2),
            ddr_add(ddr_add_lfacts(mx1, mx2), starting_lnprobv_ddr));

  // outer_sum, orig_base_likl, and orig_base_likr are expressed as multiples
  // of starting_prob.
  for (uint32_t m13_decreasing = 0; m13_decreasing != 2; ++m13_decreasing) {
    double m13 = obs_m13;
    double m23 = obs_m23;
    // left = small m11 and m22; right = large m11 and m22.
    double l11 = orig_saved_l11;
    double l12 = orig_saved_l12;
    double l21 = orig_saved_l21;
    double l22 = orig_saved_l22;
    double r11 = orig_saved_r11;
    double r12 = orig_saved_r12;
    double r21 = orig_saved_r21;
    double r22 = orig_saved_r22;
    double base_likl = orig_base_likl;
    double base_lnlikl = orig_base_lnlikl;
    double base_epsl = orig_base_epsl;
    double base_likr = orig_base_likr;
    double base_lnlikr = orig_base_lnlikr;
    double base_epsr = orig_base_epsr;
    int32_t m13_decr_last;
    if (m13_decreasing) {
      m13_decr_last = MINV(obs_m13, obs_m21 + obs_m22);
    } else {
      m13_decr_last = -MINV(obs_m23, obs_m11 + obs_m12);
    }
    for (int32_t m13_decr = 0; m13_decr != m13_decr_last; ) {
      m13_decr += m13_decreasing * 2 - 1;
      double lik_mult;
      if (m13_decreasing) {
        m23 += 1;
        // Need to be careful to not move l{11,12,21,22} to the right of the
        // mode, or r{11,12,21,22} to the left of it.
        if (l22 != 0) {
          l12 += 1;
          lik_mult = m13 * l22 / (m23 * l12);
          l22 -= 1;
        } else {
          l11 += 1;
          lik_mult = m13 * l21 / (m23 * l11);
          l21 -= 1;
        }
        m13 -= 1;
      } else {
        m13 += 1;
        if (l11 != 0) {
          l21 += 1;
          lik_mult = m23 * l11 / (m13 * l21);
          l11 -= 1;
        } else {
          l22 += 1;
          lik_mult = m23 * l12 / (m13 * l22);
          l12 -= 1;
        }
        m23 -= 1;
      }
      if (base_likl == 0.0) {
        const double lnlik_incr = log(lik_mult);
        base_lnlikl += lnlik_incr;
        base_epsl += 3 * k2m52;
        if (fabs(lnlik_incr) >= 2) {
          base_epsl += (trunc(fabs(lnlik_incr)) - 1) * k2m52;
        }
      } else {
        base_likl *= lik_mult;
        base_epsl += 2 * k2m52;
      }
      double tail_incr1 = 0.0;
      uint32_t center_is_empty;
      Fisher23LnPLeftTailsum(starting_lnprobv_ddr, obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, &base_likl, &base_lnlikl, &base_epsl, &l11, &l12, &l21, &l22, &tie_ct, &tail_incr1, &center_is_empty);
      if (center_is_empty) {
        // All tables in this rank, and all subsequent ranks, are no more
        // probable than the initial table, and we've already counted ties.
        double m11_12 = m1x - m13;
        double m21_22 = m2x - m23;
        const dd_real rank_adj_ddr = ddr_sort_and_add_4_lfacts(m13, m23, m11_12, m21_22);
        double rank_lik = exp(ddr_sub(rank_relative_lnprobf_ddr, rank_adj_ddr).x[0]);
        if (m13_decreasing) {
          while (1) {
            const double preadd = outer_sum;
            outer_sum += rank_lik;
            if (outer_sum == preadd) {
              break;
            }
            m11_12 += 1;
            m23 += 1;
            rank_lik *= m13 * m21_22 / (m23 * m11_12);
            m13 -= 1;
            m21_22 -= 1;
          }
        } else {
          while (1) {
            const double preadd = outer_sum;
            outer_sum += rank_lik;
            if (outer_sum == preadd) {
              break;
            }
            m13 += 1;
            m21_22 += 1;
            rank_lik *= m11_12 * m23 / (m13 * m21_22);
            m11_12 -= 1;
            m23 -= 1;
          }
        }
        break;
      }
      outer_sum += tail_incr1;
      if (m13_decreasing) {
        const double prev_m13 = m13 + 1;
        if (r21 != 0) {
          r11 += 1;
          lik_mult = prev_m13 * r21 / (m23 * r11);
          r21 -= 1;
        } else {
          r12 += 1;
          lik_mult = prev_m13 * r22 / (m23 * r12);
          r22 -= 1;
        }
      } else {
        const double prev_m23 = m23 + 1;
        if (r12 != 0) {
          r22 += 1;
          lik_mult = prev_m23 * r12 / (m13 * r22);
          r12 -= 1;
        } else {
          r21 += 1;
          lik_mult = prev_m23 * r11 / (m13 * r21);
          r11 -= 1;
        }
      }
      if (base_likr == 0.0) {
        const double lnlik_incr = log(lik_mult);
        base_lnlikr += lnlik_incr;
        base_epsr += 3 * k2m52;
        if (fabs(lnlik_incr) >= 2) {
          base_epsr += (trunc(fabs(lnlik_incr)) - 1) * k2m52;
        }
      } else {
        base_likr *= lik_mult;
        base_epsr += 2 * k2m52;
      }
      double tail_incr2 = 0.0;
      Fisher23LnPRightTailsum(starting_lnprobv_ddr, obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, &base_likr, &base_lnlikr, &base_epsr, &r11, &r12, &r21, &r22, &tie_ct, &tail_incr2);
      outer_sum += tail_incr2;
    }
  }
  const double starting_lnprob =
    ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
  if (midp) {
    outer_sum -= tie_ct * 0.5;
  }
  const double result = log(outer_sum) + starting_lnprob;
  if (result > -k2m35) {
    // true p-value should always be 1 here
    // (possible todo: check boundary cases with total near 2^31)
    return 0.0;
  }
  return result;
}

#ifdef __cplusplus
}
#endif
