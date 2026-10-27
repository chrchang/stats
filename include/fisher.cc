// Fisher's Exact Test library, copyright (C) 2013-2026 Christopher Chang.
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
double Fisher22TwoSidedP(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, int32_t midp, uint32_t logp) {
  // Normalize: m11 >= m22, m12 >= m21, m11*m22 <= m12*m21.
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
  // bugfix (8 Aug 2026): previous 172-steps-from-mode check doesn't prevent
  // overflow for obs_m12 huge, obs_m21 <= 172
  const double first_inward_mult = S_CAST(double, obs_m12) * S_CAST(double, obs_m21) / (S_CAST(double, obs_m11 + 1) * S_CAST(double, obs_m22 + 1));
  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx1 = obs_m11 + obs_m21;
  const double mxx = m1x + m2x;
  const double modal_m21 = m2x * mx1 / mxx;
  if ((first_inward_mult <= 1) || ((m21 - modal_m21) * log(first_inward_mult) < 672)) {
    lik = 1;
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21;
    m22 = obs_m22;
    double center_sum = 0.5 * midp;
    double one_plus_scaled_eps = 1 + k2m52;
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= (m12 * m21) / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      one_plus_scaled_eps += 2 * k2m52;
      if (lik < one_plus_scaled_eps) {
        if (lik <= 2 - one_plus_scaled_eps) {
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
        one_plus_scaled_eps = 1 + 3 * k2m52;
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
  {
    // x=modal_m21 satisfies
    //    (m2x - x) * (mx1 - x) = x * (x + m12_minus_m21)
    // -> (x - m2x) * (x - mx1) = x * (x + m12_minus_m21)
    // -> x^2 + x*(-m2x - mx1) + m2x*mx1 = x^2 + x*m12_minus_m21
    // -> m2x*mx1 = x*(m12_minus_m21 + m2x + mx1)
    // -> x = m2x*mx1 / mxx
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
  const uint32_t tdr_lnprobv_needed = use_tdr_for_hypergeom_lnprob(S_CAST(int64_t, mxx));
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
      }
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
    Fisher22TwoSidedP_rtail_instep:
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
    // bugfix (8 Aug 2026): must take one step when lik in (1 - 3*2^{-52}, 1).
    goto Fisher22TwoSidedP_rtail_instep;
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

// R or!=1 case.
double Fisher22TwoSidedPEx(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, td_real odds_tdr, int32_t midp, uint32_t logp) {
  // Normalize: m11 >= m22, m12 >= m21, m11*m22 <= m12*odds*m21
  if (obs_m11 < obs_m22) {
    swap_i64(&obs_m11, &obs_m22);
  }
  if (obs_m12 < obs_m21) {
    swap_i64(&obs_m12, &obs_m21);
  }
  // Rounding error ok, just need pmf(m11 - 1) < pmf(m11) after this.
  if (ddr_gt(ddr_mul2d(obs_m11, obs_m22), ddr_mul(ddr_mul2d(obs_m12, obs_m21), ddr_make_td(odds_tdr)))) {
    swap_i64(&obs_m11, &obs_m12);
    swap_i64(&obs_m21, &obs_m22);
    odds_tdr = tdr_accurate_div(tdr_make1(1.0), odds_tdr);
  }
  const double odds = odds_tdr.x[0];
  const double obs_m21d = obs_m21;
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21d;
  double m22 = obs_m22;
  const double first_inward_mult = (m12 * m21 * odds) / ((m11 + 1) * (m22 + 1));
  if ((!midp) && (first_inward_mult <= 1 + 3 * k2m52)) {
    // Fast path for p = 1, similar to that in BinomTwoSidedP().
    if (first_inward_mult < 1 - 3 * k2m52) {
      return logp? 0.0 : 1.0;
    }
    const td_real first_inward_mult_numer_tdr = tdr_mul(tdr_make_dd(ddr_mul2d(m12, m21)), odds_tdr);
    const td_real first_inward_mult_denom_tdr = tdr_make_dd(ddr_mul2d(m11 + 1, m22 + 1));
    if (tdr_leq(first_inward_mult_numer_tdr, tdr_mul(first_inward_mult_denom_tdr, tdr_make2(1, 8 * _tdr_eps)))) {
      return logp? 0.0 : 1.0;
    }
    // don't need to update first_inward_mult
  }

  // Iterate outward to floating-point precision limit.
  double lik = 1;
  double tail_sum = 1 - 0.5 * midp;
  while (1) {
    m12 += 1;
    m21 += 1;
    lik *= (m11 * m22) / (m12 * m21 * odds);
    m11 -= 1;
    m22 -= 1;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
  }

  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx1 = obs_m11 + obs_m21;
  const int64_t modal_m11 = ModeFNCHypergeo(S_CAST(int64_t, m1x), S_CAST(int64_t, m2x), S_CAST(int64_t, mx1), odds);
  // If our starting contingency table is close enough to the mode to plausibly
  // contribute to center_sum (instead of being lost to the precision limit),
  // just jump back to the starting table and start summing inwards.
  //
  // For a nonnegative integer n < 2^52, we want to check if the upper-bound
  //   r + r^2 + ... + r^n
  // is <= ~2^106, since 2^53 corresponds to precision limit, and log(pmf(k))
  // is asymptotically quadratic (so actual log(pmf(mode)) - log(pmf(obs_m11))
  // approaches half the tangent-line upper bound).
  // This is obviously true for r <= 1.  For r > 1,
  //   r + r^2 + ... + r^n = r(1 - r^n)/(1-r) <= 2^106
  //   r(r^n - 1)/(r-1) <= 2^106
  //   r^n - 1 <= 2^106 * (r-1)/r
  //   r^n <= 1 + 2^106 * (1 - 1/r)
  //   n log r <= log(1 + 2^106 * (1 - 1/r))
  if ((first_inward_mult <= 1) || (modal_m11 - obs_m11) * log(first_inward_mult) < log1p(((1.0 / k2m53) * (1.0 / k2m53)) * (1 - 1 / first_inward_mult))) {
    lik = 1;
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21d;
    m22 = obs_m22;
    double center_sum = 0.5 * midp;
    double one_plus_scaled_eps = 1 + k2m52;
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= (m12 * m21 * odds) / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      one_plus_scaled_eps += 3 * k2m52;
      if (lik < one_plus_scaled_eps) {
        if (lik <= 2 - one_plus_scaled_eps) {
          tail_sum += lik;
          break;
        }
        const int64_t m22_incr = S_CAST(int64_t, m22) - obs_m22;
        td_real starting_lnprobv_tdr = tdr_make1(DBL_MAX);
        td_real ln_odds_ratio_tdr = tdr_make1(DBL_MAX);
        const intptr_t cmp_result = FNCHypergeoCompare(obs_m11, obs_m12, obs_m21, obs_m22, odds_tdr, m22_incr, &starting_lnprobv_tdr, &ln_odds_ratio_tdr, &lik);
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
    // Ok if m12 or m21 is -1, lik is already 0 in that case
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= (m12 * m21 * odds) / (m11 * m22);
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

  // First jump to the center and evaluate 1/pmf(mode).
  // Then jump to the other tail, using the same approach as in
  // Fisher22TwoSidedP().

  // Save off first-tail m21 stopping point, since it hints at the location of
  // the other tail.
  const double m21_tail_prev = m21;

  m11 = modal_m11;
  m12 = m1x - m11;
  m21 = mx1 - m11;
  m22 = m2x - m21;
  // This value (and consequently starting_lnprob_ddr) is not that accurate.
  // But we can live with that since it isn't needed for near-tie handling, it
  // only contributes to the final (mid)p-value.
  const double mode_neg_log_pmf = log(CentralInvProbFNCHypergeo(m11, m12, m21, m22, odds));

  const uint32_t tdr_lnprobv_needed = use_tdr_for_nchypergeom_lnprob(S_CAST(int64_t, m1x + m2x));
  td_real starting_lnprobv_tdr;
  td_real lnodds_tdr;
  dd_real starting_lnprob_ddr;
  if (!tdr_lnprobv_needed) {
    const dd_real starting_neg_lnprobv_ddr = ddr_sort_and_add_4_lfacts(obs_m11, obs_m12, obs_m21d, obs_m22);
    const dd_real lnodds_ddr = ddr_log(ddr_make_td(odds_tdr));
    dd_real ddrs[5];
    ddrs[0] = ddr_lfact(m11);
    ddrs[1] = ddr_lfact(m12);
    ddrs[2] = ddr_lfact(m21);
    ddrs[3] = ddr_lfact(m22);
    ddrs[4] = ddr_muld(lnodds_ddr, m21 - obs_m21d);
    const dd_real mode_neg_lnprobv_ddr = ddr_sort_and_add(5, ddrs);
    starting_lnprob_ddr = ddr_subd(ddr_sub(mode_neg_lnprobv_ddr, starting_neg_lnprobv_ddr), mode_neg_log_pmf);
    starting_lnprobv_tdr = tdr_make_dd(ddr_negate(starting_neg_lnprobv_ddr));
    lnodds_tdr = tdr_make_dd(lnodds_ddr);
  } else {
    td_real tdrs[5];
    tdrs[0] = tdr_lfact(obs_m11);
    tdrs[1] = tdr_lfact(obs_m12);
    tdrs[2] = tdr_lfact(obs_m21d);
    tdrs[3] = tdr_lfact(obs_m22);
    starting_lnprobv_tdr = tdr_negate(tdr_sort_and_add(4, tdrs));
    lnodds_tdr = tdr_log(odds_tdr);
    tdrs[0] = tdr_lfact(m11);
    tdrs[1] = tdr_lfact(m12);
    tdrs[2] = tdr_lfact(m21);
    tdrs[3] = tdr_lfact(m22);
    tdrs[4] = tdr_muld(lnodds_tdr, m21 - obs_m21d);
    const td_real mode_neg_lnprobv_tdr = tdr_sort_and_add(5, tdrs);
    starting_lnprob_ddr = ddr_make_td(tdr_subd(tdr_add(mode_neg_lnprobv_tdr, starting_lnprobv_tdr), mode_neg_log_pmf));
  }

  {
    const double modal_m21 = mx1 - modal_m11;
    m21 = 2 * modal_m21 - (m21_tail_prev + obs_m21d) * 0.5;
    m21 = S_CAST(int64_t, m21);
    if (m21 < 0) {
      m21 = 0;
    }
  }
  // Extremal case: m11=mx1, m12=m12_minus_m21, m21=0, m22=m2x
  //   log((m12 + 1) * (m21 + 1) * odds / (m11 * m22))
  // = log((m12_minus_m21 + 1) * odds / (mx1 * m2x))
  const double m12_minus_m21 = obs_m12 - obs_m21;
  double lnprobv_diff_min = log((m12_minus_m21 + 1) * odds / (mx1 * m2x)) * (1 + kSmallEpsilon);
  while (1) {
    m11 = mx1 - m21;
    m12 = m12_minus_m21 + m21;
    m22 = m2x - m21;
    double lnprobv_diff;
    if (!tdr_lnprobv_needed) {
      dd_real ddrs[5];
      ddrs[0] = ddr_lfact(m11);
      ddrs[1] = ddr_lfact(m12);
      ddrs[2] = ddr_lfact(m21);
      ddrs[3] = ddr_lfact(m22);
      ddrs[4] = ddr_muld(ddr_make_td(lnodds_tdr), m21 - obs_m21d);
      const dd_real neg_lnprobv_ddr = ddr_sort_and_add(5, ddrs);
      lnprobv_diff = -ddr_add(ddr_make_td(starting_lnprobv_tdr), neg_lnprobv_ddr).x[0];
    } else {
      td_real tdrs[5];
      tdrs[0] = tdr_lfact(m11);
      tdrs[1] = tdr_lfact(m12);
      tdrs[2] = tdr_lfact(m21);
      tdrs[3] = tdr_lfact(m22);
      tdrs[4] = tdr_muld(lnodds_tdr, m21 - obs_m21d);
      const td_real neg_lnprobv_tdr = tdr_sort_and_add(5, tdrs);
      lnprobv_diff = -tdr_add(starting_lnprobv_tdr, neg_lnprobv_tdr).x[0];
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
      const double ll_deriv = log(m12 * m21 * odds / ((m11 + 1) * (m22 + 1)));
      m21 -= ceil(-lnprobv_diff / ll_deriv);
      if (m21 < 0) {
        m21 = 0;
      }
    } else {
      double ll_deriv = DBL_MAX;
      if ((lnprobv_diff_min < -53 * kLn2) && (m21 > 0)) {
        // Tighten this threshold if that lets us sum fewer terms later.
        ll_deriv = log((m12 + 1) * (m21 + 1) * odds / (m11 * m22));
        lnprobv_diff_min = ll_deriv * (1 + kSmallEpsilon);
      }
      if (lnprobv_diff > lnprobv_diff_min) {
        lik = exp(lnprobv_diff);
        break;
      }
      if (ll_deriv == DBL_MAX) {
        ll_deriv = log((m12 + 1) * (m21 + 1) * odds / (m11 * m22));
        m21 += S_CAST(int64_t, lnprobv_diff / ll_deriv);
      }
    }
  }
  // Sum toward center, until lik >= 1.
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  const double tailenter_lik = lik;
  const double tailenter_m11 = m11;
  const double tailenter_m12 = m12;
  const double tailenter_m21 = m21;
  const double tailenter_m22 = m22;
  while (1) {
    while (lik <= one_minus_scaled_eps) {
    Fisher22TwoSidedPEx_rtail_instep:
      tail_sum += lik;
      m12 += 1;
      m21 += 1;
      lik *= m11 * m22 / (m12 * m21 * odds);
      m11 -= 1;
      m22 -= 1;
      one_minus_scaled_eps -= 3 * k2m52;
    }
    if (lik >= 2 - one_minus_scaled_eps) {
      break;
    }
    const int64_t m22_incr = S_CAST(int64_t, m22) - obs_m22;
    if (!tdr_lnprobv_needed) {
      starting_lnprobv_tdr.x[2] = DBL_MAX;
      lnodds_tdr.x[2] = DBL_MAX;
    }
    const intptr_t cmp_result = FNCHypergeoCompare(obs_m11, obs_m12, obs_m21, obs_m22, odds_tdr, m22_incr, &starting_lnprobv_tdr, &lnodds_tdr, &lik);
    if (cmp_result >= 0) {
      if (cmp_result == 0) {
        tail_sum += 1 - 0.5 * midp;
      }
      break;
    }
    one_minus_scaled_eps = 1 - 3 * k2m52;
    goto Fisher22TwoSidedPEx_rtail_instep;
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
    lik *= m12 * m21 * odds / (m11 * m22);
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
    dd_real mean_ddr;
    double variance;
    MeanAndVarianceFNCHypergeo(m1, m2, n, odds, &mean_ddr, &variance);
    if (variance == 0.0) {
      // defensive, this shouldn't actually happen?
      return odds;
    }
    const double mean_delta = -ddr_subd(mean_ddr, m11d).x[0];
    const double odds_incr = mean_delta * odds / variance;
    // printf("%.17g: %.17g %.17g %.17g\n", odds, mean_delta, variance, odds_incr);
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
    const double first_inward_mult = S_CAST(double, obs_m11) * S_CAST(double, obs_m22) / ((S_CAST(double, obs_m12) + 1) * (S_CAST(double, obs_m21) + 1));
    if ((first_inward_mult <= 1) || (delta * log(first_inward_mult) < 612)) {
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
              *orig_base_epsl_ptr = 0;  // unused in log representation
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
  const double first_inward_mult = S_CAST(double, S_CAST(int64_t, obs_m12) * obs_m21) / S_CAST(double, (obs_m11 + 1LL) * (obs_m22 + 1LL));
  if ((first_inward_mult <= 1) || (delta * log(first_inward_mult) > -612)) {
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
            *orig_base_epsr_ptr = 0;  // unused in log representation
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

// Log-likelihood of the current table relative to the observed table.
static dd_real Fisher23LnLik(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m13, uint32_t obs_m21, uint32_t obs_m22, uint32_t obs_m23, uint32_t cur_m11, uint32_t cur_m12, dd_real starting_lnprobv_ddr) {
  const uint32_t cur_m13 = obs_m11 + obs_m12 + obs_m13 - cur_m11 - cur_m12;
  double lfact_args[6];
  lfact_args[0] = cur_m11;
  lfact_args[1] = cur_m12;
  lfact_args[2] = cur_m13;
  lfact_args[3] = obs_m11 + obs_m21 - cur_m11;
  lfact_args[4] = obs_m12 + obs_m22 - cur_m12;
  lfact_args[5] = obs_m13 + obs_m23 - cur_m13;
  return ddr_sub(ddr_negate(ddr_sort_and_add_lfacts(6, lfact_args)), starting_lnprobv_ddr);
}

// 'Left' = small m11 and m22.
// In the log-likelihood representation (*base_likp == 0), *base_lnlikp is
// only used to decide when to switch back to the regular representation; no
// error bound is tracked for it.  Instead, the likelihood is recomputed from
// scratch on the switch, so accumulated rounding error does not carry over.
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
    while (last_lnp >= kLnSwitchThresh) {
      const double prev_numer = m11 * m22;
      if (prev_numer == 0) {
        // lowest-likelihood table on this side is still too probable
        *base_lnlikp = last_lnp;
        *center_is_emptyp = 0;
        *saved_m11p = m11;
        *saved_m12p = m12;
        *saved_m21p = m21;
        *saved_m22p = m22;
        return;
      }
      m12 += 1;
      m21 += 1;
      last_lnp += log(prev_numer / (m12 * m21));
      m11 -= 1;
      m22 -= 1;
    }
    lik = ddr_exp(Fisher23LnLik(obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, S_CAST(int32_t, m11), S_CAST(int32_t, m12), starting_lnprobv_ddr)).x[0];
    cur_eps = 2 * k2m52;
  }
  double m11_tail;
  double m12_tail;
  double m21_tail;
  double m22_tail;
  uint32_t entry_above = (lik >= 1 + cur_eps);
  uint32_t entry_tie = 0;
  if ((!entry_above) && (lik > 1 - cur_eps)) {
    // Resolve a near-tie at the entry table before choosing a direction.  If
    // the entry table is slightly more probable than the observed table, the
    // inward-walking branch would stop immediately and then sum outward from
    // it, wrongly counting any adjacent tables that are also slightly more
    // probable.
    const intptr_t cmp_result = Fisher23Compare(obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, S_CAST(int32_t, m11), S_CAST(int32_t, m12), starting_lnprobv_ddr, &lik);
    cur_eps = 3 * k2m52;
    entry_above = (cmp_result > 0);
    entry_tie = (cmp_result == 0);
  }
  if (entry_above) {
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
          *base_lnlikp = log(lik);
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
    uint32_t tie_ct_incr = entry_tie;
    // Entry table is already known to be no more probable than the observed
    // table.
    goto Fisher23LnPLeftTailsum_inward_entry;
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
    Fisher23LnPLeftTailsum_inward_entry:
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
        base_lnlikl += log(lik_mult);
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
        base_lnlikr += log(lik_mult);
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

// 2xk generalization of Fisher23LnP().
//
// Columns are sorted in decreasing-total order.  The two largest columns
// (indexes 0 and 1) form the base 2x2 ranks, which are handled in essentially
// the same way as Fisher23LnP(): we track the positions of the left and right
// tail boundaries as we move between adjacent ranks, so we only need to visit
// the tables near them.  Column 2 plays the role of Fisher23LnP()'s third
// column.  Each additional column j>2 adds a level of recursion: for each
// possible row-1 value of column j, we evaluate the 2xj subproblem.
//
// The important differences from Fisher23LnP() are:
// 1. At every level, iteration starts from the modal rank of the current
//    subproblem, instead of from the rank containing the observed table.
//    Over a single level, the maximum relative-likelihood across a rank is a
//    log-concave function of the rank index (it's a max-convolution of
//    log-concave sequences).  So, starting from the modal rank, once we
//    encounter a rank where no table is more probable than the observed
//    table, the same is true of all subsequent ranks in that direction, and
//    those ranks can be summed in closed form via Vandermonde's identity.
//    (This also implies that ties in those ranks can only occur at the rank
//    containing the first such rank's mode(s).)
// 2. The mode of a 2xj subproblem with fixed margins can be maintained with a
//    single greedy step when the row-1 total changes by 1, since the
//    likelihood is a separable log-concave function of the row-1 entries.
// 3. Since we don't start at the observed table, the relative-likelihood of
//    each 2x3 subproblem's mode is computed from scratch (with dd_real
//    log-factorials).  This also lets us avoid accumulating error when
//    walking long distances in log-space.
// 4. The contribution of a nonempty rank only depends on the multiset of
//    row-1 values in each run of equal-total columns that have been fixed so
//    far.  So for columns with index >= 2 in such a run, a nonempty rank is
//    only evaluated when the run's row-1 values fixed so far are in
//    nonincreasing order (lower index >= higher index), and its contribution
//    is multiplied by a factor which accounts for the skipped permutations;
//    for other nonempty ranks, we only check the mode to determine whether to
//    keep going (or, for column 2, only update the tail boundaries).  Every
//    permutation of a nonempty rank is also nonempty, and all nonempty ranks
//    are visited, so this is exact.  Empty ranks are still handled normally,
//    since their permutations may be summed in closed form instead of
//    visited.  (Symmetry between columns 0 and 1 isn't exploited.)

typedef struct Fisher2xkCtxStruct {
  uint32_t col_ct;
  const int32_t* col_tots;
  const int32_t* obs_r1;
  // cum_col_tots[j] = col_tots[0] + ... + col_tots[j-1]
  const int32_t* cum_col_tots;
  // inner_const_ddrs[j] = lfact(cum_col_tots[j]) - lfact(col_tots[0]) - ... -
  //                       lfact(col_tots[j-1])
  const dd_real* inner_const_ddrs;
  // Row-1 entries of the current contingency table.
  int32_t* cur_r1;
  // Level j saves a copy of the inner subproblem's mode at
  // saved_modes[j * col_ct].
  int32_t* saved_modes;
  // If column j has index >= 2 and is in a run of equal-total columns with
  // index >= 2, and isn't the highest-index column in that run, sym_ns[j] is
  // (highest index in run) - j + 1.  Otherwise it's 0.
  const uint32_t* sym_ns;
  uint64_t* fac_args_buf;
  double* lfact_args_buf;
  dd_real starting_lnprobv_ddr;
  td_real starting_lnprobv_tdr;
  // Normally nullptr.  When Fisher2xkLevel2() is used on a subproblem of a
  // larger table (see Fisher33Base()), exact comparisons need to cover the
  // entire table: full_obs[0..(full_cell_ct-1)] are then the observed table's
  // cells, fixed_cells[] are the current table's cells outside the
  // subproblem (there are full_cell_ct - 2 * col_ct of them), and
  // starting_lnprobv_tdr corresponds to full_obs[].  (starting_lnprobv_ddr is
  // adjusted to account for fixed_cells[] instead.)
  const int32_t* full_obs;
  const int32_t* fixed_cells;
  uint32_t full_cell_ct;
  // Sum of relative-likelihoods of all tables no more probable than the
  // observed table, expressed as a multiple of the observed table's
  // probability.
  double outer_sum;
  double tie_ct;
} Fisher2xkCtx;

// Returns index of the column whose row-1 entry should be incremented to
// maximize the likelihood, setting *numer_ptr and *denom_ptr to the
// likelihood multiplier.  Returns UINT32_MAX if no row-1 entry can be
// incremented.
static uint32_t Fisher2xkBestIncr(const int32_t* col_tots, const int32_t* x, uint32_t col_ct, int64_t* numer_ptr, int64_t* denom_ptr) {
  uint32_t best_idx = UINT32_MAX;
  int64_t best_numer = 0;
  int64_t best_denom = 1;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    const int64_t numer = col_tots[col_idx] - x[col_idx];
    const int64_t denom = x[col_idx] + 1;
    if (numer * best_denom > best_numer * denom) {
      best_idx = col_idx;
      best_numer = numer;
      best_denom = denom;
    }
  }
  *numer_ptr = best_numer;
  *denom_ptr = best_denom;
  return best_idx;
}

// Decrement counterpart of Fisher2xkBestIncr().
static uint32_t Fisher2xkBestDecr(const int32_t* col_tots, const int32_t* x, uint32_t col_ct, int64_t* numer_ptr, int64_t* denom_ptr) {
  uint32_t best_idx = UINT32_MAX;
  int64_t best_numer = 0;
  int64_t best_denom = 1;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    const int64_t numer = x[col_idx];
    const int64_t denom = col_tots[col_idx] - x[col_idx] + 1;
    if (numer * best_denom > best_numer * denom) {
      best_idx = col_idx;
      best_numer = numer;
      best_denom = denom;
    }
  }
  *numer_ptr = best_numer;
  *denom_ptr = best_denom;
  return best_idx;
}

// Given a modal table, this produces a modal table with row-1 total increased
// by 1.
static inline void Fisher2xkModeIncr(const int32_t* col_tots, uint32_t col_ct, int32_t* x) {
  int64_t numer;
  int64_t denom;
  x[Fisher2xkBestIncr(col_tots, x, col_ct, &numer, &denom)] += 1;
}

static inline void Fisher2xkModeDecr(const int32_t* col_tots, uint32_t col_ct, int32_t* x) {
  int64_t numer;
  int64_t denom;
  x[Fisher2xkBestDecr(col_tots, x, col_ct, &numer, &denom)] -= 1;
}

static void Fisher2xkInitMode(const int32_t* col_tots, uint32_t col_ct, int64_t tot, int32_t row1_tot, int32_t* x) {
  // Start from (rounded-down) proportional allocation; this is within a few
  // steps of the mode.
  int32_t cur_row1_tot = 0;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    x[col_idx] = S_CAST(int32_t, (S_CAST(int64_t, col_tots[col_idx]) * row1_tot) / tot);
    cur_row1_tot += x[col_idx];
  }
  for (; cur_row1_tot != row1_tot; ++cur_row1_tot) {
    Fisher2xkModeIncr(col_tots, col_ct, x);
  }
  // Since the likelihood is a separable log-concave function, a table is
  // modal iff no single-unit exchange between columns increases it.  The
  // best exchange always pairs the best increment with the best decrement
  // (if they're in the same column, no exchange is an improvement).
  while (1) {
    int64_t incr_numer;
    int64_t incr_denom;
    int64_t decr_numer;
    int64_t decr_denom;
    const uint32_t incr_idx = Fisher2xkBestIncr(col_tots, x, col_ct, &incr_numer, &incr_denom);
    const uint32_t decr_idx = Fisher2xkBestDecr(col_tots, x, col_ct, &decr_numer, &decr_denom);
    if ((incr_idx == UINT32_MAX) || (decr_idx == UINT32_MAX) || (incr_idx == decr_idx) || (incr_numer * decr_numer <= incr_denom * decr_denom)) {
      return;
    }
    x[incr_idx] += 1;
    x[decr_idx] -= 1;
  }
}

// Number of modal tables, given that x[] is one of them.
static double Fisher2xkModeCt(const int32_t* col_tots, uint32_t col_ct, const int32_t* x) {
  int64_t incr_numer;
  int64_t incr_denom;
  int64_t decr_numer;
  int64_t decr_denom;
  const uint32_t incr_idx = Fisher2xkBestIncr(col_tots, x, col_ct, &incr_numer, &incr_denom);
  const uint32_t decr_idx = Fisher2xkBestDecr(col_tots, x, col_ct, &decr_numer, &decr_denom);
  if ((incr_idx == UINT32_MAX) || (decr_idx == UINT32_MAX) || (incr_numer * decr_numer != incr_denom * decr_denom)) {
    return 1;
  }
  // The other modes are those reachable by simultaneously incrementing s
  // best-increment columns and decrementing s best-decrement columns.  (The
  // two column sets are disjoint.)
  uint32_t incr_ct = 0;
  uint32_t decr_ct = 0;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    const int64_t cur_x = x[col_idx];
    const int64_t cur_col_tot = col_tots[col_idx];
    incr_ct += ((cur_col_tot - cur_x) * incr_denom == incr_numer * (cur_x + 1));
    decr_ct += (cur_x * decr_denom == decr_numer * (cur_col_tot - cur_x + 1));
  }
  // sum_s (incr_ct choose s) * (decr_ct choose s) = (incr_ct + decr_ct)
  // choose incr_ct.
  double result = 1;
  for (uint32_t uii = 1; uii <= incr_ct; ++uii) {
    result = (result * (decr_ct + uii)) / uii;
  }
  return result;
}

// Log-likelihood of the current table relative to the observed table.
static dd_real Fisher2xkLnLik(Fisher2xkCtx* ctxp) {
  const uint32_t col_ct = ctxp->col_ct;
  const int32_t* col_tots = ctxp->col_tots;
  const int32_t* cur_r1 = ctxp->cur_r1;
  double* lfact_args = ctxp->lfact_args_buf;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    lfact_args[2 * col_idx] = cur_r1[col_idx];
    lfact_args[2 * col_idx + 1] = col_tots[col_idx] - cur_r1[col_idx];
  }
  return ddr_sub(ddr_negate(ddr_sort_and_add_lfacts(2 * col_ct, lfact_args)), ctxp->starting_lnprobv_ddr);
}

// Exact comparison of the current table's likelihood with the observed
// table's.
static intptr_t Fisher2xkCompare(Fisher2xkCtx* ctxp, double* dbl_ptr) {
  const uint32_t col_ct = ctxp->col_ct;
  const int32_t* col_tots = ctxp->col_tots;
  const int32_t* obs_r1 = ctxp->obs_r1;
  const int32_t* cur_r1 = ctxp->cur_r1;
  const uint32_t full_cell_ct = ctxp->full_cell_ct;
  const uint32_t arg_ct = full_cell_ct? full_cell_ct : (2 * col_ct);
  uint64_t* numer_factorial_args = ctxp->fac_args_buf;
  uint64_t* denom_factorial_args = &(ctxp->fac_args_buf[arg_ct]);
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    numer_factorial_args[2 * col_idx] = obs_r1[col_idx];
    numer_factorial_args[2 * col_idx + 1] = col_tots[col_idx] - obs_r1[col_idx];
    denom_factorial_args[2 * col_idx] = cur_r1[col_idx];
    denom_factorial_args[2 * col_idx + 1] = col_tots[col_idx] - cur_r1[col_idx];
  }
  if (full_cell_ct) {
    const int32_t* full_obs = ctxp->full_obs;
    const int32_t* fixed_cells = ctxp->fixed_cells;
    for (uint32_t cell_idx = 0; cell_idx != full_cell_ct; ++cell_idx) {
      numer_factorial_args[cell_idx] = full_obs[cell_idx];
    }
    for (uint32_t cell_idx = 2 * col_ct; cell_idx != full_cell_ct; ++cell_idx) {
      denom_factorial_args[cell_idx] = fixed_cells[cell_idx - 2 * col_ct];
    }
  }
  td_real ln_odds_ratio_tdr = tdr_make1(0.0);
  return CompareFactorialProducts(arg_ct, tdr_make1(1.0), 0, 0, numer_factorial_args, denom_factorial_args, &ctxp->starting_lnprobv_tdr, &ln_odds_ratio_tdr, dbl_ptr);
}

// m11 and m12 are base-rank row-1 entries in Fisher2xkBaseTailsum()'s frame of
// reference.
static inline void Fisher2xkSetBase(Fisher2xkCtx* ctxp, uint32_t swap01, double m11, double m12) {
  ctxp->cur_r1[swap01] = S_CAST(int32_t, m11);
  ctxp->cur_r1[1 - swap01] = S_CAST(int32_t, m12);
}

// Adds the closed-form sum of ranks x_j=t, t+dir, t+2*dir, ... at level j,
// where all tables in these ranks are known to be no more probable than the
// observed table.
static void Fisher2xkClosedFormSum(Fisher2xkCtx* ctxp, uint32_t level, int32_t row1_tot, int32_t start_x, int32_t dir) {
  const uint32_t col_ct = ctxp->col_ct;
  const int32_t* col_tots = ctxp->col_tots;
  const int32_t* cur_r1 = ctxp->cur_r1;
  double* lfact_args = ctxp->lfact_args_buf;
  // Sum of likelihoods of all 2x<level> tables with column totals
  // col_tots[0..(level-1)] and row-1 total r is
  //   (C choose r) / (col_tots[0]! ... col_tots[level-1]!)
  // where C = cum_col_tots[level].
  const double inner_col_tot = ctxp->cum_col_tots[level];
  const double col_tot = col_tots[level];
  double x = start_x;
  double inner_r1 = row1_tot - start_x;
  uint32_t lfact_arg_ct = 0;
  for (uint32_t col_idx = level + 1; col_idx != col_ct; ++col_idx) {
    lfact_args[lfact_arg_ct++] = cur_r1[col_idx];
    lfact_args[lfact_arg_ct++] = col_tots[col_idx] - cur_r1[col_idx];
  }
  lfact_args[lfact_arg_ct++] = x;
  lfact_args[lfact_arg_ct++] = col_tot - x;
  lfact_args[lfact_arg_ct++] = inner_r1;
  lfact_args[lfact_arg_ct++] = inner_col_tot - inner_r1;
  const dd_real rank_lnlik_ddr =
    ddr_sub(ddr_sub(ctxp->inner_const_ddrs[level], ctxp->starting_lnprobv_ddr),
            ddr_sort_and_add_lfacts(lfact_arg_ct, lfact_args));
  double rank_lik = exp(rank_lnlik_ddr.x[0]);
  double outer_sum = ctxp->outer_sum;
  while (1) {
    const double preadd = outer_sum;
    outer_sum += rank_lik;
    // These multipliers are zero when we run out of tables.
    double rank_lik_mult;
    if (dir > 0) {
      rank_lik_mult = (col_tot - x) * inner_r1 / ((x + 1) * (inner_col_tot - inner_r1 + 1));
      x += 1;
      inner_r1 -= 1;
    } else {
      rank_lik_mult = x * (inner_col_tot - inner_r1) / ((col_tot - x + 1) * (inner_r1 + 1));
      x -= 1;
      inner_r1 += 1;
    }
    // Rank sums are log-concave, so we're done once they're decreasing and
    // negligible.
    if ((outer_sum == preadd) && (rank_lik_mult <= 1)) {
      break;
    }
    rank_lik *= rank_lik_mult;
  }
  ctxp->outer_sum = outer_sum;
}

// Generalization of Fisher23LnPLeftTailsum().  Differences:
// - Exact comparisons are performed against the full current 2xk table.
// - When swap01 is set, base columns 0 and 1 are swapped (so this computes
//   the right tail).
// - When entry_above is set, the saved table is known to be more probable
//   than the observed table.
// - When skip_sum is set, the tail boundary is updated but the tail is not
//   summed (*totalp is set to 0).
// - When a log-likelihood is converted back to a likelihood, it is recomputed
//   from scratch, so error does not accumulate.
// - When the saved table's likelihood is within epsilon of the observed
//   table's, the near-tie is resolved before choosing a direction.
static void Fisher2xkBaseTailsum(Fisher2xkCtx* ctxp, uint32_t swap01, uint32_t entry_above, uint32_t skip_sum, double* base_likp, double* base_lnlikp, double* base_epsp, double* saved_m11p, double* saved_m12p, double* saved_m21p, double* saved_m22p, double* totalp, uint32_t* center_is_emptyp) {
  double total = 0;
  double lik = *base_likp;
  double cur_eps = *base_epsp;
  double m11 = *saved_m11p;
  double m12 = *saved_m12p;
  double m21 = *saved_m21p;
  double m22 = *saved_m22p;
  *center_is_emptyp = 0;
  if (lik == 0.0) {
    double last_lnp = *base_lnlikp;
    while (last_lnp >= kLnSwitchThresh) {
      const double prev_numer = m11 * m22;
      if (prev_numer == 0) {
        // lowest-likelihood table on this side is still too probable
        *base_lnlikp = last_lnp;
        *saved_m11p = m11;
        *saved_m12p = m12;
        *saved_m21p = m21;
        *saved_m22p = m22;
        *totalp = 0.0;
        return;
      }
      m12 += 1;
      m21 += 1;
      last_lnp += log(prev_numer / (m12 * m21));
      m11 -= 1;
      m22 -= 1;
    }
    Fisher2xkSetBase(ctxp, swap01, m11, m12);
    lik = ddr_exp(Fisher2xkLnLik(ctxp)).x[0];
    cur_eps = 2 * k2m52;
  }
  uint32_t entry_tie = 0;
  if ((!entry_above) && (lik < 1 + cur_eps)) {
    if (lik > 1 - cur_eps) {
      Fisher2xkSetBase(ctxp, swap01, m11, m12);
      const intptr_t cmp_result = Fisher2xkCompare(ctxp, &lik);
      cur_eps = 3 * k2m52;
      entry_above = (cmp_result > 0);
      entry_tie = (cmp_result == 0);
    }
  } else {
    entry_above = 1;
  }
  double m11_tail;
  double m12_tail;
  double m21_tail;
  double m22_tail;
  if (entry_above) {
    while (1) {
      const double prev_numer = m11 * m22;
      if (prev_numer == 0) {
        // lowest-likelihood table on this side is still too probable
        *saved_m11p = m11;
        *saved_m12p = m12;
        *saved_m21p = m21;
        *saved_m22p = m22;
        if (lik < kSwitchThresh) {
          *base_likp = lik;
          *base_epsp = cur_eps;
        } else {
          *base_likp = 0;
          *base_lnlikp = log(lik);
        }
        *totalp = 0.0;
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
        Fisher2xkSetBase(ctxp, swap01, m11, m12);
        const intptr_t cmp_result = Fisher2xkCompare(ctxp, &lik);
        cur_eps = 3 * k2m52;
        if (cmp_result <= 0) {
          ctxp->tie_ct += (cmp_result == 0);
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
    uint32_t tie_ct_incr = entry_tie;
    while (1) {
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
          // Mode of this rank is no more probable than the observed table.
          if (tie_ct_incr) {
            ctxp->tie_ct += tie_ct_incr + (lik_mult_numer == lik_mult_denom);
          }
          *center_is_emptyp = 1;
          return;
        }
      }
      lik *= lik_mult;
      m12 -= 1;
      m21 -= 1;
      cur_eps += 2 * k2m52;
      if (lik > 1 - cur_eps) {
        if (lik >= 1 + cur_eps) {
          break;
        }
        Fisher2xkSetBase(ctxp, swap01, m11, m12);
        const intptr_t cmp_result = Fisher2xkCompare(ctxp, &lik);
        cur_eps = 3 * k2m52;
        if (cmp_result > 0) {
          break;
        }
        tie_ct_incr += (cmp_result == 0);
      }
    }
    *base_likp = lik;
    ctxp->tie_ct += tie_ct_incr;
    lik = tailenter_lik;
  }
  *base_epsp = cur_eps;
  *saved_m11p = m11;
  *saved_m12p = m12;
  *saved_m21p = m21;
  *saved_m22p = m22;
  if (skip_sum) {
    *totalp = 0.0;
    return;
  }

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
}

// Return values for Fisher2xkLevel{2}().
// - kFisher2xkNonempty: at least one table in the subproblem is more probable
//   than the observed table.  outer_sum and tie_ct have been updated.
// - kFisher2xkEmpty: the subproblem's mode is less probable than the observed
//   table.  outer_sum and tie_ct have not been updated.
// - kFisher2xkEmptyTie: the subproblem's mode is exactly as probable as the
//   observed table.  outer_sum and tie_ct have not been updated.
static const uint32_t kFisher2xkNonempty = 0;
static const uint32_t kFisher2xkEmpty = 1;
static const uint32_t kFisher2xkEmptyTie = 2;

// Compares the mode of the current subproblem (whose row-1 entries must
// already be in cur_r1[]) with the observed table, returning the status
// Fisher2xkLevel() would return for that subproblem.
static uint32_t Fisher2xkModeStatus(Fisher2xkCtx* ctxp) {
  const double lnlik = Fisher2xkLnLik(ctxp).x[0];
  if (lnlik < -kSmallEpsilon) {
    return kFisher2xkEmpty;
  }
  if (lnlik <= kSmallEpsilon) {
    double unused;
    const intptr_t cmp_result = Fisher2xkCompare(ctxp, &unused);
    if (cmp_result <= 0) {
      return (cmp_result == 0)? kFisher2xkEmptyTie : kFisher2xkEmpty;
    }
  }
  return kFisher2xkNonempty;
}

// Multiplier for the contribution of a nonempty rank at the given level,
// based on cur_r1[level] and the row-1 values of the higher-index columns in
// the same equal-total run (which are already in nonincreasing order).
// Returns 0 if cur_r1[level] is smaller than cur_r1[level + 1] (so the rank's
// contribution is accounted for by a permutation), and otherwise n/m, where n
// is the number of run columns with index >= level and m is the number of
// them with row-1 value equal to cur_r1[level].  The product of these
// multipliers over a run is the number of distinct permutations of the run's
// row-1 values.
static inline double Fisher2xkSymFactor(const Fisher2xkCtx* ctxp, uint32_t level) {
  const uint32_t n = ctxp->sym_ns[level];
  if (!n) {
    return 1;
  }
  const int32_t* x = &(ctxp->cur_r1[level]);
  if (x[0] < x[1]) {
    return 0;
  }
  uint32_t m = 1;
  while ((m != n) && (x[m] == x[0])) {
    ++m;
  }
  return S_CAST(double, n) / S_CAST(double, m);
}

// On entry, cur_r1[0..2] must be the mode of the 2x3 subproblem with row-1
// total row1_tot and the current cur_r1[3..] values.  It is restored on exit.
static uint32_t Fisher2xkLevel2(Fisher2xkCtx* ctxp, int32_t row1_tot) {
  int32_t* cur_r1 = ctxp->cur_r1;
  const int32_t* col_tots = ctxp->col_tots;
  const int32_t mode0 = cur_r1[0];
  const int32_t mode1 = cur_r1[1];
  const int32_t mode2 = cur_r1[2];
  double orig_base_lik;
  double orig_base_lnlik = 0.0;
  double orig_base_eps;
  {
    const dd_real lnlik_ddr = Fisher2xkLnLik(ctxp);
    const double lnlik = lnlik_ddr.x[0];
    if (lnlik < -kSmallEpsilon) {
      return kFisher2xkEmpty;
    }
    if (lnlik <= kSmallEpsilon) {
      const intptr_t cmp_result = Fisher2xkCompare(ctxp, &orig_base_lik);
      if (cmp_result <= 0) {
        return (cmp_result == 0)? kFisher2xkEmptyTie : kFisher2xkEmpty;
      }
      orig_base_eps = 3 * k2m52;
    } else if (lnlik < kLnSwitchThresh) {
      orig_base_lik = ddr_exp(lnlik_ddr).x[0];
      orig_base_eps = 2 * k2m52;
    } else {
      orig_base_lik = 0.0;
      orig_base_lnlik = lnlik;
      orig_base_eps = 0.0;
    }
  }
  // Mode of the starting rank is more probable than the observed table.
  // Find its tail boundaries by walking outward from the mode.
  // As in Fisher23LnP(), "left" = small m11 and m22.
  double orig_saved_l11 = mode0;
  double orig_saved_l12 = mode1;
  double orig_saved_l21 = col_tots[0] - mode0;
  double orig_saved_l22 = col_tots[1] - mode1;
  double orig_saved_r11 = orig_saved_l11;
  double orig_saved_r12 = orig_saved_l12;
  double orig_saved_r21 = orig_saved_l21;
  double orig_saved_r22 = orig_saved_l22;
  double orig_base_likl = orig_base_lik;
  double orig_base_lnlikl = orig_base_lnlik;
  double orig_base_epsl = orig_base_eps;
  double orig_base_likr = orig_base_lik;
  double orig_base_lnlikr = orig_base_lnlik;
  double orig_base_epsr = orig_base_eps;
  {
    // See Fisher2xkSymFactor().  Tail boundaries must be tracked even when
    // the rank's contribution is skipped.
    const double weight = Fisher2xkSymFactor(ctxp, 2);
    const double saved_tie_ct = ctxp->tie_ct;
    uint32_t center_is_empty;
    double tail_incr1 = 0.0;
    Fisher2xkBaseTailsum(ctxp, 0, 1, weight == 0, &orig_base_likl, &orig_base_lnlikl, &orig_base_epsl, &orig_saved_l11, &orig_saved_l12, &orig_saved_l21, &orig_saved_l22, &tail_incr1, &center_is_empty);
    double tail_incr2 = 0.0;
    Fisher2xkBaseTailsum(ctxp, 1, 1, weight == 0, &orig_base_likr, &orig_base_lnlikr, &orig_base_epsr, &orig_saved_r12, &orig_saved_r11, &orig_saved_r22, &orig_saved_r21, &tail_incr2, &center_is_empty);
    ctxp->outer_sum += weight * (tail_incr1 + tail_incr2);
    ctxp->tie_ct = saved_tie_ct + weight * (ctxp->tie_ct - saved_tie_ct);
  }

  for (uint32_t m13_decreasing = 0; m13_decreasing != 2; ++m13_decreasing) {
    double m13 = mode2;
    double m23 = col_tots[2] - mode2;
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
      m13_decr_last = MINV(mode2, col_tots[0] + col_tots[1] - mode0 - mode1);
    } else {
      m13_decr_last = -MINV(col_tots[2] - mode2, mode0 + mode1);
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
      cur_r1[2] = S_CAST(int32_t, m13);
      if (base_likl == 0.0) {
        // Exact value is recomputed before it matters.
        base_lnlikl += log(lik_mult);
      } else {
        base_likl *= lik_mult;
        base_epsl += 2 * k2m52;
      }
      const double weight = Fisher2xkSymFactor(ctxp, 2);
      const double saved_tie_ct = ctxp->tie_ct;
      double tail_incr1 = 0.0;
      uint32_t center_is_empty;
      Fisher2xkBaseTailsum(ctxp, 0, 0, weight == 0, &base_likl, &base_lnlikl, &base_epsl, &l11, &l12, &l21, &l22, &tail_incr1, &center_is_empty);
      if (center_is_empty) {
        // Empty ranks (including any ties counted for this one) are never
        // weighted.
        // All tables in this rank, and all subsequent ranks, are no more
        // probable than the observed table, and we've already counted ties.
        Fisher2xkClosedFormSum(ctxp, 2, row1_tot, S_CAST(int32_t, m13), m13_decreasing? -1 : 1);
        break;
      }
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
        base_lnlikr += log(lik_mult);
      } else {
        base_likr *= lik_mult;
        base_epsr += 2 * k2m52;
      }
      double tail_incr2 = 0.0;
      Fisher2xkBaseTailsum(ctxp, 1, 0, weight == 0, &base_likr, &base_lnlikr, &base_epsr, &r12, &r11, &r22, &r21, &tail_incr2, &center_is_empty);
      ctxp->outer_sum += weight * (tail_incr1 + tail_incr2);
      ctxp->tie_ct = saved_tie_ct + weight * (ctxp->tie_ct - saved_tie_ct);
    }
  }
  cur_r1[0] = mode0;
  cur_r1[1] = mode1;
  cur_r1[2] = mode2;
  return kFisher2xkNonempty;
}

// On entry, cur_r1[0..level] must be the mode of the 2x(level+1) subproblem
// with row-1 total row1_tot and the current cur_r1[(level+1)..] values.  It is
// restored on exit.
static uint32_t Fisher2xkLevel(Fisher2xkCtx* ctxp, uint32_t level, int32_t row1_tot) {
  if (level == 2) {
    return Fisher2xkLevel2(ctxp, row1_tot);
  }
  int32_t* cur_r1 = ctxp->cur_r1;
  const int32_t* col_tots = ctxp->col_tots;
  int32_t* saved_mode = &(ctxp->saved_modes[level * ctxp->col_ct]);
  const int32_t start_x = cur_r1[level];
  const int32_t x_min = MAXV(0, row1_tot - ctxp->cum_col_tots[level]);
  const int32_t x_max = MINV(col_tots[level], row1_tot);
  memcpy(saved_mode, cur_r1, level * sizeof(int32_t));
  for (uint32_t x_decreasing = 0; x_decreasing != 2; ++x_decreasing) {
    int32_t x = start_x;
    if (x_decreasing) {
      if (x == x_min) {
        break;
      }
      memcpy(cur_r1, saved_mode, level * sizeof(int32_t));
      --x;
      Fisher2xkModeIncr(col_tots, level, cur_r1);
    }
    while (1) {
      cur_r1[level] = x;
      uint32_t inner_status;
      const double weight = Fisher2xkSymFactor(ctxp, level);
      if (weight == 1) {
        inner_status = Fisher2xkLevel(ctxp, level - 1, row1_tot - x);
      } else {
        if (weight == 0) {
          // If this rank is nonempty, its contribution is accounted for by
          // the permutation in nonincreasing order.
          inner_status = Fisher2xkModeStatus(ctxp);
        } else {
          const double saved_outer_sum = ctxp->outer_sum;
          const double saved_tie_ct = ctxp->tie_ct;
          ctxp->outer_sum = 0.0;
          ctxp->tie_ct = 0.0;
          inner_status = Fisher2xkLevel(ctxp, level - 1, row1_tot - x);
          ctxp->outer_sum = saved_outer_sum + weight * ctxp->outer_sum;
          ctxp->tie_ct = saved_tie_ct + weight * ctxp->tie_ct;
        }
      }
      if (inner_status != kFisher2xkNonempty) {
        if ((!x_decreasing) && (x == start_x)) {
          // Modal rank is empty, so all ranks are empty.
          return inner_status;
        }
        if (inner_status == kFisher2xkEmptyTie) {
          ctxp->tie_ct += Fisher2xkModeCt(col_tots, level, cur_r1);
        }
        Fisher2xkClosedFormSum(ctxp, level, row1_tot, x, x_decreasing? -1 : 1);
        break;
      }
      if (x_decreasing) {
        if (x == x_min) {
          break;
        }
        --x;
        Fisher2xkModeIncr(col_tots, level, cur_r1);
      } else {
        if (x == x_max) {
          break;
        }
        ++x;
        Fisher2xkModeDecr(col_tots, level, cur_r1);
      }
    }
  }
  memcpy(cur_r1, saved_mode, level * sizeof(int32_t));
  cur_r1[level] = start_x;
  return kFisher2xkNonempty;
}

double Fisher2xkLnP(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff) {
  uint32_t nz_col_ct = 0;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    nz_col_ct += (obs[col_idx] + obs[col_ct + col_idx] != 0);
  }
  if (nz_col_ct <= 3) {
    int32_t small_obs[6] = {0, 0, 0, 0, 0, 0};
    uint32_t write_idx = 0;
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      if (obs[col_idx] + obs[col_ct + col_idx] != 0) {
        small_obs[write_idx] = obs[col_idx];
        small_obs[write_idx + 3] = obs[col_ct + col_idx];
        ++write_idx;
      }
    }
    if (nz_col_ct == 3) {
      return Fisher23LnP(small_obs[0], small_obs[1], small_obs[2], small_obs[3], small_obs[4], small_obs[5], midp);
    }
    return Fisher22TwoSidedP(small_obs[0], small_obs[1], small_obs[3], small_obs[4], midp, 1);
  }
  const uint32_t raw_col_ct = col_ct;
  col_ct = nz_col_ct;
  unsigned char* alloc_base;
  if (unlikely(pgl_malloc(col_ct * (sizeof(dd_real) + 4 * sizeof(int64_t) + 2 * sizeof(double) + (4 + col_ct) * sizeof(int32_t) + sizeof(int32_t)) + 4 * sizeof(double), &alloc_base))) {
    return S_CAST(double, NAN);
  }
  dd_real* inner_const_ddrs = R_CAST(dd_real*, alloc_base);
  uint64_t* fac_args_buf = R_CAST(uint64_t*, &(inner_const_ddrs[col_ct]));
  double* lfact_args_buf = R_CAST(double*, &(fac_args_buf[4 * col_ct]));
  int32_t* col_tots = R_CAST(int32_t*, &(lfact_args_buf[2 * col_ct + 4]));
  int32_t* obs_r1 = &(col_tots[col_ct]);
  int32_t* cum_col_tots = &(obs_r1[col_ct]);
  int32_t* cur_r1 = &(cum_col_tots[col_ct]);
  int32_t* saved_modes = &(cur_r1[col_ct]);
  uint32_t* sym_ns = R_CAST(uint32_t*, &(saved_modes[col_ct * col_ct]));

  // Remove empty columns, and sort the rest in decreasing-total order.
  // (Insertion sort is fine, since runtime is exponential in col_ct.)
  int32_t row1_tot = 0;
  int32_t tot = 0;
  {
    uint32_t write_idx = 0;
    for (uint32_t col_idx = 0; write_idx != col_ct; ++col_idx) {
      const int32_t cur_obs_r1 = obs[col_idx];
      const int32_t cur_col_tot = cur_obs_r1 + obs[raw_col_ct + col_idx];
      if (!cur_col_tot) {
        continue;
      }
      uint32_t insert_idx = write_idx;
      for (; insert_idx; --insert_idx) {
        if (col_tots[insert_idx - 1] >= cur_col_tot) {
          break;
        }
        col_tots[insert_idx] = col_tots[insert_idx - 1];
        obs_r1[insert_idx] = obs_r1[insert_idx - 1];
      }
      col_tots[insert_idx] = cur_col_tot;
      obs_r1[insert_idx] = cur_obs_r1;
      ++write_idx;
      row1_tot += cur_obs_r1;
      tot += cur_col_tot;
    }
  }
  // Find runs of equal-total columns, ignoring columns 0 and 1.
  for (uint32_t run_start = 0; run_start != col_ct; ) {
    uint32_t run_end = run_start + 1;
    while ((run_end != col_ct) && (col_tots[run_end] == col_tots[run_start])) {
      ++run_end;
    }
    for (uint32_t col_idx = run_start; col_idx != run_end; ++col_idx) {
      sym_ns[col_idx] = (col_idx >= 2)? (run_end - col_idx) : 0;
    }
    sym_ns[run_end - 1] = 0;
    run_start = run_end;
  }
  Fisher2xkCtx ctx;
  ctx.col_ct = col_ct;
  ctx.col_tots = col_tots;
  ctx.obs_r1 = obs_r1;
  ctx.cum_col_tots = cum_col_tots;
  ctx.inner_const_ddrs = inner_const_ddrs;
  ctx.cur_r1 = cur_r1;
  ctx.saved_modes = saved_modes;
  ctx.sym_ns = sym_ns;
  ctx.fac_args_buf = fac_args_buf;
  ctx.lfact_args_buf = lfact_args_buf;
  ctx.full_obs = nullptr;
  ctx.fixed_cells = nullptr;
  ctx.full_cell_ct = 0;
  ctx.outer_sum = 0.0;
  ctx.tie_ct = 0.0;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    lfact_args_buf[2 * col_idx] = obs_r1[col_idx];
    lfact_args_buf[2 * col_idx + 1] = col_tots[col_idx] - obs_r1[col_idx];
  }
  ctx.starting_lnprobv_ddr = ddr_negate(ddr_sort_and_add_lfacts(2 * col_ct, lfact_args_buf));
  ctx.starting_lnprobv_tdr.x[0] = ctx.starting_lnprobv_ddr.x[0];
  ctx.starting_lnprobv_tdr.x[1] = ctx.starting_lnprobv_ddr.x[1];
  ctx.starting_lnprobv_tdr.x[2] = DBL_MAX;
  // lnprobf, so we can add it to starting_lnprobv at the end and convert
  // outer_sum into the final log-[mid]p-value.
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    lfact_args_buf[col_idx] = col_tots[col_idx];
  }
  lfact_args_buf[col_ct] = row1_tot;
  lfact_args_buf[col_ct + 1] = tot - row1_tot;
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_sort_and_add_lfacts(col_ct + 2, lfact_args_buf),
            ddr_lfact(tot));
  const double starting_lnprob =
    ddr_add(lnprobf_ddr, ctx.starting_lnprobv_ddr).x[0];
  {
    // Every table counted in the p-value is no more probable than the
    // observed table, so the p-value is bounded above by (number of tables
    // with these margins) * (observed table probability).  The number of
    // tables is at most the product of the row-1 value range sizes of all
    // columns but one.
    const int32_t row2_tot = tot - row1_tot;
    double ln_table_ct = 0.0;
    double max_ln_range = 0.0;
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      const int32_t cur_col_tot = col_tots[col_idx];
      const int32_t range = MINV(cur_col_tot, row1_tot) - MAXV(0, cur_col_tot - row2_tot) + 1;
      const double ln_range = log(S_CAST(double, range));
      ln_table_ct += ln_range;
      if (ln_range > max_ln_range) {
        max_ln_range = ln_range;
      }
    }
    double ln_p_bound = starting_lnprob + ln_table_ct - max_ln_range;
    // Leave a margin for rounding error.
    ln_p_bound += fabs(ln_p_bound) * k2m35;
    if (ln_p_bound < ln_p_cutoff) {
      free(alloc_base);
      return MINV(ln_p_bound, 0.0);
    }
  }
  {
    int32_t cum_col_tot = 0;
    dd_real lfact_sum_ddr = ddr_maked(0.0);
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      cum_col_tots[col_idx] = cum_col_tot;
      inner_const_ddrs[col_idx] = ddr_sub(ddr_lfact(cum_col_tot), lfact_sum_ddr);
      cum_col_tot += col_tots[col_idx];
      lfact_sum_ddr = ddr_add(lfact_sum_ddr, ddr_lfact(col_tots[col_idx]));
    }
  }
  Fisher2xkInitMode(col_tots, col_ct, tot, row1_tot, cur_r1);
  const uint32_t status = Fisher2xkLevel(&ctx, col_ct - 1, row1_tot);
  double result;
  if (status != kFisher2xkNonempty) {
    // Observed table is modal (kFisher2xkEmpty should be impossible), so the
    // non-mid p-value is 1.
    result = 0.0;
    if (midp) {
      result = log1p(-0.5 * Fisher2xkModeCt(col_tots, col_ct, cur_r1) * exp(starting_lnprob));
    }
  } else {
    if (midp) {
      ctx.outer_sum -= ctx.tie_ct * 0.5;
    }
    result = log(ctx.outer_sum) + starting_lnprob;
    if (result > -k2m35) {
      // true p-value should always be 1 here
      result = 0.0;
    }
  }
  free(alloc_base);
  return result;
}

// 3xk generalization of Fisher2xkLnP().
//
// Each column j >= 1 contributes two levels of recursion: first over its
// row-1 entry (the "A" level), then over its row-2 entry (the "B" level).  Row
// 3 and column 0 are determined by the margins.  The overall structure at
// each level is the same as in Fisher2xkLevel(): we start from the modal rank
// of the current subproblem, iterate outward in both directions until we hit
// a rank where no table is more probable than the observed table, and then
// sum the remaining ranks in that direction in closed form.  This is
// justified by the following:
// - The log-likelihood of a column, as a function of its row-1 and row-2
//   entries a and b, is -lfact(a) - lfact(b) - lfact(c-a-b); this is
//   M-natural-concave.  Max-convolutions of M-natural-concave functions are
//   M-natural-concave, and in two dimensions, sums of M-natural-concave
//   functions are M-natural-concave.  So the maximum likelihood over a rank,
//   as a function of the (a, b) entries of the current column, is
//   M-natural-concave; in particular, it's concave along each b-line, and its
//   maximum over each b-line is concave in a.
// - The mode of a subproblem is the solution to a separable convex
//   transportation problem.  When a level variable changes by 1, the new mode
//   can be obtained from the old one by applying the best unit cycle through
//   that cell (successive shortest paths), and with 3 rows, these cycles have
//   length 4 or 6.
// - Rank sums have closed forms via the multinomial Vandermonde identity.
// The innermost two levels (columns 0 and 1, with everything else fixed) are
// handled by Fisher2xkLevel2(), since the transpose of that subproblem is a
// 2x3 table.
// This code is written for general 3xk tables, but beyond 3x3 it's usually
// slower than the network algorithm (which merges subproblems with the same
// remaining row totals), so only the 3x3 entry point is exposed for now.
//
// Cells are identified by "keys": (r, j) for j >= 1 has key 3j - r, and
// column-0 cells have key -1.  The level variable with key L is (0, L/3) when
// L is divisible by 3, and (1, (L+1)/3) when L == 2 (mod 3); when we're at
// level L, cells with key < L are free and the rest are fixed.

typedef struct Fisher3xkCtxStruct {
  uint32_t col_ct;
  const int32_t* col_tots;
  const int32_t* row_tots;
  // obs and cur are in row-major order.
  const int32_t* obs;
  int32_t* cur;
  // Level L saves a copy of cur at saved_tables[3 * col_ct * L].
  int32_t* saved_tables;
  unsigned char* flex_buf;
  double* lfact_args_buf;
  dd_real starting_lnprobv_ddr;
  td_real starting_lnprobv_tdr;
  double outer_sum;
  double tie_ct;
} Fisher3xkCtx;

static inline int32_t Fisher3xkCellKey(uint32_t row_idx, uint32_t col_idx) {
  if (!col_idx) {
    return -1;
  }
  return S_CAST(int32_t, 3 * col_idx) - S_CAST(int32_t, row_idx);
}

// Sets limbs[0..6] to the product of factors[0..(factor_ct-1)], where
// factor_ct <= 6 and each factor is less than 2^32.
static void U32ProductLimbs(const uint32_t* factors, uint32_t factor_ct, uint32_t* limbs) {
  limbs[0] = 1;
  for (uint32_t limb_idx = 1; limb_idx != 7; ++limb_idx) {
    limbs[limb_idx] = 0;
  }
  for (uint32_t factor_idx = 0; factor_idx != factor_ct; ++factor_idx) {
    const uint64_t factor = factors[factor_idx];
    uint64_t carry = 0;
    for (uint32_t limb_idx = 0; limb_idx != 7; ++limb_idx) {
      const uint64_t prod = limbs[limb_idx] * factor + carry;
      limbs[limb_idx] = S_CAST(uint32_t, prod);
      carry = prod >> 32;
    }
  }
}

// Exact comparison of two products of at most 6 factors each.
static intptr_t CompareU32Products(const uint32_t* a, uint32_t a_ct, const uint32_t* b, uint32_t b_ct) {
  uint32_t a_limbs[7];
  uint32_t b_limbs[7];
  U32ProductLimbs(a, a_ct, a_limbs);
  U32ProductLimbs(b, b_ct, b_limbs);
  for (uint32_t limb_idx = 7; limb_idx; ) {
    --limb_idx;
    if (a_limbs[limb_idx] != b_limbs[limb_idx]) {
      return (a_limbs[limb_idx] > b_limbs[limb_idx])? 1 : -1;
    }
  }
  return 0;
}

// Computes the likelihood ratio of applying the given unit changes to cur[],
// as prod(numer) / prod(denom).  Returns 0 if some decremented cell is zero.
static uint32_t Fisher3xkCycleFactors(const int32_t* cur, const uint32_t* cells, const int32_t* deltas, uint32_t cell_ct, uint32_t* numer, uint32_t* numer_ct_ptr, uint32_t* denom, uint32_t* denom_ct_ptr, double* ratio_ptr) {
  uint32_t numer_ct = 0;
  uint32_t denom_ct = 0;
  double ratio = 1.0;
  for (uint32_t uii = 0; uii != cell_ct; ++uii) {
    const int32_t val = cur[cells[uii]];
    if (deltas[uii] > 0) {
      denom[denom_ct++] = val + 1;
      ratio /= val + 1;
    } else {
      if (!val) {
        return 0;
      }
      numer[numer_ct++] = val;
      ratio *= val;
    }
  }
  *numer_ct_ptr = numer_ct;
  *denom_ct_ptr = denom_ct;
  *ratio_ptr = ratio;
  return 1;
}

// Returns positive if a_numer/a_denom > b_numer/b_denom, etc.
static intptr_t Fisher3xkCompareRatios(const uint32_t* a_numer, uint32_t a_numer_ct, const uint32_t* a_denom, uint32_t a_denom_ct, const uint32_t* b_numer, uint32_t b_numer_ct, const uint32_t* b_denom, uint32_t b_denom_ct) {
  uint32_t lhs[6];
  uint32_t rhs[6];
  memcpy(lhs, a_numer, a_numer_ct * sizeof(int32_t));
  memcpy(&(lhs[a_numer_ct]), b_denom, b_denom_ct * sizeof(int32_t));
  memcpy(rhs, b_numer, b_numer_ct * sizeof(int32_t));
  memcpy(&(rhs[b_numer_ct]), a_denom, a_denom_ct * sizeof(int32_t));
  return CompareU32Products(lhs, a_numer_ct + b_denom_ct, rhs, b_numer_ct + a_denom_ct);
}

// Changes the level variable at (row_idx, col_idx) by delta (+1 or -1), and
// updates the free cells to the mode of the new subproblem.  Returns 0
// (without changing anything) if the new value is infeasible.
static uint32_t Fisher3xkShift(Fisher3xkCtx* ctxp, uint32_t row_idx, uint32_t col_idx, int32_t delta) {
  const uint32_t col_ct = ctxp->col_ct;
  int32_t* cur = ctxp->cur;
  const uint32_t var_cell = row_idx * col_ct + col_idx;
  if ((delta < 0) && (!cur[var_cell])) {
    return 0;
  }
  uint32_t best_cells[5];
  uint32_t best_cell_ct = 0;
  uint32_t best_numer[3];
  uint32_t best_numer_ct = 0;
  uint32_t best_denom[3];
  uint32_t best_denom_ct = 0;
  double best_ratio = 0.0;
  const int32_t deltas[5] = {-delta, delta, -delta, delta, -delta};
  for (uint32_t row2_idx = row_idx + 1; row2_idx != 3; ++row2_idx) {
    const uint32_t row3_idx = 3 - row_idx - row2_idx;
    for (uint32_t col2_idx = 0; col2_idx != col_idx; ++col2_idx) {
      uint32_t cells[5];
      cells[0] = row2_idx * col_ct + col_idx;
      cells[1] = row2_idx * col_ct + col2_idx;
      for (uint32_t col3_idx = 0; col3_idx != col_idx; ++col3_idx) {
        uint32_t cell_ct;
        if (col3_idx == col2_idx) {
          // 4-cycle
          cells[2] = row_idx * col_ct + col2_idx;
          cell_ct = 3;
        } else {
          // 6-cycle
          cells[2] = row3_idx * col_ct + col2_idx;
          cells[3] = row3_idx * col_ct + col3_idx;
          cells[4] = row_idx * col_ct + col3_idx;
          cell_ct = 5;
        }
        uint32_t numer[3];
        uint32_t numer_ct;
        uint32_t denom[3];
        uint32_t denom_ct;
        double ratio;
        if (!Fisher3xkCycleFactors(cur, cells, deltas, cell_ct, numer, &numer_ct, denom, &denom_ct, &ratio)) {
          continue;
        }
        if (best_cell_ct) {
          if (ratio < best_ratio * (1 + kSmallEpsilon)) {
            if ((ratio <= best_ratio * (1 - kSmallEpsilon)) || (Fisher3xkCompareRatios(numer, numer_ct, denom, denom_ct, best_numer, best_numer_ct, best_denom, best_denom_ct) <= 0)) {
              continue;
            }
          }
        }
        memcpy(best_cells, cells, cell_ct * sizeof(int32_t));
        best_cell_ct = cell_ct;
        memcpy(best_numer, numer, numer_ct * sizeof(int32_t));
        best_numer_ct = numer_ct;
        memcpy(best_denom, denom, denom_ct * sizeof(int32_t));
        best_denom_ct = denom_ct;
        best_ratio = ratio;
      }
    }
  }
  if (!best_cell_ct) {
    return 0;
  }
  cur[var_cell] += delta;
  for (uint32_t uii = 0; uii != best_cell_ct; ++uii) {
    cur[best_cells[uii]] += deltas[uii];
  }
  return 1;
}

// With improve set, this applies the first likelihood-increasing unit cycle
// among free cells (key < key_thresh) that it finds, and returns 1 (or
// returns 0 if there is none).
// Otherwise, for every likelihood-preserving unit cycle among free cells,
// flex_buf[] bit 0 is set for the incremented cells and bit 1 is set for the
// decremented cells.
static uint32_t Fisher3xkScanCycles(Fisher3xkCtx* ctxp, int32_t key_thresh, uint32_t improve) {
  const uint32_t col_ct = ctxp->col_ct;
  int32_t* cur = ctxp->cur;
  unsigned char* flex = ctxp->flex_buf;
  const int32_t deltas[6] = {1, -1, 1, -1, 1, -1};
  // 4-cycles: (p,u)+ (q,u)- (q,v)+ (p,v)-
  // 6-cycles: (p,u)+ (q,u)- (q,v)+ (s,v)- (s,w)+ (p,w)-
  for (uint32_t row_p = 0; row_p != 3; ++row_p) {
    for (uint32_t row_q = 0; row_q != 3; ++row_q) {
      if (row_q == row_p) {
        continue;
      }
      const uint32_t row_s = 3 - row_p - row_q;
      for (uint32_t col_u = 0; col_u != col_ct; ++col_u) {
        if ((Fisher3xkCellKey(row_p, col_u) >= key_thresh) || (Fisher3xkCellKey(row_q, col_u) >= key_thresh)) {
          continue;
        }
        for (uint32_t col_v = 0; col_v != col_ct; ++col_v) {
          if ((col_v == col_u) || (Fisher3xkCellKey(row_q, col_v) >= key_thresh)) {
            continue;
          }
          for (uint32_t col_w = 0; col_w <= col_ct; ++col_w) {
            uint32_t cells[6];
            uint32_t cell_ct;
            cells[0] = row_p * col_ct + col_u;
            cells[1] = row_q * col_ct + col_u;
            cells[2] = row_q * col_ct + col_v;
            if (col_w == col_ct) {
              if (Fisher3xkCellKey(row_p, col_v) >= key_thresh) {
                continue;
              }
              cells[3] = row_p * col_ct + col_v;
              cell_ct = 4;
            } else {
              if ((col_w == col_u) || (col_w == col_v) || (Fisher3xkCellKey(row_s, col_v) >= key_thresh) || (Fisher3xkCellKey(row_s, col_w) >= key_thresh) || (Fisher3xkCellKey(row_p, col_w) >= key_thresh)) {
                continue;
              }
              cells[3] = row_s * col_ct + col_v;
              cells[4] = row_s * col_ct + col_w;
              cells[5] = row_p * col_ct + col_w;
              cell_ct = 6;
            }
            uint32_t numer[3];
            uint32_t numer_ct;
            uint32_t denom[3];
            uint32_t denom_ct;
            double ratio;
            if (!Fisher3xkCycleFactors(cur, cells, deltas, cell_ct, numer, &numer_ct, denom, &denom_ct, &ratio)) {
              continue;
            }
            if (improve) {
              if (ratio <= 1 + kSmallEpsilon) {
                if ((ratio <= 1 - kSmallEpsilon) || (CompareU32Products(numer, numer_ct, denom, denom_ct) <= 0)) {
                  continue;
                }
              }
              for (uint32_t uii = 0; uii != cell_ct; ++uii) {
                cur[cells[uii]] += deltas[uii];
              }
              return 1;
            }
            if ((fabs(ratio - 1) < kSmallEpsilon) && (!CompareU32Products(numer, numer_ct, denom, denom_ct))) {
              for (uint32_t uii = 0; uii != cell_ct; ++uii) {
                flex[cells[uii]] |= (deltas[uii] > 0)? 1 : 2;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

// Number of modal tables of the subproblem where cells with key < key_thresh
// are free, given that cur[] is one of them.
// Each free cell's value is the same across all modal tables, except for
// "flexible" cells which can take one of two adjacent values; every
// assignment of those values with the right margins is modal (by
// complementary slackness).  Cells are flexible iff they're part of a
// likelihood-preserving unit cycle.
static double Fisher3xkModeCt(Fisher3xkCtx* ctxp, int32_t key_thresh) {
  const uint32_t col_ct = ctxp->col_ct;
  unsigned char* flex = ctxp->flex_buf;
  memset(flex, 0, 3 * col_ct);
  Fisher3xkScanCycles(ctxp, key_thresh, 0);
  // Flexible cells that are currently at the higher value contribute to the
  // row/column "need" counts.
  uint32_t row_needs[3] = {0, 0, 0};
  for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      row_needs[row_idx] += (flex[row_idx * col_ct + col_idx] >> 1) & 1;
    }
  }
  if (!(row_needs[0] + row_needs[1] + row_needs[2])) {
    return 1;
  }
  const uint32_t dp_width = row_needs[1] + 1;
  const uint32_t dp_size = (row_needs[0] + 1) * dp_width;
  double* dp;
  if (unlikely(pgl_malloc(2 * dp_size * sizeof(double), &dp))) {
    return S_CAST(double, NAN);
  }
  double* dp_next = &(dp[dp_size]);
  for (uint32_t uii = 0; uii != dp_size; ++uii) {
    dp[uii] = 0.0;
  }
  dp[0] = 1.0;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    uint32_t col_flex_mask = 0;
    uint32_t col_need = 0;
    for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
      const uint32_t cur_flex = flex[row_idx * col_ct + col_idx];
      if (cur_flex) {
        col_flex_mask |= 1U << row_idx;
        col_need += cur_flex >> 1;
      }
    }
    if (!col_flex_mask) {
      continue;
    }
    for (uint32_t uii = 0; uii != dp_size; ++uii) {
      dp_next[uii] = 0.0;
    }
    for (uint32_t subset = col_flex_mask; ; subset = (subset - 1) & col_flex_mask) {
      if (S_CAST(uint32_t, __builtin_popcount(subset)) == col_need) {
        const uint32_t incr0 = subset & 1;
        const uint32_t incr1 = (subset >> 1) & 1;
        for (uint32_t n0 = incr0; n0 <= row_needs[0]; ++n0) {
          for (uint32_t n1 = incr1; n1 <= row_needs[1]; ++n1) {
            dp_next[n0 * dp_width + n1] += dp[(n0 - incr0) * dp_width + n1 - incr1];
          }
        }
      }
      if (!subset) {
        break;
      }
    }
    double* tmp = dp;
    dp = dp_next;
    dp_next = tmp;
  }
  const double result = dp[row_needs[0] * dp_width + row_needs[1]];
  free((dp < dp_next)? dp : dp_next);
  return result;
}

// Adds the closed-form sum of ranks v, v+dir, v+2*dir, ... at the given
// level, where v is the current value of the level variable, and all tables
// in these ranks are known to be no more probable than the observed table.
static void Fisher3xkClosedFormSum(Fisher3xkCtx* ctxp, int32_t lvl, int32_t dir) {
  const uint32_t col_ct = ctxp->col_ct;
  const int32_t* col_tots = ctxp->col_tots;
  const int32_t* cur = ctxp->cur;
  const uint32_t col_idx = (lvl + 2) / 3;
  const uint32_t is_a = !(lvl % 3);
  // Row totals remaining after columns > col_idx.
  int64_t rem[3];
  for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
    rem[row_idx] = ctxp->row_tots[row_idx];
    for (uint32_t col2_idx = col_idx + 1; col2_idx != col_ct; ++col2_idx) {
      rem[row_idx] -= cur[row_idx * col_ct + col2_idx];
    }
  }
  int64_t inner_col_tot = 0;
  double* neg_args = ctxp->lfact_args_buf;
  uint32_t neg_ct = 0;
  for (uint32_t col2_idx = 0; col2_idx != col_idx; ++col2_idx) {
    inner_col_tot += col_tots[col2_idx];
    neg_args[neg_ct++] = col_tots[col2_idx];
  }
  for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
    for (uint32_t col2_idx = col_idx + 1; col2_idx != col_ct; ++col2_idx) {
      neg_args[neg_ct++] = cur[row_idx * col_ct + col2_idx];
    }
  }
  const double col_tot = col_tots[col_idx];
  const double r1 = rem[0];
  const double r2 = rem[1];
  const double r3 = rem[2];
  const double inner_tot = inner_col_tot;
  double pos_args[2];
  uint32_t pos_ct;
  double x;
  double a = 0.0;
  if (is_a) {
    // Sum of 1/(cell factorials) over the rank is
    //   (C! T!) / (v! (r1-v)! (c-v)! (C-r1+v)! r2! r3! prod_{i<j} c_i!)
    // where C = sum_{i<j} c_i and T = r2 + r3.
    x = cur[col_idx];
    pos_args[0] = inner_tot;
    pos_args[1] = r2 + r3;
    pos_ct = 2;
    neg_args[neg_ct++] = x;
    neg_args[neg_ct++] = r1 - x;
    neg_args[neg_ct++] = col_tot - x;
    neg_args[neg_ct++] = inner_tot - r1 + x;
    neg_args[neg_ct++] = r2;
    neg_args[neg_ct++] = r3;
  } else {
    // Sum is
    //   C! / (a! v! (c-a-v)! (r1-a)! (r2-v)! (r3-c+a+v)! prod_{i<j} c_i!)
    a = cur[col_idx];
    x = cur[col_ct + col_idx];
    pos_args[0] = inner_tot;
    pos_ct = 1;
    neg_args[neg_ct++] = a;
    neg_args[neg_ct++] = x;
    neg_args[neg_ct++] = col_tot - a - x;
    neg_args[neg_ct++] = r1 - a;
    neg_args[neg_ct++] = r2 - x;
    neg_args[neg_ct++] = r3 - col_tot + a + x;
  }
  const dd_real rank_lnlik_ddr =
    ddr_sub(ddr_sub(ddr_sort_and_add_lfacts(pos_ct, pos_args),
                    ddr_sort_and_add_lfacts(neg_ct, neg_args)),
            ctxp->starting_lnprobv_ddr);
  double rank_lik = exp(rank_lnlik_ddr.x[0]);
  double outer_sum = ctxp->outer_sum;
  while (1) {
    const double preadd = outer_sum;
    outer_sum += rank_lik;
    double rank_lik_mult;
    if (is_a) {
      if (dir > 0) {
        rank_lik_mult = (r1 - x) * (col_tot - x) / ((x + 1) * (inner_tot - r1 + x + 1));
      } else {
        rank_lik_mult = x * (inner_tot - r1 + x) / ((r1 - x + 1) * (col_tot - x + 1));
      }
    } else {
      if (dir > 0) {
        rank_lik_mult = (col_tot - a - x) * (r2 - x) / ((x + 1) * (r3 - col_tot + a + x + 1));
      } else {
        rank_lik_mult = x * (r3 - col_tot + a + x) / ((col_tot - a - x + 1) * (r2 - x + 1));
      }
    }
    x += dir;
    // The multiplier is zero when we run out of tables.  Otherwise, rank sums
    // are log-concave, so we're done once they're decreasing and negligible.
    if ((rank_lik_mult <= 0.0) || ((outer_sum == preadd) && (rank_lik_mult <= 1))) {
      break;
    }
    rank_lik *= rank_lik_mult;
  }
  ctxp->outer_sum = outer_sum;
}

// Innermost levels: columns 0 and 1 vary, and everything else is fixed.
// The transpose of this 3x2 subproblem is a 2x3 table, so we hand it to
// Fisher2xkLevel2(), which tracks the tail boundaries incrementally.  Return
// values and cur[] requirements are as for Fisher3xkLevel().
static uint32_t Fisher33Base(Fisher3xkCtx* ctxp) {
  const uint32_t col_ct = ctxp->col_ct;
  const int32_t* cur = ctxp->cur;
  // View rows are columns 0 and 1 of the original table; view columns are
  // the original rows, in decreasing remaining-total order.
  int32_t row_rems[3];
  for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
    row_rems[row_idx] = cur[row_idx * col_ct] + cur[row_idx * col_ct + 1];
  }
  uint32_t perm[3] = {0, 1, 2};
  for (uint32_t uii = 1; uii != 3; ++uii) {
    for (uint32_t ujj = uii; ujj && (row_rems[perm[ujj - 1]] < row_rems[perm[ujj]]); --ujj) {
      const uint32_t tmp = perm[ujj - 1];
      perm[ujj - 1] = perm[ujj];
      perm[ujj] = tmp;
    }
  }
  int32_t view_col_tots[3];
  int32_t view_cur_r1[3];
  int32_t view_cum_col_tots[3];
  dd_real view_inner_const_ddrs[3];
  uint32_t view_sym_ns[3] = {0, 0, 0};
  // Buffer sizes assume col_ct == 3.
  int32_t fixed_cells[3];
  uint64_t fac_args_buf[18];
  double lfact_args_buf[16];
  for (uint32_t view_col_idx = 0; view_col_idx != 3; ++view_col_idx) {
    const uint32_t row_idx = perm[view_col_idx];
    view_col_tots[view_col_idx] = row_rems[row_idx];
    view_cur_r1[view_col_idx] = cur[row_idx * col_ct];
  }
  // Only the level-2 entries of cum_col_tots and inner_const_ddrs are used.
  view_cum_col_tots[2] = view_col_tots[0] + view_col_tots[1];
  view_inner_const_ddrs[2] = ddr_sub(ddr_lfact(view_cum_col_tots[2]), ddr_add_lfacts(view_col_tots[0], view_col_tots[1]));
  uint32_t fixed_ct = 0;
  for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
    for (uint32_t col_idx = 2; col_idx != col_ct; ++col_idx) {
      fixed_cells[fixed_ct] = cur[row_idx * col_ct + col_idx];
      lfact_args_buf[fixed_ct] = fixed_cells[fixed_ct];
      ++fixed_ct;
    }
  }
  Fisher2xkCtx view;
  view.col_ct = 3;
  view.col_tots = view_col_tots;
  // obs_r1 is only used by Fisher2xkCompare(), where it's overwritten by
  // full_obs.
  view.obs_r1 = view_cur_r1;
  view.cum_col_tots = view_cum_col_tots;
  view.inner_const_ddrs = view_inner_const_ddrs;
  view.cur_r1 = view_cur_r1;
  view.saved_modes = nullptr;
  view.sym_ns = view_sym_ns;
  view.fac_args_buf = fac_args_buf;
  view.lfact_args_buf = lfact_args_buf;
  view.starting_lnprobv_ddr = ddr_add(ctxp->starting_lnprobv_ddr, ddr_sort_and_add_lfacts(fixed_ct, lfact_args_buf));
  view.starting_lnprobv_tdr = ctxp->starting_lnprobv_tdr;
  view.full_obs = ctxp->obs;
  view.fixed_cells = fixed_cells;
  view.full_cell_ct = 3 * col_ct;
  view.outer_sum = ctxp->outer_sum;
  view.tie_ct = ctxp->tie_ct;
  const uint32_t status = Fisher2xkLevel2(&view, ctxp->col_tots[0]);
  ctxp->outer_sum = view.outer_sum;
  ctxp->tie_ct = view.tie_ct;
  return status;
}
// On entry, cur[] must be the mode of the subproblem where cells with key <=
// lvl are free.  It is restored on exit.  Return values are as for
// Fisher2xkLevel().
static uint32_t Fisher3xkLevel(Fisher3xkCtx* ctxp, int32_t lvl) {
  if (lvl == 3) {
    return Fisher33Base(ctxp);
  }
  const uint32_t col_ct = ctxp->col_ct;
  const uint32_t table_size = 3 * col_ct;
  int32_t* cur = ctxp->cur;
  const uint32_t col_idx = (lvl + 2) / 3;
  const uint32_t row_idx = (lvl % 3)? 1 : 0;
  const int32_t next_lvl = row_idx? (lvl - 2) : (lvl - 1);
  int32_t* saved = &(ctxp->saved_tables[S_CAST(uintptr_t, lvl) * table_size]);
  memcpy(saved, cur, table_size * sizeof(int32_t));
  for (uint32_t decreasing = 0; decreasing != 2; ++decreasing) {
    const int32_t dir = decreasing? -1 : 1;
    if (decreasing) {
      memcpy(cur, saved, table_size * sizeof(int32_t));
      if (!Fisher3xkShift(ctxp, row_idx, col_idx, -1)) {
        break;
      }
    }
    uint32_t is_start = !decreasing;
    while (1) {
      const uint32_t inner_status = Fisher3xkLevel(ctxp, next_lvl);
      if (inner_status != kFisher2xkNonempty) {
        if (is_start) {
          // Modal rank is empty, so all ranks are empty.
          return inner_status;
        }
        if (inner_status == kFisher2xkEmptyTie) {
          ctxp->tie_ct += Fisher3xkModeCt(ctxp, lvl);
        }
        Fisher3xkClosedFormSum(ctxp, lvl, dir);
        break;
      }
      is_start = 0;
      if (!Fisher3xkShift(ctxp, row_idx, col_idx, dir)) {
        break;
      }
    }
  }
  memcpy(cur, saved, table_size * sizeof(int32_t));
  return kFisher2xkNonempty;
}

static double Fisher3xkLnP(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff) {
  int64_t row_tots_i64[3] = {0, 0, 0};
  uint32_t nz_col_ct = 0;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    const int32_t v0 = obs[col_idx];
    const int32_t v1 = obs[col_ct + col_idx];
    const int32_t v2 = obs[2 * col_ct + col_idx];
    row_tots_i64[0] += v0;
    row_tots_i64[1] += v1;
    row_tots_i64[2] += v2;
    nz_col_ct += ((v0 + v1 + v2) != 0);
  }
  for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
    if (!row_tots_i64[row_idx]) {
      // Reduce to 2xk.
      int32_t* two_row_obs;
      if (unlikely(pgl_malloc(2 * col_ct * sizeof(int32_t), &two_row_obs))) {
        return S_CAST(double, NAN);
      }
      const uint32_t row_a = (row_idx == 0)? 1 : 0;
      const uint32_t row_b = (row_idx == 2)? 1 : 2;
      memcpy(two_row_obs, &(obs[row_a * col_ct]), col_ct * sizeof(int32_t));
      memcpy(&(two_row_obs[col_ct]), &(obs[row_b * col_ct]), col_ct * sizeof(int32_t));
      const double result = Fisher2xkLnP(two_row_obs, col_ct, midp, ln_p_cutoff);
      free(two_row_obs);
      return result;
    }
  }
  if (nz_col_ct <= 2) {
    // Transpose to 2x3.
    int32_t transposed[6] = {0, 0, 0, 0, 0, 0};
    uint32_t write_idx = 0;
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      if (obs[col_idx] + obs[col_ct + col_idx] + obs[2 * col_ct + col_idx]) {
        for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
          transposed[write_idx * 3 + row_idx] = obs[row_idx * col_ct + col_idx];
        }
        ++write_idx;
      }
    }
    return Fisher2xkLnP(transposed, 3, midp, ln_p_cutoff);
  }
  const uint32_t raw_col_ct = col_ct;
  col_ct = nz_col_ct;
  const uintptr_t table_size = 3 * col_ct;
  unsigned char* alloc_base;
  if (unlikely(pgl_malloc((table_size + 8) * sizeof(double) + (col_ct + 3 + 2 * table_size + table_size * table_size) * sizeof(int32_t) + table_size, &alloc_base))) {
    return S_CAST(double, NAN);
  }
  double* lfact_args_buf = R_CAST(double*, alloc_base);
  int32_t* col_tots = R_CAST(int32_t*, &(lfact_args_buf[table_size + 8]));
  int32_t* row_tots = &(col_tots[col_ct]);
  int32_t* obs_sorted = &(row_tots[3]);
  int32_t* cur = &(obs_sorted[table_size]);
  int32_t* saved_tables = &(cur[table_size]);
  unsigned char* flex_buf = R_CAST(unsigned char*, &(saved_tables[table_size * table_size]));
  for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
    row_tots[row_idx] = S_CAST(int32_t, row_tots_i64[row_idx]);
  }
  // Remove empty columns, and sort the rest in decreasing-total order.
  int32_t tot = 0;
  {
    uint32_t write_idx = 0;
    for (uint32_t col_idx = 0; write_idx != col_ct; ++col_idx) {
      const int32_t v0 = obs[col_idx];
      const int32_t v1 = obs[raw_col_ct + col_idx];
      const int32_t v2 = obs[2 * raw_col_ct + col_idx];
      const int32_t cur_col_tot = v0 + v1 + v2;
      if (!cur_col_tot) {
        continue;
      }
      uint32_t insert_idx = write_idx;
      for (; insert_idx; --insert_idx) {
        if (col_tots[insert_idx - 1] >= cur_col_tot) {
          break;
        }
        col_tots[insert_idx] = col_tots[insert_idx - 1];
        for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
          obs_sorted[row_idx * col_ct + insert_idx] = obs_sorted[row_idx * col_ct + insert_idx - 1];
        }
      }
      col_tots[insert_idx] = cur_col_tot;
      obs_sorted[insert_idx] = v0;
      obs_sorted[col_ct + insert_idx] = v1;
      obs_sorted[2 * col_ct + insert_idx] = v2;
      ++write_idx;
      tot += cur_col_tot;
    }
  }
  Fisher3xkCtx ctx;
  ctx.col_ct = col_ct;
  ctx.col_tots = col_tots;
  ctx.row_tots = row_tots;
  ctx.obs = obs_sorted;
  ctx.cur = cur;
  ctx.saved_tables = saved_tables;
  ctx.flex_buf = flex_buf;
  ctx.lfact_args_buf = lfact_args_buf;
  ctx.outer_sum = 0.0;
  ctx.tie_ct = 0.0;
  for (uintptr_t cell_idx = 0; cell_idx != table_size; ++cell_idx) {
    lfact_args_buf[cell_idx] = obs_sorted[cell_idx];
  }
  ctx.starting_lnprobv_ddr = ddr_negate(ddr_sort_and_add_lfacts(table_size, lfact_args_buf));
  ctx.starting_lnprobv_tdr.x[0] = ctx.starting_lnprobv_ddr.x[0];
  ctx.starting_lnprobv_tdr.x[1] = ctx.starting_lnprobv_ddr.x[1];
  ctx.starting_lnprobv_tdr.x[2] = DBL_MAX;
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    lfact_args_buf[col_idx] = col_tots[col_idx];
  }
  for (uint32_t row_idx = 0; row_idx != 3; ++row_idx) {
    lfact_args_buf[col_ct + row_idx] = row_tots[row_idx];
  }
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_sort_and_add_lfacts(col_ct + 3, lfact_args_buf),
            ddr_lfact(tot));
  const double starting_lnprob =
    ddr_add(lnprobf_ddr, ctx.starting_lnprobv_ddr).x[0];
  {
    // See Fisher2xkLnP().  The row-1 and row-2 entries of all columns but
    // one determine the table.
    double ln_table_ct = 0.0;
    double max_ln_range = 0.0;
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      const int32_t cur_col_tot = col_tots[col_idx];
      const double ln_range = log(S_CAST(double, MINV(cur_col_tot, row_tots[0]) + 1)) + log(S_CAST(double, MINV(cur_col_tot, row_tots[1]) + 1));
      ln_table_ct += ln_range;
      if (ln_range > max_ln_range) {
        max_ln_range = ln_range;
      }
    }
    double ln_p_bound = starting_lnprob + ln_table_ct - max_ln_range;
    ln_p_bound += fabs(ln_p_bound) * k2m35;
    if (ln_p_bound < ln_p_cutoff) {
      free(alloc_base);
      return MINV(ln_p_bound, 0.0);
    }
  }
  // Find the mode: start from (rounded-down) proportional allocation, fix up
  // the row-1 and row-2 totals, and then apply likelihood-increasing unit
  // cycles until there are none left.
  for (uint32_t row_idx = 0; row_idx != 2; ++row_idx) {
    int32_t row_sum = 0;
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      const int32_t val = S_CAST(int32_t, (S_CAST(int64_t, col_tots[col_idx]) * row_tots[row_idx]) / tot);
      cur[row_idx * col_ct + col_idx] = val;
      row_sum += val;
    }
    for (uint32_t col_idx = 0; row_sum != row_tots[row_idx]; ++col_idx) {
      if (cur[col_idx] + (row_idx? cur[col_ct + col_idx] : 0) < col_tots[col_idx]) {
        cur[row_idx * col_ct + col_idx] += 1;
        ++row_sum;
      }
    }
  }
  for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
    cur[2 * col_ct + col_idx] = col_tots[col_idx] - cur[col_idx] - cur[col_ct + col_idx];
  }
  while (Fisher3xkScanCycles(&ctx, INT32_MAX, 1)) {
  }
  const uint32_t status = Fisher3xkLevel(&ctx, 3 * (col_ct - 1));
  double result;
  if (status != kFisher2xkNonempty) {
    // Observed table is modal, so the non-mid p-value is 1.
    result = 0.0;
    if (midp) {
      result = log1p(-0.5 * Fisher3xkModeCt(&ctx, INT32_MAX) * exp(starting_lnprob));
    }
  } else {
    if (midp) {
      ctx.outer_sum -= ctx.tie_ct * 0.5;
    }
    result = log(ctx.outer_sum) + starting_lnprob;
    if (result > -k2m35) {
      result = 0.0;
    }
  }
  free(alloc_base);
  return result;
}

double Fisher33LnP(const int32_t* obs, uint32_t midp, double ln_p_cutoff) {
  return Fisher3xkLnP(obs, 3, midp, ln_p_cutoff);
}

#ifdef __cplusplus
}
#endif
