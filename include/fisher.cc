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

#include <assert.h>
#include <math.h>

#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

dd_real fisher22_ln_prob_internal(int64_t m11, int64_t m12, int64_t m21, int64_t m22) {
  dd_real ddrs[8];
  ddrs[0] = ddr_lfact(m11 + m12);
  ddrs[1] = ddr_lfact(m21 + m22);
  ddrs[2] = ddr_lfact(m11 + m21);
  ddrs[3] = ddr_lfact(m12 + m22);
  ddrs[4] = ddr_negate(ddr_lfact(m11));
  ddrs[5] = ddr_negate(ddr_lfact(m12));
  ddrs[6] = ddr_negate(ddr_lfact(m21));
  ddrs[7] = ddr_negate(ddr_lfact(m22));
  return ddr_sub(ddr_sort_and_add(8, ddrs), ddr_lfact(m11 + m12 + m21 + m22));
}

// Returns positive value if m22 = obs_m22 + m22_incr has higher probability
// than m22 = obs_m22, 0 if identical probability, and negative value if lower
// probability.
// If neg_numer_tdr has not been computed yet, set its x[0] to DBL_MAX; it will
// be filled in if necessary.
intptr_t Fisher22Compare(uint64_t obs_m11, uint64_t obs_m12, uint64_t obs_m21, uint64_t obs_m22, int64_t m22_incr, td_real* neg_numer_tdr_ptr, double* dbl_ptr) {
  // Fisher 2x2 probability is
  //
  //   (m11+m12)! (m21+m22)! (m11+m21)! (m12+m22)!
  //   -------------------------------------------
  //     m11! m12! m21! m22! (m11+m12+m21+m22)!
  //
  // so likelihood ratio of interest is
  //
  //           obs_m11! obs_m12! obs_m21! obs_m22!
  //   ---------------------------------------------------
  //   (obs_m11+j)! (obs_m12-j)! (obs_m21-j)! (obs_m22+j)!
  //
  // where j=m22_incr.
  //
  // Note that HWE kind of maps to this, via
  //   m11 := obs_hets*0.5
  //   m12 := obs_hom1
  //   m21 := obs_hom2
  //   m22 := (obs_hets-1)*0.5
  uint64_t numer_factorial_args[4];
  numer_factorial_args[0] = obs_m11;
  numer_factorial_args[1] = obs_m12;
  numer_factorial_args[2] = obs_m21;
  numer_factorial_args[3] = obs_m22;
  uint64_t denom_factorial_args[4];
  denom_factorial_args[0] = obs_m11 + m22_incr;
  denom_factorial_args[1] = obs_m12 - m22_incr;
  denom_factorial_args[2] = obs_m21 - m22_incr;
  denom_factorial_args[3] = obs_m22 + m22_incr;
  td_real ln_odds_ratio_tdr = tdr_make1(0.0);
  return CompareFactorialProducts(4, tdr_make1(1.0), 0, 0, numer_factorial_args, denom_factorial_args, neg_numer_tdr_ptr, &ln_odds_ratio_tdr, dbl_ptr);
}

static inline intptr_t Fisher22CompareDdr(uint64_t obs_m11, uint64_t obs_m12, uint64_t obs_m21, uint64_t obs_m22, int64_t m22_incr, dd_real neg_numer_ddr, double* dbl_ptr) {
  td_real neg_numer_tdr = {{neg_numer_ddr.x[0], neg_numer_ddr.x[1], DBL_MAX}};
  return Fisher22Compare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, &neg_numer_tdr, dbl_ptr);
}

// obs_m11 + obs_m12 + obs_m21 + obs_m22 assumed to be <2^31.
// They're defined as int32_ts instead of uint32_ts since signed int <->
// floating-point conversions are sometimes faster than the same-width unsigned
// int <-> floating-point conversions.  (Yes, uint32 -> float64 should be fine
// on 64-bit systems, but we may as well write this to also be efficient on
// 32-bit systems when there's no real drawback to doing so.)
//
// Not difficult to extend this to obs_m11 + obs_m12 + obs_m21 + obs_m22 <
// 2^52, if the rational-arithmetic backstop in CompareFactorialProducts() is
// revised to use td_reals.
//
// (Note that the odds-ratio and odds-ratio-confidence-interval reported by R
// fisher.test can also be calculated efficiently, using e.g. the approach in
// the R BiasedUrn package's meanFNCHypergeo() and pFNCHypergeo() functions.)
double Fisher22TwoSidedP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, int32_t midp, uint32_t logp) {
  // Normalize: m11 >= m22, m12 >= m21, m11*m22 < m12*m21.
  // Note that the first two are reversed from PLINK 1.9, to get rid of
  // spurious index differences between Fisher22 and Fisher23.
  if (obs_m11 < obs_m22) {
    swap_i32(&obs_m11, &obs_m22);
  }
  if (obs_m12 < obs_m21) {
    swap_i32(&obs_m12, &obs_m21);
  }
  if (S_CAST(int64_t, obs_m11) * obs_m22 > S_CAST(int64_t, obs_m12) * obs_m21) {
    swap_i32(&obs_m11, &obs_m12);
    swap_i32(&obs_m21, &obs_m22);
  }
  if (!midp) {
    // Fast path for p=1.
    if ((obs_m11 + 1LL) * (obs_m22 + 1) == S_CAST(int64_t, obs_m12) * obs_m21) {
      return logp? 0.0 : 1.0;
    }
  }
  // Iterate outward to floating-point precision limit.

  // I experimented with making more int32 -> double casts explicit, but
  // concluded that it was hurting readability more than it was helping.
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;
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
  if ((obs_m11 + 172LL) * (obs_m22 + 172LL) >= (obs_m12 - 172LL) * (obs_m21 - 172LL)) {
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
        const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
        td_real starting_lnprobv_tdr = tdr_make1(DBL_MAX);
        const intptr_t cmp_result = Fisher22Compare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, &starting_lnprobv_tdr, &lik);
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
  // todo: use quad-doubles when total > ~2^42
  dd_real starting_lnprobv_ddr =
    ddr_negate(ddr_sort_and_add_4_lfacts(obs_m11, obs_m12, obs_m21, obs_m22));

  // Now we want to jump near the other tail, without evaluating that many
  // contingency table log-likelihoods along the way.
  //
  // Each full log-likelihood evaluation requires 4 ddr_lfact() calls.  Since
  // they are now performed with extra precision, they require hundreds of
  // floating-point operations, so we want to limit ourselves to 1-2 full
  // evaluations most of the time.  (Possible todo: use lower-accuracy Lfact()
  // to jump around, followed by ddr_lfact() when exiting the loop.  Should be
  // an easy performance win, but there's a complexity cost so I'll wait until
  // I see a scenario where this branch executes frequently...)
  //
  // The current heuristic starts by reflecting (obs_m21 + m21) * 0.5 across
  // the mode, performing a full log-likelihood check at an adjacent valid
  // point.  (It is convenient to focus on m21 here, since m21=0 corresponds to
  // the outermost table on this tail.)  Hopefully we find that we're in
  // (starting_lnprob - 62 * kLn2, starting_lnprob], so we're at or near a
  // table that actually contributes to the tail-sum.  (This window is chosen
  // to be wide enough to guarantee that at least one point falls inside when
  // obs_m11 + obs_m12 + obs_m21 + obs_m22 < 2^31.)
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
    m21 = S_CAST(int32_t, m21);
    if (m21 < 0) {
      m21 = 0;
    }
  }
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_sort_and_add_4_lfacts(m1x, m2x, mx1, obs_m12 + obs_m22),
            ddr_lfact(mxx));
  const dd_real starting_lnprob_ddr = ddr_add(lnprobf_ddr, starting_lnprobv_ddr);
  while (1) {
    m11 = mx1 - m21;
    m12 = m12_minus_m21 + m21;
    m22 = m2x - m21;
    const dd_real lnprobv_ddr =
      ddr_negate(ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
    const double lnprobv_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    // Could tighten this threshold further; I haven't performed a careful
    // error analysis yet but CompareFactorialProducts() includes a plausible
    // assumption that 2^{-60} is safe.  But code is correct as long as we're
    // guaranteed to enter the "lik < 2 - one_minus_scaled_eps" branch for
    // positive lnprobv_diff.
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
    } else if (lnprobv_diff > -62 * kLn2) {
      lik = exp(lnprobv_diff);
      break;
    } else {
      const double ll_deriv = log((m12 + 1) * (m21 + 1) / (m11 * m22));
      // Round down, to guarantee we don't overshoot.
      // We're guaranteed to make progress, since lnprobv_diff <= -62 * log(2),
      // m11 * m22 < 2^62, and (m12 + 1) * (m21 + 1) >= 1.
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
    const intptr_t cmp_result = Fisher22CompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, starting_lnprobv_ddr, &lik);
    if (cmp_result <= 0) {
      tail_sum += lik;
      if (midp && (cmp_result == 0)) {
        tail_sum -= 0.5;
      }
    }
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
    if (left_sum == INFINITY) {
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
  const dd_real starting_lnprob_ddr = fisher22_ln_prob_internal(obs_m11, obs_m12, obs_m21, obs_m22);
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
// <0.6 ULP relative error unless problem instance is huge (sum > ~2^40 if logp
// true, ~2^37 if false; see binom_ln_prob_internal error analysis).
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
    // We want to compute left_sum_ddr to at most 2^{-57} relative error.  See
    // related discussion in Pbinom() implementation.
    const double min_incr_left = 0.03125 / (m22 * m22);
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
    if (!(left_sum_ddr.x[0] < INFINITY)) {
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
      if (left_sum_ddr.x[0] == INFINITY) {
        return logp? 0.0 : 1.0;
      }
    }

    // Now compute the right-sum to at most 2^{-57} relative error.
    dd_real right_sum_ddr = ddr_muld(first_right_mult_ddr, start_lik);
    m11 = obs_m11 + 1;
    m12 = obs_m12 - 1;
    m21 = obs_m21 - 1;
    m22 = obs_m22 + 1;
    lik_ddr = right_sum_ddr;
    if (m21 > 0) {
      const double min_incr_right = 0.03125 / (m21 * m21);
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
    if (!(denom_ddr.x[0] < INFINITY)) {
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
  // log-likelihood (this is the only step that could introduce >0.1 ULP error,
  // and then only for huge n), and accumulate the tail-sum from there.
  // (todo: opportunistically use the _lfact() path when that's within the
  // error budget and rates to be faster than the left_sum / (left_sum +
  // right_sum) approach.  Reasonable to wait until td_real functions are added
  // to this library, though.)
  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx2 = obs_m12 + obs_m22;
  const double mxx = m1x + m2x;
  const double modal_m22 = m2x * mx2 / mxx;
  if (modal_m22 <= m22 + 168) {
    dd_real lik_ddr = ddr_maked(1.0);
    dd_real right_sum_ddr = ddr_maked(0.0);
    if (m21 > 0) {
      const double min_incr_right = 0.03125 / (m21 * m21);
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
      const double min_incr_left = 0.03125 / (m22 * m22);
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
  dd_real ln_prob_ddr = fisher22_ln_prob_internal(obs_m11, obs_m12, obs_m21, obs_m22);
  dd_real lik_ddr = ddr_maked(1.0);
  dd_real left_sum_ddr = lik_ddr;
  if (m22 > 0) {
    const double min_incr_left = 0.03125 / (m22 * m22);
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
    return abcd + final_return_incr;
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
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_sort_and_add_4_lfacts(m1x, m2x, mxx - mx2, mx2),
            ddr_lfact(mxx));
  // guessing it's worthwhile to with lnprob instead of lnprobv since that
  // should keep x0_coeff's magnitude smaller
  const dd_real mode_lnprob_ddr = ddr_sub(lnprobf_ddr,
                                          ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
  const dd_real modem1_lnprob_incr_ddr = ddr_log(ddr_accurate_div(ddr_mul2d(m11, m22), ddr_mul2d(m12 + 1, m21 + 1)));
  const dd_real modep1_lnprob_incr_ddr = ddr_log(ddr_accurate_div(ddr_mul2d(m12, m21), ddr_mul2d(m11 + 1, m22 + 1)));

  const double x2_coeff = 0.5 * ddr_add(modem1_lnprob_incr_ddr, modep1_lnprob_incr_ddr).x[0];
  const double x1_coeff = 0.5 * ddr_sub(modep1_lnprob_incr_ddr, modem1_lnprob_incr_ddr).x[0];
  // mode^2 * x2_coeff + mode * x1_coeff + x0_coeff = mode_lnprob_ddr
  const double x0_coeff = ddr_subd(mode_lnprob_ddr, m22 * prefer_fma(m22, x2_coeff, x1_coeff)).x[0];

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
    m22 = trunc((sqrt_discrim - x1_coeff) / (2 * x2_coeff));
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
    cur_lnprob_ddr = ddr_sub(lnprobf_ddr,
                             ddr_sort_and_add_4_lfacts(m11, m12, m21, m22));
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
    const double min_incr_left = 0.125 / (m22 * m22);
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
          const intptr_t cmp_result = Fisher22CompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, *starting_lnprobv_ddr_ptr, &lik);
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
        const intptr_t cmp_result = Fisher22CompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, *starting_lnprobv_ddr_ptr, &lik);
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
        const intptr_t cmp_result = Fisher22CompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, *starting_lnprobv_ddr_ptr, &lik);
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
      const intptr_t cmp_result = Fisher22CompareDdr(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, *starting_lnprobv_ddr_ptr, &lik);
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
        if (l22) {
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
        if (l11) {
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
        if (r21) {
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
        if (r12) {
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
